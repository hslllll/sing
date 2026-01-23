#!/bin/bash
set -e
export PATH=$PATH:/home/hyunsu.lim/programs/columba/build_Vanilla/
export PATH=$PATH:/home/hyunsu/Code/columba/build_Vanilla/
echo "=== Setup ==="
if ! command -v minimap2 &> /dev/null;
then
    echo "minimap2 not found. Please install it (e.g. apt install minimap2 or brew install minimap2)"
    HAS_MINIMAP=0
else
    HAS_MINIMAP=1
fi

if ! command -v samtools &> /dev/null;
then
    echo "samtools not found. Sorting will be skipped."
fi

if ! command -v bwa-mem2 &> /dev/null;
then
    echo "bwa-mem2 not found."
    HAS_BWA=0
else
    HAS_BWA=1
fi

if ! command -v columba &> /dev/null;
then
    echo "columba not found."
    HAS_COLUMBA=0
else
    HAS_COLUMBA=1
fi

if ! command -v bowtie2 &> /dev/null;
then
    echo "bowtie2 not found."
    HAS_BT2=0
else
    HAS_BT2=1
fi

echo "=== Generating Data ==="
REF="ref.fa"
R1="r1.fq"
R2="r2.fq"

python3 -c '''

import random
import sys

num_refs = 5
ref_len = 1_000_000
refs = []
ref_names = []

with open("ref.fa", "w") as f:
    for r in range(num_refs):
        seq = "".join(random.choices("ACGT", k=ref_len))
        name = f"ref{r+1}"
        refs.append(seq)
        ref_names.append(name)
        f.write(f">{name}\n{seq}\n")

n_reads = 10000
read_len = 150
fragment_len = 300
mutation_rate = 0.01

def reverse_complement(seq):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(complement.get(base, base) for base in reversed(seq))

def mutate(seq):
    bases = list(seq)
    for i in range(len(bases)):
        if random.random() < mutation_rate:
            bases[i] = random.choice("ACGT")
    return "".join(bases)

with open("r1.fq", "w") as f1, open("r2.fq", "w") as f2, open("truth.txt", "w") as ft:
    for i in range(n_reads):
        ref_idx = random.randint(0, num_refs - 1)
        ref = refs[ref_idx]
        ref_name = ref_names[ref_idx]
        
        start = random.randint(0, ref_len - fragment_len)
        is_reverse = random.choice([True, False])

        seq_left = ref[start : start + read_len]
        seq_right = ref[start + fragment_len - read_len : start + fragment_len]

        if not is_reverse:
            r1_seq = seq_left
            r2_seq = reverse_complement(seq_right)
        else:
            r1_seq = reverse_complement(seq_right)
            r2_seq = seq_left

        r1_seq = mutate(r1_seq)
        r2_seq = mutate(r2_seq)

        qual = "I" * read_len
        f1.write(f"@read_{i}/1\n{r1_seq}\n+\n{qual}\n")
        f2.write(f"@read_{i}/2\n{r2_seq}\n+\n{qual}\n")
        
        if not is_reverse:
             true_pos = start + 1
        else:
             true_pos = start + fragment_len - read_len + 1

        ft.write(f"read_{i}\t{ref_name}\t{true_pos}\n")
'''

echo "Generated ref.fa (5MB), r1.fq, r2.fq, and truth.txt"

echo "=== Building ==="
if command -v cargo &> /dev/null; then
    cargo build --release
else
    echo "cargo not found. Skipping Rust build. Please install Rust/cargo if you want to build Sing."
fi

echo "=== Benchmarking ==="

echo "Running Sing (CPU)..."
./target/release/sing build $REF reference.idx

start_time=$(date +%s%N)
if ./target/release/sing map -t 8 reference.idx -1 $R1 -2 $R2 -o sing.sam; then
    end_time=$(date +%s%N)
    elapsed_sing=$(( (end_time - start_time) / 1000000 ))
    echo "Sing: ${elapsed_sing} ms"
else
    elapsed_sing="N/A"
    echo "Sing: Failed"
fi

if [ $HAS_MINIMAP -eq 1 ]; then
    echo "Running Minimap2..."
    start_time=$(date +%s%N)
    if minimap2 -t 8 -ax sr $REF $R1 $R2 > minimap.sam 2>/dev/null; then
        end_time=$(date +%s%N)
        elapsed_mm=$(( (end_time - start_time) / 1000000 ))
        echo "Minimap2: ${elapsed_mm} ms"
    else
        elapsed_mm="N/A"
        echo "Minimap2: Failed"
    fi
else
    elapsed_mm="N/A"
fi

if [ $HAS_COLUMBA -eq 1 ]; then
    echo "Running Columba..."
    if [ ! -f "columba_index.rev.bwt" ]; then
        columba_build -r columba_index -f $REF > /dev/null 2>&1
    fi
    
    start_time=$(date +%s%N)
    if columba -t 8 -r columba_index -f $R1 -F $R2 -o columba.sam > /dev/null 2>&1; then
        end_time=$(date +%s%N)
        elapsed_columba=$(( (end_time - start_time) / 1000000 ))
        echo "Columba: ${elapsed_columba} ms"
    else
        elapsed_columba="N/A"
        echo "Columba: Failed"
    fi
else
    elapsed_columba="N/A"
fi

if [ $HAS_BWA -eq 1 ]; then
    echo "Running bwa-mem2..."
    bwa-mem2 index $REF > /dev/null 2>&1
    
    start_time=$(date +%s%N)
    if bwa-mem2 mem -t 8 $REF $R1 $R2 > bwa.sam 2>/dev/null; then
        end_time=$(date +%s%N)
        elapsed_bwa=$(( (end_time - start_time) / 1000000 ))
        echo "bwa-mem2: ${elapsed_bwa} ms"
    else
        elapsed_bwa="N/A"
        echo "bwa-mem2: Failed"
    fi
else
    elapsed_bwa="N/A"
fi

if [ $HAS_BT2 -eq 1 ]; then
    echo "Running Bowtie2..."
    bowtie2-build $REF $REF > /dev/null 2>&1
    
    start_time=$(date +%s%N)
    if bowtie2 -p 8 -x $REF -1 $R1 -2 $R2 > bowtie2.sam 2>/dev/null; then
        end_time=$(date +%s%N)
        elapsed_bt2=$(( (end_time - start_time) / 1000000 ))
        echo "Bowtie2: ${elapsed_bt2} ms"
    else
        elapsed_bt2="N/A"
        echo "Bowtie2: Failed"
    fi
else
    elapsed_bt2="N/A"
fi


echo "=== Evaluation ==="

python3 -c '''
import sys
import os

def parse_sam(sam_file, truth_dict):
    mapped = 0
    correct = 0
    
    if not os.path.exists(sam_file):
        return -1, -1

    with open(sam_file) as f:
        for line in f:
            if line.startswith("@"): continue
            parts = line.split("\t")
            if len(parts) < 4: continue
            
            qname = parts[0]
            flag = int(parts[1])
            rname = parts[2]
            pos = int(parts[3])
            
            if (flag & 4): continue
            
            if (flag & 256) or (flag & 2048): continue
            
            if (flag & 1) and not (flag & 64):
                continue

            if qname.endswith("/1") or qname.endswith("/2"):
                qname = qname[:-2]

            if qname in truth_dict:
                t_ref, t_pos = truth_dict[qname]
                if rname == t_ref and abs(pos - t_pos) < 30:
                    correct += 1
                mapped += 1
                    
    return mapped, correct

truth_file = "truth.txt"
tools = [
    ("Minimap2", "minimap.sam"),
    ("BWA-MEM2", "bwa.sam"),
    ("Columba", "columba.sam"),
    ("Sing", "sing.sam"),
    ("Bowtie2", "bowtie2.sam")
]

truth = {}
with open(truth_file) as f:
    for line in f:
        p = line.strip().split()
        if len(p) >= 3:
            truth[p[0]] = (p[1], int(p[2]))

total_pairs = len(truth)

print(f"{str(total_pairs)} pairs in ground truth.")
header = "{:<15} {:<10} {:<10} {:<10} {:<10} {:<10}".format("Tool", "Mapped", "Correct", "Recall", "Precision", "F1")
print(f"{chr(27)}[1m{header}{chr(27)}[0m")
print("-" * 75)

for name, path in tools:
    m, c = parse_sam(path, truth)
    if m == -1:
        na = "N/A"
        print(f"{name:<15} {na:<10} {na:<10} {na:<10} {na:<10} {na:<10}")
    else:
        recall = c / total_pairs * 100 if total_pairs > 0 else 0
        precision = c / m * 100 if m > 0 else 0
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        print(f"{name:<15} {str(m):<10} {str(c):<10} {recall:.2f}%     {precision:.2f}%     {f1:.2f}")
'''

echo ""
echo "=== Execution Time ==="
echo "Minimap2: ${elapsed_mm} ms"
echo "BWA-MEM2: ${elapsed_bwa} ms"
echo "Columba:  ${elapsed_columba} ms"
echo "Sing:     ${elapsed_sing} ms"
echo "Bowtie2:  ${elapsed_bt2} ms"


rm -f minimap_mapped.txt sing_mapped.txt overlap.txt fn.txt fp.txt ref.fa.*
