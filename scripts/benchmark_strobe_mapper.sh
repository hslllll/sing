#!/bin/bash
set -e

SING_BIN="./target/release/sing"
STROBEALIGN_BIN="strobealign"
MAPPER_BIN="x-mapper" 

MODE=$1
if [ "$MODE" == "h" ]; then
    echo "Mode: Human Genome"
    REF="GRCh38_latest_genomic.fna.gz"
    REF_DECOMP="GRCh38_latest_genomic.fna"
    REF_URL="https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz"
    INDEX_PREFIX="human"
elif [ "$MODE" == "y" ]; then
    echo "Mode: Yeast Genome"
    REF="Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
    REF_DECOMP="Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
    REF_URL="http://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
    INDEX_PREFIX="yeast"
elif [ "$MODE" == "a" ]; then
    echo "Mode: Arabidopsis thaliana (TAIR10)"
    REF="Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
    REF_DECOMP="Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
    REF_URL="http://ftp.ensemblgenomes.org/pub/plants/release-57/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
    INDEX_PREFIX="arabidopsis"
else
    echo "Usage: $0 [h|y|a]"
    echo "  h: Human (GRCh38)"
    echo "  y: Yeast (Saccharomyces cerevisiae)"
    echo "  a: Arabidopsis thaliana (TAIR10)"
    exit 1
fi

COVERAGE_LIST=(1 5 10)
MUT_RATES=(0.001 0.01)
THREADS=8

OUTPUT_CSV="benchmark_results_custom.csv"

check_tool() {
    if ! command -v $1 &> /dev/null; then
        if [ -f "$1" ]; then
            return 0
        fi
        echo "Error: Tool '$1' not found. Please install it or set path."
        exit 1
    fi
}
check_tool dwgsim
check_tool /usr/bin/time

echo "Checking Mappers..."
if [ ! -f "$SING_BIN" ]; then
    echo "Building Sing..."
    cargo build --release
fi

if ! command -v $STROBEALIGN_BIN &> /dev/null; then
    echo "Warning: strobealign not in PATH. Assuming ./strobealign or failing."
fi
if ! command -v $MAPPER_BIN &> /dev/null; then
    echo "Warning: X-mapper ('mapper') not in PATH. Assuming ./mapper or failing."
fi

echo "=== Preparing Data & Indices ==="

if [ ! -f "$REF" ]; then
    echo "Downloading Genome..."
    wget -c -O "$REF" "$REF_URL"
fi

if [ ! -f "$REF_DECOMP" ]; then
    gunzip -c "$REF" > "$REF_DECOMP"
fi

echo "--- Indexing Phase ---"

if [ ! -f "${INDEX_PREFIX}.idx" ]; then
    echo "[Sing] Building Index..."
    $SING_BIN build "$REF_DECOMP" "${INDEX_PREFIX}.idx"
fi

if [ ! -f "${REF_DECOMP}.sti" ] && [ ! -f "${REF_DECOMP}.strobealign_index" ]; then
    echo "[Strobealign] Building Index..."
    $STROBEALIGN_BIN --create-index "$REF_DECOMP" -r 150
fi

cat << 'EOF' > analyze_benchmark.py
import sys
import re
import math
import os

TOLERANCE = 10

def parse_dwgsim_qname(qname):
    clean_qname = qname.split('/')[0]
    match = re.match(r'^(.*?)_(\d+)_(\d+)_([01])_([01])_', clean_qname)
    if match:
        return clean_qname, match.group(1), int(match.group(2)), int(match.group(3))
    return None

def load_sam_and_calc_stats(filepath):
    tp = 0
    fp = 0
    fn = 0
    
    if not os.path.exists(filepath):
        return 0, 0, 0, 0, 0, 0

    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('@'): continue
                parts = line.split('\t')
                qname = parts[0]
                flag = int(parts[1])
                rname = parts[2]
                pos = int(parts[3])
                
                if (flag & 256) or (flag & 2048): continue
                read_num = 2 if (flag & 128) else 1
                
                parsed = parse_dwgsim_qname(qname)
                if not parsed: continue
                _, true_ref, t_p1, t_p2 = parsed
                true_pos = t_p1 if read_num == 1 else t_p2
                
                if (flag & 4):
                    fn += 1
                else:
                    if rname == true_ref and abs(pos - true_pos) <= TOLERANCE:
                        tp += 1
                    else:
                        fp += 1
    except Exception:
        pass

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn + fp) if (tp + fn + fp) > 0 else 0.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
    return tp, fp, fn, precision * 100, recall * 100, f1 * 100

def main():
    exp_name = sys.argv[1]
    
    tools = [
        ("Sing",        sys.argv[2], sys.argv[3], "sing.sam"),
        ("Strobealign", sys.argv[4], sys.argv[5], "strobealign.sam"),
        ("X-mapper",    sys.argv[6], sys.argv[7], "mapper.sam"),
    ]

    
    for name, time_val, mem_val, sam_file in tools:
        tp, fp, fn, prec, rec, f1 = load_sam_and_calc_stats(sam_file)
        
        status = "OK"
        if time_val == "N/A": status = "FAIL"
        
        sys.stderr.write(f"{name:<12} | Time: {str(time_val):<6}s | Mem: {str(mem_val):<8}kb | F1: {f1:5.2f}%\n")
        
        print(f"{exp_name},{name},{time_val},{mem_val},{prec:.4f},{rec:.4f},{f1:.4f}")

if __name__ == "__main__":
    main()
EOF

echo "Experiment,Tool,Time_sec,Mem_kb,Precision,Recall,F1" > "$OUTPUT_CSV"

echo "=== Starting Benchmarks ==="
echo "Writing results to $OUTPUT_CSV"

for COV in "${COVERAGE_LIST[@]}"; do
    for MUT in "${MUT_RATES[@]}"; do
        
        EXP_ID="Cov_${COV}_Mut_${MUT}"
        echo "------------------------------------------------"
        echo "Generating Reads: Coverage=$COV, Mutation=$MUT"
        
        sim_prefix="sim_${EXP_ID}"
        dwgsim -1 150 -2 150 -y 0 -r $MUT -C $COV "$REF_DECOMP" "$sim_prefix" > /dev/null 2>&1
        
        R1="${sim_prefix}.bwa.read1.fastq.gz"
        R2="${sim_prefix}.bwa.read2.fastq.gz"

        
        measure() {
            local LOG=$1
            shift
            
            /usr/bin/time -f "%e %M" -o "$LOG" "$@"
            
            if [ $? -eq 0 ]; then
                read T M < "$LOG"
                RET_TIME=$T
                RET_MEM=$M
            else
                RET_TIME="N/A"
                RET_MEM="N/A"
            fi
        }

        echo "Running Sing..."
        if measure "sing.log" $SING_BIN map -t $THREADS "${INDEX_PREFIX}.idx" -1 "$R1" -2 "$R2" -o sing.sam; then
            TIME_SING=$RET_TIME
            MEM_SING=$RET_MEM
        else
            TIME_SING="N/A"
            MEM_SING="N/A"
            echo "Sing Failed"
        fi

        echo "Running Strobealign..."
        CMD="$STROBEALIGN_BIN -t $THREADS $REF_DECOMP $R1 $R2 > strobealign.sam"
        if /usr/bin/time -f "%e %M" -o strobe.log bash -c "$CMD"; then
             read T M < strobe.log
             TIME_STROBE=$T
             MEM_STROBE=$M
        else
             TIME_STROBE="N/A"
             MEM_STROBE="N/A"
             echo "Strobealign Failed"
        fi

        echo "Running X-mapper..."
        CMD="$MAPPER_BIN --reference $REF_DECOMP --queries $R1 --queries $R2 --num-threads $THREADS --out-sam mapper.sam"
        if /usr/bin/time -f "%e %M" -o mapper.log bash -c "$CMD"; then
             read T M < mapper.log
             TIME_MAP=$T
             MEM_MAP=$M
        else
             TIME_MAP="N/A"
             MEM_MAP="N/A"
             echo "X-Mapper Failed"
        fi
        
        python3 analyze_benchmark.py "$EXP_ID" \
            "$TIME_SING" "$MEM_SING" \
            "$TIME_STROBE" "$MEM_STROBE" \
            "$TIME_MAP" "$MEM_MAP" >> "$OUTPUT_CSV"

        rm sing.sam strobealign.sam mapper.sam sing.log strobe.log mapper.log
        rm "${sim_prefix}"*
    done
done

echo "Done. Results saved in $OUTPUT_CSV"
