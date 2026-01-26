#!/bin/bash
set -e

COLUMBA_PATH="/home/hyunsu.lim/programs/columba/build_Vanilla"
export PATH=$PATH:/home/hyunsu/Code/columba/build_Vanilla/
export PATH=$PATH:$COLUMBA_PATH

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
elif [ "$MODE" == "m" ]; then
    echo "Mode: Maize (Zea mays)"
    REF="Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz"
    REF_DECOMP="Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa"
    REF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz"
    INDEX_PREFIX="maize"
else
    echo "Usage: $0 [h|y|a|m]"
    echo "  h: Human (GRCh38)"
    echo "  y: Yeast (Saccharomyces cerevisiae)"
    echo "  a: Arabidopsis thaliana (TAIR10)"
    echo "  m: Maize (Zea mays)"
    exit 1
fi

COVERAGE_LIST=(0.01 0.1)
MUT_RATES=(0.001 0.01)

OUTPUT_CSV="benchmark_results_f1.${MODE}.csv"

check_tool() {
    if ! command -v $1 &> /dev/null; then
        echo "Error: $1 not found. Please install it."
        exit 1
    fi
}
check_tool dwgsim
check_tool minimap2
check_tool bowtie2
check_tool bwa-mem2
check_tool samtools

echo "=== Preparing Data & Indices ==="
if command -v cargo &> /dev/null; then
    echo "Building Sing..."
    cargo build --release
fi

if [ ! -f "$REF" ]; then
    echo "Downloading Genome..."
    wget -c -O "$REF" "$REF_URL"
fi

if [ ! -f "$REF_DECOMP" ]; then
    gunzip -c "$REF" > "$REF_DECOMP"
fi

if [ ! -f "$REF_DECOMP.0123" ]; then
    echo "Indexing BWA-MEM2..."
    bwa-mem2 index "$REF_DECOMP" > /dev/null 2>&1
fi

if [ ! -f "$REF_DECOMP.1.bt2" ]; then
    echo "Indexing Bowtie2..."
    bowtie2-build -t 8 "$REF_DECOMP" "$REF_DECOMP" > /dev/null 2>&1
fi

if [ ! -f "$REF_DECOMP.bwt" ]; then
    echo "Indexing Columba..."
    columba_build -r "$REF_DECOMP" -f "$REF_DECOMP" > /dev/null 2>&1
fi

if [ ! -f "${INDEX_PREFIX}.idx" ]; then
    echo "Indexing Sing..."
    if command -v cargo &> /dev/null; then
        cargo build --release
    fi
    ./target/release/sing build "$REF_DECOMP" "${INDEX_PREFIX}.idx"
fi

cat << 'EOF' > analyze_benchmark.py
import sys
import re
import math

TOLERANCE = 15

def parse_dwgsim_qname(qname):
    """
    Parses dwgsim read name to extract ground truth.
    Format: @RefName_Start1_Start2_Strand1_Strand2_...
    Returns: (ref, pos1, pos2) where positions are 1-based.
    """
    clean_qname = qname.split('/')[0]
    
    match = re.match(r'^(.*?)_(\d+)_(\d+)_([01])_([01])_', clean_qname)
    if match:
        ref = match.group(1)
        p1 = int(match.group(2))
        p2 = int(match.group(3))
        return clean_qname, ref, p1, p2
    return None

def load_sam_and_build_truth(filepaths):
    """
    Reads multiple SAM files to:
    1. Build a comprehensive Truth dictionary from read names.
    2. Store mapping results for each tool.
    
    Returns:
        truth_dict: { (qname_base, read_num): (true_ref, true_pos) }
        tool_results: { tool_name: { (qname_base, read_num): (mapped_ref, mapped_pos) } }
    """
    truth_dict = {}
    tool_results = {name: {} for name in filepaths}

    for tool_name, filepath in filepaths.items():
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    if line.startswith('@'): continue
                    parts = line.split('\t')
                    
                    qname = parts[0]
                    flag = int(parts[1])
                    rname = parts[2]
                    pos = int(parts[3])

                    read_num = 2 if (flag & 128) else 1
                    
                    parsed = parse_dwgsim_qname(qname)
                    if parsed:
                        q_base, true_ref, t_p1, t_p2 = parsed
                        
                        if (q_base, 1) not in truth_dict:
                            truth_dict[(q_base, 1)] = (true_ref, t_p1)
                        if (q_base, 2) not in truth_dict:
                            truth_dict[(q_base, 2)] = (true_ref, t_p2)

                    if (flag & 256) or (flag & 2048):
                        continue
                    
                    q_base = qname.split('/')[0]
                    
                    if (flag & 4):
                        tool_results[tool_name][(q_base, read_num)] = None
                    else:
                        tool_results[tool_name][(q_base, read_num)] = (rname, pos)
                        
        except FileNotFoundError:
            pass

    return truth_dict, tool_results

def calculate_metrics(truth_dict, tool_map):
    """
    Calculates TP, FP, FN, Precision, Recall, F1.
    """
    TP = 0
    FP = 0
    FN = 0
    
    for key, (true_ref, true_pos) in truth_dict.items():
        mapped_val = tool_map.get(key)
        
        if mapped_val is None:
            FN += 1
        else:
            pred_ref, pred_pos = mapped_val
            
            if pred_ref == true_ref and abs(pred_pos - true_pos) <= TOLERANCE:
                TP += 1
            else:
                FP += 1

    total_positives = len(truth_dict)
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall_fraction = TP / total_positives if total_positives > 0 else 0.0
    recall = recall_fraction * 100
    f1 = 2 * (precision * recall_fraction) / (precision + recall_fraction) if (precision + recall_fraction) > 0 else 0.0
    
    return TP, FP, FN, precision * 100, recall, f1 * 100

def main():
    exp_name = sys.argv[1]
    times = {
        'Minimap2': sys.argv[2],
        'BWA-MEM2': sys.argv[3],
        'Columba':  sys.argv[4],
        'Sing':     sys.argv[5],
        'Bowtie2':  sys.argv[6]
    }
    mode = sys.argv[7]

    file_paths = {
        'Minimap2': f'minimap.{mode}.sam',
        'BWA-MEM2': f'bwa.{mode}.sam',
        'Columba':  f'columba.{mode}.sam',
        'Sing':     f'sing.{mode}.sam',
        'Bowtie2':  f'bowtie2.{mode}.sam'
    }

    truth_dict, tool_results = load_sam_and_build_truth(file_paths)
    total_reads = len(truth_dict)

    print(f"\n--- Summary: {exp_name} (Mode: {mode}) ---")
    print("Tool        | Time_ms  | Precision | Recall   | F1        | TP       | FP       | FN")

    tools_ordered = ['Minimap2', 'BWA-MEM2', 'Columba', 'Sing', 'Bowtie2']

    for tool in tools_ordered:
        time_val = times.get(tool, "N/A")
        
        if tool in tool_results and len(tool_results[tool]) > 0:
            tp, fp, fn, prec, rec, f1 = calculate_metrics(truth_dict, tool_results[tool])
            
            print(f"{tool:<10} | {time_val:<8} | {prec:6.2f}%   | {rec:6.2f}%   | {f1:6.2f}    | {tp:<8} | {fp:<8} | {fn:<8}")
            
            sys.stderr.write(f"{exp_name},{tool},{time_val},{total_reads},{tp},{fp},{fn},{prec:.4f},{rec:.4f},{f1:.4f},{mode}\n")
        else:
            print(f"{tool:<10} | {time_val:<8} | {'FAILED':<9} | {'-':<9} | {'-':<9} | -        | -        | -")
            sys.stderr.write(f"{exp_name},{tool},{time_val},{total_reads},0,0,0,0,0,0,{mode}\n")

if __name__ == "__main__":
    main()
EOF

if [ ! -f "$OUTPUT_CSV" ]; then
    echo "Exp_ID,Tool,Time_ms,TotalReads,TP,FP,FN,Precision,Recall,F1,Mode" > "$OUTPUT_CSV"
fi

for COVERAGE in "${COVERAGE_LIST[@]}"; do
    for MUT_RATE in "${MUT_RATES[@]}"; do
        EXP_ID="Cov_${COVERAGE}_Mut_${MUT_RATE}.${MODE}"
        echo "================================================="
        echo "Running: $EXP_ID (Cov: $COVERAGE, Mut: $MUT_RATE)"
        echo "================================================="

        R1="sim_${EXP_ID}.bwa.read1.fastq.gz"
        R2="sim_${EXP_ID}.bwa.read2.fastq.gz"

        if [ -f "$R1" ] && [ -f "$R2" ]; then
            echo "Simulated reads exist. Skipping dwgsim."
        else
            dwgsim -C "$COVERAGE" -1 150 -2 150 -R 0 -X 0 -r "$MUT_RATE" -y 0 -H "$REF_DECOMP" "sim_${EXP_ID}" > /dev/null 2>&1
        fi


        echo "1. Running Sing..."
        START=$(date +%s%N)
        if ./target/release/sing map -t 8 "${INDEX_PREFIX}.idx" -1 "$R1" -2 "$R2" -o "sing.${MODE}.sam" > /dev/null 2>&1; then
            END=$(date +%s%N)
            TIME_SING=$(( (END - START) / 1000000 ))
        else
            TIME_SING="N/A"
            echo "Sing Failed"
        fi

        echo "2. Running Minimap2..."
        START=$(date +%s%N)
        if minimap2 -t 8 -ax sr "$REF_DECOMP" "$R1" "$R2" > "minimap.${MODE}.sam" 2>/dev/null; then
            END=$(date +%s%N)
            TIME_MM=$(( (END - START) / 1000000 ))
        else
            TIME_MM="N/A"
            echo "Minimap2 Failed"
        fi

        echo "3. Running Columba..."
        START=$(date +%s%N)
        if command -v columba &> /dev/null; then
             if columba -t 8 -r "$REF_DECOMP" -f "$R1" -F "$R2" -o "columba.${MODE}.sam" > /dev/null 2>&1; then
                 END=$(date +%s%N)
                 TIME_COL=$(( (END - START) / 1000000 ))
             else
                 TIME_COL="N/A"
                 echo "Columba Failed"
             fi
        else
            TIME_COL="N/A"
            touch "columba.${MODE}.sam"
        fi
        
        echo "4. Running BWA-MEM2..."
        START=$(date +%s%N)
        if bwa-mem2 mem -t 8 "$REF_DECOMP" "$R1" "$R2" > "bwa.${MODE}.sam" 2>/dev/null; then
            END=$(date +%s%N)
            TIME_BWA=$(( (END - START) / 1000000 ))
        else
            TIME_BWA="N/A"
            echo "BWA Failed"
        fi

        echo "5. Running Bowtie2..."
        START=$(date +%s%N)
        if bowtie2 -p 8 -x "$REF_DECOMP" -1 "$R1" -2 "$R2" > "bowtie2.${MODE}.sam" 2>/dev/null; then
            END=$(date +%s%N)
            TIME_BT2=$(( (END - START) / 1000000 ))
        else
            TIME_BT2="N/A"
            echo "bowtie2: Failed"
        fi

        python3 analyze_benchmark.py "$EXP_ID" "$TIME_MM" "$TIME_BWA" "$TIME_COL" "$TIME_SING" "$TIME_BT2" "$MODE" | tee -a "$OUTPUT_CSV"

        rm "minimap.${MODE}.sam" "bwa.${MODE}.sam" "columba.${MODE}.sam" "sing.${MODE}.sam" "bowtie2.${MODE}.sam"
        echo ""
    done
done

echo "Benchmark Complete. Results saved in $OUTPUT_CSV"
