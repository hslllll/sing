#!/bin/bash
set -e

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
elif [ "$MODE" == "b" ]; then
    echo "Mode: Brassica rapa"
    REF="Brassica_rapa.Brapa_1.0.dna.toplevel.fa.gz"
    REF_DECOMP="Brassica_rapa.Brapa_1.0.dna.toplevel.fa"
    REF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/brassica_rapa/dna/Brassica_rapa.Brapa_1.0.dna.toplevel.fa.gz"
    INDEX_PREFIX="brassica"
else
    echo "Usage: $0 [h|y|a|m|b]"
    echo "  h: Human (GRCh38)"
    echo "  y: Yeast (Saccharomyces cerevisiae)"
    echo "  a: Arabidopsis thaliana (TAIR10)"
    echo "  m: Maize (Zea mays)"
    echo "  b: Brassica rapa"
    exit 1
fi

READS_N=30000
MUT_RATES=(0.001 0.01)
THREADS=4

OUTPUT_CSV="benchmark_results_minimap_comparision.${MODE}.csv"

check_tool() {
    if ! command -v $1 &> /dev/null; then
        echo "Error: $1 not found. Please install it."
        exit 1
    fi
}

check_tool dwgsim
check_tool minimap2
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
    cat <(pigz -p 8 -cd "$REF") > "$REF_DECOMP"
fi

if [ ! -f "${INDEX_PREFIX}.idx" ]; then
    echo "Indexing Sing..."
    if command -v cargo &> /dev/null; then
        cargo build --release
    fi
    ./target/release/sing build "$REF_DECOMP" "${INDEX_PREFIX}.idx"
fi

cat << 'EOF' > analyze_benchmark_minimap.py
import sys
import re
import math

TOLERANCE = 10

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
                    mapq = int(parts[4])

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
                    
                    if (flag & 4) or mapq == 0:
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
        'Sing':     sys.argv[3]
    }
    mode = sys.argv[4]

    file_paths = {
        'Minimap2': f'minimap.{mode}.sam',
        'Sing':     f'sing.{mode}.sam',
    }

    truth_dict, tool_results = load_sam_and_build_truth(file_paths)
    total_reads = len(truth_dict)

    print(f"\n--- Summary: {exp_name} (Mode: {mode}) ---")
    print("Tool        | Time_ms  | Precision | Recall   | F1        | TP       | FP       | FN")

    tools_ordered = ['Minimap2', 'Sing']

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

for MUT_RATE in "${MUT_RATES[@]}"; do
        EXP_ID="N_${READS_N}_Mut_${MUT_RATE}.${MODE}"
        echo "================================================="
        echo "Running: $EXP_ID (N: $READS_N, Mut: $MUT_RATE)"
        echo "================================================="

        R1="sim_${EXP_ID}.bwa.read1.fastq.gz"
        R2="sim_${EXP_ID}.bwa.read2.fastq.gz"

        if [ -f "$R1" ] && [ -f "$R2" ]; then
            echo "Simulated reads exist. Skipping dwgsim."
        else
            dwgsim -N "$READS_N" -1 150 -2 150 -R 0 -X 0 -r "$MUT_RATE" -y 0 -H "$REF_DECOMP" "sim_${EXP_ID}"
        fi


        echo "1. Running Sing..."
        START=$(date +%s%N)
        if ./target/release/sing map -t "$THREADS" "${INDEX_PREFIX}.idx" -1 "$R1" -2 "$R2" > "sing.${MODE}.sam"; then
            END=$(date +%s%N)
            TIME_SING=$(( (END - START) / 1000000 ))
        else
            TIME_SING="N/A"
            echo "Sing Failed"
        fi

        echo "2. Running Minimap2..."
        START=$(date +%s%N)
        if minimap2 -t "$THREADS" -ax sr "$REF_DECOMP" "$R1" "$R2" > "minimap.${MODE}.sam" 2>/dev/null; then
            END=$(date +%s%N)
            TIME_MM=$(( (END - START) / 1000000 ))
        else
            TIME_MM="N/A"
            echo "Minimap2 Failed"
        fi

        python3 analyze_benchmark_minimap.py "$EXP_ID" "$TIME_MM" "$TIME_SING" "$MODE" | tee -a "$OUTPUT_CSV"

        rm "minimap.${MODE}.sam"
        echo ""
done

echo "Benchmark Complete. Results saved in $OUTPUT_CSV"
