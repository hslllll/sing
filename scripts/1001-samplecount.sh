#!/usr/bin/env bash
set -euo pipefail


THREADS=16
PAIRS_DIR="./1001"
REF_FASTA="arabidopsis.fasta"
SING_INDEX="arabidopsis.sing"
SING_BIN="./sing"
BWA_BIN="bwa-mem2"
MINIMAP_BIN="minimap2"
OUT_CSV="1001G-map.csv"
COUNTS=(1 2 4 8 16 32 64 128 256 512)

usage() {
    cat <<EOF
Usage: ./1001G-map.sh [options]
  -d, --dir DIR         FASTQ directory (default: ./1001)
  -r, --ref FASTA       Reference FASTA for bwa-mem2/minimap2 (default: arabidopsis.fasta)
  -i, --idx IDX         sing index path (default: arabidopsis.sing)
  -s, --sing-bin PATH   sing binary (default: ./sing)
  -b, --bwa-bin PATH    bwa-mem2 binary (default: bwa-mem2)
  -m, --minimap-bin PATH minimap2 binary (default: minimap2)
  -t, --threads N       Threads for all tools (default: 16)
  -o, --out FILE        Output CSV (default: 1001G-map.csv)
  -c, --counts LIST     Comma-separated counts (default: 1,2,4,8,16,32,64,128,256)
  -h, --help            Show this help

Notes:
- Uses /usr/bin/time with format: Tool,Count,Real_s,User_s,Sys_s,MaxRSS_KB.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dir) PAIRS_DIR="$2"; shift 2;;
        -r|--ref) REF_FASTA="$2"; shift 2;;
        -i|--idx) SING_INDEX="$2"; shift 2;;
        -s|--sing-bin) SING_BIN="$2"; shift 2;;
        -b|--bwa-bin) BWA_BIN="$2"; shift 2;;
        -m|--minimap-bin) MINIMAP_BIN="$2"; shift 2;;
        -t|--threads) THREADS="$2"; shift 2;;
        -o|--out) OUT_CSV="$2"; shift 2;;
        -c|--counts) IFS=',' read -r -a COUNTS <<< "$2"; shift 2;;
        -h|--help) usage; exit 0;;
        --) shift; break;;
        -*) echo "Unknown option: $1" >&2; usage; exit 1;;
        *) break;;
    esac
done

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Error: $1 not found in PATH" >&2; exit 1; }; }
need_cmd "/usr/bin/time"
need_cmd "$BWA_BIN"
need_cmd "$MINIMAP_BIN"
need_cmd "$SING_BIN"

if [[ ! -d "$PAIRS_DIR" ]]; then
    echo "Error: directory $PAIRS_DIR not found" >&2
    exit 1
fi
if [[ ! -f "$REF_FASTA" ]]; then
    echo "Error: reference FASTA $REF_FASTA not found" >&2
    exit 1
fi
if [[ ! -f "$SING_INDEX" ]]; then
    echo "Error: sing index $SING_INDEX not found" >&2
    exit 1
fi

declare -a R1_LIST
declare -a R2_LIST
while IFS= read -r -d '' r1; do
    r2="${r1/_1.fastq.gz/_2.fastq.gz}"
    if [[ -f "$r2" ]]; then
        R1_LIST+=("$r1")
        R2_LIST+=("$r2")
    fi
done < <(find "$PAIRS_DIR" -maxdepth 1 -type f -name 'SRR*_1.fastq.gz' -print0 | sort -z)

TOTAL=${#R1_LIST[@]}
if [[ $TOTAL -eq 0 ]]; then
    echo "Error: no SRR*_1.fastq.gz pairs found under $PAIRS_DIR" >&2
    exit 1
fi

echo "Found $TOTAL paired samples in $PAIRS_DIR" >&2

tmp_log=$(mktemp)
trap 'rm -f "$tmp_log"' EXIT

echo "Tool,Count,Real_s,User_s,Sys_s,MaxRSS_KB" > "$OUT_CSV"

run_time_loop() {
    local tool="$1"; shift
    local count="$1"; shift
    local pairs_file="$1"; shift
    local body="$1"; shift

    local script
    script=$(mktemp)
    cat > "$script" <<EOF
set -euo pipefail
pairs_file="$pairs_file"
while IFS=$'\t' read -r r1 r2; do
$body
done < "\$pairs_file"
EOF

    echo "[$tool] n=$count" >&2
    /usr/bin/time -f "${tool},${count},%e,%U,%S,%M" -o "$tmp_log" bash "$script"
    cat "$tmp_log" >> "$OUT_CSV"
    rm -f "$script"
}

for count in "${COUNTS[@]}"; do
    if [[ $count -gt $TOTAL ]]; then
        echo "Skipping count=$count (only $TOTAL pairs available)" >&2
        continue
    fi

    pairs_file=$(mktemp)
    for ((i=0; i<count; i++)); do
        printf "%s\t%s\n" "${R1_LIST[$i]}" "${R2_LIST[$i]}" >> "$pairs_file"
    done

    sing_body="$SING_BIN map -t $THREADS $SING_INDEX -1 \"\$r1\" -2 \"\$r2\" > /dev/null"
    bwa_body="$BWA_BIN mem -t $THREADS $REF_FASTA \"\$r1\" \"\$r2\" > /dev/null"
    minimap_body="$MINIMAP_BIN -ax sr -t $THREADS $REF_FASTA \"\$r1\" \"\$r2\" > /dev/null"

    run_time_loop "sing" "$count" "$pairs_file" "$sing_body"
    run_time_loop "bwa-mem2" "$count" "$pairs_file" "$bwa_body"
    run_time_loop "minimap2" "$count" "$pairs_file" "$minimap_body"

    rm -f "$pairs_file"

done

echo "Benchmark complete. Results: $OUT_CSV" >&2
