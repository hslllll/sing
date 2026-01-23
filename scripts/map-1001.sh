#!/bin/bash

if [ -n "$PBS_O_WORKDIR" ]; then cd "$PBS_O_WORKDIR"; fi
set -e

WORKDIR=$(pwd)
SING_BIN="${WORKDIR}/sing"
INPUT_DIR="${WORKDIR}/1001_raw" 
REF_DIR="${WORKDIR}/1001/pangenome"
OUT_BASE_DIR="${WORKDIR}/1001G-stats"
SUMMARY_CSV="${WORKDIR}/mapping_rate_matrix.csv"
RAM_DISK="/dev/shm/sing_scratch"
REF_LIST_FILE="${WORKDIR}/ref_list.txt"

REF_PARALLEL=32
THREADS_SING=4

if [ ! -x "$SING_BIN" ]; then exit 1; fi
mkdir -p "$OUT_BASE_DIR"
mkdir -p "$RAM_DISK"

if [ ! -f "$SUMMARY_CSV" ]; then
    echo "Reference,Sample,Total_Reads,Mapped_Reads,Percentage,Time_Sec" > "$SUMMARY_CSV"
fi

rm -f "$REF_LIST_FILE"
find "$REF_DIR" -name "*.fasta.gz" -o -name "*.fna.gz" | sort | while read -r REF_FASTA; do
    bname=$(basename "$REF_FASTA")
    RID="${bname%%.*}"
    INDEX_PATH="${REF_DIR}/${RID}.sing"
    
    if [ ! -f "$INDEX_PATH" ]; then
        "$SING_BIN" build "$REF_FASTA" "$INDEX_PATH"
    fi
    
    echo "$RID $INDEX_PATH" >> "$REF_LIST_FILE"
done

run_one_ref() {
    local ref_id=$1
    local index_path=$2
    
    
    local done_marker="${OUT_BASE_DIR}/${ref_id}/${SAMPLE_NAME}.done"
    
    if [ -z "$index_path" ]; then return; fi
    
    mkdir -p "${OUT_BASE_DIR}/${ref_id}"
    if [ -f "$done_marker" ]; then return; fi
    
    local start_time=$(date +%s.%N)
    
    local result
    result=$("$SING_BIN" map -t "$THREADS_SING" "$index_path" -1 "$R1_RAM" -2 "$R2_RAM" 2>/dev/null | \
    awk 'BEGIN {t=0;m=0} /^@/ {next} {t++; if(and($2,4)==0) m++} END {print t, m}')
    
    read -r total mapped <<< "$result"
    
    local end_time=$(date +%s.%N)
    local duration=$(awk "BEGIN {print $end_time - $start_time}")
    
    local pct=0
    if [ "$total" -gt 0 ]; then
        pct=$(awk "BEGIN {printf \"%.2f\", ($mapped/$total)*100}")
    fi
    
    flock -x "${SUMMARY_CSV}.lock" -c "echo '$ref_id,$SAMPLE_NAME,$total,$mapped,$pct,$duration' >> '$SUMMARY_CSV'"
    touch "$done_marker"
}
export -f run_one_ref
export SING_BIN OUT_BASE_DIR SUMMARY_CSV THREADS_SING

echo "Starting Sample-Centric Mapping (From RAW)..."

find "$INPUT_DIR" -maxdepth 1 -name "*_1.fastq" | sort | while read -r R1_RAW; do
    R2_RAW="${R1_RAW/_1.fastq/_2.fastq}"
    
    if [ ! -f "$R2_RAW" ]; then continue; fi
    
    export SAMPLE_NAME=$(basename "$R1_RAW" | sed 's/_1.fastq//')
    
    export R1_RAM="${RAM_DISK}/${SAMPLE_NAME}_1.fastq"
    export R2_RAM="${RAM_DISK}/${SAMPLE_NAME}_2.fastq"
    
    echo ">>> Processing Sample: $SAMPLE_NAME"
    
    cp "$R1_RAW" "$R1_RAM" &
    PID1=$!
    cp "$R2_RAW" "$R2_RAM" &
    PID2=$!
    wait $PID1 $PID2
    
    cat "$REF_LIST_FILE" | \
    xargs -P "$REF_PARALLEL" -n 2 bash -c 'run_one_ref "$@"' _ 
    
    rm -f "$R1_RAM" "$R2_RAM"
    
done

echo "All Done."
