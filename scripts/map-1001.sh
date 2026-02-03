#PBS -N map-1001
#PBS -l nodes=node03:ppn=128
#PBS -j oe
#PBS -o /home/hyunsu.lim/log/map-1001.log

if [ -n "$PBS_O_WORKDIR" ]; then cd "$PBS_O_WORKDIR"; fi
set -e

WORKDIR=$(pwd)
SING_BIN="${WORKDIR}/sing"
INPUT_DIR="${WORKDIR}/1001_raw" 
REF_DIR="${WORKDIR}/1001/pangenome"
OUT_BASE_DIR="${WORKDIR}/1001G-stats"
SUMMARY_CSV="${WORKDIR}/mapping_rate.csv"
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
find -L "${REF_DIR}/" -name "*.fasta.gz" -o -name "*.fna.gz" | sort | while read -r REF_FASTA; do
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

    local log_file="${OUT_BASE_DIR}/${ref_id}/${SAMPLE_NAME}.log"
    "$SING_BIN" map -t "$THREADS_SING" "$index_path" -1 "$R1_RAM" -2 "$R2_RAM" > /dev/null 2> "$log_file"
    
    local end_time=$(date +%s.%N)
    local duration=$(awk "BEGIN {print $end_time - $start_time}")

    local total=$(grep "Total reads:" "$log_file" | awk '{print $3}')
    local mapped=$(grep "Mapped reads:" "$log_file" | awk '{print $3}')
    local pct=$(grep "Mapped reads:" "$log_file" | awk -F '[()]' '{print $2}' | sed 's/%//')

    if [ -z "$total" ]; then total=0; fi
    if [ -z "$mapped" ]; then mapped=0; fi
    if [ -z "$pct" ]; then pct=0; fi
    
    rm -f "$log_file"
    
    flock -x "${SUMMARY_CSV}.lock" -c "echo '$ref_id,$SAMPLE_NAME,$total,$mapped,$pct,$duration' >> '$SUMMARY_CSV'"
    touch "$done_marker"
}
export -f run_one_ref
export SING_BIN OUT_BASE_DIR SUMMARY_CSV THREADS_SING

echo "Starting Sample-Centric Mapping (From RAW)..."

find -L "${INPUT_DIR}/" -maxdepth 1 -name "*_1.fastq" | sort | while read -r R1_RAW; do
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
