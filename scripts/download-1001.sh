#!/bin/bash
set -e

PROJECT_ACC="SRP056687"
OUT_DIR="./1001"
PARALLEL_JOBS=8

mkdir -p "$OUT_DIR"

if ! command -v fastq-dl &> /dev/null; then
    exit 1
fi

echo "[1/2] Fetching Run List for $PROJECT_ACC..."

LIST_FILE="accession_list.txt"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PROJECT_ACC}&result=read_run&fields=run_accession&format=tsv" \
    | cut -f1 \
    | tail -n +2 \
    > "$LIST_FILE"

TOTAL_SAMPLES=$(wc -l < "$LIST_FILE")
echo "total ${TOTAL_SAMPLES} samples"
echo "list file: $LIST_FILE"

echo "[2/2] Starting Parallel Download ($PARALLEL_JOBS jobs)..."
echo "Target Directory: $OUT_DIR"

cat "$LIST_FILE" | xargs -P "$PARALLEL_JOBS" -I {} \
    fastq-dl -a {} \
    --provider ena \
    --outdir "$OUT_DIR" \
    --silent

echo "Download done."
