#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: ./compare_bam_genotype.sh [-r ref.fa] [-o out_dir] [-t threads] [-m java_mem] bam1.bam bam2.bam bam3.bam

Options:
  -r, --ref        Reference FASTA (default: ref.fa)
  -o, --out        Output directory (default: compare_bam_out)
  -t, --threads    Threads (default: 8)
  -m, --memory     Java heap for GATK (default: 8g)
  -h, --help       Show this help
EOF
}

REF="ref.fa"
OUT_DIR="compare_bam_out"
THREADS=${THREADS:-8}
JAVA_MEM=${JAVA_MEM:-8g}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--ref) REF="$2"; shift 2;;
        -o|--out) OUT_DIR="$2"; shift 2;;
        -t|--threads) THREADS="$2"; shift 2;;
        -m|--memory) JAVA_MEM="$2"; shift 2;;
        -h|--help) usage; exit 0;;
        --) shift; break;;
        -*) echo "Unknown option: $1" >&2; usage; exit 1;;
        *) break;;
    esac
done

if [[ $# -ne 3 ]]; then
    echo "Error: provide exactly three BAM files." >&2
    usage
    exit 1
fi

BAMS=()
for bam in "$@"; do
    if [[ ! -f "$bam" ]]; then
        echo "Error: $bam not found." >&2
        exit 1
    fi
    BAMS+=("$bam")
done

if [[ ! -f "$REF" ]]; then
    echo "Error: reference FASTA $REF not found." >&2
    exit 1
fi

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Error: $1 not in PATH." >&2; exit 1; }; }
for bin in gatk samtools bcftools python3; do
    need_cmd "$bin"
done

if ! python3 - <<'PY' >/dev/null 2>&1; then
import importlib
for mod in ("matplotlib", "matplotlib_venn"):
    importlib.import_module(mod)
PY
    echo "Error: python3 modules matplotlib and matplotlib_venn are required." >&2
    exit 1
fi

mkdir -p "$OUT_DIR"/vcfs "$OUT_DIR"/plots "$OUT_DIR"/work

DICT="${REF%.*}.dict"
if [[ ! -f "$DICT" ]]; then
    gatk CreateSequenceDictionary -R "$REF" -O "$DICT"
fi
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

KEYS=()
LABELS=()

process_sample() {
    local in_bam="$1"; shift
    local label="$1"; shift

    local bam="$OUT_DIR/work/${label}.bam"
    local gvcf="$OUT_DIR/vcfs/${label}.g.vcf.gz"
    local raw_vcf="$OUT_DIR/vcfs/${label}.raw.vcf.gz"
    local vcf="$OUT_DIR/vcfs/${label}.vcf.gz"
    local keys="$OUT_DIR/work/${label}.keys.txt"

    samtools addreplacerg -r "ID:${label}\tSM:${label}\tPL:ILLUMINA" "$in_bam" \
        | samtools sort -@ "$THREADS" -O BAM -o "$bam"
    samtools index -@ "$THREADS" "$bam"

    gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
        -R "$REF" \
        -I "$bam" \
        -O "$gvcf" \
        -ERC GVCF

    gatk --java-options "-Xmx${JAVA_MEM}" GenotypeGVCFs \
        -R "$REF" \
        -V "$gvcf" \
        -O "$raw_vcf"

    bcftools norm -f "$REF" -m -any -O z -o "$vcf" "$raw_vcf"
    tabix -p vcf "$vcf"

    bcftools query -f '%CHROM:%POS:%REF:%ALT\n' "$vcf" > "$keys"

    rm -f "$bam" "${bam}.bai"

    KEYS+=("$keys")
    LABELS+=("$label")
}

for bam in "${BAMS[@]}"; do
    base=$(basename "$bam")
    label=${base%.bam}
    process_sample "$bam" "$label"
done

PLOT_PNG="$OUT_DIR/plots/venn.svg"
COUNTS_TSV="$OUT_DIR/plots/venn_counts.tsv"

python3 - "$PLOT_PNG" "$COUNTS_TSV" "${LABELS[@]}" "${KEYS[@]}" <<'PY'
import sys
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib_venn import venn3

if len(sys.argv) != 9:
    sys.exit("Internal error: expected png, tsv, 3 labels, 3 key files")

png_path, tsv_path, l1, l2, l3, k1, k2, k3 = sys.argv[1:]
labels = [l1, l2, l3]
key_files = [k1, k2, k3]

sets = []
for f in key_files:
    with open(f) as fh:
        sets.append(set(line.strip() for line in fh if line.strip()))

counts = {
    f"only_{labels[0]}": len(sets[0] - sets[1] - sets[2]),
    f"only_{labels[1]}": len(sets[1] - sets[0] - sets[2]),
    f"only_{labels[2]}": len(sets[2] - sets[0] - sets[1]),
    f"{labels[0]}_{labels[1]}": len((sets[0] & sets[1]) - sets[2]),
    f"{labels[0]}_{labels[2]}": len((sets[0] & sets[2]) - sets[1]),
    f"{labels[1]}_{labels[2]}": len((sets[1] & sets[2]) - sets[0]),
    "all_three": len(sets[0] & sets[1] & sets[2]),
}

with open(tsv_path, "w") as out:
    out.write("region\tcount\n")
    for k, v in counts.items():
        out.write(f"{k}\t{v}\n")

plt.figure(figsize=(6, 6))
venn3(subsets=sets, set_labels=labels)
plt.title("Genotyped Variant Overlap")
plt.tight_layout()
Path(png_path).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(png_path)
PY

rm -rf "$OUT_DIR/vcfs" "$OUT_DIR/work"

echo "Done."
echo "  Venn PNG  : $PLOT_PNG"
echo "  Venn TSV  : $COUNTS_TSV"
echo "  (Intermediate VCFs and BAMs have been removed)"