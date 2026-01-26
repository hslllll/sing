#!/usr/bin/env bash
set -e
eval "$(micromamba shell hook --shell=bash)"
micromamba activate tools

if [ -z "$1" ]; then
    echo "Usage: ./gatk-sim.sh [y|a|m|h]"
    echo "  y : Yeast"
    echo "  a : Arabidopsis"
    echo "  m : Maize"
    echo "  h : Human"
    exit 1
fi

MODE=$1
THREADS=16

COVERAGE=10
READ_LEN=150
MUT_RATE=0.001
INDEL_FRAC=0.05
ERR_RATE=0.01

OUT_DIR="${MODE}.gatk"

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi

case $MODE in
    y) SPECIES="yeast"; GENOME_SIZE=12157105
       REF_NAME="yeast.fa"
       REF_URL="http://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz" ;;
    a) SPECIES="arabidopsis"; GENOME_SIZE=135000000
       REF_NAME="arabidopsis.fa"
       REF_URL="http://ftp.ensemblgenomes.org/pub/plants/release-57/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz" ;;
    m) SPECIES="maize"; GENOME_SIZE=250000000
       REF_NAME="Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa"
       REF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz" ;;
    h) SPECIES="human"; GENOME_SIZE=3.2e9
       REF_NAME="GRCh38_latest_genomic.fna"
       REF_URL="https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz" ;;
    *) echo "Error: Unknown mode"; exit 1 ;;
esac

REF="$OUT_DIR/$REF_NAME"
SIM_PREFIX="$OUT_DIR/sim_reads"

if [ ! -f "$REF" ]; then
    wget -O "${REF}.gz" "$REF_URL"
    zcat "${REF}.gz" | sed 's/ .*//' > "$REF"
    rm "${REF}.gz"
fi

R1="${SIM_PREFIX}.bwa.read1.fastq.gz"
R2="${SIM_PREFIX}.bwa.read2.fastq.gz"
RAW_TRUTH="${SIM_PREFIX}.mutations.vcf"

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    dwgsim -C $COVERAGE -1 $READ_LEN -2 $READ_LEN \
           -r $MUT_RATE -R $INDEL_FRAC -e $ERR_RATE -E $ERR_RATE -y 0 \
           "$REF" "$SIM_PREFIX" > /dev/null
fi

TRUTH_VCF="$OUT_DIR/truth.vcf.gz"
if [ ! -f "$TRUTH_VCF" ]; then
    bgzip -c "$RAW_TRUTH" > "$TRUTH_VCF"
    tabix -p vcf "$TRUTH_VCF"
fi

TOOLS=("sing" "bwa" "mini")
declare -A SECONDS_MAP
declare -A MEM_MAP

run_timed_cmd() {
    TOOL=$1; CMD=$2; LOG="$OUT_DIR/${TOOL}.time_log"
    ./time -f "%e %M" -o "$LOG" bash -c "$CMD"
}

if [ ! -f "$OUT_DIR/index.sing" ]; then
    cargo run --release -- build "$REF" "$OUT_DIR/index.sing" > /dev/null 2>&1
fi

if [ ! -f "$OUT_DIR/sing.bam" ]; then
    run_timed_cmd "sing" "./target/release/sing map -t $THREADS $OUT_DIR/index.sing -1 $R1 -2 $R2 > $OUT_DIR/sing.sam"
    samtools sort -@ $THREADS -O BAM -o "$OUT_DIR/sing.bam" "$OUT_DIR/sing.sam"
    samtools index -@ $THREADS "$OUT_DIR/sing.bam"
fi
read REAL_SEC MAX_KB < "$OUT_DIR/sing.time_log"
SECONDS_MAP["sing"]=$REAL_SEC
MEM_MAP["sing"]=$MAX_KB

if [ ! -f "${REF}.bwt.2bit.64" ]; then 
    bwa-mem2 index -p "$REF" "$REF" > /dev/null 2>&1
fi

if [ ! -f "$OUT_DIR/bwa.bam" ]; then
    run_timed_cmd "bwa" "bwa-mem2 mem -t $THREADS -R \"@RG\tID:bwa\tSM:${SPECIES}\tPL:ILLUMINA\" $REF $R1 $R2 > $OUT_DIR/bwa.sam"
    samtools sort -@ $THREADS -O BAM -o "$OUT_DIR/bwa.bam" "$OUT_DIR/bwa.sam"
    samtools index -@ $THREADS "$OUT_DIR/bwa.bam"
fi
read REAL_SEC MAX_KB < "$OUT_DIR/bwa.time_log"
SECONDS_MAP["bwa"]=$REAL_SEC
MEM_MAP["bwa"]=$MAX_KB

if [ ! -f "$OUT_DIR/mini.bam" ]; then
    run_timed_cmd "mini" "minimap2 -t $THREADS -ax sr -R \"@RG\tID:mini\tSM:${SPECIES}\tPL:ILLUMINA\" $REF $R1 $R2 > $OUT_DIR/mini.sam"
    samtools sort -@ $THREADS -O BAM -o "$OUT_DIR/mini.bam" "$OUT_DIR/mini.sam"
    samtools index -@ $THREADS "$OUT_DIR/mini.bam"
fi
read REAL_SEC MAX_KB < "$OUT_DIR/mini.time_log"
SECONDS_MAP["mini"]=$REAL_SEC
MEM_MAP["mini"]=$MAX_KB

DICT="${REF%.fa}.dict"
if [ "${REF}" == "${REF%.fa}" ]; then DICT="${REF}.dict"; fi

if [ ! -f "$DICT" ]; then
    gatk CreateSequenceDictionary -R "$REF" -O "$DICT" > /dev/null 2>&1
fi
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "$REF"
fi

for tool in "${TOOLS[@]}"; do
    if [ ! -f "$OUT_DIR/${tool}.g.vcf.gz" ]; then
        gatk --java-options "-Xmx16g" HaplotypeCaller \
            -R "$REF" \
            -I "$OUT_DIR/${tool}.bam" \
            -O "$OUT_DIR/${tool}.g.vcf.gz" \
            -ERC GVCF > /dev/null 2>&1
    fi

    if [ ! -f "$OUT_DIR/${tool}.final.vcf.gz" ]; then
        gatk --java-options "-Xmx16g" GenotypeGVCFs \
            -R "$REF" \
            -V "$OUT_DIR/${tool}.g.vcf.gz" \
            -O "$OUT_DIR/${tool}.final.vcf.gz" > /dev/null 2>&1
    fi
done

TRUTH_NORM="$OUT_DIR/truth.norm.vcf.gz"
if [ ! -f "$TRUTH_NORM" ]; then
    bcftools norm -f "$REF" -m -any -O z -o "$TRUTH_NORM" "$TRUTH_VCF"
    tabix -p vcf "$TRUTH_NORM"
fi

declare -A TP_MAP FP_MAP FN_MAP PREC_MAP REC_MAP F1_MAP

for tool in "${TOOLS[@]}"; do
    EVAL_DIR="$OUT_DIR/eval_${tool}"
    mkdir -p "$EVAL_DIR"
    
    if [ ! -f "$EVAL_DIR/tool.norm.vcf.gz" ]; then
        bcftools norm -f "$REF" -m -any -O z -o "$EVAL_DIR/tool.norm.vcf.gz" "$OUT_DIR/${tool}.final.vcf.gz"
        tabix -p vcf "$EVAL_DIR/tool.norm.vcf.gz"
    fi
    
    bcftools isec -p "$EVAL_DIR" -Oz "$TRUTH_NORM" "$EVAL_DIR/tool.norm.vcf.gz" > /dev/null 2>&1
    
    FN=$(zcat "$EVAL_DIR/0000.vcf.gz" | grep -v "^#" | wc -l)
    FP=$(zcat "$EVAL_DIR/0001.vcf.gz" | grep -v "^#" | wc -l)
    TP=$(zcat "$EVAL_DIR/0002.vcf.gz" | grep -v "^#" | wc -l)
    
    TP_MAP[$tool]=$TP
    FP_MAP[$tool]=$FP
    FN_MAP[$tool]=$FN
    
    eval $(awk -v tp=$TP -v fp=$FP -v fn=$FN 'BEGIN {
        if (tp==0) { print "P=0.00 R=0.00 F=0.00" }
        else {
            prec = tp / (tp + fp) * 100;
            rec = tp / (tp + fn) * 100;
            f1 = 2 * tp / (2 * tp + fp + fn) * 100;
            printf "P=%.2f R=%.2f F=%.2f", prec, rec, f1
        }
    }')
    
    PREC_MAP[$tool]=$P
    REC_MAP[$tool]=$R
    F1_MAP[$tool]=$F
done

echo ""
echo "##############################################################################"
echo "            SIMULATION REPORT [$SPECIES] (Clean ${COVERAGE}X, Real - No Filter)"
echo "##############################################################################"

echo ""
echo "[1] Absolute Performance"
echo "------------------------------------------------------------------"
printf "%-10s | %-12s | %-12s | %-15s\n" "Tool" "Time(s)" "Mem(KB)" "Speedup(vs BWA)"
echo "------------------------------------------------------------------"
BWA_TIME=${SECONDS_MAP["bwa"]}
for tool in "${TOOLS[@]}"; do
    TIME=${SECONDS_MAP[$tool]}
    RATIO=$(awk -v bwa="$BWA_TIME" -v tool="$TIME" 'BEGIN {
        if (tool == 0) print "err"
        else printf "%.2f", bwa / tool
    }')
    printf "%-10s | %-12s | %-12s | %-15sx\n" $tool ${SECONDS_MAP[$tool]} ${MEM_MAP[$tool]} $RATIO
done
echo "------------------------------------------------------------------"

echo ""
echo "[2] Accuracy vs Ground Truth (Raw VCF - No Hard Filtering)"
echo "---------------------------------------------------------------------------------------"
printf "%-10s | %-8s | %-8s | %-8s | %-8s | %-8s | %-8s\n" "Tool" "TP" "FP" "FN" "Recall" "Precis" "F1"
echo "---------------------------------------------------------------------------------------"
for tool in "${TOOLS[@]}"; do
    printf "%-10s | %-8s | %-8s | %-8s | %-8s%% | %-8s%% | %-8s%%\n" \
        $tool ${TP_MAP[$tool]} ${FP_MAP[$tool]} ${FN_MAP[$tool]} \
        ${REC_MAP[$tool]} ${PREC_MAP[$tool]} ${F1_MAP[$tool]}
done
echo "---------------------------------------------------------------------------------------"
echo " * Condition: ${COVERAGE}X, 0.1% Mut, 10% Indel, 1% Error"
echo " * Filter: NONE (Raw Output from GenotypeGVCFs)"
echo "##############################################################################"
