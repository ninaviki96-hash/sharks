#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: miseq_to_sanger.sh -r <reference.fa> -1 <reads_R1.fastq.gz> -2 <reads_R2.fastq.gz> -o <output_prefix> [options]

Trim adapters, quality-filter MiSeq reads, align to a Sanger reference, and emit a duplicate-marked BAM with coverage/variant summaries plus an IUPAC consensus that preserves low-frequency alleles.

Required arguments:
  -r, --reference           Reference FASTA for the Sanger sequence(s)
  -1, --reads1              R1 FASTQ (gzipped)
  -2, --reads2              R2 FASTQ (gzipped)
  -o, --output-prefix       Prefix for outputs (e.g., results/sample1)

Optional arguments:
      --adapter             Adapter sequence to trim (default: Illumina TruSeq)
      --min-quality         Phred quality cutoff for 3' trimming (default: 20)
      --min-length          Discard reads shorter than this after trimming (default: 120)
      --threads             Number of CPU threads (default: 4)
      --mismatch-penalty    BWA-MEM mismatch penalty -B (default: 3)
      --gap-open-penalty    BWA-MEM gap open penalties -O (default: 6)
      --gap-extend-penalty  BWA-MEM gap extension penalties -E (default: 1)
      --clip-penalty        BWA-MEM clipping penalty -L (default: 5)
      --index-dir           Directory to write/read reference indexes (default: alongside reference)
  -h, --help                Show this help message

Outputs (using the provided prefix):
  <prefix>.trimmed_R1.fastq.gz / <prefix>.trimmed_R2.fastq.gz
  <prefix>.sorted.markdup.bam (+ .bai) with duplicates flagged (not removed)
  <prefix>.vcf.gz (+ .tbi) with allele depths retained
  <prefix>.coverage.txt from samtools coverage
  <prefix>.variant_summary.tsv with key variant metrics
  <prefix>.vcf.stats.txt from bcftools stats
  <prefix>.consensus.fasta (IUPAC-coded consensus so minor alleles are not collapsed)

Dependencies: cutadapt, bwa, samtools, and bcftools must be installed and in PATH.
USAGE
}

# Defaults
ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
MIN_QUAL=20
MIN_LEN=120
THREADS=4
MISMATCH=3
GAP_OPEN=6
GAP_EXT=1
CLIP=5
INDEX_DIR=""

# Parse arguments
if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--reference)
      REF="$2"; shift 2;;
    -1|--reads1)
      READ1="$2"; shift 2;;
    -2|--reads2)
      READ2="$2"; shift 2;;
    -o|--output-prefix)
      PREFIX="$2"; shift 2;;
    --adapter)
      ADAPTER="$2"; shift 2;;
    --min-quality)
      MIN_QUAL="$2"; shift 2;;
    --min-length)
      MIN_LEN="$2"; shift 2;;
    --threads)
      THREADS="$2"; shift 2;;
    --mismatch-penalty)
      MISMATCH="$2"; shift 2;;
    --gap-open-penalty)
      GAP_OPEN="$2"; shift 2;;
    --gap-extend-penalty)
      GAP_EXT="$2"; shift 2;;
    --clip-penalty)
      CLIP="$2"; shift 2;;
    --index-dir)
      INDEX_DIR="$2"; shift 2;;
    -h|--help)
      usage; exit 0;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1;;
  esac
done

if [[ -z "${REF:-}" || -z "${READ1:-}" || -z "${READ2:-}" || -z "${PREFIX:-}" ]]; then
  echo "Missing required arguments." >&2
  usage
  exit 1
fi

mkdir -p "$(dirname "$PREFIX")"

codex/create-scripts-for-illumina-miseq-data-processing-mph3pn
# Step 0: ensure reference (optionally copied to an index directory) is indexed
REF_BASENAME=$(basename "$REF")
if [[ -n "$INDEX_DIR" ]]; then
  mkdir -p "$INDEX_DIR"
  REF_FOR_ALIGN="$INDEX_DIR/$REF_BASENAME"
  if [[ ! -e "$REF_FOR_ALIGN" || "$REF_FOR_ALIGN" -ot "$REF" ]]; then
    echo "[INFO] Copying reference to index directory: $REF_FOR_ALIGN" >&2
    cp "$REF" "$REF_FOR_ALIGN"
  fi
else
  INDEX_DIR=$(dirname "$REF")
  REF_FOR_ALIGN="$REF"
fi

ensure_bwa_index() {
  local ref_path="$1"
  local missing=false
  for ext in amb ann bwt pac sa; do
    [[ -f "${ref_path}.${ext}" ]] || missing=true
  done
  if [[ "$missing" == true ]]; then
    echo "[INFO] BWA index for ${ref_path} not found. Building..." >&2
    bwa index "$ref_path"
  fi
  for ext in amb ann bwt pac sa; do
    if [[ ! -f "${ref_path}.${ext}" ]]; then
      echo "[ERROR] Failed to create BWA index (${ref_path}.${ext}). Ensure the index directory is writable." >&2
      exit 1
    fi
  done
}

ensure_fasta_index() {
  local ref_path="$1"
  if [[ ! -f "${ref_path}.fai" ]]; then
    echo "[INFO] FASTA index for ${ref_path} not found. Building..." >&2
    samtools faidx "$ref_path"
  fi
  if [[ ! -f "${ref_path}.fai" ]]; then
    echo "[ERROR] Failed to create FASTA index (${ref_path}.fai). Ensure the index directory is writable." >&2
    exit 1
  fi
}

ensure_bwa_index "$REF_FOR_ALIGN"
ensure_fasta_index "$REF_FOR_ALIGN"
=======
# Step 0: ensure reference is indexed for alignment and pileup
missing_index=false
for ext in amb ann bwt pac sa; do
  [[ -f "${REF}.${ext}" ]] || missing_index=true
done

if [[ "$missing_index" == true ]]; then
  echo "[INFO] BWA index for ${REF} not found. Building..." >&2
  bwa index "$REF"
fi

if [[ ! -f "${REF}.fai" ]]; then
  echo "[INFO] FASTA index for ${REF} not found. Building..." >&2
  samtools faidx "$REF"
fi
 main

# Step 1: adapter and quality trimming
cutadapt \
  -j "$THREADS" \
  -a "$ADAPTER" -A "$ADAPTER" \
  -q "$MIN_QUAL" \
  -m "$MIN_LEN" \
  -o "${PREFIX}.trimmed_R1.fastq.gz" \
  -p "${PREFIX}.trimmed_R2.fastq.gz" \
  "$READ1" "$READ2"

# Step 2: align with permissive penalties to retain rare variants
bwa mem \
  -t "$THREADS" \
  -B "$MISMATCH" \
  -O "${GAP_OPEN},${GAP_OPEN}" \
  -E "${GAP_EXT},${GAP_EXT}" \
  -L "$CLIP" \
  "$REF_FOR_ALIGN" "${PREFIX}.trimmed_R1.fastq.gz" "${PREFIX}.trimmed_R2.fastq.gz" \
  | samtools view -b -o "${PREFIX}.unsorted.bam" -

samtools sort -n -@ "$THREADS" -o "${PREFIX}.namesort.bam" "${PREFIX}.unsorted.bam"
samtools fixmate -m "${PREFIX}.namesort.bam" "${PREFIX}.fixmate.bam"
samtools sort -@ "$THREADS" -o "${PREFIX}.positionsort.bam" "${PREFIX}.fixmate.bam"
samtools markdup -@ "$THREADS" -s "${PREFIX}.positionsort.bam" "${PREFIX}.sorted.markdup.bam"
samtools index "${PREFIX}.sorted.markdup.bam"

rm -f "${PREFIX}.unsorted.bam" "${PREFIX}.namesort.bam" "${PREFIX}.fixmate.bam" "${PREFIX}.positionsort.bam"

# Step 3: variant calling with allele depths retained
codex/-miseq-1hofbl
bcftools mpileup -Ou -a AD,ADF,ADR,DP -f "$REF" "${PREFIX}.sorted.markdup.bam" \
=======
 codex/-miseq
bcftools mpileup -Ou -a AD,ADF,ADR,DP -f "$REF" "${PREFIX}.sorted.markdup.bam" \
=======
bcftools mpileup -Ou -a AD,ADF,ADR,DP -f "$REF_FOR_ALIGN" "${PREFIX}.sorted.bam" \
 main
 main
  | bcftools call -mv --ploidy 1 --keep-alts --multiallelic-caller -Oz -o "${PREFIX}.vcf.gz"
bcftools index "${PREFIX}.vcf.gz"

# Coverage and variant summaries to validate alignments and support phylogeny building
samtools coverage "${PREFIX}.sorted.markdup.bam" > "${PREFIX}.coverage.txt"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t[%AD]\n' "${PREFIX}.vcf.gz" > "${PREFIX}.variant_summary.tsv"
bcftools stats "${PREFIX}.vcf.gz" > "${PREFIX}.vcf.stats.txt"

# Step 4: consensus with IUPAC codes so low-frequency alleles remain represented
bcftools consensus --iupac-codes -f "$REF_FOR_ALIGN" "${PREFIX}.vcf.gz" > "${PREFIX}.consensus.fasta"

echo "Pipeline complete. Outputs written with prefix ${PREFIX}".
