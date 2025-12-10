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
  -h, --help                Show this help message

Outputs (using the provided prefix):
  <prefix>.trimmed_R1.fastq.gz / <prefix>.trimmed_R2.fastq.gz
  <prefix>.sorted.markdup.bam (+ .bai) with duplicates flagged (not removed)
  <prefix>.alignment.flagstat.txt from samtools flagstat for phylogenetic QC
  <prefix>.vcf.gz (+ .tbi) with allele depths retained
  <prefix>.coverage.tsv from samtools coverage
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

# Normalize the prefix to avoid hidden filenames when a trailing slash is supplied
PREFIX="${PREFIX%/}"
if [[ -z "$PREFIX" ]]; then
  echo "Error: output prefix cannot be empty after normalization." >&2
  exit 1
fi

# Split the prefix into directory and basename without invoking dirname (dash-safe),
# then resolve to an absolute prefix so leading dashes in the basename cannot be
# misinterpreted as options by downstream tools. Also ensure the parent directory
# exists and is writable.
if [[ "$PREFIX" == */* ]]; then
  OUT_DIR_RAW="${PREFIX%/*}"
  PREFIX_BASENAME="${PREFIX##*/}"
else
  OUT_DIR_RAW="."
  PREFIX_BASENAME="$PREFIX"
fi

# Expand a leading tilde for users who supply home-relative paths.
if [[ "$OUT_DIR_RAW" == ~* ]]; then
  OUT_DIR_RAW="${OUT_DIR_RAW/#\~/$HOME}"
fi

if [[ -z "$PREFIX_BASENAME" ]]; then
  echo "Error: output prefix basename cannot be empty." >&2
  exit 1
fi

if ! mkdir -p -- "$OUT_DIR_RAW"; then
  echo "Error: could not create output directory: $OUT_DIR_RAW" >&2
  exit 1
fi

if ! OUT_DIR="$(cd -- "$OUT_DIR_RAW" && pwd)"; then
  echo "Error: could not resolve output directory: $OUT_DIR_RAW" >&2
  exit 1
fi

PREFIX="${OUT_DIR}/${PREFIX_BASENAME}"

if ! touch "${OUT_DIR}/.write_test" 2>/dev/null; then
  echo "Error: output directory is not writable: $OUT_DIR" >&2
  exit 1
fi
rm -f "${OUT_DIR}/.write_test"

# Sanity checks for required tools and inputs to avoid cryptic downstream errors
for tool in cutadapt bwa samtools bcftools; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Error: '$tool' is not in PATH. Please install it or activate the appropriate environment." >&2
    exit 1
  fi
done

for f in "$REF" "$READ1" "$READ2"; do
  if [[ ! -f "$f" ]]; then
    echo "Error: input file not found: $f" >&2
    exit 1
  fi
done

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
  "$REF" "${PREFIX}.trimmed_R1.fastq.gz" "${PREFIX}.trimmed_R2.fastq.gz" \
  | samtools view -b -o "${PREFIX}.unsorted.bam" -

samtools sort -n -@ "$THREADS" -o "${PREFIX}.namesort.bam" "${PREFIX}.unsorted.bam"
samtools fixmate -m "${PREFIX}.namesort.bam" "${PREFIX}.fixmate.bam"
samtools sort -@ "$THREADS" -o "${PREFIX}.positionsort.bam" "${PREFIX}.fixmate.bam"
samtools markdup -@ "$THREADS" -s "${PREFIX}.positionsort.bam" "${PREFIX}.sorted.markdup.bam"
samtools index "${PREFIX}.sorted.markdup.bam"
samtools flagstat "${PREFIX}.sorted.markdup.bam" > "${PREFIX}.alignment.flagstat.txt"

rm -f "${PREFIX}.unsorted.bam" "${PREFIX}.namesort.bam" "${PREFIX}.fixmate.bam" "${PREFIX}.positionsort.bam"

# Step 3: variant calling with allele depths retained
bcftools mpileup -Ou -a AD,ADF,ADR,DP -f "$REF" "${PREFIX}.sorted.markdup.bam" \
  | bcftools call -mv --ploidy 1 --keep-alts --multiallelic-caller -Oz -o "${PREFIX}.vcf.gz"
bcftools index "${PREFIX}.vcf.gz"

# Coverage and variant summaries to validate alignments and support phylogeny building
samtools coverage "${PREFIX}.sorted.markdup.bam" > "${PREFIX}.coverage.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t[%AD]\n' "${PREFIX}.vcf.gz" > "${PREFIX}.variant_summary.tsv"
bcftools stats "${PREFIX}.vcf.gz" > "${PREFIX}.vcf.stats.txt"

# Step 4: consensus with IUPAC codes so low-frequency alleles remain represented
bcftools consensus --iupac-codes -f "$REF" "${PREFIX}.vcf.gz" > "${PREFIX}.consensus.fasta"

echo "Pipeline complete. Outputs written with prefix ${PREFIX}".
