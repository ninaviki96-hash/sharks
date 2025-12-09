#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: miseq_to_ont.sh -r <ont_contigs.fa> -1 <reads_R1.fastq.gz> -2 <reads_R2.fastq.gz> -o <output_prefix> [options]

Trim adapters, quality-filter MiSeq reads, align to ONT contigs with BWA-MEM (default) or minimap2 (sr preset), and emit
variant calls with coverage/variant summaries tuned for sensitivity in variable regions.

Required arguments:
  -r, --reference           ONT contig/reference FASTA
  -1, --reads1              R1 FASTQ (gzipped)
  -2, --reads2              R2 FASTQ (gzipped)
  -o, --output-prefix       Prefix for outputs (e.g., results/sample1)

Optional arguments:
      --adapter             Adapter sequence to trim (default: Illumina TruSeq)
      --min-quality         Phred quality cutoff for 3' trimming (default: 20)
      --min-length          Discard reads shorter than this after trimming (default: 100)
      --threads             Number of CPU threads (default: 4)
      --aligner             "bwa" or "minimap2" (default: bwa)
      --mismatch-penalty    BWA-MEM mismatch penalty -B (default: 3)
      --gap-open-penalty    BWA-MEM gap open penalties -O (default: 6)
      --gap-extend-penalty  BWA-MEM gap extension penalties -E (default: 1)
      --clip-penalty        BWA-MEM clipping penalty -L (default: 5)
      --downsample-fraction Subsample alignments with samtools view -s (0-1, skipped if >=1 or unset)
      --downsample-seed     Seed for samtools -s (integer, default: 42)
  -h, --help                Show this help message

Outputs (using the provided prefix):
  <prefix>.trimmed_R1.fastq.gz / <prefix>.trimmed_R2.fastq.gz
  <prefix>.sorted.markdup.bam (+ .bai) with duplicates flagged (not removed)
  <prefix>.vcf.gz (+ .tbi) with allele depths retained
  <prefix>.coverage.txt from samtools coverage
  <prefix>.variant_summary.tsv with key variant metrics
  <prefix>.vcf.stats.txt from bcftools stats
  <prefix>.consensus.fasta (IUPAC-coded consensus so minor alleles are not collapsed)

Dependencies: cutadapt, bwa or minimap2, samtools, and bcftools must be installed and in PATH.
USAGE
}

# Defaults
ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
MIN_QUAL=20
MIN_LEN=100
THREADS=4
ALIGNER="bwa"
MISMATCH=3
GAP_OPEN=6
GAP_EXT=1
CLIP=5
DOWNSAMPLE=""
DOWNSAMPLE_SEED=42

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
    --aligner)
      ALIGNER="$2"; shift 2;;
    --mismatch-penalty)
      MISMATCH="$2"; shift 2;;
    --gap-open-penalty)
      GAP_OPEN="$2"; shift 2;;
    --gap-extend-penalty)
      GAP_EXT="$2"; shift 2;;
    --clip-penalty)
      CLIP="$2"; shift 2;;
    --downsample-fraction)
      DOWNSAMPLE="$2"; shift 2;;
    --downsample-seed)
      DOWNSAMPLE_SEED="$2"; shift 2;;
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

for tool in cutadapt samtools bcftools; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Error: '$tool' is not in PATH. Please install it or activate the appropriate environment." >&2
    exit 1
  fi
done

if [[ "$ALIGNER" == "bwa" ]]; then
  if ! command -v bwa >/dev/null 2>&1; then
    echo "Error: 'bwa' is not in PATH but is required for --aligner bwa." >&2
    exit 1
  fi
elif [[ "$ALIGNER" == "minimap2" ]]; then
  if ! command -v minimap2 >/dev/null 2>&1; then
    echo "Error: 'minimap2' is not in PATH but is required for --aligner minimap2." >&2
    exit 1
  fi
else
  echo "Unsupported aligner: $ALIGNER (choose 'minimap2' or 'bwa')." >&2
  exit 1
fi

for f in "$REF" "$READ1" "$READ2"; do
  if [[ ! -f "$f" ]]; then
    echo "Error: input file not found: $f" >&2
    exit 1
  fi
done

if [[ -n "$DOWNSAMPLE" ]]; then
  if ! [[ "$DOWNSAMPLE" =~ ^0?\.[0-9]+$ ]]; then
    echo "--downsample-fraction must be between 0 and 1 (exclusive)." >&2
    exit 1
  fi
  python3 - <<PY
import sys
value = float("$DOWNSAMPLE")
if not (0 < value < 1):
    sys.stderr.write("--downsample-fraction must be >0 and <1.\n")
    sys.exit(1)
PY
  # Strip leading zero for samtools -s format (seed.fraction)
  DOWNSAMPLE_FMT="${DOWNSAMPLE_SEED}.${DOWNSAMPLE#0.}"
fi

OUT_DIR="$(dirname "$PREFIX")"
if ! mkdir -p "$OUT_DIR"; then
  echo "Error: could not create output directory: $OUT_DIR" >&2
  exit 1
fi
if ! touch "$OUT_DIR/.write_test" 2>/dev/null; then
  echo "Error: output directory is not writable: $OUT_DIR" >&2
  exit 1
fi
rm -f "$OUT_DIR/.write_test"

# Step 1: adapter and quality trimming
cutadapt \
  -j "$THREADS" \
  -a "$ADAPTER" -A "$ADAPTER" \
  -q "$MIN_QUAL" \
  -m "$MIN_LEN" \
  -o "${PREFIX}.trimmed_R1.fastq.gz" \
  -p "${PREFIX}.trimmed_R2.fastq.gz" \
  "$READ1" "$READ2"

# Step 2: align with high sensitivity to retain polymorphisms
if [[ "$ALIGNER" == "minimap2" ]]; then
  minimap2 -ax sr -t "$THREADS" "$REF" "${PREFIX}.trimmed_R1.fastq.gz" "${PREFIX}.trimmed_R2.fastq.gz" \
    | samtools view -b -o "${PREFIX}.unsorted.bam" -
elif [[ "$ALIGNER" == "bwa" ]]; then
  bwa mem \
    -t "$THREADS" \
    -B "$MISMATCH" \
    -O "${GAP_OPEN},${GAP_OPEN}" \
    -E "${GAP_EXT},${GAP_EXT}" \
    -L "$CLIP" \
    "$REF" "${PREFIX}.trimmed_R1.fastq.gz" "${PREFIX}.trimmed_R2.fastq.gz" \
    | samtools view -b -o "${PREFIX}.unsorted.bam" -
else
  echo "Unsupported aligner: $ALIGNER (choose 'minimap2' or 'bwa')." >&2
  exit 1
fi

ALIGN_BAM="${PREFIX}.unsorted.bam"

# Optional downsampling to control coverage while avoiding allele bias
if [[ -n "$DOWNSAMPLE" ]]; then
  samtools view -@ "$THREADS" -b -s "$DOWNSAMPLE_FMT" "$ALIGN_BAM" -o "${PREFIX}.downsampled.bam"
  ALIGN_BAM="${PREFIX}.downsampled.bam"
fi

# Step 3: mark duplicates without removal to minimize depth bias
samtools sort -n -@ "$THREADS" -o "${PREFIX}.namesort.bam" "$ALIGN_BAM"
samtools fixmate -m "${PREFIX}.namesort.bam" "${PREFIX}.fixmate.bam"
samtools sort -@ "$THREADS" -o "${PREFIX}.positionsort.bam" "${PREFIX}.fixmate.bam"
samtools markdup -@ "$THREADS" -s "${PREFIX}.positionsort.bam" "${PREFIX}.sorted.markdup.bam"
samtools index "${PREFIX}.sorted.markdup.bam"

rm -f "${PREFIX}.unsorted.bam" "${PREFIX}.downsampled.bam" "${PREFIX}.namesort.bam" "${PREFIX}.fixmate.bam" "${PREFIX}.positionsort.bam"

# Step 4: variant calling with allele depths retained for minor variants
bcftools mpileup -Ou -a AD,ADF,ADR,DP -f "$REF" "${PREFIX}.sorted.markdup.bam" \
  | bcftools call -mv --ploidy 1 --keep-alts --multiallelic-caller -Oz -o "${PREFIX}.vcf.gz"
bcftools index "${PREFIX}.vcf.gz"

# Coverage and variant summaries to validate ONT scaffolds
samtools coverage "${PREFIX}.sorted.markdup.bam" > "${PREFIX}.coverage.txt"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t[%AD]\n' "${PREFIX}.vcf.gz" > "${PREFIX}.variant_summary.tsv"
bcftools stats "${PREFIX}.vcf.gz" > "${PREFIX}.vcf.stats.txt"

# IUPAC-coded consensus so low-frequency alleles remain represented
bcftools consensus --iupac-codes -f "$REF" "${PREFIX}.vcf.gz" > "${PREFIX}.consensus.fasta"

echo "Pipeline complete. Outputs written with prefix ${PREFIX}".
