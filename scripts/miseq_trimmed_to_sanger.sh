#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: miseq_trimmed_to_sanger.sh -r <reference.fa> -1 <trimmed_R1.fastq.gz> -2 <trimmed_R2.fastq.gz> -o <output_prefix> [options]

Align pre-trimmed MiSeq reads to a Sanger reference. The script will build a BWA index for the Sanger FASTA if it is absent,
then produce a sorted BAM, variant calls, and an IUPAC consensus that keeps low-frequency alleles.

Required arguments:
  -r, --reference           Reference FASTA for the Sanger sequence(s)
  -1, --reads1              Adapter/quality-trimmed R1 FASTQ (gzipped)
  -2, --reads2              Adapter/quality-trimmed R2 FASTQ (gzipped)
  -o, --output-prefix       Prefix for outputs (e.g., results/sample1)

Optional arguments:
      --threads             Number of CPU threads (default: 4)
      --mismatch-penalty    BWA-MEM mismatch penalty -B (default: 3)
      --gap-open-penalty    BWA-MEM gap open penalties -O (default: 6)
      --gap-extend-penalty  BWA-MEM gap extension penalties -E (default: 1)
      --clip-penalty        BWA-MEM clipping penalty -L (default: 5)
  -h, --help                Show this help message

Outputs (using the provided prefix):
  <prefix>.sorted.bam (+ .bai)
  <prefix>.vcf.gz (+ .tbi) with allele depths retained
  <prefix>.consensus.fasta (IUPAC-coded consensus so minor alleles are not collapsed)

Dependencies: bwa, samtools, and bcftools must be installed and in PATH.
USAGE
}

# Defaults
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

for tool in bwa samtools bcftools; do
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

# Step 1: ensure a BWA index exists for the Sanger reference
bwa_index_missing() {
  local base="$1"
  [[ ! -f "${base}.bwt" || ! -f "${base}.pac" || ! -f "${base}.ann" || ! -f "${base}.amb" || ! -f "${base}.sa" ]]
}

if bwa_index_missing "$REF"; then
  echo "Building BWA index for ${REF}..." >&2
  bwa index "$REF"
fi

# Step 2: align trimmed reads with permissive penalties to retain rare variants
bwa mem \
  -t "$THREADS" \
  -B "$MISMATCH" \
  -O "${GAP_OPEN},${GAP_OPEN}" \
  -E "${GAP_EXT},${GAP_EXT}" \
  -L "$CLIP" \
  "$REF" "$READ1" "$READ2" \
  | samtools view -b -o "${PREFIX}.unsorted.bam" -

samtools sort -@ "$THREADS" -o "${PREFIX}.sorted.bam" "${PREFIX}.unsorted.bam"
samtools index "${PREFIX}.sorted.bam"
rm -f "${PREFIX}.unsorted.bam"

# Step 3: variant calling with allele depths retained
bcftools mpileup -Ou -a AD,ADF,ADR,DP -f "$REF" "${PREFIX}.sorted.bam" \
  | bcftools call -mv --ploidy 1 --keep-alts --multiallelic-caller -Oz -o "${PREFIX}.vcf.gz"
bcftools index "${PREFIX}.vcf.gz"

# Step 4: consensus with IUPAC codes so low-frequency alleles remain represented
bcftools consensus --iupac-codes -f "$REF" "${PREFIX}.vcf.gz" > "${PREFIX}.consensus.fasta"

echo "Alignment complete. Outputs written with prefix ${PREFIX}".
