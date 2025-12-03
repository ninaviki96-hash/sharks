#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: miseq_to_sanger.sh -r <reference.fa> -1 <reads_R1.fastq.gz> -2 <reads_R2.fastq.gz> -o <output_prefix> [options]

Trim adapters, quality-filter MiSeq reads, align to a Sanger reference, and emit a sorted BAM plus an IUPAC consensus that preserves low-frequency alleles.

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
  <prefix>.sorted.bam (+ .bai)
  <prefix>.vcf.gz (+ .tbi) with allele depths retained
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

mkdir -p "$(dirname "$PREFIX")"

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

samtools sort -@ "$THREADS" -o "${PREFIX}.sorted.bam" "${PREFIX}.unsorted.bam"
samtools index "${PREFIX}.sorted.bam"
rm -f "${PREFIX}.unsorted.bam"

# Step 3: variant calling with allele depths retained
bcftools mpileup -Ou -a AD,ADF,ADR,DP -f "$REF" "${PREFIX}.sorted.bam" \
  | bcftools call -mv --ploidy 1 --keep-alts --multiallelic-caller -Oz -o "${PREFIX}.vcf.gz"
bcftools index "${PREFIX}.vcf.gz"

# Step 4: consensus with IUPAC codes so low-frequency alleles remain represented
bcftools consensus --iupac-codes -f "$REF" "${PREFIX}.vcf.gz" > "${PREFIX}.consensus.fasta"

echo "Pipeline complete. Outputs written with prefix ${PREFIX}".
