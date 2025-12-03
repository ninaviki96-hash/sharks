#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: ont_to_sanger.sh -r <reference.fa> -i <reads.fastq[.gz]> -o <output_prefix> [options]

Align ONT reads to Sanger reference sequences with minimap2 (-ax map-ont) using scoring tuned to reduce over-clipping,
produce a sorted BAM plus an IUPAC consensus that preserves low-frequency variants, and optionally polish with racon
and/or medaka.

Required arguments:
  -r, --reference           Reference FASTA for the Sanger sequence(s)
  -i, --reads               ONT reads FASTQ (optionally gzipped)
  -o, --output-prefix       Prefix for outputs (e.g., results/sample1)

Optional arguments:
      --threads             Number of CPU threads (default: 4)
      --min-length          Minimum read length to keep with filtlong (default: 800)
      --keep-percent        Percentage of best reads to retain with filtlong (default: 90)
      --skip-filtlong       Skip read quality filtering (use raw reads)
      --mm2-mismatch        minimap2 -B mismatch penalty to discourage clipping (default: 4)
      --mm2-gap-open        minimap2 -O gap open penalties short,long (default: 6,24)
      --mm2-gap-extend      minimap2 -E gap extension penalties short,long (default: 2,1)
      --mm2-end-bonus       minimap2 --end-bonus to reward full-length alignments (default: 10)
      --mm2-zdrop           minimap2 -z soft,hard z-drop to avoid premature clipping (default: 200,400)
      --racon-rounds        Number of racon polishing rounds (default: 0)
      --medaka-model        medaka model name to run medaka_consensus after racon (e.g., r941_min_high_g360)
  -h, --help                Show this help message

Outputs (using the provided prefix):
  <prefix>.filtered.fastq.gz               (or the input reads if --skip-filtlong is used)
  <prefix>.sorted.bam (+ .bai)             Alignments to the Sanger reference with tuned clipping penalties
  <prefix>.vcf.gz (+ .tbi)                 Haploid VCF with allele depths retained to preserve minor evidence
  <prefix>.consensus.fasta                 IUPAC-coded consensus from the VCF so low-frequency alleles remain represented
  <prefix>.racon_roundN.fasta              (optional) Haploid racon consensus after each polishing iteration
  <prefix>.medaka/consensus.fasta          (optional) medaka consensus using racon/bcftools draft

Dependencies: minimap2, filtlong, samtools, and bcftools must be installed. racon and medaka are optional for polishing.
USAGE
}

THREADS=4
MIN_LEN=800
KEEP_PCT=90
USE_FILTLONG=1
MM2_MISMATCH=4
MM2_GAP_OPEN="6,24"
MM2_GAP_EXT="2,1"
MM2_END_BONUS=10
MM2_ZDROP="200,400"
RACON_ROUNDS=0
MEDAKA_MODEL=""

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--reference)
      REF="$2"; shift 2;;
    -i|--reads)
      READS="$2"; shift 2;;
    -o|--output-prefix)
      PREFIX="$2"; shift 2;;
    --threads)
      THREADS="$2"; shift 2;;
    --min-length)
      MIN_LEN="$2"; shift 2;;
    --keep-percent)
      KEEP_PCT="$2"; shift 2;;
    --skip-filtlong)
      USE_FILTLONG=0; shift 1;;
    --mm2-mismatch)
      MM2_MISMATCH="$2"; shift 2;;
    --mm2-gap-open)
      MM2_GAP_OPEN="$2"; shift 2;;
    --mm2-gap-extend)
      MM2_GAP_EXT="$2"; shift 2;;
    --mm2-end-bonus)
      MM2_END_BONUS="$2"; shift 2;;
    --mm2-zdrop)
      MM2_ZDROP="$2"; shift 2;;
    --racon-rounds)
      RACON_ROUNDS="$2"; shift 2;;
    --medaka-model)
      MEDAKA_MODEL="$2"; shift 2;;
    -h|--help)
      usage; exit 0;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1;;
  esac
done

if [[ -z "${REF:-}" || -z "${READS:-}" || -z "${PREFIX:-}" ]]; then
  echo "Missing required arguments." >&2
  usage
  exit 1
fi

mkdir -p "$(dirname "$PREFIX")"

FILTERED_READS="$READS"
if [[ "$USE_FILTLONG" -eq 1 ]]; then
  if ! command -v filtlong >/dev/null 2>&1; then
    echo "filtlong is required for filtering unless --skip-filtlong is provided." >&2
    exit 1
  fi
  FILTERED_READS="${PREFIX}.filtered.fastq.gz"
  filtlong --min_length "$MIN_LEN" --keep_percent "$KEEP_PCT" "$READS" | gzip > "$FILTERED_READS"
fi

# Align ONT reads with minimap2 tuned to reduce over-clipping at variable ends
ALIGN_BAM="${PREFIX}.unsorted.bam"
minimap2 -ax map-ont \
  -B "$MM2_MISMATCH" \
  -O "$MM2_GAP_OPEN" \
  -E "$MM2_GAP_EXT" \
  --end-bonus "$MM2_END_BONUS" \
  -z "$MM2_ZDROP" \
  -t "$THREADS" \
  "$REF" "$FILTERED_READS" \
  | samtools view -b -o "$ALIGN_BAM" -

samtools sort -@ "$THREADS" -o "${PREFIX}.sorted.bam" "$ALIGN_BAM"
samtools index "${PREFIX}.sorted.bam"
rm -f "$ALIGN_BAM"

# Variant calling retaining allele depths so minor alleles remain visible
bcftools mpileup -Ou -a AD,ADF,ADR,DP -f "$REF" "${PREFIX}.sorted.bam" \
  | bcftools call -mv --ploidy 1 --keep-alts --multiallelic-caller -Oz -o "${PREFIX}.vcf.gz"
bcftools index "${PREFIX}.vcf.gz"

bcftools consensus --iupac-codes -f "$REF" "${PREFIX}.vcf.gz" > "${PREFIX}.consensus.fasta"

# Optional racon polishing iterations starting from the bcftools consensus
CURRENT_DRAFT="${PREFIX}.consensus.fasta"
if [[ "$RACON_ROUNDS" -gt 0 ]]; then
  if ! command -v racon >/dev/null 2>&1; then
    echo "racon not found in PATH but --racon-rounds was requested." >&2
    exit 1
  fi
  for i in $(seq 1 "$RACON_ROUNDS"); do
    minimap2 -ax map-ont -t "$THREADS" "$CURRENT_DRAFT" "$FILTERED_READS" > "${PREFIX}.racon_round${i}.sam"
    racon -t "$THREADS" "$FILTERED_READS" "${PREFIX}.racon_round${i}.sam" "$CURRENT_DRAFT" > "${PREFIX}.racon_round${i}.fasta"
    CURRENT_DRAFT="${PREFIX}.racon_round${i}.fasta"
    rm -f "${PREFIX}.racon_round${i}.sam"
  done
fi

# Optional medaka polishing using the latest draft
if [[ -n "$MEDAKA_MODEL" ]]; then
  if ! command -v medaka_consensus >/dev/null 2>&1; then
    echo "medaka_consensus not found in PATH but --medaka-model was provided." >&2
    exit 1
  fi
  medaka_consensus \
    -i "$FILTERED_READS" \
    -d "$CURRENT_DRAFT" \
    -o "${PREFIX}.medaka" \
    -t "$THREADS" \
    -m "$MEDAKA_MODEL"
fi

echo "Pipeline complete. Outputs written with prefix ${PREFIX}".
