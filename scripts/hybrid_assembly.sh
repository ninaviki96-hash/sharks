#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: hybrid_assembly.sh -i <ont_reads.fastq.gz> -o <output_prefix> [options]

Run a long-read-first hybrid assembly with Flye, optionally scaffold with hybridSPAdes
when short reads are supplied, polish with short-read/Sanger alignments via Pilon, and
emit validation metrics to confirm variable immune regions are retained.

Required arguments:
  -i, --ont-reads            ONT reads (FASTQ, gzipped)
  -o, --output-prefix        Prefix for outputs (e.g., results/sample1)

Optional arguments:
  -1, --short-reads1         MiSeq/Illumina R1 FASTQ (gzipped)
  -2, --short-reads2         MiSeq/Illumina R2 FASTQ (gzipped)
  -s, --sanger               Sanger reference/reads FASTA or FASTQ
      --genome-size          Genome/target size for Flye (default: 3m)
      --threads              CPU threads (default: 8)
      --memory-gb            Memory limit passed to hybridSPAdes (default: 32)
      --pilon-rounds         Number of Pilon polishing rounds (default: 1)
      --pilon-bin            Pilon executable (default: pilon)
      --skip-hybridspades    Skip hybridSPAdes scaffolding even if short reads are supplied
  -h, --help                 Show this help message

Outputs (using the provided prefix):
  <prefix>.flye/assembly.fasta                  Flye contigs from ONT reads
  <prefix>.hybridspades/scaffolds.fasta         Hybrid scaffolds (if short reads provided)
  <prefix>.pilon_roundX/pilon_roundX.fasta      Pilon-polished assemblies
  <prefix>.alignments.short.bam                 Short-read alignment to final assembly (+ .bai)
  <prefix>.alignments.sanger.bam                Sanger alignment to final assembly (+ .bai)
  <prefix>.variants.vcf.gz (+ .tbi)             Variant calls with allele depths
  <prefix>.consensus.iupac.fasta                IUPAC consensus to retain rare alleles
  <prefix>.coverage.tsv / coverage.png          Depth table and optional plot (via gnuplot)
  <prefix>.qc.summary.txt                       Basic assembly/QC metrics

Dependencies: flye, minimap2, samtools, bcftools, java+pilon, and optionally spades.py (hybridSPAdes) and gnuplot.
USAGE
}

# Defaults
GENOME_SIZE="3m"
THREADS=8
MEMORY_GB=32
PILON_ROUNDS=1
PILON_BIN="pilon"
SKIP_HYBRIDSPADES=false

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--ont-reads)
      ONT_READS="$2"; shift 2;;
    -o|--output-prefix)
      PREFIX="$2"; shift 2;;
    -1|--short-reads1)
      READ1="$2"; shift 2;;
    -2|--short-reads2)
      READ2="$2"; shift 2;;
    -s|--sanger)
      SANGER="$2"; shift 2;;
    --genome-size)
      GENOME_SIZE="$2"; shift 2;;
    --threads)
      THREADS="$2"; shift 2;;
    --memory-gb)
      MEMORY_GB="$2"; shift 2;;
    --pilon-rounds)
      PILON_ROUNDS="$2"; shift 2;;
    --pilon-bin)
      PILON_BIN="$2"; shift 2;;
    --skip-hybridspades)
      SKIP_HYBRIDSPADES=true; shift 1;;
    -h|--help)
      usage; exit 0;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1;;
  esac
 done

if [[ -z "${ONT_READS:-}" || -z "${PREFIX:-}" ]]; then
  echo "Missing required arguments." >&2
  usage
  exit 1
fi

if { [[ -n "${READ1:-}" ]] && [[ -z "${READ2:-}" ]]; } || { [[ -n "${READ2:-}" ]] && [[ -z "${READ1:-}" ]]; }; then
  echo "Both --short-reads1 and --short-reads2 must be provided together." >&2
  exit 1
fi

mkdir -p "$(dirname "$PREFIX")"
PREFIX_BASENAME="$(basename "$PREFIX")"
export PREFIX_BASENAME

log_msg() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

align_reads() {
  local asm="$1"
  SHORT_BAM=""
  SANGER_BAM=""

  if [[ -n "${READ1:-}" ]]; then
    log_msg "Aligning short reads to ${asm} for polishing/validation."
    minimap2 -ax sr -t "$THREADS" "$asm" "$READ1" "$READ2" | \
      samtools sort -@ "$THREADS" -o "${PREFIX}.alignments.short.bam" -
    samtools index "${PREFIX}.alignments.short.bam"
    SHORT_BAM="${PREFIX}.alignments.short.bam"
  fi

  if [[ -n "${SANGER:-}" ]]; then
    log_msg "Aligning Sanger references to ${asm} for polishing/validation."
    minimap2 -ax sr -t "$THREADS" "$asm" "$SANGER" | \
      samtools sort -@ "$THREADS" -o "${PREFIX}.alignments.sanger.bam" -
    samtools index "${PREFIX}.alignments.sanger.bam"
    SANGER_BAM="${PREFIX}.alignments.sanger.bam"
  fi
}

log_msg "Starting Flye assembly from ONT reads."
FLYE_DIR="${PREFIX}.flye"
flye --nano-raw "$ONT_READS" --genome-size "$GENOME_SIZE" --threads "$THREADS" --out-dir "$FLYE_DIR"
ASM_PATH="${FLYE_DIR}/assembly.fasta"

if [[ -n "${READ1:-}" && "$SKIP_HYBRIDSPADES" == false ]]; then
  log_msg "Running hybridSPAdes scaffolding with ONT + short reads."
  HYBRID_DIR="${PREFIX}.hybridspades"
  spades.py --threads "$THREADS" --memory "$MEMORY_GB" --careful --nanopore "$ONT_READS" -1 "$READ1" -2 "$READ2" -o "$HYBRID_DIR"
  if [[ -f "${HYBRID_DIR}/scaffolds.fasta" ]]; then
    ASM_PATH="${HYBRID_DIR}/scaffolds.fasta"
  else
    log_msg "hybridSPAdes did not produce scaffolds.fasta; retaining Flye assembly."
  fi
fi

if [[ -n "${READ1:-}" || -n "${SANGER:-}" ]]; then
  for round in $(seq 1 "$PILON_ROUNDS"); do
    log_msg "Pilon polishing round ${round} using current assembly ${ASM_PATH}."
    align_reads "$ASM_PATH"

    PILON_OUTDIR="${PREFIX}.pilon_round${round}"
    mkdir -p "$PILON_OUTDIR"

    PILON_CMD=("${PILON_BIN}" --genome "$ASM_PATH" --threads "$THREADS" --output "pilon_round${round}" --outdir "$PILON_OUTDIR" --changes)
    if [[ -n "$SHORT_BAM" ]]; then
      PILON_CMD+=(--frags "$SHORT_BAM")
    fi
    if [[ -n "$SANGER_BAM" ]]; then
      PILON_CMD+=(--unpaired "$SANGER_BAM")
    fi

    "${PILON_CMD[@]}"
    ASM_PATH="${PILON_OUTDIR}/pilon_round${round}.fasta"
  done
fi

log_msg "Generating final alignments for QC to ${ASM_PATH}."
align_reads "$ASM_PATH"

log_msg "Calling variants and building IUPAC consensus to preserve minor alleles."
if [[ -n "$SHORT_BAM" ]]; then
  POLISH_BAM="$SHORT_BAM"
elif [[ -n "$SANGER_BAM" ]]; then
  POLISH_BAM="$SANGER_BAM"
else
  echo "No alignments available for variant calling; skipping QC." >&2
  exit 0
fi

bcftools mpileup -Ou -a AD,ADF,ADR,DP -f "$ASM_PATH" "$POLISH_BAM" | \
  bcftools call -mv --ploidy 1 --keep-alts --multiallelic-caller -Oz -o "${PREFIX}.variants.vcf.gz"
bcftools index "${PREFIX}.variants.vcf.gz"
bcftools consensus --iupac-codes -f "$ASM_PATH" "${PREFIX}.variants.vcf.gz" > "${PREFIX}.consensus.iupac.fasta"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t[%AD]\n' "${PREFIX}.variants.vcf.gz" > "${PREFIX}.variant_retention.tsv"
bcftools stats "${PREFIX}.variants.vcf.gz" > "${PREFIX}.variants.stats.txt"

log_msg "Computing coverage metrics."
samtools coverage "$POLISH_BAM" > "${PREFIX}.coverage.tsv"
samtools depth -a "$POLISH_BAM" > "${PREFIX}.coverage.by_pos.tsv"
if command -v gnuplot >/dev/null 2>&1; then
  gnuplot <<GNUPLOT
set terminal pngcairo size 1200,400
set output "${PREFIX}.coverage.png"
set datafile separator "\t"
set title "Depth across assembly"
set xlabel "Position"
set ylabel "Depth"
plot "${PREFIX}.coverage.by_pos.tsv" using 2:3 with lines title "Coverage"
GNUPLOT
fi

log_msg "Summarizing assembly statistics."
python3 - <<PY
import os
from pathlib import Path
from statistics import median

def load_lengths(fasta_path):
    lengths = []
    with open(fasta_path) as handle:
        seq_len = 0
        for line in handle:
            if line.startswith('>'):
                if seq_len:
                    lengths.append(seq_len)
                seq_len = 0
            else:
                seq_len += len(line.strip())
        if seq_len:
            lengths.append(seq_len)
    return lengths

asm = Path("$ASM_PATH")
prefix_base = os.environ.get("PREFIX_BASENAME", "assembly")
lengths = load_lengths(asm)
lengths.sort(reverse=True)
if not lengths:
    raise SystemExit("No contigs found in assembly")
total = sum(lengths)
med = median(lengths)

n50 = 0
cumsum = 0
for l in lengths:
    cumsum += l
    if cumsum >= total / 2:
        n50 = l
        break

with open(asm.parent / f"{prefix_base}.qc.summary.txt", "w") as out:
    out.write(f"Assembly: {asm}\n")
    out.write(f"Contigs: {len(lengths)}\n")
    out.write(f"Total_bp: {total}\n")
    out.write(f"Median_len: {med}\n")
    out.write(f"N50: {n50}\n")
PY

log_msg "Hybrid assembly + validation complete. Final assembly: ${ASM_PATH}".
