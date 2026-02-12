#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: assemble_and_map_miseq.sh -r <ont_reads.(fa|fq|gz)> -1 <reads_R1.fastq.gz> -2 <reads_R2.fastq.gz> -o <output_prefix> [options]

Treat raw ONT reads (FASTA/FASTQ) as the reference, trim MiSeq reads, align with BWA-MEM (default) or minimap2 -ax sr, and emit BAM/VCF/consensus artifacts for phylogenetics.

Required arguments:
  -r, --ont-reads         ONT reads FASTA/FASTQ (optionally gzipped) used directly as the reference
  -1, --reads1            MiSeq R1 FASTQ (gzipped)
  -2, --reads2            MiSeq R2 FASTQ (gzipped)
  -o, --output-prefix     Prefix for outputs (e.g., results/sample1)

Optional arguments:
      --adapter           Adapter sequence to trim (default: Illumina TruSeq)
      --min-quality       Phred quality cutoff for 3' trimming (default: 20)
      --min-length        Discard reads shorter than this after trimming (default: 100)
      --threads           Number of CPU threads (default: 4)
      --aligner           "bwa" (default) or "minimap2" for MiSeq alignment
      --mismatch-penalty  BWA-MEM mismatch penalty -B (default: 3)
      --gap-open-penalty  BWA-MEM gap open penalties -O (default: 6)
      --gap-extend-penalty BWA-MEM gap extension penalties -E (default: 1)
      --clip-penalty      BWA-MEM clipping penalty -L (default: 5)
      --downsample-fraction Subsample alignments with samtools view -s (0-1, skipped if >=1 or unset)
      --downsample-seed   Seed for samtools -s (integer, default: 42)
  -h, --help              Show this help message

Outputs (using the provided prefix):
  <prefix>.ont_reference.fasta (FASTA copy of ONT reads for indexing)
  <prefix>.trimmed_R1.fastq.gz / <prefix>.trimmed_R2.fastq.gz
  <prefix>.sorted.markdup.bam (+ .bai) with duplicates flagged (not removed)
  <prefix>.alignment.flagstat.txt from samtools flagstat for phylogenetic QC
  <prefix>.vcf.gz (+ .tbi) with allele depths retained
  <prefix>.coverage.tsv from samtools coverage
  <prefix>.variant_summary.tsv with key variant metrics
  <prefix>.vcf.stats.txt from bcftools stats
  <prefix>.consensus.fasta (IUPAC-coded consensus so minor alleles are not collapsed)

Dependencies: cutadapt, bwa or minimap2, samtools, bcftools, and python3 must be installed and in PATH.
USAGE
}

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

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--ont-reads)
      ONT_READS="$2"; shift 2;;
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

if [[ -z "${ONT_READS:-}" || -z "${READ1:-}" || -z "${READ2:-}" || -z "${PREFIX:-}" ]]; then
  echo "Missing required arguments." >&2
  usage
  exit 1
fi

PREFIX="${PREFIX%/}"
if [[ -z "$PREFIX" ]]; then
  echo "Error: output prefix cannot be empty after normalization." >&2
  exit 1
fi

if [[ "$PREFIX" == */* ]]; then
  OUT_DIR_RAW="${PREFIX%/*}"
  PREFIX_BASENAME="${PREFIX##*/}"
else
  OUT_DIR_RAW="."
  PREFIX_BASENAME="$PREFIX"
fi

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

for tool in cutadapt samtools bcftools python3; do
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

for f in "$ONT_READS" "$READ1" "$READ2"; do
  if [[ ! -f "$f" ]]; then
    echo "Error: input file not found: $f" >&2
    exit 1
  fi
  if [[ ! -r "$f" ]]; then
    echo "Error: input file is not readable: $f" >&2
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
  DOWNSAMPLE_FMT="${DOWNSAMPLE_SEED}.${DOWNSAMPLE#0.}"
fi

REF_FASTA="${PREFIX}.ont_reference.fasta"

if [[ ! -s "$REF_FASTA" ]]; then
  echo "Preparing ONT reads as FASTA reference: $REF_FASTA"
  python3 - "$ONT_READS" "$REF_FASTA" <<'PY'
import gzip
import sys
from pathlib import Path

src = Path(sys.argv[1])
dest = Path(sys.argv[2])

def opener(path):
    if path.suffix == '.gz':
        return gzip.open(path, 'rt')
    return open(path, 'r')

with opener(src) as fh:
    first = fh.read(1)
    if not first:
        sys.stderr.write(f"Error: input file {src} is empty.\n")
        sys.exit(1)
    fh.seek(0)
    if first == '>':
        with open(dest, 'w') as out:
            for line in fh:
                out.write(line)
    elif first == '@':
        with open(dest, 'w') as out:
            while True:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline()
                plus = fh.readline()
                qual = fh.readline()
                if not qual:
                    sys.stderr.write("Error: FASTQ appears truncated.\n")
                    sys.exit(1)
                out.write('>' + header[1:])
                out.write(seq)
    else:
        sys.stderr.write("Error: reference must be FASTA or FASTQ (optionally gzipped).\n")
        sys.exit(1)
PY
fi

REF="$REF_FASTA"

if [[ "$ALIGNER" == "bwa" ]]; then
  NEED_INDEX=false
  for ext in amb ann bwt pac sa; do
    if [[ ! -f "${REF}.${ext}" ]]; then
      NEED_INDEX=true
      break
    fi
  done
  if [[ "$NEED_INDEX" == true ]]; then
    echo "BWA index not found for $REF. Building index..."
    bwa index "$REF"
  fi
fi

cutadapt \
  -j "$THREADS" \
  -a "$ADAPTER" -A "$ADAPTER" \
  -q "$MIN_QUAL" \
  -m "$MIN_LEN" \
  -o "${PREFIX}.trimmed_R1.fastq.gz" \
  -p "${PREFIX}.trimmed_R2.fastq.gz" \
  "$READ1" "$READ2"

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

if [[ -n "$DOWNSAMPLE" ]]; then
  samtools view -@ "$THREADS" -b -s "$DOWNSAMPLE_FMT" "$ALIGN_BAM" -o "${PREFIX}.downsampled.bam"
  ALIGN_BAM="${PREFIX}.downsampled.bam"
fi

samtools sort -n -@ "$THREADS" -o "${PREFIX}.namesort.bam" "$ALIGN_BAM"
samtools fixmate -m "${PREFIX}.namesort.bam" "${PREFIX}.fixmate.bam"
samtools sort -@ "$THREADS" -o "${PREFIX}.positionsort.bam" "${PREFIX}.fixmate.bam"
samtools markdup -@ "$THREADS" -s "${PREFIX}.positionsort.bam" "${PREFIX}.sorted.markdup.bam"
samtools index "${PREFIX}.sorted.markdup.bam"
samtools flagstat "${PREFIX}.sorted.markdup.bam" > "${PREFIX}.alignment.flagstat.txt"

rm -f "${PREFIX}.unsorted.bam" "${PREFIX}.downsampled.bam" "${PREFIX}.namesort.bam" "${PREFIX}.fixmate.bam" "${PREFIX}.positionsort.bam" 2>/dev/null || true

bcftools mpileup -Ou -a AD,ADF,ADR,DP -f "$REF" "${PREFIX}.sorted.markdup.bam" \
  | bcftools call -mv --ploidy 1 --keep-alts --multiallelic-caller -Oz -o "${PREFIX}.vcf.gz"
bcftools index "${PREFIX}.vcf.gz"

samtools coverage "${PREFIX}.sorted.markdup.bam" > "${PREFIX}.coverage.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t[%AD]\n' "${PREFIX}.vcf.gz" > "${PREFIX}.variant_summary.tsv"
bcftools stats "${PREFIX}.vcf.gz" > "${PREFIX}.vcf.stats.txt"
bcftools consensus --iupac-codes -f "$REF" "${PREFIX}.vcf.gz" > "${PREFIX}.consensus.fasta"

echo "Pipeline complete. Outputs written with prefix ${PREFIX}."
