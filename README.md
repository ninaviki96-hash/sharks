# sharks

cDNA processing kit for immunology.

## MiSeq-to-Sanger alignment pipeline

Use `scripts/miseq_to_sanger.sh` to trim adapters, quality-filter MiSeq reads, align them to Sanger reference sequences, and emit both sorted BAM files and IUPAC-coded consensus FASTAs so low-frequency alleles remain represented instead of being collapsed away.

### Requirements
- `cutadapt`
- `bwa`
- `samtools`
- `bcftools`

### Recommended parameters for immunoglobulin libraries
The defaults are tuned for typical 2×250 bp MiSeq runs against single-molecule Sanger references while keeping rare variants:
- Adapter: Illumina TruSeq (`AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`)
- Quality/length filters: `--min-quality 20`, `--min-length 120` (raises stringency for long amplicons)
- Alignment penalties to avoid over-collapsing somatic variants:
  - `--mismatch-penalty 2-3` (2 is more permissive for SHM-rich IGHV/IGKV)
  - `--gap-open-penalty 6`, `--gap-extend-penalty 1`
  - `--clip-penalty 3-5` (3 for aggressive retention of terminal variability; 5 if adapters/primers remain)
- Consensus: IUPAC output from `bcftools consensus --iupac-codes` so mixed alleles are preserved in the FASTA rather than forced to a single base.

### Example command
For paired-end 2×250 bp MiSeq reads from an IGHV amplicon aligned to a Sanger reference:

```bash
./scripts/miseq_to_sanger.sh \
  -r references/IGHV_Sanger.fa \
  -1 data/sample_R1.fastq.gz \
  -2 data/sample_R2.fastq.gz \
  -o results/sample1 \
  --threads 8 \
  --mismatch-penalty 2 \
  --clip-penalty 3 \
  --min-length 120
```

Outputs will include `results/sample1.sorted.bam` (+ index), `results/sample1.vcf.gz` with allele depths, and `results/sample1.consensus.fasta` containing IUPAC symbols where minor alleles are detected.
