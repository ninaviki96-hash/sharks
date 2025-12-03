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

## MiSeq-to-ONT contig validation pipeline

Use `scripts/miseq_to_ont.sh` to trim adapters, optionally down-sample, align MiSeq reads to ONT contigs with either minimap2 (`-ax sr`) or BWA-MEM, flag duplicates without removing them, and emit coverage/variant summaries to verify scaffold accuracy in variable regions.

### Requirements
- `cutadapt`
- `minimap2` or `bwa`
- `samtools`
- `bcftools`

### Recommended parameters for high-sensitivity validation
- Adapter: Illumina TruSeq (`AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`)
- Quality/length filters: `--min-quality 20`, `--min-length 100` (keeps longer MiSeq inserts while removing short noisy tails)
- Aligner: `--aligner minimap2` (default) for tolerant short-read mapping to ONT contigs; use `--aligner bwa` with permissive penalties for problematic regions:
  - `--mismatch-penalty 3` (drop to 2 for SHM-rich segments)
  - `--gap-open-penalty 6`, `--gap-extend-penalty 1`
  - `--clip-penalty 5` (reduce to 3 to retain terminal variability)
- Coverage control: `--downsample-fraction 0.25 --downsample-seed 42` to cap extreme depth while preserving random representation of alleles.
- Duplicate handling: duplicates are marked (not removed) via `samtools markdup -s` to avoid over-amplification bias but still retain depth supporting minor alleles.
- Variant calling: haploid `bcftools call --ploidy 1 --keep-alts --multiallelic-caller` to keep alternate alleles and allele depths.

### Outputs
Using an output prefix like `results/sample1`, the pipeline writes:
- `results/sample1.sorted.markdup.bam` (+ `.bai`) with duplicates flagged
- `results/sample1.vcf.gz` (+ `.tbi`) with allele depths
- `results/sample1.coverage.txt` from `samtools coverage` to spot uneven depth
- `results/sample1.variant_summary.tsv` (chrom, pos, ref, alt, QUAL, depth, AD)
- `results/sample1.vcf.stats.txt` from `bcftools stats` for quick variant QC

### Example command

```bash
./scripts/miseq_to_ont.sh \
  -r ont_contigs/reference.fa \
  -1 data/sample_R1.fastq.gz \
  -2 data/sample_R2.fastq.gz \
  -o results/sample1 \
  --threads 8 \
  --aligner minimap2 \
  --min-length 100 \
  --downsample-fraction 0.25 \
  --clip-penalty 3
```

Outputs will include depth/variant summaries alongside BAM/VCF files to help confirm ONT scaffold correctness in variable regions without collapsing minor alleles.

## ONT-to-Sanger alignment and polishing pipeline

Use `scripts/ont_to_sanger.sh` to quality-filter ONT reads (via `filtlong`), align them to Sanger references with minimap2 `-ax map-ont` using reduced clipping penalties, call variants with allele depths preserved, and optionally polish a haploid consensus with racon and/or medaka while retaining an IUPAC consensus that reflects minor alleles.

### Requirements
- `minimap2`
- `filtlong` (skip with `--skip-filtlong`)
- `samtools`
- `bcftools`
- Optional: `racon` and `medaka` for polishing

### Parameter guidance for common ONT datasets
- Read length filtering: `--min-length 700-1000` and `--keep-percent 85-95` work well for typical ONT R9.4.1 runs with N50 around 8–15 kb while keeping high-quality subreads for amplicons of 600 bp–2 kb.
- Alignment scoring to reduce over-clipping at variable ends:
  - `--mm2-mismatch 3-4` (lower values retain more mismatches instead of clipping)
  - `--mm2-gap-open 6,18-24` and `--mm2-gap-extend 2,1-2` to discourage aggressive soft-clipping across indels
  - `--mm2-end-bonus 8-12` and `--mm2-zdrop 200,400` to reward full-length alignments on shorter amplicons.
- Polishing: `--racon-rounds 1-2` provides a fast haploid polish; add `--medaka-model r941_min_high_g360` (or a model matching your flowcell/kit) to refine homopolymers after racon. The VCF/IUPAC consensus is still derived from the alignment to preserve low-frequency alleles for interpretation.

### Example commands

**Amplicons ~900 bp with ONT reads N50 ≈ 10 kb:**
```bash
./scripts/ont_to_sanger.sh \
  -r references/IGHV_Sanger.fa \
  -i data/ont_amplicons.fastq.gz \
  -o results/sample1 \
  --threads 8 \
  --min-length 800 \
  --keep-percent 90 \
  --mm2-mismatch 3 \
  --mm2-end-bonus 12 \
  --racon-rounds 1
```

**Longer (~2 kb) amplicons with higher ONT N50 (~15 kb) and medaka polishing:**
```bash
./scripts/ont_to_sanger.sh \
  -r references/Sanger_targets.fa \
  -i data/ont_library.fastq.gz \
  -o results/sample2 \
  --threads 12 \
  --min-length 1000 \
  --keep-percent 85 \
  --mm2-gap-open 6,18 \
  --mm2-zdrop 250,500 \
  --racon-rounds 2 \
  --medaka-model r941_min_high_g360
```

Outputs include `results/sampleX.sorted.bam` (+ index) aligned with tuned clipping penalties, `results/sampleX.vcf.gz` with allele depths to retain minor variants, `results/sampleX.consensus.fasta` with IUPAC symbols, and optional racon/medaka polished FASTAs.
