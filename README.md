# sharks

cDNA processing kit for immunology.

## MiSeq-to-Sanger alignment pipeline

Use `scripts/miseq_to_sanger.sh` to trim adapters, quality-filter MiSeq reads, align them to Sanger reference sequences, and emit duplicate-marked BAM/coverage/variant summaries plus IUPAC-coded consensus FASTAs so low-frequency alleles remain represented instead of being collapsed away.

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

< codex/-miseq-t0x3ex
Use an output prefix such as `results/sample1` without a trailing slash; the scripts will normalize the prefix and create parent
directories automatically so files land beside the chosen prefix rather than as hidden dotfiles.

=======
< codex/-miseq-lvc8lw
Use an output prefix such as `results/sample1` without a trailing slash; the scripts will normalize the prefix and create parent
directories automatically so files land beside the chosen prefix rather than as hidden dotfiles.

=======
 codex/-miseq-dogg7t
Use an output prefix such as `results/sample1` without a trailing slash; the scripts will normalize the prefix and create parent
directories automatically so files land beside the chosen prefix rather than as hidden dotfiles.

=======
 codex/-miseq-r8twkl
Use an output prefix such as `results/sample1` without a trailing slash; the scripts will normalize the prefix and create parent
directories automatically so files land beside the chosen prefix rather than as hidden dotfiles.

=======
 main
 main
> main
> main
Outputs will include `results/sample1.sorted.markdup.bam` (+ index), `results/sample1.vcf.gz` with allele depths, `results/sample1.coverage.txt`, `results/sample1.variant_summary.tsv`, `results/sample1.vcf.stats.txt`, and `results/sample1.consensus.fasta` containing IUPAC symbols where minor alleles are detected.

## MiSeq-to-Sanger alignment for pre-trimmed reads

Use `scripts/miseq_trimmed_to_sanger.sh` when adapter/quality trimming has already been done. The script will build a BWA index for the provided Sanger FASTA if it is missing, align trimmed MiSeq reads with permissive penalties to retain rare variants, and output a sorted BAM, VCF with allele depths, and an IUPAC consensus FASTA.

### Requirements
- `bwa`
- `samtools`
- `bcftools`

### Example command

```bash
./scripts/miseq_trimmed_to_sanger.sh \
  -r references/IGHV_Sanger.fa \
  -1 data/sample_trimmed_R1.fastq.gz \
  -2 data/sample_trimmed_R2.fastq.gz \
  -o results/sample1 \
  --threads 8 \
  --mismatch-penalty 2 \
  --clip-penalty 3
```

Outputs will include `results/sample1.sorted.bam` (+ index), `results/sample1.vcf.gz` with allele depths, and `results/sample1.consensus.fasta` with IUPAC symbols representing minor alleles.

## MiSeq-to-ONT contig validation pipeline

Use `scripts/miseq_to_ont.sh` to trim adapters, optionally down-sample, align MiSeq reads to ONT contigs with BWA-MEM (default) or minimap2 (`-ax sr`), flag duplicates without removing them, and emit coverage/variant summaries to verify scaffold accuracy in variable regions.

### Requirements
- `cutadapt`
- `bwa` or `minimap2`
- `samtools`
- `bcftools`

### Recommended parameters for high-sensitivity validation
- Adapter: Illumina TruSeq (`AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`)
- Quality/length filters: `--min-quality 20`, `--min-length 100` (keeps longer MiSeq inserts while removing short noisy tails)
- Aligner: `--aligner bwa` (default) with permissive penalties for problematic regions:
  - `--mismatch-penalty 3` (drop to 2 for SHM-rich segments)
  - `--gap-open-penalty 6`, `--gap-extend-penalty 1`
  - `--clip-penalty 5` (reduce to 3 to retain terminal variability)
  Use `--aligner minimap2` for tolerant short-read mapping to ONT contigs when BWA-MEM struggles.
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
- `results/sample1.consensus.fasta` (IUPAC-coded so minor alleles persist rather than being collapsed)

### Example command

```bash
./scripts/miseq_to_ont.sh \
  -r ont_contigs/reference.fa \
  -1 data/sample_R1.fastq.gz \
  -2 data/sample_R2.fastq.gz \
  -o results/sample1 \
  --threads 8 \
  # --aligner minimap2 \\
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

## Hybrid ONT + MiSeq/Sanger assembly & validation pipeline

Use `scripts/hybrid_assembly.sh` to assemble ONT reads with Flye, optionally scaffold with hybridSPAdes when MiSeq pairs are available, polish against MiSeq/Sanger alignments via Pilon, and emit QC/coverage summaries that highlight depth and allele retention in immune variable regions.

```mermaid
graph LR
  A[ONT reads] -->|Flye| B[Flye contigs]
  B -->|optional| C[hybridSPAdes scaffolds]
  C -->|or| B
  C --> D[Pilon polish (MiSeq/Sanger alignments)]
  B --> D
  D --> E[IUPAC consensus + QC]
  E --> F[Coverage/variant retention reports]
```

### Requirements
- `flye`
- `minimap2`
- `samtools`
- `bcftools`
- `pilon` (Java)
- Optional: `spades.py` (hybridSPAdes), `gnuplot`

### Why this workflow preserves immune variability
- Flye tolerates ONT indels and retains long-range haplotypes across IGHV/IGKV regions.
- hybridSPAdes scaffolding uses accurate MiSeq pairs to resolve repeats without collapsing long-range ONT signal.
- Pilon polishing aligns MiSeq/Sanger reads back to the assembly while keeping alternate alleles for rare/low-frequency variants.
- Post-assembly validation emits coverage-by-position tables, bcftools variant summaries with allele depths, and optional coverage plots to flag regions where depth drops could hide rare alleles.

### Example commands

**ONT-only contigs with MiSeq+Sanger polishing/QC**
```bash
./scripts/hybrid_assembly.sh \
  -i data/ont_reads.fastq.gz \
  -1 data/miseq_R1.fastq.gz \
  -2 data/miseq_R2.fastq.gz \
  -s references/sanger_targets.fa \
  -o results/sample_hybrid \
  --genome-size 3m \
  --threads 12 \
  --pilon-rounds 2
```

**ONT+MiSeq scaffolds with hybridSPAdes followed by Pilon**
```bash
./scripts/hybrid_assembly.sh \
  -i data/ont_reads.fastq.gz \
  -1 data/miseq_R1.fastq.gz \
  -2 data/miseq_R2.fastq.gz \
  -o results/sample_scaffolds \
  --genome-size 3.2m \
  --threads 16 \
  --memory-gb 48 \
  --pilon-rounds 1
```

Key outputs live under the chosen prefix (e.g., `results/sample_hybrid.*`): Flye and hybridSPAdes assemblies, Pilon-polished FASTAs, short-read/Sanger BAMs with indexes, VCF + IUPAC consensus to retain alternate alleles, coverage TSVs/PNG, and a `*.qc.summary.txt` with contig counts/N50.
