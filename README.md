# sharks

cDNA processing kit for immunology.

## MiSeq → ONT assembly + alignment pipeline

Use `scripts/assemble_and_map_miseq.sh` when you have raw ONT long reads (FASTQ/FASTA) and paired-end MiSeq data for the same target. The script first assembles the ONT reads into contigs (Flye by default, Miniasm+Racon as an option) and then trims + aligns the MiSeq reads to those contigs with duplicate marking, coverage/variant summaries, and an IUPAC consensus ready for downstream phylogenetics. Assembly is skipped automatically if `<prefix>.assembled_contigs.fasta` already exists.

### Requirements
- `cutadapt`
- `flye` **or** `minimap2 + miniasm + racon`
- `bwa` **or** `minimap2`
- `samtools`
- `bcftools`

### Recommended parameters for high-sensitivity validation
The defaults target high-sensitivity MiSeq alignment against noisy ONT assemblies while keeping rare variants:
- Adapter: Illumina TruSeq (`AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`)
- Quality/length filters: `--min-quality 20`, `--min-length 100`
- Aligner: `--aligner bwa` (default) with permissive penalties for variable regions:
  - `--mismatch-penalty 3` (drop to 2 for highly mutated loci)
  - `--gap-open-penalty 6`, `--gap-extend-penalty 1`
  - `--clip-penalty 5` (reduce to 3 to retain terminal variability)
  Use `--aligner minimap2` for tolerant short-read mapping to long contigs when BWA-MEM struggles.
- Coverage control: `--downsample-fraction 0.25 --downsample-seed 42` to cap extreme depth while preserving random representation of alleles.
- Duplicate handling: duplicates are marked (not removed) via `samtools markdup -s` to avoid bias while preserving depth for minor alleles.
- Variant calling: haploid `bcftools call --ploidy 1 --keep-alts --multiallelic-caller` to keep alternate alleles and allele depths.

### Example command
For 2×250 bp MiSeq reads aligned to contigs assembled from ONT FASTQ:

```bash
./scripts/assemble_and_map_miseq.sh \
  -r data/ont_raw.fastq.gz \
  -1 data/sample_R1.fastq.gz \
  -2 data/sample_R2.fastq.gz \
  -o results/sample1 \
  --threads 12 \
  --assembler flye \
  --min-length 100
```

Use an output prefix such as `results/sample1` without a trailing slash; the scripts will normalize the prefix and create parent directories automatically so files land beside the chosen prefix rather than as hidden dotfiles.

Outputs include `results/sample1.assembled_contigs.fasta` (if built), `results/sample1.sorted.markdup.bam` (+ index), `results/sample1.alignment.flagstat.txt`, `results/sample1.vcf.gz` with allele depths, `results/sample1.coverage.tsv`, `results/sample1.variant_summary.tsv`, `results/sample1.vcf.stats.txt`, and `results/sample1.consensus.fasta` containing IUPAC symbols where minor alleles are detected.

## MiSeq → Sanger unified pipeline

Use `scripts/miseq_to_sanger_unified.sh` to trim adapters, quality-filter MiSeq reads, align them to Sanger reference sequences with BWA-MEM (default) or minimap2, and emit the same BAM/coverage/variant/consensus artifacts as the ONT workflow. Only the MiSeq reads are trimmed; the Sanger references are used as-is.

### Requirements
- `cutadapt`
- `bwa` or `minimap2`
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
./scripts/miseq_to_sanger_unified.sh \
  -r references/IGHV_Sanger.fa \
  -1 data/sample_R1.fastq.gz \
  -2 data/sample_R2.fastq.gz \
  -o results/sample1 \
  --threads 8 \
  --mismatch-penalty 2 \
  --clip-penalty 3 \
  --min-length 120
```

Outputs will include `results/sample1.sorted.markdup.bam` (+ index), `results/sample1.alignment.flagstat.txt`, `results/sample1.vcf.gz` with allele depths, `results/sample1.coverage.tsv`, `results/sample1.variant_summary.tsv`, `results/sample1.vcf.stats.txt`, and `results/sample1.consensus.fasta` containing IUPAC symbols where minor alleles are detected.

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
