# Changelog

This page summarizes user-visible ParaChrSNP changes reflected in the current
repository. Dates describe repository development milestones rather than a
formal semantic-version release series.

## 2026-09-28

### BAM input mode

- Added `Types: reads|bam`, with `reads` as the backward-compatible default.
- BAM mode discovers `{bam_dir}/{sample}.bam`, checks the input BAMs, links
  them into the existing duplicate-removed BAM interface, and indexes them.
- Skips FASTQ cleaning, alignment, samblaster and RastQC in BAM mode.
- Marks reads-only quality and depth metrics as unavailable in the report.

## 2026-08-18

### GLnexus workflow

- Removed the redundant post-joint-calling bcftools filter from GLnexus mode.
- GLnexus output is now separated directly into the final SNP and INDEL VCFs.
- GenomicsDB and CombineGVCFs continue to use the existing GATK hard filters.
- Retained the established `combined.*.filtered.vcf.gz` filenames for
  compatibility with downstream modules and existing configurations.

## 2026-08-03

### Documentation

- Added Read the Docs configuration and a Sphinx documentation project.
- Organized documentation into Installation, Quick start, Usage, Changelog,
  and FAQ sections.
- Added local dependency definitions for reproducible Read the Docs builds.
- Documented GLnexus-specific post-joint-calling filters and their defaults.

### GLnexus filtering

- Added a dedicated bcftools filter path for GLnexus pVCFs, which do not
  contain the annotations required by the GATK hard-filter expressions.
- Added configurable PASS status, biallelic-site, allele-quality, genotype
  depth, genotype quality, missingness, and minor-allele-count filters.
- Preserved the original GATK SNP and INDEL filtering behavior for GenomicsDB
  and CombineGVCFs joint calling.

### Alignment and duplicate removal

- Updated the supported minibwa build to `0.6-r416`.
- Integrated samblaster directly between alignment and SAMtools sorting.
- Replaced the separate GATK MarkDuplicates workflow path with streaming
  duplicate removal for the optimized alignment workflow.
- Added a GATK-compatible `PERCENT_DUPLICATION` summary derived from the
  samblaster log.

### Joint calling

- Added `glnexus` as a selectable value of `params.joint_calling.method`.
- Preserved `genomicsdb` and `combine_gvcfs` for user-selected compatibility.
- Implemented per-chromosome GLnexus calling followed by ordered bcftools
  concatenation.
- Kept `glnexus_cli` external to the container so users can download the
  executable separately into `scripts/glnexus_cli`.

### VCF compatibility

- Added `--vcf-half-call missing` to PLINK commands that read GLnexus VCFs.
- Half-called genotypes such as `./1` are treated as missing instead of causing
  PLINK 1.9 to terminate.

### Configuration handling

- Improved support for selecting an alternative file with `--configfile`.
- Propagated the active configuration file to precheck and reporting rules.
- Added a small `config.test.yaml` for workflow validation.

## 2026-06

- Added minibwa as an alignment backend alongside BWA-MEM2 and BWA.
- Refined thread allocation for alignment and sorting.
- Added benchmark workflows comparing chromosome-wise and traditional joint
  calling layouts.

## Existing workflow capabilities

- Chromosome-wise GATK HaplotypeCaller execution.
- GenomicsDB and legacy CombineGVCFs joint calling.
- Configurable GATK hard filtering for SNPs and INDELs.
- PLINK and HapMap format conversion.
- Optional Beagle imputation, CNVnator, ADMIXTURE, LD decay, Pi, SNP density,
  PCA, genetic-distance, phylogenetic-tree, and SnpEff modules.
- Integrated precheck and HTML reporting.

## Reporting changes

When reporting a problem, record the following information so changes can be
traced accurately:

```text
ParaChrSNP Git commit
Container filename and SHA256
Snakemake version
Selected aligner and version
Selected joint-calling method
Command line
Relevant rule log
```

Obtain the current Git commit with:

```bash
git rev-parse --short HEAD
```

Obtain an image checksum with:

```bash
sha256sum ParaChrSNP.sif
```
