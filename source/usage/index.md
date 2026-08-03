# Usage

ParaChrSNP is controlled through `config.yaml`. The workflow reads the sample
list and chromosome list from this file, constructs a Snakemake DAG, and runs
each rule in the configured Singularity/Apptainer container.

## Complete configuration file

Start by downloading the complete configuration template. It contains the
core workflow settings and all nine optional analysis modules.

{download}`Download config.full.example.yaml <../_static/config.full.example.yaml>`

The example deliberately sets every optional module to `enabled: true` so that
all available settings are visible. Before running ParaChrSNP, replace the
reference, sample prefixes, chromosome identifiers, population file, genome
annotation, executable paths, and resource values. Disable modules whose
inputs or software are not available.

```{literalinclude} ../_static/config.full.example.yaml
:language: yaml
:linenos:
```

After saving the edited file as `config.yaml`, continue with the input-file
requirements and execution methods below. Detailed explanations of every
optional module are provided in the
[optional-module guide](optional_modules.md).

## Input files

### Reference genome

Set the reference FASTA path:

```yaml
reference: "reference/Arabidopsis_thaliana.fasta"
```

ParaChrSNP creates or checks the required FASTA, sequence-dictionary, and
aligner indexes. Do not switch reference assemblies between workflow steps.

### Paired FASTQ files

Sample values are FASTQ prefixes without `.1.fq.gz` or `.2.fq.gz`.

```yaml
samples:
    sample1: "raw_fastq/sample1"
    sample2: "raw_fastq/sample2"
```

This configuration resolves to:

```text
raw_fastq/sample1.1.fq.gz
raw_fastq/sample1.2.fq.gz
raw_fastq/sample2.1.fq.gz
raw_fastq/sample2.2.fq.gz
```

### Chromosomes

List only chromosomes or contigs that should be called.

```yaml
chromosomes:
  - Chr01
  - Chr02
  - Chr03
```

The strings are case-sensitive and must match the FASTA identifiers exactly.

## Alignment backend

ParaChrSNP supports `minibwa`, `bwa-mem2`, and legacy `bwa` through one
configuration block.

```yaml
params:
    aligner:
        name: "minibwa"
        executable: "minibwa"
        map_threads: 8
        sort_threads: 4
        index_threads: 8
        index_extra: ""
        map_extra: ""
```

The mapping stream is passed directly through samblaster and SAMtools:

```text
minibwa map → samblaster --removeDups → samtools sort
```

The final alignment products are:

```text
duplicate_removed/{sample}.rmdup.bam
duplicate_removed/{sample}.rmdup.bam.bai
duplicate_removed/{sample}.dup.txt
```

When switching between minibwa and BWA-MEM2, reference re-indexing is expected
because the index formats are different.

## Chromosome-wise variant calling

GATK HaplotypeCaller runs independently for every sample and chromosome:

```text
gvcf/{sample}.{chromosome}.g.vcf.gz
```

This layout allows Snakemake to distribute calling jobs across CPUs or cluster
nodes without requiring one chromosome to wait for another.

## Joint-calling methods

Select the strategy with `params.joint_calling.method`.

### GenomicsDB

```yaml
params:
    joint_calling:
        method: "genomicsdb"
        reader_threads: 16
        genotype_threads: 1
        batch_size: 50
        import_java_options: '--java-options "-Xms1g -Xmx16g"'
        genotype_java_options: '--java-options "-Xms512m -Xmx128g"'
```

This method imports per-chromosome GVCFs into GenomicsDB and runs
`GenotypeGVCFs` by chromosome. It is the default GATK-oriented method for
multi-sample cohorts.

### GLnexus

```yaml
params:
    joint_calling:
        method: "glnexus"
        glnexus_executable: "scripts/glnexus_cli"
        glnexus_config: "gatk"
        glnexus_threads: 16
        glnexus_concat_threads: 32
```

GLnexus runs independently by chromosome and the resulting chromosome VCFs
are concatenated into `result_vcfs/combined.vcf.gz`. The executable must exist
at the configured path and must have execution permission.

### CombineGVCFs

```yaml
params:
    joint_calling:
        method: "combine_gvcfs"
```

This legacy mode is preserved for compatibility and small comparisons. It is
usually less attractive for large cohorts than GenomicsDB or GLnexus.

## SNP and INDEL filtering

SNPs and INDELs are selected separately and processed with GATK hard filters.
Any condition joined by `||` is sufficient to mark a site as low quality.

### SNP thresholds

```text
QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 ||
MQRankSum < -12.5 || ReadPosRankSum < -8.0
```

### INDEL thresholds

```text
QD < 2.0 || FS > 200.0 || SOR > 10.0 ||
MQRankSum < -12.5 || ReadPosRankSum < -8.0
```

`VariantFiltration` first marks failed sites. `SelectVariants
--exclude-filtered` then removes those sites from the final filtered VCFs.
These filters operate on site-level annotations; they do not perform MAF,
sample missingness, per-genotype DP, or per-genotype GQ filtering.

## Optional analyses

Optional modules are enabled independently in `config.yaml`. The table below
is a quick selection guide; the [optional-module guide](optional_modules.md)
describes each module's biological purpose, input dependencies, outputs,
parameters, and interpretation risks. A downloadable
[complete configuration example](full_configuration.md) contains the core
workflow and all optional modules.

| Module | Configuration flag | Use it to answer | Required workflow input |
| --- | --- | --- | --- |
| Beagle | `params.imputation.enabled` | Which missing genotypes can be inferred from haplotype information? | Filtered cohort SNP VCF |
| CNVnator | `params.cnv.enabled` | Which regions have read-depth evidence for deletions or duplications? | Per-sample deduplicated BAM files |
| Pi | `params.pi.enabled` | How does within-population nucleotide diversity vary across windows? | Filtered SNP VCF; optional population file |
| SNP density | `params.snp_density.enabled` | Where are SNP-rich and SNP-poor genomic windows? | Filtered SNP VCF and reference index |
| ADMIXTURE | `params.admixture.enabled` | What ancestry-component proportions describe each sample across K values? | PLINK files produced by `vcf_convert` |
| PopLDdecay | `params.ld_decay.enabled` | How quickly does linkage disequilibrium decay with physical distance? | Filtered SNP VCF; optional population file |
| PCA | `params.vcf2pca.enabled` | What are the main axes of genetic variation among samples? | Filtered SNP VCF; at least three samples |
| VCF2Dis | `params.vcf2dis.enabled` | What are the pairwise genetic distances and sample relationships? | Filtered SNP VCF; at least three samples |
| SnpEff | `params.snpeff.enabled` | Which genes, transcripts, and molecular consequences are affected? | Matching reference FASTA and GFF3/GTF annotation |

Group-aware analyses accept a two-column population file:

```text
sample1    Population_A
sample2    Population_A
sample3    Population_B
```

Sample identifiers must match those in `samples` and the VCF header.

```{toctree}
:maxdepth: 1
:caption: Detailed usage

optional_modules
full_configuration
```

## Running selected targets

Run only the configuration precheck:

```bash
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --cores 1 \
    --use-singularity \
    reports/precheck.done
```

Generate only the filtered cohort SNP VCF:

```bash
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --cores 64 \
    --use-singularity \
    result_vcfs/combined.snp.filtered.vcf.gz
```

Unlock a directory after confirming that no Snakemake process is running:

```bash
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --unlock
```

## External data directories

Bind external FASTQ and reference directories into the container:

```bash
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --cores 64 \
    --use-singularity \
    --singularity-args \
    "-B /data/fastq:/data/fastq -B /data/reference:/data/reference"
```

Paths used in `config.yaml` must be visible inside the container after binding.

## Web interface

Start the server-side launcher:

```bash
python web/parachrsnp_web.py \
    --host 0.0.0.0 \
    --port 8088 \
    --allowed-root /home/majunpeng/sda2
```

Open `http://SERVER_IP:8088` in a browser. The interface references existing
server-side data; it does not upload large sequencing files.

## Output directories

| Directory | Description |
| --- | --- |
| `qc/` | Raw-read quality-control results |
| `clean_reads/` | fastp-cleaned FASTQ files |
| `duplicate_removed/` | Deduplicated and sorted BAM files |
| `gvcf/` | Per-sample, per-chromosome GVCFs |
| `genomicsdb/` | GenomicsDB workspaces |
| `result_vcfs/` | Cohort, SNP, INDEL, and filtered VCFs |
| `missing/` | PLINK missingness tables and figures |
| `format_convert/` | PLINK and HapMap files |
| `imputation/` | Beagle output |
| `cnv/` | CNV calls, summaries, and figures |
| `pca/` | PCA coordinates and plots |
| `dis/` | Distance matrix and Newick tree |
| `admixture/` | Q matrices, CV statistics, and plots |
| `ld_decay/` | LD statistics and decay plots |
| `pi/` | Windowed Pi tables and plots |
| `snp_density/` | SNP-density tables and plots |
| `annotation/` | SnpEff database and annotated VCFs |
| `reports/` | Precheck and integrated HTML reports |

The final summary products are `reports/ParaChrSNP_report.html` and
`reports/ParaChrSNP_summary.tsv`.
