# VCF-only downstream analysis

Use `Snakefile.downstream` to run downstream analyses from an existing cohort
VCF, without repeating alignment or variant calling.

## Configure

{download}`Download the complete configuration <../../config.downstream.example.yaml>`

Edit `config.downstream.example.yaml` and set `input_vcf`, `reference` and
`container.image` to your VCF, matching reference FASTA and container image.
Run commands from the repository root. Input files must not share the workflow's
output paths under `result_vcfs/`.

The input may be `.vcf` or `.vcf.gz` and must contain GT genotypes.
Sample names are read automatically; no separate sample list is needed.
PCA and tree analysis require at least three samples.

The default analyses are missingness, PLINK/HapMap conversion, PCA,
genetic distances/tree and a summary report. Enable ADMIXTURE, LD decay,
nucleotide diversity, SNP density, imputation or annotation through the
corresponding `params.*.enabled` settings.
See [optional modules](optional_modules.md) for details.

## Run

Check the configuration and task graph.

```bash
snakemake --snakefile Snakefile.downstream --configfile config.downstream.example.yaml --cores 1 -n
# --snakefile: select the independent downstream workflow.
# --configfile: load the edited complete configuration.
# --cores 1: use one core for the check.
# -n: inspect the task graph without running jobs.
```

Run the enabled analyses.

```bash
snakemake --snakefile Snakefile.downstream --configfile config.downstream.example.yaml --cores 8 --use-singularity
# --snakefile: select the independent downstream workflow.
# --configfile: load the edited complete configuration.
# --cores 8: allow up to eight CPU threads.
# --use-singularity: use the configured container.
```

## Results

Open `reports/ParaChrSNP_report.html` for the report and figures.
Summary metrics are in `reports/ParaChrSNP_summary.tsv`.
Individual analyses write to their configured output directories.

The original VCF is preserved. SNP analyses use biallelic records with FILTER
`PASS` or `.`; multiallelic SNPs are excluded. Read-based QC and sequencing
depth are unavailable when starting from VCF.
