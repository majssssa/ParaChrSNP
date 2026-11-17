# VCF-only downstream analysis (v0.1.1)

`Snakefile.downstream` starts from an existing genotyped cohort VCF. It is an
independent workflow, not an input-mode switch in the variant-calling pipeline.
Its configuration contains no alignment, BAM-processing or variant-calling
settings. The original `Snakefile` remains the FASTQ/BAM entry point.

## Input and configuration

{download}`Download the complete downstream configuration <../../config.downstream.example.yaml>`

Use `.vcf` or `.vcf.gz` containing uniquely named samples and FORMAT/GT.
The matching reference FASTA is required. Sample IDs are inferred from the VCF
header; no `Types`, `samples` or `chromosomes` settings are needed.
Precheck rejects missing headers, duplicate sample IDs, missing GT declarations
or records, gVCF reference-block alleles and record contigs absent from the FASTA.
Contig-name validation does not prove sequence identity with the reference.

The default targets are missingness, PLINK binary/text and HapMap conversion,
PCA, genetic distances/tree, and HTML/TSV reporting. PCA/tree require at least
three samples. Other modules are disabled in the release template. Enable
`params.imputation`, `pi`, `snp_density`, `admixture`, `ld_decay` or `snpeff` with
their `enabled` settings as needed; annotation requires matching FASTA/GFF
resources. Set ADMIXTURE's K range and cross-validation for your sample count.
See [optional modules](optional_modules.md) for resource requirements.
VCF-only analysis cannot run BAM-based CNVnator.

## Keep existing results separate

Run the following setup from the cloned repository root. Keep the repository's
`rules/` and `scripts/`; the two entry/config files are not a standalone package.

Create a new working directory and copy the complete configuration.

```bash
mkdir -p downstream_run
cp config.downstream.example.yaml downstream_run/config.yaml
ln -s "$PWD/scripts" downstream_run/scripts
# mkdir -p: create the output directory if necessary.
# cp: copy the complete template for local editing.
# ln -s: expose the repository's scripts in the new working directory.
# $PWD: the absolute repository root; run setup from that directory.
```

Edit `downstream_run/config.yaml`: set the top-level `input_vcf`, `reference`
and `container.image` to **absolute paths** to your files. Paths are otherwise
resolved relative to the working directory. Keep the module output paths
relative so they remain inside `downstream_run`. Update any population-file,
annotation or Beagle resource paths if those modules are enabled.

The input can be an existing `result_vcfs/combined.snp.filtered.vcf.gz` in the
original project because outputs now go into the separate working directory.
Do not use the same file as both input and workflow output; the workflow rejects
that collision rather than overwriting the input.

## Check and run

Validate the task graph after editing the configuration.

```bash
env TMPDIR=/tmp snakemake --snakefile Snakefile.downstream \
  --configfile downstream_run/config.yaml --directory downstream_run \
  --cores 1 -n
# TMPDIR=/tmp: use an accessible temporary directory.
# --snakefile: select the independent downstream entry.
# --configfile: load the edited complete configuration.
# --directory: place outputs in the separate working directory.
# --cores 1: allocate one core for validation.
# -n: construct the DAG without executing jobs; content precheck runs later.
```

Run the enabled analyses with explicit container bindings and environment isolation.

```bash
env TMPDIR=/tmp snakemake --snakefile Snakefile.downstream \
  --configfile downstream_run/config.yaml --directory downstream_run \
  --cores 8 --use-singularity --rerun-incomplete \
  --singularity-args "--cleanenv --bind $PWD:$PWD --env TMPDIR=/tmp --env PYTHONNOUSERSITE=1 --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null --env R_LIBS_USER=/nonexistent"
# --cores 8: allow eight CPU threads across concurrent jobs.
# --use-singularity: use the configured Singularity/Apptainer image per rule.
# --rerun-incomplete: rerun interrupted/incomplete jobs.
# --singularity-args: forward the following options to the container runtime.
# --cleanenv: avoid inheriting unrelated host environment variables.
# --bind $PWD:$PWD: expose the repository at the same path inside the container.
# --env TMPDIR=/tmp: choose a container-visible temporary directory.
# --env PYTHONNOUSERSITE=1: disable host Python user packages.
# --env R_ENVIRON_USER=/dev/null: skip the host R user environment file.
# --env R_PROFILE_USER=/dev/null: skip the host R user startup script.
# --env R_LIBS_USER=/nonexistent: exclude host R user libraries.
# --snakefile, --configfile and --directory: same meanings as in the dry-run.
```

Bind external input/reference/annotation directories too, for example by adding
`--bind /data:/data` inside `--singularity-args` when those files are under
`/data`. Absolute paths alone do not make host files visible in the container.
`No such file or directory` from faidx usually indicates a missing bind;
`mkdtemp ... failed` can indicate an inaccessible temporary directory.

## Results and rerunning

The original VCF is preserved. A sorted BGZF copy is staged under
`downstream_run/result_vcfs/`. SNP analyses retain biallelic SNP records whose
FILTER is `PASS` or `.`; INDEL records with those FILTER values are exported
separately. No additional GATK hard filters are applied. A SNP-only input may
therefore produce an empty INDEL file. Multiallelic SNPs are excluded, not split.

Open `downstream_run/reports/ParaChrSNP_report.html` for results and figures;
the companion `ParaChrSNP_summary.tsv` contains summary metrics. Read QC and
read-derived sequencing depth are unavailable in VCF mode. Tree/matrix outputs
retain the full sample names, including names longer than ten characters.

Rerun the same command to resume failed jobs; completed outputs are reused.
Removing only the final report outputs rebuilds the report without rerunning
completed analyses. Start a new working directory when changing datasets.
