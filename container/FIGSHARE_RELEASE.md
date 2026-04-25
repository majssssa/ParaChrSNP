# ParaChrSNP Figshare Container Release

This document describes how users should run ParaChrSNP after downloading the
container image from Figshare.

## Required Files

Users need two things:

```text
1. The ParaChrSNP GitHub repository
2. The ParaChrSNP.sif container image downloaded from Figshare
```

The container already includes:

```text
snakemake
bwa
samtools
bcftools
gatk
picard
fastp
plink
TASSEL/run_pipeline.pl
rastqc
VCF2Dis
VCF2PCACluster
Python, pandas, numpy, matplotlib
```

## Project Inputs

The user still needs to prepare project-specific input files:

```text
raw_fastq/
reference/
pop.info
config.yaml
```

FASTQ files should follow the naming convention used in `config.yaml`:

```text
raw_fastq/{sample}.1.fq.gz
raw_fastq/{sample}.2.fq.gz
```

## Recommended Run Command

Check the workflow without running jobs.

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 4 --use-singularity -n

# snakemake: Run the ParaChrSNP Snakemake workflow.
# --snakefile Snakefile: Use the workflow entry file in the GitHub repository.
# --configfile config.yaml: Use the user-edited project config file.
# --cores 4: Allow Snakemake to use 4 CPU cores.
# --use-singularity: Run workflow jobs inside the container defined in Snakefile/config.yaml.
# -n: Dry-run mode; print planned jobs without executing them.
```

Run the complete workflow.

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 30 --use-singularity --keep-going

# --cores 30: Allow Snakemake to use up to 30 CPU cores.
# --use-singularity: Run every rule inside ParaChrSNP.sif.
# --keep-going: Continue running independent jobs after a job fails.
```

If FASTQ files or the reference genome are outside the GitHub project
directory, bind those external directories through Singularity arguments.

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 30 --use-singularity --singularity-args "-B /data/fastq:/data/fastq -B /data/reference:/data/reference"

# --singularity-args: Pass extra bind options directly to Singularity.
# -B /data/fastq:/data/fastq: Mount an external FASTQ directory into the same path inside the container.
# -B /data/reference:/data/reference: Mount an external reference directory into the same path inside the container.
```

Unlock a stale Snakemake directory.

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 1 --use-singularity --unlock

# --unlock: Remove a stale Snakemake lock after an interrupted run.
```

## Important Notes

The user needs Snakemake and Singularity/Apptainer available on the host system.
The host Snakemake process schedules jobs, and each job is executed inside
`ParaChrSNP.sif`.

The project `config.yaml` contains:

```yaml
container:
    image: "ParaChrSNP.sif"
```

Place `ParaChrSNP.sif` in the GitHub project root, or edit this path to point to
the downloaded Figshare image.

All paths in `config.yaml` should be relative to the project directory whenever
possible.

Absolute paths are also allowed, but the corresponding host directories must be
mounted with `--singularity-args "-B host_path:container_path"`.

Do not put raw FASTQ files, BAM files, GVCFs, VCFs, or result directories into
the GitHub repository. These files should stay on the user's filesystem and be
mounted at runtime.
