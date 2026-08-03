# Installation

ParaChrSNP is distributed as a Snakemake workflow and is designed to run with
Singularity or Apptainer. The container is the recommended installation method
because it keeps GATK, minibwa, samblaster, SAMtools, PLINK, R, and the other
workflow dependencies in one reproducible environment.

## Requirements

Before installing ParaChrSNP, prepare:

- A 64-bit Linux system.
- Git for downloading the workflow source code.
- Singularity or Apptainer for running the container.
- Enough storage for FASTQ, BAM, GVCF, VCF, and temporary files.
- A reference genome in FASTA format.
- Paired-end short-read data in compressed FASTQ format.

The required CPU count, memory, and storage depend on genome size, sequencing
depth, sample number, and the selected joint-calling method. Large cohorts
should use a fast temporary directory with sufficient free space.

## Download the workflow

Clone the source repository and enter its root directory.

```bash
git clone https://github.com/majssssa/ParaChrSNP.git
cd ParaChrSNP
```

- `git clone`: downloads the complete ParaChrSNP repository.
- `cd ParaChrSNP`: enters the workflow root containing `Snakefile` and
  `config.yaml`.

## Download the container

Download the released Singularity/Apptainer image.

```bash
singularity pull \
    ParaChrSNP.sif \
    http://www.majunpeng.com/ParaChrSNP/ParaChrSNP.sif
```

- `singularity pull`: downloads a container image.
- `ParaChrSNP.sif`: local output name of the container.
- The final argument is the image download address.

Apptainer can be used instead of Singularity.

```bash
apptainer pull \
    ParaChrSNP.sif \
    http://www.majunpeng.com/ParaChrSNP/ParaChrSNP.sif
```

## Configure the container path

Set the container image in `config.yaml`.

```yaml
container:
    image: "ParaChrSNP.sif"
```

An absolute path can be used when the image is stored outside the project.

```yaml
container:
    image: "/data/containers/ParaChrSNP.sif"
```

## Optional sandbox conversion

For workflows with many short jobs, converting the image to a sandbox can
reduce repeated image-mount overhead.

```bash
singularity build \
    --sandbox ParaChrSNP_sandbox \
    ParaChrSNP.sif
```

- `build`: creates a new container from an existing image.
- `--sandbox`: requests an extracted directory rather than another SIF file.
- `ParaChrSNP_sandbox`: output sandbox directory.
- `ParaChrSNP.sif`: source image.

After conversion, update the configuration.

```yaml
container:
    image: "ParaChrSNP_sandbox"
```

## Verify the installation

Check the main programs inside the image.

```bash
apptainer exec ParaChrSNP.sif minibwa version
apptainer exec ParaChrSNP.sif samblaster --help
apptainer exec ParaChrSNP.sif samtools --version
apptainer exec ParaChrSNP.sif gatk --version
```

The current optimized container is expected to provide minibwa `0.6-r416`,
samblaster `0.1.26`, and SAMtools `1.21`.

## Install GLnexus support

GLnexus is optional and is intentionally not bundled in the main container.
When `params.joint_calling.method` is set to `glnexus`, place an executable
named `glnexus_cli` at:

```text
scripts/glnexus_cli
```

Then grant execution permission.

```bash
chmod +x scripts/glnexus_cli
scripts/glnexus_cli --help
```

The executable is deliberately ignored by Git because the distributed binary
is too large for normal GitHub storage.

## Conda-based development environment

The container is recommended for routine analysis. Developers who need to
inspect or test individual rules can create the supplied Conda environment.

```bash
conda env create -f environment.yaml
conda activate parachrsnp
```

This environment does not replace the container path configured for Snakemake
unless the workflow is run without `--use-singularity`.
