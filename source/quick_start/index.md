# Quick start

This page runs the Arabidopsis example through input checking and the complete
ParaChrSNP workflow. Start in the root of a cloned ParaChrSNP repository, with
Snakemake and Singularity/Apptainer installed. Download the current
`ParaChrSNP.sif` image as described in [Installation](../installation/index.md)
and place it in the repository root; it is not included in the Git repository.
An older local file with the same name may lack samblaster, so verify the
program versions using the Installation instructions before running. The
GLnexus executable is downloaded in step 3 below.

## 1. Download and extract the example

Download the example archive.

```bash
wget http://www.majunpeng.com/ParaChrSNP/example.tar.gz
```

- `wget`: downloads a file from the supplied URL.
- `example.tar.gz`: compressed example dataset.

Extract the archive.

```bash
tar -xvf example.tar.gz
```

- `-x`: extracts files.
- `-v`: prints extracted file names.
- `-f`: reads the following archive name.

## 2. Prepare the input directories

Create the expected directory structure.

```bash
mkdir -p raw_fastq reference
```

Move the example reference and reads.

```bash
mv example/Arabidopsis_thaliana* reference/
mv example/*.fq.gz raw_fastq/
```

Paired reads must use the following names:

```text
raw_fastq/{sample}.1.fq.gz
raw_fastq/{sample}.2.fq.gz
```

For example, `ERR16804307.1.fq.gz` and `ERR16804307.2.fq.gz` are recognized as
sample `ERR16804307`.

## 3. Use the complete example configuration

Use the repository's `config.test.yaml` without copying a shortened YAML
excerpt into a new file. It contains every required workflow parameter, the
four Arabidopsis samples, the optimized container path, and GLnexus as the
joint-calling method. The complete, runnable file is shown below and can also
be downloaded directly:

{download}`Download config.test.yaml <../../config.test.yaml>`

```{literalinclude} ../../config.test.yaml
:language: yaml
:linenos:
```

The regular `config.yaml` contains tea-genome example values and must not be
used for this Arabidopsis quick start. For another species, edit a complete
configuration file such as `config.yaml` or the [full template](../usage/index.md).
Chromosome names must exactly match the reference FASTA headers.

Download GLnexus v1.4.1 to the path configured in `config.test.yaml` and make
it executable.

```bash
wget -O scripts/glnexus_cli https://github.com/dnanexus-rnd/GLnexus/releases/download/v1.4.1/glnexus_cli
chmod +x scripts/glnexus_cli

# wget: Download the official GLnexus v1.4.1 release binary.
# -O scripts/glnexus_cli: Save it at the executable path in config.test.yaml.
# chmod +x: Allow the downloaded binary to run.
```

## 4. Run the precheck

Validate the configuration and input files before starting expensive jobs.

```bash
snakemake \
    --snakefile Snakefile \
    --configfile config.test.yaml \
    --cores 1 \
    --use-singularity \
    reports/precheck.done
```

- `--snakefile Snakefile`: selects the workflow entry point.
- `--configfile config.test.yaml`: selects the complete Arabidopsis example configuration.
- `--cores 1`: allocates one core to the precheck target.
- `--use-singularity`: runs containerized rules with `container.image`.
- `reports/precheck.done`: requests only the precheck target.

Inspect these files if validation fails:

```text
reports/precheck.tsv
reports/precheck.html
logs/precheck/precheck.log
```

## 5. Perform a dry-run

Build and inspect the complete job graph without executing it.

```bash
snakemake \
    --snakefile Snakefile \
    --configfile config.test.yaml \
    --cores 64 \
    --use-singularity \
    --dry-run
```

The dry-run should finish without `MissingInputException`, configuration
errors, or ambiguous rule errors.

## 6. Run the workflow

Start the complete analysis.

```bash
snakemake \
    --snakefile Snakefile \
    --configfile config.test.yaml \
    --cores 64 \
    --use-singularity \
    --keep-going \
    --rerun-incomplete
```

- `--cores 64`: allows Snakemake to schedule up to 64 CPU cores.
- `--keep-going`: continues independent jobs after one branch fails.
- `--rerun-incomplete`: regenerates outputs interrupted during an earlier run.

## 7. Check the main results

The central outputs are:

```text
result_vcfs/combined.vcf.gz
result_vcfs/combined.snp.filtered.vcf.gz
result_vcfs/combined.indel.filtered.vcf.gz
reports/ParaChrSNP_report.html
reports/ParaChrSNP_summary.tsv
```

Snakemake can be run again with the same command. Completed outputs are reused,
and only missing, outdated, or incomplete jobs are scheduled.
