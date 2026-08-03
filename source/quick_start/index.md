# Quick start

This page runs the bundled Arabidopsis example through input checking and the
complete ParaChrSNP workflow.

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

## 3. Edit the minimal configuration

At minimum, verify the reference, container, samples, chromosomes, aligner, and
joint-calling method.

```yaml
reference: "reference/Arabidopsis_thaliana.fasta"

container:
    image: "ParaChrSNP.sif"

samples:
    ERR16804307: "raw_fastq/ERR16804307"
    ERR16805220: "raw_fastq/ERR16805220"

chromosomes:
  - NC_003070.9
  - NC_003071.7
  - NC_003074.8
  - NC_003075.7
  - NC_003076.8

params:
    aligner:
        name: "minibwa"
        executable: "minibwa"
        map_threads: 4
        sort_threads: 2
        index_threads: 4
        index_extra: ""
        map_extra: ""

    joint_calling:
        method: "genomicsdb"
```

Chromosome names must exactly match the sequence identifiers in the FASTA
header and index.

## 4. Run the precheck

Validate the configuration and input files before starting expensive jobs.

```bash
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --cores 1 \
    --use-singularity \
    reports/precheck.done
```

- `--snakefile Snakefile`: selects the workflow entry point.
- `--configfile config.yaml`: selects the analysis configuration.
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
    --configfile config.yaml \
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
    --configfile config.yaml \
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
