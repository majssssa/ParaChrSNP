# FAQ

## Why does Snakemake report a missing reference genome?

`reference` is resolved relative to the directory in which Snakemake is
started. Confirm the configured path and run from the ParaChrSNP project root.

```bash
pwd
ls -lh reference/your_genome.fasta
```

Also confirm that `--configfile` points to the intended configuration. The
reference named in `config.test.yaml` can differ from the one in `config.yaml`.

## Why does samblaster report “Missing header on input sam file”?

samblaster requires a valid SAM header from the aligner. This error usually
means the input SAM is empty, truncated, or was generated without header lines.
Prefer the integrated streaming rule:

```text
minibwa map | samblaster --removeDups | samtools sort
```

Inspect the aligner log before troubleshooting samblaster itself.

## Does samblaster have a thread option?

samblaster does not expose a normal multi-thread option. Parallelism is
provided by the aligner and `samtools sort` on either side of the pipe.

## Why are both `{sample}.0000.bam` and `{sample}.rmdup.bam` present?

`{sample}.0000.bam` is a temporary chunk created by `samtools sort -T`. It is
normally removed after a successful sort. It can remain when Snakemake or the
pipeline is interrupted. `{sample}.rmdup.bam` is the expected final output,
but it should not be trusted if Snakemake marked the job incomplete.

After confirming that no Snakemake or SAMtools process is active, remove stale
sort chunks and rerun incomplete jobs.

```bash
rm duplicate_removed/*.0000.bam
snakemake --use-singularity --cores 64 --rerun-incomplete
```

## Why does GATK MarkDuplicates reject an unsorted SAM file?

GATK MarkDuplicates requires coordinate-sorted or queryname-sorted input. The
optimized ParaChrSNP path no longer sends an unsorted SAM file to GATK; it uses
samblaster during mapping and then sorts the resulting alignments with
SAMtools.

## Why is `scripts/glnexus_cli` missing?

GLnexus is optional and its binary is not committed to Git or bundled in the
container. Download it separately, place it at `scripts/glnexus_cli`, and make
it executable.

```bash
chmod +x scripts/glnexus_cli
```

Alternatively, set `glnexus_executable` to an absolute path visible inside the
container environment.

## Why does PLINK stop on a GT half-call?

GLnexus can emit genotypes such as `./1` or `0/.`. PLINK 1.9 rejects these by
default. ParaChrSNP passes:

```text
--vcf-half-call missing
```

This treats the entire incomplete diploid genotype as missing instead of
inferring a genotype from one observed allele.

## Why does Snakemake report incomplete files?

Snakemake records outputs from interrupted or failed jobs as incomplete. Resume
with:

```bash
snakemake \
    --use-singularity \
    --configfile config.yaml \
    --cores 64 \
    --rerun-incomplete
```

Do not use `--cleanup-metadata` unless the output has been independently
validated and you intentionally want Snakemake to accept it.

## Why can the workflow directory not be locked?

First check whether another Snakemake process is using the same working
directory. If none is running, unlock it:

```bash
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --unlock
```

Never unlock a directory while another workflow process is still active.

## Why does a GATK filter expression contain an illegal newline?

The complete expression must be passed as one shell argument. Do not insert a
physical newline inside the quoted expression.

```bash
gatk VariantFiltration \
    -R reference.fasta \
    -V combined.snp.vcf.gz \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "SNP_filter" \
    -O combined.snp.marked.vcf.gz
```

## Does VariantFiltration remove failed sites?

No. `VariantFiltration` adds a value to the `FILTER` field. ParaChrSNP follows
it with `SelectVariants --exclude-filtered` to create a VCF containing only
passing sites.

## Why does a container work on the host but fail in Snakemake?

Check `container.image`, file permissions, external path bindings, and whether
the configured executable is visible inside the container.

```bash
apptainer exec ParaChrSNP.sif minibwa version
apptainer exec ParaChrSNP.sif samblaster --help
```

When data are outside the project directory, pass explicit bind mounts through
`--singularity-args`.

## Why does Read the Docs fail to import a Sphinx extension?

Every extension listed in `source/conf.py` must be installed through
`source/requirements.txt`, and the repository-root `.readthedocs.yaml` must
contain:

```yaml
python:
  install:
    - requirements: source/requirements.txt
```

The `.readthedocs.yaml` file must be at the repository root. A nested file in
`source/` is not used as the default project configuration.

## Where should I look after a failed rule?

Start with the latest file in `.snakemake/log/`, identify the first failed rule,
and then open that rule's dedicated log under `logs/`. The final Snakemake
message often reports only that one command returned a non-zero status; the
tool-specific log normally contains the actual cause.
