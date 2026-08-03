# Complete configuration example

The complete example below contains the core workflow settings and all nine
optional analysis modules. Every optional module is set to `enabled: true` to
show a full run. For a real project, enable only the analyses that are needed
and whose required inputs and executables are available.

{download}`Download config.full.example.yaml <../_static/config.full.example.yaml>`

Before using the file, replace:

- `reference/genome.fa` with the actual reference FASTA;
- the sample names and FASTQ prefixes;
- `Chr01`, `Chr02`, and `Chr03` with exact FASTA sequence identifiers;
- `pop.info` with a valid sample-population file, or use an empty string where
  the module permits it;
- `reference/genes.gff3` with an annotation from the same reference assembly;
- executable paths and CPU/memory values for the local environment.

The example has three samples because PCA and VCF2Dis require at least three.
Beagle, ADMIXTURE, and population-genomic statistics generally need more
samples for biologically reliable inference even when software can run on a
small test dataset.

```{literalinclude} ../_static/config.full.example.yaml
:language: yaml
:linenos:
```

## Validate the edited configuration

Run a dry-run to validate paths, rule dependencies, and the generated DAG
without executing analysis jobs.

```bash
snakemake \
    --snakefile Snakefile \
    --configfile config.full.yaml \
    --cores 1 \
    --use-singularity \
    --dry-run

# --snakefile Snakefile: selects the ParaChrSNP workflow entry point.
# --configfile config.full.yaml: loads the edited full configuration file.
# --cores 1: allocates one logical core for DAG validation.
# --use-singularity: applies the configured Singularity/Apptainer container.
# --dry-run: constructs and checks the workflow without running jobs.
```

Missing FASTQ, reference, annotation, population, Beagle, or GLnexus files must
be corrected before a full run. A successful YAML parse alone does not prove
that biological inputs use matching sample and chromosome identifiers.

