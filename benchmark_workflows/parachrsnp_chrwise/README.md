# ParaChrSNP chromosome-wise benchmark workflow

This independent Snakemake workflow benchmarks the ParaChrSNP-style strategy from paired-end FASTQ files to `result_vcfs/combined.vcf.gz`.

Pipeline:

```text
bwa-mem2 mem -> samtools sort -> GATK MarkDuplicates -> samtools index
-> per sample x chromosome GATK HaplotypeCaller
-> chromosome-level GenomicsDBImport
-> chromosome-level GenotypeGVCFs
-> GatherVcfs
```

Edit `config.yaml` before running. FASTQ files are expected as `{prefix}.1.fq.gz` and `{prefix}.2.fq.gz`.

```bash
snakemake --use-singularity --cores 60 --singularity-args "-B /path/to/data"
```

Benchmark files are written to `benchmarks/`.
