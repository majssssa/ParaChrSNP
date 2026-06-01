# Traditional CombineGVCFs benchmark workflow

This independent Snakemake workflow benchmarks a traditional whole-genome GVCF strategy from paired-end FASTQ files to `result_vcfs/combined.vcf.gz`.

Pipeline:

```text
bwa mem -> samtools sort -> GATK MarkDuplicates -> samtools index
-> per-sample whole-genome GATK HaplotypeCaller
-> CombineGVCFs
-> GenotypeGVCFs
```

Edit `config.yaml` before running. FASTQ files are expected as `{prefix}.1.fq.gz` and `{prefix}.2.fq.gz`.

```bash
snakemake --use-singularity --cores 60 --singularity-args "-B /path/to/data"
```

Benchmark files are written to `benchmarks/`.
