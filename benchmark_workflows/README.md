# ParaChrSNP lightweight benchmark workflows

This directory contains two minimal workflows for runtime comparison from read alignment to `combined.vcf.gz`.

- `parachrsnp_chrwise/`: ParaChrSNP-style benchmark using `bwa-mem2`, per-sample x chromosome `HaplotypeCaller`, chromosome-level `GenomicsDBImport`, chromosome-level `GenotypeGVCFs`, and `GatherVcfs`.
- `traditional_combinegvcfs/`: Traditional benchmark using `bwa mem`, per-sample whole-genome `HaplotypeCaller`, `CombineGVCFs`, and `GenotypeGVCFs`.

Edit each `config.yaml` before running on a new dataset. The final output of each workflow is:

```text
result_vcfs/combined.vcf.gz
```

Run the ParaChrSNP-style benchmark:

```bash
cd benchmark_workflows/parachrsnp_chrwise
snakemake --use-singularity --cores 60 --singularity-args "-B /home/majunpeng/sda2"
```

Run the traditional benchmark:

```bash
cd benchmark_workflows/traditional_combinegvcfs
snakemake --use-singularity --cores 60 --singularity-args "-B /home/majunpeng/sda2"
```

Benchmark files are written under each workflow's `benchmarks/` directory. To compare only the alignment-to-VCF workflow, build reference indexes once before timed runs, or exclude index-rule time from the final summary.
