# 使用指南

ParaChrSNP 根据配置中的样本、染色体、输入类型和模块开关生成 Snakemake 任务图。**必须提供完整配置文件**；只复制某个 `params` 片段会造成 `KeyError`。拟南芥示例使用 `config.test.yaml`，其他项目可从[完整模板](full_configuration.md)复制后修改。

## 输入类型与文件命名

在配置顶部选择 `Types: reads` 或 `Types: bam`。同一次运行不能混用两种输入；切换模式或替换样本时，建议使用新的工作目录，避免旧 GVCF 和报告被误认为新结果。

### 双端 FASTQ

`Types: reads` 时，`samples` 的值是去掉 `.1.fq.gz`/`.2.fq.gz` 的前缀。流程依次执行 fastp、minibwa/BWA 比对、samblaster 去重复、SAMtools 排序和索引。

```yaml
Types: reads
samples:
  sample1: "raw_fastq/sample1"
  sample2: "raw_fastq/sample2"
```

上例要求 `raw_fastq/sample1.1.fq.gz`、`raw_fastq/sample1.2.fq.gz` 等文件存在。参考基因组通过 `reference` 指定；`chromosomes` 必须与 FASTA 标题中的序列 ID 完全相同，大小写也要一致。

### 已去重复 BAM

`Types: bam` 时，每个样本文件名由键名决定，即 `{bam_dir}/{sample}.bam`，`samples` 的值可以是 `null`。

```yaml
Types: bam
bam_dir: input_bam
samples:
  sample1: null
  sample2: null
```

输入 BAM 必须已经去重复、按坐标排序，并具有 `@RG` 记录；其中 `SM` 要与样本键相同。BAM 参考序列名须与配置中的 FASTA 一致。预检查会核对可读性、排序声明、样本名和 contig，但不能证明去重复已完成，也不能证明 BAM 与 FASTA 序列完全一致。输入 BAM 不强制自带索引：流程链接到 `staged_bam/{sample}.bam` 并在该目录建立 `.bam.bai`，不改写源 BAM。超过 512 MiB 的 contig 无法使用 BAI，本流程当前不支持 CSI。

BAM 模式跳过 FASTQ 质控、fastp、比对和 samblaster；报告中不存在的 reads 指标显示为 `NA`。外部 BAM 目录需要通过容器 bind 挂载。

## 比对与染色体级检测

`params.aligner.name` 可选 `minibwa`、`bwa-mem2` 或 `bwa`。reads 模式的流式处理顺序为：

```text
minibwa map → samblaster --removeDups → samtools sort
```

主要 BAM 输出为 `duplicate_removed/{sample}.rmdup.bam`、其 `.bai` 索引以及 `{sample}.dup.txt`。GATK HaplotypeCaller 对每个样本、每条配置染色体单独生成 `gvcf/{sample}.{chrom}.g.vcf.gz`，便于并行计算。

## 联合分型模式

由 `params.joint_calling.method` 选择模式：

| 模式 | 工作方式 | 附加要求 |
| --- | --- | --- |
| `glnexus` | 逐染色体使用 GLnexus，随后按配置顺序合并 | `scripts/glnexus_cli` 可执行 |
| `genomicsdb` | 导入逐染色体 GVCF 到 GenomicsDB，再运行 GenotypeGVCFs | 容器中的 GATK |
| `combine_gvcfs` | 兼容旧流程，适合小规模对照 | 容器中的 GATK |

拟南芥 `config.test.yaml` 和项目主配置显式选择 `glnexus`；若从自建配置中省略此字段，Snakefile 当前的回退值仍是 `genomicsdb`，因此建议始终明确写出模式。

```yaml
params:
  joint_calling:
    method: "glnexus"
    glnexus_executable: "scripts/glnexus_cli"
    glnexus_config: "gatk"
```

最终联合 VCF 为 `result_vcfs/combined.vcf.gz`。GLnexus 使用其内置过滤逻辑，后续**不再**额外进行 GATK/bcftools 质量过滤，只拆分 SNP 和 INDEL。为兼容下游路径，拆分文件仍名为 `combined.snp.filtered.vcf.gz` 和 `combined.indel.filtered.vcf.gz`，文件名中的 `filtered` 不代表 GLnexus 之后又过滤了一次。

GenomicsDB 与 CombineGVCFs 则继续使用 GATK 硬过滤。SNP 条件为 `QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0`；INDEL 条件为 `QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0`。`VariantFiltration` 标记失败位点，随后 `SelectVariants --exclude-filtered` 将其排除。这些是位点级条件，不等于按 MAF、样本缺失率、基因型 DP/GQ 过滤。

## 可选分析

各模块由自己的 `enabled` 开关控制；详细原理、参数与结果见[可选模块指南](optional_modules.md)。

| 模块 | 开关 | 主要用途 |
| --- | --- | --- |
| Beagle | `params.imputation.enabled` | 推断缺失基因型 |
| CNVnator | `params.cnv.enabled` | 基于测序深度检测 CNV |
| Pi | `params.pi.enabled` | 窗口核苷酸多样性 |
| SNP density | `params.snp_density.enabled` | 统计染色体窗口 SNP 密度 |
| ADMIXTURE | `params.admixture.enabled` | 推断群体祖源成分 |
| PopLDdecay | `params.ld_decay.enabled` | 分析连锁不平衡衰减 |
| PCA | `params.vcf2pca.enabled` | 主成分分析 |
| VCF2Dis | `params.vcf2dis.enabled` | 遗传距离矩阵及树 |
| SnpEff | `params.snpeff.enabled` | SNP/INDEL 功能注释 |

涉及群体分组的模块可使用两列 `pop.info`：第一列是 VCF 样本名，第二列是群体名；样本名必须完全匹配。

```text
sample1    Population_A
sample2    Population_A
sample3    Population_B
```

```{toctree}
:maxdepth: 1
:caption: 详细使用说明

optional_modules
full_configuration
```

## 运行与检查

先只运行预检查，避免在大规模计算后才发现输入错误。

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 1 --use-singularity reports/precheck.done

# --snakefile Snakefile: 指定流程入口。
# --configfile config.yaml: 使用已按项目修改的完整配置文件。
# --cores 1: 给预检查分配一个核心。
# --use-singularity: 在配置指定的容器内运行规则。
# reports/precheck.done: 只请求输入预检查目标。
```

确认任务图，再启动全部目标。

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 64 --use-singularity -n
snakemake --snakefile Snakefile --configfile config.yaml --cores 64 --use-singularity --rerun-incomplete

# --cores 64: 全局最多分配 64 个 CPU 核心，可按机器配置调整。
# -n: 仅生成任务图，不执行规则。
# --rerun-incomplete: 重新运行上次中断且标记为未完成的任务。
```

数据在项目目录以外时，通过 `--singularity-args "-B 主机目录:容器目录"` 挂载；配置中的路径必须在容器内可见。只有确认没有其他 Snakemake 进程运行后，才可执行 `snakemake --snakefile Snakefile --configfile config.yaml --unlock` 清理陈旧锁。

## 结果目录

| 目录 | 内容 |
| --- | --- |
| `qc/`、`clean_reads/` | 原始 reads 质控和 fastp 清洗结果 |
| `duplicate_removed/`、`staged_bam/` | reads 模式去重复 BAM 或 BAM 模式链接与索引 |
| `gvcf/`、`genomicsdb/` | 单样本 GVCF 和可选的 GenomicsDB 数据库 |
| `result_vcfs/` | 联合、SNP、INDEL VCF |
| `missing/`、`format_convert/` | 缺失率统计、PLINK 与 HapMap 格式 |
| `imputation/`、`cnv/`、`pi/`、`snp_density/` | 对应可选分析结果 |
| `pca/`、`dis/`、`admixture/`、`ld_decay/`、`annotation/` | 群体分析、系统关系和注释 |
| `reports/` | 预检查和综合 HTML/TSV 报告 |

最终报告是 `reports/ParaChrSNP_report.html` 和 `reports/ParaChrSNP_summary.tsv`。报告中的变异数直接扫描 VCF 记录，不依赖可能缺失的索引计数元数据。浏览器启动与监控详见[Web 界面指南](../web_interface/index.md)。
