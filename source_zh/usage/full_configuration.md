# 完整配置示例

## 从已有 VCF 开始

使用新版入口 `Snakefile.downstream` 和完整配置模板 `config.downstream.example.yaml`，设置 `input_vcf`、匹配的参考基因组 `reference` 及容器路径。此入口和配置仅包含 VCF 整理与下游分析，无需指定输入模式。样本从 VCF 表头自动读取，无需填写 FASTQ/BAM 样本清单。建议在独立分析目录运行，避免复用其他数据的历史输出。

先检查任务依赖，再执行下游分析。

```bash
snakemake --snakefile Snakefile.downstream --configfile config.downstream.example.yaml --cores 1 -n
snakemake --snakefile Snakefile.downstream --configfile config.downstream.example.yaml --cores 8 --use-singularity
# --snakefile Snakefile.downstream: 指定支持 VCF 输入的新版工作流。
# --configfile config.downstream.example.yaml: 指定 VCF、参考及下游分析配置。
# --cores: 指定可用的总线程数。
# -n: 仅检查任务图，不执行分析。
# --use-singularity: 使用配置中指定的容器运行工具。
```

支持普通 VCF 和 gzip VCF；输入应包含 GT 基因型，不应包含 gVCF 参考区块。流程保留原文件，对副本排序、压缩和建立索引；SNP 下游使用 FILTER 为 PASS 或 `.` 的双等位 SNP，同时导出 PASS/`.` 的 INDEL。不再次执行 GATK 硬过滤。默认启用 PCA 和距离树（至少三个样本）；其他模块通过原有 `enabled` 参数开启。SnpEff 需要匹配的注释资源；此独立入口不包含依赖 BAM 的 CNV 检测模块。参考 FASTA 用于报告、π 和 SNP 密度等分析。报告中的读段质控与测序深度显示为不可用，不从旧 FASTQ 结果推算。

以下模板包含主流程和全部 9 个可选模块。模板为了展示设置项，将可选模块设为 `enabled: true`；**不能未经修改就直接运行**。请先替换参考 FASTA、样本前缀、染色体标识、群体文件、同一组装版本的 GFF3/GTF、容器路径，以及 CPU/内存参数。仅保留有真实输入且确实需要的模块。

{download}`下载 config.full.example.yaml <../../source/_static/config.full.example.yaml>`

```{literalinclude} ../../source/_static/config.full.example.yaml
:language: yaml
:linenos:
```

三样本足以演示 PCA 和 VCF2Dis 的调用，但通常不足以保证 Beagle、ADMIXTURE 等结果具有可靠的生物学解释。若只想验证拟南芥流程，请使用[快速开始](../quick_start/index.md)的 `config.test.yaml`，不要使用此“全部模块开启”的模板。

## 检查修改后的配置

把模板另存为项目配置并修改实际路径后，先进行 dry-run。

```bash
snakemake --snakefile Snakefile --configfile config.full.yaml --cores 1 --use-singularity -n

# --snakefile Snakefile: 指定 ParaChrSNP 工作流。
# --configfile config.full.yaml: 加载已按实际数据修改的完整配置。
# --cores 1: 使用一个逻辑核心检查任务图。
# --use-singularity: 按 container.image 使用容器。
# -n: 不执行分析，只验证输入和依赖关系。
```

YAML 能解析并不代表参考、注释、样本名和染色体名在生物学上匹配；正式运行前还应执行 `reports/precheck.done` 并查看报告。
