# 更新日志

这里记录当前仓库中对用户可见的重要变更。v0.1.1 之前的日期条目是开发里程碑。

## v0.1.1 - 2026-11-17

- 新增 `Snakefile.downstream` 和 `config.downstream.example.yaml`，直接从已有
  群体 VCF 执行独立下游分析，不包含上游参数。
- 增加记录级输入校验和输入/输出冲突保护；报告明确区分 VCF 分析结果与
  无法获得的读段质控和测序深度指标。
- 修复遗传距离矩阵和 Newick 树中长样本名被替换的问题。
- 增加中英文教程，覆盖独立工作目录、容器绑定、临时目录及 Python/R 隔离。
- 完成现有完整群体数据的下游测试、区域数据的 SnpEff/Beagle 测试，
  验证单独恢复报告不会重算上游结果。
- 原 FASTQ/BAM 流程和容器镜像保持不变。

## 2026-09-29

- 新增独立的中文 Sphinx/Read the Docs 文档源，与英文站点共用同一仓库。
- 快速开始页面直接展示并提供完整的拟南芥 `config.test.yaml`，避免把配置片段误用为独立文件。
- 拟南芥示例使用 GLnexus，并将镜像名称与安装页的 `ParaChrSNP.sif` 保持一致。

## 2026-09-28

### BAM 输入

- 新增 `Types: reads|bam`；默认仍是 reads。
- BAM 模式自动查找 `{bam_dir}/{sample}.bam`，核对输入后链接并建立索引；跳过 FASTQ 质控、fastp、比对和 samblaster。
- 报告中无法由 BAM 输入得到的 reads 质控指标显示为 `NA`。

### 报告计数

- HTML/TSV 报告直接扫描 VCF 记录数，避免索引没有计数元数据时把有效 VCF 错报为 0。
- 在四样本 BAM 模式实测数据上复核。

## 2026-08-18

- GLnexus 联合分型后删除重复的二次质量过滤，只将结果拆分为 SNP 与 INDEL。
- GenomicsDB 和 CombineGVCFs 仍保留 GATK 硬过滤；沿用既有 `combined.*.filtered.vcf.gz` 路径以兼容下游模块。

## 2026-08-03

- 添加 Read the Docs/Sphinx 文档、依赖文件及完整配置示例。
- 集成 minibwa `0.6-r416`、samblaster 流式去重复及其重复率统计。
- 增加 GLnexus 逐染色体联合分型，同时保留 GenomicsDB 和 CombineGVCFs。`glnexus_cli` 由用户另外下载到 `scripts/glnexus_cli`。
- 向读取 GLnexus VCF 的 PLINK 命令添加 `--vcf-half-call missing`，处理 `./1` 等半缺失基因型。
- 改进 `--configfile` 选择的配置在预检查和报告中的传递。

## 2026-06

- 在 BWA-MEM2 和 BWA 之外增加 minibwa 比对后端，调整比对、排序线程分配。
- 增加染色体级与传统联合分型布局的性能比较流程。

## 报告问题时请提供

```text
ParaChrSNP Git 提交号
容器文件名及 SHA256
Snakemake 版本
比对工具和版本
联合分型模式
执行命令
相关规则日志
```

查看当前代码版本。

```bash
git rev-parse --short HEAD

# git rev-parse: 解析 Git 版本标识。
# --short HEAD: 输出当前提交的短哈希。
```

计算容器校验值。

```bash
sha256sum ParaChrSNP.sif

# sha256sum: 计算文件 SHA256 校验和。
# ParaChrSNP.sif: 当前分析使用的镜像，应与配置一致。
```
