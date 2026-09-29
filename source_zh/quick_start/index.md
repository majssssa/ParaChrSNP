# 快速开始：拟南芥四样本

本页假定已按[安装指南](../installation/index.md)获取代码、可用的 `ParaChrSNP.sif`，并安装宿主机 Snakemake 与 Singularity/Apptainer。请在仓库根目录执行命令。此示例使用完整的 `config.test.yaml`，**不要从页面中截取少量 YAML 字段另存为配置**。

## 1. 下载并解压示例数据

下载拟南芥示例压缩包。

```bash
wget http://www.majunpeng.com/ParaChrSNP/example.tar.gz

# wget: 下载示例数据。
# URL: 示例数据包地址，下载文件名为 example.tar.gz。
```

解压下载的压缩包。

```bash
tar -xvf example.tar.gz

# tar: 处理归档文件。
# -x: 解压；-v: 显示文件名；-f example.tar.gz: 指定输入包。
```

## 2. 整理输入

创建流程预期的输入目录。

```bash
mkdir -p raw_fastq reference

# mkdir: 创建目录。
# -p: 目录已经存在时不报错，并按需创建父目录。
# raw_fastq、reference: 分别存放双端测序数据和参考基因组。
```

把示例参考基因组和测序文件移至对应目录。

```bash
mv example/Arabidopsis_thaliana* reference/
mv example/*.fq.gz raw_fastq/

# mv: 移动文件。
# example/Arabidopsis_thaliana*: 示例参考序列及其相关文件。
# example/*.fq.gz: 示例双端 FASTQ 文件。
# reference/、raw_fastq/: 与 config.test.yaml 中路径一致的目标目录。
```

四个样本须分别具有 `{sample}.1.fq.gz` 和 `{sample}.2.fq.gz`。例如 `ERR16804307.1.fq.gz` 与 `ERR16804307.2.fq.gz`。

## 3. 使用完整示例配置

仓库自带的 `config.test.yaml` 已包含所有必需参数、四个样本、五条染色体、`Types: reads`、`container.image: ParaChrSNP.sif` 和 `params.joint_calling.method: glnexus`。页面展示的文件与仓库文件相同，也可直接下载：

{download}`下载完整 config.test.yaml <../../config.test.yaml>`

```{literalinclude} ../../config.test.yaml
:language: yaml
:linenos:
```

默认 `config.yaml` 是茶树示例，**不能**用于本页拟南芥数据。若使用自己的物种，请以[完整配置模板](../usage/full_configuration.md)为起点，并将染色体名改为与 FASTA 标题完全一致的标识。

由于此配置选择 GLnexus，执行下一步前还要下载 `scripts/glnexus_cli`；具体命令见[安装指南](../installation/index.md)。

## 4. 运行输入预检查

在正式计算前验证输入文件、样本与参考信息。

```bash
snakemake --snakefile Snakefile --configfile config.test.yaml --cores 1 --use-singularity reports/precheck.done

# --snakefile Snakefile: 使用 ParaChrSNP 主工作流。
# --configfile config.test.yaml: 加载完整拟南芥示例配置。
# --cores 1: 为预检查分配一个核心。
# --use-singularity: 通过配置指定的容器执行规则。
# reports/precheck.done: 仅运行输入预检查目标。
```

若预检查失败，先查看 `reports/precheck.tsv`、`reports/precheck.html` 和 `logs/precheck/precheck.log`。

## 5. 检查任务图

使用 dry-run 检查依赖关系，不执行分析任务。

```bash
snakemake --snakefile Snakefile --configfile config.test.yaml --cores 64 --use-singularity -n

# --snakefile Snakefile: 指定主工作流。
# --configfile config.test.yaml: 使用本页完整示例配置。
# --cores 64: 告诉 Snakemake 可调度的 CPU 总数；可按机器资源调整。
# --use-singularity: 使用配置指定的镜像。
# -n: 只生成任务图，不运行分析命令。
```

## 6. 正式运行

预检查和 dry-run 均通过后，开始分析。

```bash
snakemake --snakefile Snakefile --configfile config.test.yaml --cores 64 --use-singularity --keep-going --rerun-incomplete

# --snakefile Snakefile: 指定流程入口。
# --configfile config.test.yaml: 继续使用相同的四样本配置。
# --cores 64: 同时调度任务的 CPU 上限，可按机器资源降低。
# --use-singularity: 在容器中运行工具。
# --keep-going: 某一分支失败时允许无依赖关系的任务继续。
# --rerun-incomplete: 重新执行之前中断且被标记为未完成的输出。
```

主要结果在 `result_vcfs/combined.vcf.gz`、`result_vcfs/combined.snp.filtered.vcf.gz`、`result_vcfs/combined.indel.filtered.vcf.gz` 和 `reports/ParaChrSNP_report.html`。相同命令再次运行时，Snakemake 会复用已完成结果。
