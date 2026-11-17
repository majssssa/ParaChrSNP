# 已有 VCF 的独立下游分析（v0.1.1）

`Snakefile.downstream` 是直接接收已有群体 VCF 的独立工作流，不是原始
变异检测流程中的输入模式开关。配置不包含比对、BAM 处理或变异检测参数。
原始 `Snakefile` 仍用于 FASTQ/BAM 全流程。

## 输入与配置

{download}`下载完整下游配置模板 <../../config.downstream.example.yaml>`

输入支持 `.vcf` 或 `.vcf.gz`，要求样本名唯一并包含 FORMAT/GT 基因型，
同时提供匹配的参考 FASTA。样本名自动从表头读取，不需要填写 `Types`、
`samples` 或 `chromosomes`。预检查拒绝缺失表头、重复样本名、缺失 GT、
gVCF 参考区块以及染色体名称不在参考 FASTA 中的记录。
染色体名称检查不能证明参考序列本身完全一致。

默认执行缺失率统计、PLINK 二进制/文本与 HapMap 转换、PCA、遗传距离与
系统发育树、HTML/TSV 报告。PCA 和距离树要求至少三个样本。
发布模板默认关闭其他可选模块。可根据需要开启 `params.imputation`、
`pi`、`snp_density`、`admixture`、`ld_decay`、`snpeff` 的 `enabled`。
注释模块需要匹配的 FASTA/GFF；ADMIXTURE 的 K 范围和交叉验证参数应适合
实际样本数量。资源要求见[可选模块指南](optional_modules.md)。
仅有 VCF 时不能运行依赖 BAM 的 CNVnator。

## 隔离原有结果

以下准备命令从克隆后的仓库根目录执行。必须保留仓库中的 `rules/` 和
`scripts/`，不能只复制两个入口和配置文件运行。

创建新的工作目录并复制完整配置。

```bash
mkdir -p downstream_run
cp config.downstream.example.yaml downstream_run/config.yaml
ln -s "$PWD/scripts" downstream_run/scripts
# mkdir -p: 创建输出目录，目录已存在时不报错。
# cp: 将完整模板复制为本次运行配置。
# ln -s: 在新工作目录中链接仓库的绘图和处理脚本。
# $PWD: 当前仓库根目录的绝对路径；请在仓库根目录执行。
```

编辑 `downstream_run/config.yaml`，将顶层 `input_vcf`、`reference` 和
`container.image` 改成实际文件的**绝对路径**。相对路径会按工作目录解析。
下游输出路径保持相对路径，使结果写入 `downstream_run`。
如果开启群体分组、注释或填补，还要调整相关资源路径。

输入可以直接使用原项目中的 `result_vcfs/combined.snp.filtered.vcf.gz`，
因为新结果会写入独立工作目录。输入与流程输出不能指向同一文件；
流程会拒绝这种冲突，以保护原文件。

## 检查与执行

修改配置后先检查任务图。

```bash
env TMPDIR=/tmp snakemake --snakefile Snakefile.downstream \
  --configfile downstream_run/config.yaml --directory downstream_run \
  --cores 1 -n
# TMPDIR=/tmp: 使用可访问的临时目录。
# --snakefile: 指定独立下游入口。
# --configfile: 加载本次运行的完整配置。
# --directory: 指定独立工作目录和输出位置。
# --cores 1: 为检查分配一个核心。
# -n: 仅构建任务图，不执行规则；内容预检查在实际执行时进行。
```

显式绑定目录并隔离宿主机用户环境后，执行所有已开启的分析。

```bash
env TMPDIR=/tmp snakemake --snakefile Snakefile.downstream \
  --configfile downstream_run/config.yaml --directory downstream_run \
  --cores 8 --use-singularity --rerun-incomplete \
  --singularity-args "--cleanenv --bind $PWD:$PWD --env TMPDIR=/tmp --env PYTHONNOUSERSITE=1 --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null --env R_LIBS_USER=/nonexistent"
# --cores 8: 并发任务最多分配八个 CPU 线程。
# --use-singularity: 按配置为每条规则调用 Singularity/Apptainer 容器。
# --rerun-incomplete: 重跑中断或未完成任务。
# --singularity-args: 向容器传递后续参数。
# --cleanenv: 避免继承无关的宿主机环境变量。
# --bind $PWD:$PWD: 将仓库绑定到容器中的相同路径。
# --env TMPDIR=/tmp: 显式使用容器可访问的临时目录。
# --env PYTHONNOUSERSITE=1: 不加载宿主机 Python 用户包。
# --env R_ENVIRON_USER=/dev/null: 不加载宿主机 R 用户环境文件。
# --env R_PROFILE_USER=/dev/null: 不加载宿主机 R 用户启动脚本。
# --env R_LIBS_USER=/nonexistent: 不加载宿主机 R 用户包。
# --snakefile、--configfile、--directory: 含义与上面的检查命令相同。
```

如果输入、参考或注释在仓库之外，例如 `/data`，还需在
`--singularity-args` 中添加 `--bind /data:/data`。绝对路径并不会自动挂载
宿主机目录。faidx 报文件不存在时应检查绑定；`mkdtemp ... failed` 则应
检查临时目录是否在容器内可访问。

## 输出与恢复

原始 VCF 不会被修改，排序并 BGZF 压缩的副本写入
`downstream_run/result_vcfs/`。SNP 分析保留 FILTER 为 `PASS` 或 `.` 的
双等位 SNP，同时导出相同 FILTER 条件下的 INDEL。不重复执行 GATK 硬过滤。
只有 SNP 的输入可能产生空 INDEL 文件；多等位 SNP 被排除，不会自动拆分。

综合报告为 `downstream_run/reports/ParaChrSNP_report.html`，统计表为同目录
的 `ParaChrSNP_summary.tsv`。仅凭 VCF 无法得到读段质控和基于 reads 的测序
深度，报告会显示不可用。距离矩阵和树保留完整样本名，包括超过十个字符的名称。

失败后使用相同命令继续，已完成的输出会被复用。只删除最终报告文件再运行，
可单独重建报告而不重算已完成的分析。更换数据集时请使用新的工作目录。
