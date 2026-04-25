# ParaChrSNP

`ParaChrSNP` 是一个按染色体并行进行 SNP calling 的 Snakemake 流程。流程从双端 FASTQ 开始，依次完成原始数据质控、`fastp` 清洗、`bwa mem | samtools sort` 比对、GATK 按染色体 calling、GVCF 合并、联合分型、SNP/INDEL 过滤、VCF 缺失率统计、常用格式转换、PCA 和遗传距离/系统发育树分析。

## Workflow

![ParaChrSNP workflow](figures/ParaChrSNP_workflow.png)

**Figure 1.** Overview of the ParaChrSNP workflow. ParaChrSNP is a portable Snakemake workflow packaged with Singularity/Apptainer for chromosome-wise variant discovery and downstream population genomic analysis. Paired-end FASTQ files and a reference genome are first processed through raw-read quality control, read trimming, reference indexing, alignment, duplicate removal and BAM indexing. GATK HaplotypeCaller is then executed in parallel across samples and chromosomes to generate per-chromosome GVCFs, followed by sample-level GVCF merging and cohort-level joint genotyping. The resulting VCF is split into SNP and INDEL datasets and filtered using configurable thresholds. Filtered SNPs are further used for genotype missingness assessment, PLINK/HapMap format conversion, and optional PCA, genetic distance estimation and phylogenetic tree construction.

## 项目结构

- `Snakefile`: Snakemake 主入口，读取 `config.yaml` 并加载 `rules/` 中的规则。
- `rules/`: 每个分析步骤对应一个 `.rules` 文件。
- `config.yaml`: 参考基因组、样本、染色体、容器镜像和流程参数。
- `cluster.yaml`: SLURM 集群资源配置。
- `container/`: 容器构建文件、容器运行说明和 Figshare 发布说明。
- `scripts/`: 辅助脚本，例如 PLINK missing 结果绘图脚本。
- `raw_fastq/`: 默认原始 FASTQ 输入目录。
- `reference/`: 默认参考基因组目录。
- `pop.info`: 可选的样本分组文件，格式为 `sample group`。如果不提供，PCA 和遗传距离分析不会传入分组参数。

## 推荐运行方式：Singularity 容器

ParaChrSNP 推荐使用 Singularity/Apptainer 容器运行。用户宿主机只需要安装：

```text
snakemake
singularity 或 apptainer
```

流程所需的软件已经打包在 `ParaChrSNP.sif` 中，包括：

```text
bwa
samtools
bcftools
gatk
fastp
plink
TASSEL
rastqc
VCF2Dis
VCF2PCACluster
Python/pandas/numpy/matplotlib
```

## 下载项目和容器

下载 ParaChrSNP 项目代码。

```bash
git clone https://github.com/your_name/ParaChrSNP.git
cd ParaChrSNP

# git clone: 从 GitHub 下载 ParaChrSNP 项目代码。
# cd ParaChrSNP: 进入项目目录。
```

从 Figshare 下载容器镜像。上传 Figshare 后，把下面的 `FIGSHARE_DIRECT_DOWNLOAD_URL` 替换为真实下载链接。

```bash
singularity pull ParaChrSNP.sif "FIGSHARE_DIRECT_DOWNLOAD_URL"

# singularity pull: 从远程地址下载 Singularity 镜像。
# ParaChrSNP.sif: 下载后保存到项目根目录的容器文件名。
# FIGSHARE_DIRECT_DOWNLOAD_URL: Figshare 的直接下载链接，通常形如 https://figshare.com/ndownloader/files/文件ID。
```

如果使用 Apptainer：

```bash
apptainer pull ParaChrSNP.sif "FIGSHARE_DIRECT_DOWNLOAD_URL"

# apptainer pull: 使用 Apptainer 下载容器镜像。
# ParaChrSNP.sif: 下载后的本地容器文件。
# FIGSHARE_DIRECT_DOWNLOAD_URL: Figshare 直接下载链接。
```

## 配置输入文件

默认 FASTQ 命名格式如下：

```text
raw_fastq/{sample}.1.fq.gz
raw_fastq/{sample}.2.fq.gz
```

运行前需要检查并修改 `config.yaml`：

- `container.image`: 容器镜像路径，默认是项目根目录下的 `ParaChrSNP.sif`。
- `reference`: 参考基因组 FASTA 路径。
- `samples`: 样本名和 FASTQ 文件前缀。
- `chromosomes`: 需要逐条染色体 calling 的染色体名称，必须和参考基因组 FASTA 中的序列 ID 一致。
- `params.snp_filter` 和 `params.indel_filter`: SNP/INDEL 过滤参数。
- `params.vcf2pca.enabled`: 是否在完整流程中运行 PCA 分析，`true` 表示运行，`false` 表示不运行。
- `params.vcf2dis.enabled`: 是否在完整流程中运行遗传距离和系统发育树分析，`true` 表示运行，`false` 表示不运行。
- `params.vcf2pca.sample_group`: 可选的 PCA 样本分组文件。留空或删除该参数时，不传入 `-InSampleGroup`。
- `params.vcf2dis.sample_group`: 可选的遗传距离样本分组文件。留空或删除该参数时，不传入 `-InSampleGroup`。

容器模式下，`config.yaml` 中的软件建议使用命令名，不要写本机绝对路径：

```yaml
vcf2pca:
    enabled: true
    executable: "VCF2PCACluster"

vcf2dis:
    enabled: true
    executable: "VCF2Dis"
```

## 运行流程

检查 DAG 和输入文件是否完整。

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 4 --use-singularity -n

# snakemake: 运行 Snakemake 工作流。
# --snakefile Snakefile: 指定流程入口文件。
# --configfile config.yaml: 指定流程配置文件。
# --cores 4: 允许 Snakemake 使用 4 个 CPU 核心。
# --use-singularity: 使用 Snakefile/config.yaml 中声明的 Singularity 容器执行各个 rule。
# -n: dry-run，只检查流程，不真正运行任务。
```

运行完整流程。

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 30 --use-singularity --keep-going

# --cores 30: 最多使用 30 个 CPU 核心。
# --use-singularity: 在 ParaChrSNP.sif 容器中运行各个 rule。
# --keep-going: 某个任务失败后，继续运行其他不依赖失败任务的作业。
```

如果 FASTQ 或参考基因组在项目目录外，需要额外挂载外部目录。

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 30 --use-singularity --singularity-args "-B /data/fastq:/data/fastq -B /data/reference:/data/reference"

# --singularity-args: 传递额外的 Singularity 参数。
# -B /data/fastq:/data/fastq: 将宿主机 FASTQ 目录挂载到容器内相同路径。
# -B /data/reference:/data/reference: 将宿主机参考基因组目录挂载到容器内相同路径。
```

如果项目目录中的输入文件是软链接，并且软链接指向项目目录外的位置，也必须挂载软链接的真实目标目录。例如：

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 12 --use-singularity --singularity-args "-B /home/majunpeng/sda2:/home/majunpeng/sda2"

# --cores 12: 最多使用 12 个 CPU 核心。
# --use-singularity: 使用 ParaChrSNP.sif 容器运行每个 rule。
# --singularity-args: 传递额外的 Singularity 挂载参数。
# -B /home/majunpeng/sda2:/home/majunpeng/sda2: 将宿主机的 /home/majunpeng/sda2 挂载到容器内相同路径，使 reference/ 和 raw_fastq/ 中指向该目录的软链接在容器内也能正常访问。
```

解除 Snakemake 异常退出后留下的锁。

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 1 --use-singularity --unlock

# --unlock: 删除 Snakemake 工作目录中的锁文件。
```

## 输出结果

主要结果文件包括：

- `qc/`: RastQC 原始数据质控结果。
- `clean_reads/`: `fastp` 清洗后的 FASTQ。
- `sorted_reads/`: 排序后的 BAM。
- `duplicate_removed/`: 去重复后的 BAM。
- `gvcf/`: 每个样本、每条染色体的 GVCF，以及样本级合并 GVCF。
- `result_vcfs/combined.vcf.gz`: 联合分型后的原始 VCF。
- `result_vcfs/combined.snp.filtered.vcf.gz`: 过滤后的 SNP 结果。
- `result_vcfs/combined.indel.filtered.vcf.gz`: 过滤后的 INDEL 结果。
- `missing/`: PLINK 缺失率统计结果和缺失率分布图。
- `format_convert/`: PLINK binary、PLINK text 和 HapMap 格式转换结果。
- `pca/ParaChrSNP.eigenvec`: VCF2PCACluster 输出的 PCA 坐标。
- `pca/ParaChrSNP.eigenval`: VCF2PCACluster 输出的 PCA 特征值。
- `dis/ParaChrSNP.p_dis.mat`: VCF2Dis 输出的样本遗传距离矩阵。
- `dis/ParaChrSNP.p_dis.nwk`: VCF2Dis 输出的 Newick 格式系统发育树。

## 发布说明

`ParaChrSNP.sif` 体积较大，不建议直接提交到 GitHub 仓库。推荐把容器上传到 Figshare、Zenodo 或实验室服务器，并在 README 中提供直接下载链接。

容器校验值可以用下面的命令生成：

```bash
sha256sum ParaChrSNP.sif

# sha256sum: 计算容器文件的 SHA256 校验值，用于确认下载文件是否完整。
# ParaChrSNP.sif: 需要校验的容器镜像。
```
