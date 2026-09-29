# 安装

ParaChrSNP 面向 64 位 Linux，推荐使用 Singularity/Apptainer 容器运行。开始前需准备 Git、宿主机上的 Snakemake、Singularity 或 Apptainer，以及容纳 FASTQ、BAM、GVCF、VCF 和临时文件的磁盘空间。实际 CPU、内存需求取决于基因组大小、深度和样本数。

## 获取流程代码

下载仓库并进入根目录；后续命令均在此目录运行。

```bash
git clone https://github.com/majssssa/ParaChrSNP.git
cd ParaChrSNP

# git clone: 从 GitHub 获取流程源代码。
# 仓库地址: ParaChrSNP 官方代码仓库。
# cd ParaChrSNP: 进入含 Snakefile 和配置文件的工作目录。
```

## 获取容器

下载安装页发布的镜像，并保存为配置文件要求的 `ParaChrSNP.sif`。

```bash
singularity pull ParaChrSNP.sif http://www.majunpeng.com/ParaChrSNP/ParaChrSNP.sif

# singularity pull: 下载 Singularity 镜像；也可使用 apptainer pull。
# ParaChrSNP.sif: 本地镜像文件名，必须与 container.image 一致。
# URL: 项目提供的容器下载地址。
```

如果本机已有同名旧镜像，不能只凭文件名判断是否为新版。下载后检查关键程序；预期新版 minibwa 为 `0.6-r416`，并包含 samblaster。

```bash
singularity exec ParaChrSNP.sif minibwa version
singularity exec ParaChrSNP.sif samblaster --help
singularity exec ParaChrSNP.sif samtools --version

# singularity exec: 在镜像内执行后面的程序。
# ParaChrSNP.sif: 要检查的镜像路径。
# minibwa version: 显示比对软件版本，应为 0.6-r416。
# samblaster --help: 验证去重复程序是否存在并可启动。
# samtools --version: 查看 SAMtools 版本。
```

镜像不在 Git 仓库中。如果镜像放在其他目录，修改配置中的 `container.image` 为容器内外都可见的路径。

## 安装 GLnexus

选择 `params.joint_calling.method: glnexus` 时，需另行下载 `glnexus_cli`；它不包含在主镜像和 Git 仓库中。

将官方 v1.4.1 可执行文件保存到流程默认位置。

```bash
wget -O scripts/glnexus_cli https://github.com/dnanexus-rnd/GLnexus/releases/download/v1.4.1/glnexus_cli

# wget: 从指定网址下载文件。
# -O scripts/glnexus_cli: 将文件写到配置中 glnexus_executable 指定的位置。
# URL: GLnexus v1.4.1 的发布文件地址。
```

添加执行权限并检查程序是否可启动。

```bash
chmod +x scripts/glnexus_cli
scripts/glnexus_cli --help

# chmod +x: 给文件增加可执行权限。
# scripts/glnexus_cli: 刚下载的程序；--help 用于显示帮助信息。
```

若改用 `genomicsdb` 或 `combine_gvcfs`，可以不安装 GLnexus。

## 可选：Conda 开发环境

需要单独检查规则或脚本时，可安装仓库中的 Conda 环境；正常分析仍推荐容器。

```bash
conda env create -f environment.yaml
conda activate parachrsnp

# conda env create: 创建 Conda 环境。
# -f environment.yaml: 从仓库提供的依赖文件读取软件清单。
# conda activate parachrsnp: 激活新建环境。
```

如果数据目录在项目根目录之外，运行 Snakemake 时还需通过 `--singularity-args` 将外部路径挂载进容器，详见[使用指南](../usage/index.md)。
