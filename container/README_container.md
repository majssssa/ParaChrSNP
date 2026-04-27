# ParaChrSNP Singularity Container

## Build

构建 ParaChrSNP 容器镜像。

```bash
singularity build ParaChrSNP.sif container/ParaChrSNP.def

# singularity build: 根据 definition 文件构建 Singularity/Apptainer 镜像。
# ParaChrSNP.sif: 输出的容器镜像文件。
# container/ParaChrSNP.def: 容器构建定义文件。
```

如果集群不允许普通用户直接构建，可以在支持 Singularity/Apptainer build 的服务器上构建，再把 `ParaChrSNP.sif` 拷贝到集群。

## Run

使用 Snakemake 原生 Singularity 模式执行 dry-run。

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 1 --use-singularity -n

# snakemake: 运行 Snakemake 工作流。
# --snakefile Snakefile: 指定流程入口文件。
# --configfile config.yaml: 指定用户编辑后的配置文件。
# --cores 1: 使用 1 个 CPU 核心。
# --use-singularity: 使用 Snakefile/config.yaml 中声明的容器运行各个 rule。
# -n: dry-run，只检查流程，不真正运行任务。
```

运行完整流程。

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 30 --use-singularity --keep-going

# --cores 30: 最多使用 30 个 CPU 核心。
# --use-singularity: 使用 ParaChrSNP.sif 执行各个 rule。
# --keep-going: 某个任务失败后，继续运行其他不依赖失败任务的任务。
```

## Container-Specific Config Notes

容器中已经把以下程序放入 `PATH`：

```text
snakemake
bwa
bwa-mem2
samtools
bcftools
gatk
picard
fastp
plink
run_pipeline.pl
rastqc
VCF2Dis
VCF2PCACluster
Plot2Deig
Plot3Deig
snpEff
convert
rsvg-convert
Rscript
```

容器还会设置：

```bash
export LD_LIBRARY_PATH=/opt/conda/envs/parachrsnp/lib:$LD_LIBRARY_PATH

# LD_LIBRARY_PATH: 让 VCF2PCACluster、VCF2Dis 等本地二进制程序能找到 conda 环境中的动态链接库。
```

项目配置文件中包含容器镜像路径：

```yaml
container:
    image: "ParaChrSNP.sif"
```

把 Figshare 下载得到的 `ParaChrSNP.sif` 放在项目根目录即可。如果镜像放在其他位置，修改 `config.yaml` 中的 `container.image`。

容器运行时建议使用容器内命令名：

```yaml
vcf2pca:
    executable: "VCF2PCACluster"

vcf2dis:
    executable: "VCF2Dis"
```
