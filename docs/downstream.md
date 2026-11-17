# 独立 VCF 下游分析入口

发布入口为 `Snakefile.downstream`，配套配置为
`config.downstream.example.yaml`。该入口不包含序列比对、BAM 处理、
变异检测或联合分型规则。原始全流程仍使用 `Snakefile`。

## 配置

在仓库根目录运行。修改示例配置中的 `input_vcf`、`reference` 和
`container.image`，使其指向实际文件。输入 VCF 不能与流程输出文件同路径。
样本名自动从 VCF 表头读取，不需要另设样本清单或输入模式。

默认运行缺失率、PLINK/HapMap 转换、PCA、距离矩阵和树、汇总报告。
其他下游模块使用各自的 `enabled` 开关。参考与注释必须匹配；
输入应为包含 GT 的普通群体 VCF，不能包含 gVCF 参考区块。
分析用 SNP 为 FILTER=PASS 或 `.` 的双等位 SNP。

先检查任务图。

```bash
snakemake --snakefile Snakefile.downstream --configfile config.downstream.example.yaml --cores 1 -n
# --snakefile: 指定独立下游入口。
# --configfile: 指定已修改的下游配置。
# --cores 1: 使用一个核心检查任务图。
# -n: 只检查任务依赖，不执行分析。
```

使用配置中的容器运行下游分析。

```bash
snakemake --snakefile Snakefile.downstream --configfile config.downstream.example.yaml --cores 8 --use-singularity
# --snakefile: 指定独立下游入口。
# --configfile: 指定已修改的下游配置。
# --cores 8: 设置总可用线程数上限为八。
# --use-singularity: 通过 Singularity/Apptainer 执行各条规则。
```

不要在已有其他数据结果的目录中直接复用输出。独立工作目录运行时，
需确保配置引用的 `scripts/` 路径可用，并把外部参考、输入和脚本目录
绑定到容器内的相同路径。容器受宿主机 Python/R 用户配置影响时，
可在 `--singularity-args` 中设置 `--cleanenv`，并通过
`--env PYTHONNOUSERSITE=1 --env R_ENVIRON_USER=/dev/null
--env R_PROFILE_USER=/dev/null --env R_LIBS_USER=/nonexistent`
隔离宿主机的用户包和启动配置。

从仓库根目录启动独立输出目录时，将仓库路径显式绑定进容器，
并指定容器可访问的 `/tmp` 临时目录。以下命令使用已准备好的本地配置，
且假定 `downstream_run/scripts` 已指向仓库的 `scripts/`。

```bash
env TMPDIR=/tmp snakemake \
  --snakefile Snakefile.downstream \
  --configfile config.downstream.local.yaml \
  --directory downstream_run \
  --cores 8 --use-singularity --rerun-incomplete \
  --singularity-args "--cleanenv --bind $PWD:$PWD --env TMPDIR=/tmp --env PYTHONNOUSERSITE=1 --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null --env R_LIBS_USER=/nonexistent"
# TMPDIR=/tmp: 为调度器和规则选择容器可访问的临时目录。
# --snakefile: 指定独立下游入口。
# --configfile: 指定本机已修改的配置文件。
# --directory: 将输出写入独立目录。
# --cores: 限制总可用线程数。
# --use-singularity: 使用配置中的容器执行规则。
# --rerun-incomplete: 重跑上次未完成的任务。
# --singularity-args: 向每次容器执行传递后续参数。
# --cleanenv: 清除继承的宿主机环境变量。
# --bind $PWD:$PWD: 将当前仓库目录绑定到容器中的相同位置。
# --env TMPDIR=/tmp: 显式指定容器中的临时目录。
# --env PYTHONNOUSERSITE=1: 不加载宿主机 Python 用户包。
# --env R_ENVIRON_USER=/dev/null: 不加载宿主机 R 用户环境文件。
# --env R_PROFILE_USER=/dev/null: 不加载宿主机 R 用户启动脚本。
# --env R_LIBS_USER=/nonexistent: 不加载宿主机 R 用户包。
```

如果输入、参考或注释位于仓库外，还需为这些外部目录添加 `--bind`。

## GitHub 文件依赖

这两个发布文件属于完整 ParaChrSNP 仓库，不能脱离 `rules/` 和 `scripts/`
单独执行。上传时需同步相关新文件，尤其是：

- `rules/vcf_entry_2.rules`
- `rules/vcf2dis_1.rules`
- `rules/report_vcf_1.rules`
- `scripts/parachrsnp-report_7.py`

其余下游规则和绘图脚本沿用仓库已有文件。无需上传输入 VCF、参考序列、
容器镜像或测试输出。本地历史版本不再作为推荐发布入口。
