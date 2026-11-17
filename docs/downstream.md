# 已有 VCF 的下游分析

使用 `Snakefile.downstream`，可以直接从已有群体 VCF 开展下游分析，
无需重新进行比对和变异检测。

## 配置

[下载完整配置模板](../config.downstream.example.yaml)

编辑 `config.downstream.example.yaml`，将 `input_vcf`、`reference` 和
`container.image` 分别设置为输入 VCF、匹配的参考 FASTA 和容器文件路径。
以下命令均在仓库根目录执行。输入文件不能与 `result_vcfs/` 中的流程输出同路径。

输入支持 `.vcf` 或 `.vcf.gz`，需包含 GT 基因型。样本名自动从 VCF 表头读取，
无需额外填写样本清单。PCA 和距离树要求至少三个样本。

默认执行缺失率统计、PLINK/HapMap 格式转换、PCA、遗传距离与系统发育树、
综合报告。通过对应的 `params.*.enabled` 参数，可以开启 ADMIXTURE、LD 衰减、
核苷酸多样性、SNP 密度、基因型填补和功能注释。
详细参数见[可选模块指南](https://parachrsnp.readthedocs.io/zh-cn/latest/usage/optional_modules.html)。

## 运行

先检查配置和任务图。

```bash
snakemake --snakefile Snakefile.downstream --configfile config.downstream.example.yaml --cores 1 -n
# --snakefile: 指定独立下游工作流。
# --configfile: 加载已修改的完整配置。
# --cores 1: 使用一个核心进行检查。
# -n: 仅检查任务图，不执行分析。
```

执行已开启的下游分析。

```bash
snakemake --snakefile Snakefile.downstream --configfile config.downstream.example.yaml --cores 8 --use-singularity
# --snakefile: 指定独立下游工作流。
# --configfile: 加载已修改的完整配置。
# --cores 8: 最多使用八个 CPU 线程。
# --use-singularity: 使用配置中指定的容器。
```

## 结果

综合报告及图表位于 `reports/ParaChrSNP_report.html`，
汇总统计位于 `reports/ParaChrSNP_summary.tsv`。
各项分析结果保存在对应的输出目录中。

原始 VCF 不会被修改。SNP 分析使用 FILTER 为 `PASS` 或 `.` 的双等位 SNP，
多等位 SNP 会被排除。仅凭 VCF 无法获得读段质控及基于 reads 的测序深度。
