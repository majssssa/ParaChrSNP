# 完整配置示例

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
