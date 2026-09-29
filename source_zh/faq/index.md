# 常见问题

## Snakemake 报缺少参考基因组？

`reference` 的相对路径相对于启动 Snakemake 的工作目录解析。请在项目根目录运行，并确认 `--configfile` 指向正确配置：`config.test.yaml` 是拟南芥示例，`config.yaml` 可能是其他物种。

检查当前位置和参考文件。

```bash
pwd
ls -lh reference/your_genome.fasta

# pwd: 显示当前工作目录。
# ls -lh: 以易读大小列出目标文件，确认路径和文件非空。
# reference/your_genome.fasta: 应替换为配置中真正的参考路径。
```

## 为什么复制快速开始的少量 YAML 会出现 `KeyError: genotype_gvcfs`？

Snakefile 会读取多个必需的 `params`，只包含 `aligner` 和 `joint_calling` 的片段不是可运行配置。拟南芥请直接使用完整的 `config.test.yaml`；自定义项目请下载[完整模板](../usage/full_configuration.md)后修改，不能逐项补一个键就认为已经完整。

## samblaster 报 `Missing header on input sam file`？

samblaster 的输入必须是含 `@HD/@SQ` 等头信息的有效 SAM。此错误多由上游比对失败、SAM 为空或截断引起。先查看比对日志；集成规则采用 `minibwa map | samblaster --removeDups | samtools sort` 的流式连接。

## samblaster 能设置线程数吗？

samblaster 没有常规多线程选项；流程通过上游比对工具及下游 `samtools sort` 并行。

## 为什么同时看到 `{sample}.0000.bam` 和 `{sample}.rmdup.bam`？

前者通常是 `samtools sort -T` 的临时分块，正常完成后会清除；中断时可能残留。后者是最终目标，但若 Snakemake 标记任务未完成，不能直接信任。确认所有相关进程已退出后，检查具体临时文件，再让 Snakemake 用 `--rerun-incomplete` 重跑。不要在任务仍运行时删除分块。

## GATK MarkDuplicates 为什么拒绝未排序 SAM？

该工具要求坐标或读名排序的输入。当前优化流程不再把未排序 SAM 交给 GATK MarkDuplicates，而是在比对输出流中使用 samblaster，然后使用 SAMtools 排序。

## 为什么找不到 `scripts/glnexus_cli`？

GLnexus 可执行文件不在 Git 仓库或主容器中；选择 `method: glnexus` 时须按[安装指南](../installation/index.md)下载并授权，也可在 `glnexus_executable` 中指定容器内可见的绝对路径。

## PLINK 报 GT half-call？

GLnexus 可能产生 `./1` 或 `0/.`。PLINK 1.9 默认拒绝这种不完整二倍体基因型；流程向 PLINK 传入 `--vcf-half-call missing`，将整个基因型当作缺失，不凭单个等位基因推断另一个。

## 为什么有 incomplete files？

Snakemake 会把中断或失败任务的输出标记为未完成。使用原配置重新执行，并添加 `--rerun-incomplete`。

```bash
snakemake --snakefile Snakefile --configfile config.yaml --cores 64 --use-singularity --rerun-incomplete

# --snakefile Snakefile: 使用项目主流程。
# --configfile config.yaml: 必须与原任务使用的配置相同。
# --cores 64: 最大调度核心数，可调整。
# --use-singularity: 通过容器执行。
# --rerun-incomplete: 重跑未完成输出。
```

未经独立验证，不应通过 `--cleanup-metadata` 强行接受可疑输出。

## 工作目录无法加锁？

先确认没有其他 Snakemake 进程在使用这个目录。若确定是陈旧锁，再执行解锁。

```bash
snakemake --snakefile Snakefile --configfile config.yaml --unlock

# --snakefile Snakefile: 指定工作流入口。
# --configfile config.yaml: 指向当前任务配置。
# --unlock: 只清理陈旧锁，不运行分析。
```

## GATK 过滤表达式为什么报非法换行？

整条表达式必须作为**一个** shell 参数传入，不能在引号内部插入物理换行。`VariantFiltration` 只在 `FILTER` 字段标记失败位点；流程随后由 `SelectVariants --exclude-filtered` 排除。GLnexus 模式不使用这一步 GATK 硬过滤。

## 容器能手动运行，Snakemake 却失败？

检查 `container.image` 是否指向新镜像、文件是否可读、外部目录是否 bind、所选可执行文件是否在容器中。旧的本地 `ParaChrSNP.sif` 可能仍是 minibwa 0.1 且缺少 samblaster；不要仅按文件名判断版本。

## Read the Docs 报缺少 Sphinx 扩展？

`conf.py` 中启用的 `myst_parser` 必须由 `source/requirements.txt` 安装。英文项目使用仓库根目录 `.readthedocs.yaml`；中文项目需在 Read the Docs 后台把 **Build configuration file** 明确设置为 `source_zh/.readthedocs.yaml`，不能只在仓库新增文件而不修改项目设置。

## 规则失败后看哪里？

先查看 `.snakemake/log/` 中最新日志，找第一个失败规则，再打开 `logs/` 下对应工具日志。Snakemake 最后的“非零退出码”通常不是根因；工具日志才有具体报错。
