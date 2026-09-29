# 可选分析模块

核心变异检测之外，ParaChrSNP 有 9 个可独立启用的模块。只有 `enabled: true` 的模块才会加入默认全流程目标；启用模块不会自动启用其他模块。以下配置均为**片段**，必须写入完整配置文件的 `params` 下，不能单独保存为可运行配置。

| 模块 | 核心问题 | 主要输入 | 主要输出 |
| --- | --- | --- | --- |
| Beagle | 缺失基因型能否由单倍型推断？ | 群体 SNP VCF | 填补后的 VCF |
| CNVnator | 哪些区域的读段深度提示缺失或扩增？ | 单样本 BAM | CNV 与群体汇总 |
| Pi | 群体内部多样性如何沿染色体变化？ | SNP VCF、可选群体文件 | 窗口 π 值 |
| SNP density | SNP 富集或贫乏区在哪里？ | SNP VCF | 窗口 SNP 数 |
| ADMIXTURE | 样本的祖源成分比例如何？ | PLINK 二进制数据 | 各 K 值的 Q 矩阵和 CV 误差 |
| PopLDdecay | 连锁不平衡如何随距离衰减？ | SNP VCF、可选群体文件 | LD 衰减曲线 |
| PCA | 遗传变异的主要轴是什么？ | SNP VCF | 主成分坐标和图 |
| VCF2Dis | 样本间遗传距离如何？ | SNP VCF | 距离矩阵和 Newick 树 |
| SnpEff | 变异影响哪些基因和转录本？ | VCF、对应组装的注释 | 功能注释 VCF |

## 群体信息文件

Pi、ADMIXTURE、PopLDdecay、PCA 和 VCF2Dis 可使用两列文本：第一列是 VCF 中的样本 ID，第二列是群体名。样本 ID 区分大小写。

```text
sample1    Population_A
sample2    Population_A
sample3    Population_B
```

Pi 和 PopLDdecay 的 `pop_info` 为空时会把所有样本视为一组；PCA 和 VCF2Dis 的 `sample_group` 为空时处理全部 VCF 样本。

## Beagle：缺失基因型填补

Beagle 利用样本间单倍型关联推断缺失基因型，适合需要更完整基因型矩阵的下游分析。但推断值不是直接测得的基因型；精度依赖样本数、标记密度、亲缘关系、等位基因频率和原始缺失率。应比较填补前后的缺失率，并结合概率信息评估不确定性。

```yaml
params:
  imputation:
    enabled: true
    input_vcf: "result_vcfs/combined.snp.filtered.vcf.gz"
    output_prefix: "imputation/combined.snp.filtered.beagle"
    jar: ""
    java_options: "-Xmx16g"
    threads: 8
    extra: ""
```

`input_vcf` 是待填补 VCF；`output_prefix` 决定 `.vcf.gz` 与 `.tbi` 输出名；`jar` 可明确指定 Beagle JAR，留空则搜索容器常见位置；`java_options` 控制 Java 内存；`threads` 传给 Beagle 的 `nthreads`；`extra` 用于遗传图谱、参考单倍型面板等额外参数。主要结果是 `imputation/combined.snp.filtered.beagle.vcf.gz`。

## CNVnator：测序深度 CNV

CNVnator 从 BAM 读段深度寻找可能的缺失和重复。覆盖度、比对质量、重复序列、参考组装质量和 bin 大小都会影响结果；CNV 边界是近似值，重要事件应独立验证。

```yaml
params:
  cnv:
    enabled: true
    software: "cnvnator"
    executable: "cnvnator"
    vcf_converter: "cnvnator2VCF.pl"
    bin_size: 100
    reference_dir: "reference"
    extra: ""
```

`software` 目前仅实现 `cnvnator`；`executable` 和 `vcf_converter` 指定程序；`bin_size` 是深度统计窗口，越小分辨率越高、噪声也越大且需要更高覆盖度；`reference_dir` 是参考 FASTA 所在目录；`extra` 传入附加参数。模块输出单样本 ROOT、文本、TSV、VCF，以及 `cnv/combined.cnv.tsv`、`cnv/combined.cnv.summary.tsv` 和图表。

## Pi：窗口核苷酸多样性

π 估计群体内每个位点的平均成对差异，而不是简单数 SNP。跨群体比较前，应统一样本数、缺失率、基因型过滤与可检测区域的处理方式。

```yaml
params:
  pi:
    enabled: true
    input_vcf: "result_vcfs/combined.snp.filtered.vcf.gz"
    pop_info: "pop.info"
    output_dir: "pi"
    window_size: 100000
    window_step: 10000
    extra: ""
```

`pop_info` 是可选群体文件，留空视为单一 `All` 组；`window_size` 为窗口长度（bp）；`window_step` 为相邻窗口起点间距，小于窗口长度时窗口重叠；`output_dir` 决定结果位置。主要结果为 `pi/combined.windowed.pi.tsv`、`pi/pi.summary.tsv` 和窗口图。

## SNP density：SNP 密度

按固定、不重叠的染色体窗口统计**保留的 SNP 记录数**。该数值不校正样本数、等位基因频率、缺失基因型和不可比对区，不能直接解释为核苷酸多样性。

```yaml
params:
  snp_density:
    enabled: true
    input_vcf: "result_vcfs/combined.snp.filtered.vcf.gz"
    output_dir: "snp_density"
    window_size: 1000000
```

`window_size` 为非重叠窗口宽度（bp）；模块还使用配置中的染色体列表及参考 `.fai`。输出 `snp_density/snp_density.tsv` 和 PDF/PNG/SVG/TIFF 图。

## ADMIXTURE：群体结构

ADMIXTURE 在不同祖源成分数 K 下估计各样本比例。流程先处理 PLINK 文件、按缺失率筛选并进行 LD 剪枝，再运行各 K 值。CV 误差较低可帮助比较模型，但还要考虑群体采样不均、亲缘关系及重复运行稳定性。

```yaml
params:
  admixture:
    enabled: true
    input_prefix: "format_convert/combined.snp.filtered"
    output_dir: "admixture"
    executable: "admixture"
    k_min: 1
    k_max: 10
    cv: 10
    threads: 8
    prune_window: 50
    prune_step: 10
    prune_r2: 0.2
    geno: 0.1
    pop_info: "pop.info"
    show_sample_names: true
    plink_extra: "--allow-extra-chr"
    normalize_extra: "--allow-extra-chr --set-missing-var-ids @:#"
    admixture_plink_extra: "--allow-extra-chr 0"
    extra: ""
```

`input_prefix` 不带 `.bed/.bim/.fam`；`k_min/k_max` 给出 K 范围；`cv` 是交叉验证折数；`prune_window/prune_step/prune_r2` 控制 PLINK LD 剪枝；`geno` 是位点最大缺失率，例如 `0.1` 表示剔除缺失率大于 10% 的位点；`pop_info` 决定绘图分组；`show_sample_names` 控制是否显示样本标签。输出各 K 值 Q 矩阵、`admixture/cv_errors.tsv`、群体结构图及 CV 曲线。

## PopLDdecay：LD 衰减

PopLDdecay 统计位点对之间的连锁不平衡及其随物理距离的变化。曲线同时受重组、群体历史、选择、位点筛选和样本量影响；比较群体时应保持过滤条件一致。

```yaml
params:
  ld_decay:
    enabled: true
    input_vcf: "result_vcfs/combined.snp.filtered.vcf.gz"
    pop_info: "pop.info"
    output_dir: "ld_decay"
    executable: "PopLDdecay"
    max_dist: 300
    threads: 1
    extra: ""
```

`max_dist` 传给 PopLDdecay 的 `-MaxDist`，单位为 kb；`threads` 是 Snakemake 资源声明，当前命令**没有**将它作为 PopLDdecay 线程选项；`pop_info` 留空则合并所有样本。输出每组 `.stat.gz`、`ld_decay/combined.ld_decay.tsv` 和曲线图。

## PCA：主成分分析

PCA 将相关基因型变异压缩为正交主成分，有助于识别大尺度群体结构、离群样本和批次效应。当前模块直接使用筛选后的 SNP VCF；高密度连锁位点可能主导结果，敏感分析可另行使用 LD 剪枝后的数据。至少需要 3 个样本。

```yaml
params:
  vcf2pca:
    enabled: true
    executable: "VCF2PCACluster"
    sample_group: "pop.info"
    output_prefix: "pca/ParaChrSNP"
    plot_prefix: "pca/ParaChrSNP.plot"
    plot2_executable: "Plot2Deig"
    plot3_executable: "Plot3Deig"
    plot3_preview_script: "scripts/pca3d-preview_1.R"
    threads: 8
    extra: ""
```

`sample_group` 可用于分组着色；`output_prefix` 控制特征向量/值文件；`plot_prefix` 控制图片；`plot2_executable/plot3_executable` 和 `plot3_preview_script` 指定绘图程序；`threads` 为分配线程。输出 `.eigenvec`、`.eigenval`、二维和三维图。

## VCF2Dis：遗传距离与树

VCF2Dis 根据筛选后的 SNP 计算样本两两遗传距离并构建 Newick 树。该树用于探索样本关系，不等于已充分评估支持度的系统发育推断；解释时需考虑缺失值、连锁、距离定义和分支支持。至少需要 3 个样本。

```yaml
params:
  vcf2dis:
    enabled: true
    executable: "VCF2Dis"
    sample_group: "pop.info"
    output_matrix: "dis/ParaChrSNP.p_dis.mat"
    output_tree: "dis/ParaChrSNP.p_dis.nwk"
    tree_method: 1
    extra: ""
```

`sample_group` 可限制或标记样本；`output_matrix/output_tree` 是结果路径；`tree_method` 通过 `-TreeMethod` 传给 VCF2Dis。输出距离矩阵和 Newick 树。

## SnpEff：变异功能注释

SnpEff 根据基因结构预测 SNP/INDEL 对基因、转录本及蛋白的影响。流程先从参考 FASTA 和 GFF3/GTF 构建自定义数据库。FASTA、注释和 VCF 必须来自**同一组装版本**，序列 ID 也要完全一致，否则预测会误导。

```yaml
params:
  snpeff:
    enabled: true
    annotate_snp: true
    annotate_indel: true
    executable: "snpEff"
    genome_name: "custom_genome"
    data_dir: "annotation/snpeff_data"
    config_file: "annotation/snpeff.config"
    genome_fasta: "reference/genome.fa"
    annotation_file: "reference/genes.gff3"
    annotation_format: "gff3"
    output_prefix: "annotation/combined"
    database_done: "annotation/snpeff_db.done"
    java_options: "-Xmx16g"
    build_check_options: "-noCheckCds -noCheckProtein"
    threads: 1
    extra: ""
```

`annotate_snp/annotate_indel` 选择输出类别；`genome_name` 是数据库 ID；`data_dir/config_file` 是数据库与配置路径；`annotation_format` 可取 `gff3` 或 `gtf`；`java_options` 控制内存；`build_check_options` 中的宽松检查仅在缺少 CDS/蛋白参考时使用。每类变异输出压缩且索引的注释 VCF、HTML 统计与基因汇总。

## 实践建议

- 首先运行 `reports/precheck.done`，解决所有路径、样本名和参考不一致问题。
- 位点级质量过滤不等于针对研究目的的基因型 DP/GQ、缺失率或等位基因频率过滤。
- ADMIXTURE 使用 LD 剪枝位点；PCA 和距离分析也要警惕密集连锁位点的影响。
- 比较群体间 Pi、LD 或祖源比例之前，检查样本量和缺失率是否平衡。
- 仅用与比对参考基因组相同组装版本的注释运行 SnpEff。
