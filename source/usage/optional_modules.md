# Optional analysis modules

ParaChrSNP provides nine optional analyses after the core read-mapping and
variant-calling workflow. Each module is controlled independently with an
`enabled` switch, so enabling one module does not automatically enable the
others.

## Choosing modules

| Module | Biological question | Primary input | Main result |
| --- | --- | --- | --- |
| Beagle | Can missing genotypes be statistically inferred? | Filtered cohort SNP VCF | Imputed VCF |
| CNVnator | Which genomic regions show read-depth gains or losses? | Deduplicated BAM files | Per-sample CNV calls and cohort summaries |
| Pi | How much nucleotide diversity occurs within each population? | Filtered SNP VCF and optional population file | Windowed nucleotide diversity |
| SNP density | Where are SNP-rich and SNP-poor regions? | Filtered SNP VCF | SNP counts in genomic windows |
| ADMIXTURE | How many ancestry components explain the samples? | PLINK binary files | Q matrices, CV errors, and structure plots |
| PopLDdecay | How rapidly does linkage disequilibrium decay with distance? | Filtered SNP VCF and optional population file | LD-decay statistics and curves |
| PCA | What are the main axes of genetic variation? | Filtered SNP VCF | Eigenvectors, eigenvalues, and plots |
| VCF2Dis | What are the pairwise genetic distances and sample relationships? | Filtered SNP VCF | Distance matrix and Newick tree |
| SnpEff | Which genes and molecular consequences are affected by variants? | Filtered SNP/INDEL VCF and genome annotation | Functionally annotated VCF |

An `enabled: true` setting adds that module's final files to Snakemake's
default targets. An `enabled: false` setting excludes the module from a normal
full run.

## Population information file

Pi, ADMIXTURE, PopLDdecay, PCA, and VCF2Dis can use a two-column sample-group
file. The first column is the VCF sample identifier and the second is the
population or group name.

```text
sample1    Population_A
sample2    Population_A
sample3    Population_B
```

Sample identifiers are case-sensitive and must match the VCF header. Pi and
PopLDdecay treat all samples as one group when `pop_info` is empty. PCA and
VCF2Dis analyse all VCF samples when `sample_group` is empty.

## Beagle genotype imputation

Beagle uses linkage information among samples to infer missing genotypes. It
is useful when downstream population analyses require a more complete
genotype matrix. Imputation does not create experimentally observed data:
accuracy depends on sample number, marker density, relatedness, allele
frequency, and the amount of missing data. Always inspect missingness before
and after imputation and retain genotype-probability information when judging
uncertain calls.

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

| Parameter | Meaning |
| --- | --- |
| `enabled` | Adds the imputed VCF and index to the workflow targets. |
| `input_vcf` | Cohort VCF to impute; normally the filtered SNP VCF. |
| `output_prefix` | Prefix used for the Beagle `.vcf.gz` output. |
| `jar` | Explicit Beagle JAR path. If empty, the rule searches common container installation directories. |
| `java_options` | JVM memory and other Java settings. |
| `threads` | Value passed to Beagle as `nthreads`. |
| `extra` | Additional Beagle arguments, such as a reference panel or genetic map. |

Main outputs are
`imputation/combined.snp.filtered.beagle.vcf.gz` and its `.tbi` index.

## CNVnator read-depth CNV calling

CNVnator detects deletions and duplications from read-depth changes. Unlike
the SNP-based modules, it reads
`duplicate_removed/{sample}.rmdup.bam` and the corresponding BAM index for
each sample. Coverage, mapping quality, repeat content, genome assembly
quality, and the selected bin size all affect CNV calls. CNV boundaries are
approximate and important events should be validated with an independent
method.

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

| Parameter | Meaning |
| --- | --- |
| `enabled` | Enables per-sample CNV calling and cohort summaries. |
| `software` | CNV backend label; currently `cnvnator` is implemented. |
| `executable` | CNVnator executable name or path. |
| `vcf_converter` | Path or command for `cnvnator2VCF.pl`. |
| `bin_size` | Read-depth bin size in bp. Smaller bins provide finer resolution but require higher depth and are noisier. |
| `reference_dir` | Directory containing the reference FASTA used during CNVnator processing. |
| `extra` | Additional CNVnator arguments. |

The rule creates per-sample ROOT, text, TSV, and VCF files, plus
`cnv/combined.cnv.tsv`, `cnv/combined.cnv.summary.tsv`, and cohort plots.

## Window-based nucleotide diversity

Pi estimates average pairwise nucleotide differences within each population
in sliding genomic windows. It is intended for within-population diversity,
not simple variant counting. Comparisons among groups are most meaningful
when sample sizes, genotype filtering, missingness, and callable genomic
regions are handled consistently.

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

| Parameter | Meaning |
| --- | --- |
| `enabled` | Enables the Pi table, summary, and plots. |
| `input_vcf` | Filtered cohort SNP VCF. |
| `pop_info` | Optional two-column sample/population file; empty means one `All` group. |
| `output_dir` | Directory for Pi results. |
| `window_size` | Genomic window length in bp. |
| `window_step` | Distance in bp between consecutive window starts; a value smaller than the window size creates overlapping windows. |
| `extra` | Additional arguments passed to the Pi calculation script. |

Main outputs are `pi/combined.windowed.pi.tsv`, `pi/pi.summary.tsv`, and Pi
Manhattan-style plots in PDF, PNG, SVG, and TIFF formats.

## SNP density

SNP density counts retained SNP records in fixed, non-overlapping chromosome
windows. It highlights variant-rich and variant-poor regions but does not
correct for allele frequency, sample size, missing genotypes, or inaccessible
sequence; it must therefore not be interpreted as nucleotide diversity.

```yaml
params:
    snp_density:
        enabled: true
        input_vcf: "result_vcfs/combined.snp.filtered.vcf.gz"
        output_dir: "snp_density"
        window_size: 1000000
```

| Parameter | Meaning |
| --- | --- |
| `enabled` | Enables the SNP-density table and plots. |
| `input_vcf` | Filtered cohort SNP VCF to count. |
| `output_dir` | Directory for tables and figures. |
| `window_size` | Width of each non-overlapping window in bp. |

The module uses the configured chromosome list and reference `.fai`. Its main
outputs are `snp_density/snp_density.tsv` and figures in four formats.

## ADMIXTURE population structure

ADMIXTURE estimates each sample's membership proportions across a specified
range of ancestral components (`K`). ParaChrSNP first normalizes the PLINK
dataset, filters variants by missingness, performs LD pruning, and then runs
ADMIXTURE for every K. The lowest cross-validation error is useful for model
comparison, but biological interpretation must also consider sampling design,
unequal population sizes, related individuals, and stability across repeated
runs.

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

| Parameter | Meaning |
| --- | --- |
| `input_prefix` | PLINK prefix without `.bed`, `.bim`, or `.fam`; normally produced by `vcf_convert`. |
| `k_min`, `k_max` | Inclusive range of ancestry-component counts. |
| `cv` | Number of ADMIXTURE cross-validation folds. |
| `threads` | Threads used by PLINK and ADMIXTURE. |
| `prune_window` | PLINK LD-pruning window size in variants. |
| `prune_step` | Number of variants by which the pruning window advances. |
| `prune_r2` | Pairwise LD threshold; variants above this threshold are pruned. |
| `geno` | PLINK maximum variant missing-call rate. For example, `0.1` removes variants missing in more than 10% of samples. |
| `pop_info` | Optional group file used to order and colour samples in plots. |
| `show_sample_names` | Shows sample labels when `true`; dense plots may be clearer with `false`. |
| `*_extra` | Additional backend-specific PLINK or ADMIXTURE arguments. |

Outputs include the LD-pruned PLINK dataset, a Q matrix for every K,
`admixture/cv_errors.tsv`, structure plots, and CV-error plots.

## PopLDdecay linkage-disequilibrium decay

PopLDdecay calculates pairwise LD and summarizes how LD decreases with
physical distance. The decay curve can reflect recombination, demographic
history, selection, marker ascertainment, and sample size. Compare groups
using consistent site filters and sufficient numbers of individuals.

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

| Parameter | Meaning |
| --- | --- |
| `input_vcf` | Filtered cohort SNP VCF. |
| `pop_info` | Optional two-column population file; empty analyses all samples together. |
| `output_dir` | Directory for statistics, sample lists, logs, and figures. |
| `executable` | PopLDdecay command or path. |
| `max_dist` | Maximum pairwise distance passed to PopLDdecay with `-MaxDist` (in kb for PopLDdecay). |
| `threads` | Snakemake resource allocation for this rule; the current command does not pass a thread option to PopLDdecay. |
| `extra` | Additional PopLDdecay arguments. |

Main outputs are per-group `.stat.gz` files,
`ld_decay/combined.ld_decay.tsv`, and decay plots in four formats.

## PCA

PCA summarizes correlated genotype variation into orthogonal principal
components. It is useful for identifying broad genetic structure, outliers,
and batch effects. ParaChrSNP currently runs PCA directly on the filtered SNP
VCF; for analyses sensitive to tightly linked markers, consider supplying a
suitably LD-pruned dataset through a separate analysis.

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

| Parameter | Meaning |
| --- | --- |
| `executable` | VCF2PCACluster executable name or path. |
| `sample_group` | Optional sample-group file used for labels and colours. |
| `output_prefix` | Prefix for eigenvector and eigenvalue files. |
| `plot_prefix` | Prefix for 2D and 3D figures. |
| `plot2_executable`, `plot3_executable` | VCF2PCACluster plotting commands. |
| `plot3_preview_script` | R script used to generate static 3D previews. |
| `threads` | Threads allocated to VCF2PCACluster. |
| `extra` | Additional VCF2PCACluster arguments. |

At least three configured samples are required. Outputs include `.eigenvec`,
`.eigenval`, PC1-PC2 plots, and PC1-PC2-PC3 plots.

## VCF2Dis distance matrix and tree

VCF2Dis calculates pairwise genetic distances from the filtered SNP VCF and
constructs a Newick tree. It is suitable for exploratory sample relationships;
the resulting topology is not automatically a fully specified phylogenetic
inference and should be interpreted with the distance definition, missing-data
handling, marker linkage, and branch support in mind.

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

| Parameter | Meaning |
| --- | --- |
| `executable` | VCF2Dis executable name or path. |
| `sample_group` | Optional file restricting or grouping samples. |
| `output_matrix` | Pairwise-distance matrix path. |
| `output_tree` | Newick tree path. |
| `tree_method` | Integer tree-building mode passed to VCF2Dis with `-TreeMethod`. |
| `extra` | Additional VCF2Dis arguments. |

At least three configured samples are required.

## SnpEff variant annotation

SnpEff predicts variant effects on genes, transcripts, proteins, and other
annotated features. ParaChrSNP builds a custom SnpEff database from the
configured FASTA and GFF3/GTF file before annotation. Reference sequence
identifiers, annotation sequence identifiers, and the chromosome list must
match exactly. Incorrect transcript structures or mismatched assemblies will
produce misleading effect predictions.

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

| Parameter | Meaning |
| --- | --- |
| `annotate_snp` | Produces an annotated filtered SNP VCF when `true`. |
| `annotate_indel` | Produces an annotated filtered INDEL VCF when `true`. |
| `genome_name` | Internal custom database identifier written to the SnpEff config. |
| `data_dir` | Root directory of the custom SnpEff database. |
| `config_file` | Generated local SnpEff configuration file. |
| `genome_fasta` | FASTA from the same assembly used for variant calling. |
| `annotation_file` | Matching GFF3 or GTF gene annotation. |
| `annotation_format` | `gff3` or `gtf`, controlling the SnpEff build mode. |
| `output_prefix` | Prefix for annotated VCFs and reports. |
| `database_done` | Snakemake completion marker for database construction. |
| `java_options` | JVM memory settings for database building and annotation. |
| `build_check_options` | SnpEff database-build validation options. Remove permissive options when CDS/protein references are available for strict checking. |
| `threads` | Snakemake resource allocation; SnpEff execution is primarily controlled by its own options. |
| `extra` | Additional SnpEff annotation arguments. |

For each selected variant class, the module writes a bgzip-compressed and
indexed annotated VCF, an HTML statistics report, and a gene summary table.

## Practical recommendations

- Run `reports/precheck.done` before a full analysis and resolve all reported
  path or sample-name problems.
- Apply consistent genotype-level missingness, depth, quality, and allele-
  frequency filters when required by the biological analysis; the core hard
  filters are site-level GATK filters.
- Use LD-pruned markers for ADMIXTURE and consider LD pruning for PCA and
  distance analyses when dense linked SNPs could dominate the result.
- Do not compare Pi, LD decay, or ADMIXTURE results among groups without
  checking sample-size imbalance and missingness.
- Enable SnpEff only with an annotation file belonging to the exact reference
  assembly used for mapping and calling.
