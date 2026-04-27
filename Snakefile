configfile: "config.yaml"

import os
import re

container: config.get("container", {}).get("image", "ParaChrSNP.sif")

def reference_dict(reference):
    return re.sub(r"\.(fa|fasta)$", ".dict", reference)

SAMPLE_PATTERN = "|".join(re.escape(sample) for sample in config["samples"])
CHROM_PATTERN = "|".join(re.escape(chrom) for chrom in config["chromosomes"])

OPTIONAL_TARGETS = []

SNPEFF_CONFIG = config["params"].setdefault("snpeff", {})
SNPEFF_CONFIG.setdefault("enabled", False)
SNPEFF_CONFIG.setdefault("annotate_snp", True)
SNPEFF_CONFIG.setdefault("annotate_indel", False)
SNPEFF_CONFIG.setdefault("executable", "snpEff")
SNPEFF_CONFIG.setdefault("genome_name", "custom_genome")
SNPEFF_CONFIG.setdefault("data_dir", "annotation/snpeff_data")
SNPEFF_CONFIG.setdefault("config_file", "annotation/snpeff.config")
SNPEFF_CONFIG.setdefault("genome_fasta", config["reference"])
SNPEFF_CONFIG.setdefault("annotation_file", "annotation/genes.gff3")
SNPEFF_CONFIG.setdefault("annotation_format", "gff3")
SNPEFF_CONFIG.setdefault("output_prefix", "annotation/combined")
SNPEFF_CONFIG.setdefault("database_done", "annotation/snpeff_db.done")
SNPEFF_CONFIG.setdefault("java_options", "-Xmx8g")
SNPEFF_CONFIG.setdefault("build_check_options", "-noCheckCds -noCheckProtein")
SNPEFF_CONFIG.setdefault("extra", "")

IMPUTATION_CONFIG = config["params"].setdefault("imputation", {})
IMPUTATION_CONFIG.setdefault("enabled", False)
IMPUTATION_CONFIG.setdefault("input_vcf", "result_vcfs/combined.snp.filtered.vcf.gz")
IMPUTATION_CONFIG.setdefault("output_prefix", "imputation/combined.snp.filtered.beagle")
IMPUTATION_CONFIG.setdefault("jar", "")
IMPUTATION_CONFIG.setdefault("java_options", "-Xmx4g")
IMPUTATION_CONFIG.setdefault("threads", 4)
IMPUTATION_CONFIG.setdefault("extra", "")

if config["params"].get("imputation", {}).get("enabled", False):
    IMPUTATION_PREFIX = config["params"]["imputation"].get(
        "output_prefix",
        "imputation/combined.snp.filtered.beagle",
    )
    OPTIONAL_TARGETS.extend([
        IMPUTATION_PREFIX + ".vcf.gz",
        IMPUTATION_PREFIX + ".vcf.gz.tbi",
    ])

if config["params"]["vcf2pca"].get("enabled", True) and len(config["samples"]) >= 3:
    PCA_PLOT_PREFIX = config["params"]["vcf2pca"].get(
        "plot_prefix",
        config["params"]["vcf2pca"]["output_prefix"] + ".plot",
    )
    OPTIONAL_TARGETS.extend([
        config["params"]["vcf2pca"]["output_prefix"] + ".eigenvec",
        config["params"]["vcf2pca"]["output_prefix"] + ".eigenval",
        PCA_PLOT_PREFIX + ".C.PC1_PC2.p.svg",
        PCA_PLOT_PREFIX + ".N.PC1_PC2.p.svg",
        PCA_PLOT_PREFIX + ".C.PC1_PC2.p.png",
        PCA_PLOT_PREFIX + ".N.PC1_PC2.p.png",
        PCA_PLOT_PREFIX + ".C.3DPC1PC2PC3.pdf",
        PCA_PLOT_PREFIX + ".N.3DPC1PC2PC3.pdf",
    ])

if config["params"]["vcf2dis"].get("enabled", True) and len(config["samples"]) >= 3:
    OPTIONAL_TARGETS.extend([
        config["params"]["vcf2dis"]["output_matrix"],
        config["params"]["vcf2dis"]["output_tree"],
    ])

if config["params"].get("snpeff", {}).get("enabled", False):
    SNPEFF_PREFIX = config["params"]["snpeff"].get("output_prefix", "annotation/combined")
    OPTIONAL_TARGETS.append(
        config["params"]["snpeff"].get("database_done", "annotation/snpeff_db.done")
    )
    if config["params"]["snpeff"].get("annotate_snp", True):
        OPTIONAL_TARGETS.extend([
            SNPEFF_PREFIX + ".snp.snpeff.vcf.gz",
            SNPEFF_PREFIX + ".snp.snpeff.vcf.gz.tbi",
            SNPEFF_PREFIX + ".snp.snpeff.html",
            SNPEFF_PREFIX + ".snp.snpeff.genes.txt",
        ])
    if config["params"]["snpeff"].get("annotate_indel", False):
        OPTIONAL_TARGETS.extend([
            SNPEFF_PREFIX + ".indel.snpeff.vcf.gz",
            SNPEFF_PREFIX + ".indel.snpeff.vcf.gz.tbi",
            SNPEFF_PREFIX + ".indel.snpeff.html",
            SNPEFF_PREFIX + ".indel.snpeff.genes.txt",
        ])

wildcard_constraints:
    sample=SAMPLE_PATTERN,
    chrom=CHROM_PATTERN


rule all:
    input:
        "reports/precheck.done",
        "qc/rastqc.done",
        expand("gvcf/{sample}.{chrom}.g.vcf.gz", sample=config["samples"], chrom=config["chromosomes"]),
        expand("gvcf/{sample}.g.vcf.gz", sample=config["samples"]),
        "result_vcfs/combined.vcf.gz",
        "result_vcfs/combined.indel.filtered.vcf.gz",
        "result_vcfs/combined.snp.filtered.vcf.gz",
        config["params"]["vcf_missing"]["plot_prefix"] + ".sample_missing.pdf",
        config["params"]["vcf_missing"]["plot_prefix"] + ".Variant_missing.pdf",
        config["params"]["vcf_convert"]["output_prefix"] + ".bed",
        config["params"]["vcf_convert"]["output_prefix"] + ".bim",
        config["params"]["vcf_convert"]["output_prefix"] + ".fam",
        config["params"]["vcf_convert"]["output_prefix"] + ".ped",
        config["params"]["vcf_convert"]["output_prefix"] + ".map",
        config["params"]["vcf_convert"]["output_prefix"] + ".hmp.txt",
        OPTIONAL_TARGETS,
        "reports/ParaChrSNP_report.html",
        "reports/ParaChrSNP_summary.tsv"

include: "rules/precheck.rules"
include: "rules/bam_rmdup.rules"
include: "rules/bwa_index.rules"
include: "rules/bwa_mem.rules"
include: "rules/index_rmdup.rules"
include: "rules/qc.rules"
include: "rules/clean_reads.rules"
include: "rules/combine_gvcf.rules"
include: "rules/faidx_index.rules"
include: "rules/haplo.rules"
include: "rules/indel_filter.rules"
include: "rules/indel_select.rules"
include: "rules/joint_calling.rules"
include: "rules/picard_index.rules"
include: "rules/samtools_index.rules"
include: "rules/snp_filter.rules"
include: "rules/snp_select.rules"
include: "rules/merge_sample_gvcf.rules"
include: "rules/get_chr_list.rules"
include: "rules/vcf_missing.rules"
include: "rules/vcf_convert.rules"
include: "rules/imputation.rules"
include: "rules/vcf2pca.rules"
include: "rules/vcf2dis.rules"
include: "rules/snpeff.rules"
include: "rules/report.rules"
