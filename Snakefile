configfile: "config.yaml"

import re

container: config.get("container", {}).get("image", "ParaChrSNP.sif")

def reference_dict(reference):
    return re.sub(r"\.(fa|fasta)$", ".dict", reference)

SAMPLE_PATTERN = "|".join(re.escape(sample) for sample in config["samples"])
CHROM_PATTERN = "|".join(re.escape(chrom) for chrom in config["chromosomes"])

OPTIONAL_TARGETS = []

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
        PCA_PLOT_PREFIX + ".C.3DPC1PC2PC3.pdf",
        PCA_PLOT_PREFIX + ".N.3DPC1PC2PC3.pdf",
    ])

if config["params"]["vcf2dis"].get("enabled", True) and len(config["samples"]) >= 3:
    OPTIONAL_TARGETS.extend([
        config["params"]["vcf2dis"]["output_matrix"],
        config["params"]["vcf2dis"]["output_tree"],
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
include: "rules/vcf2pca.rules"
include: "rules/vcf2dis.rules"
include: "rules/report.rules"
