configfile: "config.yaml"

import os
import re

container: config.get("container", {}).get("image", "ParaChrSNP.sif")

def reference_dict(reference):
    return re.sub(r"\.(fa|fasta)$", ".dict", reference)

def java_options_with_active_processors(java_options, threads=1):
    """Limit JVM CPU discovery unless the user already set it explicitly."""
    options = str(java_options).strip()
    if "-XX:ActiveProcessorCount" in options:
        return options

    active_option = f"-XX:ActiveProcessorCount={int(threads)}"
    match = re.search(r'--java-options\s+(["\'])(.*?)\1', options)
    if match:
        quote = match.group(1)
        current = match.group(2).strip()
        replacement = f"--java-options {quote}{current} {active_option}{quote}"
        return options[:match.start()] + replacement + options[match.end():]

    return f"{active_option} {options}".strip()

def scaled_bwa_threads(total_threads, bwa_threads, sort_threads):
    """Keep piped bwa-mem2 and samtools sort thread use within Snakemake's allocation."""
    total_threads = int(total_threads)
    bwa_threads = int(bwa_threads)
    sort_threads = int(sort_threads)
    if total_threads <= 1:
        return {"bwa": 1, "sort": 0}
    scaled_sort = min(sort_threads, total_threads - 1)
    scaled_bwa = min(bwa_threads, total_threads - scaled_sort)
    return {"bwa": max(1, scaled_bwa), "sort": max(0, scaled_sort)}

SAMPLE_PATTERN = "|".join(re.escape(sample) for sample in config["samples"])
CHROM_PATTERN = "|".join(re.escape(chrom) for chrom in config["chromosomes"])

OPTIONAL_TARGETS = []

JOINT_CALLING_CONFIG = config["params"].setdefault("joint_calling", {})
JOINT_CALLING_CONFIG.setdefault("method", "genomicsdb")
JOINT_CALLING_CONFIG.setdefault("reader_threads", 4)
JOINT_CALLING_CONFIG.setdefault("genotype_threads", 1)
JOINT_CALLING_CONFIG.setdefault("gather_threads", 1)
JOINT_CALLING_CONFIG.setdefault("batch_size", 50)
JOINT_CALLING_CONFIG.setdefault("import_java_options", '--java-options "-Xms1g -Xmx16g"')
JOINT_CALLING_CONFIG.setdefault("genotype_java_options", config["params"]["genotype_gvcfs"]["java_options"])
JOINT_CALLING_CONFIG.setdefault("gather_java_options", config["params"]["combine_vcf"]["java_options"])
JOINT_CALLING_CONFIG.setdefault("extra_import", "")
JOINT_CALLING_CONFIG.setdefault("extra_genotype", "")

JOINT_CALLING_METHOD = JOINT_CALLING_CONFIG.get("method", "genomicsdb")
if JOINT_CALLING_METHOD not in ["genomicsdb", "combine_gvcfs"]:
    raise ValueError("params.joint_calling.method must be 'genomicsdb' or 'combine_gvcfs'")

GVCF_TARGETS = expand(
    "gvcf/{sample}.{chrom}.g.vcf.gz",
    sample=config["samples"],
    chrom=config["chromosomes"],
)
if JOINT_CALLING_METHOD == "combine_gvcfs":
    GVCF_TARGETS.extend(expand("gvcf/{sample}.g.vcf.gz", sample=config["samples"]))

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

CNV_CONFIG = config["params"].setdefault("cnv", {})
CNV_CONFIG.setdefault("enabled", False)
CNV_CONFIG.setdefault("software", "cnvnator")
CNV_CONFIG.setdefault("executable", "cnvnator")
CNV_CONFIG.setdefault("vcf_converter", "cnvnator2VCF.pl")
CNV_CONFIG.setdefault("bin_size", 100)
CNV_CONFIG.setdefault("reference_dir", os.path.dirname(config["reference"]) or ".")
CNV_CONFIG.setdefault("extra", "")

PI_CONFIG = config["params"].setdefault("pi", {})
PI_CONFIG.setdefault("enabled", False)
PI_CONFIG.setdefault("input_vcf", "result_vcfs/combined.snp.filtered.vcf.gz")
PI_CONFIG.setdefault("pop_info", "")
PI_CONFIG.setdefault("output_dir", "pi")
PI_CONFIG.setdefault("window_size", 100000)
PI_CONFIG.setdefault("window_step", 10000)
PI_CONFIG.setdefault("extra", "")

SNP_DENSITY_CONFIG = config["params"].setdefault("snp_density", {})
SNP_DENSITY_CONFIG.setdefault("enabled", False)
SNP_DENSITY_CONFIG.setdefault("input_vcf", "result_vcfs/combined.snp.filtered.vcf.gz")
SNP_DENSITY_CONFIG.setdefault("output_dir", "snp_density")
SNP_DENSITY_CONFIG.setdefault("window_size", 1000000)

ADMIXTURE_CONFIG = config["params"].setdefault("admixture", {})
ADMIXTURE_CONFIG.setdefault("enabled", False)
ADMIXTURE_CONFIG.setdefault("input_prefix", config["params"]["vcf_convert"]["output_prefix"])
ADMIXTURE_CONFIG.setdefault("output_dir", "admixture")
ADMIXTURE_CONFIG.setdefault("executable", "admixture")
ADMIXTURE_CONFIG.setdefault("k_min", 1)
ADMIXTURE_CONFIG.setdefault("k_max", 10)
ADMIXTURE_CONFIG.setdefault("cv", 10)
ADMIXTURE_CONFIG.setdefault("threads", 4)
ADMIXTURE_CONFIG.setdefault("prune_window", 50)
ADMIXTURE_CONFIG.setdefault("prune_step", 10)
ADMIXTURE_CONFIG.setdefault("prune_r2", 0.2)
ADMIXTURE_CONFIG.setdefault("geno", 0.999)
ADMIXTURE_CONFIG.setdefault("pop_info", "")
ADMIXTURE_CONFIG.setdefault("show_sample_names", False)
ADMIXTURE_CONFIG.setdefault("plink_extra", config["params"]["vcf_convert"].get("plink_extra", "--allow-extra-chr"))
ADMIXTURE_CONFIG.setdefault("normalize_extra", "--allow-extra-chr --set-missing-var-ids @:#")
ADMIXTURE_CONFIG.setdefault("admixture_plink_extra", "--allow-extra-chr 0")
ADMIXTURE_CONFIG.setdefault("extra", "")
if int(ADMIXTURE_CONFIG["k_min"]) > int(ADMIXTURE_CONFIG["k_max"]):
    raise ValueError("params.admixture.k_min must be <= params.admixture.k_max")

LD_DECAY_CONFIG = config["params"].setdefault("ld_decay", {})
LD_DECAY_CONFIG.setdefault("enabled", False)
LD_DECAY_CONFIG.setdefault("input_vcf", "result_vcfs/combined.snp.filtered.vcf.gz")
LD_DECAY_CONFIG.setdefault("pop_info", "")
LD_DECAY_CONFIG.setdefault("output_dir", "ld_decay")
LD_DECAY_CONFIG.setdefault("executable", "PopLDdecay")
LD_DECAY_CONFIG.setdefault("max_dist", 300)
LD_DECAY_CONFIG.setdefault("threads", 1)
LD_DECAY_CONFIG.setdefault("extra", "")

if config["params"].get("imputation", {}).get("enabled", False):
    IMPUTATION_PREFIX = config["params"]["imputation"].get(
        "output_prefix",
        "imputation/combined.snp.filtered.beagle",
    )
    OPTIONAL_TARGETS.extend([
        IMPUTATION_PREFIX + ".vcf.gz",
        IMPUTATION_PREFIX + ".vcf.gz.tbi",
    ])

if config["params"].get("cnv", {}).get("enabled", False):
    OPTIONAL_TARGETS.extend([
        "cnv/combined.cnv.tsv",
        "cnv/combined.cnv.summary.tsv",
        "cnv/figures/cnv_count_by_sample.pdf",
        "cnv/figures/cnv_length_distribution.pdf",
        "cnv/figures/cnv_count_by_chromosome.pdf",
    ])
    OPTIONAL_TARGETS.extend(expand("cnv/{sample}.cnv.vcf", sample=config["samples"]))

if PI_CONFIG.get("enabled", False):
    PI_OUTPUT_DIR = PI_CONFIG.get("output_dir", "pi")
    OPTIONAL_TARGETS.extend([
        PI_OUTPUT_DIR + "/combined.windowed.pi.tsv",
        PI_OUTPUT_DIR + "/pi.summary.tsv",
        PI_OUTPUT_DIR + "/figures/pi_manhattan.pdf",
        PI_OUTPUT_DIR + "/figures/pi_manhattan.png",
        PI_OUTPUT_DIR + "/figures/pi_manhattan.svg",
        PI_OUTPUT_DIR + "/figures/pi_manhattan.tiff",
    ])

if SNP_DENSITY_CONFIG.get("enabled", False):
    SNP_DENSITY_OUTPUT_DIR = SNP_DENSITY_CONFIG.get("output_dir", "snp_density")
    OPTIONAL_TARGETS.extend([
        SNP_DENSITY_OUTPUT_DIR + "/snp_density.tsv",
        SNP_DENSITY_OUTPUT_DIR + "/figures/snp_density.pdf",
        SNP_DENSITY_OUTPUT_DIR + "/figures/snp_density.png",
        SNP_DENSITY_OUTPUT_DIR + "/figures/snp_density.svg",
        SNP_DENSITY_OUTPUT_DIR + "/figures/snp_density.tiff",
    ])

if ADMIXTURE_CONFIG.get("enabled", False):
    ADMIXTURE_OUTPUT_DIR = ADMIXTURE_CONFIG.get("output_dir", "admixture")
    ADMIXTURE_K_VALUES = list(range(int(ADMIXTURE_CONFIG.get("k_min", 1)), int(ADMIXTURE_CONFIG.get("k_max", 10)) + 1))
    OPTIONAL_TARGETS.extend([
        ADMIXTURE_OUTPUT_DIR + "/pruned/admixture_pruned.bed",
        ADMIXTURE_OUTPUT_DIR + "/pruned/admixture_pruned.bim",
        ADMIXTURE_OUTPUT_DIR + "/pruned/admixture_pruned.fam",
        ADMIXTURE_OUTPUT_DIR + "/cv_errors.tsv",
        ADMIXTURE_OUTPUT_DIR + "/figures/admixture_structure.pdf",
        ADMIXTURE_OUTPUT_DIR + "/figures/admixture_structure.png",
        ADMIXTURE_OUTPUT_DIR + "/figures/admixture_structure.svg",
        ADMIXTURE_OUTPUT_DIR + "/figures/admixture_structure.tiff",
        ADMIXTURE_OUTPUT_DIR + "/figures/admixture_cv_error.pdf",
        ADMIXTURE_OUTPUT_DIR + "/figures/admixture_cv_error.png",
        ADMIXTURE_OUTPUT_DIR + "/figures/admixture_cv_error.svg",
        ADMIXTURE_OUTPUT_DIR + "/figures/admixture_cv_error.tiff",
    ])
    OPTIONAL_TARGETS.extend(expand(ADMIXTURE_OUTPUT_DIR + "/admixture_pruned.{k}.Q", k=ADMIXTURE_K_VALUES))

if LD_DECAY_CONFIG.get("enabled", False):
    LD_DECAY_OUTPUT_DIR = LD_DECAY_CONFIG.get("output_dir", "ld_decay")
    OPTIONAL_TARGETS.extend([
        LD_DECAY_OUTPUT_DIR + "/combined.ld_decay.tsv",
        LD_DECAY_OUTPUT_DIR + "/figures/ld_decay.pdf",
        LD_DECAY_OUTPUT_DIR + "/figures/ld_decay.png",
        LD_DECAY_OUTPUT_DIR + "/figures/ld_decay.svg",
        LD_DECAY_OUTPUT_DIR + "/figures/ld_decay.tiff",
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
        PCA_PLOT_PREFIX + ".C.3DPC1PC2PC3.png",
        PCA_PLOT_PREFIX + ".N.3DPC1PC2PC3.png",
        PCA_PLOT_PREFIX + ".C.3DPC1PC2PC3.svg",
        PCA_PLOT_PREFIX + ".N.3DPC1PC2PC3.svg",
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
        GVCF_TARGETS,
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
include: "rules/genomicsdb.rules"
include: "rules/faidx_index.rules"
include: "rules/haplo.rules"
include: "rules/indel_filter.rules"
include: "rules/indel_select.rules"
include: "rules/index_combined_vcf.rules"
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
include: "rules/cnv.rules"
include: "rules/pi.rules"
include: "rules/snp_density.rules"
include: "rules/admixture.rules"
include: "rules/ld_decay.rules"
include: "rules/vcf2pca.rules"
include: "rules/vcf2dis.rules"
include: "rules/snpeff.rules"
include: "rules/report.rules"
