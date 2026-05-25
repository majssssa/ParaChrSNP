#!/usr/bin/env python3
"""
Generate a unified ParaChrSNP HTML report.

修改说明：
1. 原始脚本路径：scripts/parachrsnp-report_1.py
2. 本版本在统一 HTML 报告中新增流程图片展示区域。
3. PNG/SVG 图片会直接嵌入报告，PDF 图会通过网页只读输出路由进行预览。
4. 保留上一版本的每个样本测序深度统计。

Author: Junpeng Ma 1527552938@qq.com
"""

import argparse
import base64
import csv
import gzip
import html
import json
import os
import re
import subprocess
from datetime import datetime

import yaml


AUTHOR = "Junpeng Ma 1527552938@qq.com"
SCRIPT_NAME = "parachrsnp-report_2.py"
SCRIPT_FUNCTION = "Generate ParaChrSNP HTML report with embedded workflow figures."


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a unified HTML report for the ParaChrSNP workflow."
    )
    parser.add_argument("--config", required=True, help="Path to the ParaChrSNP config.yaml file.")
    parser.add_argument("--out-html", required=True, help="Output HTML report path.")
    parser.add_argument("--out-tsv", required=True, help="Output summary TSV path.")
    return parser.parse_args()


def load_config(path):
    # 读取流程配置，用于确定样本、染色体和输出前缀。
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def write_script_log():
    # 按用户要求，把脚本名称、功能、脚本位置和创建/运行时间记录到固定日志文件。
    log_path = "/home/majunpeng/script_log.txt"
    line = (
        f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\t{SCRIPT_NAME}\t"
        f"{SCRIPT_FUNCTION}\t{os.path.abspath(__file__)}\n"
    )
    try:
        with open(log_path, "a", encoding="utf-8") as handle:
            handle.write(line)
    except OSError:
        pass


def fmt_int(value):
    if value is None:
        return "NA"
    return f"{int(value):,}"


def fmt_float(value, digits=4):
    if value is None:
        return "NA"
    return f"{float(value):.{digits}f}"


def file_size(path):
    # 将文件大小转换为适合报告展示的单位。
    if not os.path.exists(path):
        return "missing"
    size = os.path.getsize(path)
    units = ["B", "KB", "MB", "GB", "TB"]
    value = float(size)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{value:.1f} {unit}"
        value /= 1024


def count_vcf_records(path):
    # 统计 VCF 非 header 行数量，用于展示 raw/filter 后变异数量。
    if not os.path.exists(path):
        return None
    if path.endswith(".gz") and (os.path.exists(path + ".tbi") or os.path.exists(path + ".csi")):
        try:
            output = subprocess.check_output(
                ["bcftools", "index", "-n", path],
                text=True,
                stderr=subprocess.DEVNULL,
            ).strip()
            return int(output)
        except (OSError, subprocess.CalledProcessError, ValueError):
            pass
    opener = gzip.open if path.endswith(".gz") else open
    mode = "rt" if path.endswith(".gz") else "r"
    count = 0
    with opener(path, mode, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("#"):
                count += 1
    return count


def read_fastp_json(path):
    # 解析 fastp json，提取 reads 数、bases 数、Q30、GC 和过滤后比例。
    if not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    before_reads = before.get("total_reads")
    after_reads = after.get("total_reads")
    retention = (after_reads / before_reads * 100) if before_reads else None
    return {
        "before_reads": before_reads,
        "after_reads": after_reads,
        "before_bases": before.get("total_bases"),
        "after_bases": after.get("total_bases"),
        "retention": retention,
        "q30_rate": after.get("q30_rate"),
        "gc_content": after.get("gc_content"),
    }


def reference_genome_size(reference):
    # 优先读取 samtools faidx 生成的 .fai 文件，避免重复扫描大型参考基因组。
    fai = reference + ".fai"
    if os.path.exists(fai):
        total = 0
        with open(fai, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    continue
                try:
                    total += int(fields[1])
                except ValueError:
                    continue
        if total > 0:
            return total

    # 如果 .fai 不存在，则直接读取 FASTA 序列长度作为兜底。
    if not os.path.exists(reference):
        return None
    opener = gzip.open if reference.endswith(".gz") else open
    mode = "rt" if reference.endswith(".gz") else "r"
    total = 0
    with opener(reference, mode, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            total += len(line.strip())
    return total or None


def parse_dup_metrics(path):
    # 解析 GATK MarkDuplicates metrics，提取重复率。
    if not os.path.exists(path):
        return {}
    header = None
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if "PERCENT_DUPLICATION" in fields:
                header = fields
                continue
            if header and len(fields) == len(header):
                data = dict(zip(header, fields))
                try:
                    return {"duplication_rate": float(data.get("PERCENT_DUPLICATION", "nan"))}
                except ValueError:
                    return {}
    return {}


def parse_plink_imiss(path):
    # 解析 PLINK imiss，统计样本缺失率均值和最大值。
    if not os.path.exists(path):
        return {}
    values = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
        for row in reader:
            try:
                values.append(float(row["F_MISS"]))
            except (KeyError, ValueError):
                continue
    if not values:
        return {}
    return {"sample_missing_mean": sum(values) / len(values), "sample_missing_max": max(values)}


def parse_plink_lmiss(path):
    # 解析 PLINK lmiss，统计位点缺失率均值和最大值。
    if not os.path.exists(path):
        return {}
    values = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
        for row in reader:
            try:
                values.append(float(row["F_MISS"]))
            except (KeyError, ValueError):
                continue
    if not values:
        return {}
    return {"variant_missing_mean": sum(values) / len(values), "variant_missing_max": max(values)}


def parse_rastqc_summary(path):
    # 解析 RastQC summary.tsv，统计 pass/warn/fail 数。
    summary = {"pass": 0, "warn": 0, "fail": 0, "rows": 0}
    if not os.path.exists(path):
        return summary
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            lower = line.lower()
            if not line.strip():
                continue
            summary["rows"] += 1
            summary["pass"] += lower.count("pass")
            summary["warn"] += lower.count("warn")
            summary["fail"] += lower.count("fail")
    return summary


def parse_log_errors(log_dir):
    # 扫描日志中的常见错误关键词，用于报告末尾提示潜在风险。
    findings = []
    if not os.path.exists(log_dir):
        return findings
    pattern = re.compile(r"(error|exception|failed|cannot open|missingoutput)", re.IGNORECASE)
    for root, _, files in os.walk(log_dir):
        for name in files:
            if not name.endswith(".log"):
                continue
            path = os.path.join(root, name)
            try:
                with open(path, "r", encoding="utf-8", errors="replace") as handle:
                    for line_number, line in enumerate(handle, start=1):
                        if line_number > 200:
                            break
                        if pattern.search(line):
                            findings.append((path, line.strip()[:220]))
                            break
            except OSError:
                continue
    return findings[:30]


def html_table(headers, rows):
    head = "".join(f"<th>{html.escape(str(header))}</th>" for header in headers)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def output_url(path):
    # 网页服务会把 /outputs/<相对路径> 映射到项目目录中的只读文件。
    return "/outputs/" + path.replace("\\", "/").lstrip("/")


def data_uri(path):
    # 将 PNG/SVG 图片编码为 data URI，使报告本身可以直接显示关键图片。
    if not os.path.exists(path):
        return ""
    lower = path.lower()
    if lower.endswith(".png"):
        mime = "image/png"
        mode = "rb"
    elif lower.endswith(".svg"):
        mime = "image/svg+xml"
        mode = "rb"
    else:
        return ""
    with open(path, mode) as handle:
        payload = base64.b64encode(handle.read()).decode("ascii")
    return f"data:{mime};base64,{payload}"


def figure_card(title, path, description=""):
    # 根据文件类型生成报告中的图片卡片。不存在的文件自动跳过。
    if not path or not os.path.exists(path):
        return ""
    title_html = html.escape(title)
    path_html = html.escape(path)
    desc_html = f"<p>{html.escape(description)}</p>" if description else ""
    lower = path.lower()
    if lower.endswith((".png", ".svg")):
        src = data_uri(path)
        media = f'<img src="{src}" alt="{title_html}">'
    elif lower.endswith(".pdf"):
        media = f'<iframe src="{html.escape(output_url(path))}" title="{title_html}"></iframe>'
    elif lower.endswith(".html"):
        media = f'<iframe src="{html.escape(output_url(path))}" title="{title_html}"></iframe>'
    else:
        media = ""
    return (
        '<div class="figure-card">'
        f"<h3>{title_html}</h3>"
        f"{desc_html}"
        f"{media}"
        f'<p class="file-link"><a href="{html.escape(output_url(path))}" target="_blank" rel="noopener">{path_html}</a></p>'
        "</div>"
    )


def first_existing(paths):
    # 从同一结果的多个候选格式中选择第一个存在的文件。
    for path in paths:
        if path and os.path.exists(path):
            return path
    return ""


def collect_figure_cards(params, samples):
    # 收集流程中常见的图形结果，存在才展示，避免可选模块关闭时报错。
    cards = []
    missing_cfg = params.get("vcf_missing", {})
    missing_plot_prefix = missing_cfg.get(
        "plot_prefix",
        missing_cfg.get("output_prefix", "missing/combined.snp.filtered") + ".miss_check",
    )
    cards.extend([
        figure_card(
            "Sample missing rate",
            missing_plot_prefix + ".sample_missing.pdf",
            "Distribution of genotype missingness across samples.",
        ),
        figure_card(
            "Variant missing rate",
            missing_plot_prefix + ".Variant_missing.pdf",
            "Distribution of genotype missingness across variant sites.",
        ),
    ])

    pca_cfg = params.get("vcf2pca", {})
    if pca_cfg.get("enabled", True) and len(samples) >= 3:
        pca_prefix = pca_cfg.get("output_prefix", "pca/ParaChrSNP")
        pca_plot_prefix = pca_cfg.get("plot_prefix", pca_prefix + ".plot")
        pca_groups = [
            (
                "PCA PC1 vs PC2, colored",
                [
                    pca_plot_prefix + ".C.PC1_PC2.p.png",
                    pca_plot_prefix + ".C.PC1_PC2.p.svg",
                    pca_prefix + ".C.PC1_PC2.p.png",
                    pca_prefix + ".C.PC1_PC2.p.svg",
                ],
            ),
            (
                "PCA PC1 vs PC2, ungrouped",
                [
                    pca_plot_prefix + ".N.PC1_PC2.p.png",
                    pca_plot_prefix + ".N.PC1_PC2.p.svg",
                    pca_prefix + ".N.PC1_PC2.p.png",
                    pca_prefix + ".N.PC1_PC2.p.svg",
                ],
            ),
            (
                "PCA PC1 vs PC3, colored",
                [
                    pca_prefix + ".C.PC1_PC3.p.png",
                    pca_prefix + ".C.PC1_PC3.p.svg",
                ],
            ),
            (
                "PCA PC1 vs PC3, ungrouped",
                [
                    pca_prefix + ".N.PC1_PC3.p.png",
                    pca_prefix + ".N.PC1_PC3.p.svg",
                ],
            ),
            (
                "PCA 3D, colored",
                [
                    pca_plot_prefix + ".C.3DPC1PC2PC3.png",
                    pca_plot_prefix + ".C.3DPC1PC2PC3.svg",
                    pca_plot_prefix + ".C.3DPC1PC2PC3.pdf",
                    pca_prefix + ".C.3DPC1PC2PC3.png",
                    pca_prefix + ".C.3DPC1PC2PC3.svg",
                    pca_prefix + ".C.3DPC1PC2PC3.pdf",
                ],
            ),
            (
                "PCA 3D, ungrouped",
                [
                    pca_plot_prefix + ".N.3DPC1PC2PC3.png",
                    pca_plot_prefix + ".N.3DPC1PC2PC3.svg",
                    pca_plot_prefix + ".N.3DPC1PC2PC3.pdf",
                    pca_prefix + ".N.3DPC1PC2PC3.png",
                    pca_prefix + ".N.3DPC1PC2PC3.svg",
                    pca_prefix + ".N.3DPC1PC2PC3.pdf",
                ],
            ),
        ]
        for title, paths in pca_groups:
            path = first_existing(paths)
            cards.append(figure_card(title, path, "Principal component analysis based on filtered SNPs."))

    return "\n".join(card for card in cards if card)


def write_tsv(path, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("section\tmetric\tvalue\n")
        for section, metric, value in rows:
            handle.write(f"{section}\t{metric}\t{value}\n")


def main():
    write_script_log()
    args = parse_args()
    config = load_config(args.config)
    samples = config.get("samples", {})
    chromosomes = config.get("chromosomes", [])
    params = config.get("params", {})
    reference = config.get("reference")
    genome_size = reference_genome_size(reference) if reference else None
    rows_for_tsv = []

    rastqc = parse_rastqc_summary("qc/summary.tsv")
    rows_for_tsv.extend([
        ("qc", "rastqc_pass", rastqc["pass"]),
        ("qc", "rastqc_warn", rastqc["warn"]),
        ("qc", "rastqc_fail", rastqc["fail"]),
    ])

    fastp_rows = []
    depth_rows = []
    for sample in samples:
        stats = read_fastp_json(f"clean_reads/{sample}.json")
        dup = parse_dup_metrics(f"duplicate_removed/{sample}.dup.txt")
        after_bases = stats.get("after_bases")
        depth = (after_bases / genome_size) if after_bases and genome_size else None
        fastp_rows.append([
            sample,
            fmt_int(stats.get("before_reads")),
            fmt_int(stats.get("after_reads")),
            fmt_float(stats.get("retention"), 2),
            fmt_float(stats.get("q30_rate"), 4),
            fmt_float(stats.get("gc_content"), 4),
            fmt_float(dup.get("duplication_rate"), 4),
        ])
        depth_rows.append([
            sample,
            fmt_int(after_bases),
            fmt_int(genome_size),
            fmt_float(depth, 2),
        ])
        rows_for_tsv.extend([
            ("depth", f"{sample}_clean_bases", after_bases if after_bases is not None else "NA"),
            ("depth", f"{sample}_reference_genome_size", genome_size if genome_size is not None else "NA"),
            ("depth", f"{sample}_clean_depth", fmt_float(depth, 4)),
        ])

    vcf_paths = [
        ("Raw cohort VCF", "result_vcfs/combined.vcf.gz"),
        ("Filtered SNP VCF", "result_vcfs/combined.snp.filtered.vcf.gz"),
        ("Filtered INDEL VCF", "result_vcfs/combined.indel.filtered.vcf.gz"),
    ]
    vcf_rows = []
    for label, path in vcf_paths:
        count = count_vcf_records(path)
        vcf_rows.append([label, path, fmt_int(count), file_size(path)])
        rows_for_tsv.append(("variant", label, count if count is not None else "NA"))

    missing_prefix = params.get("vcf_missing", {}).get("output_prefix", "missing/combined.snp.filtered")
    imiss = parse_plink_imiss(missing_prefix + ".imiss")
    lmiss = parse_plink_lmiss(missing_prefix + ".lmiss")
    rows_for_tsv.extend([
        ("missing", "sample_missing_mean", imiss.get("sample_missing_mean", "NA")),
        ("missing", "sample_missing_max", imiss.get("sample_missing_max", "NA")),
        ("missing", "variant_missing_mean", lmiss.get("variant_missing_mean", "NA")),
        ("missing", "variant_missing_max", lmiss.get("variant_missing_max", "NA")),
    ])

    convert_prefix = params.get("vcf_convert", {}).get("output_prefix", "format_convert/combined.snp.filtered")
    output_rows = [
        ["PLINK binary BED", convert_prefix + ".bed", file_size(convert_prefix + ".bed")],
        ["PLINK binary BIM", convert_prefix + ".bim", file_size(convert_prefix + ".bim")],
        ["PLINK binary FAM", convert_prefix + ".fam", file_size(convert_prefix + ".fam")],
        ["PLINK text PED", convert_prefix + ".ped", file_size(convert_prefix + ".ped")],
        ["PLINK text MAP", convert_prefix + ".map", file_size(convert_prefix + ".map")],
        ["HapMap", convert_prefix + ".hmp.txt", file_size(convert_prefix + ".hmp.txt")],
    ]
    if params.get("imputation", {}).get("enabled", False):
        imputation_prefix = params.get("imputation", {}).get(
            "output_prefix",
            "imputation/combined.snp.filtered.beagle",
        )
        output_rows.extend([
            ["Beagle imputed VCF", imputation_prefix + ".vcf.gz", file_size(imputation_prefix + ".vcf.gz")],
            ["Beagle imputed VCF index", imputation_prefix + ".vcf.gz.tbi", file_size(imputation_prefix + ".vcf.gz.tbi")],
        ])
    if params.get("vcf2pca", {}).get("enabled", True) and len(samples) >= 3:
        pca_prefix = params.get("vcf2pca", {}).get("output_prefix", "pca/ParaChrSNP")
        pca_plot_prefix = params.get("vcf2pca", {}).get("plot_prefix", pca_prefix + ".plot")
        output_rows.extend([
            ["PCA eigenvec", pca_prefix + ".eigenvec", file_size(pca_prefix + ".eigenvec")],
            ["PCA eigenval", pca_prefix + ".eigenval", file_size(pca_prefix + ".eigenval")],
            ["PCA 2D SVG", pca_plot_prefix + ".C.PC1_PC2.p.svg", file_size(pca_plot_prefix + ".C.PC1_PC2.p.svg")],
            ["PCA 3D PDF", pca_plot_prefix + ".C.3DPC1PC2PC3.pdf", file_size(pca_plot_prefix + ".C.3DPC1PC2PC3.pdf")],
        ])
    if params.get("vcf2dis", {}).get("enabled", True) and len(samples) >= 3:
        output_rows.extend([
            ["Distance matrix", params.get("vcf2dis", {}).get("output_matrix", "dis/ParaChrSNP.p_dis.mat"), file_size(params.get("vcf2dis", {}).get("output_matrix", "dis/ParaChrSNP.p_dis.mat"))],
            ["NJ tree", params.get("vcf2dis", {}).get("output_tree", "dis/ParaChrSNP.p_dis.nwk"), file_size(params.get("vcf2dis", {}).get("output_tree", "dis/ParaChrSNP.p_dis.nwk"))],
        ])
    if params.get("snpeff", {}).get("enabled", False):
        snpeff_prefix = params.get("snpeff", {}).get("output_prefix", "annotation/combined")
        if params.get("snpeff", {}).get("annotate_snp", True):
            output_rows.extend([
                ["SnpEff SNP VCF", snpeff_prefix + ".snp.snpeff.vcf.gz", file_size(snpeff_prefix + ".snp.snpeff.vcf.gz")],
                ["SnpEff SNP report", snpeff_prefix + ".snp.snpeff.html", file_size(snpeff_prefix + ".snp.snpeff.html")],
            ])
        if params.get("snpeff", {}).get("annotate_indel", False):
            output_rows.extend([
                ["SnpEff INDEL VCF", snpeff_prefix + ".indel.snpeff.vcf.gz", file_size(snpeff_prefix + ".indel.snpeff.vcf.gz")],
                ["SnpEff INDEL report", snpeff_prefix + ".indel.snpeff.html", file_size(snpeff_prefix + ".indel.snpeff.html")],
            ])

    figure_gallery = collect_figure_cards(params, samples)
    if not figure_gallery:
        figure_gallery = '<p class="muted">No figure files were detected for the enabled modules.</p>'

    log_findings = parse_log_errors("logs")
    log_rows = [[path, message] for path, message in log_findings] or [["No error-like messages detected", "NA"]]

    write_tsv(args.out_tsv, rows_for_tsv)
    os.makedirs(os.path.dirname(args.out_html), exist_ok=True)
    html_doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>ParaChrSNP Workflow Report</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 36px; color: #17212b; }}
h1 {{ margin-bottom: 6px; }}
h2 {{ margin-top: 30px; border-bottom: 2px solid #dbe4ee; padding-bottom: 6px; }}
.cards {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 14px; margin: 22px 0; }}
.card {{ background: #f8fafc; border: 1px solid #cbd5e1; border-radius: 8px; padding: 14px 16px; }}
.card strong {{ display: block; font-size: 14px; color: #52616b; margin-bottom: 8px; }}
.card span {{ font-size: 26px; font-weight: 700; }}
table {{ border-collapse: collapse; width: 100%; font-size: 14px; margin-top: 12px; }}
th, td {{ border: 1px solid #dbe4ee; padding: 8px 10px; text-align: left; vertical-align: top; }}
th {{ background: #e0f2fe; }}
code {{ background: #f1f5f9; padding: 2px 4px; border-radius: 4px; }}
a {{ color: #0369a1; }}
.muted {{ color: #64748b; }}
.figure-grid {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 18px; margin-top: 14px; }}
.figure-card {{ border: 1px solid #cbd5e1; border-radius: 8px; padding: 14px; background: #f8fafc; }}
.figure-card h3 {{ margin: 0 0 8px; font-size: 16px; }}
.figure-card p {{ margin: 6px 0 10px; color: #475569; }}
.figure-card img {{ width: 100%; max-height: 520px; object-fit: contain; background: #fff; border: 1px solid #e2e8f0; border-radius: 6px; }}
.figure-card iframe {{ width: 100%; height: 520px; border: 1px solid #e2e8f0; border-radius: 6px; background: #fff; }}
.file-link {{ font-size: 13px; word-break: break-all; }}
@media (max-width: 980px) {{
  body {{ margin: 18px; }}
  .cards, .figure-grid {{ grid-template-columns: 1fr; }}
}}
</style>
</head>
<body>
<h1>ParaChrSNP Workflow Report</h1>
<p>Author: {html.escape(AUTHOR)}</p>
<p>Generated at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
<p>Config file: <code>{html.escape(args.config)}</code></p>
<div class="cards">
  <div class="card"><strong>Samples</strong><span>{len(samples)}</span></div>
  <div class="card"><strong>Chromosomes</strong><span>{len(chromosomes)}</span></div>
  <div class="card"><strong>RastQC warnings</strong><span>{rastqc["warn"]}</span></div>
  <div class="card"><strong>RastQC failures</strong><span>{rastqc["fail"]}</span></div>
</div>

<h2>Read Processing and Duplicate Metrics</h2>
{html_table(["Sample", "Raw reads", "Clean reads", "Retention (%)", "Clean Q30 rate", "Clean GC content", "Duplicate rate"], fastp_rows)}

<h2>Sequencing Depth</h2>
<p>Depth is calculated as clean bases divided by reference genome size. Reference genome size is read from <code>{html.escape(reference + ".fai" if reference else "NA")}</code> when available.</p>
{html_table(["Sample", "Clean bases", "Reference genome size", "Clean sequencing depth (x)"], depth_rows)}

<h2>Variant Summary</h2>
{html_table(["Dataset", "File", "Variant records", "File size"], vcf_rows)}

<h2>Missingness Summary</h2>
{html_table(["Metric", "Value"], [
    ["Mean sample missing rate", fmt_float(imiss.get("sample_missing_mean"), 4)],
    ["Max sample missing rate", fmt_float(imiss.get("sample_missing_max"), 4)],
    ["Mean variant missing rate", fmt_float(lmiss.get("variant_missing_mean"), 4)],
    ["Max variant missing rate", fmt_float(lmiss.get("variant_missing_max"), 4)],
])}
<p>Sample plot: <a href="../{html.escape(missing_prefix + '.miss_check.sample_missing.pdf')}">{html.escape(missing_prefix + '.miss_check.sample_missing.pdf')}</a></p>
<p>Variant plot: <a href="../{html.escape(missing_prefix + '.miss_check.Variant_missing.pdf')}">{html.escape(missing_prefix + '.miss_check.Variant_missing.pdf')}</a></p>

<h2>Workflow Figures</h2>
<div class="figure-grid">
{figure_gallery}
</div>

<h2>Downstream Output Files</h2>
{html_table(["Output", "File", "File size"], output_rows)}

<h2>Precheck Report</h2>
<p>Precheck HTML: <a href="precheck.html">reports/precheck.html</a></p>
<p>Precheck TSV: <a href="precheck.tsv">reports/precheck.tsv</a></p>

<h2>Log Warnings and Error-like Messages</h2>
{html_table(["Log file", "First matching message"], log_rows)}
</body>
</html>
"""
    with open(args.out_html, "w", encoding="utf-8") as handle:
        handle.write(html_doc)


if __name__ == "__main__":
    main()
