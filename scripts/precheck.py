#!/usr/bin/env python3
"""
ParaChrSNP preflight validation.

Author: Junpeng Ma 1527552938@qq.com
"""

import argparse
import gzip
import html
import os
import sys
from datetime import datetime

import yaml


AUTHOR = "Junpeng Ma 1527552938@qq.com"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate ParaChrSNP input files and configuration before running the workflow."
    )
    parser.add_argument("--config", required=True, help="Path to the ParaChrSNP config.yaml file.")
    parser.add_argument("--out-tsv", required=True, help="Output TSV file recording all checks.")
    parser.add_argument("--out-html", required=True, help="Output HTML report for precheck results.")
    parser.add_argument("--done", required=True, help="Output done flag created only when no fatal error is found.")
    return parser.parse_args()


def add_result(results, level, item, status, message):
    # 统一记录检查结果，level 分为 ERROR/WARNING/INFO。
    results.append(
        {
            "level": level,
            "item": item,
            "status": status,
            "message": message,
        }
    )


def load_config(path):
    # 读取 YAML 配置文件，并检查是否是字典结构。
    with open(path, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError("config.yaml is empty or not a YAML mapping.")
    return config


def file_state(path):
    # 返回文件是否存在、是否可读、大小和软链接真实路径。
    exists = os.path.exists(path)
    readable = os.access(path, os.R_OK) if exists else False
    size = os.path.getsize(path) if exists and os.path.isfile(path) else 0
    realpath = os.path.realpath(path)
    return exists, readable, size, realpath


def is_gzip_file(path):
    # 通过 gzip magic number 快速判断压缩 FASTQ 是否像 gzip 文件。
    try:
        with open(path, "rb") as handle:
            return handle.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def fasta_headers(path, max_records=None):
    # 读取 FASTA header，用于检查 config 中的染色体名是否存在于参考基因组。
    headers = []
    opener = gzip.open if path.endswith(".gz") else open
    mode = "rt" if path.endswith(".gz") else "r"
    with opener(path, mode, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                headers.append(line[1:].strip().split()[0])
                if max_records and len(headers) >= max_records:
                    break
    return headers


def check_reference(config, results):
    reference = config.get("reference")
    if not reference:
        add_result(results, "ERROR", "reference", "failed", "Missing required key: reference")
        return

    exists, readable, size, realpath = file_state(reference)
    if not exists:
        add_result(results, "ERROR", "reference", "failed", f"Reference FASTA does not exist: {reference}")
        return
    if not readable:
        add_result(results, "ERROR", "reference", "failed", f"Reference FASTA is not readable: {reference}")
        return
    if size == 0:
        add_result(results, "ERROR", "reference", "failed", f"Reference FASTA is empty: {reference}")
        return

    if os.path.islink(reference):
        add_result(results, "INFO", "reference", "passed", f"Reference is a symlink to {realpath}")
    else:
        add_result(results, "INFO", "reference", "passed", f"Reference FASTA is readable: {reference}")

    chromosomes = config.get("chromosomes", [])
    if not chromosomes:
        add_result(results, "ERROR", "chromosomes", "failed", "No chromosome names are configured.")
        return

    try:
        headers = set(fasta_headers(reference))
    except Exception as exc:
        add_result(results, "ERROR", "reference", "failed", f"Failed to read FASTA headers: {exc}")
        return

    missing = [chrom for chrom in chromosomes if chrom not in headers]
    if missing:
        add_result(
            results,
            "ERROR",
            "chromosomes",
            "failed",
            "Chromosome(s) not found in reference FASTA: " + ", ".join(missing),
        )
    else:
        add_result(results, "INFO", "chromosomes", "passed", f"All {len(chromosomes)} configured chromosomes exist in reference.")


def check_samples(config, results):
    samples = config.get("samples")
    if not isinstance(samples, dict) or not samples:
        add_result(results, "ERROR", "samples", "failed", "Missing or empty samples mapping.")
        return

    seen_prefixes = {}
    for sample, prefix in samples.items():
        if not sample or any(char.isspace() for char in str(sample)):
            add_result(results, "ERROR", f"sample:{sample}", "failed", "Sample name is empty or contains whitespace.")
        if prefix in seen_prefixes:
            add_result(results, "ERROR", f"sample:{sample}", "failed", f"FASTQ prefix duplicated with {seen_prefixes[prefix]}: {prefix}")
        seen_prefixes[prefix] = sample

        for read in ("1", "2"):
            fq = f"{prefix}.{read}.fq.gz"
            exists, readable, size, realpath = file_state(fq)
            item = f"{sample}:R{read}"
            if not exists:
                add_result(results, "ERROR", item, "failed", f"FASTQ file does not exist: {fq}")
                continue
            if not readable:
                add_result(results, "ERROR", item, "failed", f"FASTQ file is not readable: {fq}")
                continue
            if size == 0:
                add_result(results, "ERROR", item, "failed", f"FASTQ file is empty: {fq}")
                continue
            if not is_gzip_file(fq):
                add_result(results, "ERROR", item, "failed", f"FASTQ file is not gzip-compressed or has an invalid gzip header: {fq}")
                continue
            if os.path.islink(fq):
                add_result(results, "INFO", item, "passed", f"FASTQ is readable; symlink target: {realpath}")
            else:
                add_result(results, "INFO", item, "passed", f"FASTQ is readable: {fq}")

    if len(samples) < 3:
        add_result(
            results,
            "WARNING",
            "sample_count",
            "warning",
            "Fewer than 3 samples configured. PCA and distance/tree modules will be skipped by rule all.",
        )
    else:
        add_result(results, "INFO", "sample_count", "passed", f"{len(samples)} samples configured.")


def check_optional_files(config, results):
    # 检查 pop.info 等可选分组文件是否存在，并检查样本名是否能匹配 config。
    samples = set((config.get("samples") or {}).keys())
    params = config.get("params", {})
    for module in ("vcf2pca", "vcf2dis"):
        module_cfg = params.get(module, {})
        sample_group = module_cfg.get("sample_group")
        if not sample_group:
            add_result(results, "INFO", f"{module}.sample_group", "passed", "No sample group file configured; group argument will be omitted.")
            continue
        if not os.path.exists(sample_group):
            add_result(results, "ERROR", f"{module}.sample_group", "failed", f"Configured sample group file does not exist: {sample_group}")
            continue
        group_samples = []
        with open(sample_group, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                group_samples.append(line.split()[0])
        unknown = sorted(set(group_samples) - samples)
        if unknown:
            add_result(results, "ERROR", f"{module}.sample_group", "failed", "Sample(s) in group file not found in config samples: " + ", ".join(unknown))
        else:
            add_result(results, "INFO", f"{module}.sample_group", "passed", f"Sample group file is valid: {sample_group}")

    snpeff_cfg = params.get("snpeff", {})
    if snpeff_cfg.get("enabled", False):
        genome_fasta = snpeff_cfg.get("genome_fasta")
        annotation_file = snpeff_cfg.get("annotation_file")
        for key, path in (("genome_fasta", genome_fasta), ("annotation_file", annotation_file)):
            if not path:
                add_result(results, "ERROR", f"snpeff.{key}", "failed", f"Missing required SnpEff config key: {key}")
                continue
            exists, readable, size, _ = file_state(path)
            if not exists:
                add_result(results, "ERROR", f"snpeff.{key}", "failed", f"SnpEff input file does not exist: {path}")
            elif not readable:
                add_result(results, "ERROR", f"snpeff.{key}", "failed", f"SnpEff input file is not readable: {path}")
            elif size == 0:
                add_result(results, "ERROR", f"snpeff.{key}", "failed", f"SnpEff input file is empty: {path}")
            else:
                add_result(results, "INFO", f"snpeff.{key}", "passed", f"SnpEff input file is readable: {path}")

        annotation_format = snpeff_cfg.get("annotation_format", "gff3")
        if annotation_format not in ("gff3", "gtf"):
            add_result(results, "ERROR", "snpeff.annotation_format", "failed", "SnpEff annotation_format must be gff3 or gtf.")
        else:
            add_result(results, "INFO", "snpeff.annotation_format", "passed", f"SnpEff annotation format: {annotation_format}")
    else:
        add_result(results, "INFO", "snpeff", "passed", "SnpEff annotation module is disabled.")


def check_container(config, results):
    image = (config.get("container") or {}).get("image")
    if not image:
        add_result(results, "WARNING", "container.image", "warning", "No container image is configured.")
        return
    exists, readable, size, _ = file_state(image)
    if exists and readable and size > 0:
        add_result(results, "INFO", "container.image", "passed", f"Container image is readable: {image}")
    else:
        add_result(results, "WARNING", "container.image", "warning", f"Container image is not readable from current path: {image}")


def write_tsv(results, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("level\titem\tstatus\tmessage\n")
        for row in results:
            handle.write(f"{row['level']}\t{row['item']}\t{row['status']}\t{row['message']}\n")


def write_html(results, path, config_path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    error_count = sum(1 for row in results if row["level"] == "ERROR")
    warning_count = sum(1 for row in results if row["level"] == "WARNING")
    status = "PASSED" if error_count == 0 else "FAILED"
    rows = []
    for row in results:
        css = row["level"].lower()
        rows.append(
            "<tr class='{css}'><td>{level}</td><td>{item}</td><td>{status}</td><td>{message}</td></tr>".format(
                css=css,
                level=html.escape(row["level"]),
                item=html.escape(row["item"]),
                status=html.escape(row["status"]),
                message=html.escape(row["message"]),
            )
        )
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>ParaChrSNP precheck report</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 36px; color: #17212b; }}
h1 {{ margin-bottom: 4px; }}
.summary {{ display: flex; gap: 16px; margin: 24px 0; }}
.card {{ border: 1px solid #cbd5e1; border-radius: 8px; padding: 14px 18px; background: #f8fafc; }}
table {{ border-collapse: collapse; width: 100%; font-size: 14px; }}
th, td {{ border: 1px solid #dbe4ee; padding: 8px 10px; text-align: left; vertical-align: top; }}
th {{ background: #e0f2fe; }}
tr.error td {{ background: #fee2e2; }}
tr.warning td {{ background: #fef9c3; }}
tr.info td {{ background: #f8fafc; }}
</style>
</head>
<body>
<h1>ParaChrSNP Precheck Report</h1>
<p>Author: {html.escape(AUTHOR)}</p>
<p>Config file: <code>{html.escape(config_path)}</code></p>
<p>Generated at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
<div class="summary">
  <div class="card"><strong>Status</strong><br>{status}</div>
  <div class="card"><strong>Errors</strong><br>{error_count}</div>
  <div class="card"><strong>Warnings</strong><br>{warning_count}</div>
  <div class="card"><strong>Total checks</strong><br>{len(results)}</div>
</div>
<table>
<thead><tr><th>Level</th><th>Item</th><th>Status</th><th>Message</th></tr></thead>
<tbody>
{''.join(rows)}
</tbody>
</table>
</body>
</html>
"""
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(document)


def main():
    args = parse_args()
    results = []
    try:
        config = load_config(args.config)
        check_reference(config, results)
        check_samples(config, results)
        check_optional_files(config, results)
        check_container(config, results)
    except Exception as exc:
        add_result(results, "ERROR", "precheck", "failed", str(exc))

    write_tsv(results, args.out_tsv)
    write_html(results, args.out_html, args.config)

    error_count = sum(1 for row in results if row["level"] == "ERROR")
    if error_count:
        if os.path.exists(args.done):
            os.remove(args.done)
        sys.exit(1)

    os.makedirs(os.path.dirname(args.done), exist_ok=True)
    with open(args.done, "w", encoding="utf-8") as handle:
        handle.write(datetime.now().isoformat() + "\n")


if __name__ == "__main__":
    main()
