#!/usr/bin/env python3
"""
ParaChrSNP preflight validation.

修改要求：增加 Types=reads/bam；BAM 模式检查已去重 BAM 的结构、排序、样本名和参考染色体。
原始脚本位置：scripts/precheck.py

Author: Junpeng Ma 1527552938@qq.com
"""

import argparse
import gzip
import html
import os
import shutil
import subprocess
import sys
from datetime import datetime

import yaml


AUTHOR = "Junpeng Ma 1527552938@qq.com"


def write_script_log():
    # 记录新版本脚本名称、功能、位置和执行时间。
    try:
        with open("/home/majunpeng/script_log.txt", "a", encoding="utf-8") as handle:
            handle.write(
                f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\tprecheck_1.py\t"
                f"Validate reads or deduplicated BAM input\t{os.path.realpath(__file__)}\n"
            )
    except OSError:
        # 用户在其他机器安装时没有该个人目录，不应阻止分析。
        pass


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

    input_type = str(config.get("Types", "reads")).lower()
    if input_type == "bam":
        check_bam_samples(config, results)
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


def check_bam_samples(config, results):
    # BAM 模式按样本键推导文件名；samples 的值不再作为 FASTQ 前缀使用。
    bam_dir = str(config.get("bam_dir", "input_bam"))
    chromosomes = set(config.get("chromosomes") or [])
    # 已有 FASTA 索引时，额外核对染色体长度，减少混用参考版本的风险。
    reference_lengths = {}
    fai = str(config.get("reference", "")) + ".fai"
    if os.path.isfile(fai):
        with open(fai, "r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.split("\t")
                if len(fields) >= 2 and fields[0] in chromosomes:
                    reference_lengths[fields[0]] = int(fields[1])
    if os.path.abspath(bam_dir) == os.path.abspath("staged_bam"):
        add_result(results, "ERROR", "bam_dir", "failed", "bam_dir cannot equal staged_bam.")
        return
    if not shutil.which("samtools"):
        add_result(results, "ERROR", "samtools", "failed", "samtools is required to validate BAM files.")
        return
    for sample in config["samples"]:
        # 不允许路径分隔符或通配符混入样本名，避免解析到意外文件。
        if not isinstance(sample, str) or not sample or not all(c.isalnum() or c in "._-" for c in sample):
            add_result(results, "ERROR", f"sample:{sample}", "failed", "Sample ID may only contain letters, digits, dots, underscores or hyphens.")
            continue
        path = os.path.join(bam_dir, sample + ".bam")
        exists, readable, size, realpath = file_state(path)
        if not exists or not readable or size == 0:
            add_result(results, "ERROR", f"{sample}:BAM", "failed", f"BAM is missing, unreadable or empty: {path}")
            continue
        quickcheck = subprocess.run(["samtools", "quickcheck", "-v", path], capture_output=True, text=True)
        if quickcheck.returncode:
            add_result(results, "ERROR", f"{sample}:BAM", "failed", f"BAM failed samtools quickcheck: {quickcheck.stderr.strip()}")
            continue
        header = subprocess.run(["samtools", "view", "-H", path], capture_output=True, text=True)
        if header.returncode:
            add_result(results, "ERROR", f"{sample}:BAM", "failed", f"Cannot read BAM header: {header.stderr.strip()}")
            continue
        lines = header.stdout.splitlines()
        hd = next((line for line in lines if line.startswith("@HD\t")), "")
        if "SO:coordinate" not in hd.split("\t"):
            add_result(results, "ERROR", f"{sample}:sort_order", "failed", "BAM header must declare SO:coordinate.")
        else:
            add_result(results, "INFO", f"{sample}:sort_order", "passed", "BAM declares coordinate sorting.")
        contigs = {field[3:] for line in lines if line.startswith("@SQ\t") for field in line.split("\t") if field.startswith("SN:")}
        missing = sorted(chromosomes - contigs)
        if missing:
            add_result(results, "ERROR", f"{sample}:reference", "failed", "BAM header lacks configured chromosome(s): " + ", ".join(missing))
        else:
            add_result(results, "INFO", f"{sample}:reference", "passed", "All configured chromosomes occur in BAM header.")
        bam_lengths = {}
        for line in lines:
            if line.startswith("@SQ\t"):
                fields = dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
                if "SN" in fields and "LN" in fields:
                    bam_lengths[fields["SN"]] = int(fields["LN"])
        oversized = [name for name, length in bam_lengths.items() if length > 2**29]
        if oversized:
            add_result(results, "ERROR", f"{sample}:bai_limit", "failed", "BAI indexing cannot cover contigs longer than 512 MiB: " + ", ".join(oversized))
        wrong_lengths = [name for name, length in reference_lengths.items() if bam_lengths.get(name) not in (None, length)]
        if wrong_lengths:
            add_result(results, "ERROR", f"{sample}:reference_length", "failed", "BAM @SQ length differs from FASTA index: " + ", ".join(wrong_lengths))
        elif reference_lengths:
            add_result(results, "INFO", f"{sample}:reference_length", "passed", "Configured chromosome lengths match FASTA index.")
        names = {field[3:] for line in lines if line.startswith("@RG\t") for field in line.split("\t") if field.startswith("SM:")}
        if names != {sample}:
            add_result(results, "ERROR", f"{sample}:read_group", "failed", f"BAM @RG SM must equal {sample}; found {sorted(names)}")
        else:
            add_result(results, "INFO", f"{sample}:read_group", "passed", f"BAM read-group sample is {sample}.")
        add_result(results, "INFO", f"{sample}:BAM", "passed", f"BAM is readable: {path}; real path: {realpath}")
    add_result(results, "WARNING", "deduplication", "warning", "BAM input is assumed to be duplicate-removed; this cannot be proven from its header.")
    if len(config["samples"]) < 3:
        add_result(results, "WARNING", "sample_count", "warning", "Fewer than 3 samples configured; PCA and distance/tree outputs are skipped.")


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

    pi_cfg = params.get("pi", {})
    if pi_cfg.get("enabled", False):
        pop_info = pi_cfg.get("pop_info")
        if not pop_info:
            add_result(results, "INFO", "pi.pop_info", "passed", "No pop_info configured; Pi will be calculated for all samples together.")
        elif not os.path.exists(pop_info):
            add_result(results, "ERROR", "pi.pop_info", "failed", f"Configured Pi pop_info file does not exist: {pop_info}")
        else:
            group_samples = []
            with open(pop_info, "r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    group_samples.append(line.split()[0])
            unknown = sorted(set(group_samples) - samples)
            if unknown:
                add_result(results, "ERROR", "pi.pop_info", "failed", "Sample(s) in Pi pop_info file not found in config samples: " + ", ".join(unknown))
            else:
                add_result(results, "INFO", "pi.pop_info", "passed", f"Pi pop_info file is valid: {pop_info}")

    admixture_cfg = params.get("admixture", {})
    if admixture_cfg.get("enabled", False):
        executable = admixture_cfg.get("executable", "admixture")
        if shutil.which(executable) or (os.path.exists(executable) and os.access(executable, os.X_OK)):
            add_result(results, "INFO", "admixture.executable", "passed", f"ADMIXTURE executable is available: {executable}")
        else:
            add_result(results, "ERROR", "admixture.executable", "failed", f"ADMIXTURE executable was not found: {executable}")

        pop_info = admixture_cfg.get("pop_info")
        if not pop_info:
            add_result(results, "INFO", "admixture.pop_info", "passed", "No pop_info configured; samples will be plotted by PLINK FAM order.")
        elif not os.path.exists(pop_info):
            add_result(results, "ERROR", "admixture.pop_info", "failed", f"Configured ADMIXTURE pop_info file does not exist: {pop_info}")
        else:
            group_samples = []
            with open(pop_info, "r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    group_samples.append(line.split()[0])
            unknown = sorted(set(group_samples) - samples)
            if unknown:
                add_result(results, "ERROR", "admixture.pop_info", "failed", "Sample(s) in ADMIXTURE pop_info file not found in config samples: " + ", ".join(unknown))
            else:
                add_result(results, "INFO", "admixture.pop_info", "passed", f"ADMIXTURE pop_info file is valid: {pop_info}")

        try:
            k_min = int(admixture_cfg.get("k_min", 1))
            k_max = int(admixture_cfg.get("k_max", 10))
            if k_min < 1 or k_max < k_min:
                add_result(results, "ERROR", "admixture.K", "failed", "ADMIXTURE requires k_min >= 1 and k_max >= k_min.")
            elif k_max > len(samples):
                add_result(results, "WARNING", "admixture.K", "warning", f"k_max={k_max} is greater than the sample count ({len(samples)}).")
            else:
                add_result(results, "INFO", "admixture.K", "passed", f"ADMIXTURE K range: {k_min}-{k_max}")
        except ValueError:
            add_result(results, "ERROR", "admixture.K", "failed", "ADMIXTURE k_min and k_max must be integers.")

    ld_decay_cfg = params.get("ld_decay", {})
    if ld_decay_cfg.get("enabled", False):
        executable = ld_decay_cfg.get("executable", "PopLDdecay")
        if shutil.which(executable) or (os.path.exists(executable) and os.access(executable, os.X_OK)):
            add_result(results, "INFO", "ld_decay.executable", "passed", f"PopLDdecay executable is available: {executable}")
        else:
            add_result(results, "ERROR", "ld_decay.executable", "failed", f"PopLDdecay executable was not found: {executable}")

        pop_info = ld_decay_cfg.get("pop_info")
        if not pop_info:
            add_result(results, "INFO", "ld_decay.pop_info", "passed", "No pop_info configured; LD decay will be calculated for all samples together.")
        elif not os.path.exists(pop_info):
            add_result(results, "ERROR", "ld_decay.pop_info", "failed", f"Configured LD decay pop_info file does not exist: {pop_info}")
        else:
            group_samples = []
            with open(pop_info, "r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    group_samples.append(line.split()[0])
            unknown = sorted(set(group_samples) - samples)
            if unknown:
                add_result(results, "ERROR", "ld_decay.pop_info", "failed", "Sample(s) in LD decay pop_info file not found in config samples: " + ", ".join(unknown))
            else:
                add_result(results, "INFO", "ld_decay.pop_info", "passed", f"LD decay pop_info file is valid: {pop_info}")

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
    if os.path.isdir(image) and os.access(image, os.R_OK | os.X_OK):
        add_result(results, "INFO", "container.image", "passed", f"Container sandbox is readable: {image}")
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
    write_script_log()
    args = parse_args()
    results = []
    try:
        config = load_config(args.config)
        input_type = str(config.get("Types", "reads")).lower()
        if input_type not in ("reads", "bam"):
            raise ValueError("Types must be 'reads' or 'bam'.")
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
