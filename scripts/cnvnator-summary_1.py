#!/usr/bin/env python3
"""
Parse, summarize, and plot unfiltered CNVnator results for ParaChrSNP.

修改说明：
本脚本基于 scripts/cnvnator-summary.py 生成，用于满足 CNV 模块不再过滤
CNVnator 原始 call 的需求。原始脚本位置：scripts/cnvnator-summary.py。

Author: Junpeng Ma 1527552938@qq.com
"""

import argparse
import csv
import os
import re
from collections import Counter, defaultdict
from datetime import datetime


AUTHOR = "Junpeng Ma 1527552938@qq.com"
SCRIPT_NAME = "cnvnator-summary_1.py"
SCRIPT_FUNCTION = "Normalize unfiltered CNVnator calls and merge sample CNV tables."

def write_script_log():
    # 按用户要求，把脚本名称、功能、脚本位置和运行时间记录到固定日志文件。
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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Normalize and summarize unfiltered CNVnator CNV calls for ParaChrSNP."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    filter_parser = subparsers.add_parser(
        "filter",
        help="Normalize one CNVnator raw call file without filtering and write a TSV table.",
    )
    filter_parser.add_argument("--sample", required=True, help="Sample name.")
    filter_parser.add_argument("--input", required=True, help="Input CNVnator raw call text file.")
    filter_parser.add_argument("--output", required=True, help="Output normalized CNV TSV file.")
    filter_parser.add_argument("--min-len", type=int, default=None, help=argparse.SUPPRESS)
    filter_parser.add_argument("--max-evalue", type=float, default=None, help=argparse.SUPPRESS)
    filter_parser.add_argument("--max-q0", type=float, default=None, help=argparse.SUPPRESS)

    summary_parser = subparsers.add_parser(
        "summary",
        help="Merge filtered CNV TSV files and draw summary plots.",
    )
    summary_parser.add_argument("--inputs", nargs="+", required=True, help="Per-sample CNV TSV files.")
    summary_parser.add_argument("--combined", required=True, help="Output combined CNV TSV file.")
    summary_parser.add_argument("--summary", required=True, help="Output CNV summary TSV file.")
    summary_parser.add_argument("--count-plot", required=True, help="Output CNV count by sample PDF.")
    summary_parser.add_argument("--length-plot", required=True, help="Output CNV length distribution PDF.")
    summary_parser.add_argument("--chrom-plot", required=True, help="Output CNV count by chromosome PDF.")
    return parser.parse_args()


def parse_region(region):
    # 解析 CNVnator 的 chr:start-end 坐标格式，返回染色体、起点和终点。
    match = re.match(r"^(.+):(\d+)-(\d+)$", region)
    if not match:
        return region, "", ""
    chrom, start, end = match.groups()
    return chrom, int(start), int(end)


def parse_float(value):
    # 兼容科学计数法和 NA，无法解析时返回 None。
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_cnvnator_line(line, sample):
    # CNVnator 原始输出一般为：type region length normalized_rd e1 e2 e3 e4 q0。
    fields = line.strip().split()
    if len(fields) < 4:
        return None
    cnv_type = fields[0]
    region = fields[1]
    chrom, start, end = parse_region(region)
    try:
        length = int(float(fields[2]))
    except ValueError:
        if isinstance(start, int) and isinstance(end, int):
            length = end - start + 1
        else:
            length = 0
    normalized_rd = parse_float(fields[3])
    evalues = [parse_float(value) for value in fields[4:8]]
    q0 = parse_float(fields[8]) if len(fields) > 8 else None
    valid_evalues = [value for value in evalues if value is not None]
    min_evalue = min(valid_evalues) if valid_evalues else None
    return {
        "sample": sample,
        "type": cnv_type,
        "chrom": chrom,
        "start": start,
        "end": end,
        "length": length,
        "normalized_rd": normalized_rd,
        "e_val1": evalues[0] if len(evalues) > 0 else None,
        "e_val2": evalues[1] if len(evalues) > 1 else None,
        "e_val3": evalues[2] if len(evalues) > 2 else None,
        "e_val4": evalues[3] if len(evalues) > 3 else None,
        "min_evalue": min_evalue,
        "q0": q0,
    }


def keep_cnv_type(row):
    # 不根据长度、显著性或 q0 过滤，只保留 CNVnator 识别出的 deletion/duplication 类型。
    if row["type"] not in {"deletion", "duplication"}:
        return False
    return True


def format_value(value):
    # 写出 TSV 时将缺失值统一写为 NA。
    if value is None or value == "":
        return "NA"
    return str(value)


def write_rows(path, rows):
    # 写出标准化 CNV 表格，方便后续统计和作图。
    os.makedirs(os.path.dirname(path), exist_ok=True)
    header = [
        "sample",
        "type",
        "chrom",
        "start",
        "end",
        "length",
        "normalized_rd",
        "e_val1",
        "e_val2",
        "e_val3",
        "e_val4",
        "min_evalue",
        "q0",
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_value(row.get(key)) for key in header})


def read_filtered_table(path):
    # 读取标准化后的 per-sample CNV TSV，用于合并和统计。
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                row["length"] = int(float(row.get("length", 0)))
            except ValueError:
                row["length"] = 0
            rows.append(row)
    return rows


def filter_one(args):
    # 为了兼容已有 rule 名称，命令仍叫 filter；实际行为是不做任何阈值过滤。
    rows = []
    with open(args.input, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            row = parse_cnvnator_line(line, args.sample)
            if row and keep_cnv_type(row):
                rows.append(row)
    write_rows(args.output, rows)


def summarize(args):
    all_rows = []
    for path in args.inputs:
        all_rows.extend(read_filtered_table(path))
    write_rows(args.combined, all_rows)

    sample_counts = Counter(row["sample"] for row in all_rows)
    chrom_counts = Counter(row["chrom"] for row in all_rows)
    type_counts = Counter(row["type"] for row in all_rows)
    sample_lengths = defaultdict(int)
    for row in all_rows:
        sample_lengths[row["sample"]] += row["length"]

    os.makedirs(os.path.dirname(args.summary), exist_ok=True)
    with open(args.summary, "w", encoding="utf-8") as handle:
        handle.write("section\tname\tcnv_count\ttotal_length_bp\n")
        for sample in sorted(sample_counts):
            handle.write(f"sample\t{sample}\t{sample_counts[sample]}\t{sample_lengths[sample]}\n")
        for cnv_type in sorted(type_counts):
            handle.write(f"type\t{cnv_type}\t{type_counts[cnv_type]}\tNA\n")
        for chrom in sorted(chrom_counts):
            handle.write(f"chromosome\t{chrom}\t{chrom_counts[chrom]}\tNA\n")

    # 绘图由 scripts/cnvnator-plot-nature_1.R 完成，保证项目自绘图统一使用 R 后端。


def main():
    write_script_log()
    args = parse_args()
    if args.command == "filter":
        filter_one(args)
    elif args.command == "summary":
        summarize(args)


if __name__ == "__main__":
    main()
