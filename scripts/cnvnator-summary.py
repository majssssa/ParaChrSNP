#!/usr/bin/env python3
"""
Parse, filter, summarize, and plot CNVnator results for ParaChrSNP.

Author: Junpeng Ma 1527552938@qq.com
"""

import argparse
import csv
import os
import re
from collections import Counter, defaultdict
from datetime import datetime


AUTHOR = "Junpeng Ma 1527552938@qq.com"
SCRIPT_NAME = "cnvnator-summary.py"
SCRIPT_FUNCTION = "Filter CNVnator calls, merge sample CNV tables, and draw CNV summary plots."

MPLCONFIGDIR = os.path.join(os.environ.get("TMPDIR", "/tmp"), "parachrsnp_matplotlib")
os.environ.setdefault("MPLCONFIGDIR", MPLCONFIGDIR)
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)


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
        description="Filter and summarize CNVnator CNV calls for ParaChrSNP."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    filter_parser = subparsers.add_parser(
        "filter",
        help="Filter one CNVnator raw call file and write a normalized TSV table.",
    )
    filter_parser.add_argument("--sample", required=True, help="Sample name.")
    filter_parser.add_argument("--input", required=True, help="Input CNVnator raw call text file.")
    filter_parser.add_argument("--output", required=True, help="Output filtered CNV TSV file.")
    filter_parser.add_argument("--min-len", type=int, default=1000, help="Minimum CNV length in bp.")
    filter_parser.add_argument("--max-evalue", type=float, default=0.01, help="Maximum e-value threshold.")
    filter_parser.add_argument("--max-q0", type=float, default=0.5, help="Maximum q0 threshold.")

    summary_parser = subparsers.add_parser(
        "summary",
        help="Merge filtered CNV TSV files and draw summary plots.",
    )
    summary_parser.add_argument("--inputs", nargs="+", required=True, help="Filtered per-sample CNV TSV files.")
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


def passes_filter(row, min_len, max_evalue, max_q0):
    # 根据长度、显著性和 q0 过滤 CNV，保留较可信的 deletion/duplication。
    if row["type"] not in {"deletion", "duplication"}:
        return False
    if row["length"] < min_len:
        return False
    if row["min_evalue"] is not None and row["min_evalue"] > max_evalue:
        return False
    if row["q0"] is not None and row["q0"] > max_q0:
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
    # 读取过滤后的 per-sample CNV TSV，用于合并和统计。
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
    rows = []
    with open(args.input, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            row = parse_cnvnator_line(line, args.sample)
            if row and passes_filter(row, args.min_len, args.max_evalue, args.max_q0):
                rows.append(row)
    write_rows(args.output, rows)


def save_barplot(labels, values, path, title, ylabel):
    # 绘制柱状图；样本或染色体数量较多时自动调整画布宽度。
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(path), exist_ok=True)
    width = max(7, min(24, 0.35 * max(1, len(labels))))
    fig, ax = plt.subplots(figsize=(width, 4.8))
    ax.bar(labels, values, color="#2f7ecb", edgecolor="#1f4e79", linewidth=0.4)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=60, labelsize=8)
    ax.grid(axis="y", color="#dbe4ee", linewidth=0.6)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def save_length_hist(lengths, path):
    # 绘制 CNV 长度分布图，长度跨度过大时使用 log10 坐标。
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(path), exist_ok=True)
    values = [length for length in lengths if length > 0]
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    if values:
        log_values = [__import__("math").log10(value) for value in values]
        ax.hist(log_values, bins=40, color="#2f7ecb", edgecolor="#1f4e79", linewidth=0.4)
    ax.set_title("CNV length distribution")
    ax.set_xlabel("log10(CNV length, bp)")
    ax.set_ylabel("Count")
    ax.grid(axis="y", color="#dbe4ee", linewidth=0.6)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


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

    save_barplot(
        sorted(sample_counts),
        [sample_counts[sample] for sample in sorted(sample_counts)],
        args.count_plot,
        "CNV count by sample",
        "CNV count",
    )
    save_length_hist([row["length"] for row in all_rows], args.length_plot)
    save_barplot(
        sorted(chrom_counts),
        [chrom_counts[chrom] for chrom in sorted(chrom_counts)],
        args.chrom_plot,
        "CNV count by chromosome",
        "CNV count",
    )


def main():
    write_script_log()
    args = parse_args()
    if args.command == "filter":
        filter_one(args)
    elif args.command == "summary":
        summarize(args)


if __name__ == "__main__":
    main()
