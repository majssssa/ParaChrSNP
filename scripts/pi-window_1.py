#!/usr/bin/env python3
"""
Calculate windowed nucleotide diversity (Pi) from a VCF file.

Author: Junpeng Ma 1527552938@qq.com
"""

import argparse
import gzip
import os
from collections import defaultdict
from datetime import datetime


AUTHOR = "Junpeng Ma 1527552938@qq.com"


def write_script_log():
    """记录脚本信息；容器只读时自动跳过，不影响主流程。"""
    try:
        with open("/home/majunpeng/script_log.txt", "a", encoding="utf-8") as handle:
            handle.write(
                "\t".join([
                    datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "pi-window_1.py",
                    "Calculate windowed nucleotide diversity from VCF for ParaChrSNP.",
                    os.path.abspath(__file__),
                    AUTHOR,
                ]) + "\n"
            )
    except OSError:
        pass


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate windowed nucleotide diversity (Pi) from a biallelic SNP VCF. "
            "The script accepts an optional two-column population file: sample population."
        ),
        epilog=f"Author: {AUTHOR}",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--vcf", required=True, help="Input VCF/VCF.GZ file.")
    parser.add_argument("--pop-info", default="", help="Optional sample population file. Column 1 is sample; column 2 is population.")
    parser.add_argument("--window-size", type=int, default=100000, help="Window size in bp.")
    parser.add_argument("--window-step", type=int, default=10000, help="Window step in bp.")
    parser.add_argument("--output-dir", default="pi", help="Output directory.")
    parser.add_argument("--combined", required=True, help="Combined windowed Pi output TSV.")
    parser.add_argument("--summary", required=True, help="Pi summary output TSV.")
    return parser.parse_args()


def smart_open(path):
    """根据后缀自动打开普通文本或 gzip 压缩文件。"""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_pop_info(path, vcf_samples):
    """读取群体文件；未提供时所有样本归为 All。"""
    if not path:
        return {"All": list(range(len(vcf_samples)))}

    sample_to_index = {sample: idx for idx, sample in enumerate(vcf_samples)}
    groups = defaultdict(list)
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            sample = fields[0]
            population = fields[1] if len(fields) > 1 else "All"
            if sample in sample_to_index:
                groups[population].append(sample_to_index[sample])

    if not groups:
        raise ValueError(f"No samples from pop-info were found in VCF header: {path}")
    return dict(sorted(groups.items()))


def is_biallelic_snp(ref, alt):
    """只保留双等位 SNP，跳过 INDEL、多等位和结构变异。"""
    return len(ref) == 1 and "," not in alt and len(alt) == 1 and alt not in {".", "*"}


def genotype_alleles(sample_field):
    """从样本字段中提取 GT，并返回 0/1 等位基因；缺失或复杂基因型返回空列表。"""
    gt = sample_field.split(":", 1)[0]
    if gt in {".", "./.", ".|."}:
        return []
    alleles = gt.replace("|", "/").split("/")
    parsed = []
    for allele in alleles:
        if allele == ".":
            return []
        if allele not in {"0", "1"}:
            return []
        parsed.append(int(allele))
    return parsed


def window_starts_for_site(pos, window_size, window_step):
    """返回覆盖当前位点的滑动窗口起点，窗口为 1-based 闭区间。"""
    first = ((max(1, pos - window_size + 1) - 1) // window_step) * window_step + 1
    starts = []
    start = first
    while start <= pos:
        if start <= pos <= start + window_size - 1:
            starts.append(start)
        start += window_step
    return starts


def calculate_pi(args):
    """主计算函数：逐行读取 VCF，按群体和窗口累加 Pi。"""
    write_script_log()
    os.makedirs(args.output_dir, exist_ok=True)

    samples = []
    groups = None
    windows = defaultdict(lambda: {"sum_pi": 0.0, "n_variants": 0})
    callable_sites = defaultdict(int)
    total_biallelic_snps = 0

    with smart_open(args.vcf) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                groups = parse_pop_info(args.pop_info, samples)
                continue
            if not line.strip():
                continue
            if groups is None:
                raise ValueError("VCF header line #CHROM was not found before records.")

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos_text, _vid, ref, alt = fields[:5]
            if not is_biallelic_snp(ref, alt):
                continue
            try:
                pos = int(pos_text)
            except ValueError:
                continue

            total_biallelic_snps += 1
            sample_fields = fields[9:]
            starts = window_starts_for_site(pos, args.window_size, args.window_step)
            for population, indexes in groups.items():
                ref_count = 0
                alt_count = 0
                for idx in indexes:
                    if idx >= len(sample_fields):
                        continue
                    for allele in genotype_alleles(sample_fields[idx]):
                        if allele == 0:
                            ref_count += 1
                        elif allele == 1:
                            alt_count += 1
                allele_count = ref_count + alt_count
                if allele_count < 2:
                    continue
                # 每个位点的平均成对差异数：2ab / n(n-1)。
                site_pi = (2.0 * ref_count * alt_count) / (allele_count * (allele_count - 1))
                for start in starts:
                    key = (population, chrom, start, start + args.window_size - 1)
                    windows[key]["sum_pi"] += site_pi
                    windows[key]["n_variants"] += 1
                    callable_sites[(population, chrom)] += 1

    if groups is None:
        raise ValueError("No VCF sample header was found.")

    with open(args.combined, "w", encoding="utf-8") as out:
        out.write("population\tCHROM\tBIN_START\tBIN_END\tN_VARIANTS\tSUM_PI\tPI\tMEAN_SITE_PI\n")
        for key in sorted(windows, key=lambda x: (x[0], x[1], x[2], x[3])):
            population, chrom, start, end = key
            record = windows[key]
            n_variants = record["n_variants"]
            sum_pi = record["sum_pi"]
            pi_per_bp = sum_pi / args.window_size
            mean_site_pi = sum_pi / n_variants if n_variants else 0.0
            out.write(
                f"{population}\t{chrom}\t{start}\t{end}\t{n_variants}\t"
                f"{sum_pi:.10g}\t{pi_per_bp:.10g}\t{mean_site_pi:.10g}\n"
            )

    per_population = defaultdict(lambda: {"windows": 0, "sum_pi": 0.0, "max_pi": 0.0})
    for key, record in windows.items():
        population = key[0]
        pi_per_bp = record["sum_pi"] / args.window_size
        per_population[population]["windows"] += 1
        per_population[population]["sum_pi"] += pi_per_bp
        per_population[population]["max_pi"] = max(per_population[population]["max_pi"], pi_per_bp)

    with open(args.summary, "w", encoding="utf-8") as out:
        out.write("population\tsample_count\tbiallelic_snps\twindow_count\tmean_window_pi\tmax_window_pi\n")
        for population, indexes in sorted(groups.items()):
            stats = per_population[population]
            mean_pi = stats["sum_pi"] / stats["windows"] if stats["windows"] else 0.0
            out.write(
                f"{population}\t{len(indexes)}\t{total_biallelic_snps}\t{stats['windows']}\t"
                f"{mean_pi:.10g}\t{stats['max_pi']:.10g}\n"
            )


def main():
    args = parse_args()
    if args.window_size <= 0 or args.window_step <= 0:
        raise ValueError("--window-size and --window-step must be positive integers.")
    calculate_pi(args)


if __name__ == "__main__":
    main()
