#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Junpeng Ma 1527552938@qq.com
Created: 2026-04-25 10:08:46 CST

修改要求：
将原始 R 脚本重写为 Python 脚本，保持相同的输入参数、输出文件名和绘图效果。

原始脚本位置：
/home/majunpeng/Software/ParaChrSNP/scripts/plink-missing.R
"""

import argparse
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    """解析命令行参数。"""
    parser = argparse.ArgumentParser(
        description=(
            "Plot PLINK variant-level and sample-level missing rate histograms. "
            "This script is a Python rewrite of scripts/plink-missing.R."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "lmiss",
        help=(
            "Input PLINK .lmiss file. This file is produced by PLINK --missing "
            "and contains variant-level missing rate information. The table must "
            "include a column named F_MISS."
        ),
    )
    parser.add_argument(
        "imiss",
        help=(
            "Input PLINK .imiss file. This file is produced by PLINK --missing "
            "and contains sample-level missing rate information. The table must "
            "include a column named F_MISS."
        ),
    )
    parser.add_argument(
        "out_prefix",
        help=(
            "Output prefix. The script will generate two PDF files: "
            "<out_prefix>.sample_missing.pdf and <out_prefix>.Variant_missing.pdf."
        ),
    )
    return parser.parse_args()


def read_plink_missing_table(input_file):
    """读取 PLINK missing 结果表，并检查 F_MISS 列是否存在。"""
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")

    # PLINK 的 .lmiss/.imiss 文件通常是空白字符分隔，因此使用正则空白分隔符读取。
    data = pd.read_csv(input_path, sep=r"\s+", header=0)
    if "F_MISS" not in data.columns:
        raise ValueError(f"Column F_MISS was not found in {input_file}")

    # 将 F_MISS 转换成数值，避免文本或异常字符影响统计和绘图。
    data["F_MISS"] = pd.to_numeric(data["F_MISS"], errors="coerce")
    data = data.dropna(subset=["F_MISS"])
    if data.empty:
        raise ValueError(f"No valid numeric F_MISS values were found in {input_file}")

    return data


def build_ggplot_like_bins(values, bin_count=30):
    """构建接近 ggplot2 geom_histogram(boundary = 0) 默认效果的分箱。"""
    min_value = float(np.nanmin(values))
    max_value = float(np.nanmax(values))

    if math.isclose(min_value, max_value):
        # 如果所有缺失率完全相同，给一个很窄但稳定的绘图区间，避免直方图无法分箱。
        width = 0.01 if math.isclose(min_value, 0.0) else abs(min_value) * 0.01
        return np.array([min_value - width, min_value + width])

    bin_width = (max_value - min_value) / bin_count
    start = math.floor(min_value / bin_width) * bin_width
    end = math.ceil(max_value / bin_width) * bin_width
    bins = np.arange(start, end + bin_width, bin_width)

    # 确保最大值落入最后一个 bin，避免浮点误差导致边界缺失。
    if bins[-1] < max_value:
        bins = np.append(bins, bins[-1] + bin_width)
    return bins


def plot_missing_histogram(data, y_label, output_pdf):
    """绘制缺失率直方图，并输出为 PDF。"""
    values = data["F_MISS"].to_numpy()
    total = len(values)
    mean_value = np.mean(values)

    # R 的 sd() 默认使用样本标准差，对应 numpy/pandas 中 ddof=1。
    sd_value = np.std(values, ddof=1) if total > 1 else 0.0
    label = f"Total:{total}  Mean:{mean_value:.3f}  SD:{sd_value:.3f}"

    bins = build_ggplot_like_bins(values)

    # R 脚本中 ggsave 的尺寸是 width=10, height=7，这里使用相同英寸尺寸。
    fig, ax = plt.subplots(figsize=(10, 7))

    # 复刻 ggplot2: geom_histogram(color = "black", fill = "orange", boundary = 0)。
    ax.hist(values, bins=bins, color="orange", edgecolor="black", linewidth=1.0)

    ax.set_xlabel("Missing rate")
    ax.set_ylabel(y_label)

    # 复刻 theme_bw 的主要视觉效果：白色背景、黑色边框、浅灰网格。
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")
    ax.grid(True, color="#EBEBEB", linewidth=1.0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)

    # 对应 ggplot2 annotate("text", x=Inf, y=Inf, vjust=1.5, hjust=1.1) 的右上角标注。
    ax.text(
        0.99,
        0.96,
        label,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=11,
    )

    fig.tight_layout()
    fig.savefig(output_pdf)
    plt.close(fig)


def main():
    """主函数：读取输入文件并生成两个缺失率 PDF 图。"""
    args = parse_args()

    lmiss = read_plink_missing_table(args.lmiss)
    imiss = read_plink_missing_table(args.imiss)

    sample_output = f"{args.out_prefix}.sample_missing.pdf"
    variant_output = f"{args.out_prefix}.Variant_missing.pdf"

    plot_missing_histogram(imiss, "Sample number", sample_output)
    plot_missing_histogram(lmiss, " Variant number", variant_output)

    print(f"Sample missing plot: {sample_output}")
    print(f"Variant missing plot: {variant_output}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
