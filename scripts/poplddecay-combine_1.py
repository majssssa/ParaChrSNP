#!/usr/bin/env python3
"""
Combine PopLDdecay stat.gz files into a tidy table.

Author: Junpeng Ma 1527552938@qq.com
"""

import argparse
import gzip
import os
from datetime import datetime


AUTHOR = "Junpeng Ma 1527552938@qq.com"
SCRIPT_NAME = "poplddecay-combine_1.py"
SCRIPT_PATH = os.path.abspath(__file__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine PopLDdecay .stat.gz outputs into one tidy TSV table."
    )
    parser.add_argument("--stats-dir", required=True, help="Directory containing PopLDdecay *.stat.gz files.")
    parser.add_argument("--combined", required=True, help="Output combined TSV file.")
    return parser.parse_args()


def write_script_log():
    # 按用户要求记录脚本名称、功能、路径和创建时间；失败时不影响主流程。
    try:
        with open("/home/majunpeng/script_log.txt", "a", encoding="utf-8") as handle:
            handle.write(
                f"{datetime.now().isoformat()}\t{SCRIPT_NAME}\t"
                f"Combine PopLDdecay stat files for LD decay plotting\t{SCRIPT_PATH}\n"
            )
    except OSError:
        pass


def parse_numeric_fields(parts):
    # PopLDdecay 的 stat 文件版本可能略有不同；这里取每行前两个可解析的数值列。
    values = []
    for item in parts:
        try:
            values.append(float(item))
        except ValueError:
            continue
        if len(values) == 2:
            return values
    return None


def iter_stat_rows(path):
    # 读取压缩 stat 文件，跳过注释行和表头行。
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            values = parse_numeric_fields(line.split())
            if values is None:
                continue
            yield values[0], values[1]


def population_name(path):
    name = os.path.basename(path)
    if name.endswith(".stat.gz"):
        return name[:-8]
    return os.path.splitext(name)[0]


def main():
    write_script_log()
    args = parse_args()
    stat_files = sorted(
        os.path.join(args.stats_dir, name)
        for name in os.listdir(args.stats_dir)
        if name.endswith(".stat.gz")
    )
    if not stat_files:
        raise SystemExit(f"No PopLDdecay .stat.gz files found in {args.stats_dir}")

    os.makedirs(os.path.dirname(args.combined) or ".", exist_ok=True)
    row_count = 0
    with open(args.combined, "w", encoding="utf-8") as out:
        out.write("population\tdistance_bp\tr2\n")
        for path in stat_files:
            pop = population_name(path)
            for distance, r2 in iter_stat_rows(path):
                out.write(f"{pop}\t{distance:g}\t{r2:g}\n")
                row_count += 1

    if row_count == 0:
        raise SystemExit("PopLDdecay stat files were found, but no numeric LD rows could be parsed.")


if __name__ == "__main__":
    main()
