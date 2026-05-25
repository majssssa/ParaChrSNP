#!/usr/bin/env Rscript

# 修改说明：
# 本脚本为 CNVnator 结果新增 Nature 风格绘图脚本。
# 原始绘图逻辑位于 scripts/cnvnator-summary_1.py 中，本脚本将绘图部分迁移到 R/ggplot2。

script_name <- "cnvnator-plot-nature_1.R"
script_function <- "Draw Nature-style CNVnator summary figures with R."
script_log <- "/home/majunpeng/script_log.txt"
script_path <- normalizePath(sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE)
try(
  suppressWarnings(cat(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), script_name, script_function, script_path, "\n",
    sep = "\t", file = script_log, append = TRUE
  )),
  silent = TRUE
)

usage <- paste(
  "Usage: Rscript scripts/cnvnator-plot-nature_1.R --combined <combined.cnv.tsv> --count-plot <pdf> --length-plot <pdf> --chrom-plot <pdf>",
  "",
  "Arguments:",
  "  --combined     Combined CNV table generated from CNVnator calls.",
  "  --count-plot   Output PDF path for CNV count by sample.",
  "  --length-plot  Output PDF path for CNV length distribution.",
  "  --chrom-plot   Output PDF path for CNV count by chromosome.",
  "",
  "The script also writes matching .svg and .tiff files for each PDF output.",
  sep = "\n"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1 && args[1] %in% c("-h", "--help")) {
  cat(usage, "\n")
  quit(status = 0)
}

get_arg <- function(flag, required = TRUE) {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit[1] == length(args)) {
    if (required) stop("Missing argument: ", flag)
    return(NA_character_)
  }
  args[hit[1] + 1]
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

combined_path <- get_arg("--combined")
count_plot <- get_arg("--count-plot")
length_plot <- get_arg("--length-plot")
chrom_plot <- get_arg("--chrom-plot")

cnv <- read.delim(combined_path, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("sample", "type", "chrom", "length")
missing_cols <- setdiff(required_cols, colnames(cnv))
if (length(missing_cols) > 0) {
  stop("Missing required columns in CNV table: ", paste(missing_cols, collapse = ", "))
}
cnv$length <- suppressWarnings(as.numeric(cnv$length))
cnv <- cnv[!is.na(cnv$length) & cnv$length > 0, , drop = FALSE]

palette_contract <- c(
  neutral_dark = "#272727",
  neutral_mid = "#767676",
  neutral_light = "#D8D8D8",
  deletion = "#4C78A8",
  duplication = "#E28E2C"
)

theme_nature <- function(base_size = 6.5, base_family = "Arial") {
  # Nature 风格基础主题：减少装饰，保留清晰的轴线、网格和可读标签。
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4, colour = "black"),
      plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.height = grid::unit(3.2, "mm"),
      legend.key.width = grid::unit(4.0, "mm"),
      legend.position = "top",
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "#E8E8E8"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
}

save_pub <- function(plot, pdf_path, width_mm = 89, height_mm = 66, dpi = 600) {
  # 输出 PDF/SVG/TIFF；PDF 路径由 Snakemake 作为主输出追踪。
  base <- sub("\\.pdf$", "", pdf_path, ignore.case = TRUE)
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  grDevices::svg(paste0(base, ".svg"), width = width_in, height = height_in, family = "Arial")
  print(plot)
  dev.off()

  grDevices::cairo_pdf(pdf_path, width = width_in, height = height_in, family = "Arial")
  print(plot)
  dev.off()

  grDevices::tiff(paste0(base, ".tiff"), width = width_in, height = height_in, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()
}

sample_counts <- as.data.frame(table(cnv$sample), stringsAsFactors = FALSE)
colnames(sample_counts) <- c("sample", "count")
sample_counts$count <- as.integer(sample_counts$count)
sample_counts$sample <- factor(sample_counts$sample, levels = sample_counts$sample[order(sample_counts$count)])

p_sample <- ggplot(sample_counts, aes(x = sample, y = count)) +
  geom_col(fill = "#A6CBE3", colour = palette_contract["neutral_dark"], linewidth = 0.22, width = 0.72) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(title = "CNV calls per sample", x = NULL, y = "CNV count") +
  theme_nature() +
  theme(legend.position = "none", panel.grid.major.x = element_line(linewidth = 0.2, colour = "#E8E8E8"))

p_length <- ggplot(cnv, aes(x = log10(length), fill = type)) +
  geom_histogram(bins = 40, colour = "white", linewidth = 0.12, alpha = 0.86, position = "identity") +
  scale_fill_manual(values = c(deletion = unname(palette_contract["deletion"]), duplication = unname(palette_contract["duplication"]))) +
  scale_x_continuous(labels = label_number(accuracy = 0.1), expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(title = "CNV length distribution", x = expression(log[10]("CNV length, bp")), y = "CNV count") +
  theme_nature()

chrom_counts <- as.data.frame(table(cnv$chrom, cnv$type), stringsAsFactors = FALSE)
colnames(chrom_counts) <- c("chrom", "type", "count")
chrom_totals <- aggregate(count ~ chrom, chrom_counts, sum)
chrom_counts$chrom <- factor(chrom_counts$chrom, levels = chrom_totals$chrom[order(chrom_totals$count)])

p_chrom <- ggplot(chrom_counts, aes(x = chrom, y = count, fill = type)) +
  geom_col(colour = "white", linewidth = 0.12, width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c(deletion = unname(palette_contract["deletion"]), duplication = unname(palette_contract["duplication"]))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Chromosomal distribution of CNVs", x = NULL, y = "CNV count") +
  theme_nature() +
  theme(panel.grid.major.x = element_line(linewidth = 0.2, colour = "#E8E8E8"))

save_pub(p_sample, count_plot, width_mm = 89, height_mm = max(62, 6 * length(unique(cnv$sample)) + 34))
save_pub(p_length, length_plot, width_mm = 89, height_mm = 66)
save_pub(p_chrom, chrom_plot, width_mm = 89, height_mm = max(70, 4 * length(unique(cnv$chrom)) + 36))
