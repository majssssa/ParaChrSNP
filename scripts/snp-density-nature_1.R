#!/usr/bin/env Rscript

# 功能说明：
# 本脚本统计 VCF 文件中每条染色体在固定窗口内的 SNP 位点数量，
# 并绘制适用于论文和流程报告展示的染色体 SNP 密度热图。
# 输入 VCF 通常为 ParaChrSNP 过滤后的 SNP VCF。
#
# 作者：Junpeng Ma 1527552938@qq.com

script_name <- "snp-density-nature_1.R"
script_function <- "Calculate windowed SNP density and draw a chromosome density heatmap."
script_log <- "/home/majunpeng/script_log.txt"
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_file) > 0) {
  normalizePath(sub("^--file=", "", script_file[1]), mustWork = FALSE)
} else {
  script_name
}
try(
  suppressWarnings(cat(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), script_name, script_function,
    script_path, "Junpeng Ma 1527552938@qq.com\n",
    sep = "\t", file = script_log, append = TRUE
  )),
  silent = TRUE
)

usage <- paste(
  "Usage: Rscript scripts/snp-density-nature_1.R --vcf FILE --fai FILE --window-size INT --table FILE --out FILE [--chromosomes CHR1,CHR2]",
  "",
  "Draw a chromosome-level SNP density heatmap from a VCF file.",
  "",
  "Required arguments:",
  "  --vcf FILE          Input VCF or compressed VCF.GZ file containing SNP variants.",
  "  --fai FILE          FASTA index file used to obtain chromosome lengths.",
  "  --window-size INT   Non-overlapping window size in base pairs.",
  "  --table FILE        Output tab-delimited windowed SNP density table.",
  "  --out FILE          Output PDF figure; PNG, SVG and TIFF are written alongside it.",
  "",
  "Optional arguments:",
  "  --chromosomes LIST  Comma-separated chromosome order; default uses chromosomes observed in the VCF.",
  "  -h, --help          Show this English help document and exit.",
  "",
  "Author: Junpeng Ma 1527552938@qq.com",
  sep = "\n"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1 && args[1] %in% c("-h", "--help")) {
  cat(usage, "\n")
  quit(status = 0)
}

get_arg <- function(flag, default = NA_character_, required = TRUE) {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit[1] == length(args)) {
    if (required) stop("Missing argument: ", flag)
    return(default)
  }
  args[hit[1] + 1]
}

vcf_path <- get_arg("--vcf")
fai_path <- get_arg("--fai")
window_size <- suppressWarnings(as.numeric(get_arg("--window-size")))
table_path <- get_arg("--table")
out_pdf <- get_arg("--out")
chromosome_arg <- get_arg("--chromosomes", default = "", required = FALSE)

if (!file.exists(vcf_path)) stop("VCF file does not exist: ", vcf_path)
if (!file.exists(fai_path)) stop("FASTA index file does not exist: ", fai_path)
if (!is.finite(window_size) || window_size <= 0 || window_size %% 1 != 0) {
  stop("--window-size must be a positive integer.")
}
options(scipen = 999)

# 读取参考基因组索引，并保留绘图所需的染色体长度。
fai <- read.delim(fai_path, header = FALSE, stringsAsFactors = FALSE)
if (ncol(fai) < 2) stop("The FASTA index must contain chromosome names and lengths.")
colnames(fai)[1:2] <- c("chromosome", "length")
fai$length <- suppressWarnings(as.numeric(fai$length))
fai <- fai[is.finite(fai$length) & fai$length > 0, c("chromosome", "length"), drop = FALSE]

# 按窗口记录每个 SNP 位点。本模块输入为过滤后的 SNP VCF，因此每条非 header
# 记录均计为一个 SNP 位点，避免误删带星号等位基因的已提取 SNP 记录。
counts <- new.env(hash = TRUE, parent = emptyenv())
observed_chromosomes <- character(0)
input_handle <- if (grepl("\\.gz$", vcf_path, ignore.case = TRUE)) gzfile(vcf_path, "rt") else file(vcf_path, "rt")

repeat {
  lines <- readLines(input_handle, n = 100000, warn = FALSE)
  if (length(lines) == 0) break
  records <- lines[substr(lines, 1, 1) != "#"]
  if (length(records) == 0) next
  fields <- strsplit(records, "\t", fixed = TRUE)
  for (record in fields) {
    if (length(record) < 2) next
    chrom <- record[1]
    pos <- suppressWarnings(as.numeric(record[2]))
    if (!is.finite(pos)) next
    bin <- floor((pos - 1) / window_size) + 1
    key <- paste(chrom, bin, sep = "\t")
    counts[[key]] <- if (exists(key, counts, inherits = FALSE)) counts[[key]] + 1 else 1
    if (!(chrom %in% observed_chromosomes)) observed_chromosomes <- c(observed_chromosomes, chrom)
  }
}
close(input_handle)

if (nzchar(chromosome_arg)) {
  chrom_order <- strsplit(chromosome_arg, ",", fixed = TRUE)[[1]]
  chrom_order <- chrom_order[nzchar(chrom_order)]
} else {
  chrom_order <- fai$chromosome[fai$chromosome %in% observed_chromosomes]
}
missing_reference <- setdiff(chrom_order, fai$chromosome)
if (length(missing_reference) > 0) {
  stop("Chromosome(s) are not present in the FASTA index: ", paste(missing_reference, collapse = ", "))
}
if (length(chrom_order) == 0) stop("No chromosomes were selected for SNP density plotting.")

# 构造包含零 SNP 窗口在内的完整统计表，确保低密度区域可以显示为灰色。
window_tables <- lapply(chrom_order, function(chrom) {
  chrom_length <- fai$length[match(chrom, fai$chromosome)]
  n_bins <- ceiling(chrom_length / window_size)
  bin <- seq_len(n_bins)
  values <- vapply(bin, function(index) {
    key <- paste(chrom, index, sep = "\t")
    if (exists(key, counts, inherits = FALSE)) counts[[key]] else 0
  }, numeric(1))
  data.frame(
    chromosome = chrom,
    start = (bin - 1) * window_size + 1,
    end = pmin(bin * window_size, chrom_length),
    snp_count = values,
    stringsAsFactors = FALSE
  )
})
density <- do.call(rbind, window_tables)
dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
write.table(density, file = table_path, sep = "\t", quote = FALSE, row.names = FALSE)

# 生成离散且可解释的颜色等级：灰色表示无 SNP，绿色至红色表示密度逐步增加。
max_count <- max(density$snp_count)
positive_bins <- density$snp_count[density$snp_count > 0]
n_levels <- 9
if (length(positive_bins) == 0) {
  breaks <- c(0, 1)
  density_colors <- "#006837"
} else {
  breaks <- unique(round(seq(1, max_count, length.out = n_levels + 1)))
  if (length(breaks) < 2) breaks <- c(1, max_count + 1)
  density_colors <- colorRampPalette(c("#006837", "#66A61E", "#D9EF00", "#FFB21A", "#EF3B2C"))(length(breaks) - 1)
}
get_color <- function(values) {
  result <- rep("#D9D9D9", length(values))
  positive <- values > 0
  if (any(positive)) {
    result[positive] <- density_colors[findInterval(values[positive], breaks, all.inside = TRUE)]
  }
  result
}

base <- sub("\\.pdf$", "", out_pdf, ignore.case = TRUE)
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
max_mb <- max(density$end) / 1e6
height_in <- max(2.6, 0.31 * length(chrom_order) + 1.35)
width_in <- 7.2

format_density_value <- function(value) {
  if (value >= 1000) {
    label <- format(round(value / 1000, 1), trim = TRUE, scientific = FALSE)
    return(paste0(label, "k"))
  }
  format(value, trim = TRUE, scientific = FALSE)
}

draw_plot <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  layout(matrix(c(1, 2), nrow = 1), widths = c(8.6, 1.45))

  par(mar = c(1.0, 6.1, 4.1, 0.35), family = "Arial", xaxs = "i", yaxs = "i")
  n_chrom <- length(chrom_order)
  plot(
    NA, xlim = c(0, max_mb), ylim = c(0.5, n_chrom + 0.5),
    axes = FALSE, ann = FALSE, bty = "n"
  )
  for (i in seq_along(chrom_order)) {
    chrom <- chrom_order[i]
    rows <- density[density$chromosome == chrom, , drop = FALSE]
    y <- n_chrom - i + 1
    rect(
      xleft = (rows$start - 1) / 1e6, ybottom = y - 0.31,
      xright = rows$end / 1e6, ytop = y + 0.31,
      col = get_color(rows$snp_count), border = NA
    )
  }
  axis(2, at = rev(seq_along(chrom_order)), labels = chrom_order, las = 1, tick = FALSE, line = -0.4, cex.axis = 0.61)
  ticks <- pretty(c(0, max_mb), n = 6)
  ticks <- ticks[ticks >= 0 & ticks <= max_mb]
  axis(3, at = ticks, labels = paste0(format(ticks, trim = TRUE), " Mb"), lwd = 0.8, tck = -0.025, cex.axis = 0.64, line = 0.05)
  title(main = paste0("SNP density (", format(window_size / 1e6, trim = TRUE), " Mb windows)"), line = 2.75, cex.main = 0.78, font.main = 2)

  par(mar = c(1.0, 0.35, 4.1, 0.15), family = "Arial")
  legend_labels <- if (length(density_colors) == 1) {
    if (max_count == 0) "0" else paste0("1-", max_count)
  } else {
    paste0(vapply(head(breaks, -1), format_density_value, character(1)), "-", vapply(tail(breaks, -1), format_density_value, character(1)))
  }
  legend_colors <- c("#D9D9D9", density_colors)
  legend_labels <- c("0", legend_labels)
  plot(NA, xlim = c(0, 1.7), ylim = c(0, length(legend_colors) + 2), axes = FALSE, ann = FALSE, bty = "n")
  for (j in seq_along(legend_colors)) {
    y <- length(legend_colors) - j + 1
    rect(0.04, y, 0.38, y + 0.9, col = legend_colors[j], border = NA)
    text(0.46, y + 0.45, labels = legend_labels[j], adj = 0, cex = 0.61)
  }
  text(0.04, length(legend_colors) + 1.25, labels = "SNPs", adj = 0, font = 2, cex = 0.68)
}

svg(paste0(base, ".svg"), width = width_in, height = height_in, family = "Arial")
draw_plot()
dev.off()

cairo_pdf(out_pdf, width = width_in, height = height_in, family = "Arial")
draw_plot()
dev.off()

png(paste0(base, ".png"), width = width_in, height = height_in, units = "in", res = 300, type = "cairo")
draw_plot()
dev.off()

tiff(paste0(base, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw")
draw_plot()
dev.off()
