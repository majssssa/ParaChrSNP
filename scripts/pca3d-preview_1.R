#!/usr/bin/env Rscript

# 修改说明：
# 这个脚本用于从 VCF2PCACluster 生成的 eigenvec 文件中绘制 3D PCA 报告预览图。
# 原始第三方软件 Plot3Deig 生成的 PDF 不做修改；本脚本只额外生成 PNG/SVG 预览，
# 解决 HTML 报告在部分浏览器或本地打开时无法直接显示 PDF 的问题。
# 原始结果来源：rules/vcf2pca.rules 中 Plot3Deig 输出的 3D PCA PDF。

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || is.na(x) || x == "") y else x
}

AUTHOR <- "Junpeng Ma 1527552938@qq.com"
SCRIPT_NAME <- "pca3d-preview_1.R"
SCRIPT_FUNCTION <- "Generate lightweight 3D PCA preview images from VCF2PCACluster eigenvec output for ParaChrSNP HTML reports."
SCRIPT_PATH <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% SCRIPT_NAME), mustWork = FALSE)
CREATED_TIME <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

write_script_log <- function() {
    log_file <- "/home/majunpeng/script_log.txt"
    line <- paste(
        CREATED_TIME,
        SCRIPT_NAME,
        SCRIPT_FUNCTION,
        SCRIPT_PATH,
        AUTHOR,
        sep = "\t"
    )
    try(write(line, file = log_file, append = TRUE), silent = TRUE)
}

usage <- function() {
    cat("
pca3d-preview_1.R

Description:
  Generate PNG and SVG preview images for 3D PCA results from a VCF2PCACluster eigenvec file.
  The original Plot3Deig PDF output is not modified.

Usage:
  Rscript scripts/pca3d-preview_1.R <eigenvec> <output_prefix>

Arguments:
  eigenvec:
    VCF2PCACluster eigenvec file. It must contain SampleName, Group, Cluster, PC1, PC2 and PC3 columns.

  output_prefix:
    Output prefix for preview figures. The script writes:
      <output_prefix>.C.3DPC1PC2PC3.png
      <output_prefix>.C.3DPC1PC2PC3.svg
      <output_prefix>.N.3DPC1PC2PC3.png
      <output_prefix>.N.3DPC1PC2PC3.svg

Author:
  Junpeng Ma 1527552938@qq.com
\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1 && args[1] == "-h") {
    usage()
    quit(status = 0)
}
if (length(args) != 2) {
    usage()
    quit(status = 1)
}

write_script_log()

eigenvec_file <- args[1]
output_prefix <- args[2]

if (!file.exists(eigenvec_file)) {
    stop("Input eigenvec file does not exist: ", eigenvec_file)
}

dat <- read.table(eigenvec_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("SampleName", "PC1", "PC2", "PC3")
missing_columns <- setdiff(required_columns, colnames(dat))
if (length(missing_columns) > 0) {
    stop("Missing required columns in eigenvec file: ", paste(missing_columns, collapse = ", "))
}

dat$PC1 <- as.numeric(dat$PC1)
dat$PC2 <- as.numeric(dat$PC2)
dat$PC3 <- as.numeric(dat$PC3)
dat <- dat[complete.cases(dat[, c("PC1", "PC2", "PC3")]), ]
if (nrow(dat) == 0) {
    stop("No valid PCA coordinates were found in eigenvec file.")
}

if (!"Group" %in% colnames(dat)) {
    dat$Group <- "UnGroup"
}

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

nature_palette <- c(
    "#3B6FB6", "#D95F02", "#1B9E77", "#7570B3",
    "#E7298A", "#66A61E", "#A6761D", "#4D4D4D"
)

draw_one <- function(file, colored = TRUE, device = c("png", "svg")) {
    device <- match.arg(device)
    if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
        stop("The R package 'scatterplot3d' is required for PCA 3D preview plotting.")
    }
    if (device == "png") {
        png(file, width = 1800, height = 1500, res = 300, type = "cairo")
    } else {
        svg(file, width = 6.0, height = 5.0, onefile = FALSE)
    }
    on.exit(dev.off(), add = TRUE)

    groups <- if (colored) dat$Group else rep("All samples", nrow(dat))
    group_levels <- unique(groups)
    color_map <- setNames(rep(nature_palette, length.out = length(group_levels)), group_levels)
    point_colors <- unname(color_map[groups])

    par(
        mar = c(3.6, 3.6, 1.0, 1.2),
        mgp = c(2.4, 0.65, 0),
        tcl = -0.25,
        family = "sans",
        cex.axis = 0.72,
        cex.lab = 0.82
    )

    plot3d <- scatterplot3d::scatterplot3d(
        x = dat$PC1,
        y = dat$PC2,
        z = dat$PC3,
        color = point_colors,
        pch = 19,
        cex.symbols = 1.25,
        angle = 45,
        scale.y = 0.75,
        grid = TRUE,
        box = TRUE,
        xlab = "PC1",
        ylab = "PC2",
        zlab = "PC3",
        main = "",
        axis = TRUE,
        tick.marks = TRUE,
        label.tick.marks = TRUE
    )

    if (nrow(dat) <= 30) {
        label_pos <- plot3d$xyz.convert(dat$PC1, dat$PC2, dat$PC3)
        text(label_pos$x, label_pos$y, labels = dat$SampleName, pos = 3, cex = 0.5, col = "#333333")
    }

    if (colored && length(group_levels) > 1) {
        legend(
            "topright",
            legend = group_levels,
            pt.bg = unname(color_map[group_levels]),
            pch = 21,
            col = "white",
            bty = "n",
            cex = 0.68,
            title = "Group"
        )
    }
}

for (mode in c("C", "N")) {
    colored <- mode == "C"
    png_file <- paste0(output_prefix, ".", mode, ".3DPC1PC2PC3.png")
    svg_file <- paste0(output_prefix, ".", mode, ".3DPC1PC2PC3.svg")
    draw_one(png_file, colored = colored, device = "png")
    draw_one(svg_file, colored = colored, device = "svg")
}
