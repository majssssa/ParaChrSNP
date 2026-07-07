#!/usr/bin/env Rscript

# Author: Junpeng Ma 1527552938@qq.com
# 功能：整理 ADMIXTURE 的 Q 矩阵和 CV error，并绘制群体结构图。
# 原始脚本：/home/majunpeng/Software/ParaChrSNP/scripts/admixture-plot-nature_1.R
# 修改要求：当提供 pop_info 文件时，ADMIXTURE 结构图中的样本顺序必须严格按照 pop_info 第一列出现的顺序绘制；未写入 pop_info 的 FAM 样本追加在最后，并保持 FAM 原始顺序。

script_name <- "admixture-plot-nature_2.R"
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  normalizePath("scripts/admixture-plot-nature_2.R", mustWork = FALSE)
}

write_script_log <- function() {
  line <- paste(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    script_name,
    "Draw ADMIXTURE population structure and CV error plots",
    script_path,
    sep = "\t"
  )
  suppressWarnings(try(write(line, file = "/home/majunpeng/script_log.txt", append = TRUE), silent = TRUE))
}

usage <- function() {
  cat("Usage: Rscript admixture-plot-nature_2.R --q-dir DIR --fam FILE --k-min INT --k-max INT --cv-out FILE --structure-out FILE --cv-plot-out FILE [--pop-info FILE] [--show-sample-names]\n\n")
  cat("Options:\n")
  cat("  --q-dir DIR           Directory containing admixture_pruned.K.Q and logs/K*.log files.\n")
  cat("  --fam FILE            PLINK FAM file used by ADMIXTURE.\n")
  cat("  --pop-info FILE       Optional two-column sample population file.\n")
  cat("  --k-min INT           Minimum K value.\n")
  cat("  --k-max INT           Maximum K value.\n")
  cat("  --cv-out FILE         Output TSV file for CV errors.\n")
  cat("  --structure-out FILE  Output PDF file for ADMIXTURE structure plot.\n")
  cat("  --cv-plot-out FILE    Output PDF file for CV error plot.\n")
  cat("  --show-sample-names   Show sample names on the x-axis of the structure plot.\n")
  cat("  -h, --help            Show this help message and exit.\n")
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 1 && args[1] %in% c("-h", "--help")) {
    usage()
    quit(save = "no", status = 0)
  }
  opts <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--show-sample-names") {
      opts[["show-sample-names"]] <- TRUE
      i <- i + 1
      next
    }
    if (!startsWith(key, "--") || i == length(args)) {
      stop(paste("Invalid argument:", key))
    }
    opts[[substring(key, 3)]] <- args[i + 1]
    i <- i + 2
  }
  required <- c("q-dir", "fam", "k-min", "k-max", "cv-out", "structure-out", "cv-plot-out")
  missing <- required[!required %in% names(opts)]
  if (length(missing) > 0) {
    usage()
    stop(paste("Missing required argument(s):", paste(missing, collapse = ", ")))
  }
  opts
}

read_fam_samples <- function(path) {
  fam <- read.table(path, stringsAsFactors = FALSE, fill = TRUE)
  if (ncol(fam) < 2) {
    stop("FAM file must contain at least two columns.")
  }
  fam[[2]]
}

read_pop_info <- function(path) {
  pop <- read.table(path, stringsAsFactors = FALSE, comment.char = "#")
  if (ncol(pop) < 2) {
    stop("pop_info must contain at least two columns: sample population")
  }
  pop <- pop[, 1:2]
  names(pop) <- c("sample", "population")
  pop
}

extract_cv <- function(log_path) {
  if (!file.exists(log_path)) {
    return(NA_real_)
  }
  lines <- readLines(log_path, warn = FALSE)
  hit <- grep("CV error", lines, value = TRUE)
  if (length(hit) == 0) {
    return(NA_real_)
  }
  nums <- regmatches(hit[length(hit)], gregexpr("[0-9]+\\.?[0-9]*(e[-+]?[0-9]+)?", hit[length(hit)], ignore.case = TRUE))[[1]]
  as.numeric(nums[length(nums)])
}

save_all <- function(pdf_path, plot_fun, width = 7.2, height = 4.8) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  stem <- sub("\\.pdf$", "", pdf_path, ignore.case = TRUE)
  pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  plot_fun()
  dev.off()
  png(paste0(stem, ".png"), width = width * 300, height = height * 300, res = 300, type = "cairo")
  plot_fun()
  dev.off()
  svg(paste0(stem, ".svg"), width = width, height = height)
  plot_fun()
  dev.off()
  tiff(paste0(stem, ".tiff"), width = width * 300, height = height * 300, res = 300, compression = "lzw", type = "cairo")
  plot_fun()
  dev.off()
}

draw_structure <- function(q_dir, samples, order_idx, k_values, show_sample_names = FALSE) {
  palette <- c("#3B6FB6", "#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3")
  bottom_margin <- if (show_sample_names) 4.6 else 0.4
  outer_bottom <- if (show_sample_names) 0.4 else 3.2
  par(mfrow = c(length(k_values), 1), mar = c(bottom_margin, 4.2, 0.25, 0.6), oma = c(outer_bottom, 0.2, 0.8, 0.2), xaxs = "i", yaxs = "i")
  for (k in k_values) {
    q_path <- file.path(q_dir, paste0("admixture_pruned.", k, ".Q"))
    q <- as.matrix(read.table(q_path, stringsAsFactors = FALSE))
    if (nrow(q) != length(samples)) {
      stop(paste("Sample count in", q_path, "does not match FAM file."))
    }
    q <- q[order_idx, , drop = FALSE]
    cols <- rep(palette, length.out = ncol(q))
    x_mid <- barplot(t(q), col = cols, border = NA, space = 0, axes = FALSE)
    axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), las = 1, cex.axis = 0.7, tick = FALSE)
    if (show_sample_names) {
      axis(1, at = x_mid, labels = samples[order_idx], las = 2, cex.axis = 0.65, tick = FALSE)
    }
    mtext(paste0("K=", k), side = 2, line = 2.8, las = 1, cex = 0.9)
    box(bty = "l", col = "#374151")
  }
  if (!show_sample_names) {
    mtext("Samples", side = 1, outer = TRUE, line = 1.4, cex = 0.9)
  }
}

draw_cv <- function(cv_df) {
  par(mar = c(4.4, 4.6, 1.2, 1), mgp = c(2.6, 0.8, 0), las = 1)
  plot(
    cv_df$K,
    cv_df$CV_error,
    type = "b",
    pch = 16,
    lwd = 1.8,
    col = "#3B6FB6",
    xlab = "K",
    ylab = "Cross-validation error",
    bty = "l",
    xaxt = "n"
  )
  axis(1, at = cv_df$K)
  grid(col = "#E5E7EB", lwd = 0.7)
}

write_script_log()
opts <- parse_args()
q_dir <- opts[["q-dir"]]
k_values <- seq(as.integer(opts[["k-min"]]), as.integer(opts[["k-max"]]))
samples <- read_fam_samples(opts$fam)
order_idx <- seq_along(samples)

if (!is.null(opts[["pop-info"]]) && nzchar(opts[["pop-info"]])) {
  pop <- read_pop_info(opts[["pop-info"]])
  ordered_samples <- unique(pop$sample)
  matched_idx <- match(ordered_samples, samples)
  matched_idx <- matched_idx[!is.na(matched_idx)]
  missing_idx <- setdiff(seq_along(samples), matched_idx)
  order_idx <- c(matched_idx, missing_idx)
}

cv_df <- data.frame(
  K = k_values,
  CV_error = vapply(k_values, function(k) extract_cv(file.path(q_dir, "logs", paste0("K", k, ".log"))), numeric(1))
)
dir.create(dirname(opts[["cv-out"]]), recursive = TRUE, showWarnings = FALSE)
write.table(cv_df, opts[["cv-out"]], sep = "\t", row.names = FALSE, quote = FALSE)

height <- max(3.6, 1.15 * length(k_values))
show_sample_names <- isTRUE(opts[["show-sample-names"]])
save_all(opts[["structure-out"]], function() draw_structure(q_dir, samples, order_idx, k_values, show_sample_names), width = 8, height = height)
save_all(opts[["cv-plot-out"]], function() draw_cv(cv_df), width = 5.4, height = 4.2)
