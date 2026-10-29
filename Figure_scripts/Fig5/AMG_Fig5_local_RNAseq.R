#!/usr/bin/env Rscript
# Figure 5: local RNA-seq effects and strand-specific RPKM coverage.
# Run from AMG_cleaned_Git.Rproj. No other R script or prepared coverage is needed.
# Inputs: RNAseq_data/DE/PAO1_*_All.xlsx, annotation_input/PAO1_curated.tsv and
# RNAseq_data/coverage/<sample>_fw.bigwig / <sample>_re.bigwig.
# Upstream: use the supplied complete DE workbooks and RPKM tracks directly.
# To rebuild from reads, follow RNAseq_analysis/preprocessing/README.md:
# 01_download_SRA_reads.sh -> 02_rename_SRA_reads.sh -> 03_preprocess_RNAseq_reads.sh.
# RNAseq_analysis/AMG_RNAseq_featureCounts_edgeR.R refits DE from counts in a
# separate out/ directory. Review and copy new DE workbooks and coverage tracks
# to RNAseq_data/DE/ and RNAseq_data/coverage/ before adopting them here.
# First-sheet DE columns must include PAO1_ID, logFC and FDR; annotation supplies
# genomic start/end/strand.
# Each window is divided into 50-bp bins; uncovered bases contribute zero.
# Coverage is averaged within bins, then mean +/- SD is calculated across the
# three biological replicates. Reverse-strand values are negative for display.
# All coordinates are 1-based inclusive. BigWigs must contain nonnegative RPKM
# magnitudes on each strand; read alignment and RPKM normalization are upstream.
# Outputs: combined/individual figures, regional gene effects, coverage mean/SD,
# sample-level coverage, input-file assignments and settings as TSV/TSV.gz.
#
# Run: source("Figure_scripts/Fig5/AMG_Fig5_local_RNAseq.R")
# Output folder: Figure_scripts/Fig5/out_AMG_Fig5_local_RNAseq/.
#
# R packages used/checked: GenomeInfoDb, GenomicRanges, IRanges, S4Vectors, dplyr,
# gggenes, ggplot2, ggrepel, ggtext, openxlsx, patchwork, purrr, readr, rtracklayer,
# scales, tibble, tidyr.

# 1. Input files and output directory
DE_DIR <- "RNAseq_data/DE"

ANNOTATION_FILE <- "annotation_input/PAO1_curated.tsv"

# Files are named <sample>_fw.bigwig and <sample>_re.bigwig.
BIGWIG_DIR <- "RNAseq_data/coverage"

OUTPUT_DIR <- "Figure_scripts/Fig5/out_AMG_Fig5_local_RNAseq"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# GENE-SPECIFIC SETTINGS
# ============================================================
#
# This is the main section to edit.
#
# region_start / region_end are explicit plotting ranges.
# guide_sequence is the PARENT full-match protospacer sequence
# identifying the parent sequence; coordinates are supplied in GENE_CONFIG.
#
# The guide coordinates below refer to PAO1 NC_002516.2 and are used directly;
# this script does not realign parent sequences to recalculate their positions.
#
# Gene order in the final figure:
#   atpI -> accC -> secA -> ftsZ
GENE_CONFIG <- list(atpI = list(target_gene = "atpI", target_id = "PA5561", guide_display = "g_atpI1-m2", 
  guide_name = "atpIg1m2", guide_sequence = "TCCCAGCAGTCCGGAATATC", guide_start = 6254689L, guide_end = 6254708L, 
  guide_strand = "+", region_start = 6239500L, region_end = 6264300L, silenced_condition = "atpIg1m2", 
  silenced_samples = c("16_atpIg1m2", "17_atpIg1m2", "18_atpIg1m2"), de_stem = "PAO1_atpIg1m2_vs_control1_All"), 
  accC = list(target_gene = "accC", target_id = "PA4848", guide_display = "g_accC1-m3", guide_name = "accCg1m3", 
    guide_sequence = "CGGTGACCATCTCAGATACC", guide_start = 5443685L, guide_end = 5443704L, 
    guide_strand = "-", region_start = 5438500L, region_end = 5449000L, silenced_condition = "accCg1m3", 
    silenced_samples = c("19_accCg1m3", "20_accCg1m3", "21_accCg1m3"), de_stem = "PAO1_accCg1m3_vs_control1_All"), 
  secA = list(target_gene = "secA", target_id = "PA4403", guide_display = "g_secA1-m3", guide_name = "secAg1m3", 
    guide_sequence = "CAACAGCTTGGAAACGTACT", guide_start = 4935210L, guide_end = 4935229L, 
    guide_strand = "+", region_start = 4930000L, region_end = 4940500L, silenced_condition = "secAg1m3", 
    silenced_samples = c("22_secAg1m3", "23_secAg1m3", "24_secAg1m3"), de_stem = "PAO1_secAg1m3_vs_control1_All"), 
  ftsZ = list(target_gene = "ftsZ", target_id = "PA4407", guide_display = "g_ftsZ1-m4", guide_name = "ftsZg1m4", 
    guide_sequence = "ACCGATCACTTTTATGACGG", guide_start = 4940431L, guide_end = 4940450L, 
    guide_strand = "+", region_start = 4935200L, region_end = 4945600L, silenced_condition = "ftsZg1m4", 
    silenced_samples = c("13_ftsZg1m4", "14_ftsZg1m4", "15_ftsZg1m4"), de_stem = "PAO1_ftsZg1m4_vs_control1_All"))

GENE_ORDER <- c("atpI", "accC", "secA", "ftsZ")

# ============================================================
# CONTROL SAMPLES
# ============================================================
CONTROL_CONDITION <- "control1"

CONTROL_SAMPLES <- c("10_control1", "11_control1", "12_control1")

# ============================================================
# COVERAGE / DE SETTINGS
# ============================================================
COVERAGE_BIN_SIZE_BP <- 50L

# Set only if automatic longest-chromosome detection is wrong.
BIGWIG_SEQNAME_OVERRIDE <- NA_character_

DE_SIGNIFICANCE_THRESHOLD <- 0.05

# RPKM display.
SHOW_COVERAGE_SD_RIBBON <- TRUE

# Signed RPKM is overlaid directly on the log2FC panel.
# This fraction controls how much of the available log2FC vertical
# range is occupied by the largest absolute RPKM value.
COVERAGE_OVERLAY_FRACTION <- 0.9

COVERAGE_LINE_ALPHA <- 0.78

COVERAGE_RIBBON_ALPHA <- 0.18

RPKM_DISPLAY_DIVISOR <- 1000

# Gene-order strip.
GENE_STRIP_EXTRA_SPACE <- 1.1

# Labels.
SHOW_LOCAL_GENE_LABELS <- TRUE

SHOW_DE_SIGNIFICANCE_STARS <- TRUE

# ============================================================
# VISUAL SETTINGS
# ============================================================
FONT_FAMILY <- "Arial"

BASE_SIZE <- 8.5

CONTROL_COLOR <- "#4D4D4D"

KNOCKDOWN_COLOR <- "#D95F02"

TARGET_GENE_COLOR <- "#8A0000"

TARGET_OPERON_COLOR <- "#E69F00"

NEIGHBOR_GENE_COLOR <- "#C8C8C8"

DE_UP_COLOR <- "#B2182B"

DE_DOWN_COLOR <- "#2166AC"

DE_NS_COLOR <- "#8C8C8C"

GUIDE_SITE_COLOR <- "#7B3294"

# Final 2 x 2 manuscript figure.
FINAL_WIDTH_MM <- 176

FINAL_HEIGHT_MM <- 128

SAVE_PNG <- TRUE

SAVE_PDF <- TRUE

PNG_DPI <- 300

PRINT_IN_RSTUDIO <- TRUE

# ============================================================
# HELPERS
# ============================================================
required_packages <- c("dplyr", "tidyr", "purrr", "readr", "tibble", "ggplot2", "ggrepel", "gggenes", 
  "ggtext", "patchwork", "scales", "openxlsx", "rtracklayer", "GenomicRanges",
  "IRanges", "GenomeInfoDb", "S4Vectors")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_packages) > 0L) {
  stop("Install the following packages before running:\n", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(gggenes)
  library(ggtext)
  library(patchwork)
  library(scales)
  library(openxlsx)
})

# ============================================================
# INPUT READERS
# ============================================================
find_column <- function(data, candidates, required = TRUE) {
  matched <- match(tolower(candidates), tolower(names(data)))
  matched <- matched[!is.na(matched)]
  if (length(matched) == 0L) {
    if (isTRUE(required)) {
      stop("None of the expected columns were found: ", paste(candidates, collapse = ", "))
    }
    return(NA_character_)
  }
  names(data)[matched[[1L]]]
}

read_de_table <- function(path) {
  extension <- tolower(tools::file_ext(path))
  if (extension %in% c("xlsx", "xls")) {
    return(openxlsx::read.xlsx(path, sheet = 1) %>% tibble::as_tibble())
  }
  if (extension == "csv") {
    return(readr::read_csv(path, show_col_types = FALSE))
  }
  if (extension %in% c("tsv", "txt")) {
    return(readr::read_tsv(path, show_col_types = FALSE))
  }
  stop("Unsupported DE table extension: ", extension)
}

find_de_file <- function(de_stem) {
  path <- file.path(DE_DIR, paste0(de_stem, ".xlsx"))
  if (!file.exists(path)) 
    stop("Missing DE contrast: ", path)
  path
}

significance_code <- function(adjusted_p) {
  dplyr::case_when(is.na(adjusted_p) ~ "", adjusted_p < 0.001 ~ "***", adjusted_p < 0.01 ~ "**", 
    adjusted_p < 0.05 ~ "*", TRUE ~ "")
}

format_p <- function(x) {
  if (!is.finite(x)) {
    return("NA")
  }
  if (x < 0.001) {
    return(format(x, scientific = TRUE, digits = 2))
  }
  sprintf("%.3f", x)
}

# Create the six-sample file assignment for one target: three controls and
# three knockdowns. These names are set in CONTROL_SAMPLES and GENE_CONFIG.
make_sample_manifest <- function(config) {
  sample_names <- c(CONTROL_SAMPLES, config$silenced_samples)
  if (anyNA(sample_names) || anyDuplicated(sample_names)) {
    stop(config$target_gene, ": sample names must be present and unique.")
  }
  tibble::tibble(
    sample = sample_names,
    plot_condition = c(rep("Control", length(CONTROL_SAMPLES)),
      rep("CRISPRi knockdown", length(config$silenced_samples))),
    replicate = c(seq_along(CONTROL_SAMPLES), seq_along(config$silenced_samples)),
    forward_file = file.path(BIGWIG_DIR, paste0(sample_names, "_fw.bigwig")),
    reverse_file = file.path(BIGWIG_DIR, paste0(sample_names, "_re.bigwig"))
  )
}

# Weight interval values by the number of bases overlapping each bin. Divide
# by the full bin width, including uncovered bases, rather than covered bases.
weighted_bin_means <- function(start, end, score, bins) {
  if (length(start) != length(end) || length(start) != length(score) ||
      any(!is.finite(c(start, end, score))) || any(end < start) || any(score < 0)) {
    stop("BigWigs must contain finite, nonnegative RPKM interval magnitudes.")
  }
  vapply(seq_len(nrow(bins)), function(i) {
    overlap <- pmax(0, pmin(end, bins$bin_end[i]) - pmax(start, bins$bin_start[i]) + 1)
    sum(score * overlap) / (bins$bin_end[i] - bins$bin_start[i] + 1)
  }, numeric(1))
}

# Import only this genomic window from a track. All tracks use the annotation's
# coordinate system; choose the longest sequence unless an override is given.
read_bigwig_bins <- function(path, config, bins) {
  bw <- rtracklayer::BigWigFile(path)
  lengths <- GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(bw))
  chromosome <- BIGWIG_SEQNAME_OVERRIDE
  if (is.na(chromosome) || !nzchar(chromosome)) {
    finite <- which(is.finite(lengths))
    if (!length(finite)) stop("No finite sequence lengths in BigWig: ", path)
    chromosome <- names(lengths)[finite[which.max(lengths[finite])]]
  }
  if (length(chromosome) != 1L || !chromosome %in% names(lengths) ||
      !is.finite(lengths[chromosome]) || config$region_end > lengths[chromosome]) {
    stop("BigWig sequence/window mismatch: ", path)
  }
  query <- GenomicRanges::GRanges(chromosome,
    IRanges::IRanges(config$region_start, config$region_end))
  signal <- rtracklayer::import(bw, which = query, as = "GRanges")
  score <- if (length(signal)) as.numeric(S4Vectors::mcols(signal)$score) else numeric()
  if (length(score) != length(signal)) stop("Missing BigWig score: ", path)
  weighted_bin_means(GenomicRanges::start(signal), GenomicRanges::end(signal), score, bins)
}

read_annotation <- function() {
  annotation_raw <- readr::read_tsv(ANNOTATION_FILE, show_col_types = FALSE)
  id_col <- find_column(annotation_raw, c("PAO1_ID", "locus_tag", "Locus Tag", "Locus_Tag"))
  gene_col <- find_column(annotation_raw, c("gene_name", "gene", "Gene"), FALSE)
  product_col <- find_column(annotation_raw, c("product_name", "product", "Product"), FALSE)
  start_col <- find_column(annotation_raw, c("start", "Start", "left"))
  end_col <- find_column(annotation_raw, c("end", "End", "right", "stop"))
  strand_col <- find_column(annotation_raw, c("strand", "Strand"))
  operon_col <- find_column(annotation_raw, c("operon", "predicted_operon", "predicted_transcriptional_unit", 
    "transcriptional_unit", "TU", "tu"), FALSE)
  annotation_raw %>% dplyr::transmute(PAO1_ID = as.character(.data[[id_col]]), gene_name = if (!is.na(gene_col)) {
    as.character(.data[[gene_col]])
  } else {
    NA_character_
  }, product_name = if (!is.na(product_col)) {
    as.character(.data[[product_col]])
  } else {
    NA_character_
  }, start_raw = suppressWarnings(as.numeric(.data[[start_col]])), end_raw = suppressWarnings(as.numeric(.data[[end_col]])), 
    strand_raw = as.character(.data[[strand_col]]), operon_id = if (!is.na(operon_col)) {
      as.character(.data[[operon_col]])
    } else {
      NA_character_
    }) %>% dplyr::mutate(start = pmin(start_raw, end_raw), end = pmax(start_raw, end_raw), 
    strand = dplyr::case_when(strand_raw %in% c("+", "forward", "Forward", "1") ~ "+", strand_raw %in% 
      c("-", "reverse", "Reverse", "-1") ~ "-", TRUE ~ strand_raw), gene_label = dplyr::if_else(!is.na(gene_name) & 
      gene_name != "", gene_name, PAO1_ID)) %>% dplyr::filter(!is.na(PAO1_ID), PAO1_ID != 
    "", is.finite(start), is.finite(end)) %>% dplyr::distinct(PAO1_ID, .keep_all = TRUE)
}

save_figure <- function(plot_object, stem, width_mm, height_mm) {
  if (isTRUE(SAVE_PNG)) {
    ggplot2::ggsave(filename = file.path(OUTPUT_DIR, paste0(stem, ".png")), plot = plot_object, 
      width = width_mm, height = height_mm, units = "mm", dpi = PNG_DPI, bg = "white")
  }
  if (isTRUE(SAVE_PDF)) {
    ggplot2::ggsave(filename = file.path(OUTPUT_DIR, paste0(stem, ".pdf")), plot = plot_object, 
      width = width_mm, height = height_mm, units = "mm", device = grDevices::cairo_pdf, 
      bg = "white")
  }
}

# ============================================================
# LOAD ANNOTATION ONCE
# ============================================================
if (!file.exists(ANNOTATION_FILE)) {
  stop("Annotation file not found:\n", ANNOTATION_FILE)
}

annotation <- read_annotation()

# ============================================================
# BUILD ONE COMPLETE PANEL
# ============================================================
# 2. Join regional expression, summarize replicates and construct one panel
build_gene_panel <- function(config, panel_key) {
  message("")
  message("============================================================")
  message("Building panel: ", config$target_gene)
  message("============================================================")
  if (config$region_start >= config$region_end) {
    stop(config$target_gene, ": region_start must be smaller than region_end.")
  }
  if (config$guide_start < config$region_start || config$guide_end > config$region_end) {
    stop(config$target_gene, ": guide coordinates fall outside the configured plotting range.")
  }
  target_annotation <- annotation %>% dplyr::filter(PAO1_ID == config$target_id | tolower(gene_name) == 
    tolower(config$target_gene)) %>% dplyr::slice(1L)
  if (nrow(target_annotation) == 0L) {
    stop("Target gene not found in annotation: ", config$target_gene, " [", config$target_id, "]")
  }
  target_id <- target_annotation$PAO1_ID[[1L]]
  target_start <- target_annotation$start[[1L]]
  target_end <- target_annotation$end[[1L]]
  target_operon <- target_annotation$operon_id[[1L]]
  genes_in_region <- annotation %>% dplyr::filter(end >= config$region_start, start <= config$region_end) %>% 
    dplyr::arrange(start, end) %>% dplyr::mutate(is_target = PAO1_ID == target_id, is_target_operon = !is.na(target_operon) & 
    target_operon != "" & !is.na(operon_id) & operon_id == target_operon, feature_class = dplyr::case_when(is_target ~ 
    "Target gene", is_target_operon ~ "Same annotated operon", TRUE ~ "Neighboring gene"), 
    feature_class = factor(feature_class, levels = c("Target gene", "Same annotated operon", 
      "Neighboring gene")), forward = strand == "+", plot_start = pmax(start, config$region_start), 
    plot_end = pmin(end, config$region_end), plot_midpoint = (plot_start + plot_end)/2, clipped_left = start < 
      config$region_start, clipped_right = end > config$region_end)
  de_file <- find_de_file(config$de_stem)
  message("DE table: ", de_file)
  de_raw <- read_de_table(de_file)
  de_id_col <- find_column(de_raw, c("PAO1_ID", "locus_tag", "Locus_Tag", "Geneid"))
  de_logfc_col <- find_column(de_raw, c("logFC", "log2FoldChange", "log2FC"))
  de_fdr_col <- find_column(de_raw, c("FDR", "padj", "adj.P.Val", "adjusted_p_value"))
  de_table <- de_raw %>% dplyr::transmute(PAO1_ID = as.character(.data[[de_id_col]]), logFC = suppressWarnings(as.numeric(.data[[de_logfc_col]])), 
    adjusted_p = suppressWarnings(as.numeric(.data[[de_fdr_col]]))) %>% dplyr::distinct(PAO1_ID, 
    .keep_all = TRUE)
  regional_de <- genes_in_region %>% dplyr::left_join(de_table, by = "PAO1_ID") %>% dplyr::mutate(midpoint = plot_midpoint, 
    significant = is.finite(adjusted_p) & adjusted_p < DE_SIGNIFICANCE_THRESHOLD, significance = significance_code(adjusted_p), 
    de_direction = dplyr::case_when(!is.finite(logFC) ~ "Not tested", significant & logFC > 
      0 ~ "Significant up", significant & logFC < 0 ~ "Significant down", TRUE ~ "Not significant"))
  regional_de_tested <- regional_de %>% dplyr::filter(is.finite(logFC))
  if (nrow(regional_de_tested) == 0L) {
    stop(config$target_gene, ": no finite regional logFC values.")
  }
  target_de <- regional_de %>% dplyr::filter(is_target) %>% dplyr::slice(1L)
  # 2a. Import each replicate and strand directly from its BigWig.
  coverage_manifest <- make_sample_manifest(config)
  files <- c(coverage_manifest$forward_file, coverage_manifest$reverse_file)
  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) {
    stop(config$target_gene, ": missing RPKM BigWigs:\n", paste(missing_files, collapse = "\n"))
  }
  if (!is.finite(COVERAGE_BIN_SIZE_BP) || COVERAGE_BIN_SIZE_BP < 1 ||
      COVERAGE_BIN_SIZE_BP %% 1 != 0) stop("COVERAGE_BIN_SIZE_BP must be a positive integer.")
  bins <- tibble::tibble(
    bin_start = seq(config$region_start, config$region_end, by = COVERAGE_BIN_SIZE_BP),
    bin_end = pmin(bin_start + COVERAGE_BIN_SIZE_BP - 1L, config$region_end),
    position = (bin_start + bin_end) / 2
  )
  records <- list()
  for (i in seq_len(nrow(coverage_manifest))) {
    for (strand in c("Forward", "Reverse")) {
      path <- if (strand == "Forward") coverage_manifest$forward_file[i] else coverage_manifest$reverse_file[i]
      message("Reading ", config$target_gene, ": ", basename(path))
      values <- read_bigwig_bins(path, config, bins)
      records[[length(records) + 1L]] <- bins %>%
        dplyr::mutate(
          sample = coverage_manifest$sample[i],
          plot_condition = coverage_manifest$plot_condition[i],
          replicate = as.character(coverage_manifest$replicate[i]),
          transcript_strand = strand,
          signed_rpkm = values * if (strand == "Forward") 1 else -1
        ) %>%
        dplyr::select(sample, plot_condition, replicate, transcript_strand,
          bin_start, bin_end, position, signed_rpkm)
    }
  }
  coverage_sample <- dplyr::bind_rows(records) %>%
    dplyr::mutate(
      plot_condition = factor(plot_condition, levels = c("Control", "CRISPRi knockdown")),
      transcript_strand = factor(transcript_strand, levels = c("Forward", "Reverse"))
    )
  # 2b. Summarize biological replicates, then construct the original plot.
  coverage_summary <- coverage_sample %>% dplyr::group_by(plot_condition, transcript_strand, 
    bin_start, bin_end, position) %>% dplyr::summarise(n = dplyr::n(), mean_signed_rpkm = mean(signed_rpkm, 
    na.rm = TRUE), sd_signed_rpkm = stats::sd(signed_rpkm, na.rm = TRUE), lower_signed_rpkm = mean_signed_rpkm - 
    sd_signed_rpkm, upper_signed_rpkm = mean_signed_rpkm + sd_signed_rpkm, .groups = "drop")
  max_abs_cov <- max(abs(c(coverage_sample$signed_rpkm, coverage_summary$mean_signed_rpkm, coverage_summary$lower_signed_rpkm, 
    coverage_summary$upper_signed_rpkm)), na.rm = TRUE)
  if (!is.finite(max_abs_cov) || max_abs_cov <= 0) {
    max_abs_cov <- 1
  }
  de_min_data <- min(regional_de_tested$logFC, na.rm = TRUE)
  de_max_data <- max(regional_de_tested$logFC, na.rm = TRUE)
  de_low <- min(-1, floor(de_min_data - 0.35))
  de_high <- max(1, ceiling(de_max_data + 0.35))
  quantitative_abs_limit <- max(1, abs(c(de_low, de_high)))
  coverage_display_limit <- quantitative_abs_limit * COVERAGE_OVERLAY_FRACTION
  coverage_scale <- coverage_display_limit/max_abs_cov
  coverage_summary <- coverage_summary %>% dplyr::mutate(mean_scaled = mean_signed_rpkm * coverage_scale, 
    lower_scaled = lower_signed_rpkm * coverage_scale, upper_scaled = upper_signed_rpkm * coverage_scale)
  quantitative_low <- min(de_low, -coverage_display_limit)
  quantitative_high <- max(de_high, coverage_display_limit)
  gene_track_y <- quantitative_low - GENE_STRIP_EXTRA_SPACE
  plot_bottom <- gene_track_y - 0.55
  plot_top <- quantitative_high + 0.45
  de_breaks <- pretty(c(de_low, de_high), n = 5)
  de_breaks <- de_breaks[de_breaks >= de_low & de_breaks <= de_high]
  coverage_breaks <- pretty(c(-max_abs_cov, max_abs_cov), n = 5)
  coverage_breaks <- coverage_breaks[coverage_breaks >= -max_abs_cov & coverage_breaks <= max_abs_cov]
  x_span <- config$region_end - config$region_start
  info_x <- config$region_start + 0.015 * x_span
  info_y <- plot_top - 0.04 * (plot_top - plot_bottom)
  target_logfc_text <- if (nrow(target_de) == 1L && is.finite(target_de$logFC[[1L]])) {
    sprintf("%.2f", target_de$logFC[[1L]])
  } else {
    "NA"
  }
  target_fdr_text <- if (nrow(target_de) == 1L) {
    format_p(target_de$adjusted_p[[1L]])
  } else {
    "NA"
  }
  guide_display_html <- gsub(config$target_gene, paste0("<i>", config$target_gene, "</i>"), config$guide_display, 
    fixed = TRUE)
  info_label <- paste0("<i>", config$target_gene, "</i> (", target_id, ")<br>", guide_display_html, 
    "<br>RNA-seq log<sub>2</sub>FC = ", target_logfc_text, "; FDR = ", target_fdr_text)
  info_data <- tibble::tibble(x = info_x, y = info_y, label = info_label)
  condition_colors <- c(Control = CONTROL_COLOR, `CRISPRi knockdown` = KNOCKDOWN_COLOR)
  feature_colors <- c(`Target gene` = TARGET_GENE_COLOR, `Same annotated operon` = TARGET_OPERON_COLOR, 
    `Neighboring gene` = NEIGHBOR_GENE_COLOR)
  regional_up <- regional_de_tested %>% dplyr::filter(de_direction == "Significant up")
  regional_down <- regional_de_tested %>% dplyr::filter(de_direction == "Significant down")
  regional_ns <- regional_de_tested %>% dplyr::filter(!de_direction %in% c("Significant up", 
    "Significant down"))
  guide_legend_data <- tibble::tibble(guide_x = (config$guide_start + config$guide_end)/2, legend_class = "Guide targeting site")
  p <- ggplot2::ggplot() +
  ggplot2::annotate("rect", xmin = max(target_start, config$region_start), 
    xmax = min(target_end, config$region_end), ymin = plot_bottom, ymax = plot_top, fill = TARGET_GENE_COLOR, 
    alpha = 0.045) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, color = "black") + 
    ggplot2::geom_vline(data = guide_legend_data, ggplot2::aes(xintercept = guide_x, color = legend_class), 
      linetype = "dashed", linewidth = 0.6, alpha = 0.95, inherit.aes = FALSE) + {
    if (isTRUE(SHOW_COVERAGE_SD_RIBBON)) {
      ggplot2::geom_ribbon(data = coverage_summary, ggplot2::aes(x = position, ymin = lower_scaled, 
        ymax = upper_scaled, fill = plot_condition, group = interaction(plot_condition, 
         transcript_strand)), alpha = COVERAGE_RIBBON_ALPHA, color = NA)
    }
  } +
  ggplot2::geom_line(data = coverage_summary, ggplot2::aes(x = position, y = mean_scaled, 
    color = plot_condition, group = interaction(plot_condition, transcript_strand)), linewidth = 0.8, 
    alpha = COVERAGE_LINE_ALPHA) +
  ggplot2::geom_segment(data = regional_de_tested, ggplot2::aes(x = midpoint, 
    xend = midpoint, y = 0, yend = logFC), linewidth = 0.38, color = "grey58") +
  ggplot2::geom_point(data = regional_ns, 
    ggplot2::aes(x = midpoint, y = logFC), size = 1.8, color = DE_NS_COLOR) +
  ggplot2::geom_point(data = regional_up, 
    ggplot2::aes(x = midpoint, y = logFC, color = "Upregulated DEG"), size = 2, alpha = 0.92) + 
    ggplot2::geom_point(data = regional_down, ggplot2::aes(x = midpoint, y = logFC, color = "Downregulated DEG"), 
      size = 2, alpha = 0.92) +
  ggplot2::geom_point(data = regional_de_tested %>% dplyr::filter(is_target), 
    ggplot2::aes(x = midpoint, y = logFC), shape = 21, size = 3.7, stroke = 0.85, fill = NA, 
    color = TARGET_GENE_COLOR) + {
    if (isTRUE(SHOW_DE_SIGNIFICANCE_STARS)) {
      ggplot2::geom_text(data = regional_de_tested %>% dplyr::filter(significance != ""), 
        ggplot2::aes(x = midpoint, y = logFC, label = significance, vjust = dplyr::if_else(logFC >= 
         0, -0.75, 1.45)), size = 2.25, show.legend = FALSE)
    }
  } +
  gggenes::geom_gene_arrow(data = genes_in_region, ggplot2::aes(xmin = plot_start, xmax = plot_end, 
    y = gene_track_y, forward = forward, fill = feature_class), arrowhead_height = grid::unit(2.7, 
    "mm"), arrowhead_width = grid::unit(1.3, "mm"), arrow_body_height = grid::unit(3.6, "mm"), 
    color = "black", linewidth = 0.22) + {
    if (isTRUE(SHOW_LOCAL_GENE_LABELS)) {
      ggplot2::geom_text(data = genes_in_region, ggplot2::aes(x = plot_midpoint, y = gene_track_y + 
        0.2, label = gene_label), size = 2.15, fontface = "italic", angle = 45, hjust = 0, 
        vjust = 0, show.legend = FALSE, inherit.aes = FALSE)
    }
  } +
  ggtext::geom_richtext(data = info_data, ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE, 
    hjust = 0, vjust = 1, size = 2.35, family = FONT_FAMILY, fill = "white", color = "black", 
    label.color = "black", label.size = 0.22, label.r = grid::unit(0.1, "lines"), label.padding = grid::unit(c(0.16, 
      0.2, 0.16, 0.2), "lines")) +
  ggplot2::scale_color_manual(limits = c("Control", "CRISPRi knockdown", 
    "Guide targeting site", "Downregulated DEG", "Upregulated DEG"), values = c(Control = CONTROL_COLOR, 
    `CRISPRi knockdown` = KNOCKDOWN_COLOR, `Guide targeting site` = GUIDE_SITE_COLOR, `Downregulated DEG` = DE_DOWN_COLOR, 
    `Upregulated DEG` = DE_UP_COLOR), breaks = c("Control", "CRISPRi knockdown", "Guide targeting site", 
    "Downregulated DEG", "Upregulated DEG"), labels = c("Control mean RPKM (×10³)", "CRISPRi mean RPKM (×10³)", 
    "Guide targeting site", "Downregulated DEG", "Upregulated DEG"), name = NULL, drop = FALSE) + 
    ggplot2::scale_fill_manual(values = c(condition_colors, feature_colors), breaks = c("Control", 
      "CRISPRi knockdown"), labels = c("Control (control1)", "CRISPRi knockdown"), name = NULL, 
      drop = FALSE) +
  ggplot2::scale_x_continuous(limits = c(config$region_start, config$region_end), 
    labels = scales::label_comma(), expand = c(0, 0)) +
  ggplot2::scale_y_continuous(limits = c(plot_bottom, 
    plot_top), breaks = de_breaks, name = expression("RNA-seq log"[2] * " fold change"), sec.axis = ggplot2::sec_axis(transform = ~./coverage_scale/RPKM_DISPLAY_DIVISOR, 
    breaks = coverage_breaks/RPKM_DISPLAY_DIVISOR, labels = scales::label_number(accuracy = 0.1), 
    name = "Mean RPKM coverage (×10³)"), expand = c(0, 0)) +
  ggplot2::labs(title = NULL, 
    subtitle = NULL, x = "PAO1 genomic position [bp]") +
  ggplot2::guides(fill = "none", color = ggplot2::guide_legend(nrow = 1, 
    byrow = TRUE, override.aes = list(linetype = c("solid", "solid", "dashed", "blank", "blank"), 
      shape = c(NA, NA, NA, 16, 16), linewidth = c(1, 1, 0.7, 0, 0), alpha = c(COVERAGE_LINE_ALPHA, 
        COVERAGE_LINE_ALPHA, 0.95, 0.92, 0.92)))) +
  ggplot2::theme_classic(base_family = FONT_FAMILY, 
    base_size = BASE_SIZE) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "italic", 
    size = 10.5, margin = ggplot2::margin(b = 1)), plot.subtitle = ggplot2::element_text(size = 7.3, 
    color = "grey35", margin = ggplot2::margin(b = 2)), axis.title.x = ggplot2::element_text(size = 7.5), 
    axis.title.y = ggplot2::element_text(size = 7.5), axis.text = ggplot2::element_text(size = 6.7), 
    legend.position = "top", legend.justification = "left", legend.text = ggplot2::element_text(size = 6.5), 
    legend.key.width = grid::unit(1.15, "cm"), panel.grid.major.x = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(color = "grey92", 
      linewidth = 0.25), plot.margin = ggplot2::margin(t = 3, r = 8, b = 5, l = 5))
  panel_data <- list(panel = p, regional_de = regional_de, coverage_sample = coverage_sample, 
    coverage_summary = coverage_summary, coverage_manifest = coverage_manifest, config = tibble::tibble(parameter = c("target_gene", 
      "target_id", "guide_display", "guide_name", "guide_sequence", "guide_start", "guide_end", 
      "guide_strand", "region_start", "region_end", "de_file", "coverage_bin_size_bp", "coverage_overlay_fraction", 
      "coverage_line_alpha", "coverage_ribbon_alpha", "rpkm_display_divisor"), value = c(config$target_gene, 
      target_id, config$guide_display, config$guide_name, config$guide_sequence, as.character(config$guide_start), 
      as.character(config$guide_end), config$guide_strand, as.character(config$region_start), 
      as.character(config$region_end), de_file, as.character(COVERAGE_BIN_SIZE_BP), as.character(COVERAGE_OVERLAY_FRACTION), 
      as.character(COVERAGE_LINE_ALPHA), as.character(COVERAGE_RIBBON_ALPHA), as.character(RPKM_DISPLAY_DIVISOR))))
  panel_data
}

# ============================================================
# BUILD ALL FOUR PANELS
# ============================================================
# 3. Read BigWigs and build the four panels

results <- purrr::map(GENE_ORDER, function(gene_key) {
  build_gene_panel(config = GENE_CONFIG[[gene_key]], panel_key = gene_key)
})

names(results) <- GENE_ORDER

gene_panels <- purrr::map(results, "panel")

# ============================================================
# FINAL 2 x 2 FIGURE
# ============================================================
# 4. Display and save the combined and individual figures
combined_figure <- patchwork::wrap_plots(gene_panels[GENE_ORDER], ncol = 2, nrow = 2, byrow = TRUE, 
  guides = "collect") + patchwork::plot_annotation(tag_levels = "A", theme = ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", 
  family = FONT_FAMILY, size = 11))) & ggplot2::theme(legend.position = "top")

if (isTRUE(PRINT_IN_RSTUDIO)) {
  print(combined_figure)
}

save_figure(combined_figure, "Figure5_local_RNAseq_targets_ABCD", FINAL_WIDTH_MM, FINAL_HEIGHT_MM)

# Also save individual panels for manual figure assembly.
for (gene_key in GENE_ORDER) {
  save_figure(results[[gene_key]]$panel, paste0("Figure5_local_RNAseq_", gene_key), 88, 64)
}

# Export every distinct source table once as TSV.
for (gene_key in GENE_ORDER) {
  readr::write_tsv(results[[gene_key]]$regional_de, file.path(OUTPUT_DIR, paste0("Fig5_", gene_key, 
    "_gene_effects.tsv")))
  readr::write_tsv(results[[gene_key]]$coverage_summary, file.path(OUTPUT_DIR, paste0("Fig5_", 
    gene_key, "_coverage.tsv")))
  readr::write_tsv(tibble::as_tibble(results[[gene_key]]$config), file.path(OUTPUT_DIR, paste0("Fig5_", 
    gene_key, "_settings.tsv")))
}

# Sample-level bins support the plotted mean/SD; file assignments document
# which control/knockdown tracks contributed to each regional panel.
coverage_samples <- dplyr::bind_rows(lapply(results, function(x) x$coverage_sample), .id = "region_id")
coverage_files <- dplyr::bind_rows(lapply(results, function(x) x$coverage_manifest), .id = "region_id")
readr::write_tsv(coverage_samples, file.path(OUTPUT_DIR, "Fig5_coverage_samples.tsv.gz"))
readr::write_tsv(coverage_files, file.path(OUTPUT_DIR, "Fig5_coverage_files.tsv"))

# ============================================================
# RSTUDIO OBJECTS
# ============================================================
assign("Figure5_local_RNAseq_combined", combined_figure, envir = .GlobalEnv)

for (gene_key in GENE_ORDER) {
  assign(paste0("Figure5_local_RNAseq_", gene_key), results[[gene_key]]$panel, envir = .GlobalEnv)
}

# ============================================================
# SUMMARY
# ============================================================
message("")

message("============================================================")

message("FIGURE 5 LOCAL RNA-seq PANELS COMPLETE")

message("============================================================")

message("Panels: ", paste(GENE_ORDER, collapse = ", "))

message("Layout: 2 x 2")

message("Final size: ", FINAL_WIDTH_MM, " x ", FINAL_HEIGHT_MM, " mm")

message("Primary y-axis: local gene-level RNA-seq log2FC")

message("Secondary y-axis: mean strand-aware RPKM coverage in x10^3 units ", "(positive = forward; negative = reverse)")

message("Output directory:\n  ", OUTPUT_DIR)

message("RStudio object: Figure5_local_RNAseq_combined")

message("============================================================")


figure_data_files <- c(
  file.path(OUTPUT_DIR, as.vector(outer(paste0("Fig5_", GENE_ORDER),
    c("_gene_effects.tsv", "_coverage.tsv", "_settings.tsv"), paste0))),
  file.path(OUTPUT_DIR, c("Fig5_coverage_samples.tsv.gz", "Fig5_coverage_files.tsv"))
)
message("Underlying figure data:\n",paste0("  ",normalizePath(figure_data_files),collapse="\n"))
