# ============================================================================
# Figure S12A: untreated versus antibiotic-treated knockdown depletion
# ============================================================================
# Biological question: do knockdowns that lose fitness without antibiotic
# show similar, stronger or weaker depletion when an aminoglycoside is present?
# Each point represents one gene in one batch/antibiotic/time comparison.
# More negative log2FC means stronger depletion from the pooled population.
#
# Run from the AMG_cleaned_Git.Rproj working directory:
# source("Figure_scripts/FigS12/AMG_FigS12A_baseline_vs_treated.R")
#
# Inputs:
# - Task 04: prepared Stouffer and MAGeCK interaction calls used for Table S8.
# - Task 11: gene_vulnerability_long.tsv.gz; only full-match severity is used.
# - PAO1_curated.tsv: gene annotations, matched by PAO1 locus ID.
# Existing prepared inputs are sufficient; this script does not fit edgeR,
# recalculate interaction calls, or use sensitivity/EK50 estimates.
#
# Follow the flow below:
# 1. Set the comparisons, paths and display options.
# 2. Define small local functions for reading, summarising and plotting.
# 3. Import the established interaction calls and annotations.
# 4. Pair untreated and treated effects from the same batch and time point.
# 5. Select labels and describe the trend with correlation and a smooth curve.
# 6. Draw and export the combined and individual plots.
# 7. Export seven TSV files containing the underlying figure data.
#
# How to read the scatter:
# - x: untreated versus T0; y: antibiotic-treated versus the same T0.
# - Below y = x: stronger depletion with antibiotic; above: weaker depletion.
# - Colours come from the prepared direct treated-versus-untreated Stouffer
#   calls; darker colours indicate MAGeCK support in the same direction.
# - Distance from the diagonal is descriptive. It need not equal the direct
#   interaction effect because guide filtering/aggregation can differ.
# - The purple curve summarises the pattern; it does not define hit calls.
# T2 uses its matched untreated T2 baseline. Recurrence is not calculated here.
#
# ============================================================================
#
# Input paths: preprocessing_scripts/out_04_AMG_compare_MAGeCK_Stouffer/
# All_gene_comparisons.tsv.gz; preprocessing_scripts/out_11_AMG_gene_vulnerability/
# gene_vulnerability_long.tsv.gz; annotation_input/PAO1_curated.tsv.
# Outputs: PNG/PDF and seven source-data TSVs in Figure_scripts/FigS12/out_AMG_FigS12A/.
#
# R packages used/checked: data.table, dplyr, ggplot2, ggrepel, mgcv, purrr, readr,
# stringr, tibble, tidyr, writexl.

# 1. Input files and plot settings
# All paths are relative to the project root. The first input supplies the
# interaction colours; the second supplies the two coordinates of each point.
METHOD_FILE <- "preprocessing_scripts/out_04_AMG_compare_MAGeCK_Stouffer/All_gene_comparisons.tsv.gz"
VULNERABILITY_FILE <- "preprocessing_scripts/out_11_AMG_gene_vulnerability/gene_vulnerability_long.tsv.gz"
ANNOTATION_FILE <- "annotation_input/PAO1_curated.tsv"
# EDGE_DIR is not read below; point coordinates come from the prepared task-11 table.
EDGE_DIR <- "edgeR_data"

OUTPUT_DIR <- file.path("Figure_scripts", "FigS12", "out_AMG_FigS12A")


# This comparison table keeps each treated sample paired with its own batch
# and untreated baseline. TOB/AMK dose labels retain their original values.
# Most panels use T1; the TOB T2 panel uses untreated T2, not untreated T1.
# The s9_* and mageck_* columns retain the original comparison metadata;
# interaction calls below are selected by the explicit master comparison IDs.
comparison_manifest <- tibble::tribble(~loader_label, ~display_label, ~batch, ~interaction_condition, 
  ~treated_vs_T0_condition, ~baseline_condition, ~mageck_drug, ~mageck_dose, ~mageck_timing, 
  ~comparison_order, ~s9_call_col, ~s9_medLFC_col, ~s9_p_col, "TOB 0.025", "TOB 0.025", "b4", 
  "TOB_0.025_T1 - None_0_T1", "TOB_0.025_T1 - None_0_T0", "None_0_T1 - None_0_T0", "TOB", "0.025", 
  "T1", 1L, "TOB_0.025_call", "TOB_0.025_medLFC", "TOB_0.025_p_stouffer", "TOB 0.1", "TOB 0.1", 
  "b4", "TOB_0.1_T1 - None_0_T1", "TOB_0.1_T1 - None_0_T0", "None_0_T1 - None_0_T0", "TOB", "0.1", 
  "T1", 2L, "TOB_0.1_call", "TOB_0.1_medLFC", "TOB_0.1_p_stouffer", "TOB 0.4", "TOB 0.4", "b1", 
  "TOB_0.4_T1 - None_0_T1", "TOB_0.4_T1 - None_0_T0", "None_0_T1 - None_0_T0", "TOB", "0.4", 
  "T1", 3L, "TOB_0.4_call", "TOB_0.4_medLFC", "TOB_0.4_p_stouffer", "TOB 1", "TOB 1", "b4", "TOB_1_T1 - None_0_T1", 
  "TOB_1_T1 - None_0_T0", "None_0_T1 - None_0_T0", "TOB", "1", "T1", 4L, "TOB_1_call", "TOB_1_medLFC", 
  "TOB_1_p_stouffer", "TOB T2", "TOB 0.4 T2", "b1", "TOB_0.4_T2 - None_0_T2", "TOB_0.4_T2 - None_0_T0", 
  "None_0_T2 - None_0_T0", "TOB", "0.4", "T2", 5L, "TOB_T2_call", "TOB_T2_medLFC", "TOB_T2_p_stouffer", 
  "AMK", "AMK", "b2", "AMK_2_T1 - None_0_T1", "AMK_2_T1 - None_0_T0", "None_0_T1 - None_0_T0", 
  "AMK", "2", "T1", 6L, "AMK_call", "AMK_medLFC", "AMK_p_stouffer")

CONDITIONS_TO_LOAD <- comparison_manifest$loader_label

CONDITIONS_TO_PLOT <- comparison_manifest$display_label

# Display settings: panel arrangement, point transparency and gene labels.
# These affect presentation, not the imported significance decisions.
FACET_NCOL <- 3L

POINT_SIZE <- 0.75

POINT_ALPHA_NONSIG <- 0.2

POINT_ALPHA_SIG <- 0.85

# Limit labels, not plotted genes. Only same-direction hits from both methods
# are eligible; their ranking is applied after the data are joined below.
MAX_LABELS_PER_COMPARISON <- 30L

LABEL_SIZE <- 2.15

LABEL_MAX_TIME <- 5

LABEL_MAX_ITER <- 100000L

# These limits control the displayed area, without filtering the data used
# for the correlations or smooth fits. With USE_FREE_Y = TRUE, Y_LIMITS is
# not imposed on the combined figure; each panel has its own y-axis range.
X_LIMITS <- c(-12, 3)

Y_LIMITS <- c(-12, 5)

USE_FREE_Y <- TRUE

# Require at least 50 paired genes before reporting a panel's trend.
# GAM_K controls the flexibility of the smooth curve. Its mean slope is
# summarised over untreated log2FC -5 to 0 on the grid specified here.
MIN_POINTS_FOR_FIT <- 50L

GAM_K <- 6L

GAM_SLOPE_X_MIN <- -5

GAM_SLOPE_X_MAX <- 0

GAM_SLOPE_GRID_N <- 200L

ANNOTATION_X <- -11.6

ANNOTATION_Y <- 4.6

# Blue = sensitizing, red = protective; dark shades = same-direction support
# from both methods. The 'Stouffer only' display group means no matching
# MAGeCK support: this can include a missing evaluation or a different call.
INTERACTION_COLORS <- c(`Not significant` = "grey60", `Sensitizing: Stouffer only` = "#92C5DE", 
  `Sensitizing: Stouffer + MAGeCK` = "#2166AC", `Protective: Stouffer only` = "#F4A582", `Protective: Stouffer + MAGeCK` = "#B2182B")

PRINT_IN_RSTUDIO <- TRUE

EXPORT_SEPARATE_COMPARISONS <- TRUE

# 2. PACKAGES AND LOCAL FUNCTIONS
# Functions are defined here and used in the later numbered sections.
# Check the existing package requirements without installing/updating them.

packages <- c("data.table", "dplyr", "tidyr", "tibble", "purrr", "stringr", "ggplot2", "ggrepel", 
  "mgcv", "writexl", "readr", "readxl", "R.utils")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) 
  stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(mgcv)
  library(writexl)
  library(readr)
})
OUT_DIR <- file.path(OUTPUT_DIR)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(METHOD_FILE)) {
  stop("Master Supplementary Table S8 was not found:\n", METHOD_FILE)
}
# Standardise blank annotations and missing values before matching records.
clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA", "NaN")] <- NA_character_
  x
}
# Retained utility from the original script; dose_key is not called below.
dose_key <- function(x) {
  x0 <- trimws(as.character(x))
  z <- suppressWarnings(as.numeric(x0))
  out <- x0
  ok <- is.finite(z)
  out[ok] <- format(z[ok], scientific = FALSE, trim = TRUE, digits = 15)
  out[ok] <- sub("\\.?0+$", "", out[ok])
  out
}
# When repeated records contain annotations, retain the first nonblank value.
first_nonempty <- function(x) {
  x <- clean_text(x)
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_character_)
  }
  x[[1L]]
}
# Extract one batch/condition from the shared task-11 table. Its 'severity'
# column is the median adjusted log2FC of full-match guides for each gene.
# If a locus appears more than once, summarise its imported values by median
# so the returned table contains one row per PAO1_ID. No guide fit is rerun.
# Rename the effect to baseline_LFC or treated_LFC for the later join.
load_gene_condition <- function(batch_id, condition_name, value_name) {
  result_file <- VULNERABILITY_FILE
  dat <- data.table::fread(VULNERABILITY_FILE,data.table=FALSE) %>% tibble::as_tibble() %>%
    dplyr::filter(batch==batch_id) %>% dplyr::rename(medLFC=severity)
  required_cols <- c("PAO1_ID", "condition", "medLFC")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    stop("Missing columns in ", result_file, ": ", paste(missing_cols, collapse = ", "))
  }
  dat <- dat %>% mutate(PAO1_ID = clean_text(PAO1_ID), condition = trimws(as.character(condition)), 
    medLFC = as.numeric(medLFC)) %>% filter(condition == condition_name, !is.na(PAO1_ID))
  if (nrow(dat) == 0L) {
    stop("Condition was not found in ", result_file, ":\n  ", condition_name)
  }
  annotation_cols <- intersect(c("gene_name", "gene_label", "category_collapsed"), names(dat))
  out <- dat %>% arrange(PAO1_ID) %>% group_by(PAO1_ID) %>% summarise(gene_name = if ("gene_name" %in% 
    annotation_cols) {
    first_nonempty(gene_name)
  } else {
    NA_character_
  }, gene_label = if ("gene_label" %in% annotation_cols) {
    first_nonempty(gene_label)
  } else {
    NA_character_
  }, category_collapsed = if ("category_collapsed" %in% annotation_cols) {
    first_nonempty(category_collapsed)
  } else {
    NA_character_
  }, value = stats::median(medLFC, na.rm = TRUE), .groups = "drop")
  names(out)[names(out) == "value"] <- value_name
  out
}
# Correlate finite paired effects only; return NA if too few pairs remain.
# Spearman correlation asks whether genes retain a similar rank order of
# depletion, allowing a relationship that is not a straight line.
safe_cor_test <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10L) {
    return(list(estimate = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = method, exact = FALSE))
  list(estimate = unname(test$estimate), p_value = test$p.value)
}
# Summarise each panel using ALL genes with both coordinates, regardless of
# interaction significance. The GAM is a smooth descriptive trend fitted by
# REML. Its average slope over x = -5 to 0 describes how treated depletion
# changes with baseline depletion; it is not an interaction significance test.
# Return both the displayed curve coordinates and the panel statistics.
fit_panel_models <- function(df_panel, panel_name) {
  df_fit <- df_panel %>% filter(is.finite(baseline_LFC), is.finite(treated_LFC))
  if (nrow(df_fit) < MIN_POINTS_FOR_FIT) {
    return(list(prediction = tibble(), statistics = tibble(display_label = panel_name, 
      n_genes = nrow(df_fit), spearman_rho = NA_real_, spearman_p = NA_real_, gam_mean_slope_m5_0 = NA_real_)))
  }
  spearman <- safe_cor_test(df_fit$baseline_LFC, df_fit$treated_LFC, "spearman")
  fit_gam <- mgcv::gam(treated_LFC ~ s(baseline_LFC, k = GAM_K), data = df_fit, method = "REML")
  slope_grid <- tibble(baseline_LFC = seq(GAM_SLOPE_X_MIN, GAM_SLOPE_X_MAX, length.out = GAM_SLOPE_GRID_N))
  slope_y <- as.numeric(stats::predict(fit_gam, newdata = slope_grid))
  gam_mean_slope <- mean(diff(slope_y)/diff(slope_grid$baseline_LFC), na.rm = TRUE)
  x_min <- max(X_LIMITS[[1L]], min(df_fit$baseline_LFC, na.rm = TRUE))
  x_max <- min(X_LIMITS[[2L]], max(df_fit$baseline_LFC, na.rm = TRUE))
  prediction <- tibble(baseline_LFC = seq(x_min, x_max, length.out = 300L)) %>% mutate(display_label = panel_name, 
    y_gam = as.numeric(stats::predict(fit_gam, newdata = .)))
  statistics <- tibble(display_label = panel_name, n_genes = nrow(df_fit), spearman_rho = spearman$estimate, 
    spearman_p = spearman$p_value, gam_mean_slope_m5_0 = gam_mean_slope)
  list(prediction = prediction, statistics = statistics)
}
# Format correlation p values and place the statistical summaries on panels.
format_p <- function(x) {
  if (!is.finite(x)) {
    return("NA")
  }
  if (x < 2e-16) {
    return("<2e-16")
  }
  formatC(x, format = "e", digits = 2)
}
make_panel_annotation <- function(statistics) {
  statistics %>% mutate(annotation = paste0("Spearman rho=", sprintf("%.2f", spearman_rho), 
    ", p=", vapply(spearman_p, format_p, character(1)), "\nGAM slope[-5,0]=", sprintf("%.2f", 
      gam_mean_slope_m5_0)), annotation_x = ANNOTATION_X, annotation_y = ANNOTATION_Y)
}
# Save the same plot as a 300-dpi PNG and a PDF, at the supplied dimensions.
save_plot_png_pdf <- function(plot_object, filename, width, height) {
  ggplot2::ggsave(filename = file.path(OUT_DIR, paste0(filename, ".png")), plot = plot_object, 
    width = width, height = height, units = "in", dpi = 300, bg = "white")
  ggplot2::ggsave(filename = file.path(OUT_DIR, paste0(filename, ".pdf")), plot = plot_object, 
    width = width, height = height, units = "in", device = grDevices::cairo_pdf, bg = "white")
}
# 3. IMPORT THE ESTABLISHED INTERACTION CALLS
# This table already contains direct treated-versus-untreated comparisons.
# Keep those calls/p values and attach curated annotations by PAO1 locus ID.
# Required-column and duplicate checks protect the subsequent gene matching.
master_calls <- data.table::fread(METHOD_FILE,data.table=FALSE) %>% tibble::as_tibble() %>% dplyr::rename(medLFC=medLFC_stouffer)
annotation <- data.table::fread(ANNOTATION_FILE,data.table=FALSE)
master_calls$category <- annotation$category[match(master_calls$PAO1_ID,annotation$PAO1_ID)]
required <- c("comparison_id", "PAO1_ID", "gene_name", "category", "medLFC", "p_stouffer", 
  "stouffer_call", "stouffer_evaluated", "gene_lfc", "neg_fdr", "pos_fdr", "mageck_call", 
  "mageck_evaluated")
missing <- setdiff(required, names(master_calls))
if (length(missing)) 
  stop("Master Table S8 is missing: ", paste(missing, collapse = ", "))
if (anyDuplicated(master_calls[, c("comparison_id", "PAO1_ID")])) 
  stop("Duplicate master Table S8 gene-condition records.")
# Match the six plotted panels to the exact comparison IDs in the input.
# These IDs differ from the treated-versus-T0 conditions used on the y-axis.
master_manifest <- tibble::tibble(comparison_id = c("B4_TOB_0.025_T1_vs_None_T1", "B4_TOB_0.1_T1_vs_None_T1", 
  "B1_TOB_0.4_T1_vs_None_T1", "B4_TOB_1_T1_vs_None_T1", "B1_TOB_0.4_T2_vs_None_T2", "B2_AMK_2_T1_vs_None_T1"), 
  panel_label = comparison_manifest$display_label)
if (!all(master_manifest$comparison_id %in% master_calls$comparison_id)) 
  stop("Master Table S8 lacks a required panel comparison.")
master_calls <- master_calls %>% dplyr::inner_join(master_manifest, by = "comparison_id")
# Keep Stouffer effects/calls and MAGeCK effects/calls in separate small tables.
# They will provide colours, method-support flags and labels for each gene.
interaction_calls <- master_calls %>% dplyr::transmute(display_label = panel_label, PAO1_ID = clean_text(PAO1_ID), 
  gene_name = clean_text(gene_name), gene_label = dplyr::coalesce(clean_text(gene_name), 
    clean_text(PAO1_ID)), category_collapsed = clean_text(category), interaction_direction = stouffer_call, 
  interaction_medLFC = medLFC, interaction_p_stouffer = p_stouffer, stouffer_evaluated)
mageck_calls <- master_calls %>% dplyr::transmute(display_label = panel_label, PAO1_ID = clean_text(PAO1_ID), 
  mageck_gene_lfc = gene_lfc, mageck_neg_fdr = neg_fdr, mageck_pos_fdr = pos_fdr, mageck_evaluated, 
  mageck_call = dplyr::if_else(mageck_evaluated, mageck_call, "Not evaluated"))
# 4. BUILD ONE PAIRED GENE TABLE PER PANEL
# Load untreated-versus-T0 and treated-versus-T0 effects from the same batch.
# Join by PAO1_ID, then attach the independent direct interaction calls.
# Both coordinates must be finite for a gene to appear in the scatter.
plot_list <- vector("list", nrow(comparison_manifest))
for (comparison_i in seq_len(nrow(comparison_manifest))) {
  cfg <- comparison_manifest[comparison_i, , drop = FALSE]
  baseline_tbl <- load_gene_condition(cfg$batch, cfg$baseline_condition, "baseline_LFC")
  treated_tbl <- load_gene_condition(cfg$batch, cfg$treated_vs_T0_condition, "treated_LFC")
  calls_i <- interaction_calls %>% dplyr::filter(display_label == cfg$display_label)
# Use available gene names for labels, falling back to locus IDs. A shared
# call requires BOTH methods to call the SAME interaction direction.
# Grey points have no sensitizing/protective Stouffer call in this display;
# stouffer_status preserves whether that method actually evaluated the gene.
  merged_i <- baseline_tbl %>% full_join(treated_tbl, by = "PAO1_ID", suffix = c("_baseline", 
    "_treated")) %>% left_join(calls_i, by = "PAO1_ID") %>% left_join(mageck_calls %>% 
    filter(display_label == cfg$display_label) %>% dplyr::select(-display_label), by = "PAO1_ID") %>% 
    mutate(display_label = cfg$display_label, comparison_order = cfg$comparison_order, 
      gene_name = dplyr::coalesce(gene_name, gene_name_baseline, gene_name_treated), 
      gene_label = dplyr::coalesce(gene_label, gene_label_baseline, gene_label_treated), 
      category_collapsed = dplyr::coalesce(category_collapsed, category_collapsed_baseline, 
       category_collapsed_treated), gene_display = dplyr::case_when(!is.na(gene_name) & 
       gene_name != "" ~ gene_name, !is.na(gene_label) & gene_label != "" ~ gene_label, 
       TRUE ~ PAO1_ID), interaction_direction = dplyr::case_when(interaction_direction %in% 
       c("Sensitizing", "Protective") ~ interaction_direction, TRUE ~ "Not significant"), 
      mageck_call = dplyr::coalesce(mageck_call, "Not evaluated"), stouffer_evaluated = tidyr::replace_na(stouffer_evaluated, 
       FALSE), mageck_evaluated = tidyr::replace_na(mageck_evaluated, FALSE), stouffer_status = dplyr::if_else(stouffer_evaluated, 
       interaction_direction, "Not evaluated"), shared_same_direction = (interaction_direction == 
       "Sensitizing" & mageck_call == "Sensitizing") | (interaction_direction == "Protective" & 
       mageck_call == "Protective"), interaction_display_group = dplyr::case_when(interaction_direction == 
       "Sensitizing" & shared_same_direction ~ "Sensitizing: Stouffer + MAGeCK", interaction_direction == 
       "Sensitizing" ~ "Sensitizing: Stouffer only", interaction_direction == "Protective" & 
       shared_same_direction ~ "Protective: Stouffer + MAGeCK", interaction_direction == 
       "Protective" ~ "Protective: Stouffer only", TRUE ~ "Not significant")) %>% filter(is.finite(baseline_LFC), 
    is.finite(treated_LFC))
  plot_list[[comparison_i]] <- merged_i
}
# Combine all panels and fix the comparison/group order for plotting.
plot_df <- bind_rows(plot_list) %>% mutate(display_label = factor(display_label, levels = CONDITIONS_TO_PLOT), 
  interaction_direction = factor(interaction_direction, levels = c("Not significant", "Sensitizing", 
    "Protective")), interaction_display_group = factor(interaction_display_group, levels = c("Not significant", 
    "Sensitizing: Stouffer only", "Sensitizing: Stouffer + MAGeCK", "Protective: Stouffer only", 
    "Protective: Stouffer + MAGeCK")))
# 5. SELECT LABELS AND SUMMARISE THE OVERALL RELATIONSHIP
# Label only same-direction hits supported by both methods. Within each
# panel, smaller Stouffer p values rank first, then larger absolute direct
# interaction effects. MAX_LABELS_PER_COMPARISON caps this labelled subset.
labels_df <- plot_df %>% filter(shared_same_direction == TRUE) %>% group_by(display_label) %>% 
  arrange(interaction_p_stouffer, desc(abs(interaction_medLFC)), .by_group = TRUE)
if (is.finite(MAX_LABELS_PER_COMPARISON)) {
  labels_df <- labels_df %>% slice_head(n = as.integer(MAX_LABELS_PER_COMPARISON))
}
labels_df <- labels_df %>% ungroup()
# Fit the descriptive models independently for each panel; the labelled
# subset does not determine the fit. Keep curve points/statistics for export.
fit_outputs <- lapply(CONDITIONS_TO_PLOT, function(panel_i) {
  fit_panel_models(plot_df %>% filter(as.character(display_label) == panel_i), panel_i)
})
prediction_df <- bind_rows(lapply(fit_outputs, function(x) {
  x$prediction
})) %>% mutate(display_label = factor(display_label, levels = CONDITIONS_TO_PLOT))
statistics_df <- bind_rows(lapply(fit_outputs, function(x) {
  x$statistics
}))
annotation_df <- make_panel_annotation(statistics_df) %>% mutate(display_label = factor(display_label, 
  levels = CONDITIONS_TO_PLOT))
facet_scales <- if (isTRUE(USE_FREE_Y)) {
  "free_y"
} else {
  "fixed"
}
# 6. DRAW AND EXPORT THE FIGURE
# Solid diagonal: equal untreated/treated depletion (y = x).
# Dotted diagonals: a difference of +1 or -1 log2FC, shown for orientation;
# these are not a replacement for the imported interaction-call criteria.
# The combined figure is printed in RStudio when PRINT_IN_RSTUDIO is TRUE.
p_s12 <- ggplot(plot_df, aes(x = baseline_LFC, y = treated_LFC)) + geom_abline(intercept = 0, 
  slope = 1, linewidth = 0.55, color = "grey25") + geom_abline(intercept = 1, slope = 1, 
  linetype = "dotted", linewidth = 0.4, color = "steelblue") + geom_abline(intercept = -1, 
  slope = 1, linetype = "dotted", linewidth = 0.4, color = "steelblue") + geom_point(data = plot_df %>% 
  filter(interaction_direction == "Not significant"), color = INTERACTION_COLORS[["Not significant"]], 
  alpha = POINT_ALPHA_NONSIG, size = POINT_SIZE) + geom_point(data = plot_df %>% filter(interaction_direction != 
  "Not significant"), aes(color = interaction_display_group), alpha = POINT_ALPHA_SIG, size = POINT_SIZE + 
  0.35) + geom_line(data = prediction_df, aes(x = baseline_LFC, y = y_gam, group = display_label), 
  inherit.aes = FALSE, linewidth = 0.65, color = "purple3") +
  ggrepel::geom_text_repel(data = labels_df, 
  aes(label = gene_display, color = interaction_display_group), size = LABEL_SIZE, fontface = "italic", 
  box.padding = 0.25, point.padding = 0.08, min.segment.length = 0, segment.color = "grey55", 
  segment.alpha = 0.7, segment.size = 0.25, max.overlaps = Inf, max.time = LABEL_MAX_TIME, 
  max.iter = LABEL_MAX_ITER, seed = 1, show.legend = FALSE) + geom_text(data = annotation_df, 
  aes(x = annotation_x, y = annotation_y, label = annotation), inherit.aes = FALSE, hjust = 0, 
  vjust = 1, size = 2.35, lineheight = 0.95) + facet_wrap(~display_label, ncol = FACET_NCOL, 
  scales = facet_scales) + scale_color_manual(values = INTERACTION_COLORS, breaks = c("Sensitizing: Stouffer only", 
  "Sensitizing: Stouffer + MAGeCK", "Protective: Stouffer only", "Protective: Stouffer + MAGeCK"), 
  labels = c("Sensitizing: Stouffer only", "Sensitizing: both methods", "Protective: Stouffer only", 
    "Protective: both methods"), name = "Interaction call") + scale_x_continuous(breaks = c(-12, 
  -10, -5, 0, 3)) + scale_y_continuous(breaks = c(-12, -10, -5, 0, 5)) + coord_cartesian(xlim = X_LIMITS, 
  ylim = if (isTRUE(USE_FREE_Y)) {
    NULL
  } else {
    Y_LIMITS
  }, clip = "off") + labs(x = "Untreated T1 vs T0 (median adjusted log2FC)", y = "Aminoglycoside-treated vs T0 (median adjusted log2FC)", 
  caption = paste0("Purple line: GAM fit. ", "Solid diagonal: y = x; dotted diagonals: y = x +/- 1.")) + 
  theme_minimal(base_size = 11, base_family = "Arial") + theme(strip.text = element_text(face = "bold", 
  size = 10), legend.position = "right", plot.caption = element_text(hjust = 0, size = 8), 
  plot.margin = margin(15, 30, 10, 10))
if (isTRUE(PRINT_IN_RSTUDIO)) {
  print(p_s12)
}
save_plot_png_pdf(p_s12, "Supplementary_Figure_S12A_AMG_baseline_vs_treated_MAGeCK_overlay_TableS9", 
  width = 10.5, height = 6.6)
# Optionally save each comparison separately, using the same paired data,
# colours, chosen labels and fitted curve as the combined figure.
if (isTRUE(EXPORT_SEPARATE_COMPARISONS)) {
  for (panel_i in CONDITIONS_TO_PLOT) {
    plot_i <- plot_df %>% filter(as.character(display_label) == panel_i)
    labels_i <- labels_df %>% filter(as.character(display_label) == panel_i)
    pred_i <- prediction_df %>% filter(as.character(display_label) == panel_i)
    ann_i <- annotation_df %>% filter(as.character(display_label) == panel_i)
    p_i <- ggplot(plot_i, aes(x = baseline_LFC, y = treated_LFC)) + geom_abline(intercept = 0, 
      slope = 1, linewidth = 0.55, color = "grey25") + geom_abline(intercept = 1, slope = 1, 
      linetype = "dotted", linewidth = 0.4, color = "steelblue") + geom_abline(intercept = -1, 
      slope = 1, linetype = "dotted", linewidth = 0.4, color = "steelblue") + geom_point(data = plot_i %>% 
      filter(interaction_direction == "Not significant"), color = INTERACTION_COLORS[["Not significant"]], 
      alpha = POINT_ALPHA_NONSIG, size = POINT_SIZE) + geom_point(data = plot_i %>% filter(interaction_direction != 
      "Not significant"), aes(color = interaction_display_group), alpha = POINT_ALPHA_SIG, 
      size = POINT_SIZE + 0.35) + geom_line(data = pred_i, aes(x = baseline_LFC, y = y_gam), 
      inherit.aes = FALSE, linewidth = 0.65, color = "purple3") +
  ggrepel::geom_text_repel(data = labels_i, 
      aes(label = gene_display, color = interaction_display_group), size = 2.4, fontface = "italic", 
      box.padding = 0.3, point.padding = 0.08, min.segment.length = 0, segment.color = "grey55", 
      segment.alpha = 0.7, segment.size = 0.25, max.overlaps = Inf, max.time = LABEL_MAX_TIME, 
      max.iter = LABEL_MAX_ITER, seed = 1, show.legend = FALSE) + geom_text(data = ann_i, 
      aes(x = annotation_x, y = annotation_y, label = annotation), inherit.aes = FALSE, 
      hjust = 0, vjust = 1, size = 2.6, lineheight = 0.95) + scale_color_manual(values = INTERACTION_COLORS, 
      breaks = c("Sensitizing: Stouffer only", "Sensitizing: Stouffer + MAGeCK", "Protective: Stouffer only", 
       "Protective: Stouffer + MAGeCK"), labels = c("Sensitizing: Stouffer only", "Sensitizing: both methods", 
       "Protective: Stouffer only", "Protective: both methods"), name = "Interaction call") + 
      coord_cartesian(xlim = X_LIMITS, clip = "off") + labs(title = panel_i, x = "Untreated T1 vs T0 (median adjusted log2FC)", 
      y = "Aminoglycoside-treated vs T0 (median adjusted log2FC)") + theme_minimal(base_size = 11, 
      base_family = "Arial") + theme(plot.title = element_text(face = "bold"), legend.position = "right", 
      plot.margin = margin(20, 50, 15, 15))
    safe_panel_name <- gsub("[^A-Za-z0-9]+", "_", panel_i)
    save_plot_png_pdf(p_i, paste0("Supplementary_Figure_S12A_", safe_panel_name), width = 8.5, 
      height = 7)
  }
}
# 7. EXPORT THE UNDERLYING DATA
# Count support among the Stouffer interaction hits that actually appear in
# the scatter. This summary is restricted to genes with both plotted effects.
# Distinguish absent MAGeCK evaluation from evaluated but nonmatching support.
method_support_summary <- plot_df %>% filter(interaction_direction %in% c("Sensitizing", "Protective")) %>% 
  count(display_label, interaction_direction, shared_same_direction, mageck_evaluated, name = "n_interactions") %>% 
  mutate(method_support = ifelse(shared_same_direction, "Stouffer + MAGeCK", dplyr::if_else(mageck_evaluated, 
    "Stouffer only", "MAGeCK not evaluated"))) %>% arrange(display_label, interaction_direction, 
  desc(shared_same_direction))
# Seven source TSVs allow the figure to be inspected without rerunning fits:
# - plot_data: paired effects, interaction calls and display groups per gene.
# - labels: the subset receiving gene labels.
# - model_statistics / model_predictions: panel summaries and smooth curves.
# - method_support: counts of plotted Stouffer hits by MAGeCK support.
# - MAGeCK_calls: the imported MAGeCK evidence for the selected comparisons.
# - comparisons: the original panel/batch/condition mapping.
tables <- list(FigS12A_plot_data=plot_df,FigS12A_labels=labels_df,
  FigS12A_model_statistics=statistics_df,FigS12A_model_predictions=prediction_df,
  FigS12A_method_support=method_support_summary,FigS12A_MAGeCK_calls=mageck_calls,
  FigS12A_comparisons=comparison_manifest)
for(name in names(tables))readr::write_tsv(tables[[name]],file.path(OUT_DIR,paste0(name,".tsv")))
message("Baseline/treated scatter completed: ", OUT_DIR)



# Print the full paths of the source tables at completion for easy retrieval.
figure_data_files <- file.path(OUTPUT_DIR,c("FigS12A_plot_data.tsv","FigS12A_labels.tsv","FigS12A_model_statistics.tsv","FigS12A_model_predictions.tsv","FigS12A_method_support.tsv","FigS12A_MAGeCK_calls.tsv","FigS12A_comparisons.tsv"))
message("Underlying figure data:\n",paste0("  ",normalizePath(figure_data_files,winslash="/"),collapse="\n"))
