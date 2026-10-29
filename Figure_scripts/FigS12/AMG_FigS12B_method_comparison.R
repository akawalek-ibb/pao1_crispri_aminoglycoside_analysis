# ============================================================================
# Figure S12B: agreement between Stouffer and MAGeCK interaction calls
# ============================================================================
# Biological question: how many gene-antibiotic interactions are supported
# by either analysis method, and how many have support from both?
# Each count refers to gene-comparison interactions, not a recurrence call.
#
# Run from the AMG_cleaned_Git.Rproj working directory:
# source("Figure_scripts/FigS12/AMG_FigS12B_method_comparison.R")
#
# Input: All_gene_comparisons.tsv.gz from task 04, containing the prepared
# Stouffer and primary median-normalised MAGeCK calls used for Table S8.
# Existing prepared input is sufficient. No edgeR/MAGeCK analysis, new
# significance threshold, sensitivity calculation or EK50 fit is run here.
#
# Follow the flow below:
# 1. Set paths, displayed comparisons and plot options.
# 2. Load packages and define comparisons and small local functions.
# 3. Read calls and match the selected union of genes across comparisons.
# 4. Classify support, count interactions and record exceptional cases.
# 5. Draw the three-panel comparison and individual category plots.
# 6. Export six TSV files containing the underlying figure data.
#
# How to read the bars:
# - Sensitizing interactions extend left; protective interactions right.
# - 'Both' requires the same interaction direction from both methods.
# - 'Stouffer only'/'MAGeCK only' requires no call from the other method,
#   which must have evaluated the gene in that comparison.
# Opposite directions, ambiguous MAGeCK calls and unevaluated cases are kept
# in the source data but do not enter the three plotted overlap categories.
# This selected hit union is NOT the genome-wide tested-gene denominator.
#
# Outputs: PNG/PDF plots and six TSVs in Figure_scripts/FigS12/out_AMG_FigS12B/.
# ============================================================================
#
# Input path: preprocessing_scripts/out_04_AMG_compare_MAGeCK_Stouffer/
# All_gene_comparisons.tsv.gz. Its primary MAGeCK calls are used directly.
# No normalized-count or sgRNA-summary MAGeCK files are read by this figure.
#
# R packages used/checked: data.table, ggplot2, openxlsx.

# 1. Input files and plot settings
# Paths are relative to the AMG_cleaned_Git.Rproj working directory.
# Both methods' calls are imported from the same prepared comparison table.
METHOD_FILE <- "preprocessing_scripts/out_04_AMG_compare_MAGeCK_Stouffer/All_gene_comparisons.tsv.gz"
OUTPUT_DIR <- file.path("Figure_scripts", "FigS12", "out_AMG_FigS12B")

# The final export block always writes six TSVs.
# EXPORT_AUDIT_TSV is not used to switch those source-data exports on/off.
EXPORT_AUDIT_TSV <- FALSE

# Include the later TOB T2 comparison alongside the five T1 comparisons.
# T2 is shown on its own; it is not added to a recurrent-hit calculation.
INCLUDE_T2 <- TRUE

# Save PNG/PDF images and/or display enabled combined plots in RStudio.
# These presentation settings do not change the imported interaction calls.
SAVE_PLOTS <- TRUE

PRINT_IN_RSTUDIO <- TRUE

PNG_DPI <- 300

FONT_FAMILY <- "Arial"

# Draw the three method-support panels and save each category separately.
# With SPLIT_FREE_X = TRUE, panels have different axis ranges; use the printed
# counts when comparing categories rather than comparing bar lengths.
MAKE_SCENARIO2_SPLIT <- TRUE

SHOW_GROUP_COUNTS <- TRUE

SPLIT_FREE_X <- TRUE

# Colour identifies sensitizing or protective interaction direction.
# method_colors is not used by the plotting code below.
direction_colors <- c(Sensitizing = "#3B6FB6", Protective = "#C65A4A")


method_colors <- c(Stouffer = "#4D4D4D", `MAGeCK RRA` = "#0072B2", `Any method` = "#A6A6A6", `Both methods` = "#009E73")

# 2. PACKAGES, COMPARISONS AND SMALL LOCAL FUNCTIONS
# Keep the existing package requirements; no packages are installed here.

packages <- c("data.table", "ggplot2", "openxlsx")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) 
  stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
})
OUT_DIR <- file.path(OUTPUT_DIR)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
# One row identifies each batch, drug dose and sampling time. plot_order sets
# the displayed order. The *_col fields contain mapping metadata;
# the current TSV is filtered using comparison_ids constructed below.
manifest <- data.table::data.table(analysis_id = c("B4_TOB_0.025_T1", "B4_TOB_0.1_T1", "B1_TOB_0.4_T1", 
  "B4_TOB_1_T1", "B1_TOB_0.4_T2", "B2_AMK_2_T1"), display_label = c("TOB 0.025", "TOB 0.1", 
  "TOB 0.4", "TOB 1", "TOB 0.4 T2", "AMK"), timing = c("T1", "T1", "T1", "T1", "T2", "T1"), 
  stouffer_col = c("TOB_0.025_Stouffer", "TOB_0.1_Stouffer", "TOB_0.4_Stouffer", "TOB_1_Stouffer", 
    "TOB_0.4_T2_Stouffer", "AMK_Stouffer"), mageck_col = c("TOB_0.025_MAGeCK", "TOB_0.1_MAGeCK", 
    "TOB_0.4_MAGeCK", "TOB_1_MAGeCK", "TOB_0.4_T2_MAGeCK", "AMK_MAGeCK"), comparison_col = c("TOB_0.025_comparison", 
    "TOB_0.1_comparison", "TOB_0.4_comparison", "TOB_1_comparison", "TOB_0.4_T2_comparison", 
    "AMK_comparison"), plot_order = 1:6)
if (!INCLUDE_T2) {
  manifest <- manifest[timing == "T1"]
}
DISPLAY_LEVELS <- manifest[order(plot_order), display_label]
# Standardise blank text, then translate existing call labels into a common
# spelling for counting. Keep MAGeCK's ambiguous direction separate from 'none'.
# Missing evaluations are explicitly restored after this text normalisation.
clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA", "NaN")] <- NA_character_
  x
}
normalize_stouffer_call <- function(x) {
  x <- clean_chr(x)
  data.table::fcase(x == "Sensitizing", "sensitizing", x == "Protective", "protective", default = "none")
}
normalize_mageck_call <- function(x) {
  x <- clean_chr(x)
  data.table::fcase(x == "Sensitizing", "sensitizing", x == "Protective", "protective", x == 
    "Ambiguous", "ambiguous", default = "none")
}
# This column-check utility is not called below.
check_cols <- function(x, cols, label) {
  missing <- setdiff(cols, names(x))
  if (length(missing) > 0L) {
    stop(label, " missing required columns: ", paste(missing, collapse = ", "))
  }
}
# Save the same plot in PNG and PDF formats unless SAVE_PLOTS is FALSE.
save_both <- function(p, stub, width, height) {
  if (!SAVE_PLOTS) {
    return(invisible(NULL))
  }
  ggplot2::ggsave(filename = file.path(OUT_DIR, paste0(stub, ".png")), plot = p, width = width, 
    height = height, units = "in", dpi = PNG_DPI, bg = "white")
  ggplot2::ggsave(filename = file.path(OUT_DIR, paste0(stub, ".pdf")), plot = p, width = width, 
    height = height, units = "in", device = grDevices::cairo_pdf, bg = "white")
}
# Retained workbook-formatting utility; the current export uses plain TSVs
# and does not call add_sheet or create an Excel workbook.
add_sheet <- function(wb, name, x) {
  safe_name <- substr(name, 1, 31)
  openxlsx::addWorksheet(wb, safe_name)
  openxlsx::writeData(wb, safe_name, as.data.frame(x))
  openxlsx::freezePane(wb, safe_name, firstRow = TRUE)
  if (nrow(x) > 0L && ncol(x) > 0L) {
    openxlsx::addFilter(wb, safe_name, rows = 1, cols = seq_len(ncol(x)))
    openxlsx::setColWidths(wb, safe_name, cols = seq_len(ncol(x)), widths = "auto")
  }
}
# 3. READ THE PREPARED CALLS AND MATCH GENES ACROSS COMPARISONS
# Select the exact treated-versus-matched-untreated comparisons from task 04.
# No significance calls are derived from plotted counts or recalculated here.
if (!file.exists(METHOD_FILE)) {
  stop("Master Table S8 workbook was not found:\n  ", METHOD_FILE)
}
calls <- data.table::fread(METHOD_FILE)
comparison_ids <- paste0(manifest$analysis_id,"_vs_None_",manifest$timing)
calls <- calls[comparison_id %in% comparison_ids]
# Retain every gene called by either method in at least one selected panel.
# This produces a shared gene list across panels; some of these genes have
# no interaction call in a particular panel. Never use the union size as the
# number of all genome-wide tested genes, or as a false-positive denominator.
union_ids <- sort(unique(calls$PAO1_ID[calls$any_call==TRUE]))
# For each panel, match the shared gene list to its calls by PAO1 locus ID.
# Preserve raw labels, evaluation flags and the original comparison text.
# A gene missing from that panel's input is marked unevaluated, not 'no hit'.
comparison_rows <- vector("list",nrow(manifest))
for(i in seq_len(nrow(manifest))) {
  cfg <- manifest[i]; values <- calls[comparison_id==comparison_ids[i]]
  ix <- match(union_ids,values$PAO1_ID)
  z <- data.table::data.table(analysis_id=cfg$analysis_id,display_label=cfg$display_label,
    plot_order=cfg$plot_order,PAO1_ID=union_ids,stouffer_call_raw=values$stouffer_call[ix],
    mageck_call_raw=values$mageck_call[ix],stouffer_evaluated=values$stouffer_evaluated[ix],
    mageck_evaluated=values$mageck_evaluated[ix],method_comparison_raw=values$method_comparison[ix])
  z[is.na(stouffer_evaluated),stouffer_evaluated:=FALSE]
  z[is.na(mageck_evaluated),mageck_evaluated:=FALSE]
  z[,stouffer_call:=normalize_stouffer_call(stouffer_call_raw)]
  z[,mageck_call:=normalize_mageck_call(mageck_call_raw)]
  z[stouffer_evaluated==FALSE,stouffer_call:="not_evaluated"]
  z[mageck_evaluated==FALSE,mageck_call:="not_evaluated"]
  comparison_rows[[i]] <- z
}
cmp <- data.table::rbindlist(comparison_rows, use.names = TRUE, fill = TRUE)
if (nrow(cmp) == 0L) {
  stop("No rows were recovered from the final Table S8 method-comparison sheet.")
}
# 4. CLASSIFY AGREEMENT WITHOUT CHANGING EITHER METHOD'S CALL
# meta_hit refers to a Stouffer sensitizing/protective call. Shared requires
# matching directions; opposite directions are classified as conflicts.
# A method-only call requires that the other method evaluated the gene and
# returned 'none'. Unevaluated and ambiguous cases are tracked separately.
cmp[, `:=`(meta_hit, stouffer_call %in% c("sensitizing", "protective"))]
cmp[, `:=`(mageck_hit, mageck_call %in% c("sensitizing", "protective"))]
cmp[, `:=`(shared_same_direction, meta_hit & mageck_hit & stouffer_call == mageck_call)]
cmp[, `:=`(direction_conflict, meta_hit & mageck_hit & stouffer_call != mageck_call)]
cmp[, `:=`(overlap_class, data.table::fcase(!stouffer_evaluated | !mageck_evaluated, "Not evaluated", 
  shared_same_direction, "Shared", meta_hit & !mageck_hit & mageck_call == "none", "Stouffer only", 
  !meta_hit & stouffer_call == "none" & mageck_hit, "MAGeCK only", direction_conflict, "Direction conflict", 
  mageck_call == "ambiguous", "MAGeCK ambiguous", default = "No call"))]
# Reconstruct descriptive agreement labels and compare them with the text
# already stored by task 04. A warning flags mismatches for inspection but
# does not replace the imported calls. Detailed rows remain in comparison_calls.
cmp[, `:=`(comparison_expected, data.table::fcase(stouffer_call == "sensitizing" & mageck_call == 
  "sensitizing", "Shared sensitizing", stouffer_call == "protective" & mageck_call == "protective", 
  "Shared protective", stouffer_call == "sensitizing" & mageck_call == "none", "Stouffer only sensitizing", 
  stouffer_call == "protective" & mageck_call == "none", "Stouffer only protective", stouffer_call == 
    "none" & mageck_call == "sensitizing", "MAGeCK only sensitizing", stouffer_call == 
    "none" & mageck_call == "protective", "MAGeCK only protective", direction_conflict, 
  "Direction conflict", stouffer_call == "none" & mageck_call == "ambiguous", "MAGeCK ambiguous", 
  meta_hit & mageck_call == "ambiguous", "Stouffer call; MAGeCK ambiguous", default = "No call"))]
cmp[, `:=`(method_comparison_matches_TableS9, clean_chr(method_comparison_raw) == comparison_expected)]
comparison_mismatches <- cmp[stouffer_evaluated & mageck_evaluated & !is.na(method_comparison_raw) & !method_comparison_matches_TableS9]
if (nrow(comparison_mismatches) > 0L) {
  warning(nrow(comparison_mismatches), " condition-gene rows differ between reconstructed overlap class ", 
    "and the method_comparison text already stored in Table S8. ", "They are exported for audit.")
}
# Summarise calls within the selected union. These broad counts can include
# conflicts/ambiguous cases, so they need not equal the three plotted groups.
union_summary <- cmp[, .(n_rows_in_FigS12B_union = uniqueN(PAO1_ID), n_stouffer = sum(meta_hit, 
  na.rm = TRUE), n_mageck = sum(mageck_hit, na.rm = TRUE), n_any_method = sum(meta_hit | 
  mageck_hit, na.rm = TRUE), n_both_same_direction = sum(shared_same_direction, na.rm = TRUE), 
  n_direction_conflicts = sum(direction_conflict, na.rm = TRUE), n_mageck_ambiguous = sum(mageck_call == 
    "ambiguous", na.rm = TRUE)), by = .(analysis_id, display_label, plot_order)][order(plot_order)]
# Count sensitizing and protective calls separately in each comparison.
# direction_counts records total support by each/any/both methods;
# exclusive_counts separates shared, Stouffer-only and MAGeCK-only calls.
# Reshape these summaries into one row per count for convenient inspection.
directions <- c("sensitizing", "protective")
count_rows <- list()
exclusive_rows <- list()
k <- 0L
for (i in seq_len(nrow(manifest))) {
  cfg <- manifest[i]
  d <- cmp[analysis_id == cfg$analysis_id]
  for (dir_i in directions) {
    k <- k + 1L
    stouffer_dir <- d$stouffer_call == dir_i
    mageck_dir <- d$mageck_call == dir_i
    stouffer_none <- d$stouffer_call == "none"
    mageck_none <- d$mageck_call == "none"
    count_rows[[k]] <- data.table::data.table(analysis_id = cfg$analysis_id, display_label = cfg$display_label, 
      plot_order = cfg$plot_order, direction = dir_i, Stouffer = sum(stouffer_dir, na.rm = TRUE), 
      `MAGeCK RRA` = sum(mageck_dir, na.rm = TRUE), `Any method` = sum(stouffer_dir | 
       mageck_dir, na.rm = TRUE), `Both methods` = sum(stouffer_dir & mageck_dir, na.rm = TRUE))
    exclusive_rows[[k]] <- data.table::data.table(analysis_id = cfg$analysis_id, display_label = cfg$display_label, 
      plot_order = cfg$plot_order, direction = dir_i, Shared = sum(stouffer_dir & mageck_dir, 
       na.rm = TRUE), `Stouffer only` = sum(stouffer_dir & mageck_none, na.rm = TRUE), 
      `MAGeCK only` = sum(mageck_dir & stouffer_none, na.rm = TRUE))
  }
}
direction_counts <- data.table::rbindlist(count_rows)
exclusive_counts <- data.table::rbindlist(exclusive_rows)
direction_long <- data.table::melt(direction_counts, id.vars = c("analysis_id", "display_label", 
  "plot_order", "direction"), variable.name = "call_set", value.name = "n_interactions")
exclusive_long <- data.table::melt(exclusive_counts, id.vars = c("analysis_id", "display_label", 
  "plot_order", "direction"), variable.name = "overlap_class", value.name = "n_interactions")
direction_long[, `:=`(display_label, factor(display_label, levels = DISPLAY_LEVELS))]
direction_long[, `:=`(direction, factor(direction, levels = c("sensitizing", "protective"), labels = c("Sensitizing", "Protective")))]
direction_long[, `:=`(call_set, factor(call_set, levels = c("Stouffer", "MAGeCK RRA", "Any method", "Both methods")))]
exclusive_long[, `:=`(display_label, factor(display_label, levels = DISPLAY_LEVELS))]
exclusive_long[, `:=`(direction, factor(direction, levels = c("sensitizing", "protective"), labels = c("Sensitizing", "Protective")))]
exclusive_long[, `:=`(overlap_class, factor(overlap_class, levels = c("Shared", "Stouffer only", "MAGeCK only")))]
total_long <- data.table::melt(union_summary[, .(analysis_id, display_label, plot_order, Stouffer = n_stouffer, `MAGeCK RRA` = n_mageck, `Any method` = n_any_method, `Both methods` = n_both_same_direction)], 
  id.vars = c("analysis_id", "display_label", "plot_order"), variable.name = "call_set", 
  value.name = "n_interactions")
total_long[, `:=`(display_label, factor(display_label, levels = DISPLAY_LEVELS))]
total_long[, `:=`(call_set, factor(call_set, levels = c("Stouffer", "MAGeCK RRA", "Any method", "Both methods")))]
# Build the three exclusive groups used in the final figures. A conflict or
# ambiguous/unevaluated counterpart satisfies none of these three conditions,
# and therefore contributes no bar segment. No gene is counted twice within
# the same direction and method-support category.
plot_rows <- list()
k <- 0L
for (i in seq_len(nrow(manifest))) {
  cfg <- manifest[i]
  d <- cmp[analysis_id == cfg$analysis_id]
  for (dir_i in c("sensitizing", "protective")) {
    stouffer_dir <- d$stouffer_call == dir_i
    mageck_dir <- d$mageck_call == dir_i
    stouffer_none <- d$stouffer_call == "none"
    mageck_none <- d$mageck_call == "none"
    k <- k + 1L
    plot_rows[[k]] <- data.table::data.table(analysis_id = cfg$analysis_id, display_label = cfg$display_label, 
      plot_order = cfg$plot_order, direction = dir_i, `Stouffer only` = sum(stouffer_dir & 
       mageck_none, na.rm = TRUE), `MAGeCK only` = sum(mageck_dir & stouffer_none, na.rm = TRUE), 
      Both = sum(stouffer_dir & mageck_dir, na.rm = TRUE))
  }
}
plot_counts <- data.table::rbindlist(plot_rows)
plot_long <- data.table::melt(plot_counts, id.vars = c("analysis_id", "display_label", "plot_order", 
  "direction"), variable.name = "call_category", value.name = "n_interactions")
# Negative values place sensitizing counts to the left of zero; protective
# counts are positive. These signs are drawing positions, not negative counts
# or measured log2FC values; axis labels show the absolute interaction counts.
plot_long[, `:=`(signed_count, data.table::fifelse(direction == "sensitizing", -n_interactions, n_interactions))]
plot_long[, `:=`(display_label, factor(display_label, levels = DISPLAY_LEVELS))]
plot_long[, `:=`(direction, factor(direction, levels = c("sensitizing", "protective"), labels = c("Sensitizing", "Protective")))]
plot_long[, `:=`(call_category, factor(call_category, levels = c("Both", "Stouffer only", "MAGeCK only")))]
# Record the excluded direction conflicts and Stouffer hits with ambiguous
# MAGeCK direction. Warnings explain why broad method totals may exceed the
# sum of the displayed groups. All gene-level cases stay in the source TSV.
direction_conflicts <- cmp[stouffer_call %in% c("sensitizing", "protective") & mageck_call %in% c("sensitizing", "protective") & stouffer_call != mageck_call]
mageck_ambiguous_rows <- cmp[mageck_call == "ambiguous" & stouffer_call %in% c("sensitizing", "protective")]
excluded_summary <- cmp[, .(n_direction_conflicts = sum(stouffer_call %in% c("sensitizing", 
  "protective") & mageck_call %in% c("sensitizing", "protective") & stouffer_call != mageck_call, 
  na.rm = TRUE), n_stouffer_hit_mageck_ambiguous = sum(stouffer_call %in% c("sensitizing", 
  "protective") & mageck_call == "ambiguous", na.rm = TRUE)), by = .(analysis_id, display_label, 
  plot_order)][order(plot_order)]
if (nrow(direction_conflicts) > 0L) {
  warning(nrow(direction_conflicts), " gene-condition rows have opposite interaction directions ", 
    "between final Table S8 Stouffer and MAGeCK calls. They are ", "excluded from the three manuscript call categories.")
}
if (nrow(mageck_ambiguous_rows) > 0L) {
  warning(nrow(mageck_ambiguous_rows), " final Stouffer-hit gene-condition rows have ambiguous MAGeCK ", 
    "direction and are excluded from the three manuscript call categories.")
}
# Keep the raw hit totals for each method/direction and their shared subset.
# MAGeCK-to-Stouffer ratios describe relative call counts, not accuracy or
# validation rates; a ratio is NA when the Stouffer denominator is zero.
method_direction_summary <- cmp[, .(n_Stouffer_sensitizing = sum(stouffer_call == "sensitizing", 
  na.rm = TRUE), n_MAGeCK_sensitizing = sum(mageck_call == "sensitizing", na.rm = TRUE), 
  n_shared_sensitizing = sum(stouffer_call == "sensitizing" & mageck_call == "sensitizing", 
    na.rm = TRUE), n_Stouffer_protective = sum(stouffer_call == "protective", na.rm = TRUE), 
  n_MAGeCK_protective = sum(mageck_call == "protective", na.rm = TRUE), n_shared_protective = sum(stouffer_call == 
    "protective" & mageck_call == "protective", na.rm = TRUE)), by = .(analysis_id, display_label, 
  plot_order)][order(plot_order)]
method_direction_summary[, `:=`(MAGeCK_to_Stouffer_sensitizing = data.table::fifelse(n_Stouffer_sensitizing > 
  0, n_MAGeCK_sensitizing/n_Stouffer_sensitizing, NA_real_), MAGeCK_to_Stouffer_protective = data.table::fifelse(n_Stouffer_protective > 
  0, n_MAGeCK_protective/n_Stouffer_protective, NA_real_))]
# 5. DRAW THE MANUSCRIPT PLOTS
# Separate method-support categories into panels, using the same counts.
# Here colour identifies sensitizing/protective direction instead of method.
# SPLIT_FREE_X controls whether each category uses its own count-axis range.
split_long <- data.table::copy(plot_long)
split_long[, `:=`(call_category, factor(as.character(call_category), levels = c("Stouffer only", "MAGeCK only", "Both")))]
p_split_three <- ggplot2::ggplot(split_long, ggplot2::aes(x = signed_count, y = display_label, 
  fill = direction)) +
  ggplot2::geom_col(orientation = "y", width = 0.66) +
  ggplot2::geom_vline(xintercept = 0, 
  linewidth = 0.55, color = "black") + {
  if (isTRUE(SHOW_GROUP_COUNTS)) {
    ggplot2::geom_text(data = split_long[n_interactions > 0], ggplot2::aes(label = n_interactions, 
      hjust = ifelse(signed_count < 0, 1.12, -0.12)), size = 2.8, family = FONT_FAMILY, 
      fontface = "bold", show.legend = FALSE)
  }
} +
  ggplot2::facet_wrap(~call_category, ncol = 1, scales = if (isTRUE(SPLIT_FREE_X)) {
  "free_x"
} else {
  "fixed"
}) +
  ggplot2::scale_fill_manual(values = direction_colors, name = "Interaction direction", 
  drop = FALSE) +
  ggplot2::scale_x_continuous(labels = function(x) {
  abs(x)
}, expand = ggplot2::expansion(mult = c(0.12, 0.12))) +
  ggplot2::labs(x = "Number of gene-drug interactions", 
  y = NULL, title = "Method-specific and shared final interaction calls", subtitle = if (isTRUE(SPLIT_FREE_X)) {
    "Each call-category panel uses its own x-axis scale; numbers give exact interaction counts."
  } else {
    "All call-category panels use the same x-axis scale; numbers give exact interaction counts."
  }) +
  ggplot2::theme_classic(base_family = FONT_FAMILY, base_size = 11) +
  ggplot2::theme(legend.position = "top", 
  legend.title = ggplot2::element_text(face = "bold"), plot.title = ggplot2::element_text(face = "bold"), 
  strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold", 
    size = 10), axis.text.y = ggplot2::element_text(size = 8.5), panel.grid.major.x = ggplot2::element_line(color = "grey90", 
    linewidth = 0.25), plot.margin = ggplot2::margin(5.5, 16, 5.5, 8))
if (PRINT_IN_RSTUDIO && MAKE_SCENARIO2_SPLIT) {
  print(p_split_three)
}
if (MAKE_SCENARIO2_SPLIT) {
  save_both(p_split_three, "FigS12B_MAGeCK_vs_Stouffer_SCENARIO2_three_call_category_panels", 
    8.3, 9.2)
}
# Also prepare a separate plot for each method-support category, using a
# symmetric left/right range within that plot. The range can differ between
# these separate plots; their count labels remain directly comparable.
make_one_category_plot <- function(category_name) {
  z <- data.table::copy(split_long[as.character(call_category) == category_name])
  max_abs <- max(abs(z$signed_count), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs <= 0) {
    max_abs <- 1
  }
  x_lim <- max_abs * 1.18
  ggplot2::ggplot(z, ggplot2::aes(x = signed_count, y = display_label, fill = direction)) + 
    ggplot2::geom_col(orientation = "y", width = 0.66) +
  ggplot2::geom_vline(xintercept = 0, 
    linewidth = 0.55, color = "black") + {
    if (isTRUE(SHOW_GROUP_COUNTS)) {
      ggplot2::geom_text(data = z[n_interactions > 0], ggplot2::aes(label = n_interactions, 
       hjust = ifelse(signed_count < 0, 1.12, -0.12)), size = 3, family = FONT_FAMILY, 
       fontface = "bold", show.legend = FALSE)
    }
  } +
  ggplot2::scale_fill_manual(values = direction_colors, name = "Interaction direction", 
    drop = FALSE) +
  ggplot2::scale_x_continuous(labels = function(x) {
    abs(x)
  }, limits = c(-x_lim, x_lim), expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
  ggplot2::labs(x = "Number of gene-drug interactions", 
    y = NULL, title = category_name) +
  ggplot2::theme_classic(base_family = FONT_FAMILY, 
    base_size = 11) +
  ggplot2::theme(legend.position = "top", plot.title = ggplot2::element_text(face = "bold"), 
    axis.text.y = ggplot2::element_text(size = 9), panel.grid.major.x = ggplot2::element_line(color = "grey90", 
      linewidth = 0.25), plot.margin = ggplot2::margin(5.5, 18, 5.5, 8))
}
p_stouffer_only <- make_one_category_plot("Stouffer only")
p_mageck_only <- make_one_category_plot("MAGeCK only")
p_both_methods <- make_one_category_plot("Both")
if (MAKE_SCENARIO2_SPLIT) {
  save_both(p_stouffer_only, "FigS12B_MAGeCK_vs_Stouffer_SCENARIO2A_Stouffer_only", 7.2, 
    4.8)
  save_both(p_mageck_only, "FigS12B_MAGeCK_vs_Stouffer_SCENARIO2B_MAGeCK_only", 7.2, 4.8)
  save_both(p_both_methods, "FigS12B_MAGeCK_vs_Stouffer_SCENARIO2C_Both_methods", 7.2, 4.8)
}
# 6. EXPORT COUNTS AND GENE-LEVEL SUPPORT
# Keep disagreement rows available in the R session. The exported comparison
# table contains ALL union-gene/comparison records, including these cases.
# The overlap summary separately records conflicts, ambiguous cases, no calls
# within the selected union, and unevaluated comparisons.
disagreements <- cmp[overlap_class %in% c("Stouffer only", "MAGeCK only", "Direction conflict")]
overlap_summary <- cmp[, .(n_rows_in_FigS12B_union = .N, n_shared = sum(overlap_class == "Shared"), 
  n_stouffer_only = sum(overlap_class == "Stouffer only"), n_mageck_only = sum(overlap_class == 
    "MAGeCK only"), n_direction_conflict = sum(overlap_class == "Direction conflict"), 
  n_mageck_ambiguous = sum(overlap_class == "MAGeCK ambiguous"), n_no_call_within_union = sum(overlap_class == 
    "No call"), n_not_evaluated = sum(overlap_class == "Not evaluated")), by = .(analysis_id, 
  display_label, plot_order)][order(plot_order)]
# Six source TSVs:
# - plot_counts: the exact category/direction counts and their drawing signs.
# - method_summary: each method's directional totals, overlap and count ratios.
# - overlap_summary: all agreement classes within the selected gene union.
# - excluded_summary: counts of conflicts and Stouffer-hit/MAGeCK-ambiguous rows.
# - comparison_calls: per-gene raw/normalised calls, flags and agreement checks.
# - comparisons: the panel/batch/timing mapping used to select input records.
tables <- list(FigS12B_plot_counts=plot_long,FigS12B_method_summary=method_direction_summary,
  FigS12B_overlap_summary=overlap_summary,FigS12B_excluded_summary=excluded_summary,
  FigS12B_comparison_calls=cmp,FigS12B_comparisons=manifest)
for(name in names(tables))utils::write.table(as.data.frame(tables[[name]]),
  file.path(OUT_DIR,paste0(name,".tsv")),sep="\t",quote=TRUE,row.names=FALSE,na="NA")
message("Method-comparison plots completed: ", OUT_DIR)



# Print the full paths of these source tables at completion for easy retrieval.
figure_data_files <- file.path(OUTPUT_DIR,c("FigS12B_plot_counts.tsv","FigS12B_method_summary.tsv","FigS12B_overlap_summary.tsv","FigS12B_excluded_summary.tsv","FigS12B_comparison_calls.tsv","FigS12B_comparisons.tsv"))
message("Underlying figure data:\n",paste0("  ",normalizePath(figure_data_files,winslash="/"),collapse="\n"))
