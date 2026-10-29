#!/usr/bin/env Rscript
# Figure S11: six gene-level aminoglycoside interaction volcano plots.
# Run from the project root after preprocessing_scripts/03_AMG_relevant_results.R.
# Reuse its full gene results, raw-Stouffer calls and comparison-specific cutoffs.
# No guide-level testing, p-value combination, BH filtering or recurrence calling
# is repeated. Both significant and nonsignificant evaluated genes are plotted.
# Flow: read prepared results -> check thresholds/calls -> count interactions ->
# choose up to 15 labels PER DIRECTION -> draw six panels -> save source tables.
#
# Inputs: AG_gene_results.tsv.gz, AG_thresholds.tsv, Comparison_manifest.tsv
# and Analysis_settings.tsv under preprocessing_scripts/out_03_AMG_relevant_results/.
# Run: source("Figure_scripts/FigS11/AMG_FigS11_gene_drug_volcanoes.R")
# Outputs: PNG/PDF/SVG, plotted genes, labels, panel counts, thresholds and settings
# TSVs in Figure_scripts/FigS11/out_AMG_FigS11_gene_drug_volcanoes/.
#
# R packages used/checked: data.table, ggplot2, ggrepel, svglite.

# 1. Input/output locations and display settings.
INPUT_DIR <- "preprocessing_scripts/out_03_AMG_relevant_results"
OUTPUT_DIR <- "Figure_scripts/FigS11/out_AMG_FigS11_gene_drug_volcanoes"

# Keep the source script's order: five T1 comparisons plus the matched B1 T2
# comparison. Displaying T2 does not add it to the recurrence definition.
COMPARISON_LABELS <- c(
  B4_TOB_0.025_T1_vs_None_T1 = "TOB 0.025",
  B4_TOB_0.1_T1_vs_None_T1 = "TOB 0.1",
  B1_TOB_0.4_T1_vs_None_T1 = "TOB 0.4",
  B4_TOB_1_T1_vs_None_T1 = "TOB 1",
  B1_TOB_0.4_T2_vs_None_T2 = "TOB 0.4 T2",
  B2_AMK_2_T1_vs_None_T1 = "AMK"
)

# These checks ensure the inputs still match the stated Figure S11 method.
# Changing a cutoff belongs in step 03, followed by regenerating its outputs.
EXPECTED_P_STOUFFER_CUT <- 0.01
EXPECTED_MAD_MULTIPLIER <- 5
EXPECTED_MIN_ABS_MEDLFC <- 1
MAX_LABELS_PER_DIRECTION <- 15L

# This floor only keeps -log10(p) finite for plotting; raw p-values and calls remain unchanged.
PLOT_P_FLOOR <- 1e-300
VOLCANO_COLOURS <- c(Sensitizing = "#2166AC", Protective = "#B2182B", `Not significant` = "#BDBDBD")
FONT_FAMILY <- "Arial"
WIDTH_IN <- 12
HEIGHT_IN <- 8
DPI <- 300
SEED <- 1L
PRINT_IN_RSTUDIO <- TRUE
OVERWRITE <- TRUE

# Small input checks give an informative error if a prerequisite is missing.
require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

read_table <- function(path, columns) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  con <- if (endsWith(path, ".gz")) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  x <- data.table::as.data.table(utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE,
                                                quote = "\"", comment.char = "", na.strings = c("", "NA")))
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(path, " lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  x
}

read_setting <- function(x, name) {
  value <- x$value[x$setting == name]
  if (length(value) != 1L || !is.finite(suppressWarnings(as.numeric(value)))) stop("Missing/invalid setting: ", name)
  as.numeric(value)
}

same_numbers <- function(x, y, tolerance = 1e-10) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  length(x) == length(y) && all(is.finite(x) & is.finite(y)) &&
    all(abs(x - y) <= tolerance * pmax(1, abs(y)))
}

# Write only the figure's source values and settings, without duplicate Excel exports.
write_table <- function(x, path) {
  con <- if (endsWith(path, ".gz")) gzfile(path, "wt") else file(path, "wt")
  on.exit(close(con))
  utils::write.table(as.data.frame(x), con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}

require_packages(c("data.table", "ggplot2", "ggrepel", "svglite"))
if (MAX_LABELS_PER_DIRECTION < 0L || MAX_LABELS_PER_DIRECTION != as.integer(MAX_LABELS_PER_DIRECTION)) stop("Label limit must be a nonnegative integer.")
if (!is.finite(PLOT_P_FLOOR) || PLOT_P_FLOOR <= 0 || PLOT_P_FLOOR > 1) stop("Invalid plotting p-value floor.")

# 2. Read step 03's complete gene results, cutoff table and saved settings.
# Read the complete gene table, not Stouffer_calls.tsv, which contains only called genes.
paths <- c(genes = file.path(INPUT_DIR, "AG_gene_results.tsv.gz"),
           thresholds = file.path(INPUT_DIR, "AG_thresholds.tsv"),
           settings = file.path(INPUT_DIR, "Analysis_settings.tsv"))
genes <- read_table(paths["genes"], c("comparison_id", "PAO1_ID", "gene_name", "medLFC", "p_stouffer", "stouffer_call"))
thresholds <- read_table(paths["thresholds"], c("comparison_id", "n_stouffer_evaluated", "p_stouffer_cutoff",
                                              "median_medLFC", "MAD", "sensitizing_LFC_cutoff_used", "protective_LFC_cutoff_used"))
upstream_settings <- read_table(paths["settings"], c("setting", "value"))
expected <- c(P_STOUFFER_CUT = EXPECTED_P_STOUFFER_CUT, MAD_MULTIPLIER = EXPECTED_MAD_MULTIPLIER,
              MIN_ABS_MEDLFC = EXPECTED_MIN_ABS_MEDLFC)
for (name in names(expected)) {
  if (!same_numbers(read_setting(upstream_settings, name), expected[[name]]))
    stop("Step 03 setting ", name, " differs from the Figure S11 method. Check the analysis and legend together.")
}
if (anyNA(genes[, .(comparison_id, PAO1_ID)]) || anyDuplicated(genes[, .(comparison_id, PAO1_ID)])) stop("Missing or repeated gene/comparison identifiers.")
if (anyNA(thresholds$comparison_id) || anyDuplicated(thresholds$comparison_id)) stop("Missing or repeated threshold comparison IDs.")
if (length(setdiff(names(COMPARISON_LABELS), unique(genes$comparison_id))) ||
    length(setdiff(names(COMPARISON_LABELS), thresholds$comparison_id))) stop("A required S11 comparison is absent.")

# Keep the six selected comparisons and the genes with finite effect/p-value pairs.
data <- genes[comparison_id %in% names(COMPARISON_LABELS) & is.finite(medLFC) & is.finite(p_stouffer)]
if (any(data$p_stouffer < 0 | data$p_stouffer > 1)) stop("Raw p-values must be between zero and one.")
panel <- data.table::copy(thresholds[match(names(COMPARISON_LABELS), comparison_id)])
data[, display_label := factor(unname(COMPARISON_LABELS[comparison_id]), levels = unname(COMPARISON_LABELS))]
panel[, display_label := factor(unname(COMPARISON_LABELS[comparison_id]), levels = unname(COMPARISON_LABELS))]

# 3. Check that the stored calls and cutoffs correspond to these exact input values.
# MAD is the unscaled median absolute deviation (constant = 1). Final cutoffs are
# min(-1, median - 5*MAD) and max(+1, median + 5*MAD), as in the source script.
# This check stops on disagreement; it never replaces a saved call or cutoff.
for (id in names(COMPARISON_LABELS)) {
  x <- data[comparison_id == id]
  t <- panel[comparison_id == id]
  if (!nrow(x)) stop("No evaluated genes for ", id)
  center <- stats::median(x$medLFC)
  spread <- stats::mad(x$medLFC, center = center, constant = 1)
  lower <- if (spread > 0) min(-EXPECTED_MIN_ABS_MEDLFC, center - EXPECTED_MAD_MULTIPLIER * spread) else -EXPECTED_MIN_ABS_MEDLFC
  upper <- if (spread > 0) max(EXPECTED_MIN_ABS_MEDLFC, center + EXPECTED_MAD_MULTIPLIER * spread) else EXPECTED_MIN_ABS_MEDLFC
  if (!same_numbers(c(nrow(x), center, spread, lower, upper, EXPECTED_P_STOUFFER_CUT),
                    c(t$n_stouffer_evaluated, t$median_medLFC, t$MAD, t$sensitizing_LFC_cutoff_used,
                      t$protective_LFC_cutoff_used, t$p_stouffer_cutoff))) stop("Stored thresholds disagree with gene results: ", id)
  checked_call <- ifelse(x$p_stouffer < t$p_stouffer_cutoff & x$medLFC <= t$sensitizing_LFC_cutoff_used, "Sensitizing",
                         ifelse(x$p_stouffer < t$p_stouffer_cutoff & x$medLFC >= t$protective_LFC_cutoff_used, "Protective", "Not significant"))
  if (anyNA(x$stouffer_call) || any(checked_call != x$stouffer_call)) stop("Stored calls disagree with thresholds: ", id)
}

# Use gene names where available, otherwise show the PAO1 locus identifier.
data[, gene_label := ifelse(is.na(gene_name) | !nzchar(trimws(gene_name)) | gene_name == "NA", PAO1_ID, trimws(gene_name))]
data[, `:=`(minus_log10_p = -log10(pmax(p_stouffer, PLOT_P_FLOOR)),
             interaction_class = factor(stouffer_call, levels = names(VOLCANO_COLOURS)),
             significant = stouffer_call != "Not significant", abs_medLFC = abs(medLFC))]

# 4. Count the two interaction directions and prepare each panel's annotation box.
counts <- data[, .(n_all_genes = .N, n_sensitizing = sum(stouffer_call == "Sensitizing"),
                   n_protective = sum(stouffer_call == "Protective")), by = comparison_id]
panel <- merge(panel, counts, by = "comparison_id", all.x = TRUE, sort = FALSE)
panel <- panel[match(names(COMPARISON_LABELS), comparison_id)]
panel[, `:=`(n_significant = n_sensitizing + n_protective,
              lower_5MAD = ifelse(is.finite(MAD) & MAD > 0, median_medLFC - EXPECTED_MAD_MULTIPLIER * MAD, NA_real_),
              upper_5MAD = ifelse(is.finite(MAD) & MAD > 0, median_medLFC + EXPECTED_MAD_MULTIPLIER * MAD, NA_real_))]
panel[, panel_annotation := paste0("Sens = ", n_sensitizing, "   Prot = ", n_protective,
                                    "\nMedian = ", sprintf("%.3f", median_medLFC),
                                    "\nMAD = ", sprintf("%.3f", MAD),
                                    "\n5MAD range = [", sprintf("%.3f", lower_5MAD), ", ", sprintf("%.3f", upper_5MAD), "]",
                                    "\nFinal LFC cuts: <= ", sprintf("%.3f", sensitizing_LFC_cutoff_used),
                                    " / >= ", sprintf("%.3f", protective_LFC_cutoff_used))]

# 5. Select labels separately for sensitizing and protective genes in each panel.
# Smaller raw p-values come first; ties use larger absolute effects, then PAO1_ID.
# The default therefore allows up to 30 labels per panel, not 15 in total.
labels <- data.table::copy(data[significant == TRUE])
data.table::setorderv(labels, c("comparison_id", "interaction_class", "p_stouffer", "abs_medLFC", "PAO1_ID"), c(1, 1, 1, -1, 1))
labels <- labels[, head(.SD, MAX_LABELS_PER_DIRECTION), by = .(comparison_id, display_label, interaction_class)]

# 6. Draw the source script's six-panel layout and colour classes.
# Dotted vertical line: median. Dot-dashed: raw median +/- 5MAD.
# Dashed vertical lines: final effect cutoffs after the +/-1 floor.
# Dashed horizontal line: raw Stouffer p = 0.01.
# The statistics boxes are INSIDE the panels, at their upper-left corners.
p <- ggplot2::ggplot(data, ggplot2::aes(x = medLFC, y = minus_log10_p)) +
  ggplot2::geom_point(ggplot2::aes(colour = interaction_class), size = 1.55, alpha = 0.72) +
  ggplot2::geom_vline(data = panel[is.finite(median_medLFC)], ggplot2::aes(xintercept = median_medLFC),
                      linetype = "dotted", linewidth = 0.30, colour = "grey60") +
  ggplot2::geom_vline(data = panel[is.finite(lower_5MAD)], ggplot2::aes(xintercept = lower_5MAD),
                      linetype = "dotdash", linewidth = 0.38, colour = "grey65") +
  ggplot2::geom_vline(data = panel[is.finite(upper_5MAD)], ggplot2::aes(xintercept = upper_5MAD),
                      linetype = "dotdash", linewidth = 0.38, colour = "grey65") +
  ggplot2::geom_vline(data = panel, ggplot2::aes(xintercept = sensitizing_LFC_cutoff_used),
                      linetype = "dashed", linewidth = 0.55, colour = "grey25") +
  ggplot2::geom_vline(data = panel, ggplot2::aes(xintercept = protective_LFC_cutoff_used),
                      linetype = "dashed", linewidth = 0.55, colour = "grey25") +
  ggplot2::geom_hline(yintercept = -log10(EXPECTED_P_STOUFFER_CUT), linetype = "dashed", linewidth = 0.50, colour = "grey25") +
  ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = paste0("italic('", gsub("'", "", gene_label, fixed = TRUE), "')")),
                         size = 2.15, parse = TRUE, max.overlaps = Inf, box.padding = 0.28, point.padding = 0.10,
                         min.segment.length = 0, seed = SEED, show.legend = FALSE) +
  ggplot2::geom_label(data = panel, ggplot2::aes(x = -Inf, y = Inf, label = panel_annotation),
                       hjust = -0.04, vjust = 1.05, size = 2.05, linewidth = 0.12,
                       fill = "white", alpha = 0.92, inherit.aes = FALSE) +
  ggplot2::facet_wrap(~display_label, ncol = 3, scales = "free") +
  ggplot2::scale_colour_manual(values = VOLCANO_COLOURS, breaks = names(VOLCANO_COLOURS), name = "Interaction") +
  ggplot2::labs(title = "Gene-aminoglycoside interaction severity",
                 subtitle = "Raw Stouffer p < 0.01; median +/- 5 unscaled MAD with minimum |medLFC| >= 1",
                 x = "Median adjusted log2 fold change", y = expression(-log[10] * "(raw " * p[Stouffer] * ")")) +
  ggplot2::theme_classic(base_size = 9, base_family = FONT_FAMILY) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), strip.text = ggplot2::element_text(face = "bold"),
                 legend.position = "bottom", panel.spacing = grid::unit(1.0, "lines"))

# 7. Save the figure and compact source tables. Raw p-values, stored calls and
# comparison IDs remain available alongside their plotted transformations.
outputs <- c(paste0("FigS11_gene_drug_volcanoes", c(".svg", ".pdf", ".png")),
             "FigS11_plotted_values.tsv.gz", "FigS11_panel_summary.tsv", "FigS11_labelled_genes.tsv", "Analysis_settings.tsv")
if (!OVERWRITE && any(file.exists(file.path(OUTPUT_DIR, outputs)))) stop("Output exists; set OVERWRITE <- TRUE.")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
data.table::setorder(data, display_label, PAO1_ID)
settings <- data.table::data.table(setting = c(names(expected), "MAX_LABELS_PER_DIRECTION", "PLOT_P_FLOOR", "SEED",
                                               paste0("input:", names(paths)), paste0("md5:", names(paths)), "R_version"),
                                   value = c(as.character(expected), MAX_LABELS_PER_DIRECTION, PLOT_P_FLOOR, SEED,
                                             unname(paths), unname(tools::md5sum(paths)), R.version.string))
columns <- c("comparison_id", "display_label", "PAO1_ID", "gene_name", "gene_label", "medLFC", "p_stouffer",
             "minus_log10_p", "interaction_class", "significant")
write_table(data[, ..columns], file.path(OUTPUT_DIR, "FigS11_plotted_values.tsv.gz"))
write_table(panel, file.path(OUTPUT_DIR, "FigS11_panel_summary.tsv"))
write_table(labels[, ..columns], file.path(OUTPUT_DIR, "FigS11_labelled_genes.tsv"))
write_table(settings, file.path(OUTPUT_DIR, "Analysis_settings.tsv"))
ggplot2::ggsave(file.path(OUTPUT_DIR, "FigS11_gene_drug_volcanoes.svg"), p, device = svglite::svglite, width = WIDTH_IN, height = HEIGHT_IN, bg = "white")
ggplot2::ggsave(file.path(OUTPUT_DIR, "FigS11_gene_drug_volcanoes.pdf"), p, device = grDevices::cairo_pdf, width = WIDTH_IN, height = HEIGHT_IN, bg = "white")
ggplot2::ggsave(file.path(OUTPUT_DIR, "FigS11_gene_drug_volcanoes.png"), p, width = WIDTH_IN, height = HEIGHT_IN, dpi = DPI, bg = "white")
if (PRINT_IN_RSTUDIO) print(p)
FigS11_plot <- p
message("Saved Fig. S11: ", nrow(data), " gene/comparison records across six panels; ", nrow(labels), " labels.\n", OUTPUT_DIR)
