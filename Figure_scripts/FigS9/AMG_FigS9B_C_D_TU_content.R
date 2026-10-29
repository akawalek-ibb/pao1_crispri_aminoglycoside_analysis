#!/usr/bin/env Rscript
# Figure S9B-D: severity, sensitivity and EK50 in relation to TU gene content.
# Run from the project root after preprocessing steps 11 and 12.
# Reuse the same B1-B4 untreated T1 consensus as Table S5; no guide curves are refitted.
# Negative fitness means raw Stouffer p < 0.05 and medLFC <= -1 in at least one
# of those four T1 comparisons. B1 T2 is excluded.
# Flow: read consensus and fitness calls -> count genes in each TU -> select negative
# fitness genes -> test associations at gene and TU levels -> draw and save panels.
#
# Inputs: preprocessing_scripts/out_12_AMG_MH_fitness/Fig2G_input.tsv and
# MH_fitness.tsv, plus annotation_input/PAO1_curated.tsv.
# Run: source("Figure_scripts/FigS9/AMG_FigS9B_C_D_TU_content.R")
# Outputs: PNG/PDF/SVG panels, gene/TU-level plotted data, tests and settings TSVs
# in Figure_scripts/FigS9/out_AMG_FigS9B_C_D_TU_content/.
#
# R packages used/checked: data.table, ggplot2, ggrepel, svglite.

# 1. Inputs and analysis/display settings.
INPUT_FILE <- "preprocessing_scripts/out_12_AMG_MH_fitness/Fig2G_input.tsv"
FITNESS_FILE <- "preprocessing_scripts/out_12_AMG_MH_fitness/MH_fitness.tsv"
ANNOTATION_FILE <- "annotation_input/PAO1_curated.tsv"
OUTPUT_DIR <- "Figure_scripts/FigS9/out_AMG_FigS9B_C_D_TU_content"

# Count both design groups, including fitness/AMR-associated genes.
# Set "Nr_E" for previously classified essential genes only; that changes groups and tests.
TU_COUNT_DEFINITION <- "design_call"
CONTROL_PREFIXES <- c("B1_T1_vs_T0", "B2_T1_vs_T0", "B3_T1_vs_T0", "B4_T1_vs_T0")
P_STOUFFER_CUT <- 0.05
LFC_CUT <- 1
MIN_SEVERITY_ESTIMATES <- 1L
MIN_SENSITIVITY_ESTIMATES <- 1L
MIN_EK50_ESTIMATES <- 2L
TU_BIN_CAP <- 6L
N_GENE_LABELS <- 30L
FONT_FAMILY <- "Arial"
WIDTH_IN <- 9
HEIGHT_IN <- 6
DPI <- 300
SEED <- 123L
PRINT_IN_RSTUDIO <- TRUE
OVERWRITE <- TRUE

# Small checks keep missing or mismatched input tables from producing a misleading figure.
require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

read_table <- function(path, columns) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  x <- data.table::fread(path, na.strings = c("", "NA"))
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(path, " lacks columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyNA(x$PAO1_ID) || anyDuplicated(x$PAO1_ID)) stop("Missing or repeated PAO1_ID in ", path)
  x
}

# Count genes across the FULL annotation, before selecting genes for any panel.
# The original scripts used both main_essential and other_essential design groups.
# This includes fitness/AMR-associated genes; it is not an essential-only count.
# Nr_E is available if the intended count is genes previously called essential.
attach_TU_counts <- function(x, annotation_file, definition) {
  if (!definition %in% c("design_call", "Nr_E")) stop("TU_COUNT_DEFINITION must be design_call or Nr_E.")
  annotation <- read_table(annotation_file, c("PAO1_ID", "TU", "design_call", "Nr_E"))
  annotation[, TU := trimws(as.character(TU))]
  annotation[toupper(TU) %in% c("", "NA", "NAN", "NONE", "NULL"), TU := NA_character_]
  if (definition == "design_call") {
    annotation[, counted := design_call %in% c("main_essential", "other_essential")]
  } else {
    nr <- suppressWarnings(as.numeric(annotation$Nr_E))
    annotation[, counted := is.finite(nr) & nr >= 1]
  }
  counts <- annotation[!is.na(TU), .(essential_genes_in_operon = sum(counted)), by = TU]
  index <- match(x$PAO1_ID, annotation$PAO1_ID)
  if (anyNA(index)) stop("Consensus contains genes absent from the annotation.")
  x <- data.table::copy(x)
  x[, TU := annotation$TU[index]]
  x[, essential_genes_in_operon := counts$essential_genes_in_operon[match(TU, counts$TU)]]
  x
}

# Save ordinary TSVs so every plotted point and reported statistic can be inspected.
write_tables <- function(tables, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  for (name in names(tables)) {
    utils::write.table(as.data.frame(tables[[name]]), file.path(directory, paste0(name, ".tsv")),
                      sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
  }
}

save_figure <- function(plot, stem, width, height) {
  ggplot2::ggsave(file.path(OUTPUT_DIR, paste0(stem, ".svg")), plot, device = svglite::svglite,
                  width = width, height = height, bg = "white")
  ggplot2::ggsave(file.path(OUTPUT_DIR, paste0(stem, ".pdf")), plot, device = grDevices::cairo_pdf,
                  width = width, height = height, bg = "white")
  ggplot2::ggsave(file.path(OUTPUT_DIR, paste0(stem, ".png")), plot,
                  width = width, height = height, dpi = DPI, bg = "white")
  if (PRINT_IN_RSTUDIO) print(plot)
}

# Record the input versions through file checksums, plus the settings used for this run.
run_settings <- function(settings, files) {
  data.table::data.table(setting = c(names(settings), paste0("input:", names(files)),
                                    paste0("md5:", names(files)), "R_version"),
                         value = c(as.character(settings), unname(files),
                                   unname(tools::md5sum(files)), R.version.string))
}

# Do not turn a missing metric into zero when summarising a TU.
finite_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

# Preserve the original minimum-data rules. Inadequate data produce a missing test.
association_tests <- function(x, metric, level, panel) {
  kruskal <- if (nrow(x) >= 3L && data.table::uniqueN(x$TU_bin) >= 2L)
    tryCatch(stats::kruskal.test(x[[metric]] ~ droplevels(x$TU_bin)), error = function(e) NULL) else NULL
  correlation <- if (nrow(x) >= 10L)
    tryCatch(stats::cor.test(x$essential_genes_in_operon, x[[metric]], method = "spearman", exact = FALSE), error = function(e) NULL) else NULL
  value <- function(test, field) if (is.null(test)) NA_real_ else as.numeric(test[[field]])[1L]
  data.table::data.table(panel = panel, analysis_level = level, metric = metric,
                         test = c("Kruskal-Wallis", "Spearman"), n = nrow(x),
                         statistic = c(value(kruskal, "statistic"), value(correlation, "estimate")),
                         df = c(value(kruskal, "parameter"), NA_real_),
                         p_value = c(value(kruskal, "p.value"), value(correlation, "p.value")))
}

test_label <- function(x) {
  kw <- x[test == "Kruskal-Wallis"]
  sp <- x[test == "Spearman"]
  if (!is.finite(kw$p_value)) first <- "KW: not estimable" else
    first <- paste0("KW chi2=", round(kw$statistic, 2), " (df=", kw$df, "), p=", signif(kw$p_value, 3))
  if (!is.finite(sp$p_value)) second <- "Spearman: not estimable" else
    second <- paste0("Spearman rho=", round(sp$statistic, 3), ", p=", signif(sp$p_value, 3))
  paste(first, second, sep = " | ")
}

require_packages(c("data.table", "ggplot2", "ggrepel", "svglite"))
if (TU_BIN_CAP < 1L || TU_BIN_CAP != as.integer(TU_BIN_CAP)) stop("TU_BIN_CAP must be a positive integer.")

# 2. Read the annotated T1 consensus and the per-comparison fitness results.
# Use the four individual T1 effects/p-values, not the five-comparison S5 union.
x <- read_table(INPUT_FILE, c("PAO1_ID", "gene_name", "control_median_severity", "control_median_sensitivity",
                             "control_median_EK50", "n_control_severity_estimates",
                             "n_control_sensitivity_estimates", "n_control_EK50_estimates"))
fields <- c(paste0(CONTROL_PREFIXES, "_medLFC"), paste0(CONTROL_PREFIXES, "_p_stouffer"))
fitness <- read_table(FITNESS_FILE, c("PAO1_ID", fields))
hits <- vapply(CONTROL_PREFIXES, function(prefix) {
  effect <- fitness[[paste0(prefix, "_medLFC")]]
  p <- fitness[[paste0(prefix, "_p_stouffer")]]
  is.finite(effect) & is.finite(p) & p < P_STOUFFER_CUT & effect <= -LFC_CUT
}, logical(nrow(fitness)))
negative_ids <- fitness$PAO1_ID[rowSums(hits) > 0L]
x <- attach_TU_counts(x, ANNOTATION_FILE, TU_COUNT_DEFINITION)
x[, `:=`(severity = as.numeric(control_median_severity),
          sensitivity = as.numeric(control_median_sensitivity), EK50 = as.numeric(control_median_EK50))]

# 3. Restrict all panels to negative-fitness genes with an available severity estimate.
# Missing TU identifiers are excluded; they are not combined into an artificial TU.
x <- x[PAO1_ID %in% negative_ids & n_control_severity_estimates >= MIN_SEVERITY_ESTIMATES & !is.na(TU)]
x[n_control_sensitivity_estimates < MIN_SENSITIVITY_ESTIMATES, sensitivity := NA_real_]
x[n_control_EK50_estimates < MIN_EK50_ESTIMATES, EK50 := NA_real_]
x[, label_txt := ifelse(is.na(gene_name) | !nzchar(trimws(gene_name)), PAO1_ID, gsub("'", "", gene_name, fixed = TRUE))]
data.table::setorder(x, PAO1_ID)
bin_levels <- c(as.character(0:(TU_BIN_CAP - 1L)), paste0(TU_BIN_CAP, "+"))
make_bins <- function(n) factor(ifelse(n >= TU_BIN_CAP, paste0(TU_BIN_CAP, "+"), as.character(as.integer(n))), levels = bin_levels)
x[, TU_bin := make_bins(essential_genes_in_operon)]

# 4. Make a second table with ONE median per TU, using the selected negative-fitness
# genes in that TU. These tests reduce repeated counting of genes from the same TU.
# The number of counted genes in the TU still comes from the full annotation.
tu <- x[, .(n_selected_genes = .N, essential_genes_in_operon = max(essential_genes_in_operon),
             severity = finite_median(severity), sensitivity = finite_median(sensitivity), EK50 = finite_median(EK50)), by = TU]
tu[, TU_bin := make_bins(essential_genes_in_operon)]

# 5. Analyse the three metrics separately because each has a different number of
# estimable values. Spearman uses the original count; Kruskal-Wallis uses the bins.
metric_names <- c("severity", "sensitivity", "EK50")
panels <- c("FigS9B", "FigS9C", "FigS9D")
axis_labels <- c("Severity (median full-match-guide adjusted log2 fold change)",
                 "Sensitivity (predicted efficacy at sustained log2FC <= -1)", "EK50 (predicted guide efficacy)")
count_label <- if (TU_COUNT_DEFINITION == "design_call") "Essential/fitness-associated genes per TU (binned)" else
  "Previously classified essential genes per TU (binned)"
plots <- list()
gene_values <- TU_values <- all_tests <- labels_used <- list()
for (i in seq_along(metric_names)) {
  metric <- metric_names[i]
  panel <- panels[i]
  genes <- data.table::copy(x[is.finite(get(metric)) & !is.na(TU_bin)])
  units <- data.table::copy(tu[is.finite(get(metric)) & !is.na(TU_bin)])
  if (metric == "sensitivity") {
    genes <- genes[sensitivity > 0]
    units <- units[sensitivity > 0]
  }
  if (!nrow(genes)) stop("No estimable genes for ", panel)
  gene_test <- association_tests(genes, metric, "gene", panel)
  TU_test <- association_tests(units, metric, "TU", panel)
  labels <- head(genes[order(get(metric))], N_GENE_LABELS)

  # The boxes and points show genes. The second line of statistics uses TU medians.
  # Keep the source panels' flipped coordinates and reversed metric axes.
  p <- ggplot2::ggplot(genes, ggplot2::aes(x = TU_bin, y = .data[[metric]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.4) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15, height = 0, seed = SEED), alpha = 0.20, size = 0.8) +
    ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = label_txt), size = 2,
                           max.overlaps = Inf, fontface = "italic", seed = SEED) +
    ggplot2::coord_flip() + ggplot2::scale_y_reverse() +
    ggplot2::theme_minimal(base_size = 12, base_family = FONT_FAMILY) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), axis.title = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold")) +
    ggplot2::labs(title = paste0(sub("Fig", "", panel), "  ", metric, " and TU gene content"),
                   subtitle = paste0("Gene level: ", test_label(gene_test), "\nTU level: ", test_label(TU_test)),
                   x = count_label, y = axis_labels[i])
  genes[, panel := panel]
  units[, panel := panel]
  labels[, panel := panel]
  gene_values[[metric]] <- genes
  TU_values[[metric]] <- units
  labels_used[[metric]] <- labels
  all_tests[[metric]] <- data.table::rbindlist(list(gene_test, TU_test))
  plots[[panel]] <- p
}

# 6. Save the three figures, the exact per-panel gene/TU data and both sets of tests.
# Unrelated severity-versus-sensitivity plots and exploratory regression models
# from the original combined script are not repeated here.
stems <- c("FigS9B_severity_TU_content", "FigS9C_sensitivity_TU_content", "FigS9D_EK50_TU_content")
table_names <- c("FigS9BCD_gene_values", "FigS9BCD_TU_values", "FigS9BCD_labels", "FigS9BCD_statistics", "Analysis_settings")
outputs <- c(as.vector(outer(stems, c(".svg", ".pdf", ".png"), paste0)), paste0(table_names, ".tsv"))
if (!OVERWRITE && any(file.exists(file.path(OUTPUT_DIR, outputs)))) stop("Output exists; set OVERWRITE <- TRUE.")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
settings <- run_settings(c(TU_COUNT_DEFINITION = TU_COUNT_DEFINITION, CONTROL_PREFIXES = paste(CONTROL_PREFIXES, collapse = ";"),
                           P_STOUFFER_CUT = P_STOUFFER_CUT, LFC_CUT = LFC_CUT, TU_BIN_CAP = TU_BIN_CAP,
                           MIN_SEVERITY_ESTIMATES = MIN_SEVERITY_ESTIMATES, MIN_SENSITIVITY_ESTIMATES = MIN_SENSITIVITY_ESTIMATES,
                           MIN_EK50_ESTIMATES = MIN_EK50_ESTIMATES, SEED = SEED),
                         c(consensus = INPUT_FILE, fitness = FITNESS_FILE, annotation = ANNOTATION_FILE))
write_tables(list(FigS9BCD_gene_values = data.table::rbindlist(gene_values),
                  FigS9BCD_TU_values = data.table::rbindlist(TU_values), FigS9BCD_labels = data.table::rbindlist(labels_used),
                  FigS9BCD_statistics = data.table::rbindlist(all_tests), Analysis_settings = settings), OUTPUT_DIR)
for (i in seq_along(plots)) save_figure(plots[[i]], stems[i], WIDTH_IN, HEIGHT_IN)
FigS9BCD_plots <- plots
message("Saved Fig. S9B-D to ", OUTPUT_DIR)
