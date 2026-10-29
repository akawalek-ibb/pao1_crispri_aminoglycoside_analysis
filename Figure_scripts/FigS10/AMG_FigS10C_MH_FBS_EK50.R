#!/usr/bin/env Rscript
# Compare MH and FBS EK50 in Supplementary Figure S10C.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/FigS10/AMG_FigS10C_MH_FBS_EK50.R")
# Inputs:
# - preprocessing_scripts/out_13_AMG_MH_FBS_gene_effects/MH_FBS_gene_effects.tsv and
#   annotation_input/PAO1_curated.tsv; match PAO1_ID.
# - Gene effects and vulnerability/status columns are distinct inputs within the same prepared
#   gene table.
# Upstream / execution order:
# - preprocessing_scripts/11_AMG_gene_vulnerability.R →
#   preprocessing_scripts/13_AMG_MH_FBS_gene_effects.R. Use supplied edgeR_data fits;
#   preprocessing_scripts/02_AMG_metacounter.R is needed only for a refit.
# Outputs:
# - Figure_scripts/FigS10/out_AMG_FigS10C_EK50/ — PNG/PDF/SVG and plotted
#   values/statistics/settings TSVs.
# Interpretation / operation:
# - Require negative fitness and valid numeric EK50 in both media; no boundary substitution.
# - No guide-response curve is refitted here and no per-generation scaling is applied.
#
# R packages used/checked: ggplot2, ggrepel, svglite.

INPUT_FILE <- "preprocessing_scripts/out_13_AMG_MH_FBS_gene_effects/MH_FBS_gene_effects.tsv"
ANNOTATION_FILE <- "annotation_input/PAO1_curated.tsv"
OUTPUT_DIR <- "Figure_scripts/FigS10/out_AMG_FigS10C_EK50"
PRINT_IN_RSTUDIO <- TRUE
OVERWRITE <- TRUE
DPI <- 300

P_STOUFFER_CUT <- 0.05
NEGATIVE_LFC_CUT <- -1
MAX_LABELS_ESTIMATED_BOTH <- 35L
FONT_FAMILY <- "Arial"
WIDTH_IN <- 7.2
HEIGHT_IN <- 6.0
CATEGORY_COLOURS <- c("Core essential" = "#8a0000ff", "Conditional essential" = "#0070b0ff",
  "Essential others" = "#ff6600ff", "Rest" = "#a8a8a8ff")

# 2. Read the prepared gene table and its existing annotations ---------------
packages <- c("ggplot2", "ggrepel", "svglite")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required packages: ", paste(missing_packages, collapse = ", "))
if (!file.exists(INPUT_FILE)) stop("Missing prepared gene table: ", INPUT_FILE,
  ". Use the existing preprocessing 13 output, or generate it after preprocessing 11.")
if (!file.exists(ANNOTATION_FILE)) stop("Missing annotation: ", ANNOTATION_FILE)

genes <- utils::read.delim(INPUT_FILE, check.names = FALSE, stringsAsFactors = FALSE,
  quote = "\"", comment.char = "", na.strings = c("", "NA"))
annotation <- utils::read.delim(ANNOTATION_FILE, check.names = FALSE, stringsAsFactors = FALSE,
  quote = "\"", comment.char = "", na.strings = c("", "NA"))
required <- c("PAO1_ID", "gene_name", "MH_medLFC_vs_T0", "MH_p_stouffer_vs_T0", "FBS_medLFC_vs_T0", "FBS_p_stouffer_vs_T0", "MH_EK50", "FBS_EK50", "MH_EK50_status", "FBS_EK50_status")
if (!all(required %in% names(genes))) stop("Prepared gene table lacks: ",
  paste(setdiff(required, names(genes)), collapse = ", "))
if (!all(c("PAO1_ID", "category", "design_call") %in% names(annotation))) {
  stop("Annotation must contain PAO1_ID, category and design_call.")
}
if (anyNA(genes$PAO1_ID) || anyDuplicated(genes$PAO1_ID) || any(!nzchar(genes$PAO1_ID))) {
  stop("The prepared gene table must contain one nonempty PAO1_ID per gene.")
}
if (anyNA(annotation$PAO1_ID) || anyDuplicated(annotation$PAO1_ID)) {
  stop("The annotation must contain unique, nonmissing PAO1_ID values.")
}

# The curated design annotation already encodes Poulsen core/conditional
# assignments. Genes outside that annotation are retained in the Rest category.
# Essentiality categories use the priority and palette configured for this panel.
category_text <- tolower(ifelse(is.na(annotation$category), "", annotation$category))
annotation$essential_category <- "Rest"
annotation$essential_category[annotation$design_call %in% c("main_essential", "other_essential") |
    annotation$category %in% c("E", "E_AMR", "AMR") |
    grepl("^e_", category_text)] <- "Essential others"
annotation$essential_category[grepl("conditional", category_text)] <- "Conditional essential"
annotation$essential_category[grepl("core", category_text) & !grepl("conditional", category_text)] <- "Core essential"
genes$essential_category <- annotation$essential_category[match(genes$PAO1_ID, annotation$PAO1_ID)]
genes$essential_category[is.na(genes$essential_category)] <- "Rest"
genes$essential_category <- factor(genes$essential_category, levels = names(CATEGORY_COLOURS))
genes$label_txt <- ifelse(is.na(genes$gene_name) | !nzchar(trimws(genes$gene_name)),
  genes$PAO1_ID, gsub("'", "", genes$gene_name, fixed = TRUE))

# 3. Identify negative fitness effects against the shared T0 -----------------
# These are raw Stouffer p values, not BH-adjusted values.
genes$MH_negative_fitness <- is.finite(genes$MH_p_stouffer_vs_T0) &
  genes$MH_p_stouffer_vs_T0 < P_STOUFFER_CUT & is.finite(genes$MH_medLFC_vs_T0) &
  genes$MH_medLFC_vs_T0 <= NEGATIVE_LFC_CUT
genes$FBS_negative_fitness <- is.finite(genes$FBS_p_stouffer_vs_T0) &
  genes$FBS_p_stouffer_vs_T0 < P_STOUFFER_CUT & is.finite(genes$FBS_medLFC_vs_T0) &
  genes$FBS_medLFC_vs_T0 <= NEGATIVE_LFC_CUT
# Panel C uses the INTERSECTION: a negative effect in BOTH MH and FBS.
# The old standalone script defaulted to "either" despite its header saying
# "both". Here the explicit intersection follows the current panel C legend.
genes$retained_by_fitness_filter <- genes$MH_negative_fitness & genes$FBS_negative_fitness

# A valid estimate must have status "estimated" and a finite numeric EK50.
# Failed fits, inadequate support and no half-effect crossing remain excluded.
# Unlike sensitivity panel B, EK50 has no artificial boundary points here.
for (medium in c("MH", "FBS")) {
  status <- genes[[paste0(medium, "_EK50_status")]]
  value <- genes[[paste0(medium, "_EK50")]]
  if (any(is.finite(value) & !(status %in% "estimated")) ||
      any(status %in% "estimated" & !is.finite(value))) {
    stop("Inconsistent numeric EK50/status in ", medium, "; refresh preprocessing 11/13.")
  }
  if (any(is.finite(value) & (value < 0.1 | value > 1))) stop("EK50 outside the expected fitting range.")
}
genes$valid_EK50_both <- genes$MH_EK50_status %in% "estimated" & genes$FBS_EK50_status %in% "estimated" &
  is.finite(genes$MH_EK50) & is.finite(genes$FBS_EK50)
genes$plotted <- genes$retained_by_fitness_filter & genes$valid_EK50_both
plot_genes <- genes[genes$plotted, ]
statistics_genes <- plot_genes

# 4. Statistics use numeric estimates in both conditions ---------------------
# Pair each gene with itself, using the same rows for correlations and the
# two-sided paired Wilcoxon signed-rank test. The tested difference is FBS-MH.
# No artificial boundary value enters any of these tests.
if (nrow(statistics_genes) < 3L) stop("Fewer than three genes with paired numeric EK50 estimates.")
cor_s <- suppressWarnings(stats::cor.test(statistics_genes$MH_EK50, statistics_genes$FBS_EK50,
  method = "spearman", exact = FALSE))
cor_p <- suppressWarnings(stats::cor.test(statistics_genes$MH_EK50, statistics_genes$FBS_EK50,
  method = "pearson"))
paired_test <- suppressWarnings(stats::wilcox.test(statistics_genes$FBS_EK50, statistics_genes$MH_EK50,
  paired = TRUE, exact = FALSE))
statistics <- data.frame(metric = c("n_MH_negative_fitness", "n_FBS_negative_fitness",
  "n_retained_by_fitness_filter", "n_plotted", "n_numeric_estimates_both", "n_excluded_after_fitness_filter",
  "spearman_rho", "spearman_p", "pearson_R", "pearson_p", "median_MH", "median_FBS",
  "paired_Wilcoxon_W", "paired_Wilcoxon_p"),
  value = c(sum(genes$MH_negative_fitness), sum(genes$FBS_negative_fitness), sum(genes$retained_by_fitness_filter),
    nrow(plot_genes), nrow(statistics_genes), sum(genes$retained_by_fitness_filter & !genes$plotted),
    unname(cor_s$estimate), cor_s$p.value, unname(cor_p$estimate), cor_p$p.value,
    stats::median(statistics_genes$MH_EK50), stats::median(statistics_genes$FBS_EK50),
    unname(paired_test$statistic), paired_test$p.value))
subtitle_text <- paste0("Estimated in both: n = ", nrow(statistics_genes),
  "; Spearman rho = ", format(round(unname(cor_s$estimate), 2), nsmall = 2),
  "; p ", format.pval(cor_s$p.value, digits = 3, eps = 2.2e-16))

# 5. Plot paired EK50 estimates ----------------------------------------------
# Use the configured palette, font, dimensions and largest-difference labels.
# The displayed n is calculated from the input; it never selects the gene set.
label_genes <- head(plot_genes[order(-abs(plot_genes$FBS_EK50 - plot_genes$MH_EK50)), ], MAX_LABELS_ESTIMATED_BOTH)
plot_genes$selected_for_label <- plot_genes$PAO1_ID %in% label_genes$PAO1_ID
p <- ggplot2::ggplot(plot_genes, ggplot2::aes(MH_EK50, FBS_EK50)) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linewidth = 0.5, colour = "grey45") +
  ggplot2::geom_point(ggplot2::aes(colour = essential_category), size = 1.7, alpha = 0.78, stroke = 0.45) +
  ggplot2::scale_colour_manual(values = CATEGORY_COLOURS, drop = FALSE, name = "Essential category") +
  ggrepel::geom_text_repel(data = label_genes, ggplot2::aes(label = label_txt, colour = essential_category),
    size = 2.4, max.overlaps = Inf, box.padding = 0.35, point.padding = 0.20,
    min.segment.length = 0, segment.color = "grey60", fontface = "italic", show.legend = FALSE, seed = 1) +
  ggplot2::scale_x_continuous(limits = c(0.05, 1.04), breaks = seq(0.1, 1, by = 0.2)) +
  ggplot2::scale_y_continuous(limits = c(0.05, 1.04), breaks = seq(0.1, 1, by = 0.2)) +
  ggplot2::coord_equal(clip = "off") +
  ggplot2::theme_classic(base_family = FONT_FAMILY, base_size = 9) +
  ggplot2::theme(legend.position = "right", axis.title = ggplot2::element_text(face = "bold"),
    plot.margin = ggplot2::margin(5.5, 45, 20, 5.5)) +
  ggplot2::labs(x = "EK50 in MH", y = "EK50 in FBS", subtitle = subtitle_text)

# 6. Export the figure and the values needed to reproduce it ----------------
# The all-gene table records exclusions. The plotted table records displayed
# coordinates and labels; statistical estimates never use boundary substitutes.
figure_prefix <- "FigS10C_MH_FBS_EK50"
tables <- list(all_genes = genes, plotted_genes = plot_genes,
  statistics = statistics,
  Analysis_settings = data.frame(setting = c("INPUT_FILE", "input_md5", "ANNOTATION_FILE",
    "annotation_md5", "R_version", "ggplot2_version", "ggrepel_version", "P_STOUFFER_CUT", "NEGATIVE_LFC_CUT", "fitness_filter", "EK50_filter", "statistics_subset"),
    value = c(INPUT_FILE, unname(tools::md5sum(INPUT_FILE)), ANNOTATION_FILE,
      unname(tools::md5sum(ANNOTATION_FILE)), R.version.string,
      as.character(utils::packageVersion("ggplot2")), as.character(utils::packageVersion("ggrepel")),
      P_STOUFFER_CUT, NEGATIVE_LFC_CUT, "negative in both MH and FBS", "finite estimated values in both; no boundary substitution", "all plotted genes"), stringsAsFactors = FALSE))
table_paths <- file.path(OUTPUT_DIR, paste0(figure_prefix, "_", names(tables), ".tsv"))
figure_paths <- file.path(OUTPUT_DIR, paste0(figure_prefix, c(".png", ".pdf", ".svg")))
if (!OVERWRITE && any(file.exists(c(table_paths, figure_paths)))) {
  stop("Output files already exist; set OVERWRITE <- TRUE to replace them.")
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (i in seq_along(tables)) {
  utils::write.table(tables[[i]], table_paths[i], sep = "\t", quote = TRUE,
    row.names = FALSE, na = "NA", qmethod = "double")
}

# Keep the named plot object available for redisplay in RStudio.
FigS10C_plot <- p
if (isTRUE(PRINT_IN_RSTUDIO)) print(FigS10C_plot)
ggplot2::ggsave(figure_paths[1], FigS10C_plot, width = WIDTH_IN, height = HEIGHT_IN,
  units = "in", dpi = DPI, bg = "white")
ggplot2::ggsave(figure_paths[2], FigS10C_plot, width = WIDTH_IN, height = HEIGHT_IN,
  units = "in", device = grDevices::cairo_pdf, bg = "white")
ggplot2::ggsave(figure_paths[3], FigS10C_plot, width = WIDTH_IN, height = HEIGHT_IN,
  units = "in", device = svglite::svglite, bg = "white")
message("Saved Fig. S10C:\n", paste(figure_paths, collapse = "\n"))
message("Underlying figure data and statistics:\n", paste(table_paths, collapse = "\n"))
message("Panel C gene count for the legend: n = ", nrow(plot_genes),
  ". Use this input-derived count when updating the manuscript caption.")
