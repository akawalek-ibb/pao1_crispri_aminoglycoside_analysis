#!/usr/bin/env Rscript
# Compare MH and FBS severity in Supplementary Figure S10A.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/FigS10/AMG_FigS10A_MH_FBS_severity.R")
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
# - Figure_scripts/FigS10/out_AMG_FigS10A_severity/ — PNG/PDF/SVG and plotted
#   values/statistics/settings TSVs.
# Interpretation / operation:
# - Plot genes with finite paired effects; labels use the direct FBS-versus-MH raw-p/effect
#   criteria.
# - No guide-response curve is refitted here and no per-generation scaling is applied.
#
# R packages used/checked: ggplot2, ggrepel, svglite.

INPUT_FILE <- "preprocessing_scripts/out_13_AMG_MH_FBS_gene_effects/MH_FBS_gene_effects.tsv"
ANNOTATION_FILE <- "annotation_input/PAO1_curated.tsv"
OUTPUT_DIR <- "Figure_scripts/FigS10/out_AMG_FigS10A_severity"
PRINT_IN_RSTUDIO <- TRUE
OVERWRITE <- TRUE
DPI <- 300

P_STOUFFER_CUT <- 0.05
DIRECT_LFC_CUT <- 1
LABEL_SIGNIFICANT_GENES <- TRUE
MAX_LABEL_OVERLAPS <- 30
LABEL_SIZE <- 2.2
P_STOUFFER_FLOOR <- 1e-50
MAX_NEG_LOG10_P <- 50
POINT_SIZE_RANGE <- c(0.5, 3)
FONT_FAMILY <- "sans"
WIDTH_IN <- 10
HEIGHT_IN <- 8
CATEGORY_COLOURS <- c("Core essential" = "darkred", "Conditional essential" = "#0072B2",
  "Essential others" = "#FC8D59", "Rest" = "darkgrey")

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
required <- c("PAO1_ID", "gene_name", "MH_medLFC_vs_T0", "MH_p_stouffer_vs_T0", "FBS_medLFC_vs_T0", "FBS_p_stouffer_vs_T0", "direct_FBS_vs_MH_medLFC", "direct_FBS_vs_MH_p_stouffer")
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
annotation$essential_category[annotation$category %in% c("E", "E_AMR")] <- "Essential others"
annotation$essential_category[grepl("conditional", category_text)] <- "Conditional essential"
annotation$essential_category[grepl("core", category_text) & !grepl("conditional", category_text)] <- "Core essential"
genes$essential_category <- annotation$essential_category[match(genes$PAO1_ID, annotation$PAO1_ID)]
genes$essential_category[is.na(genes$essential_category)] <- "Rest"
genes$essential_category <- factor(genes$essential_category, levels = names(CATEGORY_COLOURS))
genes$label_txt <- ifelse(is.na(genes$gene_name) | !nzchar(trimws(genes$gene_name)),
  genes$PAO1_ID, gsub("'", "", genes$gene_name, fixed = TRUE))

# 3. Select genes with measured effects in all three matched B3 contrasts ----
# Axes are MH T1/T0 and FBS T1/T0; both share the same pre-induced MH T0.
# Labels and point sizes use the separately tested FBS T1 versus MH T1 contrast.
# We do not infer a direct-comparison p value by subtracting the two axes.
genes$plotted <- is.finite(genes$MH_medLFC_vs_T0) & is.finite(genes$FBS_medLFC_vs_T0) &
  is.finite(genes$direct_FBS_vs_MH_medLFC) & is.finite(genes$direct_FBS_vs_MH_p_stouffer)
genes$direct_significant <- is.finite(genes$direct_FBS_vs_MH_p_stouffer) &
  genes$direct_FBS_vs_MH_p_stouffer < P_STOUFFER_CUT &
  is.finite(genes$direct_FBS_vs_MH_medLFC) & abs(genes$direct_FBS_vs_MH_medLFC) >= DIRECT_LFC_CUT
plot_genes <- genes[genes$plotted, ]
if (nrow(plot_genes) < 3L) stop("Fewer than three complete matched B3 gene records.")
plot_genes$label <- ifelse(LABEL_SIGNIFICANT_GENES & plot_genes$direct_significant, plot_genes$label_txt, NA_character_)
plot_genes$alpha_sig <- ifelse(plot_genes$direct_significant, 0.90, 0.25)
plot_genes$neg_log10_p_stouffer <- pmin(-log10(pmax(pmin(plot_genes$direct_FBS_vs_MH_p_stouffer, 1),
  P_STOUFFER_FLOOR)), MAX_NEG_LOG10_P)
plot_genes$observed_shift_from_MH <- plot_genes$FBS_medLFC_vs_T0 - plot_genes$MH_medLFC_vs_T0
plot_genes$direct_minus_derived <- plot_genes$direct_FBS_vs_MH_medLFC - plot_genes$observed_shift_from_MH

# 4. Summarise agreement across environments --------------------------------
cor_s <- stats::cor.test(plot_genes$MH_medLFC_vs_T0, plot_genes$FBS_medLFC_vs_T0,
  method = "spearman", exact = FALSE)
cor_p <- stats::cor.test(plot_genes$MH_medLFC_vs_T0, plot_genes$FBS_medLFC_vs_T0, method = "pearson")
statistics <- data.frame(metric = c("n_genes_plotted", "n_direct_significant",
  "n_stronger_defect_in_FBS", "n_weaker_defect_in_FBS", "spearman_rho", "spearman_p",
  "pearson_R", "pearson_p", "median_absolute_direct_minus_derived"),
  value = c(nrow(plot_genes), sum(plot_genes$direct_significant),
    sum(plot_genes$direct_significant & plot_genes$direct_FBS_vs_MH_medLFC < 0),
    sum(plot_genes$direct_significant & plot_genes$direct_FBS_vs_MH_medLFC > 0),
    unname(cor_s$estimate), cor_s$p.value, unname(cor_p$estimate), cor_p$p.value,
    stats::median(abs(plot_genes$direct_minus_derived))))

# 5. Draw severity comparison ------------------------------------------------
# Original palette, point-size/alpha rules, LOESS and font/dimensions are kept.
# The solid equality line and dashed +/-1 lines explicitly match the caption.
# LOESS is descriptive; it is not used for gene significance calls.
p <- ggplot2::ggplot(plot_genes, ggplot2::aes(MH_medLFC_vs_T0, FBS_medLFC_vs_T0,
    size = neg_log10_p_stouffer)) +
  ggplot2::geom_point(ggplot2::aes(colour = essential_category, alpha = alpha_sig)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  ggplot2::geom_smooth(method = "loess", se = TRUE, colour = "black", linewidth = 0.9) +
  ggplot2::geom_abline(intercept = 0, slope = 1, colour = "grey45", linewidth = 0.5) +
  ggplot2::geom_abline(intercept = c(-1, 1), slope = 1, linetype = "dashed", colour = "blue") +
  ggrepel::geom_text_repel(ggplot2::aes(label = label, colour = essential_category),
    fontface = "italic", max.overlaps = MAX_LABEL_OVERLAPS, size = LABEL_SIZE,
    show.legend = FALSE, min.segment.length = 0, segment.color = "black", seed = 1, na.rm = TRUE) +
  ggplot2::scale_colour_manual(values = CATEGORY_COLOURS, name = "Essential category") +
  ggplot2::scale_size_continuous(name = expression(-log[10]*"(p_stouffer)"), range = POINT_SIZE_RANGE) +
  ggplot2::scale_alpha_identity(guide = "none") +
  ggplot2::scale_x_continuous(breaks = c(-15, -10, -5, 0)) +
  ggplot2::scale_y_continuous(breaks = c(-15, -10, -5, 0)) +
  ggplot2::labs(x = "MH severity: T1 vs T0 (median log2FC)", y = "FBS severity: T1 vs T0 (median log2FC)") +
  ggplot2::theme_minimal(base_size = 14, base_family = FONT_FAMILY) +
  ggplot2::theme(legend.position = "right", legend.key = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "grey85"), panel.grid.minor = ggplot2::element_blank())

# 6. Export the figure and the values needed to reproduce it ----------------
# The all-gene table records exclusions. The plotted table records displayed
# coordinates and labels; statistical estimates never use boundary substitutes.
figure_prefix <- "FigS10A_MH_FBS_severity"
tables <- list(all_genes = genes, plotted_genes = plot_genes,
  statistics = statistics,
  Analysis_settings = data.frame(setting = c("INPUT_FILE", "input_md5", "ANNOTATION_FILE",
    "annotation_md5", "R_version", "ggplot2_version", "ggrepel_version", "P_STOUFFER_CUT", "DIRECT_LFC_CUT", "comparison"),
    value = c(INPUT_FILE, unname(tools::md5sum(INPUT_FILE)), ANNOTATION_FILE,
      unname(tools::md5sum(ANNOTATION_FILE)), R.version.string,
      as.character(utils::packageVersion("ggplot2")), as.character(utils::packageVersion("ggrepel")),
      P_STOUFFER_CUT, DIRECT_LFC_CUT, "B3 matched MH T1/T0; FBS T1/T0; FBS T1/MH T1"), stringsAsFactors = FALSE))
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
FigS10A_plot <- p
if (isTRUE(PRINT_IN_RSTUDIO)) print(FigS10A_plot)
ggplot2::ggsave(figure_paths[1], FigS10A_plot, width = WIDTH_IN, height = HEIGHT_IN,
  units = "in", dpi = DPI, bg = "white")
ggplot2::ggsave(figure_paths[2], FigS10A_plot, width = WIDTH_IN, height = HEIGHT_IN,
  units = "in", device = grDevices::cairo_pdf, bg = "white")
ggplot2::ggsave(figure_paths[3], FigS10A_plot, width = WIDTH_IN, height = HEIGHT_IN,
  units = "in", device = svglite::svglite, bg = "white")
message("Saved Fig. S10A:\n", paste(figure_paths, collapse = "\n"))
message("Underlying figure data and statistics:\n", paste(table_paths, collapse = "\n"))
message("Plotted genes: ", nrow(plot_genes), "; direct condition-dependent hits: ",
  sum(plot_genes$direct_significant), "; Pearson R = ", round(unname(cor_p$estimate), 3))
