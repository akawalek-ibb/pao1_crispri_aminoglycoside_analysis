#!/usr/bin/env Rscript
# Figure S9A: severity versus EK50 across untreated T1 control experiments.
# Run from the project root. First complete preprocessing steps 11 and 12.
# Reads the annotated consensus used for Table S5 and Fig. 2G; no guide-response
# curve is refitted. EK50 remains missing when fewer than two controls yielded it.
# The only fitted curve here describes the association between severity and EK50.
# Flow: read consensus -> add TU counts -> select finite estimates -> choose labels
# -> calculate Spearman/GAM summaries -> draw the figure -> save plotted data.
#
# Input: preprocessing_scripts/out_12_AMG_MH_fitness/Fig2G_input.tsv;
# annotation_input/PAO1_curated.tsv supplies TU membership/counts.
# Run: source("Figure_scripts/FigS9/AMG_FigS9A_severity_EK50.R")
# Outputs: PNG/PDF/SVG, plotted genes, fitted curve, statistics and settings TSVs
# in Figure_scripts/FigS9/out_AMG_FigS9A_severity_EK50/.
#
# R packages used/checked: data.table, ggnewscale, ggplot2, ggrepel, mgcv, scales,
# svglite.

# 1. Inputs, display settings and output location.
INPUT_FILE <- "preprocessing_scripts/out_12_AMG_MH_fitness/Fig2G_input.tsv"
ANNOTATION_FILE <- "annotation_input/PAO1_curated.tsv"
OUTPUT_DIR <- "Figure_scripts/FigS9/out_AMG_FigS9A_severity_EK50"

# Count both essential/fitness-associated design groups by default.
# Use "Nr_E" for previously classified essential genes only; that changes point sizes.
TU_COUNT_DEFINITION <- "design_call"
MIN_SEVERITY_ESTIMATES <- 1L
MIN_EK50_ESTIMATES <- 2L
SIZE_POINTS_BY_TU <- TRUE
SEVERITY_LABEL_CUTOFF <- -12
EK50_LABEL_CUTOFF <- 0.40
MAX_AUTOMATIC_LABELS <- 50L
MANUAL_LABELS <- c("parE", "gyrA", "gyrB", "parC", "rpoB", "murA", "ftsI", "pbpA", "pbpG",
                   "folA", "folP", "fabI", "rpsL", "rplD", "rplC", "lpxC", "bamA", "lptD", "fabF", "mreB")
CATEGORY_COLOURS <- c(`Core essential` = "darkred", `Conditional essential` = "#0072B2",
                      `Essential others` = "#FC8D59", Rest = "darkgrey")
FONT_FAMILY <- "Arial"
WIDTH_IN <- 8.5
HEIGHT_IN <- 6.5
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

require_packages(c("data.table", "ggplot2", "ggrepel", "ggnewscale", "mgcv", "scales", "svglite"))

# 2. Reuse the gene-wise medians across B1-B4 untreated T1 versus T0.
# The prepared file also contains the frozen essentiality categories used in Fig. 2G.
x <- read_table(INPUT_FILE, c("PAO1_ID", "gene_name", "category_label", "control_median_severity",
                             "control_median_EK50", "n_control_severity_estimates", "n_control_EK50_estimates"))
x <- attach_TU_counts(x, ANNOTATION_FILE, TU_COUNT_DEFINITION)
x[, `:=`(severity = as.numeric(control_median_severity), EK50 = as.numeric(control_median_EK50))]
x <- x[n_control_severity_estimates >= MIN_SEVERITY_ESTIMATES &
       n_control_EK50_estimates >= MIN_EK50_ESTIMATES & is.finite(severity) & is.finite(EK50)]
if (nrow(x) < 10L || data.table::uniqueN(x$severity) < 6L) stop("Insufficient data for the original GAM/correlation.")
if (any(!x$category_label %in% names(CATEGORY_COLOURS))) stop("Unrecognized essentiality category.")
data.table::setorder(x, PAO1_ID)
x[, category_label := factor(category_label, levels = names(CATEGORY_COLOURS))]
x[, label_txt := ifelse(is.na(gene_name) | !nzchar(trimws(gene_name)), PAO1_ID, gsub("'", "", gene_name, fixed = TRUE))]
x[, `:=`(is_known_target = label_txt %in% MANUAL_LABELS,
          size_TU = log1p(ifelse(is.finite(essential_genes_in_operon), essential_genes_in_operon, 0)))]

# 3. Label the specified antimicrobial targets and the most extreme remaining genes.
# These label cutoffs are display choices, not new significance tests.
automatic <- data.table::copy(x[!is_known_target & (severity <= SEVERITY_LABEL_CUTOFF | EK50 <= EK50_LABEL_CUTOFF)])
automatic[, label_priority := pmin(data.table::frank(severity, ties.method = "min"), data.table::frank(EK50, ties.method = "min"))]
data.table::setorder(automatic, label_priority, EK50, severity)
labels <- unique(data.table::rbindlist(list(x[is_known_target == TRUE], head(automatic, MAX_AUTOMATIC_LABELS)), fill = TRUE), by = "PAO1_ID")
labels[, label_colour := ifelse(is_known_target, "darkviolet", unname(CATEGORY_COLOURS[as.character(category_label)]))]

# 4. Describe the relationship among the plotted genes. The GAM uses the original
# six-basis smooth and REML fitting. Its smooth p-value is not a test of nonlinearity.
set.seed(SEED)
correlation <- stats::cor.test(x$severity, x$EK50, method = "spearman", exact = FALSE)
model <- mgcv::gam(EK50 ~ s(severity, k = 6), data = as.data.frame(x), method = "REML")
sm <- summary(model)
p_label <- function(p) if (p < 2.2e-16) "< 2.2e-16" else paste0("= ", format(signif(p, 3), scientific = TRUE))
annotation <- paste0("Spearman rho = ", format(round(unname(correlation$estimate), 2), nsmall = 2),
                     ", p ", p_label(correlation$p.value), "\nGAM edf = ",
                     format(round(sm$s.table[1, "edf"], 2), nsmall = 2), ", p ", p_label(sm$s.table[1, "p-value"]))

# 5. Draw the same points, labels, smooth and confidence band as the source analysis.
# The TU-size legend names the actual counted groups rather than calling all of them essential.
size_title <- if (TU_COUNT_DEFINITION == "design_call") "Essential/fitness-associated\ngenes in TU" else "Previously classified\nessential genes in TU"
p <- ggplot2::ggplot(x, ggplot2::aes(x = severity, y = EK50)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey85") +
  ggplot2::geom_hline(yintercept = EK50_LABEL_CUTOFF, linetype = "dotted", linewidth = 0.8, colour = "grey35")
if (SIZE_POINTS_BY_TU) {
  p <- p + ggplot2::geom_point(ggplot2::aes(colour = category_label, size = size_TU), alpha = 0.85) +
    ggplot2::scale_size_continuous(range = c(0.25, 1.5), name = size_title,
                                    breaks = log1p(c(0, 1, 2, 4, 6)), labels = c("0", "1", "2", "4", "6"))
} else p <- p + ggplot2::geom_point(ggplot2::aes(colour = category_label), size = 1.1, alpha = 0.85)
p <- p +
  ggplot2::scale_colour_manual(values = CATEGORY_COLOURS, breaks = names(CATEGORY_COLOURS), drop = FALSE, name = "Essentiality annotation") +
  ggnewscale::new_scale_colour() +
  ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = label_txt, colour = label_colour),
                         size = 2, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.2,
                         min.segment.length = 0, segment.color = "grey60", segment.alpha = 0.4,
                         show.legend = FALSE, fontface = "italic", seed = SEED) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(6), limits = c(0.05, 1)) +
  ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, k = 6), method.args = list(method = "REML"),
                       se = TRUE, colour = "black", fill = "grey80", linewidth = 0.8) +
  ggplot2::annotate("text", x = -Inf, y = Inf, label = annotation, hjust = -0.05, vjust = 1.15, size = 2.4, family = FONT_FAMILY) +
  ggplot2::theme_minimal(base_size = 8, base_family = FONT_FAMILY) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_line(linewidth = 0.25),
                 legend.position = "right", legend.box = "vertical", legend.title = ggplot2::element_text(face = "bold"),
                 legend.text = ggplot2::element_text(size = 7.5), axis.title = ggplot2::element_text(face = "bold", size = 8)) +
  ggplot2::labs(x = "Severity (median full-match-guide adjusted log2 fold change)",
                 y = "EK50 (predicted guide efficacy at half the fitted effect at efficacy 1)",
                 subtitle = "Untreated T1 consensus across B1-B4")

# 6. Save figures and the exact plotted data. Export the smooth from ggplot's
# own calculation so its saved curve matches the line and band drawn in the panel.
outputs <- c(paste0("FigS9A_severity_EK50", c(".svg", ".pdf", ".png")),
             paste0(c("FigS9A_plotted_genes", "FigS9A_labels", "FigS9A_statistics", "FigS9A_fitted_curve", "Analysis_settings"), ".tsv"))
if (!OVERWRITE && any(file.exists(file.path(OUTPUT_DIR, outputs)))) stop("Output exists; set OVERWRITE <- TRUE.")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
built <- ggplot2::ggplot_build(p)
smooth_layer <- which(vapply(p$layers, function(layer) inherits(layer$stat, "StatSmooth"), logical(1)))
if (length(smooth_layer) != 1L) stop("Cannot identify the fitted curve layer.")
curve <- built$data[[smooth_layer]][, intersect(c("x", "y", "ymin", "ymax", "se"), names(built$data[[smooth_layer]])), drop = FALSE]
statistics <- data.table::data.table(metric = c("n_genes", "Spearman_rho", "Spearman_p", "GAM_edf", "GAM_p"),
                                    value = c(nrow(x), unname(correlation$estimate), correlation$p.value,
                                              sm$s.table[1, "edf"], sm$s.table[1, "p-value"]))
settings <- run_settings(c(TU_COUNT_DEFINITION = TU_COUNT_DEFINITION, MIN_SEVERITY_ESTIMATES = MIN_SEVERITY_ESTIMATES,
                           MIN_EK50_ESTIMATES = MIN_EK50_ESTIMATES, consensus = "B1-B4 untreated T1 only", SEED = SEED),
                         c(consensus = INPUT_FILE, annotation = ANNOTATION_FILE))
write_tables(list(FigS9A_plotted_genes = x, FigS9A_labels = labels, FigS9A_statistics = statistics,
                  FigS9A_fitted_curve = curve, Analysis_settings = settings), OUTPUT_DIR)
save_figure(p, "FigS9A_severity_EK50", WIDTH_IN, HEIGHT_IN)
FigS9A_plot <- p
message("Saved Fig. S9A: ", nrow(x), " genes; ", OUTPUT_DIR)
