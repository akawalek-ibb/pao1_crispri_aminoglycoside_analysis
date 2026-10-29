# Plot library representation and guide/gene coverage for Figure S6.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/FigS6/AMG_FigS6_library_representation.R")
# Inputs:
# - preprocessing_scripts/out_10_AMG_library_representation/Figure_S6A.tsv through
#   Figure_S6D.tsv, Sample_key.tsv and Design_summary.tsv.
# Upstream / execution order:
# - preprocessing_scripts/09_AMG_export_GEO_raw_counts.R →
#   preprocessing_scripts/10_AMG_library_representation.R.
# Outputs:
# - Figure_scripts/FigS6/out_AMG_FigS6_library_representation/ — Separate SVG/PDF panels and
#   PNGs (enabled by default).
# Interpretation / operation:
# - Sample names/order come from the prepared sample key; coverage summaries are not recalculated.
#
# R packages used/checked: ggplot2, scales, svglite.

INPUT_DIR <- "preprocessing_scripts/out_10_AMG_library_representation"

OUTPUT_DIR <- "Figure_scripts/FigS6/out_AMG_FigS6_library_representation"

FIGURE_ID <- "S6"

FONT_FAMILY <- "sans"

ABC_HEIGHT <- 11

PANEL_WIDTHS <- c(A = 8.5, B = 13, C = 8.5, D = 11)

D_HEIGHT <- 4

EXPORT_PNG <- TRUE

PNG_DPI <- 300

PRINT_IN_RSTUDIO <- TRUE

SHOW_SEGMENT_COUNTS <- TRUE

MIN_SEGMENT_LABEL_FRACTION <- 0.04

BATCH_COLORS <- c(B1 = "#1b9e77", B2 = "#d95f02", B3 = "#7570b3", B4 = "#e7298a")

READ_BIN_COLORS <- c(`0 reads` = "grey88", `1-10` = "#66c2a5", `11-100` = "#fc8d62", `>100` = "#8da0cb")

GENE_COVERAGE_COLORS <- c(`0` = "grey88", `1` = "#66c2a5", `2` = "#fc8d62", `3-10` = "#8da0cb", 
  `>10` = "#e78ac3")

DESIGN_CALL_COLORS <- c(main_essential = "#7570b3", other_essential = "#d95f02", non_essential = "#1b9e77")

DESIGN_CALL_LABELS <- c(main_essential = "Primary essential", other_essential = "Fitness-associated", 
  non_essential = "Other genes")

# Local functions for this analysis.
# Read, prepare, plot and export in order.
packages <- c("ggplot2", "scales", "svglite")

missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) stop("Required packages: ", paste(missing, collapse = ", "))

read_input <- function(name) {
  path <- file.path(INPUT_DIR, paste0(name, ".tsv"))
  if (!file.exists(path)) 
    stop("Run 10_AMG_library_representation.R first. Missing: ", path)
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "", 
    na.strings = "NA")
}

# 1. Reuse task 10's exact sample order and four panel summaries.
# This script does not read raw count tables or apply new detection criteria.
samples <- rev(read_input("Sample_key")$label)

a <- read_input("Figure_S6A")

b <- read_input("Figure_S6B")

c <- read_input("Figure_S6C")

d <- read_input("Figure_S6D")

a$batch <- factor(a$batch, levels = names(BATCH_COLORS))

# 2. Translate guide classes and preserve the configured display order.
# The underlying proportions/counts are imported unchanged.
b$sgRNA_class <- factor(b$sgRNA_class, levels = c("perfect", "mismatch", "control"), labels = c("Full-match", 
  "Mismatch", "Non-targeting"))

b$count_bin <- factor(b$count_bin, levels = names(READ_BIN_COLORS))

c$guides_bin <- factor(c$guides_bin, levels = names(GENE_COVERAGE_COLORS))

d$design_call <- factor(d$design_call, levels = read_input("Design_summary")$design_call)

common <- ggplot2::theme_bw(base_size = 10, base_family = FONT_FAMILY) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), 
  panel.grid.minor = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 7), 
  strip.background = ggplot2::element_rect(fill = "grey95"), strip.text = ggplot2::element_text(face = "bold"), 
  legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold"))

sample_axis <- ggplot2::scale_y_discrete(limits = samples, drop = FALSE)

# 3. Panel A shows total assigned guide counts, not total sequenced read pairs.
# The division by one million only changes the axis units.
p_A <- ggplot2::ggplot(a, ggplot2::aes(x = total_reads/1e+06, y = label, fill = batch)) +
  ggplot2::geom_col(width = 0.8) + 
  ggplot2::geom_text(ggplot2::aes(label = scales::comma(total_reads)), hjust = 1.03, size = 2.3) + 
  sample_axis +
  ggplot2::scale_fill_manual(values = BATCH_COLORS) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 
  0.03))) +
  ggplot2::labs(x = "Total assigned sgRNA reads (millions)", y = NULL, title = "A. Assigned guide reads") + 
  common +
  ggplot2::theme(legend.position = "none")

# 4. Panel B shows fractions of guides in read-depth bins, separately for
# full-match, mismatch and non-targeting classes. Its denominator is guides.
p_B <- ggplot2::ggplot(b, ggplot2::aes(x = fraction, y = label, fill = count_bin)) +
  ggplot2::geom_col(width = 0.8) + 
  ggplot2::facet_grid(. ~ sgRNA_class) + sample_axis +
  ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1), 
  expand = c(0, 0)) +
  ggplot2::scale_fill_manual(name = "Reads per guide", values = READ_BIN_COLORS, 
  drop = FALSE) +
  ggplot2::labs(x = "Fraction of guides", y = NULL, title = "B. Guide read depth") + 
  common

# 5. Panel C shows fractions of genes by detected full-match-guide coverage
# in each sample. Its denominator is genes, unlike the guide fractions in B.
p_C <- ggplot2::ggplot(c, ggplot2::aes(x = fraction, y = label, fill = guides_bin)) +
  ggplot2::geom_col(width = 0.8) + 
  sample_axis +
  ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1), expand = c(0, 
  0)) +
  ggplot2::scale_fill_manual(name = "Detected full-match guides per gene", values = GENE_COVERAGE_COLORS, 
  drop = FALSE) +
  ggplot2::labs(x = "Fraction of genes", y = NULL, title = "C. Detected full-match guides per gene") + 
  common

if (SHOW_SEGMENT_COUNTS) {
  p_B <- p_B +
  ggplot2::geom_text(ggplot2::aes(label = ifelse(n_guides > 0 & fraction >= MIN_SEGMENT_LABEL_FRACTION, 
    scales::comma(n_guides), "")), position = ggplot2::position_stack(vjust = 0.5), size = 1.65)
  p_C <- p_C +
  ggplot2::geom_text(ggplot2::aes(label = ifelse(n_genes > 0 & fraction >= MIN_SEGMENT_LABEL_FRACTION, 
    scales::comma(n_genes), "")), position = ggplot2::position_stack(vjust = 0.5), size = 1.8)
}

# 6. Panel D shows represented full-match guides per gene at T0, grouped by
# library-design category. The cloning criterion and counts come from task 10;
# they do not select essential genes or impose a fitness-analysis filter here.
p_D <- ggplot2::ggplot(d, ggplot2::aes(x = n_perfect_cloned_10, y = n_genes, fill = design_call)) + 
  ggplot2::geom_col(width = 0.85, color = "black", linewidth = 0.2) +
  ggplot2::geom_text(ggplot2::aes(label = n_genes), 
  vjust = -0.25, size = 2.3) +
  ggplot2::facet_wrap(~design_call, ncol = 3, scales = "free_y", 
  labeller = ggplot2::as_labeller(DESIGN_CALL_LABELS)) +
  ggplot2::scale_fill_manual(values = DESIGN_CALL_COLORS, 
  drop = FALSE) +
  ggplot2::scale_x_continuous(breaks = scales::breaks_width(2)) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 
  0.14))) +
  ggplot2::labs(x = "Cloned full-match guides per gene", y = "Number of genes", title = "D. Full-match guide representation at T0") + 
  common +
  ggplot2::theme(legend.position = "none")

plots <- list(A = p_A, B = p_B, C = p_C, D = p_D)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 7. Display/export each panel with its configured dimensions and formats.
# The source TSVs remain in the shared preprocessing output used by Table S4.
for (panel in names(plots)) {
  plot <- plots[[panel]]
  if (PRINT_IN_RSTUDIO) 
    print(plot)
  height <- if (panel == "D") 
    D_HEIGHT
  else ABC_HEIGHT
  stem <- file.path(OUTPUT_DIR, paste0("Figure_", FIGURE_ID, panel))
  ggplot2::ggsave(paste0(stem, ".svg"), plot, device = svglite::svglite, width = PANEL_WIDTHS[[panel]], 
    height = height, bg = "white")
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, device = grDevices::cairo_pdf, width = PANEL_WIDTHS[[panel]], 
    height = height, bg = "white")
  if (EXPORT_PNG) 
    ggplot2::ggsave(paste0(stem, ".png"), plot, width = PANEL_WIDTHS[[panel]], height = height, 
      dpi = PNG_DPI, bg = "white")
}

message("Saved Figure ", FIGURE_ID, " panels A-D: ", OUTPUT_DIR)


message("Underlying figure data:\n",paste0("  ",normalizePath(file.path(INPUT_DIR,paste0("Figure_S6",c("A","B","C","D"),".tsv"))),collapse="\n"))
