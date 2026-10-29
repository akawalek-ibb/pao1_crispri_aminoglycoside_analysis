# Relate untreated T1-consensus severity to sensitivity in Figure 2G.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/Fig2/AMG_Fig2G_severity_sensitivity.R")
# Inputs:
# - preprocessing_scripts/out_12_AMG_MH_fitness/Fig2G_input.tsv: one row per locus, numeric
#   consensus estimates, contributing counts and annotations.
# Upstream / execution order:
# - preprocessing_scripts/11_AMG_gene_vulnerability.R → preprocessing_scripts/12_AMG_MH_fitness.R.
# Outputs:
# - Figure_scripts/Fig2/out_AMG_Fig2G_severity_sensitivity/ — SVG/PDF/PNG plus plotted genes,
#   labels, GAM predictions and statistics TSVs.
# Interpretation / operation:
# - No guide-response curves are refitted; the GAM describes the association between gene
#   descriptors.
#
# R packages used/checked: data.table, ggnewscale, ggplot2, ggrepel, mgcv, scales,
# svglite.

INPUT_FILE <- "preprocessing_scripts/out_12_AMG_MH_fitness/Fig2G_input.tsv"

OUTPUT_DIR <- "Figure_scripts/Fig2/out_AMG_Fig2G_severity_sensitivity"

MIN_SEVERITY_ESTIMATES <- 1L

MIN_SENSITIVITY_ESTIMATES <- 1L

SIZE_POINTS_BY_OPERON <- TRUE

SEVERITY_LABEL_CUTOFF <- -12

SENSITIVITY_LABEL_CUTOFF <- 0.4

MAX_AUTOMATIC_LABELS <- 50L

MANUAL_LABELS <- c("parE", "gyrA", "gyrB", "parC", "rpoB", "murA", "ftsI", "pbpA", "pbpG", "folA", 
  "folP", "fabI", "rpsL", "rplD", "rplC", "lpxC", "bamA", "lptD", "fabF", "mreB")

CATEGORY_COLOURS <- c(`Core essential` = "darkred", `Conditional essential` = "#0072B2", `Essential others` = "#FC8D59", 
  Rest = "darkgrey")

TARGET_LABEL_COLOUR <- "darkviolet"

FONT_FAMILY <- "Arial"

WIDTH_IN <- 8.5

HEIGHT_IN <- 6.5

DPI <- 300

SEED <- 123L

PRINT_IN_RSTUDIO <- TRUE

OVERWRITE <- TRUE

# Local functions for this analysis.
amg_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) 
    stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

amg_required <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) 
    stop(label, " lacks: ", paste(missing, collapse = ", "), call. = FALSE)
}

amg_read <- function(path) {
  if (!file.exists(path)) 
    stop("Missing input: ", path, call. = FALSE)
  con <- if (endsWith(path, ".gz")) 
    gzfile(path, "rt")
  else file(path, "rt")
  on.exit(close(con))
  data.table::as.data.table(utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE, 
    quote = "\"", comment.char = "", na.strings = c("", "NA")))
}

amg_unique <- function(x, columns, label) {
  if (anyNA(x[, columns, with = FALSE]) || anyDuplicated(x[, columns, with = FALSE])) 
    stop(label, ": missing or duplicate key: ", paste(columns, collapse = ", "), call. = FALSE)
  invisible(x)
}

amg_write <- function(tables, directory, overwrite = TRUE, gzip = character()) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(directory, paste0(names(tables), ifelse(names(tables) %in% gzip, ".tsv.gz", 
    ".tsv")))
  if (!overwrite && any(file.exists(paths))) 
    stop("Output exists; set OVERWRITE <- TRUE.")
  for (i in seq_along(tables)) {
    con <- if (endsWith(paths[i], ".gz")) 
      gzfile(paths[i], "wt")
    else file(paths[i], "wt")
    tryCatch(utils::write.table(as.data.frame(tables[[i]]), con, sep = "\t", quote = TRUE, 
      qmethod = "double", row.names = FALSE, na = "NA"), finally = close(con))
  }
  invisible(paths)
}

# Read, prepare, plot and export in order.
amg_require(c("data.table", "ggplot2", "ggnewscale", "ggrepel", "mgcv", "scales", "svglite"))

# 2. Read the annotated consensus used by Table S5
all_genes <- amg_read(INPUT_FILE)

required <- c("PAO1_ID", "gene_name", "severity", "sensitivity", "category_label", "essential_genes_in_operon", 
  "n_control_severity_estimates", "n_control_sensitivity_estimates")

amg_required(all_genes, required, INPUT_FILE)

amg_unique(all_genes, "PAO1_ID", INPUT_FILE)

# 3. Select estimable genes and their labels
# Keep genes with finite consensus measures and enough contributing batch
# estimates. Missing sensitivity is not replaced by a boundary value here;
# the correlation and smooth use this same estimable-gene subset.
x <- all_genes[n_control_severity_estimates >= MIN_SEVERITY_ESTIMATES & n_control_sensitivity_estimates >= MIN_SENSITIVITY_ESTIMATES & is.finite(severity) & is.finite(sensitivity)]

if (nrow(x) < 10L || data.table::uniqueN(x$severity) < 6L) stop("Insufficient data for the original six-basis GAM and correlation.")

if (any(!x$category_label %in% names(CATEGORY_COLOURS))) stop("Unrecognised Fig. 2G category.")

x[, `:=`(category_label, factor(category_label, levels = names(CATEGORY_COLOURS)))]

x[, `:=`(label_txt, ifelse(is.na(gene_name) | !nzchar(gene_name), PAO1_ID, gsub("'", "", gene_name, fixed = TRUE)))]

x[, `:=`(is_known_target = label_txt %in% MANUAL_LABELS, size_operon = log1p(essential_genes_in_operon))]

# Choose labels after selecting the plotted data. Manual targets and a capped
# list of extreme genes are annotated, but unlabelled genes still contribute
# to the scatter, correlation and fitted association.
automatic <- data.table::copy(x[!is_known_target & (severity <= SEVERITY_LABEL_CUTOFF | sensitivity <= SENSITIVITY_LABEL_CUTOFF)])

automatic[, `:=`(label_priority, pmin(data.table::frank(severity, ties.method = "min"), data.table::frank(sensitivity, ties.method = "min")))]

data.table::setorder(automatic, label_priority, sensitivity, severity)

labels <- unique(data.table::rbindlist(list(x[is_known_target == TRUE], head(automatic, MAX_AUTOMATIC_LABELS)), 
  fill = TRUE), by = "PAO1_ID")

labels[, `:=`(label_colour, ifelse(is_known_target, TARGET_LABEL_COLOUR, unname(CATEGORY_COLOURS[as.character(category_label)])))]

# 4. Calculate the displayed GAM and Spearman relationship
set.seed(SEED)

# The GAM describes the cross-gene association between severity and sensitivity.
# It does not refit the guide-level dose-response curves that produced these
# measures. The ribbon below is the fitted mean +/- 1.96 standard errors.
model <- mgcv::gam(sensitivity ~ s(severity, k = 6), data = as.data.frame(x), method = "REML")

model_summary <- summary(model)

correlation <- stats::cor.test(x$severity, x$sensitivity, method = "spearman", exact = FALSE)

prediction <- data.frame(severity = seq(min(x$severity), max(x$severity), length.out = 250L))

pred <- stats::predict(model, newdata = prediction, se.fit = TRUE)

prediction$fit <- as.numeric(pred$fit)

prediction$lower <- prediction$fit - 1.96 * as.numeric(pred$se.fit)

prediction$upper <- prediction$fit + 1.96 * as.numeric(pred$se.fit)

p_text <- function(value) format.pval(value, digits = 2, eps = 2e-16)

annotation <- paste0("Spearman rho = ", round(unname(correlation$estimate), 2), "; p ", ifelse(correlation$p.value < 
  2e-16, "< 2e-16", paste0("= ", p_text(correlation$p.value))), "\nGAM edf = ", round(model_summary$s.table[1, "edf"], 2), "; p ", ifelse(model_summary$s.table[1, "p-value"] < 2e-16, "< 2e-16", paste0("= ", 
  p_text(model_summary$s.table[1, "p-value"]))))

# 5. Draw points, labels, fitted line and interval
p <- ggplot2::ggplot(x, ggplot2::aes(x = severity, y = sensitivity)) +
  ggplot2::geom_vline(xintercept = 0, 
  linewidth = 0.4, colour = "grey85") +
  ggplot2::geom_hline(yintercept = SENSITIVITY_LABEL_CUTOFF, 
  linetype = "dotted", linewidth = 0.8, colour = "grey35")

if (SIZE_POINTS_BY_OPERON) {
  p <- p +
  ggplot2::geom_point(ggplot2::aes(colour = category_label, size = size_operon), alpha = 0.85) + 
    ggplot2::scale_size_continuous(range = c(0.25, 1.5), name = "Log-scaled essential\ngenes in operon", 
      breaks = log1p(c(0, 1, 2, 4, 6)), labels = c("0", "1", "2", "4", "6"))
} else p <- p +
  ggplot2::geom_point(ggplot2::aes(colour = category_label), size = 1.1, alpha = 0.85)

p <- p +
  ggplot2::scale_colour_manual(values = CATEGORY_COLOURS, breaks = names(CATEGORY_COLOURS), 
  drop = FALSE, name = "Essentiality annotation") +
  ggnewscale::new_scale_colour() +
  ggrepel::geom_text_repel(data = labels, 
  ggplot2::aes(label = label_txt, colour = label_colour), size = 2, max.overlaps = Inf, box.padding = 0.4, 
  point.padding = 0.2, min.segment.length = 0, segment.color = "grey60", segment.alpha = 0.4, 
  show.legend = FALSE, fontface = "italic", seed = SEED) +
  ggplot2::scale_colour_identity() + 
  ggplot2::geom_ribbon(data = prediction, ggplot2::aes(x = severity, ymin = lower, ymax = upper), 
    inherit.aes = FALSE, fill = "grey80", alpha = 0.4) +
  ggplot2::geom_line(data = prediction, 
  ggplot2::aes(x = severity, y = fit), inherit.aes = FALSE, colour = "black", linewidth = 0.8) + 
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(6), 
  limits = c(0.05, 1)) +
  ggplot2::annotate("text", x = -Inf, y = Inf, label = annotation, hjust = -0.05, 
  vjust = 1.2, size = 2.5, family = FONT_FAMILY) +
  ggplot2::theme_minimal(base_size = 8, base_family = FONT_FAMILY) + 
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_line(linewidth = 0.25), 
    legend.position = "right", legend.box = "vertical", legend.title = ggplot2::element_text(face = "bold"), 
    legend.text = ggplot2::element_text(size = 7.5), axis.title = ggplot2::element_text(face = "bold", 
      size = 8)) +
  ggplot2::labs(x = "Severity (median full-match-guide adjusted log2 fold change)", 
  y = "Sensitivity (predicted efficacy at sustained log2FC <= -1)", subtitle = "Untreated T1 consensus across B1-B4")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

files <- file.path(OUTPUT_DIR, paste0("Fig2G_severity_sensitivity", c(".svg", ".pdf", ".png")))

if (!OVERWRITE && any(file.exists(files))) stop("A figure output exists; set OVERWRITE <- TRUE to replace it.")

ggplot2::ggsave(files[1], p, device = svglite::svglite, width = WIDTH_IN, height = HEIGHT_IN, bg = "white")

ggplot2::ggsave(files[2], p, device = grDevices::cairo_pdf, width = WIDTH_IN, height = HEIGHT_IN, 
  bg = "white")

ggplot2::ggsave(files[3], p, width = WIDTH_IN, height = HEIGHT_IN, dpi = DPI, bg = "white")

statistics <- data.table::data.table(metric = c("n_genes", "Spearman_rho", "Spearman_p", "GAM_edf", 
  "GAM_p"), value = c(nrow(x), unname(correlation$estimate), correlation$p.value, model_summary$s.table[1, "edf"], model_summary$s.table[1, "p-value"]))

# Export the selected genes, label subset, association statistics and smooth
# predictions as separate tables; the label subset is not the analysis denominator.
tables <- list(Fig2G_plotted_genes = x, Fig2G_labels = labels, Fig2G_statistics = statistics, Fig2G_fitted_curve = prediction)

amg_write(tables, OUTPUT_DIR, OVERWRITE)

message("Saved Fig. 2G with ", nrow(x), " genes to ", OUTPUT_DIR)

if (isTRUE(PRINT_IN_RSTUDIO)) print(p)

data_paths <- normalizePath(file.path(OUTPUT_DIR, paste0(names(tables), ".tsv")), winslash = "/", 
  mustWork = TRUE)

message("Underlying figure data (Fig. 2G):\n", paste0("  ", data_paths, collapse = "\n"))

Fig2G_plot <- p

