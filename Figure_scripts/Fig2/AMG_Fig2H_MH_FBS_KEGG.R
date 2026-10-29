# Plot selected MH/FBS pathway fitness summaries for Figure 2H.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/Fig2/AMG_Fig2H_MH_FBS_KEGG.R")
# Inputs:
# - preprocessing_scripts/out_14_AMG_MH_FBS_KEGG/Fig2H_input.tsv: pathway medians, paired-gene
#   counts and NES classes.
# Upstream / execution order:
# - preprocessing_scripts/11_AMG_gene_vulnerability.R →
#   preprocessing_scripts/13_AMG_MH_FBS_gene_effects.R →
#   preprocessing_scripts/14_AMG_MH_FBS_KEGG.R.
# Outputs:
# - Figure_scripts/Fig2/out_AMG_Fig2H_MH_FBS_KEGG/ — SVG/PDF/PNG and Fig2H_plotted_pathways.tsv.
# Interpretation / operation:
# - Reuse completed pathway summaries and GSEA results; no enrichment test is rerun.
#
# R packages used/checked: data.table, ggnewscale, ggplot2, svglite.

INPUT_FILE <- "preprocessing_scripts/out_14_AMG_MH_FBS_KEGG/Fig2H_input.tsv"

OUTPUT_DIR <- "Figure_scripts/Fig2/out_AMG_Fig2H_MH_FBS_KEGG"

CONDITION_COLOURS <- c(MH = "#8E8E8E", FBS = "#D8B55A")

NES_COLOURS <- c(`Strong MH` = "#053061", `Moderate MH` = "#4393C3", Neutral = "grey90", `Moderate FBS` = "#F4A582", 
  `Strong FBS` = "#67001F")

FONT_FAMILY <- "Arial"

WIDTH_IN <- 9.5

HEIGHT_IN <- 7.5

DPI <- 300

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
amg_require(c("data.table", "ggplot2", "ggnewscale", "svglite"))

# 2. Read stored pathway effects and enrichment
selected <- amg_read(INPUT_FILE)

amg_required(selected, c("pathway_term", "n_genes_with_MH_FBS_fitness", "median_MH_medLFC", "median_FBS_medLFC", 
  "absolute_fitness_score", "NES_gsea", "NES_class", "Fig2H_rank"), INPUT_FILE)

amg_unique(selected, "pathway_term", INPUT_FILE)

if (!nrow(selected)) stop("No selected Fig. 2H pathways.")

data.table::setorder(selected, -Fig2H_rank)

selected[, `:=`(pathway_label, factor(pathway_term, levels = pathway_term))]

# Check the stored GSEA direction convention before drawing its colour key.
# The ranking was FBS-minus-MH: positive NES favours stronger negative effects
# in MH. Classes describe NES sign/magnitude, not a GSEA p-value cutoff.
expected_class <- as.character(cut(-selected$NES_gsea, breaks = c(-Inf, -1.5, -0.75, 0.75, 1.5, 
  Inf), labels = names(NES_COLOURS)))

if (!identical(is.na(expected_class), is.na(selected$NES_class)) || any(expected_class != selected$NES_class, 
  na.rm = TRUE)) stop("NES class disagrees with the FBS-minus-MH sign convention.")

# 3. Arrange paired MH/FBS values
# Place the two prepared median fitness effects side by side for each pathway.
# Point size uses the number of genes with effects in both media. Pathway
# selection and GSEA were completed in task 14; this script only displays them.
plot_data <- data.table::rbindlist(list(selected[, .(pathway_label, condition = "MH", value = median_MH_medLFC, n_genes = n_genes_with_MH_FBS_fitness)], selected[, .(pathway_label, condition = "FBS", value = median_FBS_medLFC, n_genes = n_genes_with_MH_FBS_fitness)]))

plot_data[, `:=`(condition, factor(condition, levels = names(CONDITION_COLOURS)))]

span <- range(plot_data$value, finite = TRUE)

width <- diff(span)

if (!is.finite(width) || width == 0) width <- 1

x_nes <- span[2] + 0.12 * width

selected[, `:=`(NES_class, factor(NES_class, levels = names(NES_COLOURS)))]

# 4. Draw paired effects and NES categories
p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = pathway_label)) +
  ggplot2::geom_line(ggplot2::aes(x = value, 
  group = pathway_label), colour = "grey75", linewidth = 0.35) +
  ggplot2::geom_point(ggplot2::aes(x = value, 
  fill = condition, size = n_genes), shape = 21, colour = "black", stroke = 0.35, alpha = 0.95) + 
  ggplot2::scale_fill_manual(values = CONDITION_COLOURS, name = "Condition") +
  ggplot2::scale_size_continuous(range = c(2.2, 
  6.5), name = "Genes per pathway") +
  ggnewscale::new_scale_fill() +
  ggplot2::geom_point(data = selected, 
  ggplot2::aes(x = x_nes, y = pathway_label, fill = NES_class), inherit.aes = FALSE, shape = 22, 
  size = 4, colour = "black", stroke = 0.5) +
  ggplot2::scale_fill_manual(values = NES_COLOURS, 
  drop = FALSE, na.value = "white", name = "Relative pathway vulnerability") +
  ggplot2::theme_minimal(base_size = 10, 
  base_family = FONT_FAMILY) +
  ggplot2::theme(legend.position = "right", panel.grid.major.y = ggplot2::element_blank(), 
  panel.grid.minor = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 7.5), 
  axis.text.x = ggplot2::element_text(size = 9)) +
  ggplot2::labs(x = "Median gene fitness effect (log2 fold change versus T0)", 
  y = NULL) +
  ggplot2::coord_cartesian(xlim = c(span[1], x_nes + 0.08 * width), clip = "off")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

files <- file.path(OUTPUT_DIR, paste0("Fig2H_MH_FBS_KEGG", c(".svg", ".pdf", ".png")))

if (!OVERWRITE && any(file.exists(files))) stop("A figure output exists; set OVERWRITE <- TRUE to replace it.")

ggplot2::ggsave(files[1], p, device = svglite::svglite, width = WIDTH_IN, height = HEIGHT_IN, bg = "white")

ggplot2::ggsave(files[2], p, device = grDevices::cairo_pdf, width = WIDTH_IN, height = HEIGHT_IN, 
  bg = "white")

ggplot2::ggsave(files[3], p, width = WIDTH_IN, height = HEIGHT_IN, dpi = DPI, bg = "white")

# Preserve the plotting column in the data referenced by p.
selected_export <- data.table::copy(selected)

selected_export[, `:=`(pathway_label, NULL)]

# Keep the selected pathway statistics/order as the figure's underlying data.
# No pathway medians or enrichment statistics are recalculated at export.
tables <- list(Fig2H_plotted_pathways = selected_export)

amg_write(tables, OUTPUT_DIR, OVERWRITE)

message("Saved Fig. 2H with ", nrow(selected), " pathways to ", OUTPUT_DIR)

if (isTRUE(PRINT_IN_RSTUDIO)) print(p)

data_paths <- normalizePath(file.path(OUTPUT_DIR, paste0(names(tables), ".tsv")), winslash = "/", 
  mustWork = TRUE)

message("Underlying figure data (Fig. 2H):\n", paste0("  ", data_paths, collapse = "\n"))

Fig2H_plot <- p

