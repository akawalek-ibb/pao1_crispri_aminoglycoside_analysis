# Plot selected GO overrepresentation terms for Figure 6D.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/Fig6/AMG_Fig6D_selected_GO.R")
# Inputs:
# - preprocessing_scripts/out_22_AMG_RNAseq_selected_GO/selected_GO_results.tsv and
#   Analysis_settings.tsv: both Up/Down results and display order.
# Upstream / execution order:
# - preprocessing_scripts/18_AMG_RNAseq_results.R → preprocessing_scripts/19_AMG_RNAseq_GO_ORA.R
#   → preprocessing_scripts/22_AMG_RNAseq_selected_GO.R.
# Outputs:
# - Figure_scripts/Fig6/out_AMG_Fig6D_selected_GO/ — PNG/PDF/SVG and plotted
#   term/direction/value TSVs.
# Interpretation / operation:
# - Enrichment and BH adjustment are upstream; this script only prepares the configured display.
#
# R packages used/checked: data.table, ggplot2, svglite.

INPUT_DIR <- "preprocessing_scripts/out_22_AMG_RNAseq_selected_GO"

OUTPUT_DIR <- "Figure_scripts/Fig6/out_AMG_Fig6D_selected_GO"

COMPARISON_LABELS <- c(atpI = "g_atpI1-m2", accC = "g_accC1-m3", secA = "g_secA1-m3", ftsZ = "g_ftsZ1-m4")

# Use a diamond when both Up and Down results are significant for the same term.
SHOW_BOTH_DIRECTIONS_DIAMOND <- TRUE

# Appearance and output size; these settings do not change the enrichment results.
POINT_SIZE_RANGE <- c(2.5, 8)

TERM_WRAP <- 48L

FONT_FAMILY <- "Arial"

WIDTH_IN <- 10

HEIGHT_IN <- 6.5

DPI <- 300

PRINT_IN_RSTUDIO <- TRUE

OVERWRITE <- TRUE

# Check that the required R packages are installed before starting.
amg_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) 
    stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

# Check that an input table contains the columns needed below.
amg_required <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) 
    stop(label, " lacks: ", paste(missing, collapse = ", "), call. = FALSE)
}

# Read a tab-separated table, including compressed tables ending in .gz.
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

# Stop if a gene, comparison or replicate is missing its identifier or appears twice.
amg_unique <- function(x, columns, label) {
  if (anyNA(x[, columns, with = FALSE]) || anyDuplicated(x[, columns, with = FALSE])) 
    stop(label, ": missing or duplicate key: ", paste(columns, collapse = ", "), call. = FALSE)
  invisible(x)
}

# Save the data behind the figure as tab-separated tables.
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

# Record the settings, input filenames and file checksums so the run can be traced.
amg_settings <- function(settings, inputs) {
  data.table::data.table(setting = c(names(settings), paste0("input:", names(inputs)), paste0("md5:", 
    names(inputs)), "R_version"), value = c(as.character(settings), unname(inputs), unname(tools::md5sum(inputs)), 
    R.version.string))
}

# Read the setting recorded by the upstream analysis.
rnaseq_setting <- function(directory, name, numeric = FALSE) {
  x <- amg_read(file.path(directory, "Analysis_settings.tsv"))
  value <- x$value[x$setting == name]
  if (length(value) != 1L || is.na(value)) 
    stop("Missing/ambiguous setting: ", name)
  if (numeric) 
    as.numeric(value)
  else as.character(value)
}

# Read saved TRUE/FALSE flags consistently and reject unclear entries.
rf_logical <- function(x, label) {
  x <- tolower(trimws(as.character(x)))
  if (anyNA(x) || any(!x %in% c("true", "false", "1", "0"))) 
    stop("Invalid logical values: ", label)
  x %in% c("true", "1")
}

# First build all figure files in a temporary folder, then copy the finished files
# into the output folder. Save PNG, PDF and SVG, plus any accompanying data tables.
rf_export <- function(plot, tables, output, stem, width, height, dpi, overwrite = TRUE, grob = FALSE) {
  paths <- c(paste0(stem, c(".png", ".pdf", ".svg")), if (length(tables)) paste0(names(tables), 
    ".tsv"))
  if (!overwrite && any(file.exists(file.path(output, paths)))) 
    stop("An output exists; set OVERWRITE <- TRUE.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(".rnaseq_figure_", tmpdir = output)
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE))
  if (!grob) {
    ggplot2::ggsave(file.path(stage, paths[1L]), plot, width = width, height = height, dpi = dpi, 
      bg = "white")
    ggplot2::ggsave(file.path(stage, paths[2L]), plot, device = grDevices::cairo_pdf, width = width, 
      height = height, bg = "white")
    ggplot2::ggsave(file.path(stage, paths[3L]), plot, device = svglite::svglite, width = width, 
      height = height, bg = "white")
  } else {
    grDevices::png(file.path(stage, paths[1L]), width = width, height = height, units = "in", 
      res = dpi, type = "cairo")
    tryCatch(grid::grid.draw(plot), finally = grDevices::dev.off())
    grDevices::cairo_pdf(file.path(stage, paths[2L]), width = width, height = height)
    tryCatch(grid::grid.draw(plot), finally = grDevices::dev.off())
    svglite::svglite(file.path(stage, paths[3L]), width = width, height = height)
    tryCatch(grid::grid.draw(plot), finally = grDevices::dev.off())
  }
  if (length(tables)) 
    amg_write(tables, stage, TRUE)
  for (name in paths) if (!file.copy(file.path(stage, name), file.path(output, name), overwrite = overwrite)) 
    stop("Could not save ", name)
  invisible(plot)
}

# Format guide labels so the gene name is italic and the guide number stays upright.
rf_guide_plotmath <- function(x) {
  out <- vapply(x, function(label) {
    m <- regmatches(label, regexec("^g_(.+?)([0-9]+)(-m[0-9]+)?$", label, perl = TRUE))[[1L]]
    if (length(m)) 
      paste0("paste('g_',italic('", m[2L], "'), '", m[3L], m[4L], "')")
    else paste0("'", gsub("'", "", label, fixed = TRUE), "'")
  }, character(1))
  parse(text = out)
}

# 1. Check packages and locate the prepared results.
amg_require(c("data.table", "ggplot2", "svglite"))

paths <- file.path(INPUT_DIR, c("selected_GO_results.tsv", "Analysis_settings.tsv"))

# 2. Read the selected terms and their separate results for upregulated and
# downregulated genes. Require one record per term and knockdown.
data <- amg_read(paths[1L])

amg_required(data, c("comparison_id", "ontology", "GO_ID", "GO_term", "display_order", "Up_FDR", 
  "Down_FDR", "Up_GeneRatio", "Down_GeneRatio", "displayed_direction", "displayed_signed_log10_FDR", 
  "both_directions_significant"), paths[1L])

amg_unique(data, c("comparison_id", "GO_ID"), paths[1L])

# Check that the input contains exactly the four planned knockdown comparisons.
if (!setequal(unique(data$comparison_id), names(COMPARISON_LABELS))) stop("Selected GO comparisons disagree with the four planned knockdowns.")

data[, `:=`(both_directions_significant, rf_logical(both_directions_significant, "both_directions_significant"))]

# 3. Use the direction with the smaller adjusted p-value for the bubble size.
# If the values tie, use Up. Keep both directions in the exported table.
data[, `:=`(available = is.finite(Up_FDR) | is.finite(Down_FDR), prefer_up = is.finite(Up_FDR) & (!is.finite(Down_FDR) | Up_FDR <= Down_FDR))]

data[, `:=`(dominant_FDR = ifelse(prefer_up, Up_FDR, Down_FDR), dominant_GeneRatio_text = ifelse(prefer_up, Up_GeneRatio, Down_GeneRatio))]

# Turn a stored fraction such as 12/100 into a number for bubble sizing.
# This is the fraction of GO-annotated DEGs assigned to the term in that direction.
ratio <- function(x) vapply(x, function(value) {
  if (is.na(value)) 
    return(NA_real_)
  split <- strsplit(value, "/", fixed = TRUE)[[1L]]
  if (length(split) != 2L) 
    stop("Invalid GeneRatio: ", value)
  number <- suppressWarnings(as.numeric(split))
  if (any(!is.finite(number)) || number[2L] <= 0 || number[1L] < 0 || number[1L] > number[2L]) 
    stop("Invalid GeneRatio: ", value)
  number[1L]/number[2L]
}, numeric(1))

data[, `:=`(gene_ratio, ratio(dominant_GeneRatio_text))]

if (any(data$available & !is.finite(data$gene_ratio))) stop("Returned selected GO result lacks a valid GeneRatio.")

# Use the selected-term display order and check that no term/comparison is missing.
terms <- unique(data[, .(GO_ID, GO_term, ontology, display_order)])

amg_unique(terms, "GO_ID", "Selected GO metadata")

if (nrow(data) != nrow(terms) * length(COMPARISON_LABELS)) stop("Incomplete selected GO term/comparison grid.")

data.table::setorder(terms, display_order)

# Wrap long term names for readability; their saved names stay unchanged.
terms[, `:=`(label, vapply(GO_term, function(x) paste(strwrap(x, width = TERM_WRAP), collapse = "\n"), character(1)))]

data[, `:=`(comparison = factor(comparison_id, levels = names(COMPARISON_LABELS)), term = factor(terms$label[match(GO_ID, terms$GO_ID)], levels = rev(terms$label)), ontology = factor(ontology, levels = c("BP", "MF", 
  "CC")), display_shape = ifelse(SHOW_BOTH_DIRECTIONS_DIAMOND & both_directions_significant, 
  "Both directions", "Single direction"))]

# Read the colour limit chosen in preprocessing step 22.
cap <- rnaseq_setting(INPUT_DIR, "DISPLAY_SCORE_CAP", numeric = TRUE)

# 4. Draw one bubble per term and knockdown. Bubble size shows the gene fraction.
# Colour uses the saved score: red for Up enrichment, blue for Down enrichment,
# and zero for returned nonsignificant results. A cross marks an unavailable result.
# A diamond marks significant enrichment in both directions when this option is on.
# BP, MF and CC group biological processes, molecular functions and cellular components.
p <- ggplot2::ggplot(data, ggplot2::aes(x = comparison, y = term)) +
  ggplot2::geom_blank() +
  ggplot2::geom_point(data = data[available == FALSE], shape = 4, size = 1.8, stroke = 0.4, colour = "grey45") +
  ggplot2::geom_point(data = data[available == TRUE], ggplot2::aes(size = gene_ratio, fill = displayed_signed_log10_FDR, shape = display_shape), 
  colour = "black", stroke = 0.35) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "grey92", 
  high = "#B2182B", midpoint = 0, limits = c(-cap, cap), name = "Signed -log10(FDR)\nDown                 Up") + 
  ggplot2::scale_size_continuous(range = POINT_SIZE_RANGE, name = "Gene ratio") +
  ggplot2::scale_shape_manual(values = c(`Single direction` = 21, 
  `Both directions` = 23), guide = "none") +
  ggplot2::scale_x_discrete(labels = rf_guide_plotmath(unname(COMPARISON_LABELS)), 
  drop = FALSE) +
  ggplot2::facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
  ggplot2::theme_minimal(base_size = 10, 
  base_family = FONT_FAMILY) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "grey95", 
  colour = NA), strip.text.y = ggplot2::element_text(angle = 0, face = "bold"), axis.text.x = ggplot2::element_text(angle = 45, 
  hjust = 1, colour = "black"), axis.text.y = ggplot2::element_text(colour = "black")) +
  ggplot2::labs(x = NULL, 
  y = NULL)

# 5. Save the plotted values and settings alongside PNG, PDF and SVG figures.
settings <- amg_settings(c(direction = "Smaller stored directional FDR; Up on exact ties", score = "Stored task-22 signed score; returned nonsignificant = 0, unreturned = NA", 
  SHOW_BOTH_DIRECTIONS_DIAMOND = SHOW_BOTH_DIRECTIONS_DIAMOND, DISPLAY_SCORE_CAP = cap), stats::setNames(paths, 
  basename(paths)))

rf_export(p, list(Fig6D_plotted_values = data, Analysis_settings = settings), OUTPUT_DIR, "Fig6D_selected_GO", 
  WIDTH_IN, HEIGHT_IN, DPI, OVERWRITE)

if (PRINT_IN_RSTUDIO) print(p)

message("Saved Fig. 6D: ", nrow(terms), " GO terms, ", nrow(data), " term/knockdown values.")


# Print the locations of the saved figure-data tables.
figure_data_files <- file.path(OUTPUT_DIR,c("Fig6D_plotted_values.tsv","Analysis_settings.tsv"))
message("Underlying figure data:\n",paste0("  ",normalizePath(figure_data_files,winslash="/"),collapse="\n"))
