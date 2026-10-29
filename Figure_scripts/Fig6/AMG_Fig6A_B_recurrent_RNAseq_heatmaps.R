# Show RNA-seq changes among recurrent sensitizing/protective screen genes.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/Fig6/AMG_Fig6A_B_recurrent_RNAseq_heatmaps.R")
# Inputs:
# - preprocessing_scripts/out_03_AMG_relevant_results/ — recurrence and settings;
#   preprocessing_scripts/out_18_AMG_RNAseq_results/ — complete RNA-seq results and manifest.
# Upstream / execution order:
# - preprocessing_scripts/03_AMG_relevant_results.R and
#   preprocessing_scripts/18_AMG_RNAseq_results.R.
# Outputs:
# - Figure_scripts/Fig6/out_AMG_Fig6A_B_recurrent_RNAseq_heatmaps/ — PNG/PDF/SVG heatmaps and
#   plotted gene/value/statistics TSVs.
# Interpretation / operation:
# - Default masking uses FDR < 0.05 alone; enabling FDR_and_effect additionally requires
#   |log2FC| > 1.
#
# R packages used/checked: data.table, ggplot2, pheatmap, svglite.

RNA_DIR <- "preprocessing_scripts/out_18_AMG_RNAseq_results"

INTERACTION_DIR <- "preprocessing_scripts/out_03_AMG_relevant_results"

OUTPUT_DIR <- "Figure_scripts/Fig6/out_AMG_Fig6A_B_recurrent_RNAseq_heatmaps"

# Set the displayed order and names of the five RNA-seq comparisons.
COMPARISON_LABELS <- c(TOB = "TOB alone", accC = "g_accC1-m3", secA = "g_secA1-m3", ftsZ = "g_ftsZ1-m4", 
  atpI = "g_atpI1-m2")

# Use "FDR_only" to colour changes with FDR < 0.05.
# "FDR_and_effect" additionally requires |log2FC| > RNASEQ_ABS_LOG2FC_CUT.
RNA_MASK_MODE <- "FDR_only"

RNASEQ_FDR_CUT <- 0.05

# This effect-size cutoff is used only in FDR_and_effect mode.
RNASEQ_ABS_LOG2FC_CUT <- 1

FILTER_GENES_WITH_ANY_SIGNIFICANT_RNASEQ <- TRUE

# Leave NULL for automatic selection. A supplied TSV needs direction and PAO1_ID columns.
MANUAL_GENES_FILE <- NULL

# Values outside this range share the end colour; original values are still exported.
LFC_DISPLAY_CAP <- 2.32

# Group genes by similarity of their displayed expression patterns.
CLUSTER_GENES <- TRUE

CLUSTER_DISTANCE <- "euclidean"

CLUSTER_METHOD <- "ward.D2"

# Appearance and output size. Hidden cells use this grey colour.
NS_FILL <- "#D9D9D9"

FONT_FAMILY <- "Arial"

GENE_FONT_SIZE <- 4

CONTRAST_FONT_SIZE <- 8

WIDTH_IN <- 10

HEIGHT_IN <- 3.8

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

# Use a gene name when available; otherwise show its PAO1 locus identifier.
rf_label <- function(gene, id) {
  gene <- trimws(gsub("'", "", as.character(gene), fixed = TRUE))
  gene[is.na(gene) | !nzchar(gene)] <- id[is.na(gene) | !nzchar(gene)]
  gene
}

# Load the complete RNA-seq results and their comparison information.
# Check that the tables refer to the same comparisons.
rf_rna <- function(directory) {
  paths <- file.path(directory, c("RNAseq_results.tsv.gz", "Comparison_manifest.tsv", "RNAseq_summary.tsv", 
    "Analysis_settings.tsv"))
  names(paths) <- c("results", "manifest", "summary", "settings")
  x <- lapply(paths, amg_read)
  amg_required(x$results, c("PAO1_ID", "gene_name", "comparison_id", "logFC", "FDR", "DEG_direction"), 
    paths[1L])
  amg_unique(x$results, c("PAO1_ID", "comparison_id"), paths[1L])
  amg_unique(x$manifest, "comparison_id", paths[2L])
  amg_unique(x$summary, "comparison_id", paths[3L])
  if (!setequal(unique(x$results$comparison_id), x$manifest$comparison_id)) 
    stop("Task-18 results and manifest disagree.")
  x$paths <- paths
  x
}

# Apply the chosen font to the heatmap and italicize its gene labels.
rf_font <- function(grob, family, italic_gene_labels = FALSE) {
  if (inherits(grob, "text")) {
    gp <- as.list(grob$gp)
    gp$fontfamily <- family
    grob$gp <- do.call(grid::gpar, gp)
  }
  if (inherits(grob, "gtable")) {
    grob$grobs <- lapply(grob$grobs, rf_font, family = family)
    if (italic_gene_labels) {
      i <- which(grob$layout$name == "col_names")
      if (length(i) != 1L) 
        stop("Cannot locate heatmap gene labels.")
      gp <- as.list(grob$grobs[[i]]$gp)
      gp$font <- NULL
      gp$fontface <- "italic"
      gp$col <- "black"
      grob$grobs[[i]]$gp <- do.call(grid::gpar, gp)
    }
  }
  if (inherits(grob, "gTree") && length(grob$children)) 
    grob$children <- do.call(grid::gList, lapply(grob$children, rf_font, family = family))
  grob
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

# Prepare one panel at a time: sensitizing genes or protective genes.
prepare_heatmap <- function(rna, recurrent, direction) {
  # Start with the genes already called recurrent by preprocessing step 03.
  field <- paste0("recurrent_", tolower(direction))
  meta <- data.table::copy(recurrent[recurrent[[field]]])
  n_total <- nrow(meta)
  # Keep only genes represented in the selected RNA-seq comparisons.
  present <- unique(rna$PAO1_ID[rna$comparison_id %in% names(COMPARISON_LABELS)])
  meta <- meta[PAO1_ID %in% present]
  if (!nrow(meta)) 
    stop("No recurrent ", direction, " genes represented in RNA-seq.")
  # Arrange each gene's expression changes and adjusted p-values in matching tables.
  effect <- matrix(NA_real_, nrow(meta), length(COMPARISON_LABELS), dimnames = list(meta$PAO1_ID, 
    names(COMPARISON_LABELS)))
  fdr <- effect
  z <- rna[PAO1_ID %in% meta$PAO1_ID & comparison_id %in% names(COMPARISON_LABELS)]
  ij <- cbind(match(z$PAO1_ID, meta$PAO1_ID), match(z$comparison_id, names(COMPARISON_LABELS)))
  effect[ij] <- z$logFC
  fdr[ij] <- z$FDR
  # Decide which cells will be coloured. The default uses FDR alone.
  # The optional second rule also requires an absolute log2 fold change above the cutoff.
  fdr_only <- is.finite(effect) & is.finite(fdr) & fdr < RNASEQ_FDR_CUT
  fdr_effect <- fdr_only & abs(effect) > RNASEQ_ABS_LOG2FC_CUT
  significant <- if (RNA_MASK_MODE == "FDR_only") 
    fdr_only
  else fdr_effect
  # Count how the two display rules differ, without changing the original RNA-seq results.
  audit <- data.table::data.table(direction = direction, n_recurrent = n_total, n_RNAseq_represented = nrow(meta), 
    n_genes_any_FDR = sum(rowSums(fdr_only) > 0), n_genes_any_FDR_and_effect = sum(rowSums(fdr_effect) > 
      0), n_cells_FDR_only = sum(fdr_only), n_cells_FDR_and_effect = sum(fdr_effect))
  # By default, show a gene only if at least one of the five displayed comparisons
  # passes the selected rule. This includes the TOB-alone comparison.
  keep <- if (FILTER_GENES_WITH_ANY_SIGNIFICANT_RNASEQ) 
    rowSums(significant) > 0
  else rep(TRUE, nrow(meta))
  # If a gene list is supplied, use that list instead of the automatic gene selection.
  # It must contain recurrent genes present in the RNA-seq data.
  if (!is.null(MANUAL_GENES_FILE)) {
    manual <- amg_read(MANUAL_GENES_FILE)
    amg_required(manual, c("direction", "PAO1_ID"), MANUAL_GENES_FILE)
    amg_unique(manual, c("direction", "PAO1_ID"), MANUAL_GENES_FILE)
    wanted <- manual$PAO1_ID[manual$direction == direction]
    if (!length(wanted) || length(setdiff(wanted, meta$PAO1_ID))) 
      stop("Manual heatmap selection contains absent/nonrecurrent genes for ", direction)
    keep <- meta$PAO1_ID %in% wanted
  }
  meta <- meta[keep]
  effect <- effect[keep, , drop = FALSE]
  fdr <- fdr[keep, , drop = FALSE]
  significant <- significant[keep, , drop = FALSE]
  if (!nrow(meta)) 
    stop("No genes remain for ", direction, " under ", RNA_MASK_MODE)
  # Limit extreme values for colour display and grouping only.
  # The original, uncapped fold changes are retained for export.
  clipped <- pmax(pmin(effect, LFC_DISPLAY_CAP), -LFC_DISPLAY_CAP)
  # For grouping similar genes, temporarily treat hidden or missing cells as zero.
  # This is a plotting choice, not a replacement for the measured expression change.
  cluster_values <- clipped
  cluster_values[!significant | !is.finite(cluster_values)] <- 0
  tree <- if (CLUSTER_GENES && nrow(meta) > 1L) 
    stats::hclust(stats::dist(cluster_values, method = CLUSTER_DISTANCE), method = CLUSTER_METHOD)
  else NULL
  order <- if (is.null(tree)) 
    seq_len(nrow(meta))
  else tree$order
  # Hide cells that fail the display rule; these will appear grey in the heatmap.
  shown <- clipped
  shown[!significant] <- NA_real_
  meta[, `:=`(display_order, match(seq_len(.N), order))]
  audit[, `:=`(n_plotted = nrow(meta), mask_mode = RNA_MASK_MODE)]
  list(metadata = meta, effect = effect, fdr = fdr, significant = significant, shown = shown, 
    clustering_values = cluster_values, tree = tree, order = order, audit = audit)
}

# 1. Check packages and the user settings before reading any data.
amg_require(c("data.table", "pheatmap", "svglite"))

if (!RNA_MASK_MODE %in% c("FDR_only", "FDR_and_effect") || RNASEQ_FDR_CUT <= 0 || RNASEQ_FDR_CUT > 
  1 || RNASEQ_ABS_LOG2FC_CUT < 0 || LFC_DISPLAY_CAP <= 0) stop("Invalid heatmap settings.")

# 2. Read the saved RNA-seq results and the recurrent interaction calls.
input <- rf_rna(RNA_DIR)

if (length(setdiff(names(COMPARISON_LABELS), input$manifest$comparison_id))) stop("Required RNA-seq comparisons are missing.")

recurrence_file <- file.path(INTERACTION_DIR, "AG_gene_recurrence.tsv")

recurrent <- amg_read(recurrence_file)

amg_required(recurrent, c("PAO1_ID", "gene_name", "recurrent_sensitizing", "recurrent_protective", 
  "n_independent_T1_sensitizing", "n_independent_T1_protective"), recurrence_file)

amg_unique(recurrent, "PAO1_ID", recurrence_file)

# Use the upstream recurrence cutoff and check that its counts agree with its flags.
# This script does not choose a new recurrence definition.
recurrence_min <- rnaseq_setting(INTERACTION_DIR, "RECURRENCE_MIN", numeric = TRUE)

for (direction in c("sensitizing", "protective")) {
  flag <- paste0("recurrent_", direction)
  count <- paste0("n_independent_T1_", direction)
  data.table::set(recurrent, j = flag, value = rf_logical(recurrent[[flag]], flag))
  if (any(recurrent[[flag]] != (recurrent[[count]] >= recurrence_min))) 
    stop("Stored recurrence flags/counts disagree.")
}

plots <- list()

values <- orders <- audits <- list()

# 3. Build separate panels for sensitizing and protective genes.
# The comparison order is fixed; only the gene order can be rearranged.
for (direction in c("Sensitizing", "Protective")) {
  panel <- if (direction == "Sensitizing") 
    "Fig6A"
  else "Fig6B"
  result <- prepare_heatmap(input$results, recurrent, direction)
  labels <- rf_label(result$metadata$gene_name, result$metadata$PAO1_ID)
  mat <- t(result$shown)
  rownames(mat) <- unname(COMPARISON_LABELS)
  p <- pheatmap::pheatmap(mat, cluster_rows = FALSE, cluster_cols = if (is.null(result$tree)) 
    FALSE
  else result$tree, scale = "none", color = (grDevices::colorRampPalette(c("#053061", "white", 
    "#67001F")))(100), breaks = seq(-LFC_DISPLAY_CAP, LFC_DISPLAY_CAP, length.out = 101), na_col = NS_FILL, 
    border_color = "white", labels_col = labels, fontsize_row = CONTRAST_FONT_SIZE, fontsize_col = GENE_FONT_SIZE, 
    angle_col = 90, treeheight_col = 17, fontfamily = "sans", main = paste0(direction, " recurrent genes (", 
      nrow(result$metadata), "/", result$audit$n_recurrent, ")"), legend_breaks = c(-LFC_DISPLAY_CAP, 
      0, LFC_DISPLAY_CAP), legend_labels = c(paste0("<= -", LFC_DISPLAY_CAP), "0", paste0(">= ", 
      LFC_DISPLAY_CAP)), silent = TRUE)
  p$gtable <- rf_font(p$gtable, FONT_FAMILY, italic_gene_labels = TRUE)
  ij <- cbind(rep(seq_len(nrow(result$metadata)), each = ncol(result$effect)), rep(seq_len(ncol(result$effect)), 
    times = nrow(result$metadata)))
  # Keep the original value, its FDR and its displayed value together for each cell.
  cells <- data.table::data.table(panel = panel, direction = direction, PAO1_ID = result$metadata$PAO1_ID[ij[, 1L]], gene_name = labels[ij[, 1L]], display_order = result$metadata$display_order[ij[, 1L]], comparison_id = names(COMPARISON_LABELS)[ij[, 2L]], logFC = result$effect[ij], FDR = result$fdr[ij], 
    significant = result$significant[ij], displayed_value = result$shown[ij])
  data.table::setorder(cells, display_order, comparison_id)
  source_paths <- c(input$paths, recurrence = recurrence_file)
  if (!is.null(MANUAL_GENES_FILE)) 
    source_paths <- c(source_paths, manual_genes = MANUAL_GENES_FILE)
  settings <- amg_settings(c(RNA_MASK_MODE = RNA_MASK_MODE, RNASEQ_FDR_CUT = RNASEQ_FDR_CUT, 
    RNASEQ_ABS_LOG2FC_CUT = RNASEQ_ABS_LOG2FC_CUT, FILTER_GENES_WITH_ANY_SIGNIFICANT_RNASEQ = FILTER_GENES_WITH_ANY_SIGNIFICANT_RNASEQ, 
    LFC_DISPLAY_CAP = LFC_DISPLAY_CAP, clustering = paste(CLUSTER_DISTANCE, CLUSTER_METHOD), 
    masked_values_for_clustering = 0, RECURRENCE_MIN = recurrence_min), source_paths)
  rf_export(p$gtable, list(), OUTPUT_DIR, paste0(panel, "_", tolower(direction), "_RNAseq"), 
    WIDTH_IN, HEIGHT_IN, DPI, OVERWRITE, grob = TRUE)
  values[[direction]] <- cells
  orders[[direction]] <- result$metadata[result$order]
  audits[[direction]] <- result$audit
  plots[[direction]] <- c(result, list(plot = p))
  if (PRINT_IN_RSTUDIO) {
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  }
  message(panel, ": ", nrow(result$metadata), " / ", result$audit$n_recurrent, " recurrent genes; ", 
    RNA_MASK_MODE, ".")
}

# 4. Save the plotted values, gene order, display-rule counts and run settings.
amg_write(list(Fig6AB_plotted_values = data.table::rbindlist(values), Fig6AB_gene_order = data.table::rbindlist(orders, 
  idcol = "direction"), Fig6AB_mask_comparison = data.table::rbindlist(audits), Analysis_settings = settings), 
  OUTPUT_DIR, OVERWRITE)


# Print the locations of the saved figure-data tables.
figure_data_files <- file.path(OUTPUT_DIR,c("Fig6AB_plotted_values.tsv","Fig6AB_gene_order.tsv","Fig6AB_mask_comparison.tsv","Analysis_settings.tsv"))
message("Underlying figure data:\n",paste0("  ",normalizePath(figure_data_files,winslash="/"),collapse="\n"))
