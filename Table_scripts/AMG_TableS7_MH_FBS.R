# Assemble Supplementary Table S7 from completed analysis exports.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Table_scripts/AMG_TableS7_MH_FBS.R")
# Inputs:
# - preprocessing_scripts/out_13_AMG_MH_FBS_gene_effects/;
#   preprocessing_scripts/out_14_AMG_MH_FBS_KEGG/;
#   preprocessing_scripts/out_15_AMG_MH_FBS_MAGeCK/ — completed gene/pathway/significance tables
#   and settings.
# - Table_scripts/TableS7_legend.tsv: A header row with field and description columns;
#   {{tokens}} are filled from analysis settings.
# Upstream / execution order:
# - preprocessing_scripts/11_AMG_gene_vulnerability.R →
#   preprocessing_scripts/13_AMG_MH_FBS_gene_effects.R →
#   preprocessing_scripts/14_AMG_MH_FBS_KEGG.R; preprocessing_scripts/15_AMG_MH_FBS_MAGeCK.R
#   independently imports the included MAGeCK catalogue.
# Outputs:
# - Table_scripts/out_AMG_TableS7_MH_FBS/Supplementary_Table_S7.xlsx
# Interpretation / operation:
# - Requires openxlsx and data.table.
# - This script formats the workbook; it does not refit models or redefine statistical calls.

GENE_INPUT_DIR <- "preprocessing_scripts/out_13_AMG_MH_FBS_gene_effects"

KEGG_INPUT_DIR <- "preprocessing_scripts/out_14_AMG_MH_FBS_KEGG"

MAGECK_INPUT_DIR <- "preprocessing_scripts/out_15_AMG_MH_FBS_MAGeCK"

LEGEND_FILE <- "Table_scripts/TableS7_legend.tsv"

OUTPUT_FILE <- "Table_scripts/out_AMG_TableS7_MH_FBS/Supplementary_Table_S7.xlsx"

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

amg_setting <- function(directory, name) {
  settings <- amg_read(file.path(directory, "Analysis_settings.tsv"))
  found <- settings$value[settings$setting == name]
  if (length(found) != 1L || is.na(found)) 
    stop("Missing or ambiguous preprocessing setting: ", name)
  as.character(found)
}

# Read the editable legend rather than embedding the explanatory text here.
# Replace {{placeholders}} with counts/settings from the prepared results;
# stop if a placeholder cannot be resolved, so the legend cannot silently
# report a threshold or sample count that was not supplied.
amg_legend <- function(path, tokens) {
  x <- as.data.frame(amg_read(path))
  amg_required(x, c("field", "description"), path)
  x <- x[c("field", "description")]
  x[is.na(x)] <- ""
  for (name in names(tokens)) for (column in names(x)) x[[column]] <- gsub(paste0("{{", name, 
    "}}"), as.character(tokens[[name]]), x[[column]], fixed = TRUE)
  unresolved <- unique(unlist(regmatches(unlist(x), gregexpr("\\{\\{[^}]+\\}\\}", unlist(x)))))
  if (length(unresolved)) 
    stop("Unknown legend token(s): ", paste(unresolved, collapse = ", "))
  x
}

# Workbook-writing helper: the biological/statistical results are already
# prepared before this function is called. The code below adds sheets and
# presentation formatting without refitting models or changing hit criteria.
amg_workbook <- function(tables, legend, output_file, overwrite = TRUE) {
  amg_require("openxlsx")
  if (!length(tables) || anyDuplicated(names(tables)) || any(nchar(names(tables)) > 31)) 
    stop("Invalid worksheet names.")
  if (dir.exists(output_file)) 
    stop("Workbook output path is a directory: ", output_file)
  if (file.exists(output_file) && !overwrite) 
    stop("Output exists: ", output_file)
  # A. Create the workbook and its Legend sheet. Wrapping, column widths and
  # row heights make long method descriptions readable in Excel.
  wb <- openxlsx::createWorkbook(creator = "CRISPRi manuscript")
  openxlsx::modifyBaseFont(wb, fontName = "Arial", fontSize = 10, fontColour = "#111111")
  openxlsx::addWorksheet(wb, "Legend", gridLines = FALSE)
  openxlsx::writeData(wb, "Legend", legend, colNames = FALSE, rowNames = FALSE, keepNA = FALSE)
  merged_title <- nrow(legend) > 0L && nzchar(legend$field[1L]) && !nzchar(legend$description[1L])
  if (merged_title) 
    openxlsx::mergeCells(wb, "Legend", cols = 1:2, rows = 1)
  wrap <- openxlsx::createStyle(wrapText = TRUE, valign = "top")
  openxlsx::addStyle(wb, "Legend", wrap, rows = seq_len(nrow(legend)), cols = 1:2, gridExpand = TRUE)
  openxlsx::setColWidths(wb, "Legend", cols = 1:2, widths = c(48, 112))
  openxlsx::addStyle(wb, "Legend", openxlsx::createStyle(fontSize = 13, textDecoration = "bold", 
    valign = "top", wrapText = TRUE), rows = 1, cols = 1:2, gridExpand = TRUE, stack = TRUE)
  sections <- which(legend$description == "" & nzchar(legend$field) | legend$field == "Worksheets")
  if (length(sections)) 
    openxlsx::addStyle(wb, "Legend", openxlsx::createStyle(textDecoration = "bold", valign = "top", 
      wrapText = TRUE), rows = sections, cols = 1:2, gridExpand = TRUE, stack = TRUE)
  heights <- pmax(18, 14 * pmax(ceiling(nchar(legend$field)/43), ceiling(nchar(legend$description)/105)) + 
    6)
  if (merged_title) 
    heights[1L] <- max(35, 18 * ceiling(nchar(legend$field[1L])/115) + 8)
  openxlsx::setRowHeights(wb, "Legend", rows = seq_len(nrow(legend)), heights = heights)
  openxlsx::freezePane(wb, "Legend", firstActiveRow = 3)
  header <- openxlsx::createStyle(textDecoration = "bold", wrapText = TRUE, valign = "center", 
    border = "bottom", borderColour = "#D9D9D9")
  # B. Write one worksheet for each named results table, in the supplied order.
  # Missing values appear blank. Infinite numerical values are written explicitly
  # as Inf/-Inf text rather than being silently lost during Excel export.
  for (name in names(tables)) {
    x <- as.data.frame(tables[[name]])
    infinite <- lapply(x, function(column) if (is.numeric(column)) 
      which(is.infinite(column))
    else integer())
    stored <- x
    for (column in seq_along(stored)) if (length(infinite[[column]])) 
      stored[[column]][infinite[[column]]] <- NA_real_
    openxlsx::addWorksheet(wb, name, gridLines = FALSE)
    openxlsx::writeData(wb, name, stored, headerStyle = header, keepNA = FALSE)
    for (column in seq_along(infinite)) for (row in infinite[[column]]) openxlsx::writeData(wb, 
      name, if (x[[column]][row] > 0) 
        "Inf"
      else "-Inf", startRow = row + 1L, startCol = column, colNames = FALSE)
    openxlsx::setRowHeights(wb, name, rows = 1, heights = 45)
    # C. Set display widths and frozen panes so identifiers remain visible while
    # scrolling. These are workbook layout choices, not additional data filters.
    widths <- rep(19, ncol(x))
    widths[names(x) %in% c("PAO1_ID", "PA14_ID", "gene_name")] <- c(18, 17, 18)[match(names(x)[names(x) %in% 
      c("PAO1_ID", "PA14_ID", "gene_name")], c("PAO1_ID", "PA14_ID", "gene_name"))]
    widths[names(x) %in% c("category", "design_call", "selection_direction")] <- 26
    widths[names(x) %in% c("product_name", "TU", "pathway_term")] <- 42
    widths[grepl("status", names(x))] <- 34
    widths[names(x) == "leadingEdge"] <- 80
    widths[names(x) == "condition"] <- 36
    widths[names(x) == "display_label"] <- 28
    openxlsx::setColWidths(wb, name, cols = seq_len(ncol(x)), widths = widths)
    openxlsx::freezePane(wb, name, firstActiveRow = 2, firstActiveCol = if (name == "Gene_numbers") 
      1
    else 3)
    if (nrow(x)) {
      openxlsx::addFilter(wb, name, rows = 1, cols = seq_len(ncol(x)))
      # D. Choose Excel number formats from the column names/types: scientific
      # notation for p-value-like fields, whole numbers for counts and decimals for
      # other numerical values. Display precision does not round the stored finite
      # values or repeat any multiple-testing correction.
      probability <- which(grepl("p_stouffer|fdr|p[ _.-]*value|p_gsea|^q_primary$|_score$", 
        names(x), ignore.case = TRUE))
      probability <- probability[!grepl("log10", names(x)[probability], ignore.case = TRUE)]
      count <- which(grepl("^n_|_n$|^num_guides$|_rank$|_n_conditions_|(^| )Count$|^Gene-set size$", 
        names(x), ignore.case = TRUE) & vapply(x, is.numeric, logical(1)))
      real <- which(vapply(x, is.numeric, logical(1)))
      effect <- setdiff(real, union(probability, count))
      if (length(probability)) 
        openxlsx::addStyle(wb, name, openxlsx::createStyle(numFmt = "0.000E+00"), rows = 2:(nrow(x) + 
         1), cols = probability, gridExpand = TRUE, stack = TRUE)
      if (length(effect)) 
        openxlsx::addStyle(wb, name, openxlsx::createStyle(numFmt = "0.000"), rows = 2:(nrow(x) + 
         1), cols = effect, gridExpand = TRUE, stack = TRUE)
      if (length(count)) 
        openxlsx::addStyle(wb, name, openxlsx::createStyle(numFmt = "0"), rows = 2:(nrow(x) + 
         1), cols = count, gridExpand = TRUE, stack = TRUE)
    }
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  # E. Save a temporary workbook first, then replace the requested output.
  # Protect the previous workbook during replacement and restore it if promotion
  # fails. OVERWRITE controls whether an existing completed workbook may be replaced.
  staged <- tempfile(".amg_table_", tmpdir = dirname(output_file), fileext = ".xlsx")
  on.exit(unlink(staged))
  openxlsx::saveWorkbook(wb, staged, overwrite = TRUE)
  backup <- tempfile(".previous_amg_table_", tmpdir = dirname(output_file), fileext = ".xlsx")
  protected <- FALSE
  completed <- FALSE
  on.exit({
    if (protected && !completed && !file.rename(backup, output_file)) warning("The previous workbook remains at ", 
      backup)
    if (completed) unlink(backup)
  }, add = TRUE)
  if (file.exists(output_file)) {
    if (!overwrite) 
      stop("Output exists: ", output_file)
    if (!file.rename(output_file, backup)) 
      stop("Could not protect the existing workbook.")
    protected <- TRUE
  }
  if (!file.rename(staged, output_file)) 
    stop("Could not save ", output_file)
  completed <- TRUE
  message("Saved ", output_file)
  invisible(output_file)
}

main <- function() {
  amg_require(c("data.table", "openxlsx"))
  # 1. Read four distinct prepared views of the MH/FBS comparison.
  # The Stouffer and MAGeCK sheets report direct medium-difference hits from their
  # respective methods. The vulnerability sheet compares each medium versus T0.
  # The KEGG sheet contains the full prepared pathway summary, not just the
  # subset selected for the Figure 2H plot.
  paths <- c(Stouffer_FBS_vs_MH_significant = file.path(GENE_INPUT_DIR, "Stouffer_FBS_vs_MH_significant.tsv"), 
    MAGeCK_FBS_vs_MH_significant = file.path(MAGECK_INPUT_DIR, "MAGeCK_FBS_vs_MH_significant.tsv"), 
    Sev_Sens_EK50_MH_FBS_vs_T0 = file.path(GENE_INPUT_DIR, "Sev_Sens_EK50_MH_FBS_vs_T0.tsv"), 
    KEGG_MH_vs_FBS_Fig2H = file.path(KEGG_INPUT_DIR, "KEGG_MH_vs_FBS_Fig2H.tsv"))
  tables <- lapply(paths, amg_read)
  # 2. Check the manuscript-facing column names and order for each sheet.
  # Keep the two methods and the within-medium vulnerability measures separate;
  # their p-values, effect estimates and interpretation are not interchangeable.
  expected <- list(Stouffer_FBS_vs_MH_significant = c("PAO1_ID", "gene_name", "direct_FBS_vs_MH_medLFC", 
    "direct_FBS_vs_MH_p_stouffer", "MH_medLFC_vs_T0", "MH_p_stouffer_vs_T0", "FBS_medLFC_vs_T0", 
    "FBS_p_stouffer_vs_T0"), MAGeCK_FBS_vs_MH_significant = c("PAO1_ID", "gene_name", "gene_lfc", 
    "selection_direction", "neg_fdr", "pos_fdr", "neg_pvalue", "pos_pvalue", "neg_score", "pos_score", 
    "neg_rank", "pos_rank", "num_guides"), Sev_Sens_EK50_MH_FBS_vs_T0 = c("PAO1_ID", "gene_name", 
    "MH_medLFC_vs_T0", "MH_p_stouffer_vs_T0", "FBS_medLFC_vs_T0", "FBS_p_stouffer_vs_T0", "MH_sensitivity", 
    "MH_sensitivity_status", "FBS_sensitivity", "FBS_sensitivity_status", "sensitivity_status_pair", 
    "MH_EK50", "MH_EK50_status", "FBS_EK50", "FBS_EK50_status", "EK50_status_pair"), KEGG_MH_vs_FBS_Fig2H = c("pathway_term", 
    "n_genes_in_pathway", "n_genes_with_MH_FBS_fitness", "median_MH_medLFC", "median_FBS_medLFC", 
    "pathway_delta_FBS_minus_MH", "median_direct_FBS_vs_MH_medLFC", "n_genes_direct_significant", 
    "p_gsea", "FDR_gsea", "ES", "NES_gsea", "leadingEdge"))
  # Require unique locus IDs in gene sheets and unique pathway names in the
  # pathway sheet. This catches duplicate rows before writing the workbook.
  for (name in names(tables)) {
    if (!identical(colnames(tables[[name]]), expected[[name]])) 
      stop("Unexpected columns/order in ", name, ". Re-run its preprocessing script.")
    amg_unique(tables[[name]], if (name == "KEGG_MH_vs_FBS_Fig2H") 
      "pathway_term"
    else "PAO1_ID", name)
  }
  # 3. Fill legend counts from the actual rows and import the Stouffer, MAGeCK
  # and GSEA parameters from their respective preprocessing settings.
  # These describe completed analyses; no new cutoffs are applied here.
  tokens <- list(n_stouffer = nrow(tables$Stouffer_FBS_vs_MH_significant), n_mageck = nrow(tables$MAGeCK_FBS_vs_MH_significant), 
    n_absolute = nrow(tables$Sev_Sens_EK50_MH_FBS_vs_T0), n_pathways = nrow(tables$KEGG_MH_vs_FBS_Fig2H), 
    p_stouffer = amg_setting(GENE_INPUT_DIR, "P_STOUFFER_CUT"), lfc_cut = amg_setting(GENE_INPUT_DIR, 
      "LFC_CUT"), mageck_fdr = amg_setting(MAGECK_INPUT_DIR, "MAGECK_FDR_CUT"), gsea_min = amg_setting(KEGG_INPUT_DIR, 
      "GSEA_MIN_SIZE"), gsea_max = amg_setting(KEGG_INPUT_DIR, "GSEA_MAX_SIZE"), gsea_nperm = amg_setting(KEGG_INPUT_DIR, 
      "GSEA_NPERM"))
  legend <- amg_legend(LEGEND_FILE, tokens)
  # 4. Assemble the four results sheets plus the editable Legend sheet.
  # The formatting helper preserves the prepared biological results.
  amg_workbook(tables, legend, OUTPUT_FILE, OVERWRITE)
}

main()

