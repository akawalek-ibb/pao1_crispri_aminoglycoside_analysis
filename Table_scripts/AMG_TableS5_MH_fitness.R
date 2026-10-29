# Assemble Supplementary Table S5 from completed analysis exports.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Table_scripts/AMG_TableS5_MH_fitness.R")
# Inputs:
# - preprocessing_scripts/out_12_AMG_MH_fitness/MH_fitness.tsv, Essentiality_membership.tsv,
#   Gene_numbers.tsv, MH_vulnerability.tsv and settings.
# - Table_scripts/TableS5_legend.tsv: A header row with field and description columns;
#   {{tokens}} are filled from analysis settings.
# Upstream / execution order:
# - preprocessing_scripts/02_AMG_metacounter.R (only if fitted exports must be rebuilt) →
#   preprocessing_scripts/11_AMG_gene_vulnerability.R →
#   preprocessing_scripts/12_AMG_MH_fitness.R
# Outputs:
# - Table_scripts/out_AMG_TableS5_MH_fitness/Supplementary_Table_S5.xlsx
# Interpretation / operation:
# - Requires openxlsx and data.table.
# - This script formats the workbook; it does not refit models or redefine statistical calls.

INPUT_DIR <- "preprocessing_scripts/out_12_AMG_MH_fitness"

LEGEND_FILE <- "Table_scripts/TableS5_legend.tsv"

OUTPUT_FILE <- "Table_scripts/out_AMG_TableS5_MH_fitness/Supplementary_Table_S5.xlsx"

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
  # 1. Read four completed task-12 exports, one per results worksheet.
  # MH_fitness contains individual untreated-comparison calls and their union;
  # Essentiality_membership keeps published sets separate; Gene_numbers counts
  # results by experiment; MH_vulnerability contains the task-11 T1 consensus.
  names <- c("MH_fitness", "Essentiality_membership", "Gene_numbers", "MH_vulnerability")
  tables <- stats::setNames(lapply(names, function(name) amg_read(file.path(INPUT_DIR, paste0(name, 
    ".tsv")))), names)
  # 2. Define the exact column names/order expected for the manuscript workbook.
  # This is a compatibility check, not a calculation of fitness or essentiality.
  expected <- list(MH_fitness = c("PAO1_ID", "PA14_ID", "gene_name", "product_name", "category", 
    "design_call", "TU", "CRISPRi_negative_fitness_union", "CRISPRi_positive_fitness_union", 
    "n_negative_comparisons", "n_positive_comparisons", as.vector(t(outer(c("B1_T1_vs_T0", 
      "B1_T2_vs_T0", "B2_T1_vs_T0", "B3_T1_vs_T0", "B4_T1_vs_T0"), c("_medLFC", "_p_stouffer", 
      "_call"), paste0)))), Essentiality_membership = c("PAO1_ID", "PA14_ID", "gene_name", 
    "Poulsen_core", "Zhang_2025", "This_study", "Kaczmarczyk_2025"), Gene_numbers = c("batch", 
    "condition", "display_label", "n_genes_tested", "n_negative_fitness", "n_positive_fitness", 
    "n_significant_total"), MH_vulnerability = c("PAO1_ID", "PA14_ID", "gene_name", "control_median_severity", 
    "control_median_sensitivity", "control_median_EK50", "n_control_severity_estimates", "n_control_sensitivity_estimates", 
    "n_control_EK50_estimates"))
  # Require the prepared tables to match that schema and require one row per
  # locus in the gene-level sheets. Gene_numbers instead has one row per comparison.
  for (name in names) if (!identical(colnames(tables[[name]]), expected[[name]])) 
    stop("Unexpected columns/order in ", name, ". Re-run preprocessing 12.")
  for (name in setdiff(names, "Gene_numbers")) amg_unique(tables[[name]], "PAO1_ID", name)
  # 3. Obtain legend counts from the actual tables and significance cutoffs from
  # task-12 settings. The union counts and vulnerability estimates are imported;
  # they are not recomputed from individual guides during workbook assembly.
  tokens <- list(n_MH_fitness = nrow(tables$MH_fitness), n_membership = nrow(tables$Essentiality_membership), 
    n_MH_vulnerability = nrow(tables$MH_vulnerability), n_negative_union = sum(tables$MH_fitness$CRISPRi_negative_fitness_union), 
    n_positive_union = sum(tables$MH_fitness$CRISPRi_positive_fitness_union), p_stouffer = amg_setting(INPUT_DIR, 
      "P_STOUFFER_CUT"), lfc_cut = amg_setting(INPUT_DIR, "LFC_CUT"))
  legend <- amg_legend(LEGEND_FILE, tokens)
  # 4. Pass the prepared sheets and resolved legend to the formatting helper.
  # No guide filtering, curve fitting or significance testing is performed here.
  amg_workbook(tables, legend, OUTPUT_FILE, OVERWRITE)
}

main()

