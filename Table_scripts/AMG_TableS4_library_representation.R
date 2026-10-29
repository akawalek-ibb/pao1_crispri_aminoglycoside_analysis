# Assemble Supplementary Table S4 from completed analysis exports.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Table_scripts/AMG_TableS4_library_representation.R")
# Inputs:
# - preprocessing_scripts/out_10_AMG_library_representation/Sample_key.tsv, Guide_coverage.tsv,
#   Gene_coverage.tsv, Design_summary.tsv, Library_summary.tsv and Analysis_settings.tsv.
# - Table_scripts/TableS4_legend.tsv: Two tab-separated columns without a header.
# Upstream / execution order:
# - preprocessing_scripts/09_AMG_export_GEO_raw_counts.R →
#   preprocessing_scripts/10_AMG_library_representation.R
# Outputs:
# - Table_scripts/out_AMG_TableS4_library_representation/Supplementary_Table_S4.xlsx
# Interpretation / operation:
# - Requires openxlsx.
# - This script formats the workbook; it does not refit models or redefine statistical calls.

INPUT_DIR <- "preprocessing_scripts/out_10_AMG_library_representation"

LEGEND_FILE <- "Table_scripts/TableS4_legend.tsv"

OUTPUT_FILE <- "Table_scripts/out_AMG_TableS4_library_representation/Supplementary_Table_S4.xlsx"

INCLUDE_RAW_COUNTS <- FALSE

FONT_NAME <- "Calibri"

FONT_SIZE <- 11

OVERWRITE <- TRUE

# Local functions for this analysis.
main <- function() {
  if (!requireNamespace("openxlsx", quietly = TRUE)) 
    stop("Required package: openxlsx")
  if (!file.exists(LEGEND_FILE)) 
    stop("Missing editable legend: ", LEGEND_FILE)
  if (file.exists(OUTPUT_FILE) && !OVERWRITE) 
    stop("Output exists: ", OUTPUT_FILE)
  read_input <- function(name) {
    path <- file.path(INPUT_DIR, paste0(name, ".tsv"))
    if (!file.exists(path)) 
      stop("Run 10_AMG_library_representation.R first. Missing: ", path)
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "", 
      na.strings = "NA")
  }
  # 1. Read the five prepared library-representation tables from task 10.
  # They describe samples, individual guides, target genes, design categories and
  # overall guide classes. Detection/cloning criteria were applied upstream;
  # this script does not recount reads or redefine those criteria.
  sheet_names <- c("Sample_key", "Guide_coverage", "Gene_coverage", "Design_summary", "Library_summary")
  tables <- stats::setNames(lapply(sheet_names, read_input), sheet_names)
  # 2. Read the thresholds and count-source information used to prepare those
  # tables. Use these recorded settings to fill the editable legend accurately.
  settings <- read_input("Analysis_settings")
  if (!identical(names(settings), c("parameter", "value")) || anyDuplicated(settings$parameter)) 
    stop("Invalid Analysis_settings.tsv.")
  settings <- stats::setNames(as.character(settings$value), settings$parameter)
  if (anyNA(settings[c("detection_min_reads", "cloned_mean_T0_min_reads", "count_source", "count_source_md5")])) 
    stop("Required analysis settings are missing.")
  # 3. Load the two-column, headerless legend file and substitute its placeholders.
  # Add/remove the Raw_counts legend entry according to INCLUDE_RAW_COUNTS.
  # An unresolved placeholder stops export rather than producing incomplete text.
  legend <- utils::read.delim(LEGEND_FILE, header = FALSE, sep = "\t", quote = "", comment.char = "", 
    check.names = FALSE, stringsAsFactors = FALSE, na.strings = "")
  if (ncol(legend) != 2L) 
    stop("Legend must have exactly two tab-separated columns.")
  tokens <- c(`{{n_guides}}` = format(nrow(tables$Guide_coverage), big.mark = ",", scientific = FALSE, 
    trim = TRUE), `{{detection_reads}}` = settings[["detection_min_reads"]], `{{cloning_reads}}` = settings[["cloned_mean_T0_min_reads"]])
  for (token in names(tokens)) legend[[2]] <- gsub(token, tokens[[token]], legend[[2]], fixed = TRUE)
  if (!INCLUDE_RAW_COUNTS) 
    legend <- legend[is.na(legend[[1]]) | legend[[1]] != "Raw_counts", , drop = FALSE]
  if (INCLUDE_RAW_COUNTS && !any(legend[[1]] == "Raw_counts", na.rm = TRUE)) 
    legend <- rbind(legend, c("Raw_counts", "Unfiltered guide counts from the generated count matrix, with deposited sample names and NA entries retained."))
  if (any(grepl("{{", legend[[2]], fixed = TRUE), na.rm = TRUE)) 
    stop("Unresolved placeholder in the editable legend.")
  # 4. Optionally include the full unfiltered guide-count matrix as another sheet.
  # Check its checksum against task 10, then verify guide order and sample headers.
  # This option is off by default and does not change the five coverage summaries.
  if (INCLUDE_RAW_COUNTS) {
    count_file <- settings[["count_source"]]
    if (!file.exists(count_file)) 
      stop("Missing GEO raw-count file: ", count_file)
    if (!identical(unname(tools::md5sum(count_file)), settings[["count_source_md5"]])) 
      stop("GEO count file changed; rerun preprocessing before generating the table.")
    connection <- if (endsWith(count_file, ".gz")) 
      gzfile(count_file, "rt")
    else file(count_file, "rt")
    tables$Raw_counts <- tryCatch(utils::read.delim(connection, check.names = FALSE, stringsAsFactors = FALSE, 
      quote = "", comment.char = "", na.strings = "NA"), finally = close(connection))
    if (!identical(tables$Raw_counts$spacer, tables$Guide_coverage$spacer) || !identical(names(tables$Raw_counts)[-(1:3)], 
      tables$Sample_key$label)) 
      stop("Raw counts and preprocessing TSVs do not match.")
  }
  # 5. Create the Legend sheet and apply text wrapping/row heights.
  # This affects readability only; thresholds and counts remain those imported.
  wb <- openxlsx::createWorkbook()
  openxlsx::modifyBaseFont(wb, fontSize = FONT_SIZE, fontName = FONT_NAME)
  openxlsx::addWorksheet(wb, "Legend")
  openxlsx::writeData(wb, "Legend", legend, colNames = FALSE, keepNA = FALSE)
  openxlsx::setColWidths(wb, "Legend", cols = 1:2, widths = c(31, 110))
  openxlsx::addStyle(wb, "Legend", openxlsx::createStyle(wrapText = TRUE, valign = "top"), rows = seq_len(nrow(legend)), 
    cols = 1:2, gridExpand = TRUE)
  legend_heights <- pmax(24, 15 * ceiling(pmax(nchar(ifelse(is.na(legend[[1]]), "", legend[[1]]))/29, 
    nchar(ifelse(is.na(legend[[2]]), "", legend[[2]]))/100)))
  openxlsx::setRowHeights(wb, "Legend", rows = seq_len(nrow(legend)), heights = legend_heights)
  legend_sections <- which(seq_len(nrow(legend)) == 1L | grepl("^(Worksheets|Methods and interpretation|Column definitions:)", 
    legend[[1]]))
  openxlsx::addStyle(wb, "Legend", openxlsx::createStyle(textDecoration = "bold"), rows = legend_sections, 
    cols = 1:2, gridExpand = TRUE, stack = TRUE)
  header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#E7E6E6", wrapText = TRUE, 
    valign = "top")
  # 6. Write each prepared table to its own sheet. Keep explicit NA entries in
  # Raw_counts and Guide_coverage; other missing cells are left blank. Add filters
  # and frozen panes so the large guide tables remain practical to inspect.
  for (sheet in names(tables)) {
    data <- tables[[sheet]]
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, data, headerStyle = header_style, withFilter = TRUE, keepNA = sheet %in% 
      c("Raw_counts", "Guide_coverage"), na.string = "NA")
    openxlsx::freezePane(wb, sheet, firstRow = TRUE, firstCol = sheet %in% c("Raw_counts", 
      "Guide_coverage", "Gene_coverage"))
    widths <- pmin(32, pmax(14, nchar(names(data)) + 2))
    widths[names(data) %in% c("spacer", "label", "original_sample_id")] <- 28
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = widths)
    openxlsx::setRowHeights(wb, sheet, rows = 1, heights = 34)
  }
  # 7. Display mean T0 counts to one decimal place without changing their stored
  # values, then save the completed workbook at OUTPUT_FILE.
  mean_cols <- grep("_mean_no_drug_T0$", names(tables$Guide_coverage))
  openxlsx::addStyle(wb, "Guide_coverage", openxlsx::createStyle(numFmt = "0.0"), rows = seq_len(nrow(tables$Guide_coverage)) + 
    1L, cols = mean_cols, gridExpand = TRUE, stack = TRUE)
  dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, OUTPUT_FILE, overwrite = OVERWRITE)
  message("Saved Supplementary Table S4: ", OUTPUT_FILE)
  invisible(OUTPUT_FILE)
}

main()

