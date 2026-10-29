# Assemble Supplementary Table S9 from completed analysis exports.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Table_scripts/AMG_TableS9_RNAseq.R")
# Inputs:
# - preprocessing_scripts/out_18_AMG_RNAseq_results/
# - preprocessing_scripts/out_19_AMG_RNAseq_GO_ORA/
# - preprocessing_scripts/out_20_AMG_RNAseq_GO_GSEA/
# - preprocessing_scripts/out_21_AMG_RNAseq_iModulons/
# - preprocessing_scripts/out_22_AMG_RNAseq_selected_GO/
# - Table_scripts/TableS9_legend.tsv: A header row with field and description columns;
#   {{tokens}} are filled from analysis settings.
# Upstream / execution order:
# - preprocessing_scripts/18_AMG_RNAseq_results.R → preprocessing_scripts/19_AMG_RNAseq_GO_ORA.R
#   and preprocessing_scripts/20_AMG_RNAseq_GO_GSEA.R.
# - preprocessing_scripts/22_AMG_RNAseq_selected_GO.R follows 18/19;
#   preprocessing_scripts/21_AMG_RNAseq_iModulons.R independently imports the supplied DiMA
#   results.
# Outputs:
# - Table_scripts/out_AMG_TableS9_RNAseq/Supplementary_Table_S9.xlsx
# Interpretation / operation:
# - Requires openxlsx and data.table.
# - This script formats the workbook; it does not refit models or redefine statistical calls.

RNA_DIR <- "preprocessing_scripts/out_18_AMG_RNAseq_results"

ORA_DIR <- "preprocessing_scripts/out_19_AMG_RNAseq_GO_ORA"

GSEA_DIR <- "preprocessing_scripts/out_20_AMG_RNAseq_GO_GSEA"

IMODULON_DIR <- "preprocessing_scripts/out_21_AMG_RNAseq_iModulons"

SELECTED_GO_DIR <- "preprocessing_scripts/out_22_AMG_RNAseq_selected_GO"

LEGEND_FILE <- "Table_scripts/TableS9_legend.tsv"

OUTPUT_FILE <- "Table_scripts/out_AMG_TableS9_RNAseq/Supplementary_Table_S9.xlsx"

CHECK_INTERMEDIATE_CONSISTENCY <- TRUE

OVERWRITE <- TRUE

# The schema below fixes worksheet columns/order for the manuscript table.
# RNAseq_all includes all five expression contrasts, including TOB; GO sheets
# and module/display results retain the comparison scopes set by their producers.
expected <- list(RNAseq_all = c("PAO1 locus tag", "Gene name", "atpI log2FC", "atpI FDR", "atpI P value", 
  "atpI DEG direction", "accC log2FC", "accC FDR", "accC P value", "accC DEG direction", "secA log2FC", 
  "secA FDR", "secA P value", "secA DEG direction", "ftsZ log2FC", "ftsZ FDR", "ftsZ P value", 
  "ftsZ DEG direction", "TOB log2FC", "TOB FDR", "TOB P value", "TOB DEG direction"), iModulon_all = c("module_id", 
  "contrast_label", "category", "function_description", "n_reference_members", "delta_activity", 
  "ci95_low", "ci95_high", "p_value", "q_primary", "DiMA", "selected_for_Fig6C"), GO_ORA_all = c("Comparison", 
  "DEG direction", "GO domain", "GO ID", "GO term", "FDR (BH)", "P value", "Count", "BgRatio", 
  "FoldEnrichment", "Overlapping gene IDs", "Overlapping gene names"), GO_GSEA_all = c("Comparison", 
  "GO domain", "GO ID", "GO term", "NES", "FDR (BH)", "P value", "Gene-set size", "Leading-edge gene IDs", 
  "Leading-edge gene names"), Figure6D = c("Comparison", "GO domain", "GO ID", "GO term", "Displayed direction", 
  "Displayed signed -log10(FDR)", "Up Count", "Up GeneRatio", "Up FDR", "Up result", "Up gene IDs", 
  "Up gene names", "Down Count", "Down GeneRatio", "Down FDR", "Down result", "Down gene IDs", 
  "Down gene names"))

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

rnaseq_setting <- function(directory, name, numeric = FALSE) {
  x <- amg_read(file.path(directory, "Analysis_settings.tsv"))
  value <- x$value[x$setting == name]
  if (length(value) != 1L || is.na(value)) 
    stop("Missing/ambiguous setting: ", name)
  if (numeric) 
    as.numeric(value)
  else as.character(value)
}

# Check recorded dependencies under preprocessing_scripts/out_ against their
# stored checksums. This detects stale/missing prepared results before different
# RNA-seq analysis stages are combined in one workbook.
rnaseq_check_inputs <- function(directory) {
  x <- amg_read(file.path(directory, "Analysis_settings.tsv"))
  inputs <- x[grepl("^input:", setting)]
  for (i in seq_len(nrow(inputs))) {
    path <- as.character(inputs$value[i])
    if (!startsWith(path, "preprocessing_scripts/out_")) 
      next
    hash <- x$value[x$setting == sub("^input:", "md5:", inputs$setting[i])]
    if (!file.exists(path) || length(hash) != 1L || is.na(hash) || unname(tools::md5sum(path)) != 
      hash) 
      stop("Stale/missing intermediate for ", directory, ": ", path, ". Re-run its preprocessing task.")
  }
  invisible(TRUE)
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
  # 1. Verify prepared-input consistency when enabled. Also require GO ORA and
  # GO GSEA to use the same annotation file; the two enrichment methods remain
  # separate analyses, but should refer to the same gene-to-GO catalogue.
  directories <- c(RNA_DIR, ORA_DIR, GSEA_DIR, IMODULON_DIR, SELECTED_GO_DIR)
  if (CHECK_INTERMEDIATE_CONSISTENCY) 
    for (directory in directories) rnaseq_check_inputs(directory)
  if (rnaseq_setting(ORA_DIR, "md5:annotation") != rnaseq_setting(GSEA_DIR, "md5:annotation")) 
    stop("GO ORA/GSEA used different annotation inputs; regenerate with the same annotation.")
  # 2. Read the complete long RNA-seq table and planned comparison order.
  # Keep one row per locus/comparison, using a locus ID as the display name when
  # no gene symbol is available. Do not select only significant genes.
  de <- amg_read(file.path(RNA_DIR, "RNAseq_results.tsv.gz"))
  manifest <- amg_read(file.path(RNA_DIR, "Comparison_manifest.tsv"))
  data.table::setorder(manifest, plot_order)
  amg_unique(de, c("PAO1_ID", "comparison_id"), "RNA-seq results")
  amg_required(de, c("gene_name", "logFC", "FDR", "PValue", "DEG_direction"), "RNA-seq results")
  genes <- de[, .(gene_name = if (any(!is.na(gene_name))) 
    gene_name[which(!is.na(gene_name))[1L]]
  else NA_character_), by = PAO1_ID]
  genes[is.na(gene_name), `:=`(gene_name, PAO1_ID)]
  data.table::setorder(genes, PAO1_ID)
  rna <- data.table::data.table(`PAO1 locus tag` = genes$PAO1_ID, `Gene name` = genes$gene_name)
  # 3. Convert expression results to one row per gene, with separate columns for
  # each comparison's log2FC, FDR, raw p-value and DEG direction. Align by locus ID,
  # not row position. Missing comparison rows receive Not tested as their label;
  # their numerical values remain missing, rather than becoming zero effects.
  for (id in manifest$comparison_id) {
    z <- de[comparison_id == id]
    idx <- match(genes$PAO1_ID, z$PAO1_ID)
    directions <- z$DEG_direction[idx]
    directions[is.na(directions)] <- "Not tested"
    rna[, `:=`((paste(id, c("log2FC", "FDR", "P value", "DEG direction"))), list(z$logFC[idx], 
      z$FDR[idx], z$PValue[idx], directions))]
  }
  # 4. Import all module/contrast results, including nonsignificant tests and the
  # Figure 6C selection flag. Activities, confidence intervals and global q-values
  # come from preprocessing 21; they are not estimated or corrected again here.
  imodulon <- amg_read(file.path(IMODULON_DIR, "iModulon_all.tsv"))
  amg_unique(imodulon, c("module_id", "contrast_label"), "iModulon results")
  # 5. Prepare the GO ORA sheet from Up/Down DEG-list enrichment results.
  # Retain overlap counts, background ratios and gene lists so the enrichment
  # evidence is visible. Renaming/selecting columns does not repeat ORA.
  ora <- amg_read(file.path(ORA_DIR, "GO_ORA_results.tsv.gz"))
  amg_unique(ora, c("comparison_id", "direction", "ontology", "GO_ID"), "GO ORA")
  ora[, `:=`(comparison_order, match(comparison_id, manifest$comparison_id))]
  data.table::setorder(ora, comparison_order, direction, ontology, FDR, P_value, GO_ID)
  ora_table <- ora[, .(Comparison = comparison_id, `DEG direction` = direction, `GO domain` = ontology, 
    `GO ID` = GO_ID, `GO term` = GO_term, `FDR (BH)` = FDR, `P value` = P_value, Count, BgRatio, 
    FoldEnrichment, `Overlapping gene IDs` = overlap_ids, `Overlapping gene names` = overlap_names)]
  # 6. Prepare the separate GO GSEA sheet from ranked-gene enrichment results.
  # Retain NES, corrected/raw p-values, set sizes and leading-edge genes.
  # This tests coordinated ranks rather than overlap with a thresholded DEG list.
  gsea <- amg_read(file.path(GSEA_DIR, "GO_GSEA_results.tsv.gz"))
  amg_unique(gsea, c("comparison_id", "ontology", "GO_ID"), "GO GSEA")
  gsea[, `:=`(comparison_order, match(comparison_id, manifest$comparison_id))]
  data.table::setorder(gsea, comparison_order, ontology, FDR, P_value, GO_ID)
  gsea_table <- gsea[, .(Comparison = comparison_id, `GO domain` = ontology, `GO ID` = GO_ID, 
    `GO term` = GO_term, NES, `FDR (BH)` = FDR, `P value` = P_value, `Gene-set size` = set_size, 
    `Leading-edge gene IDs` = leading_ids, `Leading-edge gene names` = leading_names)]
  # 7. Import the selected GO-term display data used for Figure 6D.
  # Keep the signed/capped display score alongside both Up/Down FDRs, result
  # statuses and overlap genes. The direction choice and display cap were already
  # applied in preprocessing 22, not chosen during workbook assembly.
  selected <- amg_read(file.path(SELECTED_GO_DIR, "selected_GO_results.tsv"))
  amg_unique(selected, c("comparison_id", "GO_ID"), "Selected GO results")
  data.table::setorder(selected, display_order, plot_order)
  selected_table <- selected[, .(Comparison = comparison_id, `GO domain` = ontology, `GO ID` = GO_ID, 
    `GO term` = GO_term, `Displayed direction` = displayed_direction, `Displayed signed -log10(FDR)` = displayed_signed_log10_FDR, 
    `Up Count` = Up_Count, `Up GeneRatio` = Up_GeneRatio, `Up FDR` = Up_FDR, `Up result` = Up_result, 
    `Up gene IDs` = Up_gene_IDs, `Up gene names` = Up_gene_names, `Down Count` = Down_Count, 
    `Down GeneRatio` = Down_GeneRatio, `Down FDR` = Down_FDR, `Down result` = Down_result, 
    `Down gene IDs` = Down_gene_IDs, `Down gene names` = Down_gene_names)]
  # 8. Collect the five data sheets and verify their exact manuscript schemas.
  # No additional DEG, GO or module-significance filters are applied here.
  tables <- list(RNAseq_all = rna, iModulon_all = imodulon, GO_ORA_all = ora_table, GO_GSEA_all = gsea_table, 
    Figure6D = selected_table)
  for (name in names(tables)) if (!identical(names(tables[[name]]), expected[[name]])) 
    stop("Unexpected named columns/order in ", name, "; check its producer.")
  # 9. Summarize ORA backgrounds and Up/Down query sizes for the editable legend.
  # When a background differs among comparisons, report those sizes separately
  # rather than implying a single shared gene universe.
  backgrounds <- amg_read(file.path(ORA_DIR, "GO_ORA_backgrounds.tsv"))
  ontologies <- strsplit(rnaseq_setting(ORA_DIR, "ONTOLOGIES"), ";", fixed = TRUE)[[1L]]
  bg_unique <- unique(backgrounds[, .(ontology, background_n)])
  bg_summary <- if (all(bg_unique[, .N, by = ontology]$N == 1L)) 
    paste(vapply(ontologies, function(o) paste0(o, " ", format(bg_unique$background_n[bg_unique$ontology == 
      o], big.mark = ",", scientific = FALSE, trim = TRUE)), character(1)), collapse = "; ")
  else paste(paste(backgrounds$comparison_id, backgrounds$ontology, backgrounds$background_n), 
    collapse = "; ")
  query_summary <- vapply(manifest$comparison_id[manifest$use_for_GO], function(id) {
    text <- vapply(ontologies, function(o) {
      x <- backgrounds[comparison_id == id & ontology == o]
      paste0(x$query_n[x$direction == "Up"], "/", x$query_n[x$direction == "Down"])
    }, character(1))
    paste0(id, ", ", paste(text, collapse = ", "))
  }, character(1))
  # 10. Retrieve DE, enrichment, module and display cutoffs from the producer
  # settings. Fill legend placeholders with those values and the recorded counts.
  tokens <- list(DE_FDR = rnaseq_setting(RNA_DIR, "DE_FDR_CUT"), DE_FC = rnaseq_setting(RNA_DIR, 
    "DE_ABS_LOG2FC_CUT"), n_modules = rnaseq_setting(IMODULON_DIR, "n_modules"), n_comparisons = rnaseq_setting(IMODULON_DIR, 
    "n_contrasts"), n_tests = rnaseq_setting(IMODULON_DIR, "n_tests"), n_selected_terms = rnaseq_setting(SELECTED_GO_DIR, 
    "n_selected_terms"), ORA_MIN = rnaseq_setting(ORA_DIR, "MIN_GS_SIZE"), ORA_MAX = rnaseq_setting(ORA_DIR, 
    "MAX_GS_SIZE"), GSEA_MIN = rnaseq_setting(GSEA_DIR, "MIN_GS_SIZE"), GSEA_MAX = rnaseq_setting(GSEA_DIR, 
    "MAX_GS_SIZE"), ORA_backgrounds = bg_summary, ORA_ontologies = paste(ontologies, collapse = ", "), 
    ORA_queries = paste(query_summary, collapse = "; "), DiMA_Q = rnaseq_setting(IMODULON_DIR, 
      "Q_CUT"), DiMA_delta = rnaseq_setting(IMODULON_DIR, "MIN_ABS_DELTA"), display_FDR = rnaseq_setting(SELECTED_GO_DIR, 
      "DISPLAY_FDR_CUT"), display_cap = rnaseq_setting(SELECTED_GO_DIR, "DISPLAY_SCORE_CAP"))
  legend <- amg_legend(LEGEND_FILE, tokens)
  # 11. Write the Legend and five data sheets with the formatting helper.
  # Expression testing, enrichment and activity modelling are not rerun.
  amg_workbook(tables, legend, OUTPUT_FILE, OVERWRITE)
  invisible(OUTPUT_FILE)
}

main()

