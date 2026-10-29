# Assemble Supplementary Table S8 from completed analysis exports.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Table_scripts/AMG_TableS8_aminoglycoside_interactions.R")
# Inputs:
# - preprocessing_scripts/out_03_AMG_relevant_results/
# - preprocessing_scripts/out_04_AMG_compare_MAGeCK_Stouffer/
# - preprocessing_scripts/out_05_AMG_TU_interactions/
# - preprocessing_scripts/out_06_AMG_TU_functional_annotation/
# - preprocessing_scripts/out_07_AMG_recurrent_gene_TU_collapse/
# - preprocessing_scripts/out_08_AMG_recurrent_KEGG_enrichment/
# - preprocessing_scripts/out_16_AMG_Murray2015_comparison/
# - preprocessing_scripts/out_17_AMG_mismatch_interaction_candidates/
# - Table_scripts/TableS8_legend.tsv: A header row with field and description columns;
#   {{tokens}} are filled from analysis settings.
# Upstream / execution order:
# - preprocessing_scripts/03_AMG_relevant_results.R; then
#   preprocessing_scripts/04_AMG_compare_MAGeCK_Stouffer.R and
#   preprocessing_scripts/05_AMG_TU_interactions.R.
# - preprocessing_scripts/06_AMG_TU_functional_annotation.R →
#   preprocessing_scripts/07_AMG_recurrent_gene_TU_collapse.R and
#   preprocessing_scripts/08_AMG_recurrent_KEGG_enrichment.R.
# - preprocessing_scripts/16_AMG_Murray2015_comparison.R and
#   preprocessing_scripts/17_AMG_mismatch_interaction_candidates.R both require completed script
#   03.
# Outputs:
# - Table_scripts/out_AMG_TableS8_aminoglycoside_interactions/Supplementary_Table_S8.xlsx
# Interpretation / operation:
# - Requires openxlsx and data.table.
# - This script formats the workbook; it does not refit models or redefine statistical calls.

GENE_DIR <- "preprocessing_scripts/out_03_AMG_relevant_results"

COMPARISON_DIR <- "preprocessing_scripts/out_04_AMG_compare_MAGeCK_Stouffer"

TU_TEST_DIR <- "preprocessing_scripts/out_05_AMG_TU_interactions"

TU_ANNOTATION_DIR <- "preprocessing_scripts/out_06_AMG_TU_functional_annotation"

TU_COLLAPSE_DIR <- "preprocessing_scripts/out_07_AMG_recurrent_gene_TU_collapse"

KEGG_DIR <- "preprocessing_scripts/out_08_AMG_recurrent_KEGG_enrichment"

MURRAY_DIR <- "preprocessing_scripts/out_16_AMG_Murray2015_comparison"

MISMATCH_DIR <- "preprocessing_scripts/out_17_AMG_mismatch_interaction_candidates"

LEGEND_FILE <- "Table_scripts/TableS8_legend.tsv"

OUTPUT_FILE <- "Table_scripts/out_AMG_TableS8_aminoglycoside_interactions/Supplementary_Table_S8.xlsx"

CHECK_INTERMEDIATE_CONSISTENCY <- TRUE

OVERWRITE <- TRUE

# The schema below specifies the required columns/order of each results sheet.
# These checks ensure the assembled workbook matches its prepared inputs;
# the long lists do not implement additional statistical tests.
expected <- list(counts_thresholds = c("comparison_id", "display_label", "n_stouffer_evaluated", 
  "n_stouffer_sensitizing", "n_stouffer_protective", "n_mageck_evaluated", "n_mageck_sensitizing", 
  "n_mageck_protective", "n_shared_sensitizing", "n_shared_protective", "n_stouffer_only_sensitizing", 
  "n_stouffer_only_protective", "n_mageck_only_sensitizing", "n_mageck_only_protective", "n_union_called", 
  "p_stouffer_cutoff", "median_medLFC", "MAD", "sensitizing_cutoff_source", "sensitizing_LFC_cutoff_used", 
  "protective_LFC_cutoff_used"), Stouffer_calls = c("PAO1_ID", "gene_name", "product_name", "essential_category_code", 
  "TU", "n_independent_T1_sensitizing", "n_independent_T1_protective", "recurrent_sensitizing", 
  "recurrent_protective", "TOB_0.025_medLFC", "TOB_0.1_medLFC", "TOB_0.4_medLFC", "TOB_1_medLFC", 
  "TOB_T2_medLFC", "AMK_medLFC", "TOB_0.025_p_stouffer", "TOB_0.1_p_stouffer", "TOB_0.4_p_stouffer", 
  "TOB_1_p_stouffer", "TOB_T2_p_stouffer", "AMK_p_stouffer", "TOB_0.025_call", "TOB_0.1_call", 
  "TOB_0.4_call", "TOB_1_call", "TOB_T2_call", "AMK_call", "TOB_0.025_p_stouffer_BH", "TOB_0.1_p_stouffer_BH", 
  "TOB_0.4_p_stouffer_BH", "TOB_1_p_stouffer_BH", "AMK_p_stouffer_BH"), MAGeCK_calls = c("comparison_id", 
  "PAO1_ID", "gene_name", "mageck_gene_lfc", "mageck_neg_fdr", "mageck_pos_fdr", "mageck_call", 
  "medLFC_stouffer", "p_stouffer", "stouffer_call", "stouffer_evaluated", "mageck_evaluated", 
  "method_comparison", "shared_same_direction", "any_call"), MAGeCK_vs_Stouffer = c("PAO1_ID", 
  "gene_name", "n_stouffer_T1_sensitizing", "n_stouffer_T1_protective", "n_mageck_T1_sensitizing", 
  "n_mageck_T1_protective", "n_shared_T1_sensitizing", "n_shared_T1_protective", "stouffer_recurrent_T1_sensitizing", 
  "stouffer_recurrent_T1_protective", "mageck_recurrent_T1_sensitizing", "mageck_recurrent_T1_protective", 
  "shared_recurrent_T1_sensitizing", "shared_recurrent_T1_protective"), Mismatch_based_evidence = c("comparison_id", 
  "PAO1_ID", "gene_name", "product_name", "direction", "interaction_medLFC", "conventional_gene_level_p_stouffer", 
  "interaction_threshold (sensitivity)", "partial_bin_extreme", "EC50", "Emax"), Gene_level_hits_by_TU = c("TU_ID", 
  "n_all_genes_in_TU", "n_recurrent_genes", "recurrent_gene_PAO1_IDs", "recurrent_gene_names", 
  "n_recurrent_sensitizing_genes", "recurrent_sensitizing_gene_PAO1_IDs", "recurrent_sensitizing_gene_names", 
  "n_recurrent_protective_genes", "recurrent_protective_gene_PAO1_IDs", "recurrent_protective_gene_names", 
  "TU_contains_recurrent_sensitizing", "TU_contains_recurrent_protective", "TU_recurrent_direction", 
  "protective_conditions_union", "sensitizing_conditions_union"), TU_level_interactions = c("TU_ID", 
  "batch", "condition", "tu_medLFC", "tu_p_stouffer", "n_perfect_guides", "n_targeted_genes", 
  "targeted_PAO1_IDs", "assigned_effect_TU"), KEGG_enrichment = c("term", "direction", "gene_overlap_n", 
  "gene_total_term_n", "gene_coverage_pct", "gene_fold_enrichment", "gene_odds_ratio", "gene_p_value", 
  "gene_FDR_BH", "gene_overlap_members", "TU_overlap_n", "TU_total_term_n", "TU_coverage_pct", 
  "TU_fold_enrichment", "TU_odds_ratio", "TU_p_value", "TU_FDR_BH", "TU_overlap_members"), Murray_comparison = c("PAO1_ID", 
  "gene_name", "recurrent_direction", "Murray_n_conditions_tested", "Murray_same_direction_conditions", 
  "Murray_consensus_same_direction_min2", "Murray_any_opposite_direction", "Murray_opposite_direction_conditions", 
  "Murray_GEN_PA14_IDs", "Murray_GEN_gene", "Murray_GEN_logFC", "Murray_NEO_logFC", "Murray_TOB_logFC", 
  "Murray_GEN_FDR", "Murray_NEO_FDR", "Murray_TOB_FDR", "Murray_GEN_direction", "Murray_NEO_direction", 
  "Murray_TOB_direction"))

# Verify the recorded checksums of upstream inputs under preprocessing_scripts/.
# A missing/changed dependency stops assembly, preventing accidental mixing of
# results from different runs. Some annotation files are needed for this check
# even though they do not become worksheets themselves.
check_intermediate_inputs <- function(directory) {
  settings <- s8_read(file.path(directory, "Analysis_settings.tsv"))
  inputs <- settings[grepl("^input:", setting)]
  for (i in seq_len(nrow(inputs))) {
    path <- as.character(inputs$value[i])
    if (!startsWith(path, "preprocessing_scripts/")) 
      next
    key <- sub("^input:", "md5:", inputs$setting[i])
    previous <- settings$value[settings$setting == key]
    if (!file.exists(path) || length(previous) != 1L || is.na(previous) || unname(tools::md5sum(path)) != 
      previous) 
      stop("Stale/missing intermediate for ", directory, ": ", path, ". Re-run its preprocessing task.")
  }
}

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

s8_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) 
    stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

s8_required <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) 
    stop(label, " lacks: ", paste(missing, collapse = ", "), call. = FALSE)
}

s8_read <- function(path, columns = NULL) {
  if (!file.exists(path)) 
    stop("Missing input: ", path, call. = FALSE)
  read_once <- function(nrows = -1L, classes = NA) {
    con <- if (endsWith(path, ".gz")) 
      gzfile(path, "rt")
    else file(path, "rt")
    on.exit(close(con))
    utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = "", 
      na.strings = c("", "NA"), nrows = nrows, colClasses = classes)
  }
  if (is.null(columns)) 
    return(data.table::as.data.table(read_once()))
  header <- names(read_once(0L))
  s8_required(as.data.frame(stats::setNames(rep(list(character()), length(header)), header)), 
    columns, path)
  data.table::as.data.table(read_once(classes = ifelse(header %in% columns, "character", "NULL")))
}

s8_unique <- function(x, keys, label) {
  if (anyNA(x[, keys, with = FALSE]) || anyDuplicated(x[, keys, with = FALSE])) 
    stop(label, ": missing or duplicate key: ", paste(keys, collapse = ", "), call. = FALSE)
  invisible(x)
}

s8_setting <- function(directory, name, numeric = FALSE) {
  x <- s8_read(file.path(directory, "Analysis_settings.tsv"))
  value <- x$value[x$setting == name]
  if (length(value) != 1L || is.na(value)) 
    stop("Missing or ambiguous setting: ", name)
  if (numeric) 
    as.numeric(value)
  else as.character(value)
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
  s8_require(c("data.table", "openxlsx"))
  # 1. Check consistency across the contributing analyses when the option is on.
  # Each producer records its own inputs; this checks those dependencies rather
  # than assuming that files found in the expected folders belong to the same run.
  directories <- c(GENE_DIR, COMPARISON_DIR, TU_TEST_DIR, TU_ANNOTATION_DIR, TU_COLLAPSE_DIR, 
    KEGG_DIR, MURRAY_DIR, MISMATCH_DIR)
  if (CHECK_INTERMEDIATE_CONSISTENCY) 
    for (directory in directories) check_intermediate_inputs(directory)
  # 2. Read the nine prepared results sheets. Their roles differ:
  # Stouffer_calls holds primary gene hits; MAGeCK_calls compares methods within
  # each gene/comparison, while MAGeCK_vs_Stouffer summarizes T1 recurrence.
  # Gene_level_hits_by_TU groups recurrent genes; TU_level_interactions tests
  # pooled guides within TUs. The remaining sheets report thresholds, exploratory
  # mismatch evidence, KEGG enrichment and the published Murray comparison.
  paths <- c(counts_thresholds = file.path(COMPARISON_DIR, "counts_thresholds.tsv"), Stouffer_calls = file.path(GENE_DIR, 
    "Stouffer_calls.tsv"), MAGeCK_calls = file.path(COMPARISON_DIR, "MAGeCK_calls.tsv"), MAGeCK_vs_Stouffer = file.path(COMPARISON_DIR, 
    "MAGeCK_vs_Stouffer.tsv"), Mismatch_based_evidence = file.path(MISMATCH_DIR, "Mismatch_based_evidence.tsv"), 
    Gene_level_hits_by_TU = file.path(TU_COLLAPSE_DIR, "Gene_level_hits_by_TU.tsv"), TU_level_interactions = file.path(TU_TEST_DIR, 
      "TU_level_interactions.tsv.gz"), KEGG_enrichment = file.path(KEGG_DIR, "KEGG_enrichment.tsv"), 
    Murray_comparison = file.path(MURRAY_DIR, "Murray_comparison.tsv"))
  tables <- lapply(paths, s8_read)
  # 3. Require the expected columns/order, then check row identities at the
  # appropriate level: gene, gene/comparison, gene/comparison/direction, TU,
  # TU/comparison, pathway/direction or published gene/direction.
  # These checks preserve the distinct units of the analyses.
  for (name in names(tables)) {
    if (!identical(names(tables[[name]]), expected[[name]])) 
      stop("Unexpected named columns/order in ", name, "; re-run its producer.")
  }
  s8_unique(tables$Stouffer_calls, "PAO1_ID", "Stouffer_calls")
  s8_unique(tables$MAGeCK_calls, c("PAO1_ID", "comparison_id"), "MAGeCK_calls")
  s8_unique(tables$MAGeCK_vs_Stouffer, "PAO1_ID", "MAGeCK_vs_Stouffer")
  s8_unique(tables$Mismatch_based_evidence, c("PAO1_ID", "comparison_id", "direction"), "Mismatch_based_evidence")
  s8_unique(tables$Gene_level_hits_by_TU, "TU_ID", "Gene_level_hits_by_TU")
  s8_unique(tables$TU_level_interactions, c("TU_ID", "batch", "condition"), "TU_level_interactions")
  s8_unique(tables$KEGG_enrichment, c("term", "direction"), "KEGG_enrichment")
  s8_unique(tables$Murray_comparison, c("PAO1_ID", "recurrent_direction"), "Murray_comparison")
  # 4. Read supporting counts/settings for the legend: recurrent genes, BH
  # support and the gene/TU enrichment query/background sizes. These support
  # the methods description but are not added as extra results worksheets.
  manifest <- s8_read(file.path(GENE_DIR, "Comparison_manifest.tsv"))
  recurrent <- s8_read(file.path(GENE_DIR, "AG_gene_recurrence.tsv"))
  bh <- s8_read(file.path(GENE_DIR, "BH_summary.tsv"))
  kegg <- s8_read(file.path(KEGG_DIR, "KEGG_settings.tsv"))
  sens <- kegg[direction == "Sensitizing"]
  prot <- kegg[direction == "Protective"]
  if (nrow(sens) != 1L || nrow(prot) != 1L) 
    stop("Invalid KEGG query-size summary.")
  # Use one legend value when sensitizing/protective profile parameters agree;
  # otherwise report both explicitly so directional differences remain visible.
  paired_setting <- function(a, b) {
    x <- s8_setting(MISMATCH_DIR, a)
    y <- s8_setting(MISMATCH_DIR, b)
    if (x == y) 
      x
    else paste0(x, " (sensitizing), ", y, " (protective)")
  }
  # 5. Populate the editable legend from the completed analyses, including their
  # actual recurrence, FDR, MAD, bin/LOESS and external-comparison parameters.
  # This does not redefine primary hits using the additional supporting analyses.
  tokens <- list(n_sensitizing = sum(recurrent$recurrent_sensitizing), n_protective = sum(recurrent$recurrent_protective), 
    n_TU_sensitizing = sum(tables$Gene_level_hits_by_TU$TU_contains_recurrent_sensitizing), 
    n_TU_protective = sum(tables$Gene_level_hits_by_TU$TU_contains_recurrent_protective), n_T1 = sum(manifest$count_for_recurrence), 
    recurrence_conditions = paste(manifest$display_label[manifest$count_for_recurrence], collapse = "; "), 
    p_stouffer = s8_setting(GENE_DIR, "P_STOUFFER_CUT"), lfc_floor = s8_setting(GENE_DIR, "MIN_ABS_MEDLFC"), 
    mad_multiplier = s8_setting(GENE_DIR, "MAD_MULTIPLIER"), bh_cut = s8_setting(GENE_DIR, 
      "P_BH_CUT"), recurrence_min = s8_setting(GENE_DIR, "RECURRENCE_MIN"), mageck_fdr = s8_setting(COMPARISON_DIR, 
      "MAGECK_FDR_CUT"), BH_sensitizing = bh$BH_supported_recurrent[bh$direction == "Sensitizing"], 
    BH_protective = bh$BH_supported_recurrent[bh$direction == "Protective"], gene_query_sensitizing = sens$gene_query_n, 
    gene_query_protective = prot$gene_query_n, TU_query_sensitizing = sens$TU_query_n, TU_query_protective = prot$TU_query_n, 
    gene_universe = sens$gene_universe_n, TU_universe = sens$TU_universe_n, n_KEGG_terms = sens$n_terms_tested, 
    min_guides = s8_setting(MISMATCH_DIR, "MIN_TOTAL_GUIDES"), min_bins = s8_setting(MISMATCH_DIR, 
      "MIN_QUALIFIED_BINS"), min_bin_guides = s8_setting(MISMATCH_DIR, "MIN_BIN_N"), partial_max = s8_setting(MISMATCH_DIR, 
      "PARTIAL_REPRESSION_MAX"), loess_min_bins = s8_setting(MISMATCH_DIR, "LOESS_MIN_BINS"), 
    loess_auc = paired_setting("LOESS_AUC_MIN_SENS", "LOESS_AUC_MIN_PROT"), loess_width = paired_setting("LOESS_MIN_WIDTH_SENS", 
      "LOESS_MIN_WIDTH_PROT"), loess_min_x = s8_setting(MISMATCH_DIR, "LOESS_AUC_X_MIN"), 
    loess_max_x = s8_setting(MISMATCH_DIR, "LOESS_AUC_X_MAX"), murray_fdr = s8_setting(MURRAY_DIR, 
      "MURRAY_FDR_CUT"), murray_lfc = s8_setting(MURRAY_DIR, "MURRAY_ABS_LFC_CUT"))
  # 6. Arrange selected worksheets for inspection, keeping recurrent/shared
  # hits prominent. Sorting changes row order only, not membership or estimates.
  data.table::setorder(tables$Stouffer_calls, -recurrent_sensitizing, -recurrent_protective, 
    PAO1_ID)
  data.table::setorder(tables$MAGeCK_calls, comparison_id, method_comparison, PAO1_ID)
  data.table::setorder(tables$MAGeCK_vs_Stouffer, -shared_recurrent_T1_sensitizing, -shared_recurrent_T1_protective, 
    PAO1_ID)
  data.table::setorder(tables$KEGG_enrichment, -direction, gene_p_value, TU_p_value, term)
  data.table::setorderv(tables$Murray_comparison, c("Murray_same_direction_conditions", "recurrent_direction", 
    "Murray_any_opposite_direction", "PAO1_ID"), na.last = TRUE)
  legend <- amg_legend(LEGEND_FILE, tokens)
  # 7. Write the resolved legend and the nine results sheets using the formatting
  # helper. No count models, enrichment tests or response curves are fitted here.
  amg_workbook(tables, legend, OUTPUT_FILE, OVERWRITE)
}

main()

