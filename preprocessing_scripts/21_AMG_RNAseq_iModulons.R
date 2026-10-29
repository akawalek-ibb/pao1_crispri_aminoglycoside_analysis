# Prepare complete iModulon results and the Figure 6C module selection.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/21_AMG_RNAseq_iModulons.R")
# Inputs:
# - RNAseq_data/iModulon_analysis/DiMA_primary_all.tsv: 116 modules × five contrasts, including
#   nonsignificant rows.
# - annotation_input/precise411-im_table.csv: matching module annotations.
# Upstream / execution order:
# - Use the supplied DiMA table, or RNAseq_analysis/AMG_RNAseq_featureCounts_edgeR.R →
#   RNAseq_data/iModulon_analysis/preprocesing/imodulon_activity_pipeline.py.
# - Python output location is set in config.json; the current config writes to recomputed/, not
#   the frozen input above.
# Outputs:
# - preprocessing_scripts/out_21_AMG_RNAseq_iModulons/ — iModulon_all.tsv,
#   iModulon_selection.tsv and Analysis_settings.tsv.
# Interpretation / operation:
# - No projection or model is refitted. Preserve global BH q-values; display flags require
#   coverage OK, q < 0.05 and |delta| > 5.

INPUT_FILE <- "RNAseq_data/iModulon_analysis/DiMA_primary_all.tsv"

ATLAS_FILE <- "annotation_input/precise411-im_table.csv"

OUTPUT_DIR <- "preprocessing_scripts/out_21_AMG_RNAseq_iModulons"

Q_CUT <- 0.05

MIN_ABS_DELTA <- 5

EXPECTED_MODULES <- 116L

CHECK_GLOBAL_BH <- TRUE

OVERWRITE <- TRUE

CONTRASTS <- data.frame(contrast_id = c("atpIg1m2_vs_control1", "accCg1m3_vs_control1", "secAg1m3_vs_control1", 
  "ftsZg1m4_vs_control1", "control1_TOB_vs_control1"), comparison_id = c("atpI", "accC", "secA", "ftsZ", 
  "TOB"), contrast_label = c("atpI silencing", "accC silencing", "secA silencing", "ftsZ silencing", 
  "TOB treatment"), stringsAsFactors = FALSE)

CONTRAST_ALIASES <- c(atpIg1m2_vs_control1 = "atpIg1m2_vs_control1", atpI_vs_control1 = "atpIg1m2_vs_control1", 
  accC_vs_control1 = "accCg1m3_vs_control1", secA_vs_control1 = "secAg1m3_vs_control1", ftsZ_vs_control1 = "ftsZg1m4_vs_control1")

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

amg_numeric <- function(x, columns, label) {
  for (column in columns) {
    original <- x[[column]]
    value <- suppressWarnings(as.numeric(original))
    if (any(!is.na(original) & (!is.finite(value)))) 
      stop(label, ": invalid numeric value in ", column, call. = FALSE)
    data.table::set(x, j = column, value = value)
  }
  invisible(x)
}

amg_settings <- function(settings, inputs) {
  data.table::data.table(setting = c(names(settings), paste0("input:", names(inputs)), paste0("md5:", 
    names(inputs)), "R_version"), value = c(as.character(settings), unname(inputs), unname(tools::md5sum(inputs)), 
    R.version.string))
}

rnaseq_clean <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

rnaseq_column <- function(x, aliases, label, required = TRUE) {
  found <- aliases[aliases %in% names(x)]
  if (!length(found)) {
    if (required) 
      stop(label, " lacks a column matching: ", paste(aliases, collapse = ", "))
    return(NA_character_)
  }
  found[1L]
}

rnaseq_read <- function(path) {
  if (!file.exists(path)) 
    stop("Missing input: ", path)
  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    amg_require("openxlsx")
    return(data.table::as.data.table(openxlsx::read.xlsx(path, sheet = 1, check.names = FALSE)))
  }
  if (grepl("\\.csv$", path, ignore.case = TRUE)) 
    return(data.table::as.data.table(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, 
      na.strings = c("", "NA"))))
  amg_read(path)
}

rnaseq_probability <- function(x, columns, label) {
  amg_numeric(x, columns, label)
  for (column in columns) if (any(x[[column]] < 0 | x[[column]] > 1, na.rm = TRUE)) 
    stop(label, ": probability outside [0,1] in ", column)
  invisible(x)
}

rnaseq_truth <- function(x, label) {
  value <- tolower(trimws(as.character(x)))
  if (anyNA(value) || any(!value %in% c("true", "false", "t", "f", "1", "0"))) 
    stop(label, ": invalid/missing logical value.")
  value %in% c("true", "t", "1")
}

rnaseq_write <- function(tables, directory, overwrite = TRUE, gzip = character()) {
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

main <- function() {
  amg_require("data.table")
  if (!is.finite(Q_CUT) || Q_CUT <= 0 || Q_CUT > 1 || !is.finite(MIN_ABS_DELTA) || MIN_ABS_DELTA < 
    0) 
    stop("Invalid iModulon cutoffs.")
  # 1. Import the complete upstream differential iModulon activity (DiMA) results.
  # This script does not project expression onto the atlas or re-estimate module
  # activities; it checks and formats the supplied activity differences.
  x <- rnaseq_read(INPUT_FILE)
  required <- c("contrast_id", "module_id", "delta_activity", "ci95_low", "ci95_high", "p_value", 
    "q_primary", "q_global", "projection_coverage_ok", "n_reference_members")
  amg_required(x, required, INPUT_FILE)
  x[, `:=`(module_id = rnaseq_clean(module_id), contrast_id = rnaseq_clean(contrast_id))]
  x[contrast_id %in% names(CONTRAST_ALIASES), `:=`(contrast_id, unname(CONTRAST_ALIASES[contrast_id]))]
  amg_unique(x, c("module_id", "contrast_id"), INPUT_FILE)
  # 2. Require all planned module-by-contrast tests, including nonsignificant
  # ones. A prefiltered significant-only table cannot reproduce the global
  # multiple-testing correction or the display-selection rule.
  if (!setequal(unique(x$contrast_id), CONTRASTS$contrast_id)) 
    stop("Expected exactly the five planned contrasts in the complete DiMA file.")
  if (data.table::uniqueN(x$module_id) != EXPECTED_MODULES || any(x[, .N, by = module_id]$N != 
    nrow(CONTRASTS))) 
    stop("Incomplete module-by-contrast grid; supply the full unfiltered DiMA table.")
  amg_numeric(x, c("delta_activity", "ci95_low", "ci95_high", "n_reference_members"), INPUT_FILE)
  rnaseq_probability(x, c("p_value", "q_primary", "q_global"), INPUT_FILE)
  if (any(!is.finite(x$delta_activity)) || any(x$ci95_low > x$ci95_high, na.rm = TRUE) || any(x$delta_activity < 
    x$ci95_low | x$delta_activity > x$ci95_high, na.rm = TRUE)) 
    stop("Invalid activity estimates/intervals.")
  if (any(!is.finite(x$n_reference_members) | x$n_reference_members < 0 | x$n_reference_members != 
    floor(x$n_reference_members))) 
    stop("Invalid reference-member counts.")
  if (any(x[, data.table::uniqueN(n_reference_members), by = module_id]$V1 != 1L)) 
    stop("Reference-member counts disagree among contrasts.")
  # 3. Verify that q_primary represents BH over the complete planned test family,
  # not a separate correction within each contrast. Optionally recalculate BH
  # for checking, while preserving the upstream p-values and activity estimates.
  if (!isTRUE(all.equal(x$q_primary, x$q_global, tolerance = 1e-10, check.attributes = FALSE))) 
    stop("S9 requires global BH q_primary values, not per-contrast correction.")
  if (CHECK_GLOBAL_BH) {
    bh <- rep(NA_real_, nrow(x))
    finite <- is.finite(x$p_value)
    bh[finite] <- stats::p.adjust(x$p_value[finite], method = "BH")
    if (any(is.finite(bh) != is.finite(x$q_primary)) || any(abs(bh[finite] - x$q_primary[finite]) > 
      1e-10)) 
      stop("Stored q_primary does not match global BH over the complete planned test family.")
  }
  # 4. Call a module/contrast DiMA result only when projection coverage passes,
  # global q is below the cutoff and the absolute activity difference exceeds
  # the configured minimum. Keep all other test rows in the exported table.
  coverage <- rnaseq_truth(x$projection_coverage_ok, "projection_coverage_ok")
  x[, `:=`(coverage_ok = coverage, DiMA_calculated = coverage & is.finite(q_primary) & q_primary < 
    Q_CUT & abs(delta_activity) > MIN_ABS_DELTA)]
  x[, `:=`(comparison_id = CONTRASTS$comparison_id[match(contrast_id, CONTRASTS$contrast_id)], 
    contrast_label = CONTRASTS$contrast_label[match(contrast_id, CONTRASTS$contrast_id)], plot_order = match(contrast_id, 
      CONTRASTS$contrast_id))]
  # 5. Attach atlas categories/descriptions from the supplied reference, or
  # use embedded annotations if no external metadata file exists. Check that
  # embedded and external labels agree when both are available.
  annotation_source <- INPUT_FILE
  if (file.exists(ATLAS_FILE)) {
    atlas <- rnaseq_read(ATLAS_FILE)
    candidates <- names(atlas)[vapply(atlas, function(v) all(unique(x$module_id) %in% trimws(as.character(v))), 
      logical(1))]
    if (!length(candidates)) 
      stop("Atlas metadata does not match all module IDs.")
    id <- candidates[1L]
    cat <- rnaseq_column(atlas, c("category", "functional_category"), ATLAS_FILE)
    fun <- rnaseq_column(atlas, c("function_description", "function", "description"), ATLAS_FILE)
    metadata <- unique(data.table::data.table(module_id = rnaseq_clean(atlas[[id]]), category = rnaseq_clean(atlas[[cat]]), 
      function_description = rnaseq_clean(atlas[[fun]])))
    annotation_source <- ATLAS_FILE
  } else {
    amg_required(x, c("category", "function_description"), "Embedded DiMA annotations")
    metadata <- unique(x[, .(module_id, category = rnaseq_clean(category), function_description = rnaseq_clean(function_description))])
  }
  amg_unique(metadata, "module_id", "PRECISE411 metadata")
  if (any(!unique(x$module_id) %in% metadata$module_id) || anyNA(metadata$category)) 
    stop("Missing official module annotations.")
  for (column in intersect(c("category", "function_description"), names(x))) {
    original <- rnaseq_clean(x[[column]])
    official <- metadata[[column]][match(x$module_id, metadata$module_id)]
    if (any(!is.na(original) & !is.na(official) & original != official)) 
      stop("Embedded and external PRECISE411 annotations disagree in ", column, "; use the metadata from the same reference release.")
  }
  drop <- intersect(c("category", "function_description"), names(x))
  if (length(drop)) 
    x[, `:=`((drop), NULL)]
  x <- merge(x, metadata, by = "module_id", all.x = TRUE, sort = FALSE)
  # 6. Select a module for Figure 6C when it passes DiMA criteria in at least
  # one contrast. Keep its other contrasts too, allowing the figure to show
  # the full pattern instead of only significant cells.
  selection <- x[, .(category = category[1L], function_description = function_description[1L], 
    n_DiMA_contrasts = sum(DiMA_calculated), max_abs_delta = max(abs(delta_activity)), minimum_q = if (any(is.finite(q_primary))) 
      min(q_primary[is.finite(q_primary)])
    else NA_real_, selected_for_Fig6C = any(DiMA_calculated)), by = module_id]
  x[, `:=`(selected_for_Fig6C, selection$selected_for_Fig6C[match(module_id, selection$module_id)])]
  data.table::setorder(x, module_id, plot_order)
  result <- x[, .(module_id, contrast_label, category, function_description, n_reference_members, 
    delta_activity, ci95_low, ci95_high, p_value, q_primary, DiMA = DiMA_calculated, selected_for_Fig6C)]
  data.table::setorder(selection, category, -max_abs_delta, module_id)
  inputs <- c(DiMA = INPUT_FILE, annotation = annotation_source)
  settings <- amg_settings(c(Q_CUT = Q_CUT, MIN_ABS_DELTA = MIN_ABS_DELTA, n_modules = data.table::uniqueN(x$module_id), 
    n_contrasts = nrow(CONTRASTS), n_tests = nrow(x), n_selected_modules = sum(selection$selected_for_Fig6C), 
    CHECK_GLOBAL_BH = CHECK_GLOBAL_BH, estimation = "upstream projection/model results preserved", 
    selection = "coverage OK AND q_primary < cutoff AND abs(delta) > cutoff in at least one comparison"), 
    inputs)
  # 7. Export all module/contrast results, module-selection metadata and cutoffs
  # for Table S9 and Figure 6C. Upstream activity estimates are unchanged.
  rnaseq_write(list(iModulon_all = result, iModulon_selection = selection, Analysis_settings = settings), 
    OUTPUT_DIR, OVERWRITE)
  message("iModulons: ", nrow(result), " results; ", sum(selection$selected_for_Fig6C), " modules selected for Fig. 6C.")
  invisible(result)
}

main()

