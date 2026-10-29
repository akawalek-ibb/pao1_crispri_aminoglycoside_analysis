# Test interactions after pooling full-match guides by transcriptional unit (TU).
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/05_AMG_TU_interactions.R")
# Inputs:
# - preprocessing_scripts/out_03_AMG_relevant_results/b*_AG_guides.tsv.gz,
#   Comparison_manifest.tsv and Analysis_settings.tsv.
# - annotation_input/PAO1_curated.tsv: PAO1_ID-to-TU mapping.
# Upstream / execution order:
# - preprocessing_scripts/03_AMG_relevant_results.R with EXPORT_GUIDE_TABLES <- TRUE; its fitted
#   screen inputs come from preprocessing_scripts/02_AMG_metacounter.R.
# Outputs:
# - preprocessing_scripts/out_05_AMG_TU_interactions/ — TU_level_interactions.tsv.gz,
#   TU_thresholds.tsv and Analysis_settings.tsv.
# Interpretation / operation:
# - This pooled-guide statistical test differs from grouping recurrent gene hits by TU.
# - TU Stouffer combination clamps guide p-values to [1e-15, 1−1e-15].

INPUT_DIR <- "preprocessing_scripts/out_03_AMG_relevant_results"

ANNOTATION_FILE <- "annotation_input/PAO1_curated.tsv"

OUTPUT_DIR <- "preprocessing_scripts/out_05_AMG_TU_interactions"

ALLOW_GUIDE_TU_FALLBACK <- TRUE

TU_P_EPS <- 1e-15

OVERWRITE <- TRUE

# Combine guide p-values directly within a transcriptional unit (TU).
# Clamp extreme values to avoid infinite normal scores; effect direction is
# assigned separately from the median adjusted guide effect.
s8_TU_stouffer <- function(p) {
  p <- p[is.finite(p)]
  if (!length(p)) 
    return(NA_real_)
  p <- pmin(pmax(p, TU_P_EPS), 1 - TU_P_EPS)
  stats::pnorm(sum(stats::qnorm(1 - p))/sqrt(length(p)), lower.tail = FALSE)
}

# Local functions for this analysis.
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

s8_clean <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | toupper(x) %in% c("", "NA", "NAN", "NONE", "NULL")] <- NA_character_
  x
}

s8_unique <- function(x, keys, label) {
  if (anyNA(x[, keys, with = FALSE]) || anyDuplicated(x[, keys, with = FALSE])) 
    stop(label, ": missing or duplicate key: ", paste(keys, collapse = ", "), call. = FALSE)
  invisible(x)
}

s8_annotation <- function(path) {
  x <- s8_read(path)
  s8_required(x, c("PAO1_ID", "gene_name", "product_name", "category", "design_call", "TU"), 
    path)
  x[, `:=`(PAO1_ID = s8_clean(PAO1_ID), TU = s8_clean(TU))]
  s8_unique(x, "PAO1_ID", path)
  x
}

s8_collapse <- function(x, sep = "; ") {
  x <- sort(unique(s8_clean(x)))
  x <- x[!is.na(x)]
  if (length(x)) 
    paste(x, collapse = sep)
  else NA_character_
}

s8_threshold <- function(lfc, p, p_cut, multiplier, floor) {
  finite <- is.finite(lfc) & is.finite(p)
  if (!any(finite)) 
    stop("A required comparison has no finite evaluated results.")
  center <- stats::median(lfc[finite])
  spread <- stats::mad(lfc[finite], center = center, constant = 1)
  lower <- if (is.finite(spread) && spread > 0) 
    center - multiplier * spread
  else NA_real_
  upper <- if (is.finite(spread) && spread > 0) 
    center + multiplier * spread
  else NA_real_
  data.table::data.table(n_stouffer_evaluated = sum(finite), p_stouffer_cutoff = p_cut, median_medLFC = center, 
    MAD = spread, sensitizing_cutoff_source = if (is.finite(lower) && lower <= -floor) 
      paste0("median - ", multiplier, "MAD")
    else paste0("fixed -", floor, " / ", floor), sensitizing_LFC_cutoff_used = if (is.finite(lower)) 
      min(-floor, lower)
    else -floor, protective_LFC_cutoff_used = if (is.finite(upper)) 
      max(floor, upper)
    else floor)
}

s8_call <- function(lfc, p, p_cut, lower, upper) {
  data.table::fcase(is.finite(p) & p < p_cut & is.finite(lfc) & lfc <= lower, "Sensitizing", 
    is.finite(p) & p < p_cut & is.finite(lfc) & lfc >= upper, "Protective", default = "Not significant")
}

s8_settings <- function(settings, inputs = character()) {
  data.table::data.table(setting = c(names(settings), paste0("input:", names(inputs)), paste0("md5:", 
    names(inputs)), "R_version"), value = c(as.character(settings), unname(inputs), unname(tools::md5sum(inputs)), 
    R.version.string))
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

s8_write <- function(tables, directory, overwrite = TRUE, gzip = character()) {
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
  s8_require("data.table")
  if (s8_setting(INPUT_DIR, "EXPORT_GUIDE_TABLES") != "TRUE") 
    stop("Run 03 with EXPORT_GUIDE_TABLES <- TRUE before TU aggregation.")
  # 1. Restrict analysis to the five independent T1 drug comparisons.
  # Read the screening cutoffs from task 03 so TU and gene calls use the same
  # parameter values, while estimating their effect distributions separately.
  manifest <- s8_read(file.path(INPUT_DIR, "Comparison_manifest.tsv"))[count_for_recurrence == 
    TRUE]
  p_cut <- s8_setting(INPUT_DIR, "P_STOUFFER_CUT", numeric = TRUE)
  multiplier <- s8_setting(INPUT_DIR, "MAD_MULTIPLIER", numeric = TRUE)
  floor <- s8_setting(INPUT_DIR, "MIN_ABS_MEDLFC", numeric = TRUE)
  if (!is.finite(TU_P_EPS) || TU_P_EPS <= 0 || TU_P_EPS >= 0.5) 
    stop("Invalid TU guide-p clamp.")
  annotation <- s8_annotation(ANNOTATION_FILE)
  paths <- file.path(INPUT_DIR, paste0(unique(manifest$batch), "_AG_guides.tsv.gz"))
  # 2. Read guide effects and map each targeted gene to its annotated TU.
  # The optional fallback uses the TU stored with the guide when the curated
  # assignment is missing. Retain full-match guides with valid, nonzero raw LFC.
  parts <- lapply(seq_along(paths), function(i) {
    x <- s8_read(paths[i])
    s8_required(x, c("PAO1_ID", "gene_name", "spacer", "type", "TU", "LFC", "LFC.adj", "pval", 
      "batch", "condition", "comparison_id"), paths[i])
    x <- x[comparison_id %in% manifest$comparison_id]
    s8_unique(x, c("batch", "condition", "spacer"), paths[i])
    x[, `:=`(PAO1_ID = s8_clean(PAO1_ID), TU_guide = s8_clean(TU))]
    idx <- match(x$PAO1_ID, annotation$PAO1_ID)
    x[, `:=`(TU_ID, annotation$TU[idx])]
    if (ALLOW_GUIDE_TU_FALLBACK) 
      x[is.na(TU_ID), `:=`(TU_ID, TU_guide)]
    x[is.na(gene_name), `:=`(gene_name, annotation$gene_name[match(PAO1_ID, annotation$PAO1_ID)])]
    x <- x[type == "perfect" & !is.na(PAO1_ID) & !is.na(TU_ID) & is.finite(LFC) & LFC != 0 & 
      is.finite(LFC.adj) & is.finite(pval)]
    x
  })
  guides <- data.table::rbindlist(parts, use.names = TRUE)
  if (!nrow(guides)) 
    stop("No eligible full-match guides for TU analysis.")
  # 3. Pool eligible individual guides across genes belonging to the same TU.
  # Calculate a TU median effect and combine guide p-values directly; this is
  # not an average of gene-level p-values or a count of significant member genes.
  tus <- guides[, .(tu_medLFC = stats::median(LFC.adj), tu_p_stouffer = s8_TU_stouffer(pval), 
    n_perfect_guides = data.table::uniqueN(spacer), n_targeted_genes = data.table::uniqueN(PAO1_ID), 
    targeted_PAO1_IDs = s8_collapse(PAO1_ID)), by = .(TU_ID, batch, condition)]
  # 4. Recalculate effect cutoffs from the TU-level distribution in each
  # comparison, then assign sensitizing/protective TU calls.
  thresholds <- tus[, s8_threshold(tu_medLFC, tu_p_stouffer, p_cut, multiplier, floor), by = .(batch, 
    condition)]
  tus <- merge(tus, thresholds[, .(batch, condition, sensitizing_LFC_cutoff_used, protective_LFC_cutoff_used)], 
    by = c("batch", "condition"), all.x = TRUE, sort = FALSE)
  tus[, `:=`(assigned_effect_TU, s8_call(tu_medLFC, tu_p_stouffer, p_cut, sensitizing_LFC_cutoff_used, 
    protective_LFC_cutoff_used))]
  data.table::setorder(tus, batch, condition, tu_p_stouffer, TU_ID)
  columns <- c("TU_ID", "batch", "condition", "tu_medLFC", "tu_p_stouffer", "n_perfect_guides", 
    "n_targeted_genes", "targeted_PAO1_IDs", "assigned_effect_TU")
  settings <- s8_settings(c(P_STOUFFER_CUT = p_cut, MAD_MULTIPLIER = multiplier, MIN_ABS_MEDLFC = floor, 
    TU_P_EPS = TU_P_EPS, ALLOW_GUIDE_TU_FALLBACK = ALLOW_GUIDE_TU_FALLBACK, contrast_scope = "five T1 AG comparisons; no T2"), 
    c(stats::setNames(paths, unique(manifest$batch)), annotation = ANNOTATION_FILE, call_settings = file.path(INPUT_DIR, 
      "Analysis_settings.tsv"), manifest = file.path(INPUT_DIR, "Comparison_manifest.tsv")))
  # 5. Export TU effects and calls, the cutoffs actually used, and provenance.
  # Targeted-gene and guide counts show how much evidence supports each TU.
  s8_write(list(TU_level_interactions = tus[, ..columns], TU_thresholds = thresholds, Analysis_settings = settings), 
    OUTPUT_DIR, OVERWRITE, gzip = "TU_level_interactions")
  message("S8 pooled-guide TU analysis: ", nrow(tus), " TU/comparison records.")
  invisible(tus)
}

main()

