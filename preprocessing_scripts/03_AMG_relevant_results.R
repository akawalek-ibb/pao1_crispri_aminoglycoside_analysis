# Prepare gene–aminoglycoside interaction calls and reusable guide subsets.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/03_AMG_relevant_results.R")
# Inputs:
# - edgeR_data/b*_median_melted_results.tsv[.gz]: PAO1_ID, condition, medLFC, p_stouffer.
# - edgeR_data/b*_melted_results.tsv[.gz]: guide effects, p-values and relative efficacies.
# - annotation_input/PAO1_curated.tsv: unique PAO1_ID, annotations and TU assignments.
# Upstream / execution order:
# - preprocessing_scripts/02_AMG_metacounter.R if fitted exports are absent; otherwise reuse the
#   supplied exports.
# Outputs:
# - preprocessing_scripts/out_03_AMG_relevant_results/ — AG_gene_results.tsv.gz,
#   AG_gene_recurrence.tsv, Stouffer_calls.tsv, thresholds, manifest, settings and
#   b*_AG_guides.tsv.gz.
# Interpretation / operation:
# - Primary calls: raw p < 0.01 and median ± 5 unscaled MAD, with an absolute LFC floor of 1.
# - Recurrence requires the same direction in ≥2 of five T1 comparisons; T2 is separate.
# - BH support is additional to these primary rules. Guide exports feed scripts 05, 17 and Fig. S13.

INPUT_DIR <- "edgeR_data"

ANNOTATION_FILE <- "annotation_input/PAO1_curated.tsv"

OUTPUT_DIR <- "preprocessing_scripts/out_03_AMG_relevant_results"

P_STOUFFER_CUT <- 0.01

MAD_MULTIPLIER <- 5

MIN_ABS_MEDLFC <- 1

P_BH_CUT <- 0.05

RECURRENCE_MIN <- 2L

EXPORT_GUIDE_TABLES <- TRUE

GUIDE_CHUNK_SIZE <- 100000L

OVERWRITE <- TRUE

COMPARISONS <- utils::read.table(text = "batch|condition|comparison_id|display_label|prefix|plot_order|count_for_recurrence|drug|dose|timing\nb4|TOB_0.025_T1 - None_0_T1|B4_TOB_0.025_T1_vs_None_T1|TOB 0.025|TOB_0.025|1|TRUE|TOB|0.025|T1\nb4|TOB_0.1_T1 - None_0_T1|B4_TOB_0.1_T1_vs_None_T1|TOB 0.1|TOB_0.1|2|TRUE|TOB|0.1|T1\nb1|TOB_0.4_T1 - None_0_T1|B1_TOB_0.4_T1_vs_None_T1|TOB 0.4|TOB_0.4|3|TRUE|TOB|0.4|T1\nb4|TOB_1_T1 - None_0_T1|B4_TOB_1_T1_vs_None_T1|TOB 1|TOB_1|4|TRUE|TOB|1|T1\nb1|TOB_0.4_T2 - None_0_T2|B1_TOB_0.4_T2_vs_None_T2|TOB 0.4 T2|TOB_T2|5|FALSE|TOB|0.4|T2\nb2|AMK_2_T1 - None_0_T1|B2_AMK_2_T1_vs_None_T1|AMK|AMK|6|TRUE|AMK|2|T1", 
  header = TRUE, sep = "|", stringsAsFactors = FALSE)

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

s8_numeric <- function(x, columns, label) {
  for (name in columns) {
    original <- x[[name]]
    value <- suppressWarnings(as.numeric(original))
    if (any(!is.na(original) & !is.finite(value))) 
      stop(label, ": invalid numeric value in ", name)
    data.table::set(x, j = name, value = value)
  }
  invisible(x)
}

s8_probability <- function(x, columns, label) {
  s8_numeric(x, columns, label)
  for (name in columns) if (any(x[[name]] < 0 | x[[name]] > 1, na.rm = TRUE)) 
    stop(label, ": probability outside [0,1] in ", name)
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

s8_file <- function(directory, batch, suffix) {
  paths <- file.path(directory, paste0(batch, suffix, c(".tsv.gz", ".tsv")))
  found <- paths[file.exists(paths)]
  if (!length(found)) 
    stop("Missing input; run 02_AMG_metacounter.R:\n", paste(paths, collapse = "\n"), call. = FALSE)
  found[1L]
}

# Estimate the centre and unscaled MAD of all finite gene effects in a comparison.
# The negative/positive cutoffs are at least as extreme as -1/+1 by default;
# genes also need the raw Stouffer p-value cutoff to receive an interaction call.
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

s8_bh <- function(p) {
  out <- rep(NA_real_, length(p))
  keep <- is.finite(p)
  out[keep] <- stats::p.adjust(p[keep], method = "BH")
  out
}

s8_settings <- function(settings, inputs = character()) {
  data.table::data.table(setting = c(names(settings), paste0("input:", names(inputs)), paste0("md5:", 
    names(inputs)), "R_version"), value = c(as.character(settings), unname(inputs), unname(tools::md5sum(inputs)), 
    R.version.string))
}

s8_stage <- function(tables, directory, gzip = character()) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  for (name in names(tables)) {
    path <- file.path(directory, paste0(name, if (name %in% gzip) 
      ".tsv.gz"
    else ".tsv"))
    con <- if (endsWith(path, ".gz")) 
      gzfile(path, "wt")
    else file(path, "wt")
    tryCatch(utils::write.table(as.data.frame(tables[[name]]), con, sep = "\t", quote = TRUE, 
      qmethod = "double", row.names = FALSE, na = "NA"), finally = close(con))
  }
}

s8_commit <- function(stage, directory, overwrite) {
  files <- list.files(stage, full.names = TRUE)
  files <- files[!dir.exists(files)]
  targets <- file.path(directory, basename(files))
  if (!overwrite && any(file.exists(targets))) 
    stop("Output exists; set OVERWRITE <- TRUE to replace it.")
  if (any(dir.exists(targets))) 
    stop("An output target is a directory.")
  backup <- file.path(stage, "backup")
  dir.create(backup)
  old <- which(file.exists(targets))
  moved_old <- moved_new <- integer()
  completed <- FALSE
  on.exit({
    if (!completed) {
      unlink(targets[moved_new])
      for (i in moved_old) if (!file.rename(file.path(backup, basename(targets[i])), targets[i])) warning("Previous output remains in ", 
        backup)
    }
  }, add = TRUE)
  for (i in old) {
    if (!file.rename(targets[i], file.path(backup, basename(targets[i])))) 
      stop("Could not protect existing output: ", targets[i])
    moved_old <- c(moved_old, i)
  }
  for (i in seq_along(files)) {
    if (!file.rename(files[i], targets[i])) 
      stop("Could not save: ", targets[i])
    moved_new <- c(moved_new, i)
  }
  completed <- TRUE
  unlink(backup, recursive = TRUE)
  invisible(targets)
}

# Extract only the requested drug comparisons from large individual-guide files.
# Read in chunks to limit memory use; keep guide effects and predicted efficacies
# for downstream TU and mismatch-profile analyses. This does not refit edgeR.
s8_stream_guides <- function(path, output, batch, manifest, chunk_size = 100000L) {
  wanted <- c("PAO1_ID", "gene_name", "spacer", "type", "y_pred", "original", "TU", "condition", 
    "LFC", "pval", "LFC.adj")
  con <- if (endsWith(path, ".gz")) 
    gzfile(path, "rt")
  else file(path, "rt")
  on.exit(close(con), add = TRUE)
  header <- names(utils::read.delim(con, nrows = 0, check.names = FALSE, quote = "\"", comment.char = ""))
  if (length(setdiff(wanted, header))) 
    stop(path, " lacks guide columns: ", paste(setdiff(wanted, header), collapse = ", "))
  close(con)
  con <- if (endsWith(path, ".gz")) 
    gzfile(path, "rt")
  else file(path, "rt")
  readLines(con, n = 1L, warn = FALSE)
  out <- gzfile(output, "wt")
  on.exit(close(out), add = TRUE)
  first <- TRUE
  retained <- 0L
  repeat {
    chunk <- tryCatch(utils::read.delim(con, header = FALSE, col.names = header, nrows = chunk_size, 
      check.names = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = "", na.strings = c("", 
        "NA"), colClasses = ifelse(header %in% wanted, "character", "NULL")), error = function(e) {
      if (grepl("no lines available in input", conditionMessage(e), fixed = TRUE)) 
        return(NULL)
      stop(e)
    })
    if (is.null(chunk) || !nrow(chunk)) 
      break
    x <- data.table::as.data.table(chunk[chunk$condition %in% manifest$condition, , drop = FALSE])
    if (!nrow(x)) 
      next
    s8_numeric(x, c("y_pred", "LFC", "pval", "LFC.adj"), path)
    x[type == "perfect", `:=`(y_pred, 1)]
    x[, `:=`(batch, batch)]
    x[, `:=`(comparison_id, manifest$comparison_id[match(condition, manifest$condition)])]
    utils::write.table(as.data.frame(x), out, sep = "\t", quote = TRUE, qmethod = "double", 
      row.names = FALSE, col.names = first, na = "NA")
    first <- FALSE
    retained <- retained + nrow(x)
  }
  if (first) 
    stop("No guide rows matched the required comparisons in ", path)
  data.table::data.table(batch = batch, retained_guide_rows = retained)
}

main <- function() {
  s8_require("data.table")
  if (P_STOUFFER_CUT <= 0 || P_STOUFFER_CUT > 1 || P_BH_CUT <= 0 || P_BH_CUT > 1 || MAD_MULTIPLIER <= 
    0 || MIN_ABS_MEDLFC <= 0 || RECURRENCE_MIN < 2 || RECURRENCE_MIN > sum(COMPARISONS$count_for_recurrence)) 
    stop("Invalid call settings.")
  if (anyDuplicated(COMPARISONS[c("batch", "condition")]) || anyDuplicated(COMPARISONS$comparison_id)) 
    stop("Duplicate comparison manifest entries.")
  # 1. Read the gene annotation and locate complete gene-level edgeR exports.
  # The comparisons selected below are drug versus matched untreated controls;
  # their effects describe drug interactions, rather than fitness versus T0.
  annotation <- s8_annotation(ANNOTATION_FILE)
  files <- stats::setNames(vapply(unique(COMPARISONS$batch), function(b) s8_file(INPUT_DIR, b, 
    "_median_melted_results"), character(1)), unique(COMPARISONS$batch))
  # 2. For each contributing batch, retain the planned comparisons and gene rows.
  # Recalculate BH across finite gene tests in each complete selected comparison,
  # before selecting hits; verify any BH values already stored by metacounter.
  parts <- lapply(names(files), function(b) {
    x <- s8_read(files[[b]])
    s8_required(x, c("PAO1_ID", "condition", "medLFC", "p_stouffer"), files[[b]])
    if ("type" %in% names(x)) 
      x <- x[type == "perfect"]
    x[, `:=`(PAO1_ID, s8_clean(PAO1_ID))]
    x <- x[!is.na(PAO1_ID) & condition %in% COMPARISONS$condition[COMPARISONS$batch == b]]
    s8_unique(x, c("PAO1_ID", "condition"), files[[b]])
    s8_numeric(x, "medLFC", files[[b]])
    s8_probability(x, "p_stouffer", files[[b]])
    x[, `:=`(p_BH_recomputed, s8_bh(p_stouffer)), by = condition]
    if ("p_stouffer_BH" %in% names(x)) {
      s8_probability(x, "p_stouffer_BH", files[[b]])
      finite <- is.finite(x$p_stouffer_BH) & is.finite(x$p_BH_recomputed)
      if (any(is.finite(x$p_stouffer_BH) != is.finite(x$p_BH_recomputed)) || any(abs(x$p_stouffer_BH[finite] - 
        x$p_BH_recomputed[finite]) > 1e-08)) 
        stop("Stored BH values differ from full-comparison recalculation in ", files[[b]])
    }
    x[, `:=`(p_stouffer_BH = p_BH_recomputed, batch = b)]
    x[, `:=`(c("p_BH_recomputed"), NULL)]
    x[, .(PAO1_ID, batch, condition, medLFC, p_stouffer, p_stouffer_BH)]
  })
  genes <- data.table::rbindlist(parts)
  manifest <- data.table::as.data.table(COMPARISONS)
  genes <- merge(genes, manifest, by = c("batch", "condition"), all.x = TRUE, sort = FALSE)
  if (length(setdiff(manifest$comparison_id, genes$comparison_id))) 
    stop("Some required AG comparisons are absent.")
  if (length(setdiff(genes$PAO1_ID, annotation$PAO1_ID))) 
    stop("Gene results contain loci absent from curated annotation.")
  genes <- merge(genes, annotation[, .(PAO1_ID, gene_name, product_name, essential_category_code = category, 
    design_call, TU)], by = "PAO1_ID", all.x = TRUE, sort = FALSE)
  # 3. Calculate separate effect cutoffs for each comparison, then call genes.
  # Negative effects meeting both criteria are sensitizing; positive effects are
  # protective. The BH-supported column adds a filter to these primary calls.
  thresholds <- genes[, s8_threshold(medLFC, p_stouffer, P_STOUFFER_CUT, MAD_MULTIPLIER, MIN_ABS_MEDLFC), 
    by = .(comparison_id, display_label, plot_order, count_for_recurrence)]
  genes <- merge(genes, thresholds[, .(comparison_id, sensitizing_LFC_cutoff_used, protective_LFC_cutoff_used)], 
    by = "comparison_id", all.x = TRUE, sort = FALSE)
  genes[, `:=`(stouffer_evaluated = is.finite(medLFC) & is.finite(p_stouffer), stouffer_call = s8_call(medLFC, 
    p_stouffer, P_STOUFFER_CUT, sensitizing_LFC_cutoff_used, protective_LFC_cutoff_used))]
  genes[, `:=`(stouffer_call_BH, ifelse(stouffer_call != "Not significant" & is.finite(p_stouffer_BH) & 
    p_stouffer_BH < P_BH_CUT, stouffer_call, "Not significant"))]
  data.table::setorder(genes, plot_order, PAO1_ID)
  # 4. Count same-direction calls across the five designated T1 comparisons.
  # The T2 follow-up is reported but does not count toward recurrence. A gene can
  # have recurrent evidence in each direction; the two counts are kept separately.
  recurrence <- genes[, .(n_independent_T1_tested = sum(count_for_recurrence & stouffer_evaluated), 
    n_independent_T1_sensitizing = sum(count_for_recurrence & stouffer_call == "Sensitizing"), 
    n_independent_T1_protective = sum(count_for_recurrence & stouffer_call == "Protective"), 
    n_BH_T1_sensitizing = sum(count_for_recurrence & stouffer_call_BH == "Sensitizing"), n_BH_T1_protective = sum(count_for_recurrence & 
      stouffer_call_BH == "Protective")), by = PAO1_ID]
  recurrence[, `:=`(recurrent_sensitizing = n_independent_T1_sensitizing >= RECURRENCE_MIN, recurrent_protective = n_independent_T1_protective >= 
    RECURRENCE_MIN, BH_recurrent_sensitizing = n_BH_T1_sensitizing >= RECURRENCE_MIN, BH_recurrent_protective = n_BH_T1_protective >= 
    RECURRENCE_MIN)]
  # 5. Build one row per gene, with separate effect/p-value/call columns for each
  # comparison. Keep a complete recurrence table and a hit table containing genes
  # called in at least one comparison; these are different views of the same calls.
  wide <- merge(annotation[, .(PAO1_ID, gene_name, product_name, essential_category_code = category, 
    TU)], recurrence, by = "PAO1_ID", all.y = TRUE, sort = FALSE)
  for (i in seq_len(nrow(manifest))) {
    z <- genes[comparison_id == manifest$comparison_id[i]]
    idx <- match(wide$PAO1_ID, z$PAO1_ID)
    wide[, `:=`((paste0(manifest$prefix[i], c("_medLFC", "_p_stouffer", "_call", "_p_stouffer_BH"))), 
      list(z$medLFC[idx], z$p_stouffer[idx], z$stouffer_call[idx], z$p_stouffer_BH[idx]))]
  }
  called_ids <- unique(genes$PAO1_ID[genes$stouffer_call != "Not significant"])
  calls <- wide[PAO1_ID %in% called_ids]
  prefixes <- manifest$prefix
  table_columns <- c("PAO1_ID", "gene_name", "product_name", "essential_category_code", "TU", 
    "n_independent_T1_sensitizing", "n_independent_T1_protective", "recurrent_sensitizing", 
    "recurrent_protective", paste0(prefixes, "_medLFC"), paste0(prefixes, "_p_stouffer"), paste0(prefixes, 
      "_call"), paste0(prefixes[manifest$count_for_recurrence], "_p_stouffer_BH"))
  calls <- calls[, ..table_columns]
  data.table::setorder(calls, -recurrent_sensitizing, -recurrent_protective, -n_independent_T1_sensitizing, 
    -n_independent_T1_protective, PAO1_ID)
  data.table::setorder(wide, PAO1_ID)
  data.table::setorder(thresholds, plot_order)
  # 6. Summarize how many primary recurrent hits also survive the additional BH
  # filter. This summary does not replace the primary raw-p-value screening rules.
  bh_summary <- data.table::data.table(direction = c("Sensitizing", "Protective"), original_recurrent = c(sum(wide$recurrent_sensitizing), 
    sum(wide$recurrent_protective)), BH_supported_recurrent = c(sum(wide$BH_recurrent_sensitizing), 
    sum(wide$BH_recurrent_protective)))
  bh_summary[, `:=`(retained_percent, ifelse(original_recurrent > 0, 100 * BH_supported_recurrent/original_recurrent, 
    NA_real_))]
  settings <- s8_settings(c(P_STOUFFER_CUT = P_STOUFFER_CUT, MAD_MULTIPLIER = MAD_MULTIPLIER, 
    MIN_ABS_MEDLFC = MIN_ABS_MEDLFC, P_BH_CUT = P_BH_CUT, RECURRENCE_MIN = RECURRENCE_MIN, 
    EXPORT_GUIDE_TABLES = EXPORT_GUIDE_TABLES, recurrence_conditions = paste(manifest$display_label[manifest$count_for_recurrence], 
      collapse = "; "), BH_family = "All finite named-gene p-values within each complete comparison"), 
    c(files, annotation = ANNOTATION_FILE))
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(".s8_", tmpdir = OUTPUT_DIR)
  dir.create(stage)
  on.exit(if (!length(list.files(file.path(stage, "backup")))) unlink(stage, recursive = TRUE))
  # 7. Optionally prepare the individual-guide subsets used downstream.
  # Stage all tables together, then install them only after export succeeds.
  # Settings and input checksums identify the files underlying these calls.
  if (EXPORT_GUIDE_TABLES) {
    guide_audit <- lapply(names(files), function(b) {
      path <- s8_file(INPUT_DIR, b, "_melted_results")
      s8_stream_guides(path, file.path(stage, paste0(b, "_AG_guides.tsv.gz")), b, manifest[batch == 
        b], GUIDE_CHUNK_SIZE)
    })
    s8_stage(list(Guide_extraction = data.table::rbindlist(guide_audit)), stage)
  }
  s8_stage(list(AG_gene_results = genes, AG_gene_recurrence = wide, Stouffer_calls = calls, AG_thresholds = thresholds, 
    BH_summary = bh_summary, Comparison_manifest = manifest, Analysis_settings = settings), 
    stage, gzip = "AG_gene_results")
  s8_commit(stage, OUTPUT_DIR, OVERWRITE)
  message("S8 Stouffer: ", nrow(calls), " genes with a call; recurrent ", sum(wide$recurrent_sensitizing), 
    " sensitizing / ", sum(wide$recurrent_protective), " protective.")
  invisible(list(genes = genes, calls = calls, recurrence = wide))
}

main()

