# Compare primary Stouffer interactions with completed MAGeCK RRA results.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/04_AMG_compare_MAGeCK_Stouffer.R")
# Inputs:
# - preprocessing_scripts/out_03_AMG_relevant_results/ — gene results, recurrence, manifest and
#   settings.
# - MAGeCK_data/mageck_rra_postprocessing_AG/MAGeCK_drug_contrasts_all_genes.tsv.gz (complete
#   all-gene catalogue).
# Upstream / execution order:
# - preprocessing_scripts/03_AMG_relevant_results.R;
#   MAGeCK_data/01_MAGeCK_preprocessing_and_run.R → generated Bash runner →
#   MAGeCK_data/02_MAGeCK_postprocessing.R.
# - The included MAGeCK catalogue can be reused without rerunning MAGeCK.
# Outputs:
# - preprocessing_scripts/out_04_AMG_compare_MAGeCK_Stouffer/ — counts_thresholds.tsv,
#   MAGeCK_calls.tsv, MAGeCK_vs_Stouffer.tsv, All_gene_comparisons.tsv.gz and settings.
# Interpretation / operation:
# - Join by PAO1 locus and the six matched drug-versus-control comparisons.
# - Directional MAGeCK FDR < 0.05; bidirectional significance is Ambiguous.
# - This comparison does not redefine the primary Stouffer recurrent sets.

STOUFFER_DIR <- "preprocessing_scripts/out_03_AMG_relevant_results"

MAGECK_FILE_CANDIDATES <- c("MAGeCK_data/MAGeCK_drug_contrasts_all_genes.tsv.gz", "MAGeCK_data/mageck_rra_postprocessing_AG/MAGeCK_drug_contrasts_all_genes.tsv.gz")

OUTPUT_DIR <- "preprocessing_scripts/out_04_AMG_compare_MAGeCK_Stouffer"

MAGECK_FDR_CUT <- 0.05

OVERWRITE <- TRUE

s8_method_comparison <- function(stouffer, mageck) {
  data.table::fcase(stouffer == "Sensitizing" & mageck == "Sensitizing", "Shared sensitizing", 
    stouffer == "Protective" & mageck == "Protective", "Shared protective", stouffer == "Sensitizing" & 
      mageck == "Not significant", "Stouffer only sensitizing", stouffer == "Protective" & 
      mageck == "Not significant", "Stouffer only protective", stouffer == "Not significant" & 
      mageck == "Sensitizing", "MAGeCK only sensitizing", stouffer == "Not significant" & 
      mageck == "Protective", "MAGeCK only protective", stouffer %in% c("Sensitizing", "Protective") & 
      mageck %in% c("Sensitizing", "Protective") & stouffer != mageck, "Direction conflict", 
    stouffer == "Not significant" & mageck == "Ambiguous", "MAGeCK ambiguous", stouffer %in% 
      c("Sensitizing", "Protective") & mageck == "Ambiguous", "Stouffer call; MAGeCK ambiguous", 
    default = "No call")
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
  if (!is.finite(MAGECK_FDR_CUT) || MAGECK_FDR_CUT <= 0 || MAGECK_FDR_CUT > 1) 
    stop("Invalid MAGeCK FDR cutoff.")
  # 1. Read the primary Stouffer calls, their effect thresholds and comparison list.
  # Use that list to match MAGeCK to the same treatments, batches and time points.
  manifest <- s8_read(file.path(STOUFFER_DIR, "Comparison_manifest.tsv"))
  st <- s8_read(file.path(STOUFFER_DIR, "AG_gene_results.tsv.gz"))
  thresholds <- s8_read(file.path(STOUFFER_DIR, "AG_thresholds.tsv"))
  recurrence_min <- s8_setting(STOUFFER_DIR, "RECURRENCE_MIN", numeric = TRUE)
  s8_required(st, c("comparison_id", "PAO1_ID", "gene_name", "medLFC", "p_stouffer", "stouffer_call", 
    "stouffer_evaluated"), "03 gene results")
  paths <- unique(unlist(lapply(MAGECK_FILE_CANDIDATES, function(p) c(p, sub("\\.gz$", "", p)))))
  found <- paths[file.exists(paths)]
  if (length(found) != 1L) 
    stop("Keep one complete MAGeCK input. Checked:\n", paste(paths, collapse = "\n"))
  raw <- s8_read(found)
  required <- c("PAO1_ID", "batch", "shift_drug", "shift_dose", "shift_timing", "base_drug", 
    "base_timing", "comparison_type", "gene_lfc", "neg_fdr", "pos_fdr")
  s8_required(raw, required, found)
  raw[, `:=`(PAO1_ID = s8_clean(PAO1_ID), batch = tolower(trimws(batch)), shift_drug = toupper(trimws(shift_drug)), 
    shift_timing = toupper(trimws(shift_timing)), dose_numeric = suppressWarnings(as.numeric(shift_dose)))]
  # 2. Select each drug-versus-matched-untreated MAGeCK comparison.
  # Use directional MAGeCK FDRs alone: negative selection is sensitizing and
  # positive selection is protective. Significant results in both directions
  # are marked ambiguous; no additional gene-LFC cutoff is applied here.
  parts <- lapply(seq_len(nrow(manifest)), function(i) {
    spec <- manifest[i]
    x <- raw[batch == spec$batch & shift_drug == spec$drug & is.finite(dose_numeric) & dose_numeric == 
      spec$dose & shift_timing == spec$timing]
    x <- x[comparison_type == "drug_vs_matched_control" & !is.na(base_drug) & toupper(trimws(base_drug)) == 
      "NONE" & !is.na(base_timing) & toupper(trimws(base_timing)) == spec$timing]
    if ("base_dose" %in% names(x)) 
      x <- x[!is.na(base_dose) & suppressWarnings(as.numeric(base_dose)) == 0]
    if (!nrow(x)) 
      stop("No matched MAGeCK rows for ", spec$comparison_id)
    s8_unique(x, "PAO1_ID", paste("MAGeCK", spec$comparison_id))
    s8_numeric(x, "gene_lfc", spec$comparison_id)
    s8_probability(x, c("neg_fdr", "pos_fdr"), spec$comparison_id)
    x[, `:=`(comparison_id = spec$comparison_id, mageck_evaluated = is.finite(neg_fdr) | is.finite(pos_fdr))]
    x[, `:=`(mageck_call, data.table::fcase(is.finite(neg_fdr) & neg_fdr < MAGECK_FDR_CUT & 
      is.finite(pos_fdr) & pos_fdr < MAGECK_FDR_CUT, "Ambiguous", is.finite(neg_fdr) & neg_fdr < 
      MAGECK_FDR_CUT, "Sensitizing", is.finite(pos_fdr) & pos_fdr < MAGECK_FDR_CUT, "Protective", 
      default = "Not significant"))]
    x[, .(comparison_id, PAO1_ID, gene_lfc, neg_fdr, pos_fdr, mageck_evaluated, mageck_call)]
  })
  mg <- data.table::rbindlist(parts)
  s8_unique(st, c("comparison_id", "PAO1_ID"), "Stouffer results")
  # 3. Join methods by locus and comparison, retaining genes found by either.
  # Keep evaluated/not-evaluated flags so missing method results remain identifiable
  # even though absent call labels are represented as Not significant.
  pair <- merge(st[, .(comparison_id, PAO1_ID, gene_name, medLFC_stouffer = medLFC, p_stouffer, 
    stouffer_call, stouffer_evaluated)], mg, by = c("comparison_id", "PAO1_ID"), all = TRUE, 
    sort = FALSE)
  annotation <- s8_annotation(s8_setting(STOUFFER_DIR, "input:annotation"))
  metadata <- annotation[, .(PAO1_ID, gene_name)]
  if (length(setdiff(mg$PAO1_ID, metadata$PAO1_ID))) 
    stop("MAGeCK loci are absent from the curated annotation.")
  pair[is.na(gene_name), `:=`(gene_name, metadata$gene_name[match(PAO1_ID, metadata$PAO1_ID)])]
  pair[is.na(stouffer_evaluated), `:=`(stouffer_evaluated, FALSE)]
  pair[is.na(mageck_evaluated), `:=`(mageck_evaluated, FALSE)]
  pair[is.na(stouffer_call), `:=`(stouffer_call, "Not significant")]
  pair[is.na(mageck_call), `:=`(mageck_call, "Not significant")]
  pair[, `:=`(method_comparison = s8_method_comparison(stouffer_call, mageck_call), shared_same_direction = stouffer_call %in% 
    c("Sensitizing", "Protective") & stouffer_call == mageck_call, any_call = stouffer_call != 
    "Not significant" | mageck_call != "Not significant")]
  pair <- merge(pair, manifest[, .(comparison_id, display_label, plot_order, count_for_recurrence)], 
    by = "comparison_id", all.x = TRUE, sort = FALSE)
  # 4. Build a complete comparison grid for genes called by at least one method.
  # Classify agreement, method-specific calls and conflicting directions without
  # changing which genes belong to the primary Stouffer hit sets.
  union_ids <- unique(pair$PAO1_ID[pair$any_call])
  grid <- data.table::CJ(PAO1_ID = sort(union_ids), comparison_id = manifest$comparison_id, unique = TRUE)
  calls <- merge(grid, pair, by = c("PAO1_ID", "comparison_id"), all.x = TRUE, sort = FALSE)
  idx <- match(calls$comparison_id, manifest$comparison_id)
  calls[, `:=`(plot_order = manifest$plot_order[idx], count_for_recurrence = manifest$count_for_recurrence[idx])]
  calls[is.na(gene_name), `:=`(gene_name, metadata$gene_name[match(PAO1_ID, metadata$PAO1_ID)])]
  for (column in c("stouffer_evaluated", "mageck_evaluated", "shared_same_direction", "any_call")) data.table::set(calls, 
    which(is.na(calls[[column]])), column, FALSE)
  for (column in c("stouffer_call", "mageck_call")) data.table::set(calls, which(is.na(calls[[column]])), 
    column, "Not significant")
  calls[is.na(method_comparison), `:=`(method_comparison, "No call")]
  data.table::setorder(calls, plot_order, PAO1_ID)
  # 5. Count T1 recurrence separately for each method and for same-direction
  # agreement in the same comparisons. T2 remains outside recurrence counts.
  method_recurrence <- calls[, .(gene_name = gene_name[which(!is.na(gene_name))[1]], n_stouffer_T1_sensitizing = sum(count_for_recurrence & 
    stouffer_call == "Sensitizing"), n_stouffer_T1_protective = sum(count_for_recurrence & 
    stouffer_call == "Protective"), n_mageck_T1_sensitizing = sum(count_for_recurrence & mageck_call == 
    "Sensitizing"), n_mageck_T1_protective = sum(count_for_recurrence & mageck_call == "Protective"), 
    n_shared_T1_sensitizing = sum(count_for_recurrence & stouffer_call == "Sensitizing" & mageck_call == 
      "Sensitizing"), n_shared_T1_protective = sum(count_for_recurrence & stouffer_call == 
      "Protective" & mageck_call == "Protective")), by = PAO1_ID]
  for (method in c("stouffer", "mageck", "shared")) for (direction in c("sensitizing", "protective")) method_recurrence[, 
    `:=`((paste0(method, "_recurrent_T1_", direction)), get(paste0("n_", method, "_T1_", direction)) >= 
      recurrence_min)]
  data.table::setorder(method_recurrence, PAO1_ID)
  # 6. Count shared and method-specific hits per comparison, retaining the
  # Stouffer effect cutoffs beside them for interpretation.
  counts <- lapply(seq_len(nrow(manifest)), function(i) {
    id <- manifest$comparison_id[i]
    a <- st[comparison_id == id]
    b <- mg[comparison_id == id]
    z <- calls[comparison_id == id]
    data.table::data.table(comparison_id = id, display_label = manifest$display_label[i], n_stouffer_evaluated = sum(a$stouffer_evaluated), 
      n_stouffer_sensitizing = sum(a$stouffer_call == "Sensitizing"), n_stouffer_protective = sum(a$stouffer_call == 
        "Protective"), n_mageck_evaluated = sum(b$mageck_evaluated), n_mageck_sensitizing = sum(b$mageck_call == 
        "Sensitizing"), n_mageck_protective = sum(b$mageck_call == "Protective"), n_shared_sensitizing = sum(z$method_comparison == 
        "Shared sensitizing"), n_shared_protective = sum(z$method_comparison == "Shared protective"), 
      n_stouffer_only_sensitizing = sum(z$method_comparison == "Stouffer only sensitizing"), 
      n_stouffer_only_protective = sum(z$method_comparison == "Stouffer only protective"), 
      n_mageck_only_sensitizing = sum(z$method_comparison == "MAGeCK only sensitizing"), 
      n_mageck_only_protective = sum(z$method_comparison == "MAGeCK only protective"), n_union_called = sum(z$any_call))
  })
  counts <- data.table::rbindlist(counts)
  threshold_columns <- c("comparison_id", "p_stouffer_cutoff", "median_medLFC", "MAD", "sensitizing_cutoff_source", 
    "sensitizing_LFC_cutoff_used", "protective_LFC_cutoff_used")
  counts <- merge(counts, thresholds[, ..threshold_columns], by = "comparison_id", all.x = TRUE, 
    sort = FALSE)
  counts <- counts[match(manifest$comparison_id, comparison_id)]
  table_columns <- c("comparison_id", "PAO1_ID", "gene_name", "mageck_gene_lfc", "mageck_neg_fdr", 
    "mageck_pos_fdr", "mageck_call", "medLFC_stouffer", "p_stouffer", "stouffer_call", "stouffer_evaluated", 
    "mageck_evaluated", "method_comparison", "shared_same_direction", "any_call")
  data.table::setnames(calls, c("gene_lfc", "neg_fdr", "pos_fdr"), c("mageck_gene_lfc", "mageck_neg_fdr", 
    "mageck_pos_fdr"))
  calls <- calls[, ..table_columns]
  settings <- s8_settings(c(MAGECK_FDR_CUT = MAGECK_FDR_CUT, RECURRENCE_MIN = recurrence_min, 
    gene_LFC_filter = "none", primary_gene_sets = "Stouffer; MAGeCK does not redefine them"), 
    c(MAGeCK = found, Stouffer = file.path(STOUFFER_DIR, "AG_gene_results.tsv.gz"), call_settings = file.path(STOUFFER_DIR, 
      "Analysis_settings.tsv"), manifest = file.path(STOUFFER_DIR, "Comparison_manifest.tsv")))
  # 7. Export the hit grid, per-gene recurrence and comparison-level counts.
  # All_gene_comparisons also retains the complete paired-method results.
  s8_write(list(MAGeCK_calls = calls, MAGeCK_vs_Stouffer = method_recurrence, counts_thresholds = counts, 
    All_gene_comparisons = pair, Analysis_settings = settings), OUTPUT_DIR, OVERWRITE, gzip = "All_gene_comparisons")
  message("S8 MAGeCK comparison: ", length(union_ids), " genes in the six-comparison publication grid.")
  invisible(list(calls = calls, recurrence = method_recurrence, counts = counts))
}

main()

