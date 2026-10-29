# Standardize RNA-seq differential expression for four knockdowns and TOB alone.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/18_AMG_RNAseq_results.R")
# Inputs:
# - RNAseq_data/DE/PAO1_{atpIg1m2,accCg1m3,secAg1m3,ftsZg1m4,control1_TOB}_vs_control1_All.xlsx.
# - First worksheet: unique named PAO1_ID, gene_name, logFC, LR, PValue and FDR; use complete
#   results, not DEG-only files.
# Upstream / execution order:
# - The supplied DE workbooks are the starting point. Rebuilding starts with
#   RNAseq_analysis/AMG_RNAseq_featureCounts_edgeR.R.
# - That script writes separate filenames under RNAseq_analysis/out/; it does not populate this
#   input directory automatically.
# Outputs:
# - preprocessing_scripts/out_18_AMG_RNAseq_results/ — RNAseq_results.tsv.gz,
#   Comparison_manifest.tsv, RNAseq_summary.tsv and settings.
# Interpretation / operation:
# - Preserve source p-values/FDR; label DEGs at FDR < 0.05 and |log2FC| > 1.
# - GO ranks are sign(logFC) × sqrt(LR). TOB is retained in expression tables but excluded from
#   GO analyses.

INPUT_DIR <- "RNAseq_data/DE"

OUTPUT_DIR <- "preprocessing_scripts/out_18_AMG_RNAseq_results"

DE_FDR_CUT <- 0.05

DE_ABS_LOG2FC_CUT <- 1

OVERWRITE <- TRUE

COMPARISONS <- data.frame(comparison_id = c("atpI", "accC", "secA", "ftsZ", "TOB"), contrast_label = c("atpI silencing", 
  "accC silencing", "secA silencing", "ftsZ silencing", "TOB treatment"), file_name = c("PAO1_atpIg1m2_vs_control1_All.xlsx", 
  "PAO1_accCg1m3_vs_control1_All.xlsx", "PAO1_secAg1m3_vs_control1_All.xlsx", "PAO1_ftsZg1m4_vs_control1_All.xlsx", 
  "PAO1_control1_TOB_vs_control1_All.xlsx"), stringsAsFactors = FALSE)

# Local functions for this analysis.
amg_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) 
    stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
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
  amg_require(c("data.table", "openxlsx"))
  if (!is.finite(DE_FDR_CUT) || DE_FDR_CUT <= 0 || DE_FDR_CUT > 1 || !is.finite(DE_ABS_LOG2FC_CUT) || 
    DE_ABS_LOG2FC_CUT <= 0) 
    stop("Invalid DEG thresholds.")
  # 1. Locate complete differential-expression outputs for each planned contrast.
  # TOB remains in the RNA-seq results but is excluded from GO analyses by use_for_GO.
  # This script standardizes existing results; it does not reanalyse read counts.
  manifest <- data.table::as.data.table(COMPARISONS)
  amg_unique(manifest, "comparison_id", "Comparison manifest")
  manifest[, `:=`(plot_order = seq_len(.N), source_file = file.path(INPUT_DIR, file_name), use_for_GO = comparison_id != 
    "TOB")]
  if (any(!file.exists(manifest$source_file))) 
    stop("Missing complete edgeR inputs:\n", paste(manifest$source_file[!file.exists(manifest$source_file)], 
      collapse = "\n"))
  # 2. Read each complete comparison and recognize supported column names.
  # Retain tested and untested genes, checking locus uniqueness and valid effects,
  # p-values, FDRs and nonnegative likelihood-ratio statistics.
  parts <- lapply(seq_len(nrow(manifest)), function(i) {
    spec <- manifest[i]
    raw <- rnaseq_read(spec$source_file)
    id <- rnaseq_column(raw, c("PAO1_ID", "locus_tag", "Locus_Tag", "Geneid", "gene_id"), spec$source_file)
    gene <- rnaseq_column(raw, c("gene_name", "Gene", "gene"), spec$source_file, required = FALSE)
    fc <- rnaseq_column(raw, c("logFC", "log2FC", "log2FoldChange"), spec$source_file)
    p <- rnaseq_column(raw, c("PValue", "pvalue", "P_value", "p_value", "P.Value"), spec$source_file)
    fdr <- rnaseq_column(raw, c("FDR", "padj", "adj.P.Val", "adjusted_p_value"), spec$source_file)
    lr <- rnaseq_column(raw, "LR", spec$source_file)
    x <- data.table::data.table(PAO1_ID = rnaseq_clean(raw[[id]]), gene_name = if (is.na(gene)) 
      NA_character_
    else rnaseq_clean(raw[[gene]]), logFC = raw[[fc]], PValue = raw[[p]], FDR = raw[[fdr]], 
      LR = raw[[lr]], source_row = seq_len(nrow(raw)))
    amg_unique(x, "PAO1_ID", spec$source_file)
    amg_numeric(x, c("logFC", "LR"), spec$source_file)
    rnaseq_probability(x, c("PValue", "FDR"), spec$source_file)
    if (any(x$LR < 0, na.rm = TRUE)) 
      stop("Negative LR in ", spec$source_file)
    # 3. Preserve the supplied p-values/FDRs. Build a signed ranking statistic
    # from sign(logFC) * sqrt(LR) for downstream GSEA: positive scores favour
    # upregulation, while negative scores favour downregulation.
    x[, `:=`(comparison_id = spec$comparison_id, contrast_label = spec$contrast_label, plot_order = spec$plot_order, 
      use_for_GO = spec$use_for_GO, tested_for_DE = is.finite(PValue) & is.finite(logFC), 
      rank_score = sign(logFC) * sqrt(LR))]
    # 4. Assign DEG direction using FDR and a STRICT fold-change cutoff.
    # Effects exactly at +/-DE_ABS_LOG2FC_CUT are not called Up/Down; count these
    # boundary cases separately in the summary. Missing tests remain Not tested.
    x[, `:=`(DEG_direction, data.table::fcase(tested_for_DE & is.finite(FDR) & FDR < DE_FDR_CUT & 
      logFC > DE_ABS_LOG2FC_CUT, "Up", tested_for_DE & is.finite(FDR) & FDR < DE_FDR_CUT & 
      logFC < -DE_ABS_LOG2FC_CUT, "Down", !tested_for_DE, "Not tested", default = "Not significant"))]
    for (column in intersect(c("product_name", "strand", "start", "end", "logCPM"), names(raw))) data.table::set(x, 
      j = column, value = raw[[column]])
    x
  })
  # 5. Combine comparisons, check consistent gene names, and export a long table
  # with the manifest, comparison-level counts and source-file checksums.
  results <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  names_check <- results[!is.na(gene_name), .(n_names = data.table::uniqueN(gene_name)), by = PAO1_ID]
  if (any(names_check$n_names > 1L)) 
    stop("Gene names disagree between the supplied DE tables; reconcile the source annotation.")
  data.table::setorder(results, plot_order, source_row)
  summary <- results[, .(n_rows = .N, n_tested = sum(tested_for_DE), n_Up = sum(DEG_direction == 
    "Up"), n_Down = sum(DEG_direction == "Down"), n_FDR_only = sum(tested_for_DE & is.finite(FDR) & 
    FDR < DE_FDR_CUT), n_at_FC_boundary = sum(abs(logFC) == DE_ABS_LOG2FC_CUT, na.rm = TRUE)), 
    by = .(comparison_id, contrast_label, plot_order)]
  settings <- amg_settings(c(DE_FDR_CUT = DE_FDR_CUT, DE_ABS_LOG2FC_CUT = DE_ABS_LOG2FC_CUT, 
    DEG_effect_operator = "strict >", rank_method = "sign(logFC) * sqrt(LR); one-df edgeR likelihood-ratio contrast", 
    p_values = "preserved; no BH recalculation"), stats::setNames(manifest$source_file, manifest$comparison_id))
  rnaseq_write(list(RNAseq_results = results, Comparison_manifest = manifest, RNAseq_summary = summary, 
    Analysis_settings = settings), OUTPUT_DIR, OVERWRITE, gzip = "RNAseq_results")
  message("RNA-seq: ", data.table::uniqueN(results$PAO1_ID), " loci across ", nrow(manifest), 
    " comparisons.")
  invisible(results)
}

main()

