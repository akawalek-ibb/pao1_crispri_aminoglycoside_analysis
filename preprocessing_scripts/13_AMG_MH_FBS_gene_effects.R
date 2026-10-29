# Purpose: Prepare matched B3 MH/FBS gene effects and vulnerability for Table S7.
# Prerequisites: 02_AMG_metacounter.R and 11_AMG_gene_vulnerability.R.
# Inputs: B3 median gene results and the shared vulnerability long table.
# Direct effects: FBS T1 versus MH T1; absolute profiles: each versus MH T0.
# No generation-number scaling is applied. No vulnerability curves are refitted.
# Outputs: MH_FBS_gene_effects.tsv, the two Stouffer/vulnerability table sheets,
# Summary.tsv and Analysis_settings.tsv in out_13_AMG_MH_FBS_gene_effects/.
# Run from AMG_cleaned_Git.Rproj:
# source("preprocessing_scripts/13_AMG_MH_FBS_gene_effects.R")
# USER SETTINGS
INPUT_DIR <- "edgeR_data"

ANNOTATION_FILE <- "annotation_input/PAO1_curated.tsv"

VULNERABILITY_FILE <- "preprocessing_scripts/out_11_AMG_gene_vulnerability/gene_vulnerability_long.tsv.gz"

OUTPUT_DIR <- "preprocessing_scripts/out_13_AMG_MH_FBS_gene_effects"

P_STOUFFER_CUT <- 0.05

LFC_CUT <- 1

OVERWRITE <- TRUE

COND_MH <- "None_0_T1 - None_0_T0"

COND_FBS <- "FBS_60_T1 - None_0_T0"

COND_DIRECT <- "FBS_60_T1 - None_0_T1"

MH_VULNERABILITY_ID <- "B3_None_T1_vs_T0"

FBS_VULNERABILITY_ID <- "B3_FBS_60_T1_vs_None_T0"

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

amg_annotation <- function(path) {
  x <- amg_read(path)
  amg_required(x, c("PAO1_ID", "PA14_ID", "gene_name", "product_name", "category", "design_call", 
    "TU"), path)
  x[, `:=`(PAO1_ID, trimws(as.character(PAO1_ID)))]
  amg_unique(x, "PAO1_ID", path)
  x
}

amg_gene_file <- function(directory, batch) {
  paths <- file.path(directory, paste0(batch, "_median_melted_results", c(".tsv.gz", ".tsv")))
  existing <- paths[file.exists(paths)]
  if (!length(existing)) 
    stop("Missing 02_AMG_metacounter gene results:\n", paste(paths, collapse = "\n"), call. = FALSE)
  existing[1L]
}

amg_gene_results <- function(path, conditions) {
  x <- amg_read(path)
  amg_required(x, c("PAO1_ID", "condition", "medLFC", "p_stouffer"), path)
  if ("type" %in% names(x)) 
    x <- x[tolower(as.character(type)) == "perfect"]
  x[, `:=`(PAO1_ID = trimws(as.character(PAO1_ID)), condition = trimws(as.character(condition)))]
  x <- x[condition %in% conditions & !is.na(PAO1_ID) & nzchar(PAO1_ID)]
  missing <- setdiff(conditions, unique(x$condition))
  if (length(missing)) 
    stop(path, " lacks comparisons:\n", paste(missing, collapse = "\n"), call. = FALSE)
  amg_unique(x, c("PAO1_ID", "condition"), path)
  amg_numeric(x, c("medLFC", "p_stouffer"), path)
  if (any(x$p_stouffer < 0 | x$p_stouffer > 1, na.rm = TRUE)) 
    stop(path, ": p_stouffer must be between 0 and 1.")
  x
}

amg_write <- function(tables, directory, overwrite = TRUE, gzip = character()) {
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

amg_settings <- function(settings, inputs) {
  data.table::data.table(setting = c(names(settings), paste0("input:", names(inputs)), paste0("md5:", 
    names(inputs)), "R_version"), value = c(as.character(settings), unname(inputs), unname(tools::md5sum(inputs)), 
    R.version.string))
}

main <- function() {
  amg_require("data.table")
  if (!is.finite(P_STOUFFER_CUT) || P_STOUFFER_CUT <= 0 || P_STOUFFER_CUT > 1 || !is.finite(LFC_CUT) || 
    LFC_CUT <= 0) 
    stop("Invalid direct-comparison thresholds.")
  # 1. Read three B3 contrasts: MH versus T0, FBS versus T0 and direct FBS versus
  # matched MH. The first two describe fitness within each medium; the third
  # tests a difference between media directly.
  path <- amg_gene_file(INPUT_DIR, "b3")
  raw <- amg_gene_results(path, c(COND_MH, COND_FBS, COND_DIRECT))
  annotation <- amg_annotation(ANNOTATION_FILE)
  ids <- unique(raw$PAO1_ID)
  if (length(setdiff(ids, annotation$PAO1_ID))) 
    stop("B3 results contain loci absent from curated annotation.")
  genes <- annotation[match(ids, PAO1_ID), .(PAO1_ID, gene_name)]
  # 2. Align all three effects and p-values by locus ID in one row per gene.
  # Keep missing tests visible rather than replacing them with zero effects.
  for (label in c("MH", "FBS", "direct_FBS_vs_MH")) {
    condition_value <- switch(label, MH = COND_MH, FBS = COND_FBS, direct_FBS_vs_MH = COND_DIRECT)
    z <- raw[condition == condition_value]
    idx <- match(genes$PAO1_ID, z$PAO1_ID)
    suffix <- if (label == "direct_FBS_vs_MH") 
      ""
    else "_vs_T0"
    genes[, `:=`((c(paste0(label, "_medLFC", suffix), paste0(label, "_p_stouffer", suffix))), 
      list(z$medLFC[idx], z$p_stouffer[idx]))]
  }
  # 3. Call a direct medium difference using its own p-value and absolute effect
  # cutoff. A positive direct FBS-minus-MH effect indicates higher relative guide
  # abundance in FBS; it need not mean a positive fitness effect versus T0.
  genes[, `:=`(direct_significant, is.finite(direct_FBS_vs_MH_p_stouffer) & direct_FBS_vs_MH_p_stouffer < 
    P_STOUFFER_CUT & is.finite(direct_FBS_vs_MH_medLFC) & abs(direct_FBS_vs_MH_medLFC) >= LFC_CUT)]
  vulnerability <- amg_read(VULNERABILITY_FILE)
  required <- c("PAO1_ID", "comparison_id", "condition", "batch", "severity", "sensitivity_y_pred", 
    "sensitivity_status", "EK50", "EK50_status")
  amg_required(vulnerability, required, VULNERABILITY_FILE)
  # 4. Import B3 sensitivity/EK50 estimates from the task-11 master separately
  # for each medium. Check contrast IDs, estimate statuses and severity agreement
  # before joining; no curves are refitted here.
  for (label in c("MH", "FBS")) {
    id <- if (label == "MH") 
      MH_VULNERABILITY_ID
    else FBS_VULNERABILITY_ID
    expected <- if (label == "MH") 
      COND_MH
    else COND_FBS
    z <- vulnerability[comparison_id == id]
    if (!nrow(z)) 
      stop("Missing vulnerability comparison: ", id)
    if (any(tolower(z$batch) != "b3") || any(z$condition != expected)) 
      stop("Vulnerability comparison does not match its expected B3 contrast: ", id)
    amg_unique(z, "PAO1_ID", id)
    amg_numeric(z, c("severity", "sensitivity_y_pred", "EK50"), id)
    if (any(!is.na(z$sensitivity_y_pred) & z$sensitivity_status != "estimated") || any(!is.na(z$EK50) & 
      z$EK50_status != "estimated")) 
      stop("Numeric vulnerability estimates have inconsistent statuses: ", id)
    idx <- match(genes$PAO1_ID, z$PAO1_ID)
    if (anyNA(idx)) 
      stop("The vulnerability master lacks some B3 loci; run 11 on the same annotation and metacounter inputs.")
    median_column <- paste0(label, "_medLFC_vs_T0")
    comparable <- is.finite(genes[[median_column]]) & is.finite(z$severity[idx])
    if (any(is.finite(genes[[median_column]]) != is.finite(z$severity[idx])) || any(abs(genes[[median_column]][comparable] - 
      z$severity[idx][comparable]) > 1e-08)) 
      stop("Severity differs between metacounter and vulnerability master: ", label, ". Re-run preprocessing 11 on the current inputs.")
    genes[, `:=`((paste0(label, c("_sensitivity", "_sensitivity_status", "_EK50", "_EK50_status"))), 
      list(z$sensitivity_y_pred[idx], z$sensitivity_status[idx], z$EK50[idx], z$EK50_status[idx]))]
  }
  # 5. Distinguish numeric estimates, evaluable profiles with no sustained
  # threshold crossing, and profiles lacking enough information in either medium.
  # These statuses allow downstream figures to treat boundary points appropriately.
  genes[, `:=`(sensitivity_status_pair, data.table::fcase(MH_sensitivity_status == "estimated" & 
    FBS_sensitivity_status == "estimated", "Estimated in both", MH_sensitivity_status == "no_sustained_crossing" & 
    FBS_sensitivity_status == "no_sustained_crossing", "No crossing in both", MH_sensitivity_status == 
    "estimated" & FBS_sensitivity_status == "no_sustained_crossing", "MH estimated; FBS no crossing", 
    MH_sensitivity_status == "no_sustained_crossing" & FBS_sensitivity_status == "estimated", 
    "MH no crossing; FBS estimated", default = "Non-evaluable in at least one condition"))]
  genes[, `:=`(EK50_status_pair, ifelse(MH_EK50_status == "estimated" & FBS_EK50_status == "estimated", 
    "Estimated in both", "Non-evaluable in at least one condition"))]
  genes[, `:=`(abs_direct, abs(direct_FBS_vs_MH_medLFC))]
  data.table::setorderv(genes, c("direct_FBS_vs_MH_p_stouffer", "abs_direct", "PAO1_ID"), c(1L, 
    -1L, 1L), na.last = TRUE)
  genes[, `:=`(abs_direct, NULL)]
  significant_columns <- c("PAO1_ID", "gene_name", "direct_FBS_vs_MH_medLFC", "direct_FBS_vs_MH_p_stouffer", 
    "MH_medLFC_vs_T0", "MH_p_stouffer_vs_T0", "FBS_medLFC_vs_T0", "FBS_p_stouffer_vs_T0")
  absolute_columns <- c("PAO1_ID", "gene_name", "MH_medLFC_vs_T0", "MH_p_stouffer_vs_T0", "FBS_medLFC_vs_T0", 
    "FBS_p_stouffer_vs_T0", "MH_sensitivity", "MH_sensitivity_status", "FBS_sensitivity", "FBS_sensitivity_status", 
    "sensitivity_status_pair", "MH_EK50", "MH_EK50_status", "FBS_EK50", "FBS_EK50_status", 
    "EK50_status_pair")
  # 6. Export both the direct significant-gene subset and the complete paired
  # severity/sensitivity/EK50 view, with provenance for the shared Table S7 inputs.
  hits <- genes[direct_significant == TRUE, ..significant_columns]
  absolute <- genes[, ..absolute_columns]
  data.table::setorderv(hits, c("direct_FBS_vs_MH_medLFC", "direct_FBS_vs_MH_p_stouffer", "PAO1_ID"), 
    na.last = TRUE)
  data.table::setorderv(absolute, c("MH_medLFC_vs_T0", "PAO1_ID"), na.last = TRUE)
  summary <- data.table::data.table(metric = c("n_direct_tested", "n_Stouffer_significant", "n_absolute_records", 
    "n_sensitivity_estimated_both", "n_EK50_estimated_both"), value = c(sum(is.finite(genes$direct_FBS_vs_MH_medLFC) & 
    is.finite(genes$direct_FBS_vs_MH_p_stouffer)), nrow(hits), nrow(absolute), sum(genes$sensitivity_status_pair == 
    "Estimated in both"), sum(genes$EK50_status_pair == "Estimated in both")))
  settings <- amg_settings(c(P_STOUFFER_CUT = P_STOUFFER_CUT, LFC_CUT = LFC_CUT, COND_MH = COND_MH, 
    COND_FBS = COND_FBS, COND_DIRECT = COND_DIRECT, MH_VULNERABILITY_ID = MH_VULNERABILITY_ID, 
    FBS_VULNERABILITY_ID = FBS_VULNERABILITY_ID), c(gene_results = path, annotation = ANNOTATION_FILE, 
    vulnerability = VULNERABILITY_FILE))
  amg_write(list(MH_FBS_gene_effects = genes, Stouffer_FBS_vs_MH_significant = hits, Sev_Sens_EK50_MH_FBS_vs_T0 = absolute, 
    Summary = summary, Analysis_settings = settings), OUTPUT_DIR, OVERWRITE)
  message("Prepared S7 gene effects: ", nrow(hits), " direct Stouffer hits; ", nrow(absolute), 
    " gene records. Saved to ", OUTPUT_DIR)
  invisible(genes)
}

main()

