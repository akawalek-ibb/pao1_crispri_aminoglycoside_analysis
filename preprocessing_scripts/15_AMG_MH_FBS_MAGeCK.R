# Import the matched FBS-versus-MH MAGeCK comparison for Table S7.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/15_AMG_MH_FBS_MAGeCK.R")
# Inputs:
# - MAGeCK_data/mageck_rra_postprocessing_AG/MAGeCK_drug_contrasts_all_genes.tsv.gz: complete
#   gene results.
# - annotation_input/PAO1_curated.tsv: unique PAO1_ID and gene annotations.
# Upstream / execution order:
# - Reuse the included catalogue, or run MAGeCK_data/01_MAGeCK_preprocessing_and_run.R →
#   generated Bash runner → MAGeCK_data/02_MAGeCK_postprocessing.R.
# Outputs:
# - preprocessing_scripts/out_15_AMG_MH_FBS_MAGeCK/ — MAGeCK_FBS_vs_MH_all_genes.tsv,
#   MAGeCK_FBS_vs_MH_significant.tsv, Summary.tsv and settings.
# Interpretation / operation:
# - Select B3 FBS 60 T1 versus None 0 T1. Directional FDR < 0.05; no LFC cutoff; both directions
#   is ambiguous.

MAGECK_FILE_CANDIDATES <- c("MAGeCK_data/MAGeCK_drug_contrasts_all_genes.tsv.gz", "MAGeCK_data/mageck_rra_postprocessing_AG/MAGeCK_drug_contrasts_all_genes.tsv.gz")

ANNOTATION_FILE <- "annotation_input/PAO1_curated.tsv"

OUTPUT_DIR <- "preprocessing_scripts/out_15_AMG_MH_FBS_MAGeCK"

MAGECK_FDR_CUT <- 0.05

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
  if (!is.finite(MAGECK_FDR_CUT) || MAGECK_FDR_CUT <= 0 || MAGECK_FDR_CUT > 1) 
    stop("Invalid MAGeCK FDR cutoff.")
  # 1. Locate one complete processed MAGeCK catalogue. This script extracts and
  # formats existing RRA results; it does not run MAGeCK or refit count models.
  paths <- unique(unlist(lapply(MAGECK_FILE_CANDIDATES, function(path) c(path, sub("\\.gz$", 
    "", path)))))
  found <- paths[file.exists(paths)]
  if (length(found) != 1L) 
    stop("Keep one MAGeCK all-gene input. Checked:\n", paste(paths, collapse = "\n"))
  mg <- amg_read(found)
  required <- c("PAO1_ID", "batch", "shift_drug", "shift_dose", "shift_timing", "gene_lfc", "neg_fdr", 
    "pos_fdr", "neg_pvalue", "pos_pvalue", "neg_score", "pos_score", "neg_rank", "pos_rank", 
    "num_guides")
  amg_required(mg, required, found)
  dose <- suppressWarnings(as.numeric(mg$shift_dose))
  # 2. Select the B3 60% FBS T1 versus untreated MH T1 comparison.
  # Check the reference group too, preventing accidental use of FBS versus T0.
  keep <- tolower(trimws(mg$batch)) == "b3" & toupper(trimws(mg$shift_drug)) == "FBS" & is.finite(dose) & 
    dose == 60 & toupper(trimws(mg$shift_timing)) == "T1"
  mg <- mg[which(keep)]
  if ("comparison_type" %in% names(mg)) 
    mg <- mg[comparison_type == "drug_vs_matched_control"]
  if ("base_drug" %in% names(mg)) 
    mg <- mg[!is.na(base_drug) & toupper(trimws(base_drug)) == "NONE"]
  if ("base_timing" %in% names(mg)) 
    mg <- mg[!is.na(base_timing) & toupper(trimws(base_timing)) == "T1"]
  if ("base_dose" %in% names(mg)) 
    mg <- mg[!is.na(base_dose) & suppressWarnings(as.numeric(base_dose)) == 0]
  if (!nrow(mg)) 
    stop("No B3 FBS 60 T1 versus MH T1 rows found in the MAGeCK input.")
  if (!"comparison_type" %in% names(mg) && !all(c("base_drug", "base_timing") %in% names(mg))) 
    stop("MAGeCK input must identify the matched control using comparison_type or base_drug/base_timing.")
  if ("contrast_id" %in% names(mg) && data.table::uniqueN(mg$contrast_id) != 1L) 
    stop("More than one FBS-versus-MH MAGeCK contrast remains.")
  mg[, `:=`(PAO1_ID, trimws(as.character(PAO1_ID)))]
  amg_unique(mg, "PAO1_ID", "Selected MAGeCK contrast")
  numeric_columns <- c("gene_lfc", "neg_fdr", "pos_fdr", "neg_pvalue", "pos_pvalue", "neg_score", 
    "pos_score", "neg_rank", "pos_rank", "num_guides")
  amg_numeric(mg, numeric_columns, "Selected MAGeCK contrast")
  for (column in c("neg_fdr", "pos_fdr", "neg_pvalue", "pos_pvalue")) if (any(mg[[column]] < 
    0 | mg[[column]] > 1, na.rm = TRUE)) 
    stop("Invalid probability: ", column)
  # 3. Call negative or positive selection from the corresponding MAGeCK FDR.
  # There is no additional fold-change filter. Results significant in both
  # directions are reported explicitly as ambiguous, not forced into one class.
  mg[, `:=`(calculated_direction, data.table::fcase(is.finite(neg_fdr) & neg_fdr < MAGECK_FDR_CUT & 
    is.finite(pos_fdr) & pos_fdr < MAGECK_FDR_CUT, "ambiguous_bidirectional", is.finite(neg_fdr) & 
    neg_fdr < MAGECK_FDR_CUT, "negative_selection", is.finite(pos_fdr) & pos_fdr < MAGECK_FDR_CUT, 
    "positive_selection", default = "not_significant"))]
  if ("selection_direction" %in% names(mg) && any(!is.na(mg$selection_direction) & mg$selection_direction != 
    mg$calculated_direction)) 
    stop("Stored MAGeCK directions disagree with the requested FDR threshold.")
  mg[, `:=`(selection_direction, calculated_direction)]
  # 4. Add/check gene names and preserve the RRA scores, ranks and guide counts
  # so the independent method's evidence remains available beside its calls.
  annotation <- amg_annotation(ANNOTATION_FILE)
  if (length(setdiff(mg$PAO1_ID, annotation$PAO1_ID))) 
    stop("MAGeCK loci are absent from the curated annotation.")
  if (!"gene_name" %in% names(mg)) 
    mg[, `:=`(gene_name, NA_character_)]
  missing_name <- is.na(mg$gene_name) | !nzchar(trimws(mg$gene_name))
  mg[missing_name, `:=`(gene_name, annotation$gene_name[match(PAO1_ID, annotation$PAO1_ID)])]
  columns <- c("PAO1_ID", "gene_name", "gene_lfc", "selection_direction", "neg_fdr", "pos_fdr", 
    "neg_pvalue", "pos_pvalue", "neg_score", "pos_score", "neg_rank", "pos_rank", "num_guides")
  all_genes <- mg[, ..columns]
  all_genes[, `:=`(minimum_fdr, pmin(neg_fdr, pos_fdr, na.rm = TRUE))]
  data.table::setorderv(all_genes, c("minimum_fdr", "gene_lfc", "PAO1_ID"), c(1L, 1L, 1L), na.last = TRUE)
  all_genes[, `:=`(minimum_fdr, NULL)]
  # 5. Keep both the full tested-gene table and its significant subset.
  # Export the significant Table S7 input and counts by direction with provenance.
  hits <- all_genes[selection_direction != "not_significant"]
  data.table::setorderv(hits, c("gene_lfc", "PAO1_ID"), na.last = TRUE)
  summary <- data.table::data.table(metric = c("n_MAGeCK_tested", "n_MAGeCK_significant", "n_negative_selection", 
    "n_positive_selection", "n_ambiguous_bidirectional"), value = c(nrow(all_genes), nrow(hits), 
    sum(hits$selection_direction == "negative_selection"), sum(hits$selection_direction == 
      "positive_selection"), sum(hits$selection_direction == "ambiguous_bidirectional")))
  settings <- amg_settings(c(MAGECK_FDR_CUT = MAGECK_FDR_CUT, LFC_filter = "none", selected_contrast = "B3 FBS 60 T1 versus None 0 T1"), 
    c(MAGeCK = found, annotation = ANNOTATION_FILE))
  amg_write(list(MAGeCK_FBS_vs_MH_all_genes = all_genes, MAGeCK_FBS_vs_MH_significant = hits, 
    Summary = summary, Analysis_settings = settings), OUTPUT_DIR, OVERWRITE)
  message("Prepared S7 MAGeCK sheet: ", nrow(hits), " significant genes. Saved to ", OUTPUT_DIR)
  invisible(hits)
}

main()

