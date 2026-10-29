# Group primary recurrent gene hits by their annotated transcriptional units.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/07_AMG_recurrent_gene_TU_collapse.R")
# Inputs:
# - preprocessing_scripts/out_03_AMG_relevant_results/AG_gene_recurrence.tsv and
#   Comparison_manifest.tsv.
# - preprocessing_scripts/out_06_AMG_TU_functional_annotation/gene_TU_map.tsv.
# Upstream / execution order:
# - preprocessing_scripts/03_AMG_relevant_results.R and
#   preprocessing_scripts/06_AMG_TU_functional_annotation.R; these independent inputs must exist
#   before this script.
# Outputs:
# - preprocessing_scripts/out_07_AMG_recurrent_gene_TU_collapse/ — Gene_level_hits_by_TU.tsv,
#   recurrent_gene_TU_membership.tsv and settings.
# Interpretation / operation:
# - This grouping does not use pooled-guide TU significance or alter gene calls.

GENE_DIR <- "preprocessing_scripts/out_03_AMG_relevant_results"

ANNOTATION_DIR <- "preprocessing_scripts/out_06_AMG_TU_functional_annotation"

OUTPUT_DIR <- "preprocessing_scripts/out_07_AMG_recurrent_gene_TU_collapse"

OVERWRITE <- TRUE

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

s8_collapse <- function(x, sep = "; ") {
  x <- sort(unique(s8_clean(x)))
  x <- x[!is.na(x)]
  if (length(x)) 
    paste(x, collapse = sep)
  else NA_character_
}

s8_conditions <- function(x, order) {
  x <- order[order %in% unique(as.character(x))]
  if (length(x)) 
    paste(x, collapse = "; ")
  else NA_character_
}

s8_settings <- function(settings, inputs = character()) {
  data.table::data.table(setting = c(names(settings), paste0("input:", names(inputs)), paste0("md5:", 
    names(inputs)), "R_version"), value = c(as.character(settings), unname(inputs), unname(tools::md5sum(inputs)), 
    R.version.string))
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
  # 1. Read the primary recurrent gene calls and the canonical TU map.
  # This groups existing gene hits; it does not repeat the pooled-guide TU test
  # in task 05 and does not recalculate gene-level significance.
  genes <- s8_read(file.path(GENE_DIR, "AG_gene_recurrence.tsv"))
  manifest <- s8_read(file.path(GENE_DIR, "Comparison_manifest.tsv"))[count_for_recurrence == 
    TRUE]
  map <- s8_read(file.path(ANNOTATION_DIR, "gene_TU_map.tsv"))
  # 2. Retain genes recurrent in either direction and require a valid TU for each.
  # For each gene, collect the T1 treatments supporting each direction.
  recurrent <- genes[recurrent_sensitizing == TRUE | recurrent_protective == TRUE]
  idx <- match(recurrent$PAO1_ID, map$PAO1_ID)
  if (anyNA(idx) || anyNA(map$TU_ID[idx])) 
    stop("A recurrent gene has no canonical TU assignment.")
  recurrent[, `:=`(TU_ID, map$TU_ID[idx])]
  for (direction in c("Sensitizing", "Protective")) {
    strings <- vapply(seq_len(nrow(recurrent)), function(i) {
      values <- vapply(manifest$prefix, function(prefix) recurrent[[paste0(prefix, "_call")]][i], 
        character(1))
      s8_conditions(manifest$display_label[values == direction], manifest$display_label)
    }, character(1))
    recurrent[, `:=`((paste0(tolower(direction), "_conditions")), strings)]
  }
  all_members <- map[!is.na(TU_ID), .(n_all_genes_in_TU = data.table::uniqueN(PAO1_ID)), by = TU_ID]
  collapse_conditions <- function(x) s8_conditions(unlist(strsplit(x[!is.na(x)], "\\s*;\\s*")), 
    manifest$display_label)
  # 3. Collapse recurrent member genes into one row per TU, retaining gene names,
  # IDs and supporting conditions. Count sensitizing and protective members
  # separately, and report TUs containing evidence in both directions.
  collapsed <- recurrent[, .(n_recurrent_genes = data.table::uniqueN(PAO1_ID), recurrent_gene_PAO1_IDs = s8_collapse(PAO1_ID), 
    recurrent_gene_names = s8_collapse(gene_name), n_recurrent_sensitizing_genes = sum(recurrent_sensitizing), 
    recurrent_sensitizing_gene_PAO1_IDs = s8_collapse(PAO1_ID[recurrent_sensitizing]), recurrent_sensitizing_gene_names = s8_collapse(gene_name[recurrent_sensitizing]), 
    n_recurrent_protective_genes = sum(recurrent_protective), recurrent_protective_gene_PAO1_IDs = s8_collapse(PAO1_ID[recurrent_protective]), 
    recurrent_protective_gene_names = s8_collapse(gene_name[recurrent_protective]), TU_contains_recurrent_sensitizing = any(recurrent_sensitizing), 
    TU_contains_recurrent_protective = any(recurrent_protective), protective_conditions_union = collapse_conditions(protective_conditions[recurrent_protective]), 
    sensitizing_conditions_union = collapse_conditions(sensitizing_conditions[recurrent_sensitizing])), 
    by = TU_ID]
  collapsed[, `:=`(TU_recurrent_direction, data.table::fcase(TU_contains_recurrent_sensitizing & 
    TU_contains_recurrent_protective, "Both directions", TU_contains_recurrent_sensitizing, 
    "Sensitizing", TU_contains_recurrent_protective, "Protective", default = "None"))]
  collapsed <- merge(collapsed, all_members, by = "TU_ID", all.x = TRUE, sort = FALSE)
  columns <- c("TU_ID", "n_all_genes_in_TU", "n_recurrent_genes", "recurrent_gene_PAO1_IDs", 
    "recurrent_gene_names", "n_recurrent_sensitizing_genes", "recurrent_sensitizing_gene_PAO1_IDs", 
    "recurrent_sensitizing_gene_names", "n_recurrent_protective_genes", "recurrent_protective_gene_PAO1_IDs", 
    "recurrent_protective_gene_names", "TU_contains_recurrent_sensitizing", "TU_contains_recurrent_protective", 
    "TU_recurrent_direction", "protective_conditions_union", "sensitizing_conditions_union")
  collapsed <- collapsed[, ..columns]
  data.table::setorder(collapsed, -TU_contains_recurrent_sensitizing, -TU_contains_recurrent_protective, 
    -n_recurrent_genes, TU_ID)
  member_columns <- c("TU_ID", "PAO1_ID", "gene_name", "recurrent_sensitizing", "recurrent_protective", 
    "sensitizing_conditions", "protective_conditions")
  # 4. Also retain a gene-level membership table so each collapsed TU row can
  # be traced back to its individual recurrent genes.
  membership <- recurrent[, ..member_columns]
  data.table::setorder(membership, TU_ID, PAO1_ID)
  settings <- s8_settings(c(definition = "Collapse of primary recurrent genes, independent of pooled-guide TU testing"), 
    c(recurrence = file.path(GENE_DIR, "AG_gene_recurrence.tsv"), TU_map = file.path(ANNOTATION_DIR, 
      "gene_TU_map.tsv"), manifest = file.path(GENE_DIR, "Comparison_manifest.tsv")))
  # 5. Export the grouped summary and its individual-gene membership evidence.
  s8_write(list(Gene_level_hits_by_TU = collapsed, recurrent_gene_TU_membership = membership, 
    Analysis_settings = settings), OUTPUT_DIR, OVERWRITE)
  message("S8 gene-hit collapse: ", nrow(collapsed), " TUs; ", sum(collapsed$TU_contains_recurrent_sensitizing), 
    " sensitizing / ", sum(collapsed$TU_contains_recurrent_protective), " protective.")
  invisible(collapsed)
}

main()

