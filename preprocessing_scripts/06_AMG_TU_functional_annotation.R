# Prepare gene and TU functional memberships for recurrent-hit enrichment.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/06_AMG_TU_functional_annotation.R")
# Inputs:
# - annotation_input/PAO1_curated.tsv: unique PAO1_ID and TU.
# - annotation_input/PAO1_pathways_all_withCOG.tsv: gene memberships with source database and term.
# Upstream / execution order:
# - No upstream analysis script is required; use the supplied annotation/reference files.
# Outputs:
# - preprocessing_scripts/out_06_AMG_TU_functional_annotation/ —
#   gene_functional_membership.tsv.gz, TU_functional_membership.tsv.gz, gene_TU_map.tsv and
#   settings.
# Interpretation / operation:
# - KEGG and eggNOG_KEGG are separate databases. This script performs no statistical tests.

ANNOTATION_FILE <- "annotation_input/PAO1_curated.tsv"

PATHWAY_FILE <- "annotation_input/PAO1_pathways_all_withCOG.tsv"

OUTPUT_DIR <- "preprocessing_scripts/out_06_AMG_TU_functional_annotation"

DATABASES <- c("KEGG", "eggNOG_KEGG", "COG", "PseudoCAP")

OVERWRITE <- TRUE

s8_clean_term <- function(x) {
  x <- gsub("<[^>]+>", "", as.character(x))
  for (name in c("beta", "gamma", "alpha")) x <- gsub(paste0("&", name, ";"), name, x, fixed = TRUE)
  x <- gsub("&amp;", "&", x, fixed = TRUE)
  s8_clean(gsub("[[:space:]]+", " ", x))
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
  # 1. Read gene-to-TU assignments and the functional annotation catalogue.
  # Recognize the supported column-name variants without changing locus IDs.
  annotation <- s8_annotation(ANNOTATION_FILE)
  raw <- s8_read(PATHWAY_FILE)
  pick <- function(candidates) {
    found <- intersect(candidates, names(raw))
    if (!length(found)) 
      stop(PATHWAY_FILE, " lacks ", paste(candidates, collapse = " or "))
    found[1L]
  }
  id <- pick(c("Locus_Tag", "Locus Tag"))
  db <- pick(c("Pathway_Database", "Pathway Database"))
  term <- pick(c("Pathway_Name", "Pathway Name"))
  # 2. Keep unique gene/database/term memberships for the requested databases.
  # One gene may legitimately belong to several terms; this is not a hit analysis.
  # Terms are selected by database label, not by exclusions based on their names.
  genes <- data.table::data.table(PAO1_ID = s8_clean(raw[[id]]), database = s8_clean(raw[[db]]), 
    term = s8_clean_term(raw[[term]]))
  genes <- unique(genes[database %in% DATABASES & !is.na(PAO1_ID) & !is.na(term)])
  map <- annotation[, .(PAO1_ID, gene_name, TU_ID = TU)]
  # 3. Give each TU the terms assigned to any annotated member gene.
  # Count each TU only once per term, even if several member genes share it.
  tus <- merge(genes, map[, .(PAO1_ID, TU_ID)], by = "PAO1_ID", all.x = TRUE, sort = FALSE)
  tus <- unique(tus[!is.na(TU_ID), .(TU_ID, database, term)])
  data.table::setorder(genes, database, term, PAO1_ID)
  data.table::setorder(tus, database, term, TU_ID)
  data.table::setorder(map, PAO1_ID)
  settings <- s8_settings(c(DATABASES = paste(DATABASES, collapse = "; "), KEGG_filter = "database label only; no term-name exclusions", 
    TU_membership = "term inherited from any annotated member gene"), c(annotation = ANNOTATION_FILE, 
    pathways = PATHWAY_FILE))
  # 4. Export gene memberships, inherited TU memberships and the gene-to-TU map.
  # Downstream enrichment scripts use these as reference annotations.
  s8_write(list(gene_functional_membership = genes, TU_functional_membership = tus, gene_TU_map = map, 
    Analysis_settings = settings), OUTPUT_DIR, OVERWRITE, gzip = c("gene_functional_membership", 
    "TU_functional_membership"))
  message("Prepared gene/TU functional annotation.")
  invisible(list(genes = genes, TUs = tus, map = map))
}

main()

