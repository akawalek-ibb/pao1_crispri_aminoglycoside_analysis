# Test KEGG overrepresentation of recurrent sensitizing and protective genes/TUs.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/08_AMG_recurrent_KEGG_enrichment.R")
# Inputs:
# - preprocessing_scripts/out_03_AMG_relevant_results/AG_gene_recurrence.tsv and
#   Analysis_settings.tsv.
# - preprocessing_scripts/out_06_AMG_TU_functional_annotation/gene_functional_membership.tsv.gz,
#   TU_functional_membership.tsv.gz and gene_TU_map.tsv.
# Upstream / execution order:
# - preprocessing_scripts/03_AMG_relevant_results.R and
#   preprocessing_scripts/06_AMG_TU_functional_annotation.R.
# Outputs:
# - preprocessing_scripts/out_08_AMG_recurrent_KEGG_enrichment/ — KEGG_enrichment.tsv,
#   KEGG_background.tsv, KEGG_settings.tsv and Analysis_settings.tsv.
# Interpretation / operation:
# - Background: genes evaluated in ≥2 T1 drug comparisons and mapped to a TU; unique mapped TUs
#   form the TU background.
# - Retain terms with 5–500 background members at both levels. Use one-sided Fisher tests and BH
#   by direction and level.

GENE_DIR <- "preprocessing_scripts/out_03_AMG_relevant_results"

ANNOTATION_DIR <- "preprocessing_scripts/out_06_AMG_TU_functional_annotation"

OUTPUT_DIR <- "preprocessing_scripts/out_08_AMG_recurrent_KEGG_enrichment"

DATABASE <- "KEGG"

MIN_TERM_SIZE <- 5L

MAX_TERM_SIZE <- 500L

OVERWRITE <- TRUE

# For each pathway, ask whether recurrent hits occur more often than expected
# among the eligible background units. Use a one-sided Fisher test and apply
# BH across pathways within this direction/level. Units are either genes or TUs.
s8_ora <- function(membership, universe, query, terms, unit, prefix) {
  universe <- sort(unique(universe))
  query <- sort(unique(query))
  if (length(setdiff(query, universe))) 
    stop("A recurrent query unit is outside its evaluable background.")
  rows <- lapply(terms, function(term_value) {
    members <- sort(unique(membership[[unit]][membership$term == term_value & membership[[unit]] %in% 
      universe]))
    overlap <- intersect(members, query)
    a <- length(overlap)
    q <- length(query)
    n <- length(members)
    u <- length(universe)
    table <- matrix(c(a, q - a, n - a, u - q - n + a), nrow = 2, byrow = TRUE)
    if (any(table < 0)) 
      stop("Invalid enrichment contingency table for ", term_value)
    fit <- stats::fisher.test(table, alternative = "greater")
    expected <- q * n/u
    data.table::data.table(term = term_value, overlap_n = a, total_term_n = n, coverage_pct = 100 * 
      a/n, fold_enrichment = if (expected > 0) 
      a/expected
    else NA_real_, odds_ratio = unname(fit$estimate), p_value = fit$p.value, overlap_members = s8_collapse(overlap))
  })
  out <- data.table::rbindlist(rows)
  out[, `:=`(FDR_BH, stats::p.adjust(p_value, method = "BH"))]
  data.table::setcolorder(out, c("term", "overlap_n", "total_term_n", "coverage_pct", "fold_enrichment", 
    "odds_ratio", "p_value", "FDR_BH", "overlap_members"))
  data.table::setnames(out, setdiff(names(out), "term"), paste0(prefix, "_", setdiff(names(out), 
    "term")))
  out
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

s8_collapse <- function(x, sep = "; ") {
  x <- sort(unique(s8_clean(x)))
  x <- x[!is.na(x)]
  if (length(x)) 
    paste(x, collapse = sep)
  else NA_character_
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
  if (MIN_TERM_SIZE < 1 || MAX_TERM_SIZE < MIN_TERM_SIZE) 
    stop("Invalid KEGG term-size range.")
  genes <- s8_read(file.path(GENE_DIR, "AG_gene_recurrence.tsv"))
  map <- s8_read(file.path(ANNOTATION_DIR, "gene_TU_map.tsv"))
  gene_terms <- s8_read(file.path(ANNOTATION_DIR, "gene_functional_membership.tsv.gz"))[database == 
    DATABASE]
  tu_terms <- s8_read(file.path(ANNOTATION_DIR, "TU_functional_membership.tsv.gz"))[database == 
    DATABASE]
  # 1. Define the background before counting pathway hits.
  # Eligible genes were evaluated in enough T1 comparisons to satisfy recurrence
  # and have a TU assignment. Their unique TUs define the TU-level background.
  minimum_tested <- s8_setting(GENE_DIR, "RECURRENCE_MIN", numeric = TRUE)
  background <- genes[n_independent_T1_tested >= minimum_tested, .(PAO1_ID, n_independent_T1_tested)]
  background <- merge(background, map[, .(PAO1_ID, TU_ID)], by = "PAO1_ID", all.x = TRUE, sort = FALSE)
  background <- background[!is.na(TU_ID)]
  gene_universe <- sort(unique(background$PAO1_ID))
  tu_universe <- sort(unique(background$TU_ID))
  if (!length(gene_universe) || !length(tu_universe)) 
    stop("No evaluable background units for KEGG enrichment.")
  # 2. Calculate pathway sizes within these experimental backgrounds.
  # Retain pathways meeting the configured size limits at both gene and TU levels
  # so the two enrichment summaries compare the same collection of pathways.
  gene_sizes <- gene_terms[PAO1_ID %in% gene_universe, .(gene_n = data.table::uniqueN(PAO1_ID)), 
    by = term]
  tu_sizes <- tu_terms[TU_ID %in% tu_universe, .(TU_n = data.table::uniqueN(TU_ID)), by = term]
  sizes <- merge(gene_sizes, tu_sizes, by = "term")
  sizes <- sizes[gene_n >= MIN_TERM_SIZE & gene_n <= MAX_TERM_SIZE & TU_n >= MIN_TERM_SIZE & 
    TU_n <= MAX_TERM_SIZE]
  data.table::setorder(sizes, term)
  if (!nrow(sizes)) 
    stop("No KEGG terms pass both gene/TU size requirements.")
  all_results <- query_sizes <- vector("list", 2L)
  # 3. Analyse recurrent sensitizing and protective genes separately.
  # At TU level, the query is the unique TUs containing those recurrent genes,
  # not the significant pooled-guide TUs from task 05. This reduces repeated
  # counting of linked member genes while retaining the primary gene-hit definition.
  for (i in seq_along(c("Sensitizing", "Protective"))) {
    direction <- c("Sensitizing", "Protective")[i]
    flag <- paste0("recurrent_", tolower(direction))
    query_gene <- genes$PAO1_ID[genes[[flag]]]
    if (length(setdiff(query_gene, gene_universe))) 
      stop("Recurrent ", direction, " genes are missing from the canonical evaluable KEGG background.")
    query_tu <- unique(map$TU_ID[match(query_gene, map$PAO1_ID)])
    a <- s8_ora(gene_terms, gene_universe, query_gene, sizes$term, "PAO1_ID", "gene")
    b <- s8_ora(tu_terms, tu_universe, query_tu, sizes$term, "TU_ID", "TU")
    z <- merge(a, b, by = "term", sort = FALSE)
    z[, `:=`(direction, direction)]
    data.table::setcolorder(z, c("term", "direction", setdiff(names(z), c("term", "direction"))))
    data.table::setorderv(z, c("gene_p_value", "gene_fold_enrichment", "term"), c(1L, -1L, 
      1L))
    all_results[[i]] <- z
    query_sizes[[i]] <- data.table::data.table(direction = direction, gene_query_n = length(unique(query_gene)), 
      TU_query_n = length(query_tu), gene_universe_n = length(gene_universe), TU_universe_n = length(tu_universe), 
      n_terms_tested = nrow(sizes))
  }
  out <- data.table::rbindlist(all_results)
  data.table::setorder(background, PAO1_ID)
  # Record thresholds and background definitions beside the enrichment results;
  # these distinguish the measured background from the whole annotated genome.
  settings <- s8_settings(c(DATABASE = DATABASE, MIN_TERM_SIZE = MIN_TERM_SIZE, MAX_TERM_SIZE = MAX_TERM_SIZE, 
    background_minimum_evaluable_T1 = minimum_tested, BH_family = "Retained terms, separately by direction and gene/TU level", 
    TU_queries = "TUs containing primary recurrent gene hits; not pooled-guide TU calls"), 
    c(recurrence = file.path(GENE_DIR, "AG_gene_recurrence.tsv"), gene_membership = file.path(ANNOTATION_DIR, 
      "gene_functional_membership.tsv.gz"), TU_membership = file.path(ANNOTATION_DIR, "TU_functional_membership.tsv.gz"), 
      TU_map = file.path(ANNOTATION_DIR, "gene_TU_map.tsv"), call_settings = file.path(GENE_DIR, 
        "Analysis_settings.tsv")))
  s8_write(list(KEGG_enrichment = out, KEGG_background = background, KEGG_settings = data.table::rbindlist(query_sizes), 
    Analysis_settings = settings), OUTPUT_DIR, OVERWRITE)
  message("S8 KEGG: ", nrow(sizes), " terms in each direction; backgrounds ", length(gene_universe), 
    " genes / ", length(tu_universe), " TUs.")
  invisible(out)
}

main()

