# Test GO overrepresentation among up- and downregulated knockdown genes.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/19_AMG_RNAseq_GO_ORA.R")
# Inputs:
# - preprocessing_scripts/out_18_AMG_RNAseq_results/RNAseq_results.tsv.gz,
#   Comparison_manifest.tsv and Analysis_settings.tsv.
# - annotation_input/PAO1_GO_unified_propagated.tsv.gz: PAO1 locus, GO ID, term and ontology
#   memberships.
# Upstream / execution order:
# - preprocessing_scripts/18_AMG_RNAseq_results.R. Use the supplied GO reference; it is not
#   rebuilt by this script.
# Outputs:
# - preprocessing_scripts/out_19_AMG_RNAseq_GO_ORA/ — GO_ORA_results.tsv.gz,
#   GO_ORA_backgrounds.tsv and Analysis_settings.tsv.
# Interpretation / operation:
# - Test four knockdowns only; use one-sided hypergeometric tails and 5–500-member terms.
# - BH families are separate by comparison, direction and ontology, over returned
#   query-overlapping terms.

RNA_DIR <- "preprocessing_scripts/out_18_AMG_RNAseq_results"

GO_FILE <- "annotation_input/PAO1_GO_unified_propagated.tsv.gz"

OUTPUT_DIR <- "preprocessing_scripts/out_19_AMG_RNAseq_GO_ORA"

MIN_GS_SIZE <- 5L

MAX_GS_SIZE <- 500L

ONTOLOGIES <- c("BP", "MF", "CC")

OVERWRITE <- TRUE

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

rnaseq_setting <- function(directory, name, numeric = FALSE) {
  x <- amg_read(file.path(directory, "Analysis_settings.tsv"))
  value <- x$value[x$setting == name]
  if (length(value) != 1L || is.na(value)) 
    stop("Missing/ambiguous setting: ", name)
  if (numeric) 
    as.numeric(value)
  else as.character(value)
}

rnaseq_check_inputs <- function(directory) {
  x <- amg_read(file.path(directory, "Analysis_settings.tsv"))
  inputs <- x[grepl("^input:", setting)]
  for (i in seq_len(nrow(inputs))) {
    path <- as.character(inputs$value[i])
    if (!startsWith(path, "preprocessing_scripts/out_")) 
      next
    hash <- x$value[x$setting == sub("^input:", "md5:", inputs$setting[i])]
    if (!file.exists(path) || length(hash) != 1L || is.na(hash) || unname(tools::md5sum(path)) != 
      hash) 
      stop("Stale/missing intermediate for ", directory, ": ", path, ". Re-run its preprocessing task.")
  }
  invisible(TRUE)
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

# Read the supplied gene-to-GO memberships, which already include propagated
# ancestor terms. This function cleans memberships; it does not build the GO
# hierarchy. Exclude NOT annotations, obsolete terms and the three root terms,
# then retain unique gene/term memberships separately for BP, MF and CC.
rnaseq_go <- function(path) {
  if (!file.exists(path) && endsWith(path, ".gz") && file.exists(sub("\\.gz$", "", path))) 
    path <- sub("\\.gz$", "", path)
  raw <- rnaseq_read(path)
  gene <- rnaseq_column(raw, c("PAO1_ID", "Locus_Tag", "locus_tag"), path)
  accession <- rnaseq_column(raw, c("GO_ID", "GO ID", "GO", "go_id", "GO_accession", "Accession"), 
    path)
  term <- rnaseq_column(raw, c("GO_Term", "GO Term", "term"), path)
  namespace <- rnaseq_column(raw, c("Namespace", "namespace", "GO_namespace"), path)
  qualifier <- rnaseq_column(raw, c("Qualifier", "Qualifiers", "QUALIFIER"), path, required = FALSE)
  if (!is.na(qualifier)) 
    raw <- raw[is.na(get(qualifier)) | !grepl("(^|\\|)NOT(\\||$)", as.character(get(qualifier)))]
  x <- data.table::data.table(PAO1_ID = rnaseq_clean(raw[[gene]]), GO_ID = rnaseq_clean(raw[[accession]]), 
    GO_term = rnaseq_clean(raw[[term]]), Namespace = rnaseq_clean(raw[[namespace]]))
  x[, `:=`(ontology, c(biological_process = "BP", molecular_function = "MF", cellular_component = "CC")[Namespace])]
  x <- x[!is.na(PAO1_ID) & !is.na(GO_ID) & !is.na(GO_term) & !is.na(ontology)]
  if (any(!grepl("^GO:[0-9]{7}$", x$GO_ID))) 
    stop("GO annotation contains an invalid accession.")
  x <- x[!grepl("^obsolete ", GO_term) & !GO_ID %in% c("GO:0008150", "GO:0003674", "GO:0005575")]
  x <- unique(x[, .(PAO1_ID, GO_ID, GO_term, ontology)])
  check <- x[, .(names = data.table::uniqueN(GO_term), domains = data.table::uniqueN(ontology)), 
    by = GO_ID]
  if (!nrow(x) || any(check$names != 1L | check$domains != 1L)) 
    stop("GO accessions have inconsistent names/domains or annotation is empty.")
  attr(x, "source_file") <- path
  x
}

rnaseq_gene_names <- function(ids, de) {
  idx <- match(ids, de$PAO1_ID)
  out <- de$gene_name[idx]
  out[is.na(out) | !nzchar(out)] <- ids[is.na(out) | !nzchar(out)]
  out
}

# For one comparison/ontology, define the background as tested genes with an
# annotation in that ontology. Analyse Up and Down DEG lists separately.
# After the term-size limits, test only terms overlapping the query using the
# one-sided hypergeometric tail. Apply BH across those returned terms within
# each direction, before any later selection of terms for display.
rnaseq_ora_one <- function(de, annotation, min_size, max_size) {
  universe <- unique(de$PAO1_ID[de$tested_for_DE & de$PAO1_ID %in% annotation$PAO1_ID])
  sets <- annotation[PAO1_ID %in% universe, .(members = list(unique(PAO1_ID)), GO_term = GO_term[1L]), 
    by = GO_ID]
  N <- length(universe)
  answer <- sizes <- vector("list", 2L)
  for (i in seq_along(answer)) {
    direction <- c("Up", "Down")[i]
    query <- unique(de$PAO1_ID[de$tested_for_DE & de$DEG_direction == direction & de$PAO1_ID %in% 
      universe])
    n <- length(query)
    retained <- sets[lengths(members) >= min_size & lengths(members) <= max_size]
    retained[, `:=`(overlap, lapply(members, function(ids) query[query %in% ids]))]
    retained <- retained[lengths(overlap) > 0L]
    sizes[[i]] <- data.table::data.table(direction = direction, background_n = N, query_n = n, 
      n_terms_tested = nrow(retained))
    empty <- data.table::data.table(direction = character(), GO_ID = character(), GO_term = character(), 
      FDR = numeric(), P_value = numeric(), Count = integer(), GeneRatio = character(), BgRatio = character(), 
      FoldEnrichment = numeric(), overlap_ids = character(), overlap_names = character())
    if (!nrow(retained)) {
      answer[[i]] <- empty
      next
    }
    counts <- lengths(retained$overlap)
    m <- lengths(retained$members)
    p <- stats::phyper(counts - 1L, m, N - m, n, lower.tail = FALSE)
    answer[[i]] <- data.table::data.table(direction = direction, GO_ID = retained$GO_ID, GO_term = retained$GO_term, 
      FDR = stats::p.adjust(p, method = "BH"), P_value = p, Count = counts, GeneRatio = paste0(counts, 
        "/", n), BgRatio = paste0(m, "/", N), FoldEnrichment = (counts/m) * N/n, overlap_ids = vapply(retained$overlap, 
        paste, character(1), collapse = "; "), overlap_names = vapply(retained$overlap, 
        function(ids) paste(rnaseq_gene_names(ids, de), collapse = "; "), character(1)))
  }
  list(results = data.table::rbindlist(answer), sizes = data.table::rbindlist(sizes))
}

main <- function() {
  amg_require("data.table")
  if (MIN_GS_SIZE < 1L || MAX_GS_SIZE < MIN_GS_SIZE || !length(ONTOLOGIES) || any(!ONTOLOGIES %in% 
    c("BP", "MF", "CC"))) 
    stop("Invalid GO ORA settings.")
  # 1. Confirm task-18 provenance and read the complete expression results.
  # Use only comparisons flagged for GO; TOB is not in that analysis subset.
  rnaseq_check_inputs(RNA_DIR)
  de <- amg_read(file.path(RNA_DIR, "RNAseq_results.tsv.gz"))
  manifest <- amg_read(file.path(RNA_DIR, "Comparison_manifest.tsv"))[use_for_GO == TRUE]
  annotation <- rnaseq_go(GO_FILE)
  go_path <- attr(annotation, "source_file")
  result <- sizes <- list()
  k <- 0L
  # 2. Repeat overrepresentation testing for every ontology and comparison.
  # Each call uses its own tested/annotated background and the Up/Down labels
  # already assigned in task 18. No count normalization or DE testing is repeated.
  for (ontology_i in ONTOLOGIES) for (id in manifest$comparison_id) {
    k <- k + 1L
    z <- de[comparison_id == id]
    data.table::setorder(z, source_row)
    x <- rnaseq_ora_one(z, annotation[ontology == ontology_i], MIN_GS_SIZE, MAX_GS_SIZE)
    x$results[, `:=`(comparison_id = id, ontology = ontology_i)]
    x$sizes[, `:=`(comparison_id = id, ontology = ontology_i)]
    result[[k]] <- x$results
    sizes[[k]] <- x$sizes
  }
  # 3. Combine term-level results and record background/query sizes separately.
  # Terms with no query overlap are not returned; absence of a row is not an
  # estimated FDR of 1. Gene IDs and names expose each observed overlap.
  results <- data.table::rbindlist(result)
  background <- data.table::rbindlist(sizes)
  data.table::setorderv(results, c("comparison_id", "direction", "ontology", "FDR", "P_value", 
    "GO_ID"))
  display <- results[, .(Comparison = comparison_id, `DEG direction` = direction, `GO domain` = ontology, 
    `GO ID` = GO_ID, `GO term` = GO_term, `FDR (BH)` = FDR, `P value` = P_value, Count, BgRatio, 
    FoldEnrichment, `Overlapping gene IDs` = overlap_ids, `Overlapping gene names` = overlap_names)]
  settings <- amg_settings(c(MIN_GS_SIZE = MIN_GS_SIZE, MAX_GS_SIZE = MAX_GS_SIZE, ONTOLOGIES = paste(ONTOLOGIES, 
    collapse = ";"), method = "One-sided hypergeometric tail; same query-overlapping term family as clusterProfiler::enricher", 
    adjustment = "BH separately per comparison/direction/ontology", zero_overlap = "not returned", 
    DE_FDR_CUT = rnaseq_setting(RNA_DIR, "DE_FDR_CUT"), DE_ABS_LOG2FC_CUT = rnaseq_setting(RNA_DIR, 
      "DE_ABS_LOG2FC_CUT")), c(results = file.path(RNA_DIR, "RNAseq_results.tsv.gz"), manifest = file.path(RNA_DIR, 
    "Comparison_manifest.tsv"), RNA_settings = file.path(RNA_DIR, "Analysis_settings.tsv"), 
    annotation = go_path))
  # 4. Export the full returned ORA results and the test-family/background
  # descriptions for Table S9 and the selected-term analysis in task 22.
  rnaseq_write(list(GO_ORA_results = results, GO_ORA_backgrounds = background, Analysis_settings = settings), 
    OUTPUT_DIR, OVERWRITE, gzip = "GO_ORA_results")
  message("GO ORA: ", nrow(display), " returned term/comparison/direction records.")
  invisible(display)
}

main()

