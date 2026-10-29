# Prepare selected GO enrichment terms for Figure 6D and Table S9.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/22_AMG_RNAseq_selected_GO.R")
# Inputs:
# - preprocessing_scripts/out_18_AMG_RNAseq_results/RNAseq_results.tsv.gz and manifest/settings;
#   preprocessing_scripts/out_19_AMG_RNAseq_GO_ORA/GO_ORA_results.tsv.gz and settings.
# - preprocessing_scripts/input_AMG_RNAseq/selected_GO_terms.tsv: GO IDs/domain and display order.
# Upstream / execution order:
# - preprocessing_scripts/18_AMG_RNAseq_results.R → preprocessing_scripts/19_AMG_RNAseq_GO_ORA.R
#   → this script.
# Outputs:
# - preprocessing_scripts/out_22_AMG_RNAseq_selected_GO/ — selected_GO_results.tsv,
#   Selected_GO_backgrounds.tsv and Analysis_settings.tsv.
# Interpretation / operation:
# - Use a common tested/annotated background across knockdowns within each ontology.
# - If that differs from script 19, retest all eligible terms before selecting displayed terms,
#   preserving the BH family.

RNA_DIR <- "preprocessing_scripts/out_18_AMG_RNAseq_results"

ORA_DIR <- "preprocessing_scripts/out_19_AMG_RNAseq_GO_ORA"

SELECTED_TERMS_FILE <- "preprocessing_scripts/input_AMG_RNAseq/selected_GO_terms.tsv"

OUTPUT_DIR <- "preprocessing_scripts/out_22_AMG_RNAseq_selected_GO"

DISPLAY_FDR_CUT <- 0.05

DISPLAY_SCORE_CAP <- 12

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
  if (DISPLAY_FDR_CUT <= 0 || DISPLAY_FDR_CUT > 1 || DISPLAY_SCORE_CAP <= 0) 
    stop("Invalid selected-GO display settings.")
  rnaseq_check_inputs(ORA_DIR)
  # 1. Read the explicit display-term list and order. This selects which GO
  # terms to show; it does not define a smaller multiple-testing correction family.
  selected <- rnaseq_read(SELECTED_TERMS_FILE)
  amg_required(selected, c("display_order", "GO_ID", "ontology", "fallback_description"), SELECTED_TERMS_FILE)
  amg_unique(selected, "GO_ID", SELECTED_TERMS_FILE)
  amg_numeric(selected, "display_order", SELECTED_TERMS_FILE)
  if (anyNA(selected$display_order) || anyDuplicated(selected$display_order) || any(!selected$ontology %in% 
    c("BP", "MF", "CC")) || any(!grepl("^GO:[0-9]{7}$", selected$GO_ID))) 
    stop("Invalid selected GO terms/order.")
  data.table::setorder(selected, display_order)
  # 2. Read the complete expression and ORA results, then verify they use the
  # same inputs/annotation. This prevents mixing enrichments from different runs.
  de <- amg_read(file.path(RNA_DIR, "RNAseq_results.tsv.gz"))
  manifest <- amg_read(file.path(RNA_DIR, "Comparison_manifest.tsv"))[use_for_GO == TRUE]
  ora <- amg_read(file.path(ORA_DIR, "GO_ORA_results.tsv.gz"))
  annotation_path <- rnaseq_setting(ORA_DIR, "input:annotation")
  check_paths <- c(results = file.path(RNA_DIR, "RNAseq_results.tsv.gz"), manifest = file.path(RNA_DIR, 
    "Comparison_manifest.tsv"), annotation = annotation_path)
  for (name in names(check_paths)) if (!file.exists(check_paths[[name]]) || unname(tools::md5sum(check_paths[[name]])) != 
    rnaseq_setting(ORA_DIR, paste0("md5:", name))) 
    stop("Selected GO inputs differ from 19: ", name, ". Re-run with the same RNA-seq results and GO annotation.")
  annotation <- rnaseq_go(annotation_path)
  min_size <- rnaseq_setting(ORA_DIR, "MIN_GS_SIZE", numeric = TRUE)
  max_size <- rnaseq_setting(ORA_DIR, "MAX_GS_SIZE", numeric = TRUE)
  sources <- backgrounds <- list()
  # 3. Define a common tested/annotated background across the four knockdowns
  # within each ontology. If the individual backgrounds already match, reuse
  # task-19 ORA. Otherwise repeat ORA on the common intersection, testing the
  # eligible overlapping GO terms before selecting the display terms.
  for (ontology_i in unique(selected$ontology)) {
    ann <- annotation[ontology == ontology_i]
    universes <- lapply(manifest$comparison_id, function(id) unique(de$PAO1_ID[de$comparison_id == 
      id & de$tested_for_DE & de$PAO1_ID %in% ann$PAO1_ID]))
    common <- Reduce(intersect, universes)
    same <- all(vapply(universes, function(v) setequal(v, common), logical(1)))
    if (!length(common)) 
      stop("Empty common GO background for ", ontology_i)
    if (same) {
      sources[[ontology_i]] <- ora[ontology == ontology_i]
    } else {
      sources[[ontology_i]] <- data.table::rbindlist(lapply(manifest$comparison_id, function(id) {
        z <- de[comparison_id == id & PAO1_ID %in% common]
        data.table::setorder(z, source_row)
        x <- rnaseq_ora_one(z, ann, min_size, max_size)$results
        x[, `:=`(comparison_id = id, ontology = ontology_i)]
        x
      }))
    }
    backgrounds[[ontology_i]] <- data.table::data.table(ontology = ontology_i, n_common_background_genes = length(common), 
      reused_individual_ORA = same)
  }
  source <- data.table::rbindlist(sources, use.names = TRUE)
  amg_unique(source, c("comparison_id", "ontology", "GO_ID", "direction"), "Selected-GO ORA source")
  result <- list()
  k <- 0L
  # 4. Build every selected-term/comparison row and retain Up and Down evidence
  # separately. A term not returned by ORA remains missing/Not returned rather
  # than being assigned a fabricated p-value or a nonsignificant test result.
  for (i in seq_len(nrow(selected))) for (id in manifest$comparison_id) {
    k <- k + 1L
    term <- selected$GO_ID[i]
    ontology_i <- selected$ontology[i]
    found <- annotation[GO_ID == term & ontology == ontology_i, unique(GO_term)]
    if (!length(found)) 
      found <- selected$fallback_description[i]
    x <- data.table::data.table(comparison_id = id, ontology = ontology_i, GO_ID = term, GO_term = found[1L], 
      display_order = selected$display_order[i])
    for (direction_i in c("Up", "Down")) {
      z <- source[comparison_id == id & ontology == ontology_i & GO_ID == term & direction == 
        direction_i]
      if (nrow(z) > 1L) 
        stop("Duplicate selected-term directional results.")
      returned <- nrow(z) == 1L && is.finite(z$FDR)
      x[, `:=`((paste0(direction_i, c("_Count", "_GeneRatio", "_FDR", "_result", "_gene_IDs", 
        "_gene_names"))), list(if (returned) z$Count else NA_integer_, if (returned) z$GeneRatio else NA_character_, 
        if (returned) z$FDR else NA_real_, if (!returned) "Not returned" else if (z$FDR < 
         DISPLAY_FDR_CUT) "Significant" else "Returned, not significant", if (returned) z$overlap_ids else NA_character_, 
        if (returned) z$overlap_names else NA_character_))]
    }
    up <- is.finite(x$Up_FDR)
    down <- is.finite(x$Down_FDR)
    # 5. For a single display score, choose the direction with the smaller FDR
    # (Up wins a tie). Significant Up enrichment receives positive -log10(FDR);
    # Down receives a negative score. Returned nonsignificant results score zero;
    # missing results remain NA. Keep a flag if both directions are significant.
    prefer_up <- up && (!down || x$Up_FDR <= x$Down_FDR)
    dominant <- if (prefer_up) 
      x$Up_FDR
    else x$Down_FDR
    direction <- if (!up && !down) 
      "Not returned"
    else if (!is.finite(dominant) || dominant >= DISPLAY_FDR_CUT) 
      "Not significant"
    else if (prefer_up) 
      "Up"
    else "Down"
    score <- if (direction == "Up") 
      -log10(max(dominant, .Machine$double.xmin))
    else if (direction == "Down") 
      log10(max(dominant, .Machine$double.xmin))
    else if (direction == "Not significant") 
      0
    else NA_real_
    # Cap only the displayed score so extreme enrichment does not dominate the
    # colour scale. Retain the uncapped score and both directional FDRs as evidence.
    x[, `:=`(displayed_direction = direction, signed_log10_FDR = score, displayed_signed_log10_FDR = pmin(DISPLAY_SCORE_CAP, 
      pmax(-DISPLAY_SCORE_CAP, score)), both_directions_significant = up & down & Up_FDR < 
      DISPLAY_FDR_CUT & Down_FDR < DISPLAY_FDR_CUT)]
    result[[k]] <- x
  }
  results <- data.table::rbindlist(result)
  results[, `:=`(plot_order, match(comparison_id, manifest$comparison_id))]
  data.table::setorder(results, display_order, plot_order)
  inputs <- c(ORA = file.path(ORA_DIR, "GO_ORA_results.tsv.gz"), ORA_settings = file.path(ORA_DIR, 
    "Analysis_settings.tsv"), RNA = file.path(RNA_DIR, "RNAseq_results.tsv.gz"), manifest = file.path(RNA_DIR, 
    "Comparison_manifest.tsv"), selected_terms = SELECTED_TERMS_FILE, annotation = attr(annotation, 
    "source_file"))
  settings <- amg_settings(c(DISPLAY_FDR_CUT = DISPLAY_FDR_CUT, DISPLAY_SCORE_CAP = DISPLAY_SCORE_CAP, 
    n_selected_terms = nrow(selected), n_comparisons = nrow(manifest), universe = "intersection of tested annotated genes across the four knockdowns", 
    BH_selection = "all eligible overlapping GO terms before editorial term selection", tie_rule = "Up when directional FDR values tie"), 
    inputs)
  # 6. Export the ordered display data and common-background audit for Figure 6D
  # and Table S9. This is an ORA-based display, separate from GO GSEA in task 20.
  rnaseq_write(list(selected_GO_results = results, Selected_GO_backgrounds = data.table::rbindlist(backgrounds), 
    Analysis_settings = settings), OUTPUT_DIR, OVERWRITE)
  message("Selected GO terms: ", nrow(results), " term/comparison records for Fig. 6D.")
  invisible(results)
}

main()

