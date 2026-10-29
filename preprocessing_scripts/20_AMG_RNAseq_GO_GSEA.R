# Run ranked GO enrichment for the four CRISPRi knockdowns.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/20_AMG_RNAseq_GO_GSEA.R")
# Inputs:
# - preprocessing_scripts/out_18_AMG_RNAseq_results/RNAseq_results.tsv.gz,
#   Comparison_manifest.tsv and Analysis_settings.tsv.
# - annotation_input/PAO1_GO_unified_propagated.tsv.gz: gene–GO memberships and BP/MF/CC labels.
# Upstream / execution order:
# - preprocessing_scripts/18_AMG_RNAseq_results.R. Script 19 is a separate analysis and is not a
#   prerequisite.
# Outputs:
# - preprocessing_scripts/out_20_AMG_RNAseq_GO_GSEA/ — GO_GSEA_results.tsv.gz,
#   GO_GSEA_ranks.tsv.gz, GO_GSEA_diagnostics.tsv and settings.
# Interpretation / operation:
# - Requires data.table, fgsea and BiocParallel. Rank tested annotated genes by sign(logFC) ×
#   sqrt(LR).
# - Use the configured fgsea multilevel procedure and 5–500-member sets; BH is per
#   comparison/ontology.

RNA_DIR <- "preprocessing_scripts/out_18_AMG_RNAseq_results"

GO_FILE <- "annotation_input/PAO1_GO_unified_propagated.tsv.gz"

OUTPUT_DIR <- "preprocessing_scripts/out_20_AMG_RNAseq_GO_GSEA"

MIN_GS_SIZE <- 5L

MAX_GS_SIZE <- 500L

ONTOLOGIES <- c("BP", "MF", "CC")

GSEA_RANDOM_SEED <- 20261007L

RESET_AS_DOSE_WRAPPER <- TRUE

GSEA_EXPONENT <- 1

GSEA_EPS <- 1e-10

GSEA_NPERM_SIMPLE <- 1000L

GSEA_SAMPLE_SIZE <- 101L

GSEA_WORKERS <- 0L

OVERWRITE <- TRUE

# Recover the pathway members contributing to the dominant excursion of the
# ranked enrichment curve. This uses the implemented hit-slice convention;
# the exported IDs preserve its ordering for reproducibility.
rnaseq_leading_edge <- function(scores, members, exponent = 1) {
  hits <- names(scores) %in% members
  if (!any(hits) || all(hits)) 
    return(character())
  weights <- abs(scores)^exponent * hits
  if (!is.finite(sum(weights)) || sum(weights) <= 0) 
    return(character())
  run <- cumsum(weights/sum(weights)) - cumsum(!hits)/sum(!hits)
  hit_ids <- names(scores)[hits]
  hit_run <- run[hits]
  if (abs(max(run)) > abs(min(run))) 
    return(hit_ids[seq_len(which.max(hit_run))])
  i <- which.min(hit_run)
  hit_ids[-seq.int(1L, i - 1L)]
}

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

main <- function() {
  amg_require(c("data.table", "fgsea"))
  if (MIN_GS_SIZE < 1L || MAX_GS_SIZE < MIN_GS_SIZE || GSEA_EPS < 0 || GSEA_NPERM_SIMPLE < 1L || 
    GSEA_SAMPLE_SIZE < 3L || GSEA_WORKERS < 0L || GSEA_EXPONENT <= 0 || any(!ONTOLOGIES %in% 
    c("BP", "MF", "CC"))) 
    stop("Invalid GSEA settings.")
  # 1. Read complete tested expression results and the propagated GO memberships.
  # GSEA uses all eligible ranked genes, not just significant DEGs. Only the
  # comparisons flagged for GO in task 18 enter this analysis.
  rnaseq_check_inputs(RNA_DIR)
  de <- amg_read(file.path(RNA_DIR, "RNAseq_results.tsv.gz"))
  manifest <- amg_read(file.path(RNA_DIR, "Comparison_manifest.tsv"))[use_for_GO == TRUE]
  annotation <- rnaseq_go(GO_FILE)
  go_path <- attr(annotation, "source_file")
  # 2. Set the configured random-number generator and seed for enrichment.
  # Restore the caller's previous random-number state when the script finishes.
  previous_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) 
    previous_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(previous_kind))
    if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv) else if (exists(".Random.seed", 
      envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(GSEA_RANDOM_SEED)
  results <- ranks <- diagnostics <- list()
  k <- 0L
  # 3. For each comparison/ontology, rank tested annotated genes by the signed
  # likelihood-ratio score from task 18. Highest scores favour upregulation;
  # lowest scores favour downregulation. Locus IDs break exact ties in sorting.
  for (ontology_i in ONTOLOGIES) for (id in manifest$comparison_id) {
    k <- k + 1L
    ann <- annotation[ontology == ontology_i]
    z <- de[comparison_id == id & tested_for_DE & PAO1_ID %in% ann$PAO1_ID]
    if (!nrow(z) || any(!is.finite(z$rank_score))) 
      stop("Incomplete signed LR ranks for ", id, " / ", ontology_i, ". Do not silently substitute logFC ranks.")
    data.table::setorderv(z, c("rank_score", "PAO1_ID"), c(-1L, 1L))
    scores <- stats::setNames(z$rank_score, z$PAO1_ID)
    sets <- split(ann$PAO1_ID, ann$GO_ID)
    ranks[[k]] <- z[, .(comparison_id, ontology = ontology_i, PAO1_ID, gene_name, logFC, LR, 
      rank_score)]
    if (RESET_AS_DOSE_WRAPPER) 
      set.seed(get(".Random.seed", envir = .GlobalEnv))
    warnings <- character()
    fit <- NULL
    # 4. Test whether each eligible GO set is concentrated near either end of
    # that ranking using fgseaMultilevel. Keep warnings as diagnostic information;
    # use the configured set sizes, enrichment exponent and precision settings.
    invisible(utils::capture.output(fit <- withCallingHandlers(fgsea::fgseaMultilevel(pathways = sets, 
      stats = scores, minSize = MIN_GS_SIZE, maxSize = MAX_GS_SIZE, eps = GSEA_EPS, gseaParam = GSEA_EXPONENT, 
      nPermSimple = GSEA_NPERM_SIMPLE, sampleSize = GSEA_SAMPLE_SIZE, scoreType = "std", 
      nproc = GSEA_WORKERS), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })))
    fit <- data.table::as.data.table(fit)
    # 5. Apply BH separately within this comparison/ontology's returned tests.
    # Positive NES indicates enrichment toward upregulated genes; negative NES
    # indicates enrichment toward downregulated genes. Preserve missing-test counts
    # in diagnostics before dropping nonfinite results from the final term table.
    fit[, `:=`(FDR, stats::p.adjust(pval, method = "BH"))]
    diagnostics[[k]] <- data.table::data.table(comparison_id = id, ontology = ontology_i, n_ranked_genes = length(scores), 
      n_terms_returned = nrow(fit), n_missing_p = sum(!is.finite(fit$pval)), n_tied_rank_values = sum(duplicated(scores)), 
      warnings = if (length(warnings)) 
        paste(unique(warnings), collapse = "; ")
      else NA_character_)
    fit <- fit[is.finite(pval) & is.finite(FDR)]
    # 6. Add contributing leading-edge gene IDs/names to each returned GO term.
    # Keep the complete ranking separately so the enrichment inputs are inspectable.
    ids <- lapply(fit$pathway, function(term) rnaseq_leading_edge(scores, sets[[term]], GSEA_EXPONENT))
    result <- data.table::data.table(comparison_id = id, ontology = ontology_i, GO_ID = fit$pathway, 
      GO_term = ann$GO_term[match(fit$pathway, ann$GO_ID)], NES = fit$NES, FDR = fit$FDR, 
      P_value = fit$pval, set_size = fit$size, ES = fit$ES, log2err = if ("log2err" %in% 
        names(fit)) 
        fit$log2err
      else NA_real_, leading_ids = vapply(ids, paste, character(1), collapse = "; "), leading_names = vapply(ids, 
        function(v) paste(rnaseq_gene_names(v, z), collapse = "; "), character(1)))
    data.table::setorder(result, FDR, -NES, GO_ID)
    results[[k]] <- result
    message("GO GSEA: ", id, " / ", ontology_i, "; ", nrow(result), " terms.")
  }
  all_results <- data.table::rbindlist(results)
  settings <- amg_settings(c(MIN_GS_SIZE = MIN_GS_SIZE, MAX_GS_SIZE = MAX_GS_SIZE, ONTOLOGIES = paste(ONTOLOGIES, 
    collapse = ";"), GSEA_RANDOM_SEED = GSEA_RANDOM_SEED, RESET_AS_DOSE_WRAPPER = RESET_AS_DOSE_WRAPPER, 
    GSEA_EXPONENT = GSEA_EXPONENT, GSEA_EPS = GSEA_EPS, GSEA_NPERM_SIMPLE = GSEA_NPERM_SIMPLE, 
    GSEA_SAMPLE_SIZE = GSEA_SAMPLE_SIZE, GSEA_WORKERS = GSEA_WORKERS, workers_used = if (GSEA_WORKERS == 
      0L) BiocParallel::bpnworkers(BiocParallel::bpparam()) else GSEA_WORKERS, RNG_kind = "Mersenne-Twister / Inversion / Rejection", 
    fgsea_version = as.character(utils::packageVersion("fgsea")), BiocParallel_version = as.character(utils::packageVersion("BiocParallel")), 
    rank_method = rnaseq_setting(RNA_DIR, "rank_method"), adjustment = "BH separately per comparison/ontology", 
    leading_edge = "Original DOSE hit-slice rule and gene order"), c(results = file.path(RNA_DIR, 
    "RNAseq_results.tsv.gz"), manifest = file.path(RNA_DIR, "Comparison_manifest.tsv"), RNA_settings = file.path(RNA_DIR, 
    "Analysis_settings.tsv"), annotation = go_path))
  # 7. Export enrichment, ranked input genes, diagnostics and software/seed
  # metadata. This is a separate ranked analysis from the DEG-list ORA in task 19.
  rnaseq_write(list(GO_GSEA_results = all_results, GO_GSEA_ranks = data.table::rbindlist(ranks), 
    GO_GSEA_diagnostics = data.table::rbindlist(diagnostics), Analysis_settings = settings), 
    OUTPUT_DIR, OVERWRITE, gzip = c("GO_GSEA_results", "GO_GSEA_ranks"))
  invisible(all_results)
}

main()

