# Purpose: Summarise MH/FBS gene fitness effects by KEGG pathway and run GSEA.
# Shared analysis for Table S7 and Fig. 2H. Their separate scripts import these
# results to assemble the workbook or draw the figure, without repeating GSEA.
# Prerequisite: 13_AMG_MH_FBS_gene_effects.R.
# Run from AMG_cleaned_Git.Rproj:
# source("preprocessing_scripts/14_AMG_MH_FBS_KEGG.R")
#
# Data flow:
# 1. Read all gene effects and the KEGG memberships in the annotation file.
# 2. Join genes to pathways; calculate pathway medians and representation counts.
# 3. Rank represented gene symbols by FBS-minus-MH effects and run GSEA.
# 4. Select Fig. 2H pathways: >=10 paired genes, then the 20 largest absolute
#    pathway effects. Selection does not use GSEA p values or NES magnitude.
# 5. Export full pathway results, selected figure rows and GSEA source data.
#
# Keep the two summaries distinct: MH/FBS pathway medians describe the overall
# fitness effects in each medium; GSEA tests whether a pathway's genes tend to
# occur toward either end of the ranked FBS-minus-MH gene-effect list.
# Positive NES indicates enrichment toward higher FBS-minus-MH ranks; negative
# NES indicates enrichment toward lower ranks. MH/FBS labels describe relative
# tendencies, rather than every gene necessarily changing in the same direction.
# These are observed screen effects; no scaling by culture generations is applied.
# No edgeR modelling, vulnerability-curve fitting or EK50 calculation is run.
#
# Outputs in OUTPUT_DIR (plain TSV):
# - KEGG_MH_vs_FBS_Fig2H: complete pathway table imported by Table S7.
# - Fig2H_input: selected pathway rows imported by the Fig. 2H plotting script.
# - GSEA_ranked_genes / GSEA_membership: the ranking and represented gene sets.
# - Summary / Analysis_settings: result counts, parameters and input provenance.
# Pathways outside GSEA size limits remain in the full table with NA GSEA.
# Only KEGG memberships are retained; no extra pathway exclusions are introduced.

# 1. USER SETTINGS -------------------------------------------------------------
# Full gene-level output from task 13, not only its significant gene subset.
GENE_EFFECTS_FILE <- "preprocessing_scripts/out_13_AMG_MH_FBS_gene_effects/MH_FBS_gene_effects.tsv"

# Annotation rows link PAO1 locus IDs to pathway names and their source database.
PATHWAY_FILE <- "annotation_input/PAO1_pathways_all_withCOG.tsv"

# Shared output folder used by the Table S7 and Fig. 2H scripts.
OUTPUT_DIR <- "preprocessing_scripts/out_14_AMG_MH_FBS_KEGG"

# The annotation file can contain several databases; use only KEGG here.
PATHWAY_DATABASE <- "KEGG"

# GSEA size limits count represented gene symbols in the ranked list,
# rather than every gene in the original pathway annotation.
GSEA_MIN_SIZE <- 10L

GSEA_MAX_SIZE <- 500L

# Permutation count and seed retain the existing GSEA calculation.
GSEA_NPERM <- 10000L

GSEA_SEED <- 1L

# Single-process execution; retain the seed, package version and inputs
# together when reproducing permutation-based results.
GSEA_NPROC <- 1L

# Figure eligibility requires this many genes measured in BOTH media.
FIGURE_MIN_PAIRED_GENES <- 10L

# Select the largest absolute pathway effects, not the smallest GSEA p values.
FIGURE_TOP_N <- 20L

# TRUE replaces existing TSVs; FALSE stops if any target already exists.
OVERWRITE <- TRUE

# 2. LOCAL FUNCTIONS -----------------------------------------------------------
# Clean HTML tags/entities and repeated spaces in the supplied pathway labels.
clean_pathway_name <- function(x) {
  x <- gsub("<[^>]+>", "", as.character(x))
  x <- gsub("&beta;", "beta", x, fixed = TRUE)
  x <- gsub("&gamma;", "gamma", x, fixed = TRUE)
  x <- gsub("&alpha;", "alpha", x, fixed = TRUE)
  x <- gsub("&amp;", "&", x, fixed = TRUE)
  trimws(gsub("[[:space:]]+", " ", x))
}

# Check packages without installing or updating them.
amg_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) 
    stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

# Stop early if an input lacks columns needed by the analysis.
amg_required <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) 
    stop(label, " lacks: ", paste(missing, collapse = ", "), call. = FALSE)
}

# Read plain/gzipped TSV while preserving its original column names.
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

# Require one unambiguous record for each key, such as a PAO1 locus ID.
amg_unique <- function(x, columns, label) {
  if (anyNA(x[, columns, with = FALSE]) || anyDuplicated(x[, columns, with = FALSE])) 
    stop(label, ": missing or duplicate key: ", paste(columns, collapse = ", "), call. = FALSE)
  invisible(x)
}

# Convert numeric input columns and reject malformed/non-finite nonmissing values.
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

# Write each named result as a TSV, respecting the existing overwrite setting.
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

# Record settings, input paths/checksums and the R version for traceability.
amg_settings <- function(settings, inputs) {
  data.table::data.table(setting = c(names(settings), paste0("input:", names(inputs)), paste0("md5:", 
    names(inputs)), "R_version"), value = c(as.character(settings), unname(inputs), unname(tools::md5sum(inputs)), 
    R.version.string))
}

# 3. ANALYSIS ------------------------------------------------------------------
# main() follows the data flow below. All input/output paths are set above.
main <- function() {
  # 3.1 Check parameters before reading the input files.
  amg_require(c("data.table", "fgsea"))
  integer_settings <- c(GSEA_MIN_SIZE, GSEA_MAX_SIZE, GSEA_NPERM, GSEA_NPROC, FIGURE_MIN_PAIRED_GENES, 
    FIGURE_TOP_N)
  if (any(!is.finite(integer_settings) | integer_settings < 1 | integer_settings%%1 != 0) || 
    GSEA_MAX_SIZE < GSEA_MIN_SIZE) 
    stop("Invalid pathway/GSEA settings.")

  # Read the complete gene-effect table.
  # MH_medLFC_vs_T0 and FBS_medLFC_vs_T0 describe fitness within each medium.
  # direct_FBS_vs_MH_medLFC is the separately analysed direct medium contrast.
  # Its significance flag comes from task 13 and is not recalculated here.
  genes <- amg_read(GENE_EFFECTS_FILE)
  amg_required(genes, c("PAO1_ID", "gene_name", "MH_medLFC_vs_T0", "FBS_medLFC_vs_T0", "direct_FBS_vs_MH_medLFC", 
    "direct_significant"), GENE_EFFECTS_FILE)
  amg_unique(genes, "PAO1_ID", GENE_EFFECTS_FILE)
  amg_numeric(genes, c("MH_medLFC_vs_T0", "FBS_medLFC_vs_T0", "direct_FBS_vs_MH_medLFC"), GENE_EFFECTS_FILE)
  if (!is.logical(genes$direct_significant) || anyNA(genes$direct_significant)) 
    stop("direct_significant must contain TRUE/FALSE calls from preprocessing 13.")

  # 3.2 Read and standardise pathway memberships.
  # Accept the existing underscore/space annotation-column naming conventions.
  pathways <- amg_read(PATHWAY_FILE)
  pick <- function(candidates) {
    found <- intersect(candidates, names(pathways))
    if (!length(found)) 
      stop(PATHWAY_FILE, " lacks ", paste(candidates, collapse = " or "))
    found[1L]
  }
  locus_col <- pick(c("Locus_Tag", "Locus Tag"))
  term_col <- pick(c("Pathway_Name", "Pathway Name"))
  db_col <- pick(c("Pathway_Database", "Pathway Database"))
  pathways <- data.table::data.table(PAO1_ID = trimws(as.character(pathways[[locus_col]])), pathway_term = clean_pathway_name(pathways[[term_col]]), 
    pathway_db = trimws(as.character(pathways[[db_col]])))
  # Retain KEGG rows with usable locus IDs and pathway names; remove duplicate
  # membership rows. One gene can legitimately belong to several pathways.
  pathways <- unique(pathways[pathway_db == PATHWAY_DATABASE & !is.na(PAO1_ID) & nzchar(PAO1_ID) & 
    !is.na(pathway_term) & nzchar(pathway_term)])
  if (!nrow(pathways)) 
    stop("No ", PATHWAY_DATABASE, " memberships found.")
  data.table::setorder(pathways, pathway_db, pathway_term, PAO1_ID)

  # 3.3 Attach experimental effects to each annotated pathway membership.
  # all.x = TRUE retains annotation members with no experimental gene record,
  # keeping the full pathway size distinguishable from experimental coverage.
  members <- merge(pathways, genes, by = "PAO1_ID", all.x = TRUE, sort = FALSE)

  # Ignore missing/non-finite values in a median; return NA if none exist.
  safe_median <- function(x) if (any(is.finite(x))) 
    stats::median(x[is.finite(x)])
  else NA_real_

  # Calculate one row per pathway:
  # - n_genes_in_pathway: all distinct annotated locus IDs.
  # - n_genes_with_MH_FBS_fitness: genes with finite effects in BOTH media.
  # - median_MH/FBS_medLFC: available gene-effect medians in each medium.
  # - median_direct_FBS_vs_MH_medLFC: median of the direct gene-contrast effects.
  # - n_genes_direct_significant: count of inherited significant gene calls.
  # The medium-specific medians use available values separately; the paired
  # count does not impose an additional paired-only filter on these medians.
  summary <- members[, .(n_genes_in_pathway = data.table::uniqueN(PAO1_ID), n_genes_with_MH_FBS_fitness = sum(is.finite(MH_medLFC_vs_T0) & 
    is.finite(FBS_medLFC_vs_T0)), median_MH_medLFC = safe_median(MH_medLFC_vs_T0), median_FBS_medLFC = safe_median(FBS_medLFC_vs_T0), 
    median_direct_FBS_vs_MH_medLFC = safe_median(direct_FBS_vs_MH_medLFC), n_genes_direct_significant = sum(direct_significant, 
      na.rm = TRUE)), by = .(pathway_db, pathway_term)]

  # Difference between the two pathway medians. This is a different summary
  # from the median of the independently analysed direct gene contrasts.
  summary[, `:=`(pathway_delta_FBS_minus_MH, median_FBS_medLFC - median_MH_medLFC)]

  # Figure-selection score: the larger absolute pathway median in either medium.
  # Pathways can score highly when their fitness effects are strong in BOTH.
  summary[, `:=`(absolute_fitness_score, pmax(abs(median_MH_medLFC), abs(median_FBS_medLFC)))]

  # 3.4 Prepare the GSEA gene ranking and represented gene sets.
  # GSEA uses gene symbols. If several loci share a symbol within a pathway,
  # collapse their effects to one median per medium before taking the difference.
  # Only named KEGG symbols with finite values in BOTH media enter the ranking.
  gsea_members <- members[!is.na(gene_name) & nzchar(gene_name), .(MH = safe_median(MH_medLFC_vs_T0), 
    FBS = safe_median(FBS_medLFC_vs_T0)), by = .(gene_name, pathway_db, pathway_term)]
  gsea_members <- gsea_members[is.finite(MH) & is.finite(FBS)]

  # Ranking score = FBS effect versus T0 minus MH effect versus T0.
  # Positive: more negative effects in MH. Negative: more negative in FBS.
  # Use all represented symbols, without filtering to significant gene contrasts.
  gsea_members[, `:=`(delta_gene, FBS - MH)]

  # Each symbol needs one score even if it belongs to several pathways.
  # Stop if a repeated symbol would acquire conflicting pathway-specific scores.
  ranked <- unique(gsea_members[, .(gene_name, delta_gene)])
  if (anyDuplicated(ranked$gene_name)) 
    stop("Some gene symbols have conflicting pathway-specific ranks. Resolve the symbol mapping explicitly before GSEA.")
  # Sort scores from largest to smallest; names give a reproducible ordering
  # among equal values without changing the numerical ranking scores.
  data.table::setorderv(ranked, c("delta_gene", "gene_name"), c(-1L, 1L))
  if (nrow(ranked) < GSEA_MIN_SIZE) 
    stop("Too few paired, named KEGG genes for GSEA.")

  # Build one unique-symbol set per pathway and a numeric vector named by the
  # same symbols. fgsea can then match pathway members to measured gene effects.
  data.table::setorder(gsea_members, pathway_db, pathway_term, gene_name)
  sets <- gsea_members[, .(genes = list(unique(gene_name))), by = .(pathway_db, pathway_term)]
  sets[, `:=`(set_id, paste(pathway_db, pathway_term, sep = "||"))]
  pathway_list <- stats::setNames(sets$genes, sets$set_id)
  ranked_vector <- stats::setNames(ranked$delta_gene, ranked$gene_name)

  # 3.5 Run permutation-based GSEA with the existing size limits and seed.
  # This tests concentration of pathway genes toward the high/low ends of the
  # ranked effects, rather than testing whether a pathway median equals zero.
  set.seed(GSEA_SEED)
  enrichment <- data.table::as.data.table(fgsea::fgseaSimple(pathways = pathway_list, stats = ranked_vector, 
    minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE, nperm = GSEA_NPERM, nproc = GSEA_NPROC))
  amg_required(enrichment, c("pathway", "pval", "padj", "ES", "NES", "leadingEdge"), "fgsea result")

  # Recover names from the temporary database||pathway IDs. Store leading-edge
  # gene symbols as semicolon-separated text for ordinary TSV use. p_gsea,
  # FDR_gsea, ES and NES_gsea retain the values returned by fgsea.
  enrichment[, `:=`(pathway_db = sub("\\|\\|.*$", "", pathway), pathway_term = sub("^.*?\\|\\|", 
    "", pathway), leadingEdge = vapply(leadingEdge, paste, collapse = "; ", FUN.VALUE = character(1)))]
  enrichment <- enrichment[, .(pathway_db, pathway_term, p_gsea = pval, FDR_gsea = padj, ES, 
    NES_gsea = NES, leadingEdge)]

  # Retain the complete pathway summary, including pathways outside GSEA size
  # limits. Their descriptive statistics remain, with NA enrichment results.
  summary <- merge(summary, enrichment, by = c("pathway_db", "pathway_term"), all.x = TRUE, sort = FALSE)

  # 3.6 Select Fig. 2H pathways independently of GSEA significance.
  # First rank by absolute_fitness_score (alphabetical name breaks score ties).
  # Then require >= FIGURE_MIN_PAIRED_GENES and take the first FIGURE_TOP_N rows.
  data.table::setorderv(summary, c("absolute_fitness_score", "pathway_term"), c(-1L, 1L), na.last = TRUE)
  selected <- head(summary[n_genes_with_MH_FBS_fitness >= FIGURE_MIN_PAIRED_GENES & is.finite(absolute_fitness_score)], 
    FIGURE_TOP_N)
  if (!nrow(selected)) 
    stop("No pathways meet the Fig. 2H representation requirement.")

  # Record figure order in both the selected subset and the full summary.
  selected[, `:=`(Fig2H_rank, seq_len(.N))]
  summary[, `:=`(Fig2H_rank, selected$Fig2H_rank[match(pathway_term, selected$pathway_term)])]

  # Display classes use NES magnitude/sign, without a GSEA p-value filter.
  # cut(-NES) aligns positive NES with MH labels and negative NES with FBS.
  # Original NES values are preserved; missing NES remains unclassified.
  selected[, `:=`(NES_class, as.character(cut(-NES_gsea, breaks = c(-Inf, -1.5, -0.75, 0.75, 
    1.5, Inf), labels = c("Strong MH", "Moderate MH", "Neutral", "Moderate FBS", "Strong FBS"))))]

  # 3.7 Arrange the COMPLETE pathway table for the existing Table S7 sheet.
  # Preserve its column names/order and sorting for the workbook assembler.
  # The separate selected object retains the additional Fig. 2H plotting fields.
  columns <- c("pathway_term", "n_genes_in_pathway", "n_genes_with_MH_FBS_fitness", "median_MH_medLFC", 
    "median_FBS_medLFC", "pathway_delta_FBS_minus_MH", "median_direct_FBS_vs_MH_medLFC", "n_genes_direct_significant", 
    "p_gsea", "FDR_gsea", "ES", "NES_gsea", "leadingEdge")
  table <- summary[, ..columns]
  data.table::setorderv(table, c("median_MH_medLFC", "median_FBS_medLFC", "pathway_term"), na.last = TRUE)

  # Store parameters, input checksums and fgsea version so results can be traced
  # to the supplied files and the software used for the permutation analysis.
  settings <- amg_settings(c(PATHWAY_DATABASE = PATHWAY_DATABASE, GSEA_MIN_SIZE = GSEA_MIN_SIZE, 
    GSEA_MAX_SIZE = GSEA_MAX_SIZE, GSEA_NPERM = GSEA_NPERM, GSEA_SEED = GSEA_SEED, GSEA_NPROC = GSEA_NPROC, 
    GSEA_ranking = "FBS-minus-MH medLFC; paired KEGG-annotated gene symbols", FIGURE_MIN_PAIRED_GENES = FIGURE_MIN_PAIRED_GENES, 
    FIGURE_TOP_N = FIGURE_TOP_N, fgsea_version = as.character(utils::packageVersion("fgsea"))), 
    c(gene_effects = GENE_EFFECTS_FILE, pathways = PATHWAY_FILE))

  # Report full pathway count, tested GSEA sets, selected pathways and the number
  # of unique gene symbols in the GSEA ranking as separate coverage measures.
  counts <- data.table::data.table(metric = c("n_pathways", "n_GSEA_pathways", "n_Fig2H_pathways", 
    "n_GSEA_ranked_genes"), value = c(nrow(summary), nrow(enrichment), nrow(selected), nrow(ranked)))

  # 3.8 Export shared results once. Fig. 2H reads Fig2H_input; Table S7 reads
  # KEGG_MH_vs_FBS_Fig2H. The other TSVs expose the exact ranked input and
  # represented memberships for inspection and reproducibility.
  amg_write(list(KEGG_MH_vs_FBS_Fig2H = table, Fig2H_input = selected, GSEA_ranked_genes = ranked, 
    GSEA_membership = unique(gsea_members[, .(gene_name, pathway_db, pathway_term)]), Summary = counts, 
    Analysis_settings = settings), OUTPUT_DIR, OVERWRITE)
  message("Prepared S7 KEGG: ", nrow(summary), " pathways; ", nrow(enrichment), " with GSEA; ", 
    nrow(selected), " for Fig. 2H. Saved to ", OUTPUT_DIR)
  invisible(summary)
}

# Execute the analysis when this complete file is sourced in the project.
main()

