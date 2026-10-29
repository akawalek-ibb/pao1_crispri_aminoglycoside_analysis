# Purpose: Prepare untreated MH fitness, external-set membership and gene counts
# for Supplementary Table S5. Also prepare the annotated Fig. 2G input.
# Prerequisites: 02_AMG_metacounter.R and 11_AMG_gene_vulnerability.R.
# Inputs: four *_median_melted_results.tsv[.gz] files, curated annotation,
# frozen published reference memberships and the MH vulnerability consensus.
# Rules: raw Stouffer p < 0.05 and medLFC <= -1 / >= 1; no new BH filter.
# Fitness union includes B1 T1/T2 and B2-B4 T1; vulnerability uses T1 only.
# Outputs: plain TSV files in out_12_AMG_MH_fitness/. No fitting is repeated.
# Run from AMG_cleaned_Git.Rproj:
# source("preprocessing_scripts/12_AMG_MH_fitness.R")
# USER SETTINGS
INPUT_DIR <- "edgeR_data"

ANNOTATION_FILE <- "annotation_input/PAO1_curated.tsv"

EXTERNAL_SET_FILE <- "annotation_input/AMG_published_fitness_membership.tsv"

VULNERABILITY_FILE <- "preprocessing_scripts/out_11_AMG_gene_vulnerability/MH_vulnerability.tsv"

OUTPUT_DIR <- "preprocessing_scripts/out_12_AMG_MH_fitness"

P_STOUFFER_CUT <- 0.05

LFC_CUT <- 1

OVERWRITE <- TRUE

COMPARISONS <- utils::read.table(text = "batch|condition|prefix|display_label\nb1|None_0_T1 - None_0_T0|B1_T1_vs_T0|B1 control T1 vs T0\nb1|None_0_T2 - None_0_T0|B1_T2_vs_T0|B1 control T2 vs T0\nb2|None_0_T1 - None_0_T0|B2_T1_vs_T0|B2 control T1 vs T0\nb3|None_0_T1 - None_0_T0|B3_T1_vs_T0|B3 control T1 vs T0\nb4|None_0_T1 - None_0_T0|B4_T1_vs_T0|B4 control T1 vs T0", 
  header = TRUE, sep = "|", stringsAsFactors = FALSE)

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
    stop("Invalid fitness thresholds.")
  # 1. Read the gene annotation and published essentiality memberships.
  # Published sets are reference labels for overlap comparisons; they do not
  # determine whether this study calls a negative-fitness effect.
  annotation <- amg_annotation(ANNOTATION_FILE)
  amg_required(annotation, "Nr_E", ANNOTATION_FILE)
  external <- amg_read(EXTERNAL_SET_FILE)
  external_columns <- c("Poulsen_core", "Zhang_2025", "Kaczmarczyk_2025")
  amg_required(external, c("PAO1_ID", external_columns), EXTERNAL_SET_FILE)
  amg_unique(external, "PAO1_ID", EXTERNAL_SET_FILE)
  for (column in external_columns) {
    value <- toupper(as.character(external[[column]]))
    if (anyNA(value) || any(!value %in% c("TRUE", "FALSE"))) 
      stop(EXTERNAL_SET_FILE, ": invalid membership in ", column)
    data.table::set(external, j = column, value = value == "TRUE")
  }
  files <- stats::setNames(vapply(unique(COMPARISONS$batch), function(batch) amg_gene_file(INPUT_DIR, 
    batch), character(1)), unique(COMPARISONS$batch))
  # 2. Read the five untreated comparisons (four T1 batches plus B1 T2).
  # Use the stored gene medians and Stouffer p-values; no guide-level model is
  # refitted. Finite results are tested, even when they do not meet hit cutoffs.
  raw <- data.table::rbindlist(lapply(names(files), function(batch) {
    x <- amg_gene_results(files[[batch]], COMPARISONS$condition[COMPARISONS$batch == batch])
    x[, `:=`(batch, batch)]
    x
  }), use.names = TRUE, fill = TRUE)
  raw[, `:=`(tested, is.finite(medLFC) & is.finite(p_stouffer))]
  # 3. Call negative/positive fitness only when both the raw-p-value and effect
  # criteria pass. The union flags below mean a call in at least one comparison;
  # they are not a requirement for agreement across all batches.
  raw[, `:=`(call, data.table::fcase(!tested, NA_character_, tested & p_stouffer < P_STOUFFER_CUT & 
    medLFC <= -LFC_CUT, "negative_fitness", tested & p_stouffer < P_STOUFFER_CUT & medLFC >= 
    LFC_CUT, "positive_fitness", default = "not_significant"))]
  ids <- unique(raw$PAO1_ID)
  if (length(setdiff(ids, annotation$PAO1_ID))) 
    stop("Gene results contain loci absent from curated annotation: ", paste(head(setdiff(ids, 
      annotation$PAO1_ID), 8), collapse = ", "))
  annotation_columns <- c("PAO1_ID", "PA14_ID", "gene_name", "product_name", "category", "design_call", 
    "TU")
  fitness <- annotation[match(ids, PAO1_ID), ..annotation_columns]
  calls <- raw[, .(CRISPRi_negative_fitness_union = any(call == "negative_fitness", na.rm = TRUE), 
    CRISPRi_positive_fitness_union = any(call == "positive_fitness", na.rm = TRUE), n_negative_comparisons = sum(call == 
      "negative_fitness", na.rm = TRUE), n_positive_comparisons = sum(call == "positive_fitness", 
      na.rm = TRUE)), by = PAO1_ID]
  fitness <- merge(fitness, calls, by = "PAO1_ID", all.x = TRUE, sort = FALSE)
  # 4. Add each comparison's effect, p-value and call to the gene-level table,
  # and count tested/negative/positive genes separately for each experiment.
  numbers <- vector("list", nrow(COMPARISONS))
  for (i in seq_len(nrow(COMPARISONS))) {
    spec <- COMPARISONS[i, , drop = FALSE]
    z <- raw[batch == spec$batch & condition == spec$condition]
    idx <- match(fitness$PAO1_ID, z$PAO1_ID)
    fitness[, `:=`((paste0(spec$prefix, c("_medLFC", "_p_stouffer", "_call"))), list(z$medLFC[idx], 
      z$p_stouffer[idx], z$call[idx]))]
    numbers[[i]] <- data.table::data.table(batch = toupper(spec$batch), condition = spec$condition, 
      display_label = spec$display_label, n_genes_tested = sum(z$tested), n_negative_fitness = sum(z$call == 
        "negative_fitness", na.rm = TRUE), n_positive_fitness = sum(z$call == "positive_fitness", 
        na.rm = TRUE), n_significant_total = sum(z$call %in% c("negative_fitness", "positive_fitness")))
  }
  data.table::setorder(fitness, -CRISPRi_negative_fitness_union, -CRISPRi_positive_fitness_union, 
    -n_negative_comparisons, -n_positive_comparisons, PAO1_ID)
  # 5. Build the essentiality membership table from the union of published sets
  # and this study's negative-fitness genes. Keep separate TRUE/FALSE flags for
  # each reference, allowing overlap figures without merging their definitions.
  membership_ids <- union(external$PAO1_ID[rowSums(as.data.frame(external[, ..external_columns])) > 
    0], calls$PAO1_ID[calls$CRISPRi_negative_fitness_union])
  if (length(setdiff(membership_ids, annotation$PAO1_ID))) 
    stop("Published memberships include unannotated loci; update the curated mapping explicitly.")
  membership <- annotation[match(membership_ids, PAO1_ID), .(PAO1_ID, PA14_ID, gene_name)]
  membership <- merge(membership, external[, c("PAO1_ID", external_columns), with = FALSE], by = "PAO1_ID", 
    all.x = TRUE, sort = FALSE)
  for (column in external_columns) data.table::set(membership, which(is.na(membership[[column]])), 
    column, FALSE)
  membership[, `:=`(This_study, PAO1_ID %in% calls$PAO1_ID[calls$CRISPRi_negative_fitness_union])]
  data.table::setcolorder(membership, c("PAO1_ID", "PA14_ID", "gene_name", "Poulsen_core", "Zhang_2025", 
    "This_study", "Kaczmarczyk_2025"))
  data.table::setorder(membership, PAO1_ID)
  # 6. Import the MH vulnerability consensus from task 11; do not recalculate
  # sensitivity or EK50 here. Verify that consensus severity matches the median
  # of the four T1 batch effects, so mixed input runs cannot silently pass.
  vulnerability <- amg_read(VULNERABILITY_FILE)
  vulnerability_columns <- c("PAO1_ID", "PA14_ID", "gene_name", "control_median_severity", "control_median_sensitivity", 
    "control_median_EK50", "n_control_severity_estimates", "n_control_sensitivity_estimates", 
    "n_control_EK50_estimates")
  amg_required(vulnerability, vulnerability_columns, VULNERABILITY_FILE)
  amg_unique(vulnerability, "PAO1_ID", VULNERABILITY_FILE)
  amg_numeric(vulnerability, vulnerability_columns[4:9], VULNERABILITY_FILE)
  if (!setequal(vulnerability$PAO1_ID, annotation$PAO1_ID)) 
    stop("The vulnerability consensus and curated annotation have different gene universes.")
  expected_consensus <- raw[condition == "None_0_T1 - None_0_T0", .(expected_severity = if (any(is.finite(medLFC))) 
    stats::median(medLFC[is.finite(medLFC)])
  else NA_real_, expected_n = sum(is.finite(medLFC))), by = PAO1_ID]
  idx <- match(vulnerability$PAO1_ID, expected_consensus$PAO1_ID)
  expected_severity <- expected_consensus$expected_severity[idx]
  expected_n <- expected_consensus$expected_n[idx]
  expected_n[is.na(expected_n)] <- 0L
  finite <- is.finite(expected_severity) & is.finite(vulnerability$control_median_severity)
  if (any(is.finite(expected_severity) != is.finite(vulnerability$control_median_severity)) || 
    any(abs(expected_severity[finite] - vulnerability$control_median_severity[finite]) > 1e-08) || 
    any(expected_n != vulnerability$n_control_severity_estimates)) 
    stop("The MH consensus and metacounter fitness results disagree. Re-run 11 on the current inputs.")
  vulnerability <- vulnerability[, ..vulnerability_columns]
  data.table::setorderv(vulnerability, c("control_median_EK50", "control_median_sensitivity", 
    "control_median_severity", "PAO1_ID"), na.last = TRUE)
  # 7. Add essentiality display categories and the number of annotated essential
  # genes in each TU to the Figure 2G input. These annotations describe the
  # consensus measures; they do not change their values or call genes essential.
  fig2g <- merge(vulnerability, annotation[, .(PAO1_ID, category, design_call, TU, Nr_E = as.numeric(Nr_E))], 
    by = "PAO1_ID", all.x = TRUE, sort = FALSE)
  core <- union(external$PAO1_ID[external$Poulsen_core], annotation$PAO1_ID[grepl("core", tolower(annotation$category)) & 
    !grepl("conditional", tolower(annotation$category))])
  conditional <- setdiff(annotation$PAO1_ID[grepl("conditional", tolower(annotation$category))], 
    core)
  fig2g[, `:=`(category_label, data.table::fcase(PAO1_ID %in% core, "Core essential", PAO1_ID %in% 
    conditional, "Conditional essential", design_call %in% c("main_essential", "other_essential") | 
    (!is.na(Nr_E) & Nr_E >= 1) | category %in% c("E", "E_AMR", "AMR"), "Essential others", 
    default = "Rest"))]
  essential_counts <- annotation[!is.na(TU) & nzchar(TU), .(essential_genes_in_operon = sum(design_call %in% 
    c("main_essential", "other_essential"))), by = TU]
  fig2g <- merge(fig2g, essential_counts, by = "TU", all.x = TRUE, sort = FALSE)
  fig2g[is.na(essential_genes_in_operon), `:=`(essential_genes_in_operon, 0L)]
  fig2g[, `:=`(severity = control_median_severity, sensitivity = control_median_sensitivity)]
  data.table::setorder(fig2g, PAO1_ID)
  summary <- data.table::data.table(metric = c("n_MH_fitness", "n_negative_union", "n_positive_union", 
    "n_membership", "n_MH_vulnerability", "n_Kaczmarczyk_reference"), value = c(nrow(fitness), 
    sum(fitness$CRISPRi_negative_fitness_union), sum(fitness$CRISPRi_positive_fitness_union), 
    nrow(membership), nrow(vulnerability), sum(external$Kaczmarczyk_2025)))
  settings <- amg_settings(c(P_STOUFFER_CUT = P_STOUFFER_CUT, LFC_CUT = LFC_CUT, union_comparisons = nrow(COMPARISONS), 
    vulnerability_consensus = "B1-B4 untreated T1 only; supplied by 11_AMG_gene_vulnerability"), 
    c(files, annotation = ANNOTATION_FILE, external = EXTERNAL_SET_FILE, vulnerability = VULNERABILITY_FILE))
  # 8. Export shared Table S5/figure inputs, summary counts and input provenance.
  amg_write(list(MH_fitness = fitness, Essentiality_membership = membership, Gene_numbers = data.table::rbindlist(numbers), 
    MH_vulnerability = vulnerability, Fig2G_input = fig2g, Summary = summary, Analysis_settings = settings), 
    OUTPUT_DIR, OVERWRITE)
  message("Prepared S5: ", sum(fitness$CRISPRi_negative_fitness_union), " negative-fitness genes; ", 
    nrow(membership), " membership records. Saved to ", OUTPUT_DIR)
  invisible(list(fitness = fitness, membership = membership, vulnerability = vulnerability))
}

main()

