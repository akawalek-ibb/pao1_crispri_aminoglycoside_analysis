# ============================================================
# Batch-specific edgeR metacounter — all ordered group contrasts
#
# Purpose: Fit each batch independently and export guide-level
# melted_results and perfect-guide median_melted_results.
#
# Inputs:
# - Barcode_counting/count_tables/b*_overall_stats.tsv
# - Barcode_counting/count_tables/b*_experimental_design.tsv
# - annotation_input/PAO1_key.tsv, PAO1_curated.tsv, PAO1_orfs.bed
# Prerequisite: Review QC and edit experimental designs yourself.
# Packages: edgeR, limma, data.table, dplyr, tidyr, poolr.
#
# Outputs: edgeR_data/b*_melted_results.tsv.gz and
# edgeR_data/b*_median_melted_results.tsv.gz.
# Gene tables include p_stouffer_BH for Supplementary Tables S5 and S8.
# Optional fits/diagnostics: preprocessing_scripts/out_02_AMG_metacounter/.
# Run: source("preprocessing_scripts/02_AMG_metacounter.R")
#
# Analytical rules:
# - Mean untreated-T0 count >=100; replicate counts modeled directly.
# - filterByExpr, TMM, estimateGLMRobustDisp, robust glmQLFit.
# - All ordered group contrasts; control-guide median LFC centering.
# - Gene medians/Stouffer use perfect guides with raw LFC !=0.
# - Stouffer combination excludes nonfinite or p<=0 guide p-values.
# - BH adjustment across all finite gene tests within each batch/contrast.
# - Raw p_stouffer values are preserved; no interaction calls or sample exclusions.
# Shared table exports are staged until the full run succeeds.
# Failed file promotion attempts restore previous tables.
# ============================================================
# USER SETTINGS
BATCHES <- c("b1", "b2", "b3", "b4")

COUNT_DIR <- file.path("Barcode_counting", "count_tables")

ANNOTATION_DIR <- "annotation_input"

RESULT_DIR <- "edgeR_data"

OUTPUT_DIR <- file.path("preprocessing_scripts", "out_02_AMG_metacounter")

KEY_FILENAME <- "PAO1_key.tsv"

CURATED_FILENAME <- "PAO1_curated.tsv"

BED_FILENAME <- "PAO1_orfs.bed"

# Minimum mean read count across untreated T0 replicates in a batch.
MIN_MEAN_T0 <- 100

# Optional tables help inspect intermediate steps; they do not change fits.
EXPORT_DIAGNOSTICS <- FALSE

# Print each comparison name while testing it.
VERBOSE_CONTRASTS <- FALSE

# Completed outputs are replaced only after all selected batches succeed.
OVERWRITE_RESULTS <- TRUE

RANDOM_SEED <- 123

# Local functions for this analysis.
run_analysis <- function() {
  # Output bookkeeping: write temporary tables while analysing batches.
  # commit_exports() installs them at their final paths only at the end.
  # A failed run removes its pending tables; it does not save fit checkpoints.
  export_state <- new.env(parent = emptyenv())
  export_state$files <- list()
  on.exit({
    pending <- unlist(export_state$files, use.names = FALSE)
    if (length(pending)) unlink(pending)
  }, add = TRUE)
  packages <- c("edgeR", "limma", "data.table", "dplyr", "tidyr", "poolr")
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) 
    stop("Install required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(poolr)
  })
  check_columns <- function(data, required, label) {
    missing <- setdiff(required, names(data))
    if (length(missing)) 
      stop(label, " lacks columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  register_export <- function(temporary, path) {
    if (!file.exists(temporary)) 
      stop("Staged export is missing: ", temporary)
    if (path %in% names(export_state$files)) 
      stop("Output registered twice: ", path)
    export_state$files[[path]] <- temporary
  }
  write_result <- function(data, path) {
    if (file.exists(path) && !OVERWRITE_RESULTS) 
      stop("Output exists: ", path)
    temporary <- tempfile(pattern = "pending_result_", tmpdir = dirname(path), fileext = if (endsWith(path, 
      ".gz")) 
      ".tsv.gz"
    else ".tsv")
    registered <- FALSE
    on.exit(if (!registered) unlink(temporary), add = TRUE)
    data.table::fwrite(data, temporary, sep = "\t", na = "NA")
    register_export(temporary, path)
    registered <- TRUE
  }
  commit_exports <- function() {
    files <- export_state$files
    if (!length(files)) 
      stop("No outputs were staged.")
    targets <- names(files)
    if (any(!file.exists(unlist(files, use.names = FALSE)))) 
      stop("A staged output is missing; completed results were not changed.")
    if (any(dir.exists(targets))) 
      stop("An output target is a directory; completed results were not changed.")
    if (!OVERWRITE_RESULTS && any(file.exists(targets))) 
      stop("Outputs already exist; completed results were not changed.")
    backups <- list()
    promoted <- character()
    completed <- FALSE
    on.exit({
      if (!completed) {
        if (length(promoted)) unlink(promoted)
        for (target in names(backups)) {
         if (file.exists(backups[[target]]) && !file.rename(backups[[target]], target)) {
          warning("Could not restore ", target, ". Previous file remains at: ", backups[[target]], 
           call. = FALSE)
         }
        }
      } else if (length(backups)) {
        unlink(unlist(backups, use.names = FALSE))
      }
    }, add = TRUE)
    for (target in targets[file.exists(targets)]) {
      backup <- tempfile(pattern = "previous_result_", tmpdir = dirname(target))
      if (!file.rename(target, backup)) 
        stop("Could not protect the previous result: ", target)
      backups[[target]] <- backup
    }
    for (target in targets) {
      if (!file.rename(files[[target]], target)) 
        stop("Could not promote staged output: ", target)
      promoted <- c(promoted, target)
    }
    completed <- TRUE
    export_state$files <- list()
  }
  normalize_dose <- function(x) {
    x_chr <- trimws(as.character(x))
    x_num <- suppressWarnings(as.numeric(x_chr))
    out <- x_chr
    numeric_idx <- is.finite(x_num) & !is.na(x_num)
    out[numeric_idx] <- vapply(x_num[numeric_idx], function(z) {
      if (z == 0) 
        return("0")
      format(z, scientific = FALSE, trim = TRUE, digits = 15)
    }, character(1))
    out
  }
  # Combine guide p-values into one gene-level p-value. Effect direction
  # comes from the median fold change, not from this p-value combination.
  # Missing, nonfinite and zero p-values are excluded here.
  stouffer_safe <- function(p) {
    p <- as.numeric(p)
    p <- p[is.finite(p) & !is.na(p) & p > 0 & p <= 1]
    if (length(p) == 0) 
      return(NA_real_)
    p <- pmax(p, .Machine$double.xmin)
    poolr::stouffer(p)$p
  }
  # Add multiple-testing-adjusted gene p-values without replacing raw ones.
  bh_adjust_complete <- function(p) {
    p <- as.numeric(p)
    if (any(is.finite(p) & (p < 0 | p > 1))) 
      stop("Invalid gene-level Stouffer p-value.")
    adjusted <- rep(NA_real_, length(p))
    eligible <- is.finite(p) & p >= 0 & p <= 1
    adjusted[eligible] <- stats::p.adjust(p[eligible], method = "BH")
    adjusted
  }
  # Read sample metadata and form groups from drug, dose and time point.
  # For example, untreated T1 replicates share the group None_0_T1.
  # Here "condition" means a sample ID; in output tables it means a contrast.
  standardize_design <- function(design_file, batch_id) {
    if (!file.exists(design_file)) {
      stop("Missing design file: ", design_file)
    }
    design <- data.table::fread(design_file, sep = "\t", header = TRUE, fill = TRUE, na.strings = c("", 
      "NA"), data.table = TRUE)
    required_base <- c("condition", "experiment", "drug", "dose", "timing", "rep")
    missing_required <- setdiff(required_base, colnames(design))
    if (length(missing_required) > 0) {
      stop("Design ", batch_id, " is missing required columns: ", paste(missing_required, 
        collapse = ", "))
    }
    design[, `:=`(condition, trimws(as.character(condition)))]
    design[, `:=`(experiment, trimws(as.character(experiment)))]
    design[, `:=`(drug, trimws(as.character(drug)))]
    design[, `:=`(timing, toupper(trimws(as.character(timing))))]
    design[, `:=`(rep, as.character(rep))]
    design[is.na(drug) | drug == "" | tolower(drug) == "none", `:=`(drug, "None")]
    design[, `:=`(dose, normalize_dose(dose))]
    design[drug == "None", `:=`(dose, "0")]
    design[, `:=`(group, paste(drug, dose, timing, sep = "_"))]
    if (anyDuplicated(design$condition)) {
      duplicated_conditions <- unique(design$condition[duplicated(design$condition)])
      stop("Duplicated sample names in ", batch_id, " design: ", paste(duplicated_conditions, 
        collapse = ", "))
    }
    if (any(is.na(design$condition) | design$condition == "")) {
      stop("Empty condition name detected in ", batch_id, " design")
    }
    design <- design[, .(condition, experiment, drug, dose, timing, rep, group)]
    design[, `:=`(batch, batch_id)]
    data.table::setcolorder(design, c("batch", "condition", "experiment", "drug", "dose", "timing", 
      "rep", "group"))
    design
  }
  # Check the working directory, settings and required files before fitting.
  if (!file.exists("AMG_cleaned_Git.Rproj")) 
    stop("Open AMG_cleaned_Git.Rproj and run from its working directory.", call. = FALSE)
  INPUT_DIR <- file.path(COUNT_DIR)
  ANNOT_DIR <- file.path(ANNOTATION_DIR)
  RESULT_DIR <- file.path(RESULT_DIR)
  OUTPUT_DIR <- file.path(OUTPUT_DIR)
  if (!length(BATCHES) || anyDuplicated(BATCHES) || any(!BATCHES %in% c("b1", "b2", "b3", "b4"))) 
    stop("Choose unique BATCHES from b1–b4.")
  if (!is.finite(MIN_MEAN_T0) || MIN_MEAN_T0 < 0) 
    stop("MIN_MEAN_T0 must be finite and nonnegative.")
  for (setting in c("EXPORT_DIAGNOSTICS", "VERBOSE_CONTRASTS", "OVERWRITE_RESULTS")) {
    value <- get(setting)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) 
      stop(setting, " must be TRUE or FALSE.")
  }
  BATCH_CONFIG <- setNames(lapply(BATCHES, function(batch) list(design_file = file.path(INPUT_DIR, 
    paste0(batch, "_experimental_design.tsv")), counts_file = file.path(INPUT_DIR, paste0(batch, 
    "_overall_stats.tsv")))), BATCHES)
  annotation_files <- file.path(ANNOT_DIR, c(KEY_FILENAME, CURATED_FILENAME, BED_FILENAME))
  required_files <- c(annotation_files, unlist(lapply(BATCH_CONFIG, unlist), use.names = FALSE))
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files)) 
    stop("Missing inputs:\n", paste(missing_files, collapse = "\n"), call. = FALSE)
  for (directory in c(RESULT_DIR, OUTPUT_DIR)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  if (!OVERWRITE_RESULTS) {
    planned_outputs <- file.path(RESULT_DIR, c(paste0(BATCHES, "_melted_results.tsv.gz"), paste0(BATCHES, 
      "_median_melted_results.tsv.gz")))
    if (any(file.exists(planned_outputs))) 
      stop("Result files already exist; choose another RESULT_DIR or allow overwriting.")
  }
  # Read annotations once for reuse across batches. These identify guides,
  # genes and transcriptional units; they are not model covariates.
  set.seed(RANDOM_SEED)
  PAO1 <- data.table::fread(file.path(ANNOT_DIR, CURATED_FILENAME), sep = "\t")
  check_columns(PAO1, c("PAO1_ID", "PA14_ID", "gene_name", "category", "TU"), CURATED_FILENAME)
  curated_names <- PAO1[, .(PAO1_ID, PA14_ID, gene_name, category, TU)]
  pao_bed <- data.table::fread(file.path(ANNOT_DIR, BED_FILENAME), col.names = c("chromosome", 
    "left", "right", "locus_tag", "random", "strand", "gene_name", "product", "alt_gene_name"))
  pao_key <- data.table::fread(file.path(ANNOT_DIR, KEY_FILENAME), sep = "\t")
  check_columns(pao_key, c("spacer", "type", "locus_tag", "y_pred", "original", "offset"), KEY_FILENAME)
  pao_key[, `:=`(spacer, toupper(gsub("\\s+", "", as.character(spacer))))]
  if (!any(pao_key$type == "control", na.rm = TRUE)) 
    stop("The full key needs type='control' entries for LFC centering.")
  if (!any(pao_key$type == "perfect", na.rm = TRUE)) 
    stop("The full key needs type='perfect' entries for gene aggregation.")
  pao_genome <- pao_bed[pao_key[, .(spacer, type, locus_tag, y_pred, original, offset)], on = .(locus_tag)]
  annotation_md5 <- tools::md5sum(annotation_files)

  # ============================================================
  # ANALYSE ONE BATCH: repeat these steps independently for B1–B4
  # ============================================================
  process_batch <- function(batch_id, config) {
    message("Processing ", batch_id)
    # Record the batch input checksums for the run summary; fitting is unchanged.
    current_md5 <- c(
      design = unname(tools::md5sum(config$design_file)),
      counts = unname(tools::md5sum(config$counts_file))
    )

    # 1. Identify samples and their biological groups.
    # Sample names link the count table to treatment/time-point metadata.
    pao_design <- standardize_design(config$design_file, batch_id)
    if (!nrow(pao_design) || anyNA(pao_design$group)) 
      stop(batch_id, ": empty or incomplete design.")
    if (any(is.na(pao_design$timing) | !nzchar(pao_design$timing)) || any(is.na(pao_design$dose) | 
      pao_design$dose %in% c("", "NA", "NaN"))) {
      stop(batch_id, ": timing or dose metadata are missing. Review the experimental design.")
    }
    pao_design_model <- data.table::copy(pao_design)

    # 2. Read raw guide counts: spacer sequence, read count and sample ID.
    # Sum repeated records for the same spacer/sample. Samples missing from
    # the design are ignored; missing required samples stop the analysis.
    pao_counts <- data.table::fread(config$counts_file, sep = "\t", header = FALSE, fill = TRUE, 
      data.table = TRUE)
    if (ncol(pao_counts) < 3) {
      stop("Count file for ", batch_id, " must have at least three columns: spacer count condition")
    }
    pao_counts <- pao_counts[, 1:3]
    data.table::setnames(pao_counts, c("spacer", "count", "condition"))
    pao_counts <- pao_counts[!(spacer %in% c("barcode", "spacer") & count %in% c("count", "counts") & 
      condition %in% c("condition", "sample"))]
    pao_counts[, `:=`(spacer, as.character(spacer))]
    pao_counts[, `:=`(condition, as.character(condition))]
    pao_counts[, `:=`(count, as.numeric(count))]
    if (anyNA(pao_counts$count)) {
      stop("Non-numeric barcode counts detected in ", config$counts_file)
    }
    pao_counts <- pao_counts[, .(count = sum(count, na.rm = TRUE)), by = .(spacer, condition)]
    missing_in_counts <- setdiff(pao_design_model$condition, unique(pao_counts$condition))
    extra_in_counts <- setdiff(unique(pao_counts$condition), pao_design_model$condition)
    if (length(missing_in_counts) > 0) {
      stop("These design samples are missing from ", batch_id, " count table:\n", paste(missing_in_counts, 
        collapse = "\n"))
    }
    if (length(extra_in_counts) > 0) {
      warning("These count-table samples are not in the ", batch_id, " design and will be ignored:\n", 
        paste(extra_in_counts, collapse = "\n"))
    }
    pao_counts <- pao_counts[condition %in% pao_design_model$condition]
    if (any(!is.finite(pao_counts$count)) || any(pao_counts$count < 0)) 
      stop(batch_id, ": counts must be finite and nonnegative.")

    # 3. Arrange counts as guides (rows) by biological samples (columns).
    # A missing guide/sample record becomes zero. Column order follows the
    # design so that each count column receives the correct group label.
    pao_grid <- data.table::dcast(pao_counts, spacer ~ factor(condition, levels = pao_design_model$condition), 
      value.var = "count", fill = 0)
    colnames(pao_grid) <- gsub("^factor\\(condition, levels = pao_design_model\\$condition\\)", 
      "", colnames(pao_grid), fixed = FALSE)
    sample_cols <- pao_design_model$condition
    missing_cols <- setdiff(sample_cols, colnames(pao_grid))
    if (length(missing_cols) > 0) {
      stop("Missing sample columns after dcast for ", batch_id, ":\n", paste(missing_cols, 
        collapse = "\n"))
    }
    pao_grid <- pao_grid[, c("spacer", sample_cols), with = FALSE]

    # 4. Require adequate representation before the competition starts.
    # Keep guides with mean untreated-T0 count >= MIN_MEAN_T0 (default 100).
    # This filters guides, not samples, and uses raw counts. T0 counts are
    # not used as a divisor: group differences are estimated by the model.
    t0_samples <- pao_design_model[drug == "None" & timing == "T0", condition]
    if (length(t0_samples) == 0) {
      stop("No no-drug T0 samples were identified for ", batch_id)
    }
    keep_t0 <- rowMeans(as.matrix(pao_grid[, ..t0_samples]), na.rm = TRUE) >= MIN_MEAN_T0
    pao_grid <- pao_grid[keep_t0]
    if (!nrow(pao_grid)) 
      stop(batch_id, ": no guides survive the untreated-T0 filter.")
    if (any(colSums(as.matrix(pao_grid[, ..sample_cols])) <= 0)) 
      stop(batch_id, ": a sample has zero counts after the T0 guide filter.")
    if (EXPORT_DIAGNOSTICS) {
      write_result(pao_design, file.path(OUTPUT_DIR, paste0(batch_id, "_experimental_design_standardized.tsv")))
      write_result(pao_grid, file.path(OUTPUT_DIR, paste0(batch_id, "_count_matrix_T0filtered.tsv.gz")))
    }

    # 5. Tell edgeR which samples belong to each group.
    # Replicates remain separate observations. The design matrix has one
    # column per group; replicate numbers are metadata, not pairing terms.
    pao_grid_matrix <- data.matrix(pao_grid[, !names(pao_grid) %in% "spacer", with = FALSE])
    rownames(pao_grid_matrix) <- pao_grid$spacer
    pao_group <- factor(pao_design_model$group, levels = unique(pao_design_model$group))
    pao_permut <- model.matrix(~0 + pao_group)
    colnames(pao_permut) <- levels(pao_group)
    if (ncol(pao_permut) < 2L) 
      stop(batch_id, ": at least two groups are needed for contrasts.")
    if (nrow(pao_permut) <= qr(pao_permut)$rank) 
      stop(batch_id, ": no residual degrees of freedom; review the replicate design.")
    fit_source <- "fitted_de_novo"
    pao_y0 <- edgeR::DGEList(counts = pao_grid_matrix, group = pao_group, genes = rownames(pao_grid_matrix))
    # Apply edgeR's additional abundance filter using the group design.
    # A guide passing the T0 cutoff can still fail this filter.
    pao_keep <- edgeR::filterByExpr(pao_y0, design = pao_permut, group = pao_group)
    pao_y <- pao_y0[pao_keep, , keep.lib.sizes = FALSE]
    message(batch_id, ": guides retained after filterByExpr: ", nrow(pao_y$counts))

    # 6. Normalize libraries and fit the guide-abundance model.
    # TMM normalization accounts for library size/composition differences.
    # Dispersion estimates describe count variability between replicates;
    # robust estimation limits the influence of outliers. The quasi-
    # likelihood model then estimates group effects for each retained guide.
    pao_y <- pao_y %>% edgeR::calcNormFactors() %>% edgeR::estimateGLMRobustDisp(pao_permut)
    pao_fit <- edgeR::glmQLFit(pao_y, pao_permut, robust = TRUE)
    if (!nrow(pao_y$counts)) 
      stop(batch_id, ": no guides remain in the fitted object.")
    if (!nrow(pao_y$counts)) 
      stop(batch_id, ": fitted object is empty.")

    # 7. Test every ordered pair of groups within this batch.
    # A - B and B - A are both exported. For example,
    # None_0_T1 - None_0_T0 measures the change from untreated T0 to T1.
    # Negative log2FC means lower relative abundance in the first group.
    contrast_levels <- data.table::CJ(level2 = colnames(pao_permut), level1 = colnames(pao_permut))[level2 != 
      level1, paste(level2, level1, sep = " - ")]
    pao_contrast <- limma::makeContrasts(contrasts = contrast_levels, levels = pao_permut)
    message(batch_id, ": ", nrow(pao_y$counts), " fitted guides; ", length(contrast_levels), 
      " ordered contrasts.")
    results_pval <- pao_key[, .(genes = unique(spacer))]
    results_LFC <- pao_key[, .(genes = unique(spacer))]
    for (ii in seq_len(ncol(pao_contrast))) {
      fit_i <- edgeR::glmQLFTest(pao_fit, contrast = pao_contrast[, ii])
      tt <- edgeR::topTags(fit_i, n = Inf)
      res <- data.table::data.table(tt$table)
      res[, `:=`(pval, PValue)]
      if (VERBOSE_CONTRASTS) 
        message(batch_id, ": ", contrast_levels[ii])
      # Store each guide's raw PValue and log2FC, not the topTags FDR.
      # Guides absent from the fitted model remain missing in this key-based
      # scaffold and are removed later if they have no valid result.
      results_pval <- res[, .(genes, pval)][results_pval, on = .(genes)]
      data.table::setnames(results_pval, "pval", contrast_levels[ii])
      results_LFC <- res[, .(genes, logFC)][results_LFC, on = .(genes)]
      data.table::setnames(results_LFC, "logFC", contrast_levels[ii])
    }

    # 8. Attach gene names, guide classes and genomic annotations to results.
    # These joins label the fitted guide effects; they do not refit models.
    # Keep separate wide tables for p-values and fold changes at this stage.
    results_pval <- pao_genome[results_pval, on = .(spacer == genes)]
    results_pval <- dplyr::right_join(as.data.frame(curated_names), as.data.frame(results_pval), 
      by = c(PAO1_ID = "locus_tag"))
    results_pval <- data.table::as.data.table(results_pval)
    if ("gene_name.x" %in% colnames(results_pval)) {
      data.table::setnames(results_pval, "gene_name.x", "gene_name")
    }
    if ("gene_name.y" %in% colnames(results_pval)) {
      results_pval[, `:=`(gene_name.y, NULL)]
    }
    core_cols <- c("PAO1_ID", "PA14_ID", "gene_name", "spacer", "type", "y_pred", "original", 
      "offset", "category", "TU")
    missing_core_pval <- setdiff(core_cols, colnames(results_pval))
    if (length(missing_core_pval) > 0) {
      stop(batch_id, ": missing expected p-value annotation columns: ", paste(missing_core_pval, 
        collapse = ", "))
    }
    results_pval <- results_pval[, c(core_cols, contrast_levels), with = FALSE]
    results_LFC <- pao_genome[results_LFC, on = .(spacer == genes)]
    results_LFC <- dplyr::right_join(as.data.frame(curated_names), as.data.frame(results_LFC), 
      by = c(PAO1_ID = "locus_tag"))
    results_LFC <- data.table::as.data.table(results_LFC)
    if ("gene_name.x" %in% colnames(results_LFC)) {
      data.table::setnames(results_LFC, "gene_name.x", "gene_name")
    }
    if ("gene_name.y" %in% colnames(results_LFC)) {
      results_LFC[, `:=`(gene_name.y, NULL)]
    }
    missing_core_lfc <- setdiff(core_cols, colnames(results_LFC))
    if (length(missing_core_lfc) > 0) {
      stop(batch_id, ": missing expected LFC annotation columns: ", paste(missing_core_lfc, 
        collapse = ", "))
    }
    results_LFC <- results_LFC[, c(core_cols, contrast_levels), with = FALSE]
    if (EXPORT_DIAGNOSTICS) {
      write_result(results_pval, file.path(OUTPUT_DIR, paste0(batch_id, "_results_pval_wide.tsv.gz")))
      write_result(results_LFC, file.path(OUTPUT_DIR, paste0(batch_id, "_results_LFC_wide.tsv.gz")))
    }

    # 9. Convert the wide tables to one row per guide assignment/comparison.
    # Join the matching p-value and fold change. "condition" now contains
    # a comparison such as None_0_T1 - None_0_T0, rather than a sample ID.
    # Discard rows without both a p-value and a fold-change estimate.
    melted_results_pval <- results_pval %>% data.table::melt(id.vars = c("PAO1_ID", "PA14_ID", 
      "gene_name", "category", "spacer", "type", "y_pred", "original", "offset", "TU"), variable.name = "condition", 
      value.name = "pval", measure.vars = contrast_levels)
    melted_results_LFC <- results_LFC %>% data.table::melt(id.vars = c("PAO1_ID", "PA14_ID", 
      "gene_name", "category", "spacer", "type", "y_pred", "original", "offset", "TU"), variable.name = "condition", 
      value.name = "LFC", measure.vars = contrast_levels)
    melted_results <- dplyr::left_join(melted_results_LFC, melted_results_pval, by = c("PAO1_ID", 
      "PA14_ID", "gene_name", "category", "spacer", "type", "y_pred", "original", "offset", 
      "condition", "TU"))
    melted_results <- melted_results[!is.na(melted_results$pval) & !is.na(melted_results$LFC), 
      ]
    melted_results <- data.table::data.table(melted_results)
    data.table::setkey(melted_results, condition)

    # 10. Centre guide effects on non-targeting controls, then summarize genes.
    # Within each comparison, subtract the median control-guide LFC from
    # every guide LFC. This places the median non-targeting effect at zero.
    # It changes LFC.adj, but does not recalculate model p-values.
    melted_results[, `:=`(LFC.adj, LFC - median(LFC[type == "control"], na.rm = TRUE)), by = condition]
    # Gene summaries use full-match (perfect) guides only and omit guides
    # whose raw LFC is exactly zero. medLFC is the median adjusted effect;
    # p_stouffer combines the eligible raw guide p-values. Mismatch and
    # control guides remain available in the individual-guide output.
    median_melted_results <- melted_results %>% dplyr::filter(LFC != 0) %>% dplyr::filter(type == 
      "perfect") %>% data.table::as.data.table()
    median_melted_results <- median_melted_results[, .(medLFC = median(LFC.adj, na.rm = TRUE), 
      p_stouffer = stouffer_safe(pval)), by = .(PAO1_ID, PA14_ID, gene_name, type, condition)]

    # 11. Prepare identifiers and gene-level multiple-testing adjustment.
    # Split each comparison into shift (first group) and base (reference).
    # Full-match y_pred is nominally 1; it is not a measured knockdown level.
    melted_results[, `:=`(PAO1_ID, ifelse(PAO1_ID == "", NA, PAO1_ID))]
    melted_results[, `:=`(c("shift", "base"), data.table::tstrsplit(condition, " - "))]
    melted_results[type == "perfect", `:=`(y_pred, 1)]
    data.table::setorder(melted_results, pval)
    median_melted_results[, `:=`(PAO1_ID, ifelse(PAO1_ID == "", NA, PAO1_ID))]
    median_melted_results[, `:=`(c("shift", "base"), data.table::tstrsplit(condition, " - "))]
    named_gene <- !is.na(median_melted_results$PAO1_ID) & nzchar(trimws(median_melted_results$PAO1_ID))
    if (anyDuplicated(median_melted_results[named_gene, .(condition, PAO1_ID)])) {
      stop(batch_id, ": duplicate locus/contrast rows prevent unambiguous gene-level BH adjustment.")
    }
    # Adjust across finite gene p-values separately within each comparison.
    # Keep both raw p_stouffer and adjusted p_stouffer_BH. This script does
    # not select significant hits; downstream analyses apply their cutoffs.
    median_melted_results[, `:=`(p_stouffer_BH, NA_real_)]
    median_melted_results[named_gene, `:=`(p_stouffer_BH, bh_adjust_complete(p_stouffer)), 
      by = condition]
    data.table::setorder(median_melted_results, p_stouffer)
    melted_results <- melted_results[, .(PAO1_ID, PA14_ID, gene_name, category, spacer, type, 
      y_pred, original, offset, TU, condition, LFC, pval, LFC.adj, shift, base)]
    median_melted_results <- median_melted_results[, .(PAO1_ID, PA14_ID, gene_name, type, condition, 
      medLFC, p_stouffer, shift, base, p_stouffer_BH)]
    if (!nrow(melted_results) || !nrow(median_melted_results)) 
      stop(batch_id, ": guide/gene results are empty.")
    if (any(!is.finite(melted_results$LFC.adj))) 
      stop(batch_id, ": control-guide LFC centering failed; inspect control retention.")

    # 12. Stage the two result tables and return this batch's run summary.
    # Guide-level output retains individual-guide effects; gene-level output
    # contains full-match medians and combined p-values. "Prepared" means
    # staged temporarily, not yet installed at the final output filenames.
    guide_file <- file.path(RESULT_DIR, paste0(batch_id, "_melted_results.tsv.gz"))
    gene_file <- file.path(RESULT_DIR, paste0(batch_id, "_median_melted_results.tsv.gz"))
    write_result(melted_results, guide_file)
    write_result(median_melted_results, gene_file)
    message("Prepared ", basename(guide_file), " and ", basename(gene_file))
    # Record sample/guide counts, input checksums and software versions.
    # This is provenance information, not another biological calculation.
    data.frame(batch = batch_id, fit_source = fit_source, n_samples = nrow(pao_design), n_groups = ncol(pao_permut), 
      n_guides_T0_filtered = nrow(pao_grid), n_guides_fitted = nrow(pao_y$counts), n_ordered_contrasts = length(contrast_levels), 
      n_guide_result_rows = nrow(melted_results), n_gene_result_rows = nrow(median_melted_results), 
      design_md5 = current_md5[["design"]], counts_md5 = current_md5[["counts"]], key_md5 = unname(annotation_md5[1]), 
      curated_md5 = unname(annotation_md5[2]), bed_md5 = unname(annotation_md5[3]), min_mean_T0 = MIN_MEAN_T0, 
      R_version = R.version.string, edgeR_version = as.character(utils::packageVersion("edgeR")), 
      guide_results = file.path(RESULT_DIR, basename(guide_file)), gene_results = file.path(RESULT_DIR, 
        basename(gene_file)), stringsAsFactors = FALSE)
  }

  # Run the complete analysis above once per selected batch. Each call fits
  # its own model; gc() releases unused memory between calls. No batch
  # checkpoint is saved, so an interrupted run must start again.
  run_summary <- lapply(BATCHES, function(batch) {
    result <- process_batch(batch, BATCH_CONFIG[[batch]])
    gc(verbose = FALSE)
    result
  })
  write_result(data.table::rbindlist(run_summary), file.path(OUTPUT_DIR, "metacounter_run_summary.tsv"))
  # All batches and the summary are ready. Promote the staged files now;
  # existing completed results are protected during replacement.
  commit_exports()
  message("Metacounter completed. Shared guide/gene results: ", RESULT_DIR)
}

run_analysis()

