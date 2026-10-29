# ============================================================
# Sample QC and exploratory PCA — CRISPRi batches B1–B4
#
# Purpose: Review depth, guide representation, dominance and replicate
# similarity BEFORE differential fitness analysis. No sample is removed
# because it receives a QC flag. Only samples in the current design enter
# PCA; missing or zero-reference-count samples cannot enter normalization.
#
# Inputs:
# - Barcode_counting/count_tables/b*_overall_stats.tsv
# - Barcode_counting/count_tables/b*_experimental_design.tsv
# - Barcode_counting/Lib1_barcodes.fasta (the supplied counting reference).
# QC denominators follow that FASTA universe; they need not equal the canonical
# PAO1_key.tsv universe used in downstream key-based coverage tables.
# - Optional cutadapt reports from out_count_Lib1_barcodes/b*/
#
# Prerequisites: R packages data.table, edgeR, ggplot2, ggrepel, openxlsx.
# Run from RStudio: source("preprocessing_scripts/01_AMG_sample_QC_PCA.R")
# Outputs: preprocessing_scripts/out_01_AMG_sample_QC_PCA/
#   sample_QC_report.xlsx, QC/PCA PNG+PDF, optional TSVs.
#
# Interpretation: assigned counts are NOT raw sequencing-read counts.
# Guide representation after induction/treatment can change biologically.
# All numerical review cutoffs below are configurable heuristics, not
# validated pass/fail standards. This script does not modify the designs.
# ============================================================
# USER SETTINGS
BATCHES <- c("b1", "b2", "b3", "b4")

COUNT_ID_MODE <- "auto"

COUNT_DIR_RELATIVE <- file.path("Barcode_counting", "count_tables")

REFERENCE_FILE_RELATIVE <- file.path("Barcode_counting", "Lib1_barcodes.fasta")

CUTADAPT_REPORT_DIR_RELATIVE <- file.path("Barcode_counting", "out_count_Lib1_barcodes")

OUTPUT_DIR_RELATIVE <- file.path("preprocessing_scripts", "out_01_AMG_sample_QC_PCA")

# Count-based review thresholds (in the units supplied by overall_stats).
MIN_ASSIGNED_COUNTS <- 1e+06

MIN_T0_GUIDE_COVERAGE_PERCENT <- 80

MAX_T0_TOP_GUIDE_PERCENT <- 1

MAX_UNKNOWN_COUNT_PERCENT <- 1

MIN_WITHIN_GROUP_COUNT_RATIO <- 0.5

MIN_WITHIN_GROUP_COVERAGE_RATIO <- 0.8

MIN_REPLICATE_SPEARMAN <- 0.8

# Descriptive complexity metrics.
REPRESENTATION_COUNT_LEVELS <- c(1, 10, 100)

RAREFACTION_DEPTH <- 1e+06

# Rarefaction is NA for noninteger counts or depth exceeding sample counts.
# PCA settings. Guide filtering here affects diagnostics only, not edgeR fits.
PCA_MIN_COUNT <- 10

PCA_MIN_SAMPLES <- 2

PCA_TOP_VARIABLE_GUIDES <- 10000

PCA_LOGCPM_PRIOR_COUNT <- 1

# Each analysis is normalized independently with TMM. No batch correction.
# Features are centered; no unit-variance scaling or sample outlier removal.
# Plot/export settings.
LABEL_ALL_SAMPLES <- TRUE

FONT_FAMILY <- "sans"

PLOT_WIDTH <- 12

PLOT_HEIGHT <- 8

DPI <- 300

EXPORT_TSV <- FALSE

SHOW_PLOTS <- TRUE

PLOT_SEED <- 123

# LIBRARIES
REQUIRED_PACKAGES <- c("data.table", "edgeR", "ggplot2", "ggrepel", "openxlsx")

missing_packages <- REQUIRED_PACKAGES[!vapply(REQUIRED_PACKAGES, requireNamespace, logical(1), 
  quietly = TRUE)]

if (length(missing_packages)) {
  stop("Install required packages before running: ", paste(missing_packages, collapse = ", "), 
    call. = FALSE)
}

# HELPERS
normalize_dose <- function(x) {
  text <- trimws(as.character(x))
  numeric_value <- suppressWarnings(as.numeric(text))
  ok <- is.finite(numeric_value)
  text[ok] <- vapply(numeric_value[ok], function(z) format(z, scientific = FALSE, trim = TRUE, 
    digits = 15), character(1))
  text
}

read_reference <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  headers <- which(startsWith(lines, ">"))
  if (!length(headers) || headers[1] != 1L) 
    stop("Invalid barcode FASTA: ", path)
  ends <- c(headers[-1] - 1L, length(lines))
  sequences <- vapply(seq_along(headers), function(i) {
    if (ends[i] <= headers[i]) 
      return("")
    toupper(paste0(lines[seq.int(headers[i] + 1L, ends[i])], collapse = ""))
  }, character(1))
  ids <- sub("[[:space:]].*$", "", substring(lines[headers], 2))
  if (any(!nzchar(ids)) || anyDuplicated(ids) || any(!nzchar(sequences))) {
    stop("The barcode FASTA needs unique nonempty IDs and nonempty sequences.")
  }
  data.frame(fasta_id = ids, sequence = sequences, stringsAsFactors = FALSE)
}

read_design <- function(path, batch) {
  design <- as.data.frame(data.table::fread(path, sep = "\t", na.strings = c("", "NA")))
  required <- c("condition", "experiment", "drug", "dose", "timing", "rep")
  absent <- setdiff(required, names(design))
  if (length(absent)) 
    stop(batch, " design lacks columns: ", paste(absent, collapse = ", "))
  for (column in required) design[[column]] <- trimws(as.character(design[[column]]))
  if (anyNA(design$condition) || any(!nzchar(design$condition)) || anyDuplicated(design$condition)) {
    stop(batch, " design must have unique, nonempty condition values.")
  }
  design$drug[is.na(design$drug) | tolower(design$drug) %in% c("none", "na")] <- "None"
  design$drug[design$drug != "None"] <- toupper(design$drug[design$drug != "None"])
  design$dose <- normalize_dose(design$dose)
  design$dose[design$drug == "None"] <- "0"
  design$timing <- toupper(design$timing)
  if (anyNA(design$timing) || anyNA(design$dose) || any(!nzchar(design$timing)) || any(!nzchar(design$dose))) {
    stop(batch, " design needs nonempty timing and dose values.")
  }
  design$batch <- batch
  design$group <- paste(design$drug, design$dose, design$timing, sep = "_")
  # The design already contains the unique batch-prefixed sample ID.
  design$sample_key <- design$condition
  design$treatment <- ifelse(design$drug == "None", "None", paste(design$drug, design$dose))
  design
}

read_counts <- function(path) {
  counts <- data.table::fread(path, sep = "\t", header = FALSE, colClasses = "character")
  if (ncol(counts) != 3L) 
    stop("Expected exactly three columns: barcode, count, condition in ", path)
  data.table::setnames(counts, c("barcode", "count", "condition"))
  for (column in c("barcode", "count", "condition")) data.table::set(counts, j = column, value = trimws(counts[[column]]))
  header <- counts$barcode %in% c("barcode", "spacer") & counts$count %in% c("count", "counts") & 
    counts$condition %in% c("condition", "sample")
  counts <- counts[!header]
  counts$count <- suppressWarnings(as.numeric(counts$count))
  if (!nrow(counts) || anyNA(counts) || any(!is.finite(counts$count)) || any(counts$count < 0) || 
    any(!nzchar(counts$barcode)) || any(!nzchar(counts$condition))) {
    stop("Count file contains empty identifiers or invalid/nonfinite/negative counts: ", path)
  }
  counts
}

expected_rarefied_guides <- function(counts, depth) {
  if (any(abs(counts - round(counts)) > 1e-08)) 
    return(NA_real_)
  total <- sum(counts)
  if (total < depth || total == 0) 
    return(NA_real_)
  positive <- counts[counts > 0]
  certain <- total - positive < depth
  probabilities <- rep(1, length(positive))
  probabilities[!certain] <- -expm1(lchoose(total - positive[!certain], depth) - lchoose(total, 
    depth))
  sum(pmin(1, pmax(0, probabilities)))
}

sample_complexity <- function(counts) {
  total <- sum(counts)
  positive <- counts[counts > 0]
  p <- if (total > 0) 
    positive/total
  else numeric()
  entropy <- if (length(p)) 
    -sum(p * log(p))
  else NA_real_
  sorted <- sort(counts)
  gini <- if (total > 0) 
    sum((2 * seq_along(sorted) - length(sorted) - 1) * sorted)/(length(sorted) * total)
  else NA_real_
  result <- list(reference_assigned_counts = total, n_detected_guides = length(positive), detected_guides_percent = 100 * 
    length(positive)/length(counts), shannon_effective_guides = if (length(p)) exp(entropy) else NA_real_, 
    inverse_simpson_guides = if (length(p)) 1/sum(p^2) else NA_real_, gini = gini, top_guide_percent = if (total > 
      0) 100 * max(counts)/total else NA_real_, top_guide_id = if (total > 0) names(counts)[which.max(counts)] else NA_character_, 
    top10_guides_percent = if (total > 0) 100 * sum(head(sort(counts, decreasing = TRUE), 10))/total else NA_real_, 
    noninteger_reference_counts = sum(abs(counts - round(counts)) > 1e-08), rarefied_detected_guides = expected_rarefied_guides(counts, 
      RAREFACTION_DEPTH))
  for (level in REPRESENTATION_COUNT_LEVELS) result[[paste0("guides_ge_", level, "_percent")]] <- 100 * 
    mean(counts >= level)
  result
}

read_cutadapt_pairs <- function(batch, sample) {
  path <- file.path(CUTADAPT_REPORT_DIR_RELATIVE, batch, paste0(sample, "_cutadapt.txt"))
  result <- list(raw_read_pairs = NA_real_, retained_read_pairs = NA_real_)
  if (!file.exists(path)) 
    return(result)
  lines <- readLines(path, warn = FALSE)
  extract <- function(pattern) {
    matches <- grep(pattern, lines, value = TRUE)
    if (length(matches) != 1L) 
      return(NA_real_)
    token <- sub("^[^:]+:[[:space:]]*", "", matches)
    token <- sub("[[:space:]].*$", "", token)
    suppressWarnings(as.numeric(gsub(",", "", token, fixed = TRUE)))
  }
  result$raw_read_pairs <- extract("^[[:space:]]*Total read pairs processed:")
  result$retained_read_pairs <- extract("^[[:space:]]*Pairs written \\(passing filters\\):")
  result
}

save_plot <- function(plot, stub, width = PLOT_WIDTH, height = PLOT_HEIGHT) {
  if (SHOW_PLOTS) 
    print(plot)
  for (extension in c("png", "pdf")) {
    ggplot2::ggsave(file.path(OUTPUT_DIR, paste0(stub, ".", extension)), plot = plot, width = width, 
      height = height, units = "in", dpi = DPI, bg = "white", limitsize = FALSE)
  }
}

run_pca <- function(counts, metadata, analysis) {
  unlink(file.path(OUTPUT_DIR, paste0("PCA_", analysis, c(".png", ".pdf"))))
  empty <- list(scores = NULL, variance = NULL, logcpm = NULL, normalization = NULL, status = data.frame(analysis = analysis, 
    status = "not computed", reason = "fewer than 3 valid samples", n_samples = ncol(counts), 
    n_guides = 0L))
  if (ncol(counts) < 2L) 
    return(empty)
  keep <- rowSums(counts >= PCA_MIN_COUNT) >= min(PCA_MIN_SAMPLES, ncol(counts))
  if (sum(keep) < 2L) {
    empty$status$reason <- "fewer than 2 guides passing diagnostic count filter"
    return(empty)
  }
  if (any(colSums(counts[keep, , drop = FALSE]) <= 0)) {
    empty$status$reason <- "a sample has zero counts after diagnostic feature filtering; review its depth/representation"
    return(empty)
  }
  y <- edgeR::DGEList(counts = counts[keep, , drop = FALSE])
  y <- edgeR::calcNormFactors(y, method = "TMM")
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = PCA_LOGCPM_PRIOR_COUNT, normalized.lib.sizes = TRUE)
  variances <- apply(logcpm, 1, stats::var)
  finite <- apply(logcpm, 1, function(x) all(is.finite(x))) & is.finite(variances) & variances > 
    0
  logcpm <- logcpm[finite, , drop = FALSE]
  variances <- variances[finite]
  if (nrow(logcpm) < 2L) {
    empty$status$reason <- "fewer than 2 finite, variable guides"
    return(empty)
  }
  if (ncol(counts) < 3L) {
    empty$logcpm <- logcpm
    return(empty)
  }
  selected <- order(-variances, rownames(logcpm))[seq_len(min(PCA_TOP_VARIABLE_GUIDES, nrow(logcpm)))]
  fit <- stats::prcomp(t(logcpm[selected, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = min(5L, 
    ncol(logcpm) - 1L, length(selected)))
  variance <- 100 * fit$sdev^2/sum(fit$sdev^2)
  score_columns <- c("sample_key", "batch", "condition", "experiment", "drug", "dose", "timing", 
    "rep", "group", "treatment")
  scores <- cbind(metadata[, score_columns, drop = FALSE], as.data.frame(fit$x), analysis = analysis)
  table <- data.frame(analysis = analysis, PC = paste0("PC", seq_along(variance)), variance_percent = variance)
  normalization <- data.frame(sample_key = colnames(counts), analysis = analysis, filtered_library_size = y$samples$lib.size, 
    TMM_factor = y$samples$norm.factors)
  color <- if (startsWith(analysis, "within_")) 
    "treatment"
  else "batch"
  plot <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, color = .data[[color]], shape = timing)) + 
    ggplot2::geom_point(size = 3) + ggplot2::scale_shape_manual(values = c(T0 = 15, T1 = 16, 
    T2 = 17), na.value = 4) + ggplot2::labs(title = gsub("_", " ", analysis), subtitle = paste(nrow(logcpm), 
    "eligible guides;", length(selected), "variable guides; independent TMM logCPM"), x = sprintf("PC1 (%.1f%%)", 
    variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2])) + ggplot2::theme_bw(base_size = 10, 
    base_family = FONT_FAMILY)
  if (LABEL_ALL_SAMPLES) {
    plot <- plot + ggrepel::geom_text_repel(ggplot2::aes(label = sample_key), size = 2.7, max.overlaps = Inf, 
      seed = PLOT_SEED, show.legend = FALSE)
  }
  save_plot(plot, paste0("PCA_", analysis))
  list(scores = scores, variance = table, logcpm = logcpm, normalization = normalization, status = data.frame(analysis = analysis, 
    status = "computed", reason = "", n_samples = ncol(counts), n_guides = length(selected)))
}

# LOAD AND VALIDATE FILES
if (!file.exists("AMG_cleaned_Git.Rproj")) stop("Open AMG_cleaned_Git.Rproj and run from its working directory.", 
  call. = FALSE)

COUNT_DIR <- file.path(COUNT_DIR_RELATIVE)

REFERENCE_FILE <- file.path(REFERENCE_FILE_RELATIVE)

OUTPUT_DIR <- file.path(OUTPUT_DIR_RELATIVE)

COUNT_FILES <- setNames(file.path(COUNT_DIR, paste0(BATCHES, "_overall_stats.tsv")), BATCHES)

DESIGN_FILES <- setNames(file.path(COUNT_DIR, paste0(BATCHES, "_experimental_design.tsv")), BATCHES)

required_files <- c(REFERENCE_FILE, COUNT_FILES, DESIGN_FILES)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files)) stop("Missing inputs:\n", paste(missing_files, collapse = "\n"), call. = FALSE)

if (!COUNT_ID_MODE %in% c("auto", "fasta_id", "sequence")) stop("Invalid COUNT_ID_MODE.")

positive_settings <- c(MIN_ASSIGNED_COUNTS, RAREFACTION_DEPTH, PCA_MIN_COUNT, PCA_MIN_SAMPLES, 
  PCA_TOP_VARIABLE_GUIDES, PCA_LOGCPM_PRIOR_COUNT, PLOT_WIDTH, PLOT_HEIGHT, DPI, REPRESENTATION_COUNT_LEVELS)

if (any(!is.finite(positive_settings)) || any(positive_settings <= 0)) stop("Count, depth and plot settings must be finite and positive.")

if (any(c(RAREFACTION_DEPTH, PCA_MIN_SAMPLES, PCA_TOP_VARIABLE_GUIDES, REPRESENTATION_COUNT_LEVELS)%%1 != 
  0)) stop("Rarefaction depth, guide cutoffs and PCA feature/sample counts must be integers.")

if (!is.finite(MIN_REPLICATE_SPEARMAN) || MIN_REPLICATE_SPEARMAN < -1 || MIN_REPLICATE_SPEARMAN > 
  1) stop("Spearman review cutoff must be between -1 and 1.")

if (PCA_TOP_VARIABLE_GUIDES < 2L) stop("PCA_TOP_VARIABLE_GUIDES must be at least 2 for a two-axis PCA plot.")

percent_settings <- c(MIN_T0_GUIDE_COVERAGE_PERCENT, MAX_T0_TOP_GUIDE_PERCENT, MAX_UNKNOWN_COUNT_PERCENT)

if (any(!is.finite(percent_settings)) || any(percent_settings < 0 | percent_settings > 100)) stop("Percentage review cutoffs must be between 0 and 100.")

ratio_settings <- c(MIN_WITHIN_GROUP_COUNT_RATIO, MIN_WITHIN_GROUP_COVERAGE_RATIO)

if (any(!is.finite(ratio_settings)) || any(ratio_settings < 0 | ratio_settings > 1)) stop("Relative review cutoffs must be between 0 and 1.")

for (setting in c("LABEL_ALL_SAMPLES", "EXPORT_TSV", "SHOW_PLOTS")) {
  value <- get(setting)
  if (!is.logical(value) || length(value) != 1L || is.na(value)) 
    stop(setting, " must be TRUE or FALSE.")
}

if (anyDuplicated(BATCHES) || !length(BATCHES)) stop("BATCHES must contain unique identifiers.")

reference <- read_reference(REFERENCE_FILE)

designs <- setNames(lapply(BATCHES, function(batch) read_design(DESIGN_FILES[[batch]], batch)), 
  BATCHES)

qc_list <- matrices <- input_audit <- list()

chosen_mode <- NULL

for (batch in BATCHES) {
  message("Reading raw counts and computing sample QC: ", batch)
  counts <- read_counts(COUNT_FILES[[batch]])
  design <- designs[[batch]]
  observed <- unique(counts$barcode)
  mode <- COUNT_ID_MODE
  if (mode == "auto") {
    id_hits <- mean(observed %in% reference$fasta_id)
    sequence_hits <- mean(toupper(observed) %in% reference$sequence)
    if (max(id_hits, sequence_hits) == 0) 
      stop(batch, ": count IDs match neither FASTA IDs nor sequences. Supply the matching full-library reference.")
    mode <- if (identical(reference$fasta_id, reference$sequence)) 
      "sequence"
    else if (id_hits >= sequence_hits) 
      "fasta_id"
    else "sequence"
  }
  if (!is.null(chosen_mode) && chosen_mode != mode) 
    stop("Count identifier modes differ between batches. Resolve IDs before cross-batch QC.")
  chosen_mode <- mode
  reference_ids <- if (mode == "fasta_id") 
    reference$fasta_id
  else unique(reference$sequence)
  if (mode == "sequence") 
    counts$barcode <- toupper(counts$barcode)
  duplicate_rows <- counts[, .(duplicate_rows = .N - 1L), by = .(barcode, condition)]
  duplicate_rows <- duplicate_rows[, .(duplicate_rows = sum(duplicate_rows)), by = condition]
  counts <- counts[, .(count = sum(count)), by = .(barcode, condition)]
  counts$known_reference <- counts$barcode %in% reference_ids
  samples <- c(design$condition, setdiff(unique(counts$condition), design$condition))
  matrix <- matrix(0, nrow = length(reference_ids), ncol = length(samples), dimnames = list(reference_ids, 
    samples))
  known <- counts[counts$known_reference]
  matrix[cbind(match(known$barcode, reference_ids), match(known$condition, samples))] <- known$count
  totals <- counts[, .(table_assigned_counts = sum(count), unknown_assigned_counts = sum(count[!known_reference])), 
    by = condition]
  metadata <- design[match(samples, design$condition), , drop = FALSE]
  metadata$condition <- samples
  metadata$batch <- batch
  metadata$sample_key <- colnames(matrix)
  metadata$in_current_design <- samples %in% design$condition
  rows <- lapply(seq_along(samples), function(i) {
    sample <- samples[i]
    present <- sample %in% counts$condition
    total_i <- totals[match(sample, totals$condition), , drop = FALSE]
    metrics <- sample_complexity(matrix[, i])
    if (!present) {
      metrics <- lapply(metrics, function(x) NA_real_)
      metrics$top_guide_id <- NA_character_
    }
    processed <- read_cutadapt_pairs(batch, sample)
    data.frame(metadata[i, , drop = FALSE], count_rows_present = present, table_assigned_counts = total_i$table_assigned_counts, 
      unknown_assigned_counts = total_i$unknown_assigned_counts, unknown_count_percent = if (is.finite(total_i$table_assigned_counts) && 
        total_i$table_assigned_counts > 0) 
        100 * total_i$unknown_assigned_counts/total_i$table_assigned_counts
      else NA_real_, n_reference_guides = length(reference_ids), as.data.frame(metrics), 
      as.data.frame(processed), duplicate_count_rows = duplicate_rows$duplicate_rows[match(sample, 
        duplicate_rows$condition)], count_id_mode = mode, stringsAsFactors = FALSE)
  })
  qc_list[[batch]] <- data.table::rbindlist(rows, fill = TRUE)
  matrices[[batch]] <- matrix
  input_audit[[batch]] <- data.frame(batch = batch, count_id_mode = mode, n_reference_guides = length(reference_ids), 
    n_FASTA_records = nrow(reference), n_unique_reference_sequences = length(unique(reference$sequence)), 
    n_samples_in_design = nrow(design), n_samples_in_counts = length(unique(counts$condition)), 
    n_extra_count_samples = sum(!samples %in% design$condition), count_file_md5 = unname(tools::md5sum(COUNT_FILES[[batch]])), 
    design_file_md5 = unname(tools::md5sum(DESIGN_FILES[[batch]])))
  rm(counts, known)
}

qc <- as.data.frame(data.table::rbindlist(qc_list, fill = TRUE))

# Ordering of matrix columns must always match the QC rows.
raw_matrix <- do.call(cbind, matrices)

stopifnot(identical(colnames(raw_matrix), qc$sample_key))

qc$pca_eligible <- qc$in_current_design & qc$count_rows_present & is.finite(qc$reference_assigned_counts) & 
  qc$reference_assigned_counts > 0

qc$pca_note <- ifelse(!qc$in_current_design, "Outside current design", ifelse(!qc$count_rows_present, 
  "Missing count rows", ifelse(!qc$pca_eligible, "No reference-assigned counts", "Included regardless of review flags")))

qc$median_same_group_spearman <- NA_real_

qc$n_same_group_peers <- 0L

qc$count_ratio_to_group_median <- NA_real_

qc$coverage_ratio_to_group_median <- NA_real_

qc$group_sample_count <- NA_integer_

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# PCA AND WITHIN-CONDITION SIMILARITY
pca_results <- list()

for (batch in BATCHES) {
  unlink(file.path(OUTPUT_DIR, paste0("QC_similarity_", batch, c(".png", ".pdf"))))
  rows <- which(qc$batch == batch & qc$pca_eligible)
  result <- run_pca(raw_matrix[, rows, drop = FALSE], qc[rows, , drop = FALSE], paste0("within_", 
    batch))
  if (!is.null(result$logcpm)) {
    correlations <- stats::cor(result$logcpm, method = "spearman")
    for (i in seq_along(rows)) {
      peers <- which(qc$group[rows] == qc$group[rows[i]] & seq_along(rows) != i)
      qc$n_same_group_peers[rows[i]] <- length(peers)
      values <- correlations[i, peers]
      if (any(is.finite(values))) 
        qc$median_same_group_spearman[rows[i]] <- median(values[is.finite(values)])
    }
    heatmap <- as.data.frame(as.table(correlations), stringsAsFactors = FALSE)
    names(heatmap) <- c("sample_x", "sample_y", "spearman")
    plot <- ggplot2::ggplot(heatmap, ggplot2::aes(sample_x, sample_y, fill = spearman)) + ggplot2::geom_tile() + 
      ggplot2::scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, 
        limits = c(-1, 1), na.value = "grey80") + ggplot2::coord_equal() + ggplot2::labs(title = paste("Sample similarity:", 
      toupper(batch)), subtitle = "Spearman correlation on variable, count-filtered TMM logCPM guides", 
      x = NULL, y = NULL) + ggplot2::theme_minimal(base_size = 10, base_family = FONT_FAMILY) + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 7), 
        axis.text.y = ggplot2::element_text(size = 7))
    save_plot(plot, paste0("QC_similarity_", batch), width = max(PLOT_WIDTH, 4 + 0.2 * length(rows)), 
      height = max(PLOT_HEIGHT, 4 + 0.2 * length(rows)))
  }
  result$logcpm <- NULL
  pca_results[[paste0("within_", batch)]] <- result
}

subsets <- list(global_all_samples = which(qc$pca_eligible), global_untreated_T0 = which(qc$pca_eligible & 
  qc$drug == "None" & qc$timing == "T0"), global_untreated_postinduction = which(qc$pca_eligible & 
  qc$drug == "None" & qc$timing %in% c("T1", "T2")))

for (analysis in names(subsets)) {
  rows <- subsets[[analysis]]
  result <- run_pca(raw_matrix[, rows, drop = FALSE], qc[rows, , drop = FALSE], analysis)
  result$logcpm <- NULL
  pca_results[[analysis]] <- result
}

# REVIEW FLAGS (NO AUTOMATIC SAMPLE REMOVAL)
flag_rows <- list()

add_flag <- function(row, severity, issue, detail) {
  flag_rows[[length(flag_rows) + 1L]] <<- data.frame(sample_key = qc$sample_key[row], batch = qc$batch[row], 
    condition = qc$condition[row], severity = severity, issue = issue, detail = detail, stringsAsFactors = FALSE)
}

for (i in seq_len(nrow(qc))) {
  if (!qc$count_rows_present[i]) {
    add_flag(i, "Missing data", "missing_count_rows", "Design sample has no rows in the count table; no measured library metrics or PCA coordinates.")
    next
  }
  if (qc$reference_assigned_counts[i] == 0) 
    add_flag(i, "Review", "zero_reference_counts", "No counts match the reference; this sample cannot enter normalization.")
  if (qc$table_assigned_counts[i] < MIN_ASSIGNED_COUNTS) 
    add_flag(i, "Review", "low_assigned_counts", paste("Assigned counts below configured threshold", 
      MIN_ASSIGNED_COUNTS))
  if (is.finite(qc$unknown_count_percent[i]) && qc$unknown_count_percent[i] > MAX_UNKNOWN_COUNT_PERCENT) 
    add_flag(i, "Review", "unknown_count_ids", "An excessive fraction of assigned counts does not match the library reference.")
  if (qc$duplicate_count_rows[i] > 0) 
    add_flag(i, "Review", "duplicated_count_rows", "Repeated barcode/sample rows were summed, as in the metacounter. Check concatenation or lane pooling.")
  if (qc$noninteger_reference_counts[i] > 0) 
    add_flag(i, "Information", "noninteger_counts", "Fractional supplied counts preserved; integer-count rarefaction was not calculated.")
  if (is.na(qc$raw_read_pairs[i])) 
    add_flag(i, "Information", "raw_read_depth_unavailable", "Raw sequencing-read pairs cannot be inferred from overall_stats; no usable cutadapt report was found.")
  if (!qc$in_current_design[i]) {
    add_flag(i, "Information", "outside_current_design", "Count-table sample retained in the QC inventory; omitted from PCA because it has no current design row.")
    next
  }
  if (qc$drug[i] == "None" && qc$timing[i] == "T0") {
    if (qc$detected_guides_percent[i] < MIN_T0_GUIDE_COVERAGE_PERCENT) 
      add_flag(i, "Review", "low_T0_representation", "Untreated T0 detected-guide fraction is below the configured review threshold.")
    if (is.finite(qc$top_guide_percent[i]) && qc$top_guide_percent[i] > MAX_T0_TOP_GUIDE_PERCENT) 
      add_flag(i, "Review", "dominant_T0_guide", "One guide dominates the untreated T0 library; review its identifier and counts.")
  }
  peers <- which(qc$batch == qc$batch[i] & qc$group == qc$group[i] & qc$pca_eligible)
  qc$group_sample_count[i] <- length(peers)
  if (length(peers) >= 3L && qc$pca_eligible[i]) {
    depth_median <- median(qc$reference_assigned_counts[peers])
    coverage_median <- median(qc$detected_guides_percent[peers])
    qc$count_ratio_to_group_median[i] <- qc$reference_assigned_counts[i]/depth_median
    if (coverage_median > 0) 
      qc$coverage_ratio_to_group_median[i] <- qc$detected_guides_percent[i]/coverage_median
    if (qc$count_ratio_to_group_median[i] < MIN_WITHIN_GROUP_COUNT_RATIO) 
      add_flag(i, "Review", "low_depth_vs_group", "Assigned depth is low relative to same-batch, same-treatment, same-timepoint samples.")
    if (is.finite(qc$coverage_ratio_to_group_median[i]) && qc$coverage_ratio_to_group_median[i] < 
      MIN_WITHIN_GROUP_COVERAGE_RATIO) 
      add_flag(i, "Review", "low_representation_vs_group", "Guide detection is low relative to the same biological group; review depth before attributing a bottleneck.")
  }
  if (is.finite(qc$median_same_group_spearman[i]) && qc$median_same_group_spearman[i] < MIN_REPLICATE_SPEARMAN) 
    add_flag(i, "Review", "low_same_group_similarity", "Low same-condition similarity. With only two replicates, both are flagged; this cannot identify which sample is responsible.")
  if (qc$n_same_group_peers[i] == 0L) 
    add_flag(i, "Information", "replicate_similarity_unavailable", "No usable same-condition peer correlation was available; consult PCA status and the design.")
}

flags <- if (length(flag_rows)) as.data.frame(data.table::rbindlist(flag_rows)) else data.frame(sample_key = character(), 
  batch = character(), condition = character(), severity = character(), issue = character(), 
  detail = character())

qc$review_flag_count <- vapply(qc$sample_key, function(key) sum(flags$sample_key == key & flags$severity %in% 
  c("Review", "Missing data")), integer(1))

qc$review_status <- ifelse(!qc$count_rows_present, "Missing data", ifelse(qc$review_flag_count > 
  0L, "Review", "No configured review flag"))

qc$review_reasons <- vapply(qc$sample_key, function(key) paste(unique(flags$issue[flags$sample_key == 
  key & flags$severity %in% c("Review", "Missing data")]), collapse = "; "), character(1))

front_columns <- c("sample_key", "batch", "condition", "in_current_design", "review_status", "review_reasons", 
  "table_assigned_counts", "reference_assigned_counts", "raw_read_pairs", "retained_read_pairs", 
  "detected_guides_percent", "shannon_effective_guides", "top_guide_percent", "top_guide_id", 
  "median_same_group_spearman")

qc <- qc[, c(front_columns, setdiff(names(qc), front_columns)), drop = FALSE]

# QC PLOTS AND REPORT EXPORT
qc$plot_label <- factor(qc$sample_key, levels = rev(qc$sample_key))

qc$sample_scope <- ifelse(qc$in_current_design, "In design", "Outside design")

metric_specs <- list(assigned_depth = list(column = "table_assigned_counts", label = "Assigned barcode counts (supplied units)", 
  log_scale = TRUE), detected_guides = list(column = "detected_guides_percent", label = "Detected reference guides (%)", 
  log_scale = FALSE), effective_guides = list(column = "shannon_effective_guides", label = "Shannon effective number of guides", 
  log_scale = FALSE), dominant_guide = list(column = "top_guide_percent", label = "Counts assigned to the top reference guide (%)", 
  log_scale = FALSE))

for (stub in names(metric_specs)) {
  spec <- metric_specs[[stub]]
  plot <- ggplot2::ggplot(qc, ggplot2::aes(x = .data[[spec$column]], y = plot_label, color = review_status, 
    shape = sample_scope)) + ggplot2::geom_point(size = 2, na.rm = TRUE) + ggplot2::facet_grid(batch ~ 
    ., scales = "free_y", space = "free_y") + ggplot2::scale_color_manual(values = c(`Missing data` = "#000000", 
    Review = "#D95F02", `No configured review flag` = "#1B9E77")) + ggplot2::labs(x = spec$label, 
    y = NULL, color = "Review status", shape = "Design membership", subtitle = "Review flags guide inspection; the script does not exclude flagged samples") + 
    ggplot2::theme_bw(base_size = 10, base_family = FONT_FAMILY) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
  if (spec$log_scale) 
    plot <- plot + ggplot2::scale_x_continuous(trans = "log1p")
  save_plot(plot, paste0("QC_", stub), height = max(PLOT_HEIGHT, 3 + 0.18 * nrow(qc)))
}

qc$plot_label <- NULL

qc$sample_scope <- NULL

bind_component <- function(component) {
  pieces <- lapply(pca_results, function(result) result[[component]])
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (!length(pieces)) 
    return(data.frame(note = "No analysis available; see PCA_status."))
  as.data.frame(data.table::rbindlist(pieces, fill = TRUE))
}

notes <- data.frame(topic = c("Purpose", "Assigned counts", "Raw read pairs", "Reference universe", 
  "Complexity", "Rarefaction", "PCA", "Review flags", "Exclusion", "Validation"), explanation = c("Pre-model sample review using unfiltered barcode-count tables.", 
  "Totals retain supplied count units; they are not raw sequencing depth.", "Populated from cutadapt reports where available; otherwise NA. No assigned-count/raw-read mapping fraction is inferred.", 
  "Uses the full barcode FASTA. FASTA ID mode counts records; sequence mode collapses duplicate reference sequences. Inspect Input_audit.", 
  "Detection and diversity describe guide representation, not true cellular clone numbers. Postinduction/treatment changes can be biological.", 
  "Expected guides detected at a fixed draw depth, without replacement. NA for fractional counts or samples below that depth.", 
  "Separate global, baseline-control, postinduction-control and within-batch PCA. Count-filtered TMM logCPM; top-variable centered unscaled features. No batch correction.", 
  "Default numerical cutoffs are editable screening heuristics, not established quality standards. Compare same-condition replicates.", 
  "No designs are edited. Review flags never remove samples. Samples outside the current design remain in the QC inventory but do not enter PCA.", 
  "New diagnostic workflow. Actual dataset/package execution must be validated before treating the report as authoritative."), 
  stringsAsFactors = FALSE)

settings_names <- c("COUNT_ID_MODE", "MIN_ASSIGNED_COUNTS", "MIN_T0_GUIDE_COVERAGE_PERCENT", "MAX_T0_TOP_GUIDE_PERCENT", 
  "MAX_UNKNOWN_COUNT_PERCENT", "MIN_WITHIN_GROUP_COUNT_RATIO", "MIN_WITHIN_GROUP_COVERAGE_RATIO", 
  "MIN_REPLICATE_SPEARMAN", "REPRESENTATION_COUNT_LEVELS", "RAREFACTION_DEPTH", "PCA_MIN_COUNT", 
  "PCA_MIN_SAMPLES", "PCA_TOP_VARIABLE_GUIDES", "PCA_LOGCPM_PRIOR_COUNT", "LABEL_ALL_SAMPLES", 
  "EXPORT_TSV", "PLOT_SEED")

settings <- data.frame(setting = c(settings_names, "reference_md5", "R_version", REQUIRED_PACKAGES), 
  value = c(vapply(settings_names, function(name) paste(get(name), collapse = ", "), character(1)), 
    unname(tools::md5sum(REFERENCE_FILE)), R.version.string, vapply(REQUIRED_PACKAGES, function(package) as.character(utils::packageVersion(package)), 
      character(1))))

sheets <- list(Readme = notes, Samples = qc, Flags = flags, PCA_status = bind_component("status"), 
  PCA_scores = bind_component("scores"), PCA_variance = bind_component("variance"), PCA_normalization = bind_component("normalization"), 
  Input_audit = as.data.frame(data.table::rbindlist(input_audit)), Settings = settings)

workbook <- openxlsx::createWorkbook()

header_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold", 
  fgFill = "#1F4E78", fontColour = "#FFFFFF", wrapText = TRUE)

wrap_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, wrapText = TRUE, valign = "top")

for (sheet in names(sheets)) {
  openxlsx::addWorksheet(workbook, sheet)
  openxlsx::writeData(workbook, sheet, sheets[[sheet]], keepNA = TRUE, headerStyle = header_style)
  openxlsx::freezePane(workbook, sheet, firstRow = TRUE, firstCol = TRUE)
  openxlsx::setRowHeights(workbook, sheet, rows = 1, heights = 40)
  openxlsx::setColWidths(workbook, sheet, cols = seq_len(ncol(sheets[[sheet]])), widths = 18)
  if (nrow(sheets[[sheet]])) 
    openxlsx::addFilter(workbook, sheet, rows = 1, cols = seq_len(ncol(sheets[[sheet]])))
  if (EXPORT_TSV) 
    data.table::fwrite(sheets[[sheet]], file.path(OUTPUT_DIR, paste0("QC_", sheet, ".tsv")), 
      sep = "\t", na = "NA")
}

openxlsx::setColWidths(workbook, "Readme", cols = 2, widths = 100)

openxlsx::setColWidths(workbook, "Flags", cols = 6, widths = 100)

openxlsx::setColWidths(workbook, "Samples", cols = which(names(qc) %in% c("sample_key", "condition", 
  "review_reasons")), widths = 40)

openxlsx::addStyle(workbook, "Readme", wrap_style, rows = seq_len(nrow(notes)) + 1L, cols = 2, 
  gridExpand = TRUE)

openxlsx::setRowHeights(workbook, "Readme", rows = seq_len(nrow(notes)) + 1L, heights = 45)

if (nrow(flags)) {
  openxlsx::addStyle(workbook, "Flags", wrap_style, rows = seq_len(nrow(flags)) + 1L, cols = 6, 
    gridExpand = TRUE)
  openxlsx::setRowHeights(workbook, "Flags", rows = seq_len(nrow(flags)) + 1L, heights = 32)
}

openxlsx::saveWorkbook(workbook, file.path(OUTPUT_DIR, "sample_QC_report.xlsx"), overwrite = TRUE)

message("QC/PCA completed. Review samples: ", sum(qc$review_flag_count > 0L), "/", nrow(qc), ". Designs were not changed.\nReport: ", 
  file.path(OUTPUT_DIR, "sample_QC_report.xlsx"))

