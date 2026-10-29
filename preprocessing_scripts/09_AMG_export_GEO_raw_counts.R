# Purpose: Assemble a GEO-format guide-count matrix from the same unfiltered
# barcode-count tables and guide key used by 02_AMG_metacounter.R.
# Inputs: Barcode_counting/count_tables/b*_overall_stats.tsv and
# b*_experimental_design.tsv and annotation_input/PAO1_key.tsv. Sample names
# use batch-prefixed IDs (e.g. B1_1-T0_MH). GEO column labels are generated
# from the designs; reference headers define the manuscript subset
# and order when the comparison matrix is available. No crosswalk is required.
# Counts are aggregated by spacer/sample, without T0, fitness or read filters.
# All unique key spacers are retained; first key assignments define primary loci.
# Within a batch, missing sample counts for an observed spacer become zero.
# A spacer absent from all selected samples in that batch retains NA, matching
# the left-joined layout of the deposited matrix.
# Package: data.table (also required by the metacounter).
# Metadata TSVs may contain quoted or unquoted fields.
# Run from AMG_cleaned_Git.Rproj:
# source("preprocessing_scripts/09_AMG_export_GEO_raw_counts.R")
# Outputs: A compressed TSV matrix, Sample_key.tsv and a small export audit.
# USER SETTINGS
BATCHES <- c("b1", "b2", "b3", "b4")

COUNT_DIR <- "Barcode_counting/count_tables"

KEY_FILE <- "annotation_input/PAO1_key.tsv"

OUTPUT_DIR <- "preprocessing_scripts/out_09_AMG_export_GEO_raw_counts"

OUTPUT_FILENAME <- "CRISPRi_guide_raw_counts.tsv.gz"

# The deposited file is an optional independent check, never a source of counts.
REFERENCE_GEO_FILE <- "Barcode_counting/count_tables/GSE349013_CRISPRi_guide_raw_counts.tsv.gz"

REQUIRE_REFERENCE_COUNT_MATCH <- TRUE

# Local functions for this analysis.
main <- function() {
  if (!requireNamespace("data.table", quietly = TRUE)) 
    stop("Required package: data.table")
  paths <- c(KEY_FILE, file.path(COUNT_DIR, paste0(BATCHES, "_overall_stats.tsv")), file.path(COUNT_DIR, 
    paste0(BATCHES, "_experimental_design.tsv")))
  if (any(!file.exists(paths))) 
    stop("Missing metacounter inputs:\n", paste(paths[!file.exists(paths)], collapse = "\n"))
  read_tsv <- function(path) {
    connection <- if (endsWith(path, ".gz")) 
      gzfile(path, "rt")
    else file(path, "rt")
    on.exit(close(connection))
    utils::read.delim(connection, check.names = FALSE, stringsAsFactors = FALSE, quote = "\"", 
      comment.char = "", na.strings = c("", "NA"))
  }
  required <- function(data, columns, label) {
    if (length(setdiff(columns, names(data)))) 
      stop(label, " lacks: ", paste(setdiff(columns, names(data)), collapse = ", "))
  }
  # 1. Define the exported guide universe from the canonical key.
  # Require valid 20-nt spacers and guide classes. For repeated locus assignments,
  # retain the first assignment so each spacer occupies exactly one matrix row.
  key <- read_tsv(KEY_FILE)
  required(key, c("spacer", "type", "locus_tag"), "Guide key")
  key <- unique(key[c("spacer", "type", "locus_tag")])
  key$spacer <- toupper(gsub("\\s+", "", key$spacer))
  if (anyNA(key$spacer) || any(!grepl("^[ACGT]{20}$", key$spacer)) || anyNA(key$type) || any(!key$type %in% 
    c("perfect", "mismatch", "control"))) 
    stop("Invalid guide key.")
  type_counts <- tapply(key$type, key$spacer, function(value) length(unique(value)))
  if (any(type_counts != 1L)) 
    stop("A spacer has conflicting guide types in the key.")
  if (anyDuplicated(key$spacer)) 
    message("Retaining the first primary assignment for ", sum(duplicated(key$spacer)), " additional key rows.")
  key <- key[!duplicated(key$spacer), , drop = FALSE]
  if (anyNA(key$locus_tag[key$type != "control"])) 
    stop("A targeted guide lacks its primary locus.")
  # 2. Translate each design row into a GEO column label using batch, sample
  # number, time point and treatment/dose. Retain the exact input sample ID too,
  # so count-table records can be matched without guessing names.
  sample_list <- lapply(BATCHES, function(batch) {
    design <- read_tsv(file.path(COUNT_DIR, paste0(batch, "_experimental_design.tsv")))
    required(design, c("condition", "drug", "dose", "timing", "rep"), paste(batch, "design"))
    condition <- trimws(as.character(design$condition))
    if (anyNA(condition) || anyDuplicated(condition) || any(!nzchar(condition))) 
      stop("Missing or duplicate design conditions in ", batch)
    short_condition <- sub("^[bB][0-9]+[_-]", "", condition)
    if (any(!grepl("^[0-9]+", short_condition))) 
      stop("Cannot extract sample numbers from the ", batch, " design conditions.")
    number <- as.integer(sub("^([0-9]+).*$", "\\1", short_condition))
    drug <- toupper(trimws(as.character(design$drug)))
    drug[is.na(drug) | drug == ""] <- "NONE"
    dose <- suppressWarnings(as.numeric(gsub(",", ".", as.character(design$dose), fixed = TRUE)))
    dose[drug == "NONE"] <- 0
    timing <- toupper(trimws(as.character(design$timing)))
    if (anyNA(dose) || any(!is.finite(dose)) || anyNA(timing) || any(!grepl("^T[0-9]+$", timing)) || 
      anyNA(design$rep)) 
      stop("Incomplete treatment metadata in ", batch)
    dose_label <- vapply(dose, function(value) format(value, scientific = FALSE, trim = TRUE), 
      character(1))
    treatment <- paste0(drug, dose_label)
    treatment[drug == "NONE"] <- "MH"
    treatment[drug == "FBS"] <- "FBS"
    drug[drug == "NONE"] <- "None"
    data.frame(label = paste(toupper(batch), sprintf("%02d", number), timing, treatment, sep = "_"), 
      batch = toupper(batch), original_condition = condition, original_sample_id = condition,
      drug = drug, dose = dose, timing = timing, rep = design$rep, stringsAsFactors = FALSE)
  })
  samples <- do.call(rbind, sample_list)
  if (anyDuplicated(samples$label) || anyDuplicated(samples$original_sample_id)) 
    stop("Designs produce ambiguous sample labels.")
  # 3. If a deposited comparison matrix is available, use its sample subset
  # and order. It is a reference for checking, not the source of exported counts.
  reference <- if (file.exists(REFERENCE_GEO_FILE)) 
    read_tsv(REFERENCE_GEO_FILE)
  else NULL
  if (!is.null(reference)) {
    if (!identical(names(reference)[1:3], c("spacer", "type", "locus_tag"))) 
      stop("Invalid comparison-matrix metadata columns.")
    reference_labels <- names(reference)[-(1:3)]
    if (!length(reference_labels) || anyDuplicated(reference_labels)) 
      stop("Invalid comparison-matrix sample headers.")
    selected <- match(reference_labels, samples$label)
    if (anyNA(selected)) 
      stop("Experimental designs do not reproduce these comparison-file headers:\n", paste(reference_labels[is.na(selected)], 
        collapse = "\n"), "\nReview timing, treatment and dose metadata; labels were not guessed.")
    samples <- samples[selected, , drop = FALSE]
  }
  if (!setequal(tolower(samples$batch), BATCHES)) 
    stop("Selected samples do not cover the requested batches.")
  matrix <- matrix(NA_real_, nrow(key), nrow(samples), dimnames = list(NULL, samples$label))
  audit <- vector("list", nrow(samples))
  # 4. Fill the guide-by-sample matrix separately for each batch from raw counts.
  # Sum repeated spacer/sample records; no T0 abundance or fitness filter is used.
  for (batch in BATCHES) {
    indices <- which(tolower(samples$batch) == batch)
    selected <- samples[indices, , drop = FALSE]
    counts <- data.table::fread(file.path(COUNT_DIR, paste0(batch, "_overall_stats.tsv")), 
      sep = "\t", header = FALSE, colClasses = "character", fill = TRUE, data.table = TRUE)
    if (ncol(counts) < 3L) 
      stop("Expected spacer, count and condition columns for ", batch)
    counts <- counts[, 1:3]
    data.table::setnames(counts, c("spacer", "count", "condition"))
    counts <- counts[!(spacer %in% c("barcode", "spacer") & count %in% c("count", "counts") & 
      condition %in% c("condition", "sample"))]
    counts$spacer <- toupper(gsub("\\s+", "", counts$spacer))
    counts$condition <- trimws(counts$condition)
    counts$count <- suppressWarnings(as.numeric(counts$count))
    if (anyNA(counts) || any(!is.finite(counts$count)) || any(counts$count < 0) || any(counts$count%%1 != 
      0)) 
      stop("Invalid raw barcode counts in ", batch)
    if (!all(selected$original_condition %in% counts$condition)) 
      stop("A selected sample has no records in the ", batch, " raw count table.")
    counts <- counts[counts$condition %in% selected$original_condition]
    original_n <- table(factor(counts$condition, levels = selected$original_condition))
    counts <- counts[, .(count = sum(count)), by = .(spacer, condition)]
    unique_n <- table(factor(counts$condition, levels = selected$original_condition))
    rows <- match(counts$spacer, key$spacer)
    known <- !is.na(rows)
    present <- unique(rows[known])
    if (!length(present)) 
      stop("No raw spacers match the guide key in ", batch)
    # A spacer seen anywhere in this batch receives zero for missing sample records.
    # A key spacer absent throughout the batch stays NA, preserving the deposited
    # left-join convention. Reads outside the canonical key are counted in the audit.
    matrix[present, indices] <- 0
    matrix[cbind(rows[known], indices[match(counts$condition[known], selected$original_condition)])] <- counts$count[known]
    for (i in seq_along(indices)) {
      sample <- selected$original_condition[i]
      audit[[indices[i]]] <- data.frame(label = selected$label[i], batch = toupper(batch), 
        original_sample_id = selected$original_sample_id[i], assigned_reads_in_key = sum(counts$count[counts$condition == 
         sample & known]), assigned_reads_outside_key = sum(counts$count[counts$condition == 
         sample & !known]), guides_observed_in_batch = length(present), missing_guide_entries = sum(is.na(matrix[, 
         indices[i]])), duplicate_count_rows_merged = as.integer(original_n[i] - unique_n[i]))
    }
  }
  audit <- do.call(rbind, audit)
  result <- cbind(key, as.data.frame(matrix, check.names = FALSE))
  audit$reference_checked <- !is.null(reference)
  audit$reference_count_differences <- NA_integer_
  audit$reference_NA_pattern_differences <- NA_integer_
  # 5. Compare guide annotations, counts and missing-value patterns with the
  # reference. Numerical comparison treats NA as zero; NA-pattern differences
  # are reported separately in the audit.
  if (!is.null(reference)) {
    required(reference, c("spacer", "type", "locus_tag", samples$label), "GEO reference")
    if (anyDuplicated(reference$spacer) || !setequal(reference$spacer, key$spacer) || !identical(names(reference), 
      names(result))) 
      stop("Guide universe or headers differ from the GEO reference.")
    reference <- reference[match(key$spacer, reference$spacer), , drop = FALSE]
    same_locus <- (is.na(key$locus_tag) & is.na(reference$locus_tag)) | (!is.na(key$locus_tag) & 
      !is.na(reference$locus_tag) & key$locus_tag == reference$locus_tag)
    if (any(key$type != reference$type) || any(!same_locus)) 
      stop("Primary guide annotations differ from the GEO reference.")
    other <- as.matrix(reference[samples$label])
    audit$reference_NA_pattern_differences <- colSums(is.na(matrix) != is.na(other))
    generated <- matrix
    generated[is.na(generated)] <- 0
    other[is.na(other)] <- 0
    audit$reference_count_differences <- colSums(generated != other)
  }
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  # Stop before replacing completed outputs if numerical counts disagree and
  # the reference-match requirement is enabled. Retain an audit of that mismatch.
  if (REQUIRE_REFERENCE_COUNT_MATCH && any(audit$reference_count_differences > 0, na.rm = TRUE)) {
    failure_audit <- file.path(OUTPUT_DIR, "raw_count_export_reference_mismatch.tsv")
    utils::write.table(audit, failure_audit, sep = "\t", quote = FALSE, row.names = FALSE, 
      na = "NA")
    stop("Generated counts differ from the comparison matrix. Previous completed outputs were retained. Review ", 
      failure_audit)
  }
  # 6. Export the count matrix, sample metadata and per-sample audit.
  # These are raw guide counts; no library-size normalization is applied.
  output <- file.path(OUTPUT_DIR, OUTPUT_FILENAME)
  connection <- if (endsWith(output, ".gz")) 
    gzfile(output, "wt")
  else file(output, "wt")
  tryCatch(utils::write.table(result, connection, sep = "\t", quote = FALSE, row.names = FALSE, 
    na = "NA"), finally = close(connection))
  sample_key <- samples[c("label", "batch", "drug", "dose", "timing", "rep", "original_sample_id")]
  names(sample_key)[names(sample_key) == "drug"] <- "drug or treatment"
  utils::write.table(sample_key, file.path(OUTPUT_DIR, "Sample_key.tsv"), sep = "\t", quote = FALSE, 
    row.names = FALSE, na = "NA")
  utils::write.table(audit, file.path(OUTPUT_DIR, "raw_count_export_audit.tsv"), sep = "\t", 
    quote = FALSE, row.names = FALSE, na = "NA")
  message("Saved GEO-format raw counts: ", output, "\n", nrow(result), " unique guides; ", nrow(samples), 
    " samples. No T0 abundance or fitness filter was applied.")
  invisible(result)
}

main()

