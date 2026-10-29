# Join unfiltered barcode counts to guide context for Figure S8 inputs.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/25_AMG_guide_efficacy_inputs.R")
# Inputs:
# - Barcode_counting/count_tables/b*_overall_stats.tsv: headerless spacer, count, sample columns.
# - Matching b*_experimental_design.tsv and annotation_input/PAO1_key.tsv: unique spacer with
#   locus_tag, PAM, downstream context, offset and off_target_n.
# Upstream / execution order:
# - Use the supplied counts/designs. Raw-read rebuilding uses
#   Barcode_counting/count_Lib1_barcodes.sh followed by a manual count-table/sample-name
#   hand-off.
# Outputs:
# - preprocessing_scripts/out_25_AMG_guide_efficacy_inputs/ — b*_pao_join_raw_counts.tsv.gz and
#   standardized_designs/b*_experimental_design_standardized.tsv.
# Interpretation / operation:
# - No fitness test, T0 gate or score prediction is performed. Duplicate spacer/sample counts
#   are summed.
# - Missing counts for a batch-observed spacer become zero; guides absent throughout a batch
#   retain NA.
# - S8A–B reads this output layout. S8C–F also accepts frozen copies in edgeR_data/raw_counts/.

BATCHES <- c("b1", "b2", "b3", "b4")
COUNT_DIR <- file.path("Barcode_counting", "count_tables")
KEY_FILE <- file.path("annotation_input", "PAO1_key.tsv")
OUTPUT_DIR <- file.path("preprocessing_scripts", "out_25_AMG_guide_efficacy_inputs")

library(data.table)

# The key supplies the sequence context and existing off-target annotations.
key <- data.table::fread(KEY_FILE)
required <- c("spacer", "type", "locus_tag", "offset", "pam", "downstream", "off_target_n")
if (!all(required %in% names(key))) stop("Guide key lacks: ", paste(setdiff(required, names(key)), collapse = ", "))
if (anyNA(key$spacer) || anyDuplicated(key$spacer)) stop("The guide key needs one row per spacer.")
data.table::setnames(key, "locus_tag", "PAO1_ID")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
design_dir <- file.path(OUTPUT_DIR, "standardized_designs")
dir.create(design_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Repeat input preparation independently for each batch.
# Use the design's exact sample IDs and column order so downstream T0 lookups
# refer to the correct samples. No statistical comparison is performed here.
for (batch in BATCHES) {
  message("Preparing original guide/count input layout for ", batch)
  design <- data.table::fread(file.path(COUNT_DIR, paste0(batch, "_experimental_design.tsv")))
  if (!all(c("condition", "drug", "timing") %in% names(design))) stop("Incomplete design for ", batch)
  sample_columns <- as.character(design$condition)
  if (anyNA(sample_columns) || anyDuplicated(sample_columns)) stop("Invalid sample IDs in ", batch)

  # Raw barcode inputs have three columns: spacer, read count and sample ID.
  counts <- data.table::fread(file.path(COUNT_DIR, paste0(batch, "_overall_stats.tsv")),
                             header = FALSE, fill = TRUE)
  if (ncol(counts) < 3L) stop("Expected three count-table columns for ", batch)
  counts <- counts[, 1:3]
  data.table::setnames(counts, c("spacer", "count", "condition"))
  counts <- counts[!(spacer %in% c("barcode", "spacer") &
                      count %in% c("count", "counts") & condition %in% c("condition", "sample"))]
  counts[, count := as.numeric(count)]
  if (any(!is.finite(counts$count)) || any(counts$count < 0)) stop("Invalid raw counts in ", batch)
  if (length(setdiff(sample_columns, counts$condition))) stop("Design samples missing from counts for ", batch)
  # 2. Keep design-listed samples after checking that all are present.
  # Counts remain raw: no TMM normalization, T0 cutoff or fitness selection.
  counts <- counts[condition %in% sample_columns]

  # Sum duplicate count records, as in metacounter, then spread samples into
  # columns. A spacer absent from a sample has zero reads, not a missing value.
  counts <- counts[, .(count = sum(count)), by = .(spacer, condition)]
  wide <- data.table::dcast(counts, spacer ~ condition, value.var = "count", fill = 0)
  wide <- wide[, c("spacer", sample_columns), with = FALSE]

  # Preserve the canonical key's guide order and annotations. A guide never
  # observed anywhere in this batch stays NA across all sample columns, as in
  # the historical raw joins; a missing sample count for an observed guide is 0.
  # The model scripts apply their original abundance/profile filters themselves.
  wide <- wide[match(key$spacer, wide$spacer), ]
  # 3. Attach raw counts to each canonical guide's sequence/context annotations.
  # The key determines the guide universe; matching by spacer avoids relying
  # on row order in the count tables. No efficacy scores are fitted or updated.
  result <- cbind(data.table::copy(key), wide[, sample_columns, with = FALSE])
  data.table::fwrite(result, file.path(OUTPUT_DIR, paste0(batch, "_pao_join_raw_counts.tsv.gz")),
                     sep = "\t", na = "NA")

  # These fields are already sufficient for the original T0 lookup; copy the
  # provided design without changing sample names, groups, treatment or timing.
  # 4. Save the matching design beside the matrices so model scripts can select
  # untreated T0 columns. They apply their own abundance and training filters.
  data.table::fwrite(design, file.path(design_dir, paste0(batch, "_experimental_design_standardized.tsv")),
                     sep = "\t", na = "NA")
}

message("Prepared sequence/context/raw-count inputs: ", OUTPUT_DIR)
message("The statistical model inputs remain the existing edgeR_data/b*_melted_results.tsv.gz and b*_median_melted_results.tsv.gz files.")
