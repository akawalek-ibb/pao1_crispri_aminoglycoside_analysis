#!/usr/bin/env Rscript
############################################################
## PAO1 CRISPRi-seq -> INDEPENDENT AMINOGLYCOSIDE MAGeCK RRA
##
## Run from the project root: Source this complete file in RStudio.
## Inputs: Barcode_counting/count_tables/b*_overall_stats.tsv (spacer, count,
## sample; no header), b*_experimental_design.tsv and annotation_input/PAO1_key.tsv.
## sample_name_map.tsv maps the original and GEO-style sample IDs in both inputs.
## Source count/design files are not rewritten. The supplied inputs
## need no upstream R script. Raw FASTQs require barcode counting plus hand-off.
## Outputs: MAGeCK_data/mageck_rra_raw_AG/ contains four full-match count matrices,
## per-batch and combined manifests, QC/settings tables and Bash command files.
## This script prepares files only. Next run the generated master Bash script
## in your MAGeCK environment, then MAGeCK_data/02_MAGeCK_postprocessing.R.
##
## ANALYSIS SET:
##   B4 TOB 0.025 T1 vs matched None T1
##   B4 TOB 0.1   T1 vs matched None T1
##   B1 TOB 0.4   T1 vs matched None T1
##   B4 TOB 1     T1 vs matched None T1
##   B1 TOB 0.4   T2 vs matched None T2
##      -> longitudinal 14-generation follow-up; reported separately
##   B2 AMK 2     T1 vs matched None T1
##
## PRIMARY INPUT:
##   perfect/full-match targeting guides ONLY.
##   Mismatch guides are completely excluded from MAGeCK input.
##   Non-targeting controls are also excluded from MAGeCK input.
##
## MAGeCK:
##   raw annotated sgRNA counts
##     -> median normalization
##     -> MAGeCK NB sgRNA testing
##     -> alpha-RRA gene ranking
##     -> median gene-LFC
##
## No edgeR filterByExpr, no TMM, no p_stouffer, no mismatch guides.
## No mean-T0 abundance prefilter in the primary analysis.
##
## THREE CONTRAST FAMILIES:
##   1) drug_vs_matched_no_drug
##      Drug Tx vs untreated Tx
##      -> PRIMARY gene-drug interaction contrast
##
##   2) untreated_post_vs_T0
##      Untreated T1/T2 vs untreated T0
##      Manuscript control set: B1 T1, B1 T2, B2 T1, B3 T1, B4 T1
##      -> baseline CRISPRi fitness
##
##   3) drug_vs_T0
##      Drug T1/T2 vs untreated T0
##      -> total treated fitness relative to T0
##
## Expected totals:
##   7 direct manuscript treatment contrasts
##   5 untreated-post-vs-T0 contrasts
##   7 treatment-vs-T0 contrasts
##   19 contrasts total
##
## Input folder:
##
## Generated .sh files use relative paths and LF line endings for WSL.
## This R script prepares input tables and commands only. At completion it
## prints two terminal commands: activate mageckenv, then run the master .sh.
## RStudio does not need to find Conda or MAGeCK.
############################################################
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required. Install it with install.packages('data.table').")
}

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 1. PROJECT SETTINGS
# ============================================================
# Batch design/count inputs are stored in:
INPUT_DIR <- "Barcode_counting/count_tables"

# Guide annotation is stored in:
ANNOT_DIR <- file.path("annotation_input")

GUIDE_ANNOTATION_FILE <- file.path(ANNOT_DIR, "PAO1_key.tsv")

SAMPLE_NAME_MAP_FILE <- file.path("MAGeCK_data", "sample_name_map.tsv")

message("\nInput folder:\n  ", INPUT_DIR, "\nGuide annotation:\n  ", GUIDE_ANNOTATION_FILE)

OUTPUT_DIR <- "MAGeCK_data/mageck_rra_raw_AG"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. BATCH INPUTS
# ============================================================
BATCH_CONFIG <- list(b1 = list(design_file = file.path(INPUT_DIR, "b1_experimental_design.tsv"), 
  counts_file = file.path(INPUT_DIR, "b1_overall_stats.tsv")), b2 = list(design_file = file.path(INPUT_DIR, 
  "b2_experimental_design.tsv"), counts_file = file.path(INPUT_DIR, "b2_overall_stats.tsv")), 
  b3 = list(design_file = file.path(INPUT_DIR, "b3_experimental_design.tsv"), counts_file = file.path(INPUT_DIR, 
    "b3_overall_stats.tsv")), b4 = list(design_file = file.path(INPUT_DIR, "b4_experimental_design.tsv"), 
    counts_file = file.path(INPUT_DIR, "b4_overall_stats.tsv")))

BATCHES_TO_PROCESS <- c("b1", "b2", "b3", "b4")

# No batch is control-only in this analysis.
# B3 contains the FBS 60 T1 condition in addition to None T0/T1.
CONTROL_ONLY_BATCHES <- character(0)

CONTROL_ONLY_TIMINGS <- list()

# Exact manuscript treatment groups retained from the larger screen.
# These include the aminoglycoside conditions plus FBS.
AG_TARGETS <- data.table::data.table(batch = c("b4", "b4", "b1", "b4", "b1", "b2", "b3"), drug = c("TOB", 
  "TOB", "TOB", "TOB", "TOB", "AMK", "FBS"), dose = c("0.025", "0.1", "0.4", "1", "0.4", "2", 
  "60"), timing = c("T1", "T1", "T1", "T1", "T2", "T1", "T1"), analysis_role = c("recurrence_T1", 
  "recurrence_T1", "recurrence_T1", "recurrence_T1", "longitudinal_T2_14_generations", "recurrence_T1", 
  "FBS_nonrecurrence"))

EXPECTED_TOTAL_CONTRASTS <- 19L

EXPECTED_DRUG_CONTRASTS <- 7L

EXPECTED_BASELINE_CONTRASTS <- 5L

EXPECTED_TREATED_VS_T0_CONTRASTS <- 7L

# ============================================================
# 3. MAGeCK SETTINGS
# ============================================================
MAGECK_EXECUTABLE <- "mageck"

# Activate your MAGeCK environment in a terminal before running the generated
# shell script. Environment activation and MAGeCK execution are not done in R.

MAGECK_NORM_METHOD <- "median"

MAGECK_GENE_LFC_METHOD <- "median"

WRITE_MAGECK_NORMALIZED_COUNTS <- TRUE

GENERATE_DRUG_CONTRASTS <- TRUE

GENERATE_BASELINE_CRISPRI_CONTRASTS <- TRUE

GENERATE_TREATED_VS_T0_CONTRASTS <- TRUE

# ============================================================
# 4. OPTIONAL ROBUSTNESS FILTER
# ============================================================
#
# Recommended PRIMARY setting: FALSE
#
# When TRUE, this applies only your historical abundance gate.
# edgeR::filterByExpr is deliberately never used here.
APPLY_MEAN_T0_FILTER <- FALSE

MIN_MEAN_T0 <- 100

PRIMARY_GUIDE_TYPE <- "perfect"

# ============================================================
# 5. HELPERS
# ============================================================
normalize_dose <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.numeric(x_chr))
  numeric_idx <- is.finite(x_num) & !is.na(x_num)
  x_chr[numeric_idx] <- vapply(x_num[numeric_idx], function(z) {
    if (z == 0) {
      return("0")
    }
    format(z, scientific = FALSE, trim = TRUE, digits = 15)
  }, character(1))
  x_chr
}

normalize_guide_type <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[x %in% c("perfect", "full_match", "full-match", "fullmatch")] <- "perfect"
  out[x %in% c("mismatch", "mismatched", "mm")] <- "mismatch"
  out[x %in% c("control", "non-targeting", "nontargeting", "nt")] <- "control"
  unknown <- is.na(out)
  out[unknown] <- x[unknown]
  out
}

safe_filename <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(x == "", "unnamed", x)
}

quote_shell <- function(x) {
  paste0("\"", gsub("\"", "\\\"", as.character(x), fixed = TRUE), "\"")
}

write_lf_lines <- function(lines, path) {
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  payload <- paste0(paste(lines, collapse = "\n"), "\n")
  writeBin(charToRaw(payload), con)
  invisible(path)
}

# Match sample aliases to the deposited GEO labels without changing sample order.
standardize_sample_ids <- function(x, batch_id, strict = TRUE) {
  if (!file.exists(SAMPLE_NAME_MAP_FILE)) {
    stop("Missing sample-name mapping: ", SAMPLE_NAME_MAP_FILE)
  }
  sample_map <- data.table::fread(SAMPLE_NAME_MAP_FILE, sep = "\t", header = TRUE,
                                 colClasses = "character")
  required <- c("batch", "original_sample_id", "sample_id")
  if (!all(required %in% names(sample_map))) {
    stop("Sample-name mapping lacks: ", paste(setdiff(required, names(sample_map)), collapse = ", "))
  }
  sample_map <- sample_map[tolower(batch) == tolower(batch_id)]
  if (nrow(sample_map) == 0L || anyNA(sample_map) ||
      anyDuplicated(sample_map$original_sample_id) || anyDuplicated(sample_map$sample_id)) {
    stop("Missing or ambiguous sample-name mapping for ", batch_id)
  }
  aliases <- unique(data.table::rbindlist(list(
    sample_map[, .(input_id = original_sample_id, sample_id)],
    sample_map[, .(input_id = sub("^[bB][0-9]+_", "", original_sample_id), sample_id)],
    sample_map[, .(input_id = sample_id, sample_id)],
    sample_map[, .(input_id = sub("^[bB][0-9]+_", "", sample_id), sample_id)]
  )))
  aliases[, input_id := toupper(input_id)]
  aliases <- unique(aliases)
  if (anyDuplicated(aliases$input_id)) stop("Ambiguous sample aliases for ", batch_id)
  x <- trimws(as.character(x))
  idx <- match(toupper(x), aliases$input_id)
  if (isTRUE(strict) && anyNA(idx)) {
    stop("Samples absent from the GEO name map for ", batch_id, ": ",
         paste(unique(x[is.na(idx)]), collapse = ", "))
  }
  matched <- !is.na(idx)
  x[matched] <- aliases$sample_id[idx[matched]]
  x
}

standardize_design <- function(design_file, batch_id) {
  if (!file.exists(design_file)) {
    stop("Missing design file: ", design_file)
  }
  design <- data.table::fread(design_file, sep = "\t", header = TRUE, fill = TRUE, na.strings = c("", 
    "NA"), data.table = TRUE)
  required <- c("condition", "experiment", "drug", "dose", "timing", "rep")
  missing_required <- setdiff(required, names(design))
  if (length(missing_required) > 0L) {
    stop("Design ", batch_id, " is missing required columns: ", paste(missing_required, collapse = ", "))
  }
  design[, `:=`(condition, standardize_sample_ids(condition, batch_id, strict = FALSE))]
  design[, `:=`(experiment, trimws(as.character(experiment)))]
  design[, `:=`(drug, toupper(trimws(as.character(drug))))]
  design[, `:=`(timing, toupper(trimws(as.character(timing))))]
  design[, `:=`(rep, as.character(rep))]
  design[is.na(drug) | drug == "" | drug == "NONE", `:=`(drug, "None")]
  design[, `:=`(dose, normalize_dose(dose))]
  design[drug == "None", `:=`(dose, "0")]
  if (any(is.na(design$condition) | design$condition == "")) {
    stop("Empty condition name in ", batch_id)
  }
  if (anyDuplicated(design$condition)) {
    stop("Duplicated condition names in ", batch_id)
  }
  if (any(is.na(design$timing) | design$timing == "")) {
    stop("Missing timing value in ", batch_id)
  }
  if (any(design$drug != "None" & (is.na(design$dose) | design$dose == ""))) {
    stop("Missing treatment dose in ", batch_id)
  }
  if (any(is.na(design$rep) | design$rep == "")) {
    stop("Missing replicate label in ", batch_id)
  }
  design[, `:=`(group, paste(drug, dose, timing, sep = "_"))]
  design[, `:=`(batch, batch_id)]
  design[, `:=`(original_sample_order, seq_len(.N))]
  design[]
}

filter_design_to_aminoglycosides <- function(design, batch_id) {
  if (batch_id %in% CONTROL_ONLY_BATCHES) {
    timings_to_keep <- CONTROL_ONLY_TIMINGS[[batch_id]]
    if (is.null(timings_to_keep)) {
      stop("No CONTROL_ONLY_TIMINGS definition for ", batch_id)
    }
    filtered <- design[drug == "None" & timing %in% timings_to_keep]
    for (timing_i in timings_to_keep) {
      if (filtered[timing == timing_i, .N] == 0L) {
        stop(batch_id, ": required untreated ", timing_i, " samples are absent from the design.")
      }
    }
    message(batch_id, ": control-only batch; retaining untreated ", paste(timings_to_keep, 
      collapse = " and "), " samples only.")
    return(filtered[])
  }
  targets <- AG_TARGETS[batch == batch_id]
  if (nrow(targets) == 0L) {
    stop("No aminoglycoside target definition for ", batch_id, ". If this batch is intentionally control-only, add it to ", 
      "CONTROL_ONLY_BATCHES and CONTROL_ONLY_TIMINGS.")
  }
  required_timings <- unique(targets$timing)
  keep_treatment <- rep(FALSE, nrow(design))
  for (target_i in seq_len(nrow(targets))) {
    keep_treatment <- keep_treatment | (design$drug == targets$drug[[target_i]] & design$dose == 
      targets$dose[[target_i]] & design$timing == targets$timing[[target_i]])
  }
  keep_control <- design$drug == "None" & design$timing %in% c("T0", required_timings)
  filtered <- design[keep_treatment | keep_control]
  for (target_i in seq_len(nrow(targets))) {
    n_target <- filtered[drug == targets$drug[[target_i]] & dose == targets$dose[[target_i]] & 
      timing == targets$timing[[target_i]], .N]
    if (n_target == 0L) {
      stop(batch_id, ": requested aminoglycoside group is absent from design: ", targets$drug[[target_i]], 
        " ", targets$dose[[target_i]], " ", targets$timing[[target_i]])
    }
  }
  if (filtered[drug == "None" & timing == "T0", .N] == 0L) {
    stop(batch_id, ": no untreated T0 samples found after aminoglycoside filtering.")
  }
  for (timing_i in required_timings) {
    if (filtered[drug == "None" & timing == timing_i, .N] == 0L) {
      stop(batch_id, ": no matched untreated ", timing_i, " control found.")
    }
  }
  filtered[]
}

read_barcode_counts <- function(counts_file, batch_id) {
  if (!file.exists(counts_file)) {
    stop("Missing count file for ", batch_id, ": ", counts_file)
  }
  counts <- data.table::fread(counts_file, sep = "\t", header = FALSE, fill = TRUE, data.table = TRUE)
  if (ncol(counts) < 3L) {
    stop("Count file for ", batch_id, " must contain spacer, count, condition")
  }
  counts <- counts[, 1:3]
  data.table::setnames(counts, c("spacer", "count", "condition"))
  counts <- counts[!(tolower(as.character(spacer)) %in% c("barcode", "spacer") & tolower(as.character(count)) %in% 
    c("count", "counts") & tolower(as.character(condition)) %in% c("condition", "sample"))]
  counts[, `:=`(spacer, trimws(as.character(spacer)))]
  input_samples <- unique(trimws(as.character(counts$condition)))
  mapped_samples <- standardize_sample_ids(input_samples, batch_id, strict = FALSE)
  if (anyDuplicated(mapped_samples)) {
    stop("Multiple input labels map to the same sample in ", batch_id)
  }
  counts[, `:=`(condition, mapped_samples[match(trimws(as.character(condition)), input_samples)])]
  counts[, `:=`(count, suppressWarnings(as.numeric(count)))]
  if (anyNA(counts$count) || any(!is.finite(counts$count))) {
    stop("Non-numeric or non-finite counts in ", counts_file)
  }
  if (any(counts$count < 0)) {
    stop("Negative counts in ", counts_file)
  }
  if (any(abs(counts$count - round(counts$count)) > 1e-08)) {
    stop("Non-integer counts in ", counts_file)
  }
  counts[, `:=`(count, as.numeric(round(count)))]
  counts[, .(count = sum(count)), by = .(spacer, condition)]
}

read_guide_annotation <- function(path) {
  if (!file.exists(path)) {
    stop("Missing guide annotation: ", path)
  }
  key <- data.table::fread(path, sep = "\t", header = TRUE, fill = TRUE, data.table = TRUE)
  required <- c("spacer", "type", "locus_tag")
  missing_required <- setdiff(required, names(key))
  if (length(missing_required) > 0L) {
    stop("PAO1_key.tsv is missing required columns: ", paste(missing_required, collapse = ", "))
  }
  key[, `:=`(spacer, trimws(as.character(spacer)))]
  key[, `:=`(type, normalize_guide_type(type))]
  key[, `:=`(locus_tag, trimws(as.character(locus_tag)))]
  key <- key[!is.na(spacer) & spacer != ""]
  if (anyDuplicated(key$spacer)) {
    duplicate_count <- key[, .N, by = spacer][N > 1L, .N]
    warning(duplicate_count, " duplicated spacers in PAO1_key.tsv; first annotation retained")
    key <- unique(key, by = "spacer")
  }
  key[]
}

validate_design_vs_counts <- function(counts, design, batch_id) {
  missing_in_counts <- setdiff(design$condition, unique(counts$condition))
  extra_in_counts <- setdiff(unique(counts$condition), design$condition)
  if (length(missing_in_counts) > 0L) {
    stop(batch_id, ": design samples absent from counts: ", paste(missing_in_counts, collapse = ", "))
  }
  if (length(extra_in_counts) > 0L) {
    stop(batch_id, ": count-table samples absent from design: ", paste(extra_in_counts, collapse = ", "))
  }
  invisible(TRUE)
}

make_wide_raw_counts <- function(counts, design, guide_key, batch_id) {
  validate_design_vs_counts(counts, design, batch_id)
  unknown_spacers <- setdiff(unique(counts$spacer), guide_key$spacer)
  if (length(unknown_spacers) > 0L) {
    warning(batch_id, ": ", length(unknown_spacers), " count-table spacer(s) are absent from PAO1_key.tsv and will be excluded.")
  }
  counts_known <- counts[spacer %in% guide_key$spacer]
  scaffold <- data.table::CJ(spacer = unique(guide_key$spacer), condition = design$condition, 
    unique = TRUE)
  expanded <- merge(scaffold, counts_known, by = c("spacer", "condition"), all.x = TRUE, sort = FALSE)
  expanded[is.na(count), `:=`(count, 0)]
  wide <- data.table::dcast(expanded, spacer ~ condition, value.var = "count", fill = 0)
  desired_columns <- c("spacer", design$condition)
  wide <- wide[, ..desired_columns]
  wide[]
}

apply_optional_t0_filter <- function(wide, design, batch_id) {
  t0_samples <- design[drug == "None" & timing == "T0", condition]
  if (length(t0_samples) == 0L) {
    stop(batch_id, ": no untreated T0 samples found.")
  }
  mean_t0 <- rowMeans(as.matrix(wide[, ..t0_samples]), na.rm = TRUE)
  filter_status <- data.table::data.table(spacer = wide$spacer, mean_T0_count = mean_t0, passed_optional_T0_filter = mean_t0 >= 
    MIN_MEAN_T0)
  if (isTRUE(APPLY_MEAN_T0_FILTER)) {
    keep_spacers <- filter_status[passed_optional_T0_filter == TRUE, spacer]
    wide <- wide[spacer %in% keep_spacers]
  }
  list(wide = wide, filter_status = filter_status)
}

prepare_mageck_matrix <- function(wide, guide_key, design, guide_type, batch_id) {
  guide_map <- guide_key[type == guide_type]
  if (nrow(guide_map) == 0L) {
    return(data.table::data.table())
  }
  if (guide_type %in% c("perfect", "mismatch")) {
    bad_target <- guide_map[is.na(locus_tag) | locus_tag == ""]
    if (nrow(bad_target) > 0L) {
      stop(batch_id, ": ", nrow(bad_target), " ", guide_type, " guide(s) have no locus_tag in PAO1_key.tsv.")
    }
  }
  sample_columns <- design$condition
  subset_wide <- wide[spacer %in% guide_map$spacer]
  subset_wide <- merge(guide_map[, .(spacer, type, locus_tag)], subset_wide, by = "spacer", all.y = TRUE, 
    sort = FALSE)
  subset_wide[, `:=`(guide_order__, match(spacer, guide_map$spacer))]
  data.table::setorder(subset_wide, guide_order__)
  subset_wide[, `:=`(guide_order__, NULL)]
  if (guide_type == "control") {
    subset_wide[, `:=`(Gene, paste0("NTC_", seq_len(.N)))]
  } else {
    subset_wide[, `:=`(Gene, locus_tag)]
  }
  out <- subset_wide[, c(list(sgRNA = spacer, Gene = Gene), mget(sample_columns))]
  out[]
}

build_contrast_manifest <- function(design, batch_id) {
  contrast_rows <- list()
  row_i <- 0L
  treatment_groups <- unique(design[drug != "None", .(drug, dose, timing, group)])
  t0_samples <- design[drug == "None" & timing == "T0", condition]
  if ((isTRUE(GENERATE_BASELINE_CRISPRI_CONTRASTS) || isTRUE(GENERATE_TREATED_VS_T0_CONTRASTS)) && 
    length(t0_samples) == 0L) {
    stop(batch_id, ": untreated T0 samples are required for baseline and/or drug-vs-T0 contrasts.")
  }
  if (isTRUE(GENERATE_DRUG_CONTRASTS)) {
    if (nrow(treatment_groups) > 0L) {
      for (group_i in seq_len(nrow(treatment_groups))) {
        treatment_row <- treatment_groups[group_i]
        treatment_samples <- design[group == treatment_row$group, condition]
        control_samples <- design[drug == "None" & timing == treatment_row$timing, condition]
        if (length(control_samples) == 0L) {
         stop(batch_id, ": no matched no-drug control at timing ", treatment_row$timing, 
          " for ", treatment_row$group, ". Direct drug-vs-matched-control analysis cannot proceed.")
        }
        row_i <- row_i + 1L
        contrast_rows[[row_i]] <- data.table::data.table(batch = batch_id, contrast_type = "drug_vs_matched_no_drug", 
         contrast = paste0(treatment_row$drug, "_", treatment_row$dose, "_", treatment_row$timing, 
          "_vs_None_0_", treatment_row$timing), drug = treatment_row$drug, dose = treatment_row$dose, 
         timing = treatment_row$timing, treatment_samples = paste(treatment_samples, collapse = ","), 
         control_samples = paste(control_samples, collapse = ","))
      }
    }
  }
  if (isTRUE(GENERATE_BASELINE_CRISPRI_CONTRASTS)) {
    if (length(t0_samples) == 0L) {
      warning(batch_id, ": no untreated T0 samples; baseline CRISPRi contrasts skipped.")
    } else {
      post_t0_timings <- unique(design[drug == "None" & timing != "T0", timing])
      for (timing_i in post_t0_timings) {
        treatment_samples <- design[drug == "None" & timing == timing_i, condition]
        if (length(treatment_samples) == 0L) {
         next
        }
        row_i <- row_i + 1L
        contrast_rows[[row_i]] <- data.table::data.table(batch = batch_id, contrast_type = "untreated_post_vs_T0", 
         contrast = paste0("None_0_", timing_i, "_vs_None_0_T0"), drug = "None", dose = "0", 
         timing = timing_i, treatment_samples = paste(treatment_samples, collapse = ","), 
         control_samples = paste(t0_samples, collapse = ","))
      }
    }
  }
  if (isTRUE(GENERATE_TREATED_VS_T0_CONTRASTS)) {
    if (length(t0_samples) == 0L) {
      warning(batch_id, ": no untreated T0 samples; drug-vs-T0 contrasts skipped.")
    } else if (nrow(treatment_groups) > 0L) {
      for (group_i in seq_len(nrow(treatment_groups))) {
        treatment_row <- treatment_groups[group_i]
        treatment_samples <- design[group == treatment_row$group, condition]
        if (length(treatment_samples) == 0L) {
         next
        }
        row_i <- row_i + 1L
        contrast_rows[[row_i]] <- data.table::data.table(batch = batch_id, contrast_type = "drug_vs_T0", 
         contrast = paste0(treatment_row$drug, "_", treatment_row$dose, "_", treatment_row$timing, 
          "_vs_None_0_T0"), drug = treatment_row$drug, dose = treatment_row$dose, timing = treatment_row$timing, 
         treatment_samples = paste(treatment_samples, collapse = ","), control_samples = paste(t0_samples, 
          collapse = ","))
      }
    }
  }
  if (length(contrast_rows) == 0L) {
    return(data.table::data.table())
  }
  contrasts <- data.table::rbindlist(contrast_rows, fill = TRUE)
  duplicate_contrasts <- contrasts[, .N, by = .(batch, contrast)][N > 1L]
  if (nrow(duplicate_contrasts) > 0L) {
    stop(batch_id, ": duplicated MAGeCK contrast names were generated:\n", paste(duplicate_contrasts$contrast, 
      collapse = "\n"))
  }
  bad_sample_sets <- contrasts[is.na(treatment_samples) | treatment_samples == "" | is.na(control_samples) | 
    control_samples == ""]
  if (nrow(bad_sample_sets) > 0L) {
    stop(batch_id, ": one or more MAGeCK contrasts have an empty treatment/control sample set.")
  }
  contrasts[]
}

# Quote each argument for the Bash runners, including executable paths with
# spaces. MAGeCK uses the environment activated in the calling terminal.
# Count files, sample names and all MAGeCK analysis options are retained.
mageck_shell_command <- function(arguments) {
  parts <- c(MAGECK_EXECUTABLE, arguments)
  paste(vapply(parts, shQuote, character(1), type = "sh"), collapse = " ")
}

build_mageck_command <- function(count_table, output_prefix, treatment_samples, control_samples) {
  args <- c("test", "-k", count_table, "-t", treatment_samples,
    "-c", control_samples, "--norm-method", MAGECK_NORM_METHOD, "--gene-lfc-method",
    MAGECK_GENE_LFC_METHOD, "-n", output_prefix)
  if (isTRUE(WRITE_MAGECK_NORMALIZED_COUNTS)) {
    args <- c(args, "--normcounts-to-file")
  }
  mageck_shell_command(args)
}

# ============================================================
# 6. LOAD GUIDE ANNOTATION
# ============================================================
message("\nGuide annotation file:\n  ", GUIDE_ANNOTATION_FILE)

guide_key <- read_guide_annotation(GUIDE_ANNOTATION_FILE)

guide_type_summary <- guide_key[, .(n_guides = .N), by = type][order(type)]

message("\nGuide annotation summary:")

print(guide_type_summary)

unknown_types <- setdiff(unique(guide_key$type), c("perfect", "mismatch", "control"))

if (length(unknown_types) > 0L) {
  warning("Unrecognized guide type(s) in PAO1_key.tsv: ", paste(unknown_types, collapse = ", "), 
    "\nThey will be excluded from the primary MAGeCK analysis.")
}

# ============================================================
# 7. PROCESS BATCHES
# ============================================================
all_contrasts <- list()

all_batch_qc <- list()

all_commands <- character()

for (batch_id in BATCHES_TO_PROCESS) {
  message("\n============================================================", "\nProcessing ", 
    batch_id, "\n============================================================")
  config <- BATCH_CONFIG[[batch_id]]
  if (is.null(config)) {
    stop("No BATCH_CONFIG entry for ", batch_id)
  }
  batch_out <- file.path(OUTPUT_DIR, batch_id)
  dir.create(batch_out, recursive = TRUE, showWarnings = FALSE)
  # Read the biological assignments and map count/design sample aliases to
  # the same IDs. This changes the in-memory labels, not the source files.
  design <- standardize_design(config$design_file, batch_id)
  design <- filter_design_to_aminoglycosides(design, batch_id)
  design[, condition := standardize_sample_ids(condition, batch_id)]
  # Read spacer counts and keep only samples retained in this batch design.
  # These are raw counts; normalization and statistical testing happen in MAGeCK.
  counts <- read_barcode_counts(config$counts_file, batch_id)
  n_count_rows_before_AG_filter <- nrow(counts)
  counts <- counts[condition %in% design$condition]
  if (nrow(counts) == 0L) {
    stop(batch_id, ": no count rows remain after aminoglycoside-only sample filtering.")
  }
  message(batch_id, ": retained ", nrow(design), " sample columns for aminoglycoside analysis; count rows ", 
    n_count_rows_before_AG_filter, " -> ", nrow(counts), ".")
  # Record observed spacers absent from the canonical key so unassigned reads
  # remain auditable rather than being interpreted as gene-targeting guides.
  unknown_spacers_i <- sort(setdiff(unique(counts$spacer), guide_key$spacer))
  data.table::fwrite(data.table::data.table(spacer = unknown_spacers_i), file.path(batch_out, 
    paste0(batch_id, "_spacers_absent_from_PAO1_key.tsv")), sep = "\t", quote = FALSE)
  message(batch_id, ": spacers absent from PAO1_key.tsv = ", length(unknown_spacers_i), ". Audit file: ", 
    file.path(batch_out, paste0(batch_id, "_spacers_absent_from_PAO1_key.tsv")))
  # Arrange counts as guides by samples, then apply the optional T0 gate only
  # if enabled above. Select full-match targeting guides for the primary matrix.
  wide_raw <- make_wide_raw_counts(counts = counts, design = design, guide_key = guide_key, batch_id = batch_id)
  filter_result <- apply_optional_t0_filter(wide = wide_raw, design = design, batch_id = batch_id)
  wide_for_mageck <- filter_result$wide
  primary_matrix <- prepare_mageck_matrix(wide = wide_for_mageck, guide_key = guide_key, design = design, 
    guide_type = PRIMARY_GUIDE_TYPE, batch_id = batch_id)
  if (nrow(primary_matrix) == 0L) {
    stop(batch_id, ": no perfect/full-match guides available for MAGeCK.")
  }
  # Save the exact matrix and design used to build commands. These copies
  # are also inputs to postprocessing QC, so retain them with the gene summaries.
  primary_file <- file.path(batch_out, paste0(batch_id, "_MAGeCK_fullmatch_raw_counts.tsv"))
  data.table::fwrite(primary_matrix, primary_file, sep = "\t", quote = FALSE)
  data.table::fwrite(filter_result$filter_status, file.path(batch_out, paste0(batch_id, "_optional_T0_filter_status.tsv")), 
    sep = "\t")
  data.table::fwrite(design, file.path(batch_out, paste0(batch_id, "_design_used.tsv")), sep = "\t")
  # Define treatment/control sample groups for each comparison family.
  # Drug/control measures interaction; drug/T0 also contains baseline repression effects.
  contrasts <- build_contrast_manifest(design, batch_id)
  if (nrow(contrasts) == 0L) {
    warning(batch_id, ": no MAGeCK contrasts generated.")
    next
  }
  result_dir <- file.path(batch_out, "results")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  contrasts[, `:=`(count_table, primary_file)]
  contrasts[, `:=`(output_prefix, file.path(result_dir, safe_filename(paste0(batch, "__", contrast))))]
  contrasts[, `:=`(count_table_shell, basename(primary_file))]
  contrasts[, `:=`(output_prefix_shell, file.path("results", safe_filename(paste0(batch, "__", 
    contrast))))]
  # Generate shell commands without executing them. The batch runner changes
  # into its own folder before using the relative matrix and result paths.
  contrasts[, `:=`(command, mapply(FUN = build_mageck_command, count_table = count_table_shell, 
    output_prefix = output_prefix_shell, treatment_samples = treatment_samples, control_samples = control_samples, 
    USE.NAMES = FALSE))]
  data.table::fwrite(contrasts, file.path(batch_out, paste0(batch_id, "_MAGeCK_contrast_manifest.tsv")), 
    sep = "\t", quote = FALSE)
  command_file <- file.path(batch_out, paste0(batch_id, "_run_MAGeCK_RRA.sh"))
  write_lf_lines(c("#!/usr/bin/env bash", "set -euo pipefail", "", "SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"", 
    "cd \"$SCRIPT_DIR\"", "", contrasts$command), command_file)
  all_commands <- c(all_commands, paste0("bash \"$SCRIPT_DIR/", batch_id, "/", batch_id, "_run_MAGeCK_RRA.sh\""))
  # Collect batch records for the combined manifest and check how many
  # comparisons and guides were prepared. These are input QC, not hit calls.
  all_contrasts[[batch_id]] <- contrasts
  qc_row <- data.table::data.table(batch = batch_id, n_samples = nrow(design), n_raw_count_rows = nrow(counts), 
    n_annotated_spacers = nrow(guide_key), n_perfect_guides_annotation = guide_key[type == 
      "perfect", .N], n_mismatch_guides_annotation = guide_key[type == "mismatch", .N], n_control_guides_annotation = guide_key[type == 
      "control", .N], n_count_table_spacers_absent_from_PAO1_key = length(unknown_spacers_i), 
    n_guides_before_optional_T0_filter = nrow(wide_raw), n_guides_after_optional_T0_filter = nrow(wide_for_mageck), 
    n_primary_fullmatch_guides = nrow(primary_matrix), primary_analysis_uses_fullmatch_only = TRUE, 
    mismatches_in_primary_analysis = FALSE, controls_in_primary_analysis = FALSE, optional_T0_filter_applied = APPLY_MEAN_T0_FILTER, 
    min_mean_T0_if_applied = MIN_MEAN_T0, mageck_norm_method = MAGECK_NORM_METHOD, mageck_gene_lfc_method = MAGECK_GENE_LFC_METHOD, 
    n_contrasts = nrow(contrasts), n_drug_vs_matched_no_drug = contrasts[contrast_type == "drug_vs_matched_no_drug", 
      .N], n_untreated_post_vs_T0 = contrasts[contrast_type == "untreated_post_vs_T0", .N], 
    n_drug_vs_T0 = contrasts[contrast_type == "drug_vs_T0", .N])
  all_batch_qc[[batch_id]] <- qc_row
  message(batch_id, ":", "\n  full-match guides in primary MAGeCK matrix: ", nrow(primary_matrix), 
    "\n  mismatch guides in MAGeCK input: 0", "\n  non-targeting controls in MAGeCK input: 0", 
    "\n  count-table spacers absent from PAO1_key.tsv: ", length(unknown_spacers_i), "\n  total contrasts generated: ", 
    nrow(contrasts), "\n    drug vs matched no-drug: ", contrasts[contrast_type == "drug_vs_matched_no_drug", 
      .N], "\n    untreated post-T0 vs T0: ", contrasts[contrast_type == "untreated_post_vs_T0", 
      .N], "\n    drug vs T0: ", contrasts[contrast_type == "drug_vs_T0", .N])
}

# ============================================================
# 8. CROSS-BATCH OUTPUTS
# ============================================================
if (length(all_contrasts) > 0L) {
  contrast_manifest_all <- data.table::rbindlist(all_contrasts, fill = TRUE)
  data.table::fwrite(contrast_manifest_all, file.path(OUTPUT_DIR, "all_batches_MAGeCK_contrast_manifest.tsv"), 
    sep = "\t", quote = FALSE)
  contrast_family_summary <- contrast_manifest_all[, .(n_contrasts = .N), by = .(batch, contrast_type)]
  data.table::setorder(contrast_family_summary, batch, contrast_type)
  data.table::fwrite(contrast_family_summary, file.path(OUTPUT_DIR, "all_batches_MAGeCK_contrast_family_summary.tsv"), 
    sep = "\t", quote = FALSE)
  n_total <- nrow(contrast_manifest_all)
  n_direct <- contrast_manifest_all[contrast_type == "drug_vs_matched_no_drug", .N]
  n_baseline <- contrast_manifest_all[contrast_type == "untreated_post_vs_T0", .N]
  n_treated_t0 <- contrast_manifest_all[contrast_type == "drug_vs_T0", .N]
  if (n_total != EXPECTED_TOTAL_CONTRASTS || n_direct != EXPECTED_DRUG_CONTRASTS || n_baseline != 
    EXPECTED_BASELINE_CONTRASTS || n_treated_t0 != EXPECTED_TREATED_VS_T0_CONTRASTS) {
    stop("Unexpected aminoglycoside MAGeCK contrast counts.\n", "Expected total/direct/baseline/treated-vs-T0 = ", 
      EXPECTED_TOTAL_CONTRASTS, "/", EXPECTED_DRUG_CONTRASTS, "/", EXPECTED_BASELINE_CONTRASTS, 
      "/", EXPECTED_TREATED_VS_T0_CONTRASTS, "; observed = ", n_total, "/", n_direct, "/", 
      n_baseline, "/", n_treated_t0, ".")
  }
  expected_direct_labels <- c("TOB_0.025_T1_vs_None_0_T1", "TOB_0.1_T1_vs_None_0_T1", "TOB_0.4_T1_vs_None_0_T1", 
    "TOB_1_T1_vs_None_0_T1", "TOB_0.4_T2_vs_None_0_T2", "AMK_2_T1_vs_None_0_T1", "FBS_60_T1_vs_None_0_T1")
  observed_direct_labels <- contrast_manifest_all[contrast_type == "drug_vs_matched_no_drug", 
    contrast]
  missing_direct <- setdiff(expected_direct_labels, observed_direct_labels)
  extra_direct <- setdiff(observed_direct_labels, expected_direct_labels)
  if (length(missing_direct) > 0L || length(extra_direct) > 0L) {
    stop("Direct manuscript contrast manifest does not match the expected seven conditions.\n", 
      "Missing: ", paste(missing_direct, collapse = ", "), "\nExtra: ", paste(extra_direct, 
        collapse = ", "))
  }
}

if (length(all_batch_qc) > 0L) {
  batch_qc_all <- data.table::rbindlist(all_batch_qc, fill = TRUE)
  data.table::fwrite(batch_qc_all, file.path(OUTPUT_DIR, "all_batches_MAGeCK_input_QC_summary.tsv"), 
    sep = "\t", quote = FALSE)
}

if (length(all_commands) > 0L) {
  master_command_file <- file.path(OUTPUT_DIR, "run_all_MAGeCK_RRA.sh")
  # A manual master-script run records the executable version as well.
  version_command <- paste0(mageck_shell_command("--version"),
                            ' 2>&1 | tee "$SCRIPT_DIR/MAGeCK_version.txt"')
  write_lf_lines(c("#!/usr/bin/env bash", "set -euo pipefail", "", "SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"",
    "", version_command, "", all_commands), master_command_file)
}

data.table::fwrite(guide_type_summary, file.path(OUTPUT_DIR, "guide_type_summary_from_PAO1_key.tsv"), 
  sep = "\t")

settings_table <- data.table::data.table(setting = c("MAGECK_EXECUTABLE", "MAGECK_NORM_METHOD", 
  "MAGECK_GENE_LFC_METHOD", "WRITE_MAGECK_NORMALIZED_COUNTS", "GENERATE_DRUG_CONTRASTS", "GENERATE_BASELINE_CRISPRI_CONTRASTS", 
  "GENERATE_TREATED_VS_T0_CONTRASTS", "APPLY_MEAN_T0_FILTER", "MIN_MEAN_T0", 
  "PRIMARY_GUIDE_TYPE"), value = c(MAGECK_EXECUTABLE, MAGECK_NORM_METHOD, MAGECK_GENE_LFC_METHOD, 
  as.character(WRITE_MAGECK_NORMALIZED_COUNTS), as.character(GENERATE_DRUG_CONTRASTS), as.character(GENERATE_BASELINE_CRISPRI_CONTRASTS), 
  as.character(GENERATE_TREATED_VS_T0_CONTRASTS), as.character(APPLY_MEAN_T0_FILTER), 
  as.character(MIN_MEAN_T0), PRIMARY_GUIDE_TYPE))

data.table::fwrite(settings_table, file.path(OUTPUT_DIR, "MAGeCK_analysis_settings.tsv"), sep = "\t")

message("\n============================================================", "\nMAGeCK RAW-COUNT INPUT PREPARATION COMPLETE", 
  "\n============================================================", "\nOutput directory:", "\n  ", 
  OUTPUT_DIR, "\n", "\nPRIMARY RRA input:", "\n  full-match/perfect guides only", "\n", "\nMismatch guides:", 
  "\n  EXCLUDED completely from MAGeCK input", "\n", "\nNon-targeting controls:", "\n  EXCLUDED from MAGeCK input", 
  "\n", "\nStatistical prefilter before MAGeCK:", "\n  ", if (isTRUE(APPLY_MEAN_T0_FILTER)) {
    paste0("mean untreated T0 >= ", MIN_MEAN_T0, " (SECONDARY robustness mode)")
  } else {
    "NONE — raw annotated guide counts are passed to MAGeCK"
  }, "\n", "\nMAGeCK normalization:", "\n  ", MAGECK_NORM_METHOD, "\n", "\nMAGeCK gene-LFC method:", 
  "\n  ", MAGECK_GENE_LFC_METHOD, "\n", "\nContrast families:", "\n  1) drug_vs_matched_no_drug = primary gene-drug interaction", 
  "\n  2) untreated_post_vs_T0 = baseline CRISPRi fitness", "\n  3) drug_vs_T0 = total treated fitness relative to T0", 
  "\n", "\nGenerated master command file:", "\n  ", file.path(OUTPUT_DIR, "run_all_MAGeCK_RRA.sh"), 
  "\n", "\nNEXT STEPS — in a terminal opened at the project root:",
  "\nconda activate mageckenv",
  "\nbash ", paste0(file.path(OUTPUT_DIR, "run_all_MAGeCK_RRA.sh")),
  "\n\nChange mageckenv above if your environment has a different name.",
  "\nIf the correct environment is already active, run only the bash line.",
  "\nAfter MAGeCK finishes, run MAGeCK_data/02_MAGeCK_postprocessing.R in R.",
  "\nVersion record when executed: ", file.path(OUTPUT_DIR, "MAGeCK_version.txt"), "\n============================================================")

