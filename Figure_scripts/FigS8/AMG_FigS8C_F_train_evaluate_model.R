#!/usr/bin/env Rscript
# Train a Lasso model of relative guide-associated fitness and evaluate held-out genes.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/FigS8/AMG_FigS8C_F_train_evaluate_model.R")
# Inputs:
# - edgeR_data/b*_melted_results.tsv.gz and b*_median_melted_results.tsv.gz.
# - edgeR_data/raw_counts/b*_pao_join_raw_counts.tsv.gz (preferred); script 25 outputs are
#   fallback inputs.
# - Barcode_counting/count_tables/b*_experimental_design.tsv and annotation_input/PAO1_key.tsv.
# - Optional externally supplied catalogue:
#   Figure_scripts/FigS8/input_AMG_FigS8/PAO1_107_relaxed_unfiletered.csv.
# Upstream / execution order:
# - Use the supplied fitted screen and frozen count/context files; no upstream rerun is required
#   for those inputs.
# - preprocessing_scripts/02_AMG_metacounter.R rebuilds fitted exports;
#   preprocessing_scripts/25_AMG_guide_efficacy_inputs.R rebuilds count/context inputs.
# Outputs:
# - Figure_scripts/FigS8/out_AMG_FigS8C_F_model/ — model RDS, training/scores, coefficients,
#   CV/validation data, an XLSX results workbook and PNG/TIFF/PDF.
# Interpretation / operation:
# - Training requires raw gene p < 0.05, medLFC ≤ −1, mean T0 ≥100, ≥2 batches/profiles and ≥5
#   eligible guides per gene.
# - The response is within-gene-centered fitness efficacy, not measured transcript knockdown.
# - 100 outer gene-held-out splits and gene-grouped inner CV are used. Feature-set selection
#   reuses outer results, so it is not an independent final test.
# - D–F illustrate the split closest to median outer Spearman performance. New model files do
#   not replace the preserved scoring RDS.
#
# R packages used/checked: data.table, dplyr, ggplot2, ggseqlogo, glmnet, openxlsx, purrr,
# readr, readxl, scales, stringr, tibble, tidyr.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tibble)
  library(purrr)
  library(glmnet)
  library(ggplot2)
  library(ggseqlogo)
  library(scales)
  library(openxlsx)
})

# ============================================================
# USER SETTINGS
# ============================================================

BASE_SEED <- 42L
set.seed(BASE_SEED)

PROJECT_DIR <- "."

EDGE_DIR <- "edgeR_data"

# Frozen annotated raw-count files from 25_AMG_guide_efficacy_inputs.R.
# This is the preferred location for b1-b4_pao_join_raw_counts.tsv.gz.
RAW_COUNTS_DIR <- file.path(EDGE_DIR, "raw_counts")

# Original preprocessing output remains available as a fallback.
# No edgeR fit is rerun when reading either location.
PREPARED_INPUT_DIR <- file.path("preprocessing_scripts", "out_25_AMG_guide_efficacy_inputs")

OUT_DIR <- file.path("Figure_scripts", "FigS8", "out_AMG_FigS8C_F_model")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Gene/guide filtering
P_STOUFFER_CUT <- 0.05
NEGATIVE_LFC_CUT <- -1
MIN_MEAN_T0 <- 100
MIN_BATCHES_PER_GUIDE <- 2L
MIN_PROFILES_PER_GUIDE <- 2L
MIN_GUIDES_PER_GENE <- 5L

# Sequence context
SPACER_LENGTH <- 20L
DOWNSTREAM_LENGTH <- 15L

# Model
ALPHA_LASSO <- 1
INNER_NFOLDS <- 10L
OUTER_VALIDATION_FRACTION <- 0.20
N_OUTER_REPEATS <- 100L
LAMBDA_CHOICE <- "lambda.min"  # lambda.min or lambda.1se

# Automatic deployed-model selection:
# choose the simplest model whose mean outer-validation MSE
# is within one standard error of the model with minimum MSE.
DEPLOYED_MODEL_SELECTION <- "auto_one_se"
FORCE_DEPLOYED_MODEL <- NA_character_
# Valid forced values:
#   "sequence_only"
#   "sequence_offset"
#   "sequence_offset_offtarget"

# Plot/export
BASE_FAMILY <- "Arial"
PLOT_DPI <- 300

CONTROL_MANIFEST <- tibble::tribble(
  ~batch, ~condition,                    ~profile_label, ~profile_order,
  "b1",   "None_0_T1 - None_0_T0",      "B1 T1",       1L,
  "b1",   "None_0_T2 - None_0_T0",      "B1 T2",       2L,
  "b2",   "None_0_T1 - None_0_T0",      "B2 T1",       3L,
  "b3",   "None_0_T1 - None_0_T0",      "B3 T1",       4L,
  "b4",   "None_0_T1 - None_0_T0",      "B4 T1",       5L
)

# Canonical table containing every guide used in Lib1.
# Use the frozen B1 annotated raw-count file when present; retain the
# original preprocessing/metacounter locations as fallbacks. The selected
# file must already contain off_target_n; no separate guide tables are used.
LIB1_GUIDE_CANDIDATES <- c(
  file.path(RAW_COUNTS_DIR, "b1_pao_join_raw_counts.tsv.gz"),
  file.path(PREPARED_INPUT_DIR, "b1_pao_join_raw_counts.tsv.gz"),
  file.path(EDGE_DIR, "b1", "b1_pao_join_raw_counts.tsv.gz")
)

# Require complete off-target annotation for every perfect Lib1 guide.
# Mismatch and control guides are allowed to have NA because their actual
# altered/non-targeting sequences require separate off-target calculations.
STRICT_LIB1_OFFTARGET_COVERAGE <- TRUE

# Complete PAO1 guide-design catalogue(s) identified in the
# original machine-learning/design scripts. Every existing file
# is scored and exported separately.
PAO1_DESIGN_CANDIDATES <- c(file.path("Figure_scripts", "FigS8", "input_AMG_FigS8", "PAO1_107_relaxed_unfiletered.csv"))


OFFTARGET_SOURCE_CANDIDATES <- c(file.path("annotation_input", "PAO1_key.tsv"))

# ============================================================
# GENERAL HELPERS
# ============================================================

first_existing <- function(paths, description, required = TRUE) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) {
    if (required) {
      stop(
        "Could not locate ", description, ". Checked:\n",
        paste(paths, collapse = "\n")
      )
    }
    return(NA_character_)
  }
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

find_batch_file <- function(batch, suffix) {
  candidates <- c(
    if (identical(suffix, "pao_join_raw_counts.tsv.gz"))
      file.path(RAW_COUNTS_DIR, paste0(batch, "_", suffix)),
    file.path(PREPARED_INPUT_DIR, paste0(batch, "_", suffix)),
    file.path(EDGE_DIR, batch, paste0(batch, "_", suffix)),
    file.path(EDGE_DIR, batch, suffix),
    file.path(EDGE_DIR, paste0(batch, "_", suffix))
  )
  first_existing(candidates, paste0(batch, " ", suffix))
}

find_standardized_design <- function(batch) {
  candidates <- c(
    file.path(PREPARED_INPUT_DIR, "standardized_designs",
              paste0(batch, "_experimental_design_standardized.tsv")),
    file.path(
      EDGE_DIR,
      "standardized_designs",
      paste0(batch, "_experimental_design_standardized.tsv")
    ),
    file.path(
      EDGE_DIR,
      "standardized_designs",
      paste0(batch, "_experimental_design.tsv")
    ),
    file.path("Barcode_counting", "count_tables", paste0(batch, "_experimental_design.tsv"))
  )
  first_existing(candidates, paste0(batch, " experimental design"))
}

standardize_gene_column <- function(df) {
  if (!"PAO1_ID" %in% names(df)) {
    candidate <- intersect(
      c("locus_tag", "Locus Tag", "gene_id"),
      names(df)
    )[1]
    if (is.na(candidate)) stop("No PAO1_ID/locus_tag column found.")
    data.table::setnames(df, candidate, "PAO1_ID")
  }
  df[, PAO1_ID := as.character(PAO1_ID)]
  df
}

standardize_guide_columns <- function(df) {
  df <- standardize_gene_column(df)
  
  rename_first <- function(target, candidates) {
    if (target %in% names(df)) return(invisible(NULL))
    hit <- intersect(candidates, names(df))[1]
    if (!is.na(hit)) data.table::setnames(df, hit, target)
  }
  
  rename_first("spacer", c("sgrna", "guide", "target"))
  rename_first("pam", c("PAM"))
  rename_first("downstream", c("downstream_sequence", "downstream.seq"))
  rename_first("offset", c("distance", "target_offset"))
  rename_first(
    "off_target_n",
    c("offtarget", "off_target", "offtargets", "off_target_count")
  )
  rename_first("type", c("guide_type"))
  rename_first("original", c("parent", "parent_spacer", "full_sgrna"))
  
  df
}


load_offtarget_map <- local({
  cached_map <- NULL
  
  function() {
    if (!is.null(cached_map)) return(data.table::copy(cached_map))
    
    source_file <- first_existing(
      OFFTARGET_SOURCE_CANDIDATES,
      "PAO1 genome-wide guide-design table containing off_target_n"
    )
    
    source_dt <- data.table::fread(source_file, header = TRUE)
    
    spacer_col <- intersect(
      c("sgrna", "spacer", "guide", "target"),
      names(source_dt)
    )[1]
    offtarget_col <- intersect(
      c("off_target_n", "offtarget", "off_target", "offtargets", "off_target_count"),
      names(source_dt)
    )[1]
    
    if (is.na(spacer_col) || is.na(offtarget_col)) {
      stop(
        "Could not identify guide/off-target columns in ", source_file,
        ". Available columns: ", paste(names(source_dt), collapse = ", ")
      )
    }
    
    map <- source_dt[
      ,
      .(
        spacer = toupper(gsub("\\s+", "", as.character(get(spacer_col)))),
        off_target_n_lookup = suppressWarnings(as.numeric(get(offtarget_col)))
      )
    ][
      !is.na(spacer) & spacer != "" & is.finite(off_target_n_lookup)
    ]
    
    conflicts <- map[
      ,
      .(n_values = uniqueN(off_target_n_lookup)),
      by = spacer
    ][n_values > 1L]
    
    if (nrow(conflicts) > 0L) {
      stop(
        "Conflicting off-target counts were found for ",
        nrow(conflicts), " spacer sequence(s) in ", source_file
      )
    }
    
    cached_map <<- map[
      ,
      .(off_target_n_lookup = off_target_n_lookup[1]),
      by = spacer
    ]
    
    message("Loaded off-target annotation from: ", source_file)
    data.table::copy(cached_map)
  }
})

attach_offtarget_annotation <- function(df) {
  if (!"spacer" %in% names(df)) {
    stop("Cannot add off_target_n because the input has no spacer column.")
  }
  
  df[, row_order_offtarget__ := .I]
  df[, spacer := toupper(gsub("\\s+", "", as.character(spacer)))]
  
  if ("off_target_n" %in% names(df)) {
    df[, off_target_n := suppressWarnings(as.numeric(off_target_n))]
  }
  
  df <- merge(
    df,
    load_offtarget_map(),
    by = "spacer",
    all.x = TRUE,
    sort = FALSE
  )
  
  if (!"off_target_n" %in% names(df)) {
    df[, off_target_n := NA_real_]
  }
  
  conflicts <- df[
    is.finite(off_target_n) &
      is.finite(off_target_n_lookup) &
      off_target_n != off_target_n_lookup
  ]
  
  if (nrow(conflicts) > 0L) {
    warning(
      nrow(conflicts),
      " row(s) contain an off_target_n value that disagrees with the design table. ",
      "The existing input value is retained."
    )
  }
  
  df[
    !is.finite(off_target_n),
    off_target_n := off_target_n_lookup
  ]
  
  data.table::setorder(df, row_order_offtarget__)
  df[, c("row_order_offtarget__", "off_target_n_lookup") := NULL]
  
  message(
    "off_target_n coverage: ",
    sum(is.finite(df$off_target_n)), " / ", nrow(df),
    " rows (", round(100 * mean(is.finite(df$off_target_n)), 2), "%)"
  )
  
  df
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  if (dplyr::n_distinct(x[ok]) < 2L) return(NA_real_)
  if (dplyr::n_distinct(y[ok]) < 2L) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

add_sheet <- function(wb, sheet_name, x) {
  openxlsx::addWorksheet(wb, sheet_name)
  
  x_df <- as.data.frame(x)
  if (ncol(x_df) == 0L) {
    x_df <- data.frame(note = "No rows were generated.")
  }
  
  openxlsx::writeDataTable(
    wb,
    sheet = sheet_name,
    x = x_df,
    tableStyle = "TableStyleMedium2"
  )
  openxlsx::freezePane(wb, sheet_name, firstRow = TRUE)
}

export_plot <- function(plot_object, stem, width, height) {
  print(plot_object)
  
  ggplot2::ggsave(
    file.path(OUT_DIR, paste0(stem, ".png")),
    plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = PLOT_DPI
  )
  
  ggplot2::ggsave(
    file.path(OUT_DIR, paste0(stem, ".tiff")),
    plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = PLOT_DPI,
    compression = "lzw"
  )
  
  ggplot2::ggsave(
    file.path(OUT_DIR, paste0(stem, ".pdf")),
    plot_object,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )
}

sanitize_stem <- function(x) {
  x %>%
    basename() %>%
    tools::file_path_sans_ext() %>%
    stringr::str_replace_all("[^A-Za-z0-9._-]+", "_")
}

# ============================================================
# SEQUENCE FEATURE HELPERS
# ============================================================

valid_dna <- function(x, expected_length) {
  x <- toupper(as.character(x))
  !is.na(x) &
    nchar(x) == expected_length &
    stringr::str_detect(x, paste0("^[ACGT]{", expected_length, "}$"))
}

one_hot_acg <- function(seqs, prefix) {
  seqs <- toupper(as.character(seqs))
  lengths <- unique(nchar(seqs))
  if (length(lengths) != 1L) {
    stop("Sequences supplied to one_hot_acg have unequal lengths.")
  }
  
  npos <- lengths[1]
  X <- matrix(
    0,
    nrow = length(seqs),
    ncol = npos * 3L
  )
  
  colnames(X) <- paste0(
    prefix,
    rep(seq_len(npos), each = 3L),
    "_",
    rep(c("A", "C", "G"), times = npos)
  )
  
  for (i in seq_len(npos)) {
    b <- substr(seqs, i, i)
    j <- (i - 1L) * 3L
    X[, j + 1L] <- as.numeric(b == "A")
    X[, j + 2L] <- as.numeric(b == "C")
    X[, j + 3L] <- as.numeric(b == "G")
  }
  
  as.data.frame(X, check.names = FALSE)
}

build_full_design_matrix <- function(df) {
  required <- c(
    "spacer", "pam", "downstream_model", "offset", "off_target_n"
  )
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing feature columns: ", paste(missing, collapse = ", "))
  }
  
  spacer_X <- one_hot_acg(df$spacer, "spacer_")
  down_X <- one_hot_acg(df$downstream_model, "down_")
  
  pam_n <- substr(df$pam, 1, 1)
  pam_X <- tibble::tibble(
    pam_1_A = as.numeric(pam_n == "A"),
    pam_1_C = as.numeric(pam_n == "C"),
    pam_1_G = as.numeric(pam_n == "G")
  )
  
  meta_X <- tibble::tibble(
    offset_log1p = log1p(pmax(as.numeric(df$offset), 0)),
    offtarget_log1p = log1p(pmax(as.numeric(df$off_target_n), 0))
  )
  
  dplyr::bind_cols(spacer_X, pam_X, down_X, meta_X) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))
}

MODEL_FEATURE_SETS <- list(
  sequence_only = function(X) {
    X[, !names(X) %in% c("offset_log1p", "offtarget_log1p"), drop = FALSE]
  },
  sequence_offset = function(X) {
    X[, names(X) != "offtarget_log1p", drop = FALSE]
  },
  sequence_offset_offtarget = function(X) X
)

MODEL_COMPLEXITY_ORDER <- c(
  "sequence_only",
  "sequence_offset",
  "sequence_offset_offtarget"
)

# ============================================================
# GROUPED CROSS-VALIDATION HELPERS
# ============================================================

make_group_folds <- function(groups, nfolds, seed) {
  groups <- as.character(groups)
  unique_groups <- unique(groups)
  if (length(unique_groups) < 3L) {
    stop("At least three genes are required for grouped CV.")
  }
  
  nfolds_use <- min(as.integer(nfolds), length(unique_groups))
  set.seed(seed)
  shuffled <- sample(unique_groups)
  group_fold <- rep(seq_len(nfolds_use), length.out = length(shuffled))
  names(group_fold) <- shuffled
  as.integer(group_fold[groups])
}

split_genes <- function(groups, validation_fraction, seed) {
  genes <- unique(as.character(groups))
  if (length(genes) < 10L) {
    stop("Too few genes for a held-out validation split.")
  }
  
  set.seed(seed)
  n_val <- max(1L, floor(length(genes) * validation_fraction))
  val_genes <- sample(genes, n_val)
  
  list(
    train_idx = which(!groups %in% val_genes),
    validation_idx = which(groups %in% val_genes),
    validation_genes = val_genes
  )
}

fit_grouped_cv_lasso <- function(
    X,
    y,
    groups,
    seed,
    lambda_choice = "lambda.min") {
  
  foldid <- make_group_folds(groups, INNER_NFOLDS, seed)
  
  fit <- glmnet::cv.glmnet(
    x = as.matrix(X),
    y = as.numeric(y),
    family = "gaussian",
    alpha = ALPHA_LASSO,
    standardize = TRUE,
    intercept = TRUE,
    foldid = foldid,
    keep = FALSE
  )
  
  lambda <- switch(
    lambda_choice,
    "lambda.min" = fit$lambda.min,
    "lambda.1se" = fit$lambda.1se,
    stop("Unknown LAMBDA_CHOICE: ", lambda_choice)
  )
  
  list(
    cvfit = fit,
    lambda = lambda,
    foldid = foldid
  )
}

predict_lasso <- function(fit_object, X) {
  as.numeric(
    predict(
      fit_object$cvfit,
      newx = as.matrix(X),
      s = fit_object$lambda,
      type = "link"
    )
  )
}

evaluate_predictions <- function(df, pred, model_name, repeat_id) {
  eval_df <- df %>%
    dplyr::mutate(predicted = as.numeric(pred))
  
  gene_cor <- eval_df %>%
    dplyr::group_by(PAO1_ID) %>%
    dplyr::filter(dplyr::n() >= 3L) %>%
    dplyr::summarise(
      n_guides = dplyr::n(),
      spearman_r = safe_cor(predicted, y, "spearman"),
      pearson_r = safe_cor(predicted, y, "pearson"),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(spearman_r))
  
  tibble::tibble(
    repeat_id = repeat_id,
    model = model_name,
    n_guides = nrow(eval_df),
    n_genes = dplyr::n_distinct(eval_df$PAO1_ID),
    pearson = safe_cor(eval_df$predicted, eval_df$y, "pearson"),
    spearman = safe_cor(eval_df$predicted, eval_df$y, "spearman"),
    mse = mean((eval_df$predicted - eval_df$y)^2, na.rm = TRUE),
    median_gene_spearman = stats::median(gene_cor$spearman_r, na.rm = TRUE),
    n_gene_correlations = nrow(gene_cor)
  )
}

# ============================================================
# LOAD AND AGGREGATE B1-B4 CONTROL DATA
# ============================================================

profile_list <- vector("list", nrow(CONTROL_MANIFEST))
gene_call_list <- vector("list", nrow(CONTROL_MANIFEST))

for (i in seq_len(nrow(CONTROL_MANIFEST))) {
  batch_id <- CONTROL_MANIFEST$batch[i]
  condition_id <- CONTROL_MANIFEST$condition[i]
  profile_label <- CONTROL_MANIFEST$profile_label[i]
  
  message("Loading model profile ", profile_label, ": ", condition_id)
  
  melted_file <- find_batch_file(batch_id, "melted_results.tsv.gz")
  median_file <- find_batch_file(batch_id, "median_melted_results.tsv.gz")
  raw_file <- find_batch_file(batch_id, "pao_join_raw_counts.tsv.gz")
  design_file <- find_standardized_design(batch_id)
  
  melted <- standardize_gene_column(data.table::fread(melted_file))
  median_tbl <- standardize_gene_column(data.table::fread(median_file))
  raw_join <- standardize_guide_columns(data.table::fread(raw_file))
  raw_join <- attach_offtarget_annotation(raw_join)
  design <- data.table::fread(design_file)
  
  required_melted <- c(
    "PAO1_ID", "spacer", "type", "condition", "LFC.adj"
  )
  missing_melted <- setdiff(required_melted, names(melted))
  if (length(missing_melted) > 0) {
    stop(
      "Missing melted-result columns for ", batch_id, ": ",
      paste(missing_melted, collapse = ", ")
    )
  }
  
  required_features <- c(
    "PAO1_ID", "spacer", "pam", "downstream", "offset", "off_target_n"
  )
  missing_features <- setdiff(required_features, names(raw_join))
  if (length(missing_features) > 0) {
    stop(
      "Missing guide-feature columns for ", batch_id, ": ",
      paste(missing_features, collapse = ", ")
    )
  }
  
  if (!all(c("condition", "drug", "timing") %in% names(design))) {
    stop(
      "Design for ", batch_id,
      " must contain condition, drug and timing."
    )
  }
  
  t0_samples <- design[
    toupper(as.character(drug)) == "NONE" &
      toupper(as.character(timing)) == "T0",
    as.character(condition)
  ]
  
  missing_t0 <- setdiff(t0_samples, names(raw_join))
  if (length(missing_t0) > 0) {
    stop(
      "Missing T0 columns for ", batch_id, ":\n",
      paste(missing_t0, collapse = "\n")
    )
  }
  
  raw_join[
    ,
    mean_T0_count := rowMeans(as.matrix(.SD), na.rm = TRUE),
    .SDcols = t0_samples
  ]
  
  feature_tbl <- unique(
    raw_join[
      ,
      .(
        spacer = as.character(spacer),
        PAO1_ID_feature = as.character(PAO1_ID),
        pam = toupper(as.character(pam)),
        downstream = toupper(as.character(downstream)),
        offset_feature = suppressWarnings(as.numeric(offset)),
        off_target_n = suppressWarnings(as.numeric(off_target_n)),
        mean_T0_count = suppressWarnings(as.numeric(mean_T0_count))
      )
    ]
  )
  
  profile_tbl <- melted[
    condition == condition_id &
      type == "perfect" &
      !is.na(PAO1_ID) &
      PAO1_ID != "",
    .(
      PAO1_ID = as.character(PAO1_ID),
      gene_name = as.character(gene_name),
      spacer = toupper(as.character(spacer)),
      offset_melted = suppressWarnings(as.numeric(offset)),
      LFC.adj = suppressWarnings(as.numeric(LFC.adj))
    )
  ]
  
  profile_tbl <- merge(
    profile_tbl,
    feature_tbl,
    by = "spacer",
    all.x = TRUE
  )
  
  profile_tbl[
    is.na(PAO1_ID) | PAO1_ID == "",
    PAO1_ID := PAO1_ID_feature
  ]
  
  profile_tbl[
    ,
    offset := data.table::fcoalesce(offset_melted, offset_feature)
  ]
  
  profile_tbl[
    ,
    `:=`(
      batch = batch_id,
      profile_label = profile_label,
      profile_order = CONTROL_MANIFEST$profile_order[i],
      efficacy = -LFC.adj
    )
  ]
  
  profile_tbl <- profile_tbl[
    is.finite(LFC.adj) &
      is.finite(mean_T0_count) &
      mean_T0_count >= MIN_MEAN_T0 &
      is.finite(offset) &
      is.finite(off_target_n)
  ]
  
  profile_tbl[
    ,
    y_profile := efficacy - stats::median(efficacy, na.rm = TRUE),
    by = .(profile_label, PAO1_ID)
  ]
  
  profile_list[[i]] <- profile_tbl
  
  gene_calls <- median_tbl[
    condition == condition_id &
      !is.na(PAO1_ID) &
      PAO1_ID != "",
    .(
      PAO1_ID = as.character(PAO1_ID),
      gene_name = as.character(gene_name),
      batch = batch_id,
      profile_label = profile_label,
      profile_order = CONTROL_MANIFEST$profile_order[i],
      medLFC = suppressWarnings(as.numeric(medLFC)),
      p_stouffer = suppressWarnings(as.numeric(p_stouffer))
    )
  ]
  
  gene_calls[
    ,
    significant_negative :=
      is.finite(p_stouffer) &
      p_stouffer < P_STOUFFER_CUT &
      is.finite(medLFC) &
      medLFC <= NEGATIVE_LFC_CUT
  ]
  
  gene_call_list[[i]] <- gene_calls
  
  rm(melted, median_tbl, raw_join, design, feature_tbl)
  invisible(gc())
}

profile_long <- data.table::rbindlist(profile_list, use.names = TRUE, fill = TRUE)
gene_calls <- data.table::rbindlist(gene_call_list, use.names = TRUE, fill = TRUE)

# Select genes once using the original untreated negative-fitness calls; all eligible guide profiles are then used.
training_gene_ids <- gene_calls[
  significant_negative == TRUE,
  unique(PAO1_ID)
]

if (length(training_gene_ids) == 0) {
  stop("No screen-defined negative-fitness genes were identified.")
}

message("Candidate training genes: ", length(training_gene_ids))

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & as.character(x) != ""]
  if (length(x) == 0) return(NA_character_)
  as.character(x[1])
}

# Hierarchical aggregation prevents B1 from receiving double weight
# because it has both T1 and T2. Profile-centered outcomes are
# summarized within batch first, then across independent batches.
# B1 has two timepoints. Collapse within batch before taking the across-batch median.
Lib_batch <- profile_long[
  PAO1_ID %in% training_gene_ids,
  .(
    gene_name = first_nonmissing(gene_name),
    pam = first_nonmissing(pam),
    downstream = first_nonmissing(downstream),
    offset = stats::median(offset, na.rm = TRUE),
    off_target_n = stats::median(off_target_n, na.rm = TRUE),
    y_batch = stats::median(y_profile, na.rm = TRUE),
    LFC_batch = stats::median(LFC.adj, na.rm = TRUE),
    efficacy_batch = stats::median(efficacy, na.rm = TRUE),
    n_profiles_batch = data.table::uniqueN(profile_label),
    min_mean_T0_count_batch = min(mean_T0_count, na.rm = TRUE),
    median_mean_T0_count_batch = stats::median(mean_T0_count, na.rm = TRUE)
  ),
  by = .(batch, PAO1_ID, spacer)
]

Lib <- Lib_batch[
  ,
  .(
    gene_name = first_nonmissing(gene_name),
    pam = first_nonmissing(pam),
    downstream = first_nonmissing(downstream),
    offset = stats::median(offset, na.rm = TRUE),
    off_target_n = stats::median(off_target_n, na.rm = TRUE),
    y = stats::median(y_batch, na.rm = TRUE),
    LFC_agg = stats::median(LFC_batch, na.rm = TRUE),
    efficacy_agg = stats::median(efficacy_batch, na.rm = TRUE),
    n_batches = data.table::uniqueN(batch),
    n_profiles = sum(n_profiles_batch, na.rm = TRUE),
    min_mean_T0_count = min(min_mean_T0_count_batch, na.rm = TRUE),
    median_mean_T0_count = stats::median(
      median_mean_T0_count_batch,
      na.rm = TRUE
    )
  ),
  by = .(PAO1_ID, spacer)
]

Lib[
  ,
  `:=`(
    spacer = toupper(as.character(spacer)),
    pam = toupper(as.character(pam)),
    downstream = toupper(as.character(downstream)),
    downstream_model = substr(
      toupper(as.character(downstream)),
      1,
      DOWNSTREAM_LENGTH
    )
  )
]

Lib[
  ,
  valid_sequence_context :=
    valid_dna(spacer, SPACER_LENGTH) &
    !is.na(pam) &
    nchar(pam) == 3L &
    stringr::str_detect(pam, "^[ACGT]GG$") &
    valid_dna(downstream_model, DOWNSTREAM_LENGTH)
]

Lib <- Lib[
  n_batches >= MIN_BATCHES_PER_GUIDE &
    n_profiles >= MIN_PROFILES_PER_GUIDE &
    valid_sequence_context == TRUE &
    is.finite(y) &
    is.finite(offset) &
    is.finite(off_target_n)
]

Lib[
  ,
  n_guides_gene := .N,
  by = PAO1_ID
]

# The final per-gene minimum is applied after guide abundance, reproducibility and context checks.
Lib <- Lib[n_guides_gene >= MIN_GUIDES_PER_GENE]

if (nrow(Lib) == 0) stop("No guides remained after model filtering.")
if (data.table::uniqueN(Lib$PAO1_ID) < 20L) {
  stop("Too few genes remained for robust gene-held-out validation.")
}

message("Training guides: ", nrow(Lib))
message("Training genes: ", data.table::uniqueN(Lib$PAO1_ID))
message("Median profiles per guide: ", stats::median(Lib$n_profiles))

# ============================================================
# BUILD DESIGN MATRICES
# ============================================================

X_full <- build_full_design_matrix(Lib)
X_models <- lapply(MODEL_FEATURE_SETS, function(fun) fun(X_full))

# ============================================================
# REPEATED OUTER GENE-HELD-OUT VALIDATION / ABLATION
# ============================================================

outer_metric_list <- vector(
  "list",
  N_OUTER_REPEATS * length(X_models)
)
metric_index <- 1L

# Keep every guide of a gene on the same side of each train/validation split.
for (repeat_id in seq_len(N_OUTER_REPEATS)) {
  if (repeat_id %% 10L == 0L || repeat_id == 1L) {
    message("Outer validation repeat ", repeat_id, " / ", N_OUTER_REPEATS)
  }
  
  split <- split_genes(
    groups = Lib$PAO1_ID,
    validation_fraction = OUTER_VALIDATION_FRACTION,
    seed = BASE_SEED + repeat_id
  )
  
  for (model_name in names(X_models)) {
    X_current <- X_models[[model_name]]
    
    fit <- fit_grouped_cv_lasso(
      X = X_current[split$train_idx, , drop = FALSE],
      y = Lib$y[split$train_idx],
      groups = Lib$PAO1_ID[split$train_idx],
      seed = BASE_SEED + 10000L + repeat_id,
      lambda_choice = LAMBDA_CHOICE
    )
    
    pred <- predict_lasso(
      fit,
      X_current[split$validation_idx, , drop = FALSE]
    )
    
    outer_metric_list[[metric_index]] <- evaluate_predictions(
      df = as.data.frame(Lib[split$validation_idx]),
      pred = pred,
      model_name = model_name,
      repeat_id = repeat_id
    )
    
    metric_index <- metric_index + 1L
  }
}

outer_metrics <- dplyr::bind_rows(outer_metric_list)

model_summary <- outer_metrics %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    repeats = dplyr::n(),
    mean_mse = mean(mse, na.rm = TRUE),
    sd_mse = stats::sd(mse, na.rm = TRUE),
    se_mse = sd_mse / sqrt(repeats),
    mean_spearman = mean(spearman, na.rm = TRUE),
    sd_spearman = stats::sd(spearman, na.rm = TRUE),
    mean_pearson = mean(pearson, na.rm = TRUE),
    mean_median_gene_spearman = mean(
      median_gene_spearman,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    model = factor(model, levels = MODEL_COMPLEXITY_ORDER)
  ) %>%
  dplyr::arrange(model)

# ============================================================
# CHOOSE DEPLOYED FEATURE SET
# ============================================================

# Feature-set selection uses the repeated outer results below. Report this as model selection, not an untouched final test.
if (!is.na(FORCE_DEPLOYED_MODEL)) {
  if (!FORCE_DEPLOYED_MODEL %in% names(X_models)) {
    stop("Invalid FORCE_DEPLOYED_MODEL: ", FORCE_DEPLOYED_MODEL)
  }
  deployed_model <- FORCE_DEPLOYED_MODEL
} else if (DEPLOYED_MODEL_SELECTION == "auto_one_se") {
  best_row <- model_summary %>%
    dplyr::filter(mean_mse == min(mean_mse, na.rm = TRUE)) %>%
    dplyr::slice(1)
  
  mse_limit <- best_row$mean_mse + best_row$se_mse
  
  deployed_model <- model_summary %>%
    dplyr::filter(mean_mse <= mse_limit) %>%
    dplyr::arrange(model) %>%
    dplyr::slice(1) %>%
    dplyr::pull(model) %>%
    as.character()
} else {
  stop("Unknown DEPLOYED_MODEL_SELECTION: ", DEPLOYED_MODEL_SELECTION)
}

message("Deployed feature set: ", deployed_model)

X_deployed <- X_models[[deployed_model]]

# ============================================================
# REPRESENTATIVE OUTER SPLIT FOR VALIDATION FIGURES
# ============================================================

deployed_repeat_metrics <- outer_metrics %>%
  dplyr::filter(model == deployed_model)

median_outer_spearman <- stats::median(
  deployed_repeat_metrics$spearman,
  na.rm = TRUE
)

# This plotted split is chosen near the median validation correlation, rather than the best split.
representative_repeat <- deployed_repeat_metrics %>%
  dplyr::mutate(distance_to_median = abs(spearman - median_outer_spearman)) %>%
  dplyr::arrange(distance_to_median, repeat_id) %>%
  dplyr::slice(1) %>%
  dplyr::pull(repeat_id)

representative_split <- split_genes(
  groups = Lib$PAO1_ID,
  validation_fraction = OUTER_VALIDATION_FRACTION,
  seed = BASE_SEED + representative_repeat
)

representative_fit <- fit_grouped_cv_lasso(
  X = X_deployed[representative_split$train_idx, , drop = FALSE],
  y = Lib$y[representative_split$train_idx],
  groups = Lib$PAO1_ID[representative_split$train_idx],
  seed = BASE_SEED + 10000L + representative_repeat,
  lambda_choice = LAMBDA_CHOICE
)

representative_pred <- predict_lasso(
  representative_fit,
  X_deployed[representative_split$validation_idx, , drop = FALSE]
)

validation_df <- as.data.frame(Lib[representative_split$validation_idx]) %>%
  dplyr::mutate(
    predicted_efficacy = representative_pred,
    data_split = "validation"
  )

validation_metrics <- evaluate_predictions(
  validation_df,
  representative_pred,
  deployed_model,
  representative_repeat
)

validation_gene_cor <- validation_df %>%
  dplyr::group_by(PAO1_ID, gene_name) %>%
  dplyr::filter(dplyr::n() >= 3L) %>%
  dplyr::summarise(
    n_guides = dplyr::n(),
    spearman_r = safe_cor(predicted_efficacy, y, "spearman"),
    pearson_r = safe_cor(predicted_efficacy, y, "pearson"),
    .groups = "drop"
  ) %>%
  dplyr::filter(is.finite(spearman_r))

# ============================================================
# FIT FINAL MODEL ON ALL TRAINING DATA
# ============================================================

# Fit the deployment model only after evaluation; its in-sample scores are not held-out performance.
final_fit <- fit_grouped_cv_lasso(
  X = X_deployed,
  y = Lib$y,
  groups = Lib$PAO1_ID,
  seed = BASE_SEED + 50000L,
  lambda_choice = LAMBDA_CHOICE
)

final_predictions <- predict_lasso(final_fit, X_deployed)
Lib[, predicted_efficacy_score := final_predictions]
Lib[
  ,
  predicted_rank_within_gene := data.table::frank(
    -predicted_efficacy_score,
    ties.method = "min"
  ),
  by = PAO1_ID
]

# ============================================================
# FIGURE: SEQUENCE COEFFICIENT LOGO
# ============================================================

coef_matrix <- as.matrix(
  coef(final_fit$cvfit, s = final_fit$lambda)
)

coefficient_tbl <- tibble::tibble(
  feature = rownames(coef_matrix),
  coefficient = as.numeric(coef_matrix[, 1])
) %>%
  dplyr::mutate(abs_coefficient = abs(coefficient)) %>%
  dplyr::arrange(dplyr::desc(abs_coefficient))

beta <- coefficient_tbl$coefficient
names(beta) <- coefficient_tbl$feature

# Convert A/C/G dummy coefficients with T as the reference into
# centered four-base effects. This is the correct sum-to-zero
# display: subtract the mean of A, C, G and T(=0) from all four.
reference_to_centered_logo <- function(beta_vec, prefix, npos) {
  M <- matrix(
    0,
    nrow = 4,
    ncol = npos,
    dimnames = list(c("A", "C", "G", "T"), NULL)
  )
  
  for (i in seq_len(npos)) {
    raw_effects <- c(
      A = ifelse(
        paste0(prefix, i, "_A") %in% names(beta_vec),
        beta_vec[paste0(prefix, i, "_A")],
        0
      ),
      C = ifelse(
        paste0(prefix, i, "_C") %in% names(beta_vec),
        beta_vec[paste0(prefix, i, "_C")],
        0
      ),
      G = ifelse(
        paste0(prefix, i, "_G") %in% names(beta_vec),
        beta_vec[paste0(prefix, i, "_G")],
        0
      ),
      T = 0
    )
    
    M[, i] <- raw_effects - mean(raw_effects)
  }
  
  M
}

spacer_M <- reference_to_centered_logo(
  beta,
  "spacer_",
  SPACER_LENGTH
)

downstream_M <- reference_to_centered_logo(
  beta,
  "down_",
  DOWNSTREAM_LENGTH
)

pam_raw <- c(
  A = ifelse("pam_1_A" %in% names(beta), beta["pam_1_A"], 0),
  C = ifelse("pam_1_C" %in% names(beta), beta["pam_1_C"], 0),
  G = ifelse("pam_1_G" %in% names(beta), beta["pam_1_G"], 0),
  T = 0
)

pam_M <- matrix(
  0,
  nrow = 4,
  ncol = 3,
  dimnames = list(c("A", "C", "G", "T"), NULL)
)
pam_M[, 1] <- pam_raw - mean(pam_raw)

logo_M <- cbind(spacer_M, pam_M, downstream_M)

x_labels <- c(
  as.character(seq_len(SPACER_LENGTH)),
  "N", "G", "G",
  paste0("+", seq_len(DOWNSTREAM_LENGTH))
)

pam_start <- SPACER_LENGTH + 1L
pam_end <- SPACER_LENGTH + 3L

p_logo <- ggseqlogo::ggseqlogo(
  logo_M,
  method = "custom"
) +
  ggplot2::annotate(
    "rect",
    xmin = pam_start - 0.5,
    xmax = pam_end + 0.5,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.12
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq_len(ncol(logo_M)),
    labels = x_labels
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Lasso coefficient"
  ) +
  ggplot2::theme_minimal(
    base_size = 11,
    base_family = BASE_FAMILY
  ) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 7),
    axis.title.y = ggplot2::element_text(face = "bold")
  )

export_plot(p_logo, "sequence_feature_coefficient_logo", 11.0, 4.2)

# ============================================================
# FIGURE: FINAL GROUPED CV CURVE
# ============================================================

cv_curve_df <- tibble::tibble(
  lambda = final_fit$cvfit$lambda,
  log_lambda = log(final_fit$cvfit$lambda),
  mse = final_fit$cvfit$cvm,
  sd = final_fit$cvfit$cvsd
)

p_cv <- ggplot2::ggplot(
  cv_curve_df,
  ggplot2::aes(x = log_lambda, y = mse)
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mse - sd, ymax = mse + sd),
    width = 0.03,
    alpha = 0.45
  ) +
  ggplot2::geom_point(size = 1.2) +
  ggplot2::geom_line(linewidth = 0.45) +
  ggplot2::geom_vline(
    xintercept = log(final_fit$lambda),
    linetype = "dashed",
    linewidth = 0.7
  ) +
  ggplot2::labs(
    x = expression(log(lambda)),
    y = "Grouped cross-validation mean squared error"
  ) +
  ggplot2::theme_minimal(
    base_size = 12,
    base_family = BASE_FAMILY
  ) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(face = "bold")
  )

export_plot(p_cv, "lasso_grouped_cross_validation_curve", 6.5, 4.8)

# ============================================================
# FIGURE: HELD-OUT VALIDATION SCATTER
# ============================================================

val_pearson <- validation_metrics$pearson[1]
val_spearman <- validation_metrics$spearman[1]
val_mse <- validation_metrics$mse[1]

p_validation <- ggplot2::ggplot(
  validation_df,
  ggplot2::aes(x = predicted_efficacy, y = y)
) +
  ggplot2::geom_point(alpha = 0.25, size = 0.8) +
  ggplot2::geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 0.8
  ) +
  ggplot2::annotate(
    "label",
    x = -Inf,
    y = Inf,
    hjust = -0.05,
    vjust = 1.05,
    size = 3.1,
    label = paste0(
      "Validation genes = ", dplyr::n_distinct(validation_df$PAO1_ID),
      "\nGuides = ", nrow(validation_df),
      "\nPearson r = ", sprintf("%.3f", val_pearson),
      "\nSpearman rho = ", sprintf("%.3f", val_spearman),
      "\nMSE = ", signif(val_mse, 3)
    )
  ) +
  ggplot2::labs(
    x = "Predicted within-gene centered efficacy",
    y = "Observed within-gene centered efficacy"
  ) +
  ggplot2::theme_minimal(
    base_size = 12,
    base_family = BASE_FAMILY
  ) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(face = "bold")
  )

export_plot(p_validation, "heldout_gene_validation_scatter", 6.4, 5.2)

# ============================================================
# FIGURE: PER-GENE SPEARMAN DISTRIBUTION
# ============================================================

median_validation_gene_rho <- stats::median(
  validation_gene_cor$spearman_r,
  na.rm = TRUE
)

p_gene_cor <- ggplot2::ggplot(
  validation_gene_cor,
  ggplot2::aes(x = spearman_r)
) +
  ggplot2::geom_histogram(
    bins = 30,
    fill = "grey75",
    color = "black"
  ) +
  ggplot2::geom_vline(
    xintercept = median_validation_gene_rho,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  ggplot2::annotate(
    "label",
    x = Inf,
    y = Inf,
    hjust = 1.05,
    vjust = 1.1,
    size = 3.1,
    label = paste0(
      "Genes = ", nrow(validation_gene_cor),
      "\nMedian rho = ",
      sprintf("%.3f", median_validation_gene_rho)
    )
  ) +
  ggplot2::labs(
    x = "Per-gene Spearman correlation",
    y = "Number of validation genes"
  ) +
  ggplot2::theme_minimal(
    base_size = 12,
    base_family = BASE_FAMILY
  ) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(face = "bold")
  )

export_plot(p_gene_cor, "heldout_gene_spearman_distribution", 6.4, 4.8)

# ============================================================
# FIGURE: TOP-GUIDE ENRICHMENT
# ============================================================

# The expected observed LFC of a randomly selected guide is the
# within-gene mean. This removes Monte Carlo noise while matching
# the expectation of repeated random guide selection.
top_guide_enrichment <- validation_df %>%
  dplyr::group_by(PAO1_ID, gene_name) %>%
  dplyr::filter(dplyr::n() >= 3L) %>%
  dplyr::summarise(
    n_guides = dplyr::n(),
    selected_spacer = spacer[which.max(predicted_efficacy)],
    selected_predicted_score = max(predicted_efficacy, na.rm = TRUE),
    selected_observed_LFC = LFC_agg[which.max(predicted_efficacy)],
    expected_random_observed_LFC = mean(LFC_agg, na.rm = TRUE),
    delta_LFC_selected_minus_random =
      selected_observed_LFC - expected_random_observed_LFC,
    .groups = "drop"
  )

median_top_delta <- stats::median(
  top_guide_enrichment$delta_LFC_selected_minus_random,
  na.rm = TRUE
)

wilcox_top <- tryCatch(
  stats::wilcox.test(
    top_guide_enrichment$delta_LFC_selected_minus_random,
    mu = 0,
    alternative = "less",
    exact = FALSE
  ),
  error = function(e) NULL
)

p_top <- ggplot2::ggplot(
  top_guide_enrichment,
  ggplot2::aes(x = delta_LFC_selected_minus_random)
) +
  ggplot2::geom_histogram(
    bins = 35,
    fill = "grey75",
    color = "black"
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    color = "grey40",
    linewidth = 0.45
  ) +
  ggplot2::geom_vline(
    xintercept = median_top_delta,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  ggplot2::annotate(
    "label",
    x = Inf,
    y = Inf,
    hjust = 1.05,
    vjust = 1.1,
    size = 3.1,
    label = paste0(
      "Genes = ", nrow(top_guide_enrichment),
      "\nMedian delta LFC = ", sprintf("%.3f", median_top_delta),
      "\nOne-sided Wilcoxon p = ",
      ifelse(
        is.null(wilcox_top),
        "NA",
        format.pval(wilcox_top$p.value, digits = 3, eps = 1e-300)
      )
    )
  ) +
  ggplot2::labs(
    x = "Observed LFC: model-selected guide minus random-guide expectation",
    y = "Number of validation genes"
  ) +
  ggplot2::theme_minimal(
    base_size = 12,
    base_family = BASE_FAMILY
  ) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(face = "bold")
  )

export_plot(p_top, "top_guide_enrichment_validation", 6.6, 4.8)

# ============================================================
# DIAGNOSTIC ABLATION PLOT
# ============================================================

model_summary_plot <- model_summary %>%
  dplyr::mutate(
    model_label = dplyr::recode(
      as.character(model),
      sequence_only = "Sequence only",
      sequence_offset = "Sequence + offset",
      sequence_offset_offtarget = "Sequence + offset + off-target"
    ),
    model_label = factor(
      model_label,
      levels = c(
        "Sequence only",
        "Sequence + offset",
        "Sequence + offset + off-target"
      )
    )
  )

p_ablation <- ggplot2::ggplot(
  model_summary_plot,
  ggplot2::aes(x = model_label, y = mean_mse)
) +
  ggplot2::geom_col(fill = "grey70", color = "black") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean_mse - se_mse, ymax = mean_mse + se_mse),
    width = 0.15
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Mean held-out validation MSE"
  ) +
  ggplot2::theme_minimal(
    base_size = 12,
    base_family = BASE_FAMILY
  ) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
    axis.title.y = ggplot2::element_text(face = "bold")
  )

export_plot(p_ablation, "model_feature_ablation", 7.0, 4.8)

# ============================================================
# SAVE MODEL BUNDLE
# ============================================================

model_feature_names <- colnames(X_deployed)

model_bundle <- list(
  model = list(
    cvfit = final_fit$cvfit,
    lambda = final_fit$lambda,
    lambda_choice = LAMBDA_CHOICE,
    alpha = ALPHA_LASSO,
    deployed_feature_set = deployed_model,
    feature_names = model_feature_names,
    spacer_length = SPACER_LENGTH,
    downstream_length = DOWNSTREAM_LENGTH,
    use_offset = "offset_log1p" %in% model_feature_names,
    use_offtarget = "offtarget_log1p" %in% model_feature_names,
    target_definition = paste0(
      "median across control profiles of within-gene centered ",
      "efficacy, where efficacy=-LFC.adj"
    )
  ),
  training = list(
    control_manifest = CONTROL_MANIFEST,
    p_stouffer_cut = P_STOUFFER_CUT,
    negative_lfc_cut = NEGATIVE_LFC_CUT,
    min_mean_T0 = MIN_MEAN_T0,
    min_batches_per_guide = MIN_BATCHES_PER_GUIDE,
    min_profiles_per_guide = MIN_PROFILES_PER_GUIDE,
    min_guides_per_gene = MIN_GUIDES_PER_GENE,
    n_guides = nrow(Lib),
    n_genes = data.table::uniqueN(Lib$PAO1_ID),
    training_gene_ids = sort(unique(Lib$PAO1_ID))
  ),
  validation = list(
    outer_metrics = outer_metrics,
    model_summary = model_summary,
    representative_repeat = representative_repeat,
    representative_validation_genes = sort(
      representative_split$validation_genes
    ),
    representative_metrics = validation_metrics
  ),
  encoding = list(
    reference_base = "T",
    continuous_features = c(
      "offset_log1p = log1p(max(offset,0))",
      "offtarget_log1p = log1p(max(off_target_n,0))"
    )
  )
)

MODEL_FILE <- file.path(OUT_DIR, "CRISPRi_guide_efficacy_model.rds")
saveRDS(model_bundle, MODEL_FILE, compress = "xz")

# ============================================================
# GENERIC SCORING FUNCTION
# ============================================================

score_guides_with_bundle <- function(df, bundle) {
  dt <- data.table::as.data.table(data.table::copy(df))
  dt <- standardize_guide_columns(dt)
  
  required <- c("PAO1_ID", "spacer", "pam", "downstream")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("Guide table missing: ", paste(missing, collapse = ", "))
  }
  
  if (!"offset" %in% names(dt)) {
    dt[, offset := NA_real_]
  }
  if (!"off_target_n" %in% names(dt)) {
    dt[, off_target_n := NA_real_]
  }
  
  dt[
    ,
    `:=`(
      spacer = toupper(as.character(spacer)),
      pam = toupper(as.character(pam)),
      downstream = toupper(as.character(downstream)),
      offset = suppressWarnings(as.numeric(offset)),
      off_target_n = suppressWarnings(as.numeric(off_target_n)),
      downstream_model = substr(
        toupper(as.character(downstream)),
        1,
        bundle$model$downstream_length
      )
    )
  ]
  
  dt[
    ,
    score_status := data.table::fcase(
      !valid_dna(spacer, bundle$model$spacer_length),
      "invalid_spacer",
      is.na(pam) | !stringr::str_detect(pam, "^[ACGT]GG$"),
      "invalid_pam",
      !valid_dna(downstream_model, bundle$model$downstream_length),
      "invalid_downstream_context",
      !is.finite(offset) & bundle$model$use_offset,
      "missing_offset",
      !is.finite(off_target_n) & bundle$model$use_offtarget,
      "missing_offtarget",
      default = "scorable"
    )
  ]
  
  dt[, predicted_efficacy_score := NA_real_]
  
  score_idx <- which(dt$score_status == "scorable")
  if (length(score_idx) > 0) {
    score_input <- data.table::copy(dt[score_idx])
    
    # Unused metadata are set to zero only to construct the full
    # candidate matrix; the corresponding columns are removed
    # before prediction for models that do not use them.
    if (!bundle$model$use_offset) {
      score_input[, offset := 0]
    }
    if (!bundle$model$use_offtarget) {
      score_input[, off_target_n := 0]
    }
    
    X_all_score <- build_full_design_matrix(score_input)
    
    X_score <- switch(
      bundle$model$deployed_feature_set,
      sequence_only = MODEL_FEATURE_SETS$sequence_only(X_all_score),
      sequence_offset = MODEL_FEATURE_SETS$sequence_offset(X_all_score),
      sequence_offset_offtarget =
        MODEL_FEATURE_SETS$sequence_offset_offtarget(X_all_score),
      stop(
        "Unknown deployed feature set in bundle: ",
        bundle$model$deployed_feature_set
      )
    )
    
    missing_features <- setdiff(
      bundle$model$feature_names,
      colnames(X_score)
    )
    if (length(missing_features) > 0) {
      X_score[missing_features] <- 0
    }
    
    X_score <- X_score[, bundle$model$feature_names, drop = FALSE]
    
    dt$predicted_efficacy_score[score_idx] <- as.numeric(
      predict(
        bundle$model$cvfit,
        newx = as.matrix(X_score),
        s = bundle$model$lambda,
        type = "link"
      )
    )
  }
  
  dt[
    ,
    predicted_score_centered_within_gene :=
      predicted_efficacy_score -
      stats::median(predicted_efficacy_score, na.rm = TRUE),
    by = PAO1_ID
  ]
  
  dt[
    ,
    predicted_rank_within_gene := data.table::frank(
      -predicted_efficacy_score,
      ties.method = "min",
      na.last = "keep"
    ),
    by = PAO1_ID
  ]
  
  dt
}

# ============================================================
# SCORE ALL GUIDES USED IN LIB1
# ============================================================

lib1_file <- first_existing(
  LIB1_GUIDE_CANDIDATES,
  "canonical B1 Lib1 guide table",
  required = TRUE
)

lib1_scored <- NULL

if (!is.na(lib1_file)) {
  message("Scoring Lib1 guides from: ", lib1_file)
  
  lib1_all <- standardize_guide_columns(data.table::fread(lib1_file))
  
  # The selected annotated raw-count table is the canonical source. Fail early
  # rather than silently falling back to a separate design-table lookup.
  if (!"off_target_n" %in% names(lib1_all)) {
    stop(
      "The canonical Lib1 input is missing off_target_n: ", lib1_file,
      "\nRerun the metacounter with the enriched PAO1 guide key before scoring."
    )
  }
  lib1_all[, off_target_n := suppressWarnings(as.numeric(off_target_n))]
  
  if (!"type" %in% names(lib1_all)) {
    lib1_all[, type := "perfect"]
  }
  if (!"original" %in% names(lib1_all)) {
    lib1_all[, original := NA_character_]
  }
  
  perfect_missing_offtarget <- lib1_all[
    tolower(as.character(type)) == "perfect" &
      !is.finite(off_target_n)
  ]
  
  message(
    "Lib1 perfect-guide off_target_n coverage: ",
    sum(
      tolower(as.character(lib1_all$type)) == "perfect" &
        is.finite(lib1_all$off_target_n)
    ),
    " / ",
    sum(tolower(as.character(lib1_all$type)) == "perfect")
  )
  
  if (isTRUE(STRICT_LIB1_OFFTARGET_COVERAGE) &&
      nrow(perfect_missing_offtarget) > 0L) {
    missing_examples <- unique(
      perfect_missing_offtarget[, .(PAO1_ID, spacer)]
    )
    missing_examples <- utils::head(missing_examples, 10L)
    
    stop(
      nrow(perfect_missing_offtarget),
      " perfect Lib1 guide row(s) lack off_target_n in the canonical B1 input. ",
      "Example guide(s):\n",
      paste(
        paste0(missing_examples$PAO1_ID, " : ", missing_examples$spacer),
        collapse = "\n"
      )
    )
  }
  
  # Score perfect guides directly. Mismatch guides are not scored
  # as independent perfect guides; they receive the score of their
  # full-match parent in a separate column.
  perfect_idx <- which(tolower(as.character(lib1_all$type)) == "perfect")
  
  lib1_all[, predicted_efficacy_score := NA_real_]
  lib1_all[, score_status := "not_scored_nonperfect"]
  lib1_all[, predicted_score_centered_within_gene := NA_real_]
  lib1_all[, predicted_rank_within_gene := NA_real_]
  
  if (length(perfect_idx) > 0) {
    perfect_scored <- score_guides_with_bundle(
      lib1_all[perfect_idx],
      model_bundle
    )
    
    lib1_all$predicted_efficacy_score[perfect_idx] <-
      perfect_scored$predicted_efficacy_score
    lib1_all$score_status[perfect_idx] <- perfect_scored$score_status
    lib1_all$predicted_score_centered_within_gene[perfect_idx] <-
      perfect_scored$predicted_score_centered_within_gene
    lib1_all$predicted_rank_within_gene[perfect_idx] <-
      perfect_scored$predicted_rank_within_gene
  }
  
  parent_lookup <- unique(lib1_all[
    tolower(as.character(type)) == "perfect" &
      !is.na(predicted_efficacy_score),
    .(
      PAO1_ID,
      parent_spacer = as.character(spacer),
      parent_predicted_efficacy_score = predicted_efficacy_score,
      parent_predicted_rank_within_gene = predicted_rank_within_gene
    )
  ])
  
  lib1_all[
    ,
    parent_spacer := data.table::fcase(
      tolower(as.character(type)) == "perfect",
      as.character(spacer),
      !is.na(original) & as.character(original) != "",
      as.character(original),
      default = NA_character_
    )
  ]
  
  lib1_scored <- merge(
    lib1_all,
    parent_lookup,
    by = c("PAO1_ID", "parent_spacer"),
    all.x = TRUE
  )
  
  data.table::fwrite(
    lib1_scored,
    file.path(OUT_DIR, "Lib1_all_guides_with_model_scores.tsv.gz"),
    sep = "\t"
  )
  
  data.table::fwrite(
    lib1_scored,
    file.path(OUT_DIR, "Lib1_all_guides_with_model_scores.csv.gz"),
    sep = ","
  )
} else {
  warning(
    "No complete Lib1 guide table was found; Lib1 scoring was skipped."
  )
}

# ============================================================
# SCORE COMPLETE PAO1 GUIDE-DESIGN CATALOGUE(S)
# ============================================================

existing_design_files <- unique(
  normalizePath(
    PAO1_DESIGN_CANDIDATES[file.exists(PAO1_DESIGN_CANDIDATES)],
    winslash = "/",
    mustWork = TRUE
  )
)

pao1_score_summary <- list()

if (length(existing_design_files) == 0) {
  warning(
    "No PAO1 all-guide design catalogue was found. Checked:\n",
    paste(PAO1_DESIGN_CANDIDATES, collapse = "\n")
  )
} else {
  for (design_file in existing_design_files) {
    message("Scoring PAO1 design catalogue: ", design_file)
    
    design_tbl <- standardize_guide_columns(
      data.table::fread(design_file)
    )
    
    scored_tbl <- score_guides_with_bundle(design_tbl, model_bundle)
    
    stem <- sanitize_stem(design_file)
    output_file <- file.path(
      OUT_DIR,
      paste0(stem, "_with_model_scores.tsv.gz")
    )
    
    data.table::fwrite(scored_tbl, output_file, sep = "\t")
    
    pao1_score_summary[[stem]] <- data.table::data.table(
      source_file = design_file,
      output_file = output_file,
      n_rows = nrow(scored_tbl),
      n_genes = data.table::uniqueN(scored_tbl$PAO1_ID),
      n_scorable = sum(scored_tbl$score_status == "scorable"),
      n_scored = sum(is.finite(scored_tbl$predicted_efficacy_score))
    )
  }
}

pao1_score_summary_tbl <- data.table::rbindlist(
  pao1_score_summary,
  use.names = TRUE,
  fill = TRUE
)

# ============================================================
# EXPORT REUSABLE TABLES
# ============================================================

data.table::fwrite(
  profile_long,
  file.path(OUT_DIR, "model_control_profiles_long.tsv.gz"),
  sep = "\t"
)

data.table::fwrite(
  gene_calls,
  file.path(OUT_DIR, "model_training_gene_calls.tsv.gz"),
  sep = "\t"
)

data.table::fwrite(
  Lib,
  file.path(OUT_DIR, "model_training_guides_with_scores.tsv.gz"),
  sep = "\t"
)

data.table::fwrite(
  outer_metrics,
  file.path(OUT_DIR, "repeated_gene_heldout_validation_metrics.tsv"),
  sep = "\t"
)

data.table::fwrite(
  model_summary,
  file.path(OUT_DIR, "model_ablation_summary.tsv"),
  sep = "\t"
)

data.table::fwrite(
  validation_df,
  file.path(OUT_DIR, "representative_validation_guide_data.tsv.gz"),
  sep = "\t"
)

data.table::fwrite(
  validation_gene_cor,
  file.path(OUT_DIR, "representative_validation_gene_correlations.tsv"),
  sep = "\t"
)

data.table::fwrite(
  top_guide_enrichment,
  file.path(OUT_DIR, "top_guide_enrichment_validation.tsv"),
  sep = "\t"
)

data.table::fwrite(
  coefficient_tbl,
  file.path(OUT_DIR, "final_model_coefficients.tsv"),
  sep = "\t"
)

data.table::fwrite(
  cv_curve_df,
  file.path(OUT_DIR, "final_model_cv_curve.tsv"),
  sep = "\t"
)

if (nrow(pao1_score_summary_tbl) > 0) {
  data.table::fwrite(
    pao1_score_summary_tbl,
    file.path(OUT_DIR, "PAO1_design_scoring_summary.tsv"),
    sep = "\t"
  )
}

settings_tbl <- data.table::data.table(
  setting = c(
    "P_STOUFFER_CUT",
    "NEGATIVE_LFC_CUT",
    "MIN_MEAN_T0",
    "MIN_BATCHES_PER_GUIDE",
    "MIN_PROFILES_PER_GUIDE",
    "MIN_GUIDES_PER_GENE",
    "SPACER_LENGTH",
    "DOWNSTREAM_LENGTH",
    "ALPHA_LASSO",
    "INNER_NFOLDS",
    "OUTER_VALIDATION_FRACTION",
    "N_OUTER_REPEATS",
    "LAMBDA_CHOICE",
    "DEPLOYED_MODEL"
  ),
  value = as.character(c(
    P_STOUFFER_CUT,
    NEGATIVE_LFC_CUT,
    MIN_MEAN_T0,
    MIN_BATCHES_PER_GUIDE,
    MIN_PROFILES_PER_GUIDE,
    MIN_GUIDES_PER_GENE,
    SPACER_LENGTH,
    DOWNSTREAM_LENGTH,
    ALPHA_LASSO,
    INNER_NFOLDS,
    OUTER_VALIDATION_FRACTION,
    N_OUTER_REPEATS,
    LAMBDA_CHOICE,
    deployed_model
  ))
)

wb <- openxlsx::createWorkbook()
add_sheet(wb, "settings", settings_tbl)
add_sheet(wb, "control_manifest", CONTROL_MANIFEST)
add_sheet(wb, "model_summary", model_summary)
add_sheet(wb, "validation_metrics", validation_metrics)
add_sheet(wb, "validation_gene_cor", validation_gene_cor)
add_sheet(wb, "top_guide_enrichment", top_guide_enrichment)
add_sheet(wb, "coefficients", coefficient_tbl)
add_sheet(wb, "training_gene_calls", gene_calls)
add_sheet(wb, "PAO1_scoring_summary", pao1_score_summary_tbl)

openxlsx::saveWorkbook(
  wb,
  file.path(OUT_DIR, "guide_efficacy_model_results.xlsx"),
  overwrite = TRUE
)

message("\nGuide-efficacy model workflow completed.")
message("Deployed model: ", deployed_model)
message("Model bundle: ", MODEL_FILE)
message("Output directory: ", normalizePath(OUT_DIR, winslash = "/"))

# Reprint manuscript-facing plots for RStudio navigation.
print(p_logo)
print(p_cv)
print(p_validation)
print(p_gene_cor)
print(p_top)

# Complete reusable inputs and underlying panel data produced by this run.
message("Training guide data: ", file.path(OUT_DIR, "model_training_guides_with_scores.tsv.gz"))
message("Panel C coefficients: ", file.path(OUT_DIR, "final_model_coefficients.tsv"))
message("Panel C CV curve: ", file.path(OUT_DIR, "final_model_cv_curve.tsv"))
message("Panels D-E guide data: ", file.path(OUT_DIR, "representative_validation_guide_data.tsv.gz"))
message("Panel F data: ", file.path(OUT_DIR, "top_guide_enrichment_validation.tsv"))
message("Repeated evaluation: ", file.path(OUT_DIR, "repeated_gene_heldout_validation_metrics.tsv"))
