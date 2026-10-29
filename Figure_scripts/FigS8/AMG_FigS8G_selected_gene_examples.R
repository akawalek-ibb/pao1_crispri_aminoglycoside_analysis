#!/usr/bin/env Rscript
# Compare three guide models for the selected efp, lptF and parE examples.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/FigS8/AMG_FigS8G_selected_gene_examples.R")
# Inputs:
# - guide_efficacy_prediction/this_study/CRISPRi_guide_efficacy_model.rds and
#   model_training_guides_with_scores.tsv.gz: matched model/training population.
# - Figure_scripts/FigS8/input_AMG_FigS8/PAO1_sgRNA_with_Pasteur_scores.tsv.gz and
#   PAO1_all_sgRNA_predictions_CIAO.csv.
# Upstream / execution order:
# - Use the supplied matched model/training pair and external prediction tables.
# - To use a new study model, run AMG_FigS8C_F_train_evaluate_model.R and explicitly select its
#   output with MODEL_DIR.
# Outputs:
# - Figure_scripts/FigS8/out_AMG_FigS8G_selected_genes/ — common-guide predictions, correlations
#   and PNG/PDF.
# Interpretation / operation:
# - Exclude each selected gene before refitting and retuning lambda; compare all models on its
#   same common guides.
# - Scores are oriented so larger values indicate stronger predicted activity. This is separate
#   from deployment scoring.
#
# R packages used/checked: data.table, ggplot2, glmnet, patchwork.

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(ggplot2)
  library(patchwork)
})

# ============================================================
# USER SETTINGS
# ============================================================

PROJECT_DIR <- "."

MODEL_DIR <- "guide_efficacy_prediction/this_study"

EXPERIMENTAL_GUIDES_FILE <- file.path(
  MODEL_DIR,
  "model_training_guides_with_scores.tsv.gz"
)

MODEL_BUNDLE_FILE <- file.path(
  MODEL_DIR,
  "CRISPRi_guide_efficacy_model.rds"
)

PASTEUR_FILE <- file.path("Figure_scripts", "FigS8", "input_AMG_FigS8", "PAO1_sgRNA_with_Pasteur_scores.tsv.gz")

CIAO_FILE <- file.path("Figure_scripts", "FigS8", "input_AMG_FigS8", "PAO1_all_sgRNA_predictions_CIAO.csv")

# Gene names OR PAO1 locus tags can be supplied here.
# Matching of gene names is case-insensitive.
SELECTED_GENES <- c(
  "efp",
  "lptF",
  "parE"
)

OUT_DIR <- file.path("Figure_scripts", "FigS8", "out_AMG_FigS8G_selected_genes")

# Inner grouped CV used when refitting the model after excluding
# each selected gene.
INNER_NFOLDS <- 10L
BASE_SEED <- 4201L

# Plot settings.
BASE_FAMILY <- "Arial"
BASE_SIZE <- 11
PLOT_DPI <- 300
POINT_SIZE <- 2.0
POINT_ALPHA <- 0.75
FACET_SCALES <- "free"

# TRUE is recommended for the manuscript preview: every model
# row then contains exactly the same guides within each gene.
REQUIRE_COMMON_GUIDES_ALL_MODELS <- TRUE

MODEL_ORDER <- c(
  "P. aeruginosa-trained model",
  "CRISPRi-Pasteur",
  "CIAO"
)

# ============================================================
# SETUP
# ============================================================

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

required_files <- c(
  EXPERIMENTAL_GUIDES_FILE,
  MODEL_BUNDLE_FILE,
  PASTEUR_FILE,
  CIAO_FILE
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required file(s) were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ============================================================
# HELPERS
# ============================================================

normalize_dna <- function(x) {
  toupper(gsub("\\s+", "", as.character(x)))
}

normalize_text <- function(x) {
  trimws(as.character(x))
}

as_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

first_existing_column <- function(
    dt,
    candidates,
    required = TRUE,
    description = "column"
) {
  hit <- candidates[candidates %in% names(dt)]
  
  if (length(hit) == 0L) {
    if (isTRUE(required)) {
      stop(
        "Could not identify ", description, ". Checked: ",
        paste(candidates, collapse = ", "),
        "\nAvailable columns:\n",
        paste(names(dt), collapse = ", ")
      )
    }
    return(NULL)
  }
  
  hit[[1L]]
}

is_scorable_status <- function(x) {
  x <- tolower(trimws(as.character(x)))
  
  is.na(x) |
    x == "" |
    x %in% c("scorable", "scored", "ok", "success")
}

make_guide_key <- function(locus_tag, spacer, pam) {
  paste(
    normalize_text(locus_tag),
    normalize_dna(spacer),
    normalize_dna(pam),
    sep = "|"
  )
}

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  
  if (sum(keep) < 3L) {
    return(NA_real_)
  }
  
  x <- x[keep]
  y <- y[keep]
  
  if (
    length(unique(x)) < 2L ||
    length(unique(y)) < 2L
  ) {
    return(NA_real_)
  }
  
  suppressWarnings(
    stats::cor(
      x,
      y,
      method = method,
      use = "complete.obs"
    )
  )
}

valid_dna <- function(x, expected_length) {
  x <- normalize_dna(x)
  
  !is.na(x) &
    nchar(x) == expected_length &
    grepl(
      paste0("^[ACGT]{", expected_length, "}$"),
      x
    )
}

one_hot_acg <- function(seqs, prefix) {
  seqs <- normalize_dna(seqs)
  lengths <- unique(nchar(seqs))
  
  if (length(lengths) != 1L) {
    stop("Sequences supplied to one_hot_acg() have unequal lengths.")
  }
  
  npos <- lengths[[1L]]
  
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
  
  for (position in seq_len(npos)) {
    base <- substr(seqs, position, position)
    first_column <- (position - 1L) * 3L
    
    X[, first_column + 1L] <- as.numeric(base == "A")
    X[, first_column + 2L] <- as.numeric(base == "C")
    X[, first_column + 3L] <- as.numeric(base == "G")
  }
  
  as.data.frame(
    X,
    check.names = FALSE
  )
}

build_full_design_matrix <- function(df) {
  required <- c(
    "spacer",
    "pam",
    "downstream_model",
    "offset",
    "off_target_n"
  )
  
  missing <- setdiff(required, names(df))
  
  if (length(missing) > 0L) {
    stop(
      "Missing model feature columns: ",
      paste(missing, collapse = ", ")
    )
  }
  
  spacer_X <- one_hot_acg(
    df$spacer,
    "spacer_"
  )
  
  downstream_X <- one_hot_acg(
    df$downstream_model,
    "down_"
  )
  
  pam_n <- substr(df$pam, 1L, 1L)
  
  pam_X <- data.frame(
    pam_1_A = as.numeric(pam_n == "A"),
    pam_1_C = as.numeric(pam_n == "C"),
    pam_1_G = as.numeric(pam_n == "G"),
    check.names = FALSE
  )
  
  meta_X <- data.frame(
    offset_log1p = log1p(
      pmax(as.numeric(df$offset), 0)
    ),
    offtarget_log1p = log1p(
      pmax(as.numeric(df$off_target_n), 0)
    ),
    check.names = FALSE
  )
  
  data.frame(
    spacer_X,
    pam_X,
    downstream_X,
    meta_X,
    check.names = FALSE
  )
}

select_deployed_features <- function(
    X,
    deployed_feature_set
) {
  switch(
    deployed_feature_set,
    
    sequence_only = X[
      ,
      !names(X) %in% c(
        "offset_log1p",
        "offtarget_log1p"
      ),
      drop = FALSE
    ],
    
    sequence_offset = X[
      ,
      names(X) != "offtarget_log1p",
      drop = FALSE
    ],
    
    sequence_offset_offtarget = X,
    
    stop(
      "Unknown deployed feature set: ",
      deployed_feature_set
    )
  )
}

align_feature_matrix <- function(
    X,
    feature_names
) {
  missing_features <- setdiff(
    feature_names,
    colnames(X)
  )
  
  if (length(missing_features) > 0L) {
    for (feature in missing_features) {
      X[[feature]] <- 0
    }
  }
  
  extra_features <- setdiff(
    colnames(X),
    feature_names
  )
  
  if (length(extra_features) > 0L) {
    X <- X[
      ,
      setdiff(colnames(X), extra_features),
      drop = FALSE
    ]
  }
  
  X[
    ,
    feature_names,
    drop = FALSE
  ]
}

make_group_folds <- function(
    groups,
    nfolds,
    seed
) {
  groups <- as.character(groups)
  unique_groups <- unique(groups)
  
  if (length(unique_groups) < 3L) {
    stop(
      "At least three training genes are required for grouped CV."
    )
  }
  
  nfolds_use <- min(
    as.integer(nfolds),
    length(unique_groups)
  )
  
  set.seed(seed)
  shuffled <- sample(unique_groups)
  
  group_fold <- rep(
    seq_len(nfolds_use),
    length.out = length(shuffled)
  )
  
  names(group_fold) <- shuffled
  as.integer(group_fold[groups])
}

fit_gene_heldout_model <- function(
    X_train,
    y_train,
    training_groups,
    bundle,
    seed
) {
  foldid <- make_group_folds(
    groups = training_groups,
    nfolds = INNER_NFOLDS,
    seed = seed
  )
  
  cvfit <- glmnet::cv.glmnet(
    x = as.matrix(X_train),
    y = as.numeric(y_train),
    family = "gaussian",
    alpha = bundle$model$alpha,
    standardize = TRUE,
    intercept = TRUE,
    foldid = foldid,
    keep = FALSE
  )
  
  lambda <- switch(
    bundle$model$lambda_choice,
    "lambda.min" = cvfit$lambda.min,
    "lambda.1se" = cvfit$lambda.1se,
    stop(
      "Unknown lambda choice in model bundle: ",
      bundle$model$lambda_choice
    )
  )
  
  list(
    cvfit = cvfit,
    lambda = lambda
  )
}

format_number <- function(x, digits = 2L) {
  if (!is.finite(x)) {
    return("NA")
  }
  
  sprintf(
    paste0("%.", digits, "f"),
    x
  )
}

save_plot <- function(
    plot_object,
    stem,
    width,
    height
) {
  ggplot2::ggsave(
    filename = file.path(
      OUT_DIR,
      paste0(stem, ".png")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = PLOT_DPI,
    bg = "white"
  )
  
  ggplot2::ggsave(
    filename = file.path(
      OUT_DIR,
      paste0(stem, ".pdf")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}

# ============================================================
# LOAD TRAINING GUIDE TABLE AND MODEL BUNDLE
# ============================================================

message(
  "Reading experimentally measured guide data:\n  ",
  EXPERIMENTAL_GUIDES_FILE
)

# 1. Read the matched experimental guide population and saved model bundle.
# The bundle supplies feature definitions and tuning choices; deployment scores
# from a model fitted to the selected gene are not used as its held-out prediction.
guides <- data.table::fread(
  EXPERIMENTAL_GUIDES_FILE,
  header = TRUE,
  showProgress = TRUE
)

bundle <- readRDS(
  MODEL_BUNDLE_FILE
)

required_columns <- c(
  "PAO1_ID",
  "gene_name",
  "spacer",
  "pam",
  "downstream",
  "offset",
  "off_target_n",
  "y",
  "LFC_agg"
)

missing_columns <- setdiff(
  required_columns,
  names(guides)
)

if (length(missing_columns) > 0L) {
  stop(
    "Experimental guide table is missing required column(s):\n",
    paste(missing_columns, collapse = ", ")
  )
}

guides[
  ,
  `:=`(
    PAO1_ID = trimws(as.character(PAO1_ID)),
    gene_name = trimws(as.character(gene_name)),
    spacer = normalize_dna(spacer),
    pam = normalize_dna(pam),
    downstream = normalize_dna(downstream),
    offset = suppressWarnings(as.numeric(offset)),
    off_target_n = suppressWarnings(as.numeric(off_target_n)),
    y = suppressWarnings(as.numeric(y)),
    LFC_agg = suppressWarnings(as.numeric(LFC_agg)),
    guide_key = make_guide_key(PAO1_ID, spacer, pam)
  )
]

guides[
  ,
  downstream_model := substr(
    downstream,
    1L,
    bundle$model$downstream_length
  )
]

valid_rows <-
  valid_dna(
    guides$spacer,
    bundle$model$spacer_length
  ) &
  grepl("^[ACGT]GG$", guides$pam) &
  valid_dna(
    guides$downstream_model,
    bundle$model$downstream_length
  ) &
  is.finite(guides$y) &
  is.finite(guides$LFC_agg) &
  (
    !isTRUE(bundle$model$use_offset) |
      is.finite(guides$offset)
  ) &
  (
    !isTRUE(bundle$model$use_offtarget) |
      is.finite(guides$off_target_n)
  )

guides <- guides[valid_rows]
guides[, experimental_row_id__ := .I]

if (nrow(guides) == 0L) {
  stop("No valid guide rows remained after model compatibility checks.")
}

# ============================================================
# LOAD PASTEUR AND CIAO SCORES
# ============================================================

message("Reading CRISPRi-Pasteur scores:\n  ", PASTEUR_FILE)

pasteur_raw <- data.table::fread(
  PASTEUR_FILE,
  header = TRUE,
  showProgress = TRUE,
  check.names = FALSE
)

pasteur_locus_col <- first_existing_column(
  pasteur_raw,
  c("locus_tag", "PAO1_ID"),
  description = "Pasteur locus-tag column"
)

pasteur_guide_col <- first_existing_column(
  pasteur_raw,
  c("sgrna", "guide", "spacer", "gRNA sequence (20nt)"),
  description = "Pasteur guide column"
)

pasteur_pam_col <- first_existing_column(
  pasteur_raw,
  c("pam", "PAM"),
  description = "Pasteur PAM column"
)

pasteur_score_col <- first_existing_column(
  pasteur_raw,
  c("pasteur_activity_score", "Pasteur_activity_score", "Prediction", "score"),
  description = "Pasteur activity-score column"
)

pasteur_status_col <- first_existing_column(
  pasteur_raw,
  c("pasteur_score_status", "score_status"),
  required = FALSE
)

pasteur <- data.table(
  PAO1_ID = normalize_text(pasteur_raw[[pasteur_locus_col]]),
  spacer = normalize_dna(pasteur_raw[[pasteur_guide_col]]),
  pam = normalize_dna(pasteur_raw[[pasteur_pam_col]]),
  pasteur_activity_score = as_numeric_safe(pasteur_raw[[pasteur_score_col]]),
  pasteur_score_status = if (is.null(pasteur_status_col)) {
    rep("scorable", nrow(pasteur_raw))
  } else {
    as.character(pasteur_raw[[pasteur_status_col]])
  }
)

pasteur[, guide_key := make_guide_key(PAO1_ID, spacer, pam)]
pasteur[, duplicate_n__ := .N, by = guide_key]
pasteur[, status_priority__ := as.integer(!is_scorable_status(pasteur_score_status))]
pasteur[, missing_score_priority__ := as.integer(!is.finite(pasteur_activity_score))]
setorder(
  pasteur,
  guide_key,
  status_priority__,
  missing_score_priority__
)
pasteur <- pasteur[, .SD[1L], by = guide_key]
pasteur[, c("status_priority__", "missing_score_priority__") := NULL]

message("Reading CIAO scores:\n  ", CIAO_FILE)

ciao_raw <- data.table::fread(
  CIAO_FILE,
  header = TRUE,
  showProgress = TRUE,
  check.names = FALSE
)

ciao_locus_col <- first_existing_column(
  ciao_raw,
  c("PAO1_ID", "locus_tag"),
  description = "CIAO locus-tag column"
)

ciao_guide_col <- first_existing_column(
  ciao_raw,
  c("gRNA sequence (20nt)", "sgrna", "spacer"),
  description = "CIAO guide column"
)

ciao_pam_col <- first_existing_column(
  ciao_raw,
  c("PAM", "pam"),
  description = "CIAO PAM column"
)

ciao_score_col <- first_existing_column(
  ciao_raw,
  c("Activity score", "activity_score"),
  description = "CIAO activity-score column"
)

ciao <- data.table(
  PAO1_ID = normalize_text(ciao_raw[[ciao_locus_col]]),
  spacer = normalize_dna(ciao_raw[[ciao_guide_col]]),
  pam = normalize_dna(ciao_raw[[ciao_pam_col]]),
  ciao_activity_score = as_numeric_safe(ciao_raw[[ciao_score_col]])
)

ciao[, guide_key := make_guide_key(PAO1_ID, spacer, pam)]
ciao[, duplicate_n__ := .N, by = guide_key]
ciao[, missing_score_priority__ := as.integer(!is.finite(ciao_activity_score))]
setorder(ciao, guide_key, missing_score_priority__, ciao_activity_score)
ciao <- ciao[, .SD[1L], by = guide_key]
ciao[, missing_score_priority__ := NULL]

guides <- merge(
  guides,
  pasteur[, .(guide_key, pasteur_activity_score, pasteur_score_status)],
  by = "guide_key",
  all.x = TRUE,
  sort = FALSE
)

guides <- merge(
  guides,
  ciao[, .(guide_key, ciao_activity_score)],
  by = "guide_key",
  all.x = TRUE,
  sort = FALSE
)
setorder(guides, experimental_row_id__)

guides[
  ,
  `:=`(
    pasteur_available = is.finite(pasteur_activity_score) &
      is_scorable_status(pasteur_score_status),
    ciao_available = is.finite(ciao_activity_score)
  )
]

message(
  "External-score coverage among eligible experimental guides: ",
  "Pasteur ", sum(guides$pasteur_available), "/", nrow(guides),
  "; CIAO ", sum(guides$ciao_available), "/", nrow(guides)
)

# ============================================================
# RESOLVE REQUESTED GENE NAMES / LOCUS TAGS
# ============================================================

available_gene_map <- unique(
  guides[
    ,
    .(
      PAO1_ID,
      gene_name
    )
  ]
)

resolve_one_gene <- function(query) {
  query_clean <- trimws(as.character(query))
  
  locus_hit <- available_gene_map[
    toupper(PAO1_ID) == toupper(query_clean)
  ]
  
  name_hit <- available_gene_map[
    !is.na(gene_name) &
      gene_name != "" &
      tolower(gene_name) == tolower(query_clean)
  ]
  
  hit <- unique(
    rbindlist(
      list(locus_hit, name_hit),
      use.names = TRUE,
      fill = TRUE
    )
  )
  
  if (nrow(hit) == 0L) {
    partial <- available_gene_map[
      grepl(
        tolower(query_clean),
        tolower(paste(PAO1_ID, gene_name)),
        fixed = TRUE
      )
    ]
    
    suggestion_text <- if (nrow(partial) == 0L) {
      "none"
    } else {
      paste(
        utils::head(
          paste0(
            partial$gene_name,
            " (",
            partial$PAO1_ID,
            ")"
          ),
          10L
        ),
        collapse = ", "
      )
    }
    
    stop(
      "Could not resolve selected gene: ",
      query_clean,
      "\nPossible partial matches: ",
      suggestion_text
    )
  }
  
  if (data.table::uniqueN(hit$PAO1_ID) > 1L) {
    stop(
      "Selected gene name maps to more than one locus: ",
      query_clean,
      "\nMatches:\n",
      paste(
        paste0(
          hit$gene_name,
          " (",
          hit$PAO1_ID,
          ")"
        ),
        collapse = "\n"
      ),
      "\nUse the PAO1 locus tag explicitly."
    )
  }
  
  data.table(
    requested_gene = query_clean,
    PAO1_ID = hit$PAO1_ID[[1L]],
    gene_name = hit$gene_name[[1L]]
  )
}

resolved_genes <- rbindlist(
  lapply(
    SELECTED_GENES,
    resolve_one_gene
  ),
  use.names = TRUE,
  fill = TRUE
)

resolved_genes <- unique(
  resolved_genes,
  by = "PAO1_ID"
)

message("\nResolved selected genes:")
print(resolved_genes)

# ============================================================
# BUILD MODEL MATRIX ONCE
# ============================================================

# 2. Encode the guide/context features exactly as specified by the saved
# bundle. Match feature names/order before predicting to avoid mixing columns.
X_full <- build_full_design_matrix(
  guides
)

X_deployed <- select_deployed_features(
  X_full,
  bundle$model$deployed_feature_set
)

X_deployed <- align_feature_matrix(
  X_deployed,
  bundle$model$feature_names
)

# ============================================================
# LEAVE-EACH-SELECTED-GENE-OUT PREDICTIONS
# ============================================================

prediction_list <- vector(
  "list",
  nrow(resolved_genes)
)

# 3. Exclude every guide of the selected gene, refit using the remaining genes
# and retune lambda with gene-grouped inner cross-validation. Repeat independently
# for each selected gene, then predict the excluded gene's guides.
for (gene_i in seq_len(nrow(resolved_genes))) {
  gene_id <- resolved_genes$PAO1_ID[[gene_i]]
  gene_name_value <- resolved_genes$gene_name[[gene_i]]
  
  message(
    "\nPredicting selected gene ",
    gene_name_value,
    " (",
    gene_id,
    ") using a model trained on all other genes..."
  )
  
  test_index <- which(
    guides$PAO1_ID == gene_id
  )
  
  train_index <- which(
    guides$PAO1_ID != gene_id
  )
  
  if (length(test_index) < 3L) {
    warning(
      gene_name_value,
      " (",
      gene_id,
      ") has only ",
      length(test_index),
      " eligible guides; correlation estimates may be unavailable."
    )
  }
  
  heldout_fit <- fit_gene_heldout_model(
    X_train = X_deployed[
      train_index,
      ,
      drop = FALSE
    ],
    y_train = guides$y[train_index],
    training_groups = guides$PAO1_ID[train_index],
    bundle = bundle,
    seed = BASE_SEED + gene_i
  )
  
  predicted <- as.numeric(
    predict(
      heldout_fit$cvfit,
      newx = as.matrix(
        X_deployed[
          test_index,
          ,
          drop = FALSE
        ]
      ),
      s = heldout_fit$lambda,
      type = "link"
    )
  )
  
  gene_data <- copy(
    guides[test_index]
  )
  
  gene_data[
    ,
    `:=`(
      requested_gene = resolved_genes$requested_gene[[gene_i]],
      predicted_efficacy_heldout = predicted,
      heldout_lambda = heldout_fit$lambda,
      observed_efficacy = -LFC_agg
    )
  ]
  
  prediction_list[[gene_i]] <- gene_data
}

plot_data <- rbindlist(
  prediction_list,
  use.names = TRUE,
  fill = TRUE
)

plot_data[
  ,
  gene_label := ifelse(
    is.na(gene_name) | gene_name == "",
    PAO1_ID,
    paste0(gene_name, " (", PAO1_ID, ")")
  )
]

# Preserve requested order in facets.
facet_levels <- resolved_genes[
  ,
  ifelse(
    is.na(gene_name) | gene_name == "",
    PAO1_ID,
    paste0(gene_name, " (", PAO1_ID, ")")
  )
]

plot_data[
  ,
  gene_label := factor(
    gene_label,
    levels = facet_levels
  )
]

# ============================================================
# BUILD THREE-MODEL LONG TABLE
# ============================================================

plot_data[
  ,
  ciao_strength_score := -ciao_activity_score
]

if (isTRUE(REQUIRE_COMMON_GUIDES_ALL_MODELS)) {
  plot_data_common <- plot_data[
    is.finite(predicted_efficacy_heldout) &
      pasteur_available == TRUE &
      ciao_available == TRUE
  ]
} else {
  plot_data_common <- copy(plot_data)
}

if (nrow(plot_data_common) == 0L) {
  stop("No selected-gene guides had usable scores for plotting.")
}

# 4. Compare the study, Pasteur and CIAO predictions on the selected common
# guides when that option is enabled. Orient scores so larger values mean
# stronger activity; observed efficacy here is -LFC_agg, not transcript repression.
model_long <- rbindlist(
  list(
    plot_data_common[
      is.finite(predicted_efficacy_heldout),
      .(
        PAO1_ID,
        gene_name,
        gene_label,
        spacer,
        pam,
        downstream,
        offset,
        off_target_n,
        LFC_agg,
        observed_efficacy,
        y,
        n_batches,
        n_profiles,
        min_mean_T0_count,
        median_mean_T0_count,
        model = "P. aeruginosa-trained model",
        predicted_efficacy = predicted_efficacy_heldout,
        native_model_score = predicted_efficacy_heldout,
        native_score_direction = "higher = stronger",
        heldout_lambda = heldout_lambda
      )
    ],
    
    plot_data_common[
      pasteur_available == TRUE,
      .(
        PAO1_ID,
        gene_name,
        gene_label,
        spacer,
        pam,
        downstream,
        offset,
        off_target_n,
        LFC_agg,
        observed_efficacy,
        y,
        n_batches,
        n_profiles,
        min_mean_T0_count,
        median_mean_T0_count,
        model = "CRISPRi-Pasteur",
        predicted_efficacy = pasteur_activity_score,
        native_model_score = pasteur_activity_score,
        native_score_direction = "higher = stronger",
        heldout_lambda = NA_real_
      )
    ],
    
    plot_data_common[
      ciao_available == TRUE,
      .(
        PAO1_ID,
        gene_name,
        gene_label,
        spacer,
        pam,
        downstream,
        offset,
        off_target_n,
        LFC_agg,
        observed_efficacy,
        y,
        n_batches,
        n_profiles,
        min_mean_T0_count,
        median_mean_T0_count,
        model = "CIAO",
        predicted_efficacy = ciao_strength_score,
        native_model_score = ciao_activity_score,
        native_score_direction = "lower / more negative = stronger; sign inverted for plot",
        heldout_lambda = NA_real_
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

model_long[
  ,
  model := factor(
    model,
    levels = MODEL_ORDER
  )
]

# ============================================================
# PER-MODEL, PER-GENE CORRELATIONS
# ============================================================

correlation_stats <- model_long[
  ,
  .(
    n_guides = .N,
    
    # Raw-LFC orientation is retained in the exported table for
    # reference. Because stronger guides have more-negative LFC,
    # these correlations have the opposite sign from the plotted
    # efficacy-oriented correlations.
    spearman_predicted_efficacy_vs_LFC = safe_cor(
      predicted_efficacy,
      LFC_agg,
      "spearman"
    ),
    pearson_predicted_efficacy_vs_LFC = safe_cor(
      predicted_efficacy,
      LFC_agg,
      "pearson"
    ),
    
    # Efficacy orientation used in the plot and in S8G:
    # observed_efficacy = -LFC_agg, so positive correlation means
    # agreement between predicted and observed guide efficacy.
    spearman_predicted_vs_observed_efficacy = safe_cor(
      predicted_efficacy,
      observed_efficacy,
      "spearman"
    ),
    pearson_predicted_vs_observed_efficacy = safe_cor(
      predicted_efficacy,
      observed_efficacy,
      "pearson"
    ),
    
    min_observed_LFC = min(LFC_agg, na.rm = TRUE),
    max_observed_LFC = max(LFC_agg, na.rm = TRUE),
    min_predicted_efficacy = min(predicted_efficacy, na.rm = TRUE),
    max_predicted_efficacy = max(predicted_efficacy, na.rm = TRUE)
  ),
  by = .(
    model,
    PAO1_ID,
    gene_name,
    gene_label
  )
]

correlation_stats[
  ,
  annotation := paste0(
    "Spearman rho = ",
    vapply(
      spearman_predicted_vs_observed_efficacy,
      format_number,
      character(1),
      digits = 2L
    ),
    "\nPearson r = ",
    vapply(
      pearson_predicted_vs_observed_efficacy,
      format_number,
      character(1),
      digits = 2L
    ),
    "\nn = ",
    n_guides
  )
]

message("\nSelected-gene three-model correlations:")
print(
  correlation_stats[
    ,
    .(
      model,
      PAO1_ID,
      gene_name,
      n_guides,
      spearman_predicted_efficacy_vs_LFC,
      pearson_predicted_efficacy_vs_LFC,
      spearman_predicted_vs_observed_efficacy,
      pearson_predicted_vs_observed_efficacy
    )
  ]
)

# ============================================================
# COMBINED 3 x N-GENE PLOT
# ============================================================

make_model_row_plot <- function(model_name) {
  row_data <- model_long[
    as.character(model) == model_name
  ]
  
  row_stats <- correlation_stats[
    as.character(model) == model_name
  ]
  
  ggplot2::ggplot(
    row_data,
    ggplot2::aes(
      x = predicted_efficacy,
      y = observed_efficacy
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.35
    ) +
    ggplot2::geom_point(
      size = POINT_SIZE,
      alpha = POINT_ALPHA
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 0.65
    ) +
    ggplot2::geom_label(
      data = row_stats,
      ggplot2::aes(
        x = -Inf,
        y = Inf,
        label = annotation
      ),
      inherit.aes = FALSE,
      hjust = -0.04,
      vjust = 1.06,
      size = 2.7,
      label.size = 0.18,
      fill = "white"
    ) +
    ggplot2::facet_wrap(
      ~gene_label,
      nrow = 1,
      scales = FACET_SCALES
    ) +
    ggplot2::labs(
      title = model_name,
      x = "Predicted guide efficacy / activity",
      y = "Observed guide efficacy (-aggregated adjusted log2FC)"
    ) +
    ggplot2::theme_classic(
      base_family = BASE_FAMILY,
      base_size = BASE_SIZE
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = BASE_SIZE + 1
      ),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = BASE_SIZE - 1
      ),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(1.0, "lines")
    )
}

model_row_plots <- lapply(
  MODEL_ORDER,
  make_model_row_plot
)

p_selected <- patchwork::wrap_plots(
  model_row_plots,
  ncol = 1,
  heights = rep(1, length(MODEL_ORDER))
) +
  patchwork::plot_annotation(
    subtitle = paste0(
      "Columns show the same selected genes. Higher values on both axes indicate stronger guide activity; ",
      "positive correlation therefore indicates agreement between predicted and observed efficacy."
    ),
    theme = ggplot2::theme(
      plot.subtitle = ggplot2::element_text(
        family = BASE_FAMILY,
        size = BASE_SIZE
      )
    )
  )

combined_width <- max(
  10.5,
  3.8 * length(facet_levels)
)
combined_height <- 12.5

save_plot(
  p_selected,
  "selected_genes_observed_efficacy_vs_three_model_predicted_efficacy",
  width = combined_width,
  height = combined_height
)

# ============================================================
# EXPORT TABLES
# ============================================================

model_long_export_columns <- c(
  "model",
  "PAO1_ID",
  "gene_name",
  "spacer",
  "pam",
  "downstream",
  "offset",
  "off_target_n",
  "LFC_agg",
  "observed_efficacy",
  "y",
  "predicted_efficacy",
  "native_model_score",
  "native_score_direction",
  "heldout_lambda",
  "n_batches",
  "n_profiles",
  "min_mean_T0_count",
  "median_mean_T0_count"
)

model_long_export_columns <- model_long_export_columns[
  model_long_export_columns %in% names(model_long)
]

data.table::fwrite(
  model_long[
    ,
    ..model_long_export_columns
  ],
  file.path(
    OUT_DIR,
    "selected_genes_observed_efficacy_vs_three_model_predicted_efficacy_guide_data.tsv.gz"
  ),
  sep = "\t"
)

correlation_export_columns <- setdiff(
  names(correlation_stats),
  "annotation"
)

data.table::fwrite(
  correlation_stats[
    ,
    ..correlation_export_columns
  ],
  file.path(
    OUT_DIR,
    "selected_genes_observed_efficacy_vs_three_model_predicted_efficacy_correlations.tsv"
  ),
  sep = "\t"
)

data.table::fwrite(
  resolved_genes,
  file.path(
    OUT_DIR,
    "selected_genes_resolved.tsv"
  ),
  sep = "\t"
)

cat("\nSelected-gene three-model preview completed.\n")
cat("Output directory:\n  ", OUT_DIR, "\n", sep = "")
cat(
  "Combined figure:\n  selected_genes_observed_efficacy_vs_three_model_predicted_efficacy.png/pdf\n",
  sep = ""
)
cat(
  "Correlation table:\n  selected_genes_observed_efficacy_vs_three_model_predicted_efficacy_correlations.tsv\n",
  sep = ""
)

# Keep useful objects available after Source() in RStudio.
assign(
  "selected_gene_three_model_data",
  model_long,
  envir = .GlobalEnv
)
assign(
  "selected_gene_three_model_correlation_stats",
  correlation_stats,
  envir = .GlobalEnv
)
assign(
  "selected_gene_three_model_plot",
  p_selected,
  envir = .GlobalEnv
)

print(p_selected)


message("Underlying panel G guide data: ", file.path(OUT_DIR,
  "selected_genes_observed_efficacy_vs_three_model_predicted_efficacy_guide_data.tsv.gz"))
message("Panel G statistics: ", file.path(OUT_DIR,
  "selected_genes_observed_efficacy_vs_three_model_predicted_efficacy_correlations.tsv"))
