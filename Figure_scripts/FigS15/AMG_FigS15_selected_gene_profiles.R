#!/usr/bin/env Rscript
# Plot selected guide-level repression–fitness profiles for FigS15.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/FigS15/AMG_FigS15_selected_gene_profiles.R")
# Inputs:
# - edgeR_data/b*_melted_results.tsv[.gz]: PAO1_ID/gene_name, spacer, type, original parent,
#   y_pred, condition and LFC.adj.
# - guide_efficacy_prediction/this_study/Lib1_all_guides_with_model_scores.tsv.gz: full-match parent
#   sequence scores for point colours.
# Upstream / execution order:
# - Use the supplied guide exports. Rebuilding them requires
#   preprocessing_scripts/02_AMG_metacounter.R.
# - The saved sequence-score table is separate; scoring a new table is described in
#   guide_efficacy_prediction/this_study/README.md.
# Outputs:
# - Figure_scripts/FigS15/out_AMG_FigS15_selected_gene_profiles/ — SVG/PDF/PNG plots.
# - Guide points, bins, metrics and fitted curves remain as df_goi, df_bins, gene_stats and
#   predicted_all in R; they are not exported as TSVs.
# Interpretation / operation:
# - Source the complete file. Edit genes_ex and COMPARISON_IDS_TO_PLOT to select
#   genes/conditions; commented alternatives remain available.
# - Full-match guides have nominal relative efficacy 1. Mismatch y_pred is distinct from the
#   sequence score used for colours.
# - Selected profiles are binned and fitted here. The drc objective check maintains the
#   configured anchor influence across weighting conventions.
#
# R packages used/checked: data.table, dplyr, drc, ggplot2, ggrepel, hrbrthemes, purrr,
# readr, stringr, svglite, tibble.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(hrbrthemes)
  library(drc)
  library(data.table)
})

# -----------------------------
# USER SETTINGS: input paths, genes and comparisons
# -----------------------------

# Standard output folder of preprocessing_scripts/02_AMG_metacounter.R.
# Change this one path if your existing exports are stored elsewhere.
EDGE_DIR <- "edgeR_data"
GUIDE_SCORE_FILE <- "guide_efficacy_prediction/this_study/Lib1_all_guides_with_model_scores.tsv.gz"

# Save both plot versions: parent-family colours and parent efficacy-score colours.
# Set FALSE to display plots in RStudio without writing files.
SAVE_PLOTS <- TRUE
OUTPUT_DIR <- "Figure_scripts/FigS15/out_AMG_FigS15_selected_gene_profiles"
OUTPUT_PREFIX <- "FigS15_selected_gene_profiles"
SAVE_FORMATS <- c("svg", "pdf", "png")
# Export dimensions in inches; adjust these for your intended layout.
SAVE_WIDTH_IN <- 16
SAVE_HEIGHT_IN <- 16
SAVE_DPI <- 300  # PNG resolution; SVG and PDF remain vector graphics.

# Select one or several active comparison IDs from comparison_manifest below.
# Columns from left to right in Figure S15.
COMPARISON_IDS_TO_PLOT <- c(
  "B4_TOB_0.025_T1_vs_None_T1",
  "B4_TOB_0.1_T1_vs_None_T1",
  "B1_TOB_0.4_T1_vs_None_T1",
  "B2_AMK_2_T1_vs_None_T1"
)

# Rows from top to bottom in Figure S15.
# A: atpI; B: secA; C: ftsZ; D: accA, accB, accC, accD.
# Gene names or PAO1 locus IDs are accepted.
genes_ex <- c("atpI", "secA", "ftsZ", "accA", "accB", "accC", "accD")

gene_facet_labeller <- ggplot2::as_labeller(
  stats::setNames(
    paste0(
      "italic('",
      genes_ex,
      "')"
    ),
    genes_ex
  ),
  default = ggplot2::label_parsed
)

# Allow untreated-control comparisons only when their IDs are selected above.
# Figure S15 uses antibiotic-versus-untreated comparisons, so this remains FALSE.
INCLUDE_CONTROL_COMPARISONS <- FALSE

# Fixed treatment ordering used throughout this project:
# TOB 0.025, TOB 0.1, TOB 0.4, TOB 1,
# TOB T2, AMK, then FBS.
comparison_manifest <- tibble::tribble(
  ~comparison_id,                    ~batch, ~raw_condition,                         ~display_label,          ~comparison_class, ~comparison_order,
  "B4_TOB_0.025_T1_vs_None_T1",      "b4",   "TOB_0.025_T1 - None_0_T1",             "TOB 0.025\nvs untreated T1", "treatment",         1,
  "B4_TOB_0.1_T1_vs_None_T1",        "b4",   "TOB_0.1_T1 - None_0_T1",               "TOB 0.1\nvs untreated T1",   "treatment",         2,
  "B1_TOB_0.4_T1_vs_None_T1",        "b1",   "TOB_0.4_T1 - None_0_T1",               "TOB 0.4\nvs untreated T1",   "treatment",         3,
  "B4_TOB_1_T1_vs_None_T1",          "b4",   "TOB_1_T1 - None_0_T1",                 "TOB 1\nvs untreated T1",     "treatment",         4,
  "B1_TOB_0.4_T2_vs_None_T2",        "b1",   "TOB_0.4_T2 - None_0_T2",               "TOB 0.4 T2\nvs untreated T2", "treatment",         5,
  "B2_AMK_2_T1_vs_None_T1",          "b2",   "AMK_2_T1 - None_0_T1",                 "AMK 2\nvs untreated T1",     "treatment",         6,
  # "B3_FBS_60_T1_vs_None_T1",         "b3",   "FBS_60_T1 - None_0_T1",                "FBS 60",                "treatment",         7,
  #"B3_FBS_60_T1_vs_None_T0",         "b3",   "FBS_60_T1 - None_0_T0",                "FBS 60 vs T0",                "treatment",         7,
  "B1_None_T1_vs_T0",                "b1",   "None_0_T1 - None_0_T0",                "B1 control T1 vs T0",   "control",         101
  # "B1_None_T2_vs_T0",                "b1",   "None_0_T2 - None_0_T0",                "B1 control T2 vs T0",   "control",         102,
  # "B2_None_T1_vs_T0",                "b2",   "None_0_T1 - None_0_T0",                "B2 control T1 vs T0",   "control",         103,
  # "B3_None_T1_vs_T0",                "b3",   "None_0_T1 - None_0_T0",                "B3 control T1 vs T0",   "control",         104,
  # "B4_None_T1_vs_T0",                "b4",   "None_0_T1 - None_0_T0",                "B4 control T1 vs T0",   "control",         105
)

comparison_manifest_selected <- comparison_manifest %>%
  dplyr::filter(
    comparison_class == "treatment" |
      isTRUE(INCLUDE_CONTROL_COMPARISONS)
  ) %>%
  dplyr::filter(comparison_id %in% COMPARISON_IDS_TO_PLOT) %>%
  dplyr::arrange(comparison_order)

missing_requested <- setdiff(COMPARISON_IDS_TO_PLOT, comparison_manifest_selected$comparison_id)
if (length(missing_requested)) {
  stop("Enable/uncomment these manifest rows (and control comparisons when needed): ",
    paste(missing_requested, collapse = ", "))
}

conds_ex <- comparison_manifest_selected$display_label

find_melted_file <- function(batch_id) {
  candidates <- c(
    file.path(EDGE_DIR, paste0(batch_id, "_melted_results.tsv.gz")),
    file.path(EDGE_DIR, batch_id, paste0(batch_id, "_melted_results.tsv.gz")),
    file.path(EDGE_DIR, batch_id, "melted_results.tsv.gz")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop(
      "No melted_results file found for ", batch_id,
      ". Checked:\n", paste(candidates, collapse = "\n")
    )
  }
  existing[1]
}

melted_paths <- stats::setNames(
  vapply(unique(comparison_manifest_selected$batch), find_melted_file, character(1)),
  unique(comparison_manifest_selected$batch)
)

# No merged control condition is used.
# Downstream code checks this variable when handling combined profiles.
merged_condition_name <- "__not_used_for_manuscript_subset__"

# -----------------------------
# Analysis settings
# -----------------------------

do_sigmoid <- TRUE   # TRUE = slow; needed for EC50

# Intended anchor strength: legacy residual multiplier 3, giving nine times
# one guide's influence in the least-squares objective. See the drc check below.
anchor_residual_multiplier <- 3

# EC50 quality control
emax_threshold <- 1  # require |Emax| > this value

# binning
# must match loop
y_floor   <- 0.05
y_cap     <- 1.00
bin_width <- 0.05
bin_breaks <- unique(c(seq(y_floor, y_cap, by = bin_width), y_cap))
bin_mids   <- (head(bin_breaks, -1) + tail(bin_breaks, -1)) / 2

# thresholds in LFC.adj space
cutoff_sens <- -2
cutoff_prot <-  2

# --- SENS: exact loop parameters ---
recover_thr_sens <- -1.5
K_sens           <- 2
R_sens           <- 2
min_bin_n_sens   <- 3
min_lastbin_n    <- 8
strong_lastbin   <- -2

# --- PROT: exact loop parameters ---
K_prot         <- 2
min_bin_n_prot <- 3

# Less strict settings used in the recent versions
cutoff_sens <- -1
recover_thr_sens <- -1
strong_lastbin <- -1
emax_threshold <- 1

# Kept only for compatibility with old merged-condition logic.
merged_cutoff_sens <- -1
merged_recover_thr_sens <- -1
merged_strong_lastbin <- -1
merged_emax_threshold <- 1

# -----------------------------
# LOAD
# -----------------------------

melted_results <- purrr::imap_dfr(
  melted_paths,
  function(path, batch_id) {
    requested <- comparison_manifest_selected %>%
      dplyr::filter(batch == batch_id)
    
    dat <- readr::read_tsv(path, show_col_types = FALSE) %>%
      dplyr::mutate(
        batch = batch_id,
        raw_condition = as.character(condition)
      ) %>%
      dplyr::inner_join(
        requested %>%
          dplyr::select(
            batch,
            raw_condition,
            comparison_id,
            display_label,
            comparison_class,
            comparison_order
          ),
        by = c("batch", "raw_condition")
      ) %>%
      dplyr::mutate(
        condition = display_label
      )
    
    dat
  }
)

req_cols <- c(
  "PAO1_ID", "gene_name", "type", "y_pred",
  "LFC", "LFC.adj", "condition", "batch",
  "comparison_id", "raw_condition"
)
miss <- setdiff(req_cols, names(melted_results))
if (length(miss) > 0) {
  stop(
    "melted_results is missing: ",
    paste(miss, collapse = ", ")
  )
}

missing_conds <- setdiff(
  conds_ex,
  unique(melted_results$condition)
)

if (length(missing_conds) > 0) {
  warning(
    "These requested comparisons are not present and will be dropped:\n",
    paste(missing_conds, collapse = "\n")
  )
  conds_ex <- intersect(
    conds_ex,
    unique(melted_results$condition)
  )
}

if (length(conds_ex) == 0) {
  stop("No requested comparisons are present in the manuscript-subset results.")
}

melted_results <- melted_results %>%
  dplyr::mutate(
    condition = factor(condition, levels = conds_ex)
  )

message("Genes selected: ", paste(genes_ex, collapse = ", "))
message("Conditions selected:")
message(paste(conds_ex, collapse = "\n"))

# -----------------------------
# Load current Lib1 guide scores for dot coloring
# -----------------------------

# Canonical output from the current multibatch Kawalek-model workflow.
# No legacy guide_abundance_*_with_scores_v2.csv fallback is used.
if (!file.exists(GUIDE_SCORE_FILE)) {
  stop(
    "The current Lib1 guide-score file was not found:\n",
    GUIDE_SCORE_FILE,
    "\nRun the finalized multibatch guide-efficacy model script first."
  )
}

message("Using current Lib1 guide-score file: ", GUIDE_SCORE_FILE)
Lib1_guide_scores <- data.table::fread(
  GUIDE_SCORE_FILE,
  header = TRUE,
  showProgress = TRUE
)

required_score_columns <- c("type", "PAO1_ID", "spacer")
missing_score_columns <- setdiff(
  required_score_columns,
  names(Lib1_guide_scores)
)

if (length(missing_score_columns) > 0L) {
  stop(
    "The current Lib1 guide-score file is missing: ",
    paste(missing_score_columns, collapse = ", "),
    "\nAvailable columns: ",
    paste(names(Lib1_guide_scores), collapse = ", ")
  )
}

# The canonical Lib1 output uses predicted_efficacy_score.
# kawalek_efficacy_score is accepted only to remain compatible with
# files produced by the standalone RStudio scorer.
score_column_candidates <- c(
  "predicted_efficacy_score",
  "kawalek_efficacy_score"
)

score_column <- score_column_candidates[
  score_column_candidates %in% names(Lib1_guide_scores)
][1]

if (is.na(score_column)) {
  stop(
    "No Kawalek efficacy-score column was found in:\n",
    GUIDE_SCORE_FILE,
    "\nExpected one of: ",
    paste(score_column_candidates, collapse = ", ")
  )
}

status_column_candidates <- c(
  "score_status",
  "kawalek_score_status"
)

status_column <- status_column_candidates[
  status_column_candidates %in% names(Lib1_guide_scores)
][1]

Lib1_guide_scores <- Lib1_guide_scores %>%
  dplyr::mutate(
    type = tolower(trimws(as.character(type))),
    PAO1_ID = trimws(as.character(PAO1_ID)),
    spacer = toupper(gsub("\\s+", "", as.character(spacer))),
    kawalek_score__ = suppressWarnings(
      as.numeric(.data[[score_column]])
    ),
    score_status__ = if (is.na(status_column)) {
      "scorable"
    } else {
      tolower(trimws(as.character(.data[[status_column]])))
    }
  )

# Map only experimentally present perfect guides:
# parent spacer -> current Kawalek predicted efficacy.
parent_kawalek_map <- Lib1_guide_scores %>%
  dplyr::filter(
    type == "perfect",
    is.finite(kawalek_score__),
    is.na(score_status__) |
      score_status__ == "" |
      score_status__ %in% c("scorable", "scored", "ok", "success")
  ) %>%
  dplyr::transmute(
    PAO1_ID,
    parent_spacer = spacer,
    parent_Kawalek_score = kawalek_score__
  ) %>%
  dplyr::group_by(PAO1_ID, parent_spacer) %>%
  dplyr::summarise(
    parent_Kawalek_score = dplyr::first(parent_Kawalek_score),
    n_score_rows = dplyr::n(),
    .groups = "drop"
  )

if (nrow(parent_kawalek_map) == 0L) {
  stop(
    "No finite, scorable perfect-guide Kawalek scores were recovered from:\n",
    GUIDE_SCORE_FILE
  )
}

duplicate_score_keys <- parent_kawalek_map %>%
  dplyr::filter(n_score_rows > 1L)

if (nrow(duplicate_score_keys) > 0L) {
  warning(
    nrow(duplicate_score_keys),
    " PAO1_ID + spacer keys occurred more than once in the Lib1 score file. ",
    "The first identical-key score was retained."
  )
}

parent_kawalek_map <- parent_kawalek_map %>%
  dplyr::select(
    PAO1_ID,
    parent_spacer,
    parent_Kawalek_score
  )

message(
  "Perfect parent guides available for Kawalek coloring: ",
  format(nrow(parent_kawalek_map), big.mark = ",")
)

# -----------------------------
# Helpers
# -----------------------------
interp_at_cutoff_safe <- function(x, y, i_cross, cutoff, x_floor = 0.05, x_cap = 1.00) {
  if (is.na(i_cross)) return(NA_real_)
  if (i_cross <= 1)   return(pmin(pmax(x[i_cross], x_floor), x_cap))
  
  i1 <- i_cross - 1
  x1 <- x[i1]; y1 <- y[i1]
  x2 <- x[i_cross]; y2 <- y[i_cross]
  
  if (!is.finite(x1) || !is.finite(x2) || !is.finite(y1) || !is.finite(y2)) return(NA_real_)
  
  if ((y1 <= cutoff && y2 <= cutoff) || (y1 > cutoff && y2 > cutoff) || (y1 == y2)) {
    return(pmin(pmax(x2, x_floor), x_cap))
  }
  
  x_hat <- x1 + (cutoff - y1) * (x2 - x1) / (y2 - y1)
  x_hat <- max(min(x_hat, max(x1, x2)), min(x1, x2))
  x_hat <- pmin(pmax(x_hat, x_floor), x_cap)
  x_hat
}

find_sustained_cross_idx_le <- function(bin_i, yv, nv,
                                        cutoff,
                                        recover_thr,
                                        K = 2,
                                        R = 1,
                                        min_bin_n = 1,
                                        min_lastbin_n = 12,
                                        strong_lastbin = -1.0) {
  n <- length(yv)
  if (n < 2) return(NA_integer_)
  
  ord <- order(bin_i)
  bin_i <- bin_i[ord]; yv <- yv[ord]; nv <- nv[ord]
  
  idx_ok <- which((yv <= cutoff) & (nv >= min_bin_n))
  if (length(idx_ok) == 0) return(NA_integer_)
  
  cand <- vapply(idx_ok, function(i) {
    tail_y <- yv[i:n]
    tail_n <- nv[i:n]
    
    use <- tail_n >= min_bin_n
    tail_y <- tail_y[use]
    if (length(tail_y) == 0) return(FALSE)
    
    sustain_n <- sum(tail_y <= cutoff, na.rm = TRUE)
    recover_n <- sum(tail_y > recover_thr, na.rm = TRUE)
    (sustain_n >= K) && (recover_n <= R)
  }, logical(1))
  
  if (any(cand)) return(idx_ok[min(which(cand))])
  
  i_last <- idx_ok[length(idx_ok)]
  if (i_last == n && is.finite(yv[n]) && is.finite(nv[n]) &&
      nv[n] >= min_lastbin_n && yv[n] <= strong_lastbin) {
    return(n)
  }
  
  NA_integer_
}

find_run_start_adj_ge <- function(bin_i, yv, nv, cutoff, K = 2, min_bin_n = 1) {
  n <- length(yv)
  if (n < K) return(NA_integer_)
  
  ord <- order(bin_i)
  bin_i <- bin_i[ord]; yv <- yv[ord]; nv <- nv[ord]
  
  idx_ok <- which((yv >= cutoff) & (nv >= min_bin_n))
  if (length(idx_ok) == 0) return(NA_integer_)
  
  for (i in idx_ok) {
    target_bins <- bin_i[i] + 0:(K - 1)
    pos <- match(target_bins, bin_i)
    if (any(is.na(pos))) next
    if (all(yv[pos] >= cutoff, na.rm = TRUE) && all(nv[pos] >= min_bin_n, na.rm = TRUE)) return(i)
  }
  NA_integer_
}

interp_at_thr_ge <- function(x, y, i_run_start, thr, x_floor = 0.05, x_cap = 1.00) {
  if (is.na(i_run_start)) return(NA_real_)
  if (i_run_start <= 1) return(pmin(pmax(x[1], x_floor), x_cap))
  
  i1 <- i_run_start - 1
  x1 <- x[i1]; y1 <- y[i1]
  x2 <- x[i_run_start]; y2 <- y[i_run_start]
  
  if (!is.finite(x1) || !is.finite(x2) || !is.finite(y1) || !is.finite(y2) || y1 == y2) {
    return(pmin(pmax(x2, x_floor), x_cap))
  }
  
  x_hat <- x1 + (thr - y1) * (x2 - x1) / (y2 - y1)
  x_hat <- max(min(x_hat, max(x1, x2)), min(x1, x2))
  x_hat <- pmin(pmax(x_hat, x_floor), x_cap)
  x_hat
}

# -----------------------------
# Sigmoid fit helper: returns predicted p(dose) + CI
# -----------------------------

# ---- global fit settings for preview sigmoid overlay ----
# Reproducibility: drc 3.0-1 multiplies continuous-response residuals by weights,
# whereas drc 4.0-0 divides them by weights before squaring. Thus argument 3
# upweights the anchor in the older release but downweights it in the newer one.
# Check the installed objective directly: a unit residual with argument 3 gives
# 9 under multiplication or 1/9 under division. This preserves the intended
# anchor strength across both conventions without changing any real guide weight.
anchor_weight_probe <- getFromNamespace("drmEMls", "drc")(
  dose = 1, resp = 1, multCurves = function(dose, parm) rep(0, length(dose)),
  startVec = 0, robustFct = function(residual) residual^2,
  weights = 3, rmNA = FALSE
)$opfct(0)

if (isTRUE(all.equal(as.numeric(anchor_weight_probe), 9, tolerance = 1e-10))) {
  anchor_w <- anchor_residual_multiplier
} else if (isTRUE(all.equal(as.numeric(anchor_weight_probe), 1/9, tolerance = 1e-10))) {
  anchor_w <- 1 / anchor_residual_multiplier
} else {
  stop("Unrecognised drc weighting convention; check the installed package before fitting.")
}
message("drc ", utils::packageVersion("drc"), ": anchor weights argument = ",
        signif(anchor_w, 6), "; relative least-squares weight = ", anchor_residual_multiplier^2)
dose_grid <- seq(0.1, 1, 0.01) # prediction grid (optional; keep if already defined)

generate_predicted_df_2 <- function(temp) {
  temp$condition2 <- paste(temp$condition, temp$gene_name, sep = "_")
  
  model_list <- by(temp, temp$condition2, function(sub_df) {
    tryCatch({
      # real points weight
      sub_df <- sub_df %>% dplyr::mutate(w = 1)
      
      # anchor at low dose to stabilize (upweight via global anchor_w)
      anchor_point <- tibble::tibble(
        y_pred     = 0.1,
        LFC.adj    = 0,
        w          = anchor_w,
        condition  = unique(sub_df$condition),
        gene_name  = unique(sub_df$gene_name),
        condition2 = unique(sub_df$condition2)
      )
      
      sub_df_aug <- dplyr::bind_rows(sub_df, anchor_point)
      
      drc::drm(
        LFC.adj ~ y_pred,
        data    = sub_df_aug,
        weights = w,
        fct     = drc::LL2.4()
      )
    }, error = function(e) NULL)
  })
  
  model_list <- Filter(Negate(is.null), model_list)
  if (length(model_list) == 0) return(tibble::tibble())
  
  # prediction grid from global dose_grid
  model_data <- data.frame(dose = dose_grid)
  
  predicted_data <- lapply(model_list, function(model) {
    tryCatch({
      pm <- predict(model, newdata = data.frame(y_pred = model_data$dose), interval = "confidence")
      data.frame(dose = model_data$dose, p = pm[, 1], pmin = pm[, 2], pmax = pm[, 3])
    }, error = function(e) NULL)
  })
  
  predicted_data <- Filter(Negate(is.null), predicted_data)
  if (length(predicted_data) == 0) return(tibble::tibble())
  
  tibble::as_tibble(do.call(rbind, predicted_data))
}
# generate_predicted_df_2 <- function(temp) {
#   temp$condition2 <- paste(temp$condition, temp$gene_name, sep = "_")
#   
#   model_list <- by(temp, temp$condition2, function(sub_df) {
#     tryCatch({
#       # anchor at low dose to stabilize
#       anchor_point <- tibble(
#         y_pred = 0.1,
#         LFC.adj = 0,
#         condition = unique(sub_df$condition),
#         gene_name = unique(sub_df$gene_name),
#         condition2 = unique(sub_df$condition2)
#       )
#       sub_df_aug <- bind_rows(sub_df, anchor_point)
#       drm(LFC.adj ~ y_pred, data = sub_df_aug, fct = LL2.4())
#     }, error = function(e) NULL)
#   })
#   
#   model_list <- Filter(Negate(is.null), model_list)
#   if (length(model_list) == 0) return(tibble())
#   
#   model_data <- data.frame(dose = seq(0.1, 1, by = 0.01))
#   
#   predicted_data <- lapply(model_list, function(model) {
#     tryCatch({
#       pm <- predict(model, newdata = model_data, interval = "confidence")
#       data.frame(dose = model_data$dose, p = pm[,1], pmin = pm[,2], pmax = pm[,3])
#     }, error = function(e) NULL)
#   })
#   predicted_data <- Filter(Negate(is.null), predicted_data)
#   if (length(predicted_data) == 0) return(tibble())
#   
#   tibble::as_tibble(do.call(rbind, predicted_data))
# }

# -----------------------------
# EC50 from predicted curve
# - Normalize so last fitted point is 100% effect
# - EC50 is dose at 50% effect (linear interpolation)
# -----------------------------
# -----------------------------
# EC50 from predicted curve (robust)
# - Uses pmin=min(p), pmax=max(p)  (NOT first/last point)
# - EC50 = dose where curve crosses halfway between (pmin,pmax)
# - Computes only if |Emax| > threshold
# -----------------------------
# -----------------------------
# EC50 from predicted curve (endpoint-defined)
# - E0   = predicted p at lowest dose in pred_df (usually 0.1)
# - Emax = predicted p at highest dose in pred_df (usually 1.0)
# - EC50 = dose where p crosses halfway between (E0, Emax)
# - Only compute EC50 if |(Emax - E0)| >= emax_thr  (dynamic range threshold)
# -----------------------------
calc_ec50_from_pred <- function(pred_df, emax_thr = 1) {
  
  if (is.null(pred_df) || nrow(pred_df) < 5)
    return(tibble(EC50 = NA_real_, Emax = NA_real_, E0 = 0, delta = NA_real_))
  
  pred_df <- pred_df %>% arrange(dose)
  
  p <- pred_df$p
  d <- pred_df$dose
  
  if (any(!is.finite(p)) || any(!is.finite(d)))
    return(tibble(EC50 = NA_real_, Emax = NA_real_, E0 = 0, delta = NA_real_))
  
  # Endpoint-only, baseline=0
  E0   <- 0
  Emax <- p[length(p)]
  delta <- Emax
  
  if (!is.finite(delta) || abs(delta) < emax_thr)
    return(tibble(EC50 = NA_real_, Emax = Emax, E0 = E0, delta = delta))
  
  half <- 0.5 * delta
  
  idx <- which((p[-length(p)] - half) * (p[-1] - half) <= 0)
  if (length(idx) == 0)
    return(tibble(EC50 = NA_real_, Emax = Emax, E0 = E0, delta = delta))
  
  i <- idx[1]
  p1 <- p[i];   p2 <- p[i + 1]
  d1 <- d[i];   d2 <- d[i + 1]
  
  if (!is.finite(p1) || !is.finite(p2) || p1 == p2)
    return(tibble(EC50 = NA_real_, Emax = Emax, E0 = E0, delta = delta))
  
  ec50 <- d1 + (half - p1) * (d2 - d1) / (p2 - p1)
  
  tibble(EC50 = ec50, Emax = Emax, E0 = E0, delta = delta)
}
# -----------------------------
# Build df0 for selected conditions
# -----------------------------
df0 <- melted_results %>%
  filter(condition %in% conds_ex, type %in% c("mismatch","perfect")) %>%
  mutate(
    PAO1_ID   = as.character(PAO1_ID),
    gene_name = gsub("'", "", as.character(gene_name)),
    condition = as.character(condition),
    y_pred = as.numeric(y_pred),
    y_pred = ifelse(type == "perfect", 1, y_pred),
    y_pred = pmin(pmax(y_pred, y_floor), y_cap),
    LFC       = as.numeric(LFC),
    LFC.adj   = as.numeric(LFC.adj),
    original  = if ("original" %in% names(.)) as.character(original) else "NA_parent"
  ) %>%
  filter(!is.na(PAO1_ID), !is.na(gene_name), !is.na(y_pred), !is.na(LFC.adj))

df_goi <- df0 %>%
  filter(gene_name %in% genes_ex | PAO1_ID %in% genes_ex) %>%
  mutate(
    gene_ex = dplyr::case_when(
      gene_name %in% genes_ex ~ gene_name,
      PAO1_ID   %in% genes_ex ~ PAO1_ID,
      TRUE ~ gene_name
    ),
    gene_ex = factor(gene_ex, levels = genes_ex),
    condition = factor(condition, levels = conds_ex)
  )
# Define parent spacer for each guide point:
# - mismatch: parent = original
# - perfect:  parent = itself (if spacer exists), otherwise original (if already set upstream)
if ("spacer" %in% names(df_goi)) {
  df_goi <- df_goi %>%
    mutate(parent_spacer = ifelse(type == "mismatch", original, spacer))
} else {
  # fallback: assume 'original' already stores parent for both mismatch & perfect
  df_goi <- df_goi %>%
    mutate(parent_spacer = original)
}

df_goi <- df_goi %>%
  mutate(
    PAO1_ID = trimws(as.character(PAO1_ID)),
    parent_spacer = toupper(
      gsub("\\s+", "", as.character(parent_spacer))
    ),
    parent_spacer = dplyr::na_if(parent_spacer, "NA_PARENT")
  ) %>%
  left_join(
    parent_kawalek_map,
    by = c("PAO1_ID", "parent_spacer")
  )

kawalek_color_coverage <- df_goi %>%
  summarise(
    n_points = dplyr::n(),
    n_points_with_score = sum(
      is.finite(parent_Kawalek_score)
    ),
    score_coverage_pct = 100 * mean(
      is.finite(parent_Kawalek_score)
    ),
    n_unique_parent_guides = dplyr::n_distinct(
      paste(PAO1_ID, parent_spacer)
    ),
    n_unique_parent_guides_with_score = dplyr::n_distinct(
      paste(
        PAO1_ID[is.finite(parent_Kawalek_score)],
        parent_spacer[is.finite(parent_Kawalek_score)]
      )
    )
  )

message(
  "Kawalek dot-color coverage: ",
  kawalek_color_coverage$n_points_with_score,
  " / ",
  kawalek_color_coverage$n_points,
  " points (",
  sprintf("%.1f", kawalek_color_coverage$score_coverage_pct),
  "%)"
)

if (kawalek_color_coverage$n_points_with_score == 0L) {
  stop(
    "None of the plotted parent guides matched the current Lib1 score file. ",
    "Check PAO1_ID and spacer conventions."
  )
}

# per (gene_ex, condition): do we have ANY mismatch guides?
has_mismatch_tbl <- df_goi %>%
  group_by(gene_ex, condition) %>%
  summarise(has_mismatch = any(type == "mismatch"), .groups = "drop")

if (nrow(df_goi) == 0) stop("No rows after filtering: check genes_ex / conds_ex vs your melted_results.")

df_bins <- df_goi %>%
  mutate(
    bin   = cut(y_pred, breaks = bin_breaks, include.lowest = TRUE, right = TRUE),
    bin_i = as.integer(bin),
    y_bin_mid = bin_mids[bin_i]
  ) %>%
  filter(!is.na(bin_i), is.finite(y_bin_mid)) %>%
  group_by(gene_ex, condition, bin_i, y_bin_mid) %>%
  summarise(med_LFC_bin = median(LFC.adj, na.rm = TRUE),
            n_bin = n(),
            .groups = "drop") %>%
  arrange(gene_ex, condition, bin_i)

# -----------------------------
# Optional sigmoid fits per (gene_ex, condition) + EC50
# -----------------------------
predicted_all <- NULL
ec50_tbl <- NULL

if (isTRUE(do_sigmoid)) {
  
  # Only keep groups that have mismatch data
  df_fit <- df_goi %>%
    left_join(has_mismatch_tbl, by = c("gene_ex","condition")) %>%
    filter(has_mismatch)
  
  predicted_all <- df_fit %>%
    group_by(gene_ex, condition) %>%
    group_modify(~{
      out <- tryCatch(generate_predicted_df_2(.x), error = function(e) NULL)
      if (is.null(out) || nrow(out) == 0) return(tibble())
      tibble::as_tibble(out)
    }) %>%
    ungroup()
  
  ec50_tbl <- predicted_all %>%
    group_by(gene_ex, condition) %>%
    group_modify(~ calc_ec50_from_pred(.x, emax_threshold)) %>%
    ungroup()
  #filtering ec lower than available mismatch
  ec50_tbl <- ec50_tbl %>%
    left_join(
      df_fit %>%
        filter(
          type == "mismatch",
          is.finite(y_pred)
        ) %>%
        group_by(gene_ex, condition) %>%
        summarise(
          min_mismatch_y_pred = min(y_pred),
          .groups = "drop"
        ),
      by = c("gene_ex", "condition")
    ) %>%
    mutate(
      EC50 = ifelse(
        is.finite(EC50) &
          is.finite(min_mismatch_y_pred) &
          EC50 < min_mismatch_y_pred,
        NA_real_,
        EC50
      )
    )
  
} else {
  ec50_tbl <- df_goi %>%
    distinct(gene_ex, condition) %>%
    mutate(EC50 = NA_real_, Emax = NA_real_, pmin = NA_real_, pmax = NA_real_)
}

# -----------------------------
# Loop-based sens/prot crossings (still computed if you want vlines)
# -----------------------------
gene_stats <- df_goi %>%
  group_by(gene_ex, condition) %>%
  group_modify(function(dat, key) {
    
    # ---- severity (same as before) ----
    sev_val <- dat %>%
      filter(type == "perfect") %>%
      { if ("LFC" %in% names(.)) filter(., !is.na(LFC), LFC != 0) else . } %>%
      summarise(sev = median(LFC.adj, na.rm = TRUE)) %>%
      pull(sev)
    
    #   sev_val <- dat %>%
    # filter(
    #   type == "perfect",
    #   is.finite(LFC.adj)
    # ) %>%
    # summarise(
    #   sev = median(LFC.adj, na.rm = TRUE)
    # ) %>%
    # pull(sev)
    
    bins_g <- df_bins %>% filter(gene_ex == key$gene_ex, condition == key$condition)
    if (nrow(bins_g) < 2) {
      return(tibble(severity = sev_val, y_sens = NA_real_, y_prot = NA_real_))
    }
    
    # ---- SENS (loop rule) ----
    i_sens <- find_sustained_cross_idx_le(
      bin_i = bins_g$bin_i, yv = bins_g$med_LFC_bin, nv = bins_g$n_bin,
      cutoff = cutoff_sens, recover_thr = recover_thr_sens,
      K = K_sens, R = R_sens,
      min_bin_n = min_bin_n_sens,
      min_lastbin_n = min_lastbin_n,
      strong_lastbin = strong_lastbin
    )
    
    y_sens <- interp_at_cutoff_safe(
      x = bins_g$y_bin_mid, y = bins_g$med_LFC_bin, i_cross = i_sens,
      cutoff = cutoff_sens, x_floor = y_floor, x_cap = y_cap
    )
    
    # ---- PROT (loop rule) ----
    i_prot <- find_run_start_adj_ge(
      bin_i = bins_g$bin_i, yv = bins_g$med_LFC_bin, nv = bins_g$n_bin,
      cutoff = cutoff_prot, K = K_prot, min_bin_n = min_bin_n_prot
    )
    
    y_prot <- interp_at_thr_ge(
      x = bins_g$y_bin_mid, y = bins_g$med_LFC_bin, i_run_start = i_prot,
      thr = cutoff_prot, x_floor = y_floor, x_cap = y_cap
    )
    
    tibble(severity = sev_val, y_sens = y_sens, y_prot = y_prot)
  }) %>%
  ungroup() %>%
  left_join(ec50_tbl, by = c("gene_ex","condition")) %>%
  mutate(
    lab_txt = paste0(
      "Sev=",  ifelse(is.finite(severity), sprintf("%.2f", severity), "NA"),
      "\nEmax=", ifelse(is.finite(Emax), sprintf("%.2f", Emax), "NA"),
      "\nEK50=", ifelse(is.finite(EC50), sprintf("%.2f", EC50), "NA"),
      "\nSens=", ifelse(is.finite(y_sens), sprintf("%.2f", y_sens), "NA"),
      "\nProt=", ifelse(is.finite(y_prot), sprintf("%.2f", y_prot), "NA")
    )
  )

lab_df <- gene_stats %>%
  transmute(gene_ex, condition, x = -Inf, y = -Inf, lab_txt)

### dynamic limits 
# dynamic y-limits: only cap if data exceed bounds
ymin_data <- min(df_goi$LFC.adj, df_bins$med_LFC_bin, na.rm = TRUE)
ymax_data <- max(df_goi$LFC.adj, df_bins$med_LFC_bin, na.rm = TRUE)

ymin_lim <- if (is.finite(ymin_data) && ymin_data < -14) -14 else ymin_data
ymax_lim <- if (is.finite(ymax_data) && ymax_data > 5)   5   else ymax_data


# -----------------------------
# Plot
# -----------------------------
p <- ggplot(df_goi, aes(x = y_pred, y = LFC.adj)) +
  geom_point(aes(color = original), alpha = 0.47, size = 1.5) +
  
  # (optional) keep thresholds
  geom_hline(yintercept = cutoff_sens, linewidth = 0.5, linetype = "dashed", color = "grey30") +
  geom_hline(yintercept = cutoff_prot, linewidth = 0.5, linetype = "dashed", color = "grey60") +
  # severity line per facet (dotdash)
  geom_hline(
    data = gene_stats %>% filter(is.finite(severity)),
    aes(yintercept = severity),
    linewidth = 0.6, linetype = "dotdash", color = "black",
    inherit.aes = FALSE
  ) +
  
  # binned medians
  geom_line(
    data = df_bins,
    aes(x = y_bin_mid, y = med_LFC_bin, group = interaction(gene_ex, condition)),
    linewidth = 0.95, color = "black",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = df_bins,
    aes(x = y_bin_mid, y = med_LFC_bin, size = n_bin),
    shape = 21, fill = "white", color = "black", stroke = 0.7, alpha = 0.95,
    inherit.aes = FALSE
  ) +
  
  # sigmoid overlay
  { if (isTRUE(do_sigmoid) && !is.null(predicted_all) && nrow(predicted_all) > 0)
    geom_ribbon(
      data = predicted_all,
      aes(x = dose, ymin = pmin, ymax = pmax),
      alpha = 0.18,
      inherit.aes = FALSE
    )
  } +
  { if (isTRUE(do_sigmoid) && !is.null(predicted_all) && nrow(predicted_all) > 0)
    geom_line(
      data = predicted_all,
      aes(x = dose, y = p),
      linewidth = 0.85,
      inherit.aes = FALSE
    )
  } +
  
  # EC50 vline (fit-based)
  geom_vline(
    data = gene_stats %>% filter(is.finite(EC50)),
    aes(xintercept = EC50),
    linewidth = 0.7, linetype = "solid", color = "grey10",
    inherit.aes = FALSE
  ) +
  
  # (optional) keep sens/prot vlines
  geom_vline(
    data = gene_stats %>% filter(is.finite(y_sens)),
    aes(xintercept = y_sens),
    linewidth = 0.6, linetype = "dotted", color = "grey20",
    inherit.aes = FALSE
  ) +
  geom_vline(
    data = gene_stats %>% filter(is.finite(y_prot)),
    aes(xintercept = y_prot),
    linewidth = 0.6, linetype = "dotted", color = "grey40",
    inherit.aes = FALSE
  ) +
  
  # label box (bottom-left)
  geom_label(
    data = lab_df,
    aes(x = x, y = y, label = lab_txt),
    hjust = -0.05, vjust = -0.10,
    size = 2.5,
    label.size = 0.1,
    fill = "white",
    inherit.aes = FALSE
  ) +
  
  guides(color = "none") +
  scale_size_continuous(name = "Guides per bin", range = c(1.2, 3.8)) +
  facet_grid(
    gene_ex ~ condition,
    scales = "free_y",
    labeller = ggplot2::labeller(
      gene_ex = gene_facet_labeller
    )
  ) +
  theme_ipsum(base_size = 12, axis_title_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Predicted silencing efficacy",
    y = expression("Adjusted log"[2]*" fold change")
  ) + 
  coord_cartesian(ylim = c(ymin_lim, ymax_lim))

if (interactive()) print(p)



# Increase contrast by winsorizing the current Kawalek parent-guide scores
finite_kawalek_scores <- df_goi$parent_Kawalek_score[
  is.finite(df_goi$parent_Kawalek_score)
]

if (length(finite_kawalek_scores) < 2L) {
  warning(
    "Fewer than two finite Kawalek parent-guide scores are available; ",
    "the color scale will not be winsorized."
  )
  
  df_goi <- df_goi %>%
    mutate(
      parent_Kawalek_score_clip = parent_Kawalek_score
    )
} else {
  kawalek_score_q <- quantile(
    finite_kawalek_scores,
    probs = c(0.05, 0.95),
    na.rm = TRUE,
    names = FALSE
  )
  
  df_goi <- df_goi %>%
    mutate(
      parent_Kawalek_score_clip = pmin(
        pmax(
          parent_Kawalek_score,
          kawalek_score_q[1]
        ),
        kawalek_score_q[2]
      )
    )
}

p_EF <- ggplot(df_goi, aes(x = y_pred, y = LFC.adj)) +
  geom_point(aes(color = parent_Kawalek_score_clip), alpha = 0.47, size = 1.5, na.rm = TRUE) +
  
  geom_hline(yintercept = cutoff_sens, linewidth = 0.5, linetype = "dashed", color = "grey30") +
  geom_hline(yintercept = cutoff_prot, linewidth = 0.5, linetype = "dashed", color = "grey60") +
  
  geom_hline(
    data = gene_stats %>% filter(is.finite(severity)),
    aes(yintercept = severity),
    linewidth = 0.6, linetype = "dotdash", color = "black",
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = df_bins,
    aes(x = y_bin_mid, y = med_LFC_bin, group = interaction(gene_ex, condition)),
    linewidth = 0.95, color = "black",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = df_bins,
    aes(x = y_bin_mid, y = med_LFC_bin, size = n_bin),
    shape = 21, fill = "white", color = "black", stroke = 0.7, alpha = 0.95,
    inherit.aes = FALSE
  ) +
  
  { if (isTRUE(do_sigmoid) && !is.null(predicted_all) && nrow(predicted_all) > 0)
    geom_ribbon(
      data = predicted_all,
      aes(x = dose, ymin = pmin, ymax = pmax),
      alpha = 0.18,
      inherit.aes = FALSE
    )
  } +
  { if (isTRUE(do_sigmoid) && !is.null(predicted_all) && nrow(predicted_all) > 0)
    geom_line(
      data = predicted_all,
      aes(x = dose, y = p),
      linewidth = 0.85,
      inherit.aes = FALSE
    )
  } +
  
  geom_vline(
    data = gene_stats %>% filter(is.finite(EC50)),
    aes(xintercept = EC50),
    linewidth = 0.7, linetype = "solid", color = "grey10",
    inherit.aes = FALSE
  ) +
  
  geom_vline(
    data = gene_stats %>% filter(is.finite(y_sens)),
    aes(xintercept = y_sens),
    linewidth = 0.6, linetype = "dotted", color = "grey20",
    inherit.aes = FALSE
  ) +
  geom_vline(
    data = gene_stats %>% filter(is.finite(y_prot)),
    aes(xintercept = y_prot),
    linewidth = 0.6, linetype = "dotted", color = "grey40",
    inherit.aes = FALSE
  ) +
  
  geom_label(
    data = lab_df,
    aes(x = x, y = y, label = lab_txt),
    hjust = -0.05, vjust = -0.10,
    size = 2.5,
    label.size = 0.1,
    fill = "white",
    inherit.aes = FALSE
  ) +
  
  scale_color_gradient2(
    low = "#053061",      # deeper blue
    mid = "grey90",       # lighter mid
    high = "#67001F",     # deeper red
    midpoint = 0,
    na.value = "grey85",
    name = "Kawalek predicted\nefficacy"
  )+
  
  scale_size_continuous(name = "Guides per bin", range = c(1.2, 3.8)) +
  facet_grid(
    gene_ex ~ condition,
    scales = "free_y",
    labeller = ggplot2::labeller(
      gene_ex = gene_facet_labeller
    )
  ) +
  theme_ipsum(base_size = 12, axis_title_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Predicted silencing efficacy",
    y = expression("Adjusted log"[2]*" fold change")
  ) +
  coord_cartesian(ylim = c(ymin_lim, ymax_lim))

if (interactive()) print(p_EF)


# Save both colour versions of the plots when requested. Plotting and fitting
# use the settings above; the following block exports the plotted objects.
if (isTRUE(SAVE_PLOTS)) {
  formats <- unique(tolower(SAVE_FORMATS))
  if (!length(formats) || anyNA(formats) || any(!formats %in% c("svg", "pdf", "png"))) {
    stop("SAVE_FORMATS must contain svg, pdf and/or png.")
  }
  dimensions <- c(SAVE_WIDTH_IN, SAVE_HEIGHT_IN, SAVE_DPI)
  if (length(dimensions) != 3L || any(!is.finite(dimensions)) || any(dimensions <= 0)) {
    stop("Save width, height and resolution must be positive numbers.")
  }
  if ("svg" %in% formats && !requireNamespace("svglite", quietly = TRUE)) {
    stop("Install svglite to save SVG plots: install.packages('svglite').")
  }
  if ("pdf" %in% formats && !capabilities("cairo")) {
    stop("PDF export requires Cairo support in this R installation.")
  }
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(OUTPUT_DIR)) stop("Could not create output folder: ", OUTPUT_DIR)
  plots_to_save <- list(parent_families = p, parent_efficacy = p_EF)
  saved_plot_files <- character()
  for (plot_name in names(plots_to_save)) {
    for (format in formats) {
      output_file <- file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_", plot_name, ".", format))
      device <- switch(format, svg = svglite::svglite, pdf = grDevices::cairo_pdf, png = "png")
      ggplot2::ggsave(filename = output_file, plot = plots_to_save[[plot_name]], device = device,
                      width = SAVE_WIDTH_IN, height = SAVE_HEIGHT_IN, units = "in", dpi = SAVE_DPI, bg = "white")
      saved_plot_files <- c(saved_plot_files, output_file)
    }
  }
  message("Saved plots to: ", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = TRUE),
          "\n  ", paste(basename(saved_plot_files), collapse = "\n  "))
} else {
  message("Plot saving is disabled (SAVE_PLOTS <- FALSE). Plot objects: p and p_EF.")
}
