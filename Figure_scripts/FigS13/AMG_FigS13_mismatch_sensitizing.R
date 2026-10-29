#!/usr/bin/env Rscript
# Figure S13: additional sensitizing interactions supported by mismatch-guide profiles.
# Run from the project root AFTER preprocessing step 03 with guide export enabled.
# This script performs its own binning, threshold selection, LOESS filtering and
# profile fitting. Step 17 does not need to be run first and is never overwritten.
# Protective candidate selection and protective-candidate plots are omitted.
# Prot is retained only as a diagnostic threshold on the sensitizing profiles.
# Flow: read guides/calls -> bin effects -> find thresholds -> select new candidates
# -> check LOESS support -> fit/display candidate profiles -> save source tables.
#
# Inputs: task-03 b*_AG_guides.tsv.gz, AG_gene_results.tsv.gz,
# AG_gene_recurrence.tsv, Comparison_manifest.tsv and Analysis_settings.tsv.
# guide_efficacy_prediction/this_study/Lib1_all_guides_with_model_scores.tsv.gz is optional
# point-colour input; it does not select candidates or set their efficacy axis.
# Run: source("Figure_scripts/FigS13/AMG_FigS13_mismatch_sensitizing.R")
# Outputs: PNG/PDF/SVG and profile/selection/LOESS/fit/settings TSVs in
# Figure_scripts/FigS13/out_AMG_FigS13_mismatch_sensitizing/.
#
# R packages used/checked: data.table, dplyr, drc, ggplot2, ggrepel, svglite, tibble.

# 1. Files and comparisons. Five T1 drug-versus-matched-control comparisons only.
INPUT_DIR <- "preprocessing_scripts/out_03_AMG_relevant_results"
GUIDE_SCORE_FILE <- "guide_efficacy_prediction/this_study/Lib1_all_guides_with_model_scores.tsv.gz"
OUTPUT_DIR <- "Figure_scripts/FigS13/out_AMG_FigS13_mismatch_sensitizing"
COMPARISON_LABELS <- c(B4_TOB_0.025_T1_vs_None_T1 = "TOB 0.025",
                      B4_TOB_0.1_T1_vs_None_T1 = "TOB 0.1",
                      B1_TOB_0.4_T1_vs_None_T1 = "TOB 0.4",
                      B4_TOB_1_T1_vs_None_T1 = "TOB 1",
                      B2_AMK_2_T1_vs_None_T1 = "AMK")

# Keep the source rule by default: exclude genes recurrent in EITHER direction.
# "sensitizing_only" instead excludes only the recurrent sensitizing set; this
# is a change to candidate selection and may yield a different candidate list.
RECURRENCE_EXCLUSION <- "either_direction"

# 2. Guide efficacy and binning. Full-match guides are assigned efficacy 1.
# This relative mismatch efficacy is distinct from the sequence model used for point colour.
Y_FLOOR <- 0.05
Y_CAP <- 1
BIN_WIDTH <- 0.05
MIN_BIN_N <- 3L
MIN_TOTAL_GUIDES <- 10L
MIN_QUALIFIED_BINS <- 4L
PARTIAL_REPRESSION_MAX <- 0.8

# A sensitizing crossing uses qualified bins at/below -1. Starting from a
# candidate bin, require at least K_SENS such bins and allow at most R_SENS
# qualified bins above RECOVER_THR_SENS. A qualifying final bin can also provide a fallback.
# These fixed +/-1 profile thresholds are NOT the adaptive full-match MAD cutoffs.
CUTOFF_SENS <- -1
RECOVER_THR_SENS <- -1
K_SENS <- 2L
R_SENS <- 2L
MIN_LASTBIN_N <- 8L
STRONG_LASTBIN <- -1
CUTOFF_PROT <- 1
K_PROT <- 2L

# 3. LOESS controls. Changing these settings recomputes candidate support here.
# Fit qualified bin medians, weighted by the number of guides in each bin.
# LOESS needs five qualified bins by default, in addition to the four-bin QC above.
LOESS_SPAN <- 0.75
LOESS_DEGREE <- 2L
LOESS_FAMILY <- "gaussian"
LOESS_GRID_N <- 401L
LOESS_MIN_BINS <- 5L
LOESS_AUC_X_MIN <- 0.10
LOESS_AUC_X_MAX <- PARTIAL_REPRESSION_MAX
USE_LOESS_AUC_FILTER <- TRUE
LOESS_AUC_MIN_SENS <- 0.05
USE_LOESS_WIDTH_FILTER <- TRUE
LOESS_MIN_WIDTH_SENS <- 0.2

# AUC measures fitted area BELOW CUTOFF_SENS; Width measures the total efficacy
# width below it. Both are restricted to the observed-bin range within [0.1,0.8].
# Width may include separate segments; it is not necessarily one continuous run.

# 4. Separate sigmoid fit for Emax/EK50 in panel B. It does not select candidates.
# The anchor is at predicted efficacy 0.1 and effect 0. The compatibility check
# below maintains the configured anchor strength across drc weight conventions.
DO_SIGMOID <- TRUE
ANCHOR_EFFICACY <- 0.1
ANCHOR_RESIDUAL_MULTIPLIER <- 3
EMAX_THRESHOLD <- 1
DOSE_GRID <- seq(0.1, 1, by = 0.01)

# 5. Plot settings. Only sensitizing candidates receive labels and profile panels.
MAX_NEW_LABELS_PER_PANEL <- 1000L
PROFILE_PANELS_PER_PAGE <- 25L
PROFILE_NCOL <- 5L
SHOW_LOESS_IN_PROFILES <- FALSE
SAVE_LOESS_DIAGNOSTIC <- TRUE
# Display A followed by all B pages in RStudio; use the Plots pane arrows to browse.
PRINT_IN_RSTUDIO <- TRUE
PRINT_PROFILE_PAGES_IN_RSTUDIO <- TRUE
FONT_FAMILY <- "Arial"
DPI <- 300
SCATTER_WIDTH_IN <- 17
SCATTER_HEIGHT_IN <- 4.8
PROFILE_WIDTH_IN <- 16
PROFILE_ROW_HEIGHT_IN <- 3.1
PROFILE_Y_LIMITS <- c(-14, 5)
SCORE_CLIP_QUANTILES <- c(0.05, 0.95)
SEED <- 1L
OVERWRITE <- TRUE

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

# Read only requested columns when loading a large guide-score table.
read_table <- function(path, columns = NULL) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  read_once <- function(nrows = -1L, classes = NA) {
    con <- if (endsWith(path, ".gz")) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con))
    utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = "",
                    na.strings = c("", "NA"), nrows = nrows, colClasses = classes)
  }
  if (is.null(columns)) return(data.table::as.data.table(read_once()))
  fields <- names(read_once(0L))
  missing <- setdiff(columns, fields)
  if (length(missing)) stop(path, " lacks: ", paste(missing, collapse = ", "))
  data.table::as.data.table(read_once(classes = ifelse(fields %in% columns, "character", "NULL")))
}

check_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(label, " lacks: ", paste(missing, collapse = ", "))
}

flag <- function(x, label) {
  text <- tolower(trimws(as.character(x)))
  if (anyNA(text) || any(!text %in% c("true", "false", "1", "0"))) stop("Invalid flag: ", label)
  text %in% c("true", "1")
}

unique_key <- function(x, keys, label) {
  if (anyNA(x[, keys, with = FALSE]) || anyDuplicated(x[, keys, with = FALSE])) stop("Missing or repeated identifier in ", label)
}

write_table <- function(x, path) {
  con <- if (endsWith(path, ".gz")) gzfile(path, "wt") else file(path, "wt")
  on.exit(close(con))
  utils::write.table(as.data.frame(x), con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}

save_figure <- function(plot, stem, width, height) {
  ggplot2::ggsave(file.path(OUTPUT_DIR, paste0(stem, ".svg")), plot, device = svglite::svglite, width = width, height = height, bg = "white")
  ggplot2::ggsave(file.path(OUTPUT_DIR, paste0(stem, ".pdf")), plot, device = grDevices::cairo_pdf, width = width, height = height, bg = "white")
  ggplot2::ggsave(file.path(OUTPUT_DIR, paste0(stem, ".png")), plot, width = width, height = height, dpi = DPI, bg = "white")
}

# The next helpers retain the source threshold/interpolation/integration rules.
find_sustained_cross_idx_le <- function(bin_i, yv, nv, cutoff, recover_thr, K = 2L, R = 2L, min_bin_n = 3L, 
  min_lastbin_n = 8L, strong_lastbin = -1) {
  n <- length(yv)
  if (n < 2L) {
    return(NA_integer_)
  }
  ord <- order(bin_i)
  bin_i <- bin_i[ord]
  yv <- yv[ord]
  nv <- nv[ord]
  idx_ok <- which(yv <= cutoff & nv >= min_bin_n)
  if (length(idx_ok) == 0L) {
    return(NA_integer_)
  }
  candidate <- vapply(idx_ok, function(i) {
    tail_y <- yv[i:n]
    tail_n <- nv[i:n]
    keep <- tail_n >= min_bin_n
    tail_y <- tail_y[keep]
    if (length(tail_y) == 0L) {
      return(FALSE)
    }
    sum(tail_y <= cutoff, na.rm = TRUE) >= K && sum(tail_y > recover_thr, na.rm = TRUE) <= 
      R
  }, logical(1))
  if (any(candidate)) {
    return(idx_ok[min(which(candidate))])
  }
  last_index <- idx_ok[length(idx_ok)]
  if (last_index == n && is.finite(yv[n]) && is.finite(nv[n]) && nv[n] >= min_lastbin_n && yv[n] <= 
    strong_lastbin) {
    return(n)
  }
  NA_integer_
}

find_run_start_adj_ge <- function(bin_i, yv, nv, cutoff, K = 2L, min_bin_n = 3L) {
  n <- length(yv)
  if (n < K) {
    return(NA_integer_)
  }
  ord <- order(bin_i)
  bin_i <- bin_i[ord]
  yv <- yv[ord]
  nv <- nv[ord]
  idx_ok <- which(yv >= cutoff & nv >= min_bin_n)
  if (length(idx_ok) == 0L) {
    return(NA_integer_)
  }
  for (i in idx_ok) {
    positions <- match(bin_i[i] + 0:(K - 1L), bin_i)
    if (any(is.na(positions))) {
      next
    }
    if (all(yv[positions] >= cutoff, na.rm = TRUE) && all(nv[positions] >= min_bin_n, na.rm = TRUE)) {
      return(i)
    }
  }
  NA_integer_
}

interpolate_threshold <- function(x, y, index, threshold) {
  if (is.na(index)) {
    return(NA_real_)
  }
  if (index <= 1L) {
    return(pmin(pmax(x[index], Y_FLOOR), Y_CAP))
  }
  i1 <- index - 1L
  x1 <- x[i1]
  y1 <- y[i1]
  x2 <- x[index]
  y2 <- y[index]
  if (!all(is.finite(c(x1, y1, x2, y2))) || y1 == y2) {
    return(pmin(pmax(x2, Y_FLOOR), Y_CAP))
  }
  estimate <- x1 + (threshold - y1) * (x2 - x1)/(y2 - y1)
  estimate <- max(min(estimate, max(x1, x2)), min(x1, x2))
  pmin(pmax(estimate, Y_FLOOR), Y_CAP)
}

trapz_numeric <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 2L) {
    return(NA_real_)
  }
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L))/2, na.rm = TRUE)
}


# Fit and half-effect helpers calculate candidate response curves and their half-effect crossing.
generate_predicted_profile <- function(temp) {
  if (is.null(temp) || nrow(temp) == 0L) {
    return(tibble::tibble())
  }
  temp <- temp %>% dplyr::mutate(w = 1)
  anchor_point <- tibble::tibble(y_pred = ANCHOR_EFFICACY, LFC.adj = 0, w = ANCHOR_W)
  temp_aug <- dplyr::bind_rows(temp %>% dplyr::select(y_pred, LFC.adj, w), anchor_point)
  model <- tryCatch(drc::drm(LFC.adj ~ y_pred, data = temp_aug, weights = w, fct = drc::LL2.4()), 
    error = function(e) {
      NULL
    })
  if (is.null(model)) {
    return(tibble::tibble())
  }
  pred <- tryCatch(predict(model, newdata = data.frame(y_pred = DOSE_GRID), interval = "confidence"), 
    error = function(e) {
      NULL
    })
  if (is.null(pred)) {
    return(tibble::tibble())
  }
  tibble::tibble(dose = DOSE_GRID, p = pred[, 1], pmin = pred[, 2], pmax = pred[, 3])
}

calc_ec50_from_pred <- function(pred_df, emax_thr = 1) {
  if (is.null(pred_df) || nrow(pred_df) < 5L) {
    return(tibble::tibble(EC50 = NA_real_, Emax = NA_real_, E0 = 0, delta = NA_real_))
  }
  pred_df <- pred_df %>% dplyr::arrange(dose)
  p <- pred_df$p
  d <- pred_df$dose
  if (any(!is.finite(p)) || any(!is.finite(d))) {
    return(tibble::tibble(EC50 = NA_real_, Emax = NA_real_, E0 = 0, delta = NA_real_))
  }
  E0 <- 0
  Emax <- p[length(p)]
  delta <- Emax
  if (!is.finite(delta) || abs(delta) < emax_thr) {
    return(tibble::tibble(EC50 = NA_real_, Emax = Emax, E0 = E0, delta = delta))
  }
  half <- 0.5 * delta
  idx <- which((p[-length(p)] - half) * (p[-1] - half) <= 0)
  if (length(idx) == 0L) {
    return(tibble::tibble(EC50 = NA_real_, Emax = Emax, E0 = E0, delta = delta))
  }
  i <- idx[[1L]]
  p1 <- p[i]
  p2 <- p[i + 1L]
  d1 <- d[i]
  d2 <- d[i + 1L]
  if (!is.finite(p1) || !is.finite(p2) || p1 == p2) {
    return(tibble::tibble(EC50 = NA_real_, Emax = Emax, E0 = E0, delta = delta))
  }
  ec50 <- d1 + (half - p1) * (d2 - d1)/(p2 - p1)
  tibble::tibble(EC50 = ec50, Emax = Emax, E0 = E0, delta = delta)
}


# Fit LOESS for one pre-candidate. No extrapolation is used for AUC or Width.
loess_support <- function(bins) {
  q <- bins[n_bin >= MIN_BIN_N]
  empty_curve <- data.table::data.table(y_pred = numeric(), loess_LFC = numeric(), in_auc_domain = logical())
  out <- data.table::data.table(loess_status = "insufficient_bins", loess_auc_sens = NA_real_,
                                loess_width_sens = NA_real_, loess_warning = NA_character_)
  if (nrow(q) < LOESS_MIN_BINS || data.table::uniqueN(q$y_bin_mid) < LOESS_MIN_BINS) return(list(metrics = out, curve = empty_curve))
  messages <- character()
  fit <- withCallingHandlers(tryCatch(stats::loess(med_LFC_bin ~ y_bin_mid, data = as.data.frame(q), weights = n_bin,
                                span = LOESS_SPAN, degree = LOESS_DEGREE, family = LOESS_FAMILY,
                                na.action = stats::na.exclude,
                                control = stats::loess.control(surface = "direct", statistics = "approximate")), error = function(e) NULL),
                              warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (is.null(fit)) { out$loess_status <- "fit_failed"; return(list(metrics = out, curve = empty_curve)) }
  grid <- seq(max(Y_FLOOR, min(q$y_bin_mid)), min(Y_CAP, max(q$y_bin_mid)), length.out = LOESS_GRID_N)
  prediction <- tryCatch(suppressWarnings(stats::predict(fit, newdata = data.frame(y_bin_mid = grid))), error = function(e) rep(NA_real_, length(grid)))
  curve <- data.table::data.table(y_pred = grid, loess_LFC = as.numeric(prediction))
  curve <- curve[is.finite(loess_LFC)]
  curve[, in_auc_domain := y_pred >= LOESS_AUC_X_MIN & y_pred <= LOESS_AUC_X_MAX]
  domain <- curve[in_auc_domain == TRUE]
  if (nrow(domain) < 2L) { out$loess_status <- "no_auc_domain"; return(list(metrics = out, curve = curve)) }
  out[, `:=`(loess_status = "estimated",
              loess_auc_sens = trapz_numeric(domain$y_pred, pmax(0, CUTOFF_SENS - domain$loess_LFC)),
              loess_width_sens = trapz_numeric(domain$y_pred, as.numeric(domain$loess_LFC <= CUTOFF_SENS)),
              loess_warning = if (length(messages)) paste(unique(messages), collapse = "; ") else NA_character_)]
  list(metrics = out, curve = curve)
}

require_packages(c("data.table", "dplyr", "tibble", "ggplot2", "ggrepel", "svglite", if (DO_SIGMOID) "drc"))
suppressPackageStartupMessages(library(dplyr))
if (!RECURRENCE_EXCLUSION %in% c("either_direction", "sensitizing_only")) stop("Invalid RECURRENCE_EXCLUSION.")
integers <- c(MIN_BIN_N, MIN_TOTAL_GUIDES, MIN_QUALIFIED_BINS, K_SENS, MIN_LASTBIN_N, K_PROT,
              LOESS_GRID_N, LOESS_MIN_BINS, PROFILE_PANELS_PER_PAGE, PROFILE_NCOL)
if (any(!is.finite(integers)) || any(integers < 1 | integers != as.integer(integers)) ||
    R_SENS < 0 || R_SENS != as.integer(R_SENS) || !LOESS_DEGREE %in% 0:2 ||
    Y_FLOOR <= 0 || Y_CAP <= Y_FLOOR || BIN_WIDTH <= 0 || PARTIAL_REPRESSION_MAX > Y_CAP ||
    LOESS_AUC_X_MIN < Y_FLOOR || LOESS_AUC_X_MAX > Y_CAP || LOESS_AUC_X_MAX <= LOESS_AUC_X_MIN ||
    LOESS_SPAN <= 0 || !LOESS_FAMILY %in% c("gaussian", "symmetric") ||
    CUTOFF_SENS >= 0 || CUTOFF_PROT <= 0 || LOESS_AUC_MIN_SENS < 0 || LOESS_MIN_WIDTH_SENS < 0)
  stop("Invalid guide, threshold, LOESS or panel settings.")
if (length(DOSE_GRID) < 5L || any(!is.finite(DOSE_GRID)) || any(diff(DOSE_GRID) <= 0) || tail(DOSE_GRID, 1L) != 1)
  stop("DOSE_GRID must be increasing and end at efficacy 1.")

# 6. Preserve the intended anchor influence across drc versions. This is the
# current step-17 check: older drc multiplies residuals by its weights argument;
# newer drc divides by it. Probe the installed objective rather than guessing
# from the version number. All real-guide residual multipliers remain 1.
ANCHOR_W <- NA_real_
anchor_weight_probe <- NA_real_
if (DO_SIGMOID) {
  if (!is.finite(ANCHOR_RESIDUAL_MULTIPLIER) || ANCHOR_RESIDUAL_MULTIPLIER <= 0) stop("Invalid anchor strength.")
  anchor_weight_probe <- getFromNamespace("drmEMls", "drc")(
    dose = 1, resp = 1, multCurves = function(dose, parm) rep(0, length(dose)),
    startVec = 0, robustFct = function(residual) residual^2, weights = 3, rmNA = FALSE)$opfct(0)
  if (isTRUE(all.equal(as.numeric(anchor_weight_probe), 9, tolerance = 1e-10))) {
    ANCHOR_W <- ANCHOR_RESIDUAL_MULTIPLIER
  } else if (isTRUE(all.equal(as.numeric(anchor_weight_probe), 1/9, tolerance = 1e-10))) {
    ANCHOR_W <- 1 / ANCHOR_RESIDUAL_MULTIPLIER
  } else stop("Unrecognized drc weighting convention; check the installed package before fitting.")
}

# 7. Read the saved full-match calls, recurrence flags and matched guide effects.
# These calls are not recalculated when the LOESS settings are changed.
paths <- c(genes = file.path(INPUT_DIR, "AG_gene_results.tsv.gz"), recurrence = file.path(INPUT_DIR, "AG_gene_recurrence.tsv"),
           manifest = file.path(INPUT_DIR, "Comparison_manifest.tsv"))
genes <- read_table(paths["genes"])
recurrence <- read_table(paths["recurrence"])
manifest <- read_table(paths["manifest"])
check_columns(genes, c("PAO1_ID", "comparison_id", "gene_name", "medLFC", "p_stouffer", "stouffer_call"), paths["genes"])
check_columns(recurrence, c("PAO1_ID", "recurrent_sensitizing", "recurrent_protective"), paths["recurrence"])
check_columns(manifest, c("comparison_id", "batch", "count_for_recurrence"), paths["manifest"])
unique_key(genes, c("PAO1_ID", "comparison_id"), "gene results")
unique_key(recurrence, "PAO1_ID", "recurrence table")
unique_key(manifest, "comparison_id", "comparison manifest")
if (length(setdiff(names(COMPARISON_LABELS), manifest$comparison_id))) stop("A required T1 comparison is missing.")
manifest <- manifest[match(names(COMPARISON_LABELS), comparison_id)]
if (any(!flag(manifest$count_for_recurrence, "count_for_recurrence"))) stop("This figure uses T1 comparisons only.")
genes <- genes[comparison_id %in% names(COMPARISON_LABELS)]
if (length(setdiff(names(COMPARISON_LABELS), unique(genes$comparison_id)))) stop("Required gene results are missing.")
recurrence[, recurrent_sensitizing := flag(recurrent_sensitizing, "recurrent_sensitizing")]
recurrence[, recurrent_protective := flag(recurrent_protective, "recurrent_protective")]
paths_guides <- stats::setNames(file.path(INPUT_DIR, paste0(unique(manifest$batch), "_AG_guides.tsv.gz")), unique(manifest$batch))
guides <- data.table::rbindlist(lapply(paths_guides, read_table), use.names = TRUE, fill = TRUE)
check_columns(guides, c("PAO1_ID", "comparison_id", "spacer", "type", "y_pred", "LFC", "LFC.adj", "original"), "AG guide files")
guides[, type := tolower(trimws(type))]
guides <- guides[comparison_id %in% names(COMPARISON_LABELS) & !is.na(PAO1_ID) & type %in% c("perfect", "mismatch") &
                 is.finite(y_pred) & is.finite(LFC.adj) & is.finite(LFC) & LFC != 0]
unique_key(guides, c("comparison_id", "spacer"), "guide results")
guides[type == "perfect", y_pred := 1]
guides[, y_pred := pmin(pmax(y_pred, Y_FLOOR), Y_CAP)]

# 8. Summarise guide effects in efficacy bins, including full-match guides in the
# final bin. Retain all bins for display; only sufficiently populated bins support thresholds.
breaks <- unique(c(seq(Y_FLOOR, Y_CAP, by = BIN_WIDTH), Y_CAP))
mids <- (head(breaks, -1L) + tail(breaks, -1L))/2
guides[, bin_i := as.integer(cut(y_pred, breaks = breaks, include.lowest = TRUE, right = TRUE))]
bins <- guides[!is.na(bin_i), .(med_LFC_bin = stats::median(LFC.adj), n_bin = .N,
                              n_perfect_bin = sum(type == "perfect"), n_mismatch_bin = sum(type == "mismatch")),
                by = .(PAO1_ID, comparison_id, bin_i)]
bins[, y_bin_mid := mids[bin_i]]
data.table::setorder(bins, PAO1_ID, comparison_id, bin_i)
threshold_metrics <- bins[, {
  si <- find_sustained_cross_idx_le(bin_i, med_LFC_bin, n_bin, CUTOFF_SENS, RECOVER_THR_SENS,
                                    K_SENS, R_SENS, MIN_BIN_N, MIN_LASTBIN_N, STRONG_LASTBIN)
  pi <- find_run_start_adj_ge(bin_i, med_LFC_bin, n_bin, CUTOFF_PROT, K_PROT, MIN_BIN_N)
  list(sensitizing_threshold = interpolate_threshold(y_bin_mid, med_LFC_bin, si, CUTOFF_SENS),
       protective_threshold = interpolate_threshold(y_bin_mid, med_LFC_bin, pi, CUTOFF_PROT),
       n_guides_in_bins = sum(n_bin), n_qualified_bins = sum(n_bin >= MIN_BIN_N),
       threshold_qc = sum(n_bin) >= MIN_TOTAL_GUIDES && sum(n_bin >= MIN_BIN_N) >= MIN_QUALIFIED_BINS)
}, by = .(PAO1_ID, comparison_id)]
metrics <- merge(genes[, .(PAO1_ID, comparison_id, gene_name, medLFC, p_stouffer, stouffer_call)], threshold_metrics,
                  by = c("PAO1_ID", "comparison_id"), all.x = TRUE, sort = FALSE)
ri <- match(metrics$PAO1_ID, recurrence$PAO1_ID)
if (anyNA(ri)) stop("Gene results and the recurrence table have different gene coverage.")
metrics[, `:=`(recurrent_sensitizing = recurrence$recurrent_sensitizing[ri], recurrent_protective = recurrence$recurrent_protective[ri])]
metrics[, recurrent_any := recurrent_sensitizing | recurrent_protective]
metrics[, excluded_recurrent := if (RECURRENCE_EXCLUSION == "either_direction") recurrent_any else recurrent_sensitizing]
metrics[, `:=`(sensitizing_supported = !is.na(threshold_qc) & threshold_qc & is.finite(sensitizing_threshold),
               current_sensitizing_hit = stouffer_call == "Sensitizing")]
metrics[, pre_candidate := sensitizing_supported & sensitizing_threshold <= PARTIAL_REPRESSION_MAX &
                           !is.na(current_sensitizing_hit) & !current_sensitizing_hit & !excluded_recurrent]
metrics[is.na(pre_candidate), pre_candidate := FALSE]
metrics[, `:=`(loess_status = "not_a_pre_candidate", loess_auc_sens = NA_real_, loess_width_sens = NA_real_, loess_warning = NA_character_)]

# 9. Apply LOESS only to profiles that could become new sensitizing candidates.
# Skipping already-excluded profiles saves time without changing any selection.
keys <- metrics[pre_candidate == TRUE, .(PAO1_ID, comparison_id)]
loess_curves_list <- vector("list", nrow(keys))
for (i in seq_len(nrow(keys))) {
  id <- keys$PAO1_ID[i]; comparison <- keys$comparison_id[i]
  result <- loess_support(bins[PAO1_ID == id & comparison_id == comparison])
  index <- which(metrics$PAO1_ID == id & metrics$comparison_id == comparison)
  for (column in names(result$metrics)) data.table::set(metrics, index, column, result$metrics[[column]])
  if (nrow(result$curve)) loess_curves_list[[i]] <- data.table::data.table(PAO1_ID = id, comparison_id = comparison, result$curve)
}
loess_curves <- if (any(vapply(loess_curves_list, function(x) !is.null(x), logical(1))))
  data.table::rbindlist(loess_curves_list) else data.table::data.table(PAO1_ID = character(), comparison_id = character(),
                                                                  y_pred = numeric(), loess_LFC = numeric(), in_auc_domain = logical())
metrics[, `:=`(auc_pass = !USE_LOESS_AUC_FILTER | (is.finite(loess_auc_sens) & loess_auc_sens >= LOESS_AUC_MIN_SENS),
               width_pass = !USE_LOESS_WIDTH_FILTER | (is.finite(loess_width_sens) & loess_width_sens >= LOESS_MIN_WIDTH_SENS))]
metrics[, new_sensitizing_candidate := pre_candidate & auc_pass & width_pass]
metrics[, `:=`(display_label = factor(unname(COMPARISON_LABELS[comparison_id]), levels = unname(COMPARISON_LABELS)),
               gene_label = ifelse(is.na(gene_name) | !nzchar(trimws(gene_name)), PAO1_ID, gsub("'", "", gene_name, fixed = TRUE)))]
candidates <- data.table::copy(metrics[new_sensitizing_candidate == TRUE])
candidates[, condition_order := match(comparison_id, names(COMPARISON_LABELS))]
data.table::setorder(candidates, condition_order, sensitizing_threshold, PAO1_ID)
candidates[, panel_label := paste(gene_label, as.character(display_label), sep = "\n")]
candidate_keys <- candidates[, .(PAO1_ID, comparison_id)]
profile_guides <- merge(guides, candidate_keys, by = c("PAO1_ID", "comparison_id"), all = FALSE, sort = FALSE)
profile_bins <- merge(bins, candidate_keys, by = c("PAO1_ID", "comparison_id"), all = FALSE, sort = FALSE)

# 10. Attach parent-guide sequence scores for POINT COLOUR only. This optional
# file never changes mismatch efficacy, binning, LOESS support or candidate selection.
profile_guides[, parent_spacer := ifelse(type == "mismatch" & !is.na(original) & nzchar(original), original, spacer)]
profile_guides[, parent_spacer := toupper(gsub("\\s+", "", parent_spacer))]
profile_guides[, parent_score := NA_real_]
if (file.exists(GUIDE_SCORE_FILE) && nrow(profile_guides)) {
  scores <- read_table(GUIDE_SCORE_FILE, c("PAO1_ID", "spacer", "type", "predicted_efficacy_score", "score_status"))
  scores[, `:=`(type = tolower(trimws(type)), score_status = tolower(trimws(score_status)),
                 predicted_efficacy_score = suppressWarnings(as.numeric(predicted_efficacy_score)),
                 spacer = toupper(gsub("\\s+", "", spacer)))]
  scores <- scores[type == "perfect" & is.finite(predicted_efficacy_score) &
                   (is.na(score_status) | score_status %in% c("", "scorable", "scored", "ok", "success"))]
  score_map <- unique(scores[, .(PAO1_ID, parent_spacer = spacer, parent_score = predicted_efficacy_score)])
  unique_key(score_map, c("PAO1_ID", "parent_spacer"), "parent-guide scores")
  profile_guides[, parent_score := NULL]
  profile_guides <- merge(profile_guides, score_map, by = c("PAO1_ID", "parent_spacer"), all.x = TRUE, sort = FALSE)
} else if (nrow(profile_guides)) message("Guide-score file absent: profile points will be grey; candidate selection is unchanged.")
score_values <- profile_guides$parent_score[is.finite(profile_guides$parent_score)]
profile_guides[, parent_score_display := parent_score]
if (length(score_values) >= 2L) {
  limits <- stats::quantile(score_values, probs = SCORE_CLIP_QUANTILES, names = FALSE)
  profile_guides[, parent_score_display := pmin(pmax(parent_score, limits[1L]), limits[2L])]
}

# 11. Fit the selected profiles separately for Emax and EK50. Retain candidates
# even if a sigmoid fit fails; their evidence is the bin/LOESS criteria above.
fit_rows <- curve_rows <- vector("list", nrow(candidates))
for (i in seq_len(nrow(candidates))) {
  id <- candidates$PAO1_ID[i]; comparison <- candidates$comparison_id[i]
  dat <- profile_guides[PAO1_ID == id & comparison_id == comparison]
  severity <- if (any(dat$type == "perfect")) stats::median(dat$LFC.adj[dat$type == "perfect"]) else NA_real_
  warnings <- character()
  predicted <- if (DO_SIGMOID) withCallingHandlers(generate_predicted_profile(as.data.frame(dat)),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }) else tibble::tibble()
  estimated <- calc_ec50_from_pred(predicted, EMAX_THRESHOLD)
  minimum_mismatch <- if (any(dat$type == "mismatch")) min(dat$y_pred[dat$type == "mismatch"]) else NA_real_
  fit_status <- if (!DO_SIGMOID) "not_requested" else if (!nrow(predicted)) "fit_or_prediction_failed" else
    if (is.finite(estimated$EC50)) "estimated" else "no_valid_half_effect_crossing"
  if (is.finite(estimated$EC50) && is.finite(minimum_mismatch) && estimated$EC50 < minimum_mismatch) {
    estimated$EC50 <- NA_real_; fit_status <- "below_mismatch_support"
  }
  fit_rows[[i]] <- data.table::data.table(PAO1_ID = id, comparison_id = comparison, severity_profile = severity,
    EK50 = estimated$EC50, Emax = estimated$Emax, fit_status = fit_status,
    fit_warning = if (length(warnings)) paste(unique(warnings), collapse = "; ") else NA_character_)
  if (nrow(predicted)) curve_rows[[i]] <- data.table::data.table(PAO1_ID = id, comparison_id = comparison, predicted)
}
fits <- if (length(fit_rows)) data.table::rbindlist(fit_rows) else data.table::data.table(PAO1_ID = character(), comparison_id = character(),
  severity_profile = numeric(), EK50 = numeric(), Emax = numeric(), fit_status = character(), fit_warning = character())
sigmoid_curves <- if (any(vapply(curve_rows, function(x) !is.null(x), logical(1)))) data.table::rbindlist(curve_rows) else
  data.table::data.table(PAO1_ID = character(), comparison_id = character(), dose = numeric(), p = numeric(), pmin = numeric(), pmax = numeric())
candidates <- merge(candidates, fits, by = c("PAO1_ID", "comparison_id"), all.x = TRUE, sort = FALSE)
data.table::setorder(candidates, condition_order, sensitizing_threshold, PAO1_ID)
number <- function(x, digits = 2L) ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
candidates[, lab_txt := paste0("Sev=", number(severity_profile), "\nEmax=", number(Emax), "\nEK50=", number(EK50),
                                "\nSens=", number(sensitizing_threshold), "\nProt=", number(protective_threshold),
                                "\nAUC=", number(loess_auc_sens, 3L), "\nWidth=", number(loess_width_sens))]

# 12. Panel A shows sensitizing thresholds, with the other threshold-supported
# genes retained as reference points. Only the NEW sensitizing candidates are labelled.
scatter <- data.table::copy(metrics[sensitizing_supported == TRUE & is.finite(medLFC)])
scatter[, point_class := data.table::fcase(new_sensitizing_candidate, "Additional sensitizing candidate",
                                          recurrent_any, "Recurrent interaction hit", current_sensitizing_hit,
                                          "Condition-specific sensitizing hit", default = "Other genes")]
colours <- c(`Other genes` = "#BDBDBD", `Condition-specific sensitizing hit` = "#1B7837",
             `Recurrent interaction hit` = "#762A83", `Additional sensitizing candidate` = "#2166AC")
scatter_labels <- scatter[new_sensitizing_candidate == TRUE]
data.table::setorder(scatter_labels, display_label, sensitizing_threshold, medLFC)
scatter_labels <- scatter_labels[, head(.SD, MAX_NEW_LABELS_PER_PANEL), by = .(display_label)]
pA <- ggplot2::ggplot() +
  ggplot2::geom_vline(xintercept = -1, linetype = "dashed", linewidth = 0.4, colour = "grey45") +
  ggplot2::geom_hline(yintercept = PARTIAL_REPRESSION_MAX, linetype = "dotted", linewidth = 0.4, colour = "grey45") +
  ggplot2::geom_point(data = scatter[point_class == "Other genes"], ggplot2::aes(medLFC, sensitizing_threshold), colour = colours[["Other genes"]], alpha = 0.42, size = 1.45) +
  ggplot2::geom_point(data = scatter[point_class != "Other genes"], ggplot2::aes(medLFC, sensitizing_threshold, colour = point_class), alpha = 0.92, size = 2.15) +
  ggrepel::geom_text_repel(data = scatter_labels, ggplot2::aes(medLFC, sensitizing_threshold, label = gene_label, colour = point_class),
                         size = 2.4, fontface = "italic", max.overlaps = Inf, box.padding = 0.30, point.padding = 0.10,
                         min.segment.length = 0, seed = SEED, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = colours, breaks = names(colours)[-1L], drop = FALSE, name = "Highlighted genes") +
  ggplot2::scale_y_continuous(limits = c(Y_FLOOR, Y_CAP), breaks = seq(0.2, 1, by = 0.2)) +
  ggplot2::facet_wrap(~display_label, ncol = 5, scales = "free_x", drop = FALSE) +
  ggplot2::theme_classic(base_family = FONT_FAMILY, base_size = 9) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), strip.text = ggplot2::element_text(face = "bold"), legend.position = "bottom") +
  ggplot2::labs(title = paste0("Additional sensitizing interactions: ", nrow(candidates), " gene-condition records"),
                 x = "Gene-level interaction severity (median log2 fold change)", y = "Sensitizing threshold\n(predicted guide efficacy)")

# 13. Build profile pages from exactly the candidates identified in A. Keep full-
# match points at efficacy 1 visible. The black bin medians and sigmoid band/line
# follow the source display; an optional blue LOESS line helps inspect selection.
panel_meta <- candidates[, .(PAO1_ID, comparison_id, panel_label)]
profile_guides <- merge(profile_guides, panel_meta, by = c("PAO1_ID", "comparison_id"), all.x = TRUE, sort = FALSE)
profile_bins <- merge(profile_bins, panel_meta, by = c("PAO1_ID", "comparison_id"), all.x = TRUE, sort = FALSE)
sigmoid_curves <- merge(sigmoid_curves, panel_meta, by = c("PAO1_ID", "comparison_id"), all.x = TRUE, sort = FALSE)
profile_loess <- merge(loess_curves, panel_meta, by = c("PAO1_ID", "comparison_id"), all = FALSE, sort = FALSE)
profile_pages <- list()
if (nrow(candidates)) {
  pages <- split(seq_len(nrow(candidates)), ceiling(seq_len(nrow(candidates))/PROFILE_PANELS_PER_PAGE))
  for (page in seq_along(pages)) {
    stat <- data.table::copy(candidates[pages[[page]]])
    levels <- stat$panel_label
    dat <- data.table::copy(profile_guides[panel_label %in% levels])
    bin <- data.table::copy(profile_bins[panel_label %in% levels])
    pred <- data.table::copy(sigmoid_curves[panel_label %in% levels])
    smooth <- data.table::copy(profile_loess[panel_label %in% levels])
    for (object in list(dat, bin, pred, smooth, stat)) object[, panel_label := factor(panel_label, levels = levels)]
    pB <- ggplot2::ggplot(dat, ggplot2::aes(y_pred, LFC.adj)) +
      ggplot2::geom_point(ggplot2::aes(colour = parent_score_display), alpha = 0.47, size = 1.5) +
      ggplot2::geom_hline(yintercept = CUTOFF_SENS, linewidth = 0.5, linetype = "dashed", colour = "grey30") +
      ggplot2::geom_hline(yintercept = CUTOFF_PROT, linewidth = 0.5, linetype = "dashed", colour = "grey60") +
      ggplot2::geom_hline(data = stat[is.finite(severity_profile)], ggplot2::aes(yintercept = severity_profile), linewidth = 0.6, linetype = "dotdash", colour = "black") +
      ggplot2::geom_line(data = bin, ggplot2::aes(y_bin_mid, med_LFC_bin, group = panel_label), linewidth = 0.95, colour = "black", inherit.aes = FALSE) +
      ggplot2::geom_point(data = bin, ggplot2::aes(y_bin_mid, med_LFC_bin, size = n_bin), shape = 21, fill = "white", colour = "black", stroke = 0.7, inherit.aes = FALSE)
    if (nrow(pred)) pB <- pB +
      ggplot2::geom_ribbon(data = pred, ggplot2::aes(x = dose, ymin = pmin, ymax = pmax), alpha = 0.18, inherit.aes = FALSE) +
      ggplot2::geom_line(data = pred, ggplot2::aes(dose, p), linewidth = 0.85, inherit.aes = FALSE)
    if (SHOW_LOESS_IN_PROFILES && nrow(smooth)) pB <- pB +
      ggplot2::geom_line(data = smooth, ggplot2::aes(y_pred, loess_LFC), colour = "#2166AC", linetype = "longdash", linewidth = 0.65, inherit.aes = FALSE)
    pB <- pB +
      ggplot2::geom_vline(data = stat[is.finite(EK50)], ggplot2::aes(xintercept = EK50), linewidth = 0.7, colour = "grey10") +
      ggplot2::geom_vline(data = stat[is.finite(sensitizing_threshold)], ggplot2::aes(xintercept = sensitizing_threshold), linewidth = 0.6, linetype = "dotted", colour = "grey20") +
      ggplot2::geom_vline(data = stat[is.finite(protective_threshold)], ggplot2::aes(xintercept = protective_threshold), linewidth = 0.6, linetype = "dotted", colour = "grey40") +
      ggplot2::geom_label(data = stat, ggplot2::aes(x = -Inf, y = -Inf, label = lab_txt), hjust = -0.05, vjust = -0.10,
                           size = 2.5, label.size = 0.1, fill = "white", inherit.aes = FALSE) +
      ggplot2::scale_colour_gradient2(low = "#053061", mid = "grey90", high = "#67001F", midpoint = 0, na.value = "grey85", name = "Parent-guide\nsequence score") +
      ggplot2::scale_size_continuous(name = "Guides per bin", range = c(1.2, 3.8)) +
      ggplot2::facet_wrap(~panel_label, ncol = PROFILE_NCOL, scales = "free_y") +
      ggplot2::theme_minimal(base_size = 12, base_family = FONT_FAMILY) +
      ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"), panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::labs(title = "Additional sensitizing candidates", x = "Predicted guide efficacy", y = "Adjusted interaction log2 fold change") +
      ggplot2::coord_cartesian(xlim = c(Y_FLOOR, Y_CAP + 0.02), ylim = PROFILE_Y_LIMITS)
    profile_pages[[page]] <- pB
  }
} else {
  profile_pages[[1L]] <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "No additional sensitizing candidates pass the selected settings.") + ggplot2::theme_void()
}

# 14. A small diagnostic shows accepted/rejected pre-candidates while tuning LOESS.
diagnostic <- metrics[pre_candidate == TRUE]
pD <- if (!nrow(diagnostic)) {
  ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "No profiles pass the pre-candidate criteria.") + ggplot2::theme_void()
} else ggplot2::ggplot(diagnostic, ggplot2::aes(loess_width_sens, loess_auc_sens, colour = new_sensitizing_candidate)) +
  ggplot2::geom_vline(xintercept = LOESS_MIN_WIDTH_SENS, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_hline(yintercept = LOESS_AUC_MIN_SENS, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_point(na.rm = TRUE) + ggplot2::facet_wrap(~display_label, scales = "free") +
  ggplot2::scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#2166AC"), name = "Selected") +
  ggplot2::theme_classic(base_size = 10, base_family = FONT_FAMILY) +
  ggplot2::labs(x = "Total LOESS width below sensitizing threshold", y = "LOESS area below sensitizing threshold")

# 15. Save this run's settings and data without touching step 03 or step 17.
# Remove only this script's older profile-page images when overwriting, so a run
# with fewer candidates does not leave stale pages suggesting extra candidates.
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
own_outputs <- list.files(OUTPUT_DIR, pattern = "^(FigS13|Analysis_settings)[_.]", full.names = TRUE)
if (!OVERWRITE && length(own_outputs)) stop("Output exists; set OVERWRITE <- TRUE or change OUTPUT_DIR.")
if (OVERWRITE) unlink(list.files(OUTPUT_DIR, pattern = "^FigS13B_profiles_page[0-9]+\\.(svg|pdf|png)$", full.names = TRUE))
setting_names <- c("RECURRENCE_EXCLUSION", "Y_FLOOR", "Y_CAP", "BIN_WIDTH", "MIN_BIN_N", "MIN_TOTAL_GUIDES", "MIN_QUALIFIED_BINS",
 "PARTIAL_REPRESSION_MAX", "CUTOFF_SENS", "RECOVER_THR_SENS", "K_SENS", "R_SENS", "MIN_LASTBIN_N", "STRONG_LASTBIN", "CUTOFF_PROT", "K_PROT",
 "LOESS_SPAN", "LOESS_DEGREE", "LOESS_FAMILY", "LOESS_GRID_N", "LOESS_MIN_BINS", "LOESS_AUC_X_MIN", "LOESS_AUC_X_MAX",
 "USE_LOESS_AUC_FILTER", "LOESS_AUC_MIN_SENS", "USE_LOESS_WIDTH_FILTER", "LOESS_MIN_WIDTH_SENS", "DO_SIGMOID", "ANCHOR_EFFICACY",
 "ANCHOR_RESIDUAL_MULTIPLIER", "ANCHOR_W", "EMAX_THRESHOLD", "SHOW_LOESS_IN_PROFILES", "PROFILE_PANELS_PER_PAGE", "PROFILE_NCOL", "SEED")
setting_values <- vapply(setting_names, function(name) as.character(get(name)), character(1))
inputs <- c(paths, paths_guides)
if (file.exists(GUIDE_SCORE_FILE)) inputs <- c(inputs, guide_scores = GUIDE_SCORE_FILE)
settings <- data.table::data.table(setting = c(setting_names, "DOSE_GRID", "drc_version", "drc_probe", "input_zero_LFC", "SCORE_CLIP_QUANTILES",
                                               paste0("input:", names(inputs)), paste0("md5:", names(inputs)), "R_version"),
  value = c(setting_values, paste(DOSE_GRID, collapse = ";"), if (DO_SIGMOID) as.character(utils::packageVersion("drc")) else "not used",
            as.character(anchor_weight_probe), "raw LFC == 0 excluded", paste(SCORE_CLIP_QUANTILES, collapse = ";"),
            unname(inputs), unname(tools::md5sum(inputs)), R.version.string))
write_table(metrics, file.path(OUTPUT_DIR, "FigS13_selection_audit.tsv.gz"))
write_table(scatter, file.path(OUTPUT_DIR, "FigS13A_points.tsv.gz"))
write_table(candidates, file.path(OUTPUT_DIR, "FigS13_candidates.tsv"))
write_table(profile_guides, file.path(OUTPUT_DIR, "FigS13B_guides.tsv.gz"))
write_table(profile_bins, file.path(OUTPUT_DIR, "FigS13B_bins.tsv.gz"))
write_table(loess_curves, file.path(OUTPUT_DIR, "FigS13_LOESS_curves.tsv.gz"))
write_table(sigmoid_curves, file.path(OUTPUT_DIR, "FigS13_sigmoid_curves.tsv.gz"))
write_table(settings, file.path(OUTPUT_DIR, "Analysis_settings.tsv"))
save_figure(pA, "FigS13A_additional_sensitizing", SCATTER_WIDTH_IN, SCATTER_HEIGHT_IN)
for (page in seq_along(profile_pages)) {
  n_panels <- if (nrow(candidates)) min(PROFILE_PANELS_PER_PAGE, nrow(candidates) - (page - 1L)*PROFILE_PANELS_PER_PAGE) else 1L
  save_figure(profile_pages[[page]], sprintf("FigS13B_profiles_page%02d", page), PROFILE_WIDTH_IN,
              max(4, ceiling(n_panels/PROFILE_NCOL)*PROFILE_ROW_HEIGHT_IN + 1.2))
}
if (SAVE_LOESS_DIAGNOSTIC) save_figure(pD, "FigS13_LOESS_diagnostic", 10, 6)
FigS13A_plot <- pA
FigS13B_plots <- profile_pages
FigS13_LOESS_diagnostic_plot <- pD
message("Saved Fig. S13: ", nrow(candidates), " additional sensitizing gene-condition records (",
        data.table::uniqueN(candidates$PAO1_ID), " genes); ", nrow(keys), " profiles checked by LOESS.\n", OUTPUT_DIR)


# Report the files saved by this run, then display the main panels in RStudio.
plot_stems <- c("FigS13A_additional_sensitizing", sprintf("FigS13B_profiles_page%02d", seq_along(profile_pages)),
                if (SAVE_LOESS_DIAGNOSTIC) "FigS13_LOESS_diagnostic")
plot_files <- unlist(lapply(plot_stems, function(stem) paste0(stem, c(".svg", ".pdf", ".png"))), use.names = FALSE)
data_files <- c("FigS13_selection_audit.tsv.gz", "FigS13A_points.tsv.gz", "FigS13_candidates.tsv",
                "FigS13B_guides.tsv.gz", "FigS13B_bins.tsv.gz", "FigS13_LOESS_curves.tsv.gz",
                "FigS13_sigmoid_curves.tsv.gz", "Analysis_settings.tsv")
message("Output folder: ", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = TRUE),
        "\nPlot files:\n  ", paste(plot_files, collapse = "\n  "),
        "\nSource tables and settings:\n  ", paste(data_files, collapse = "\n  "))
if (interactive()) {
  if (PRINT_IN_RSTUDIO) print(FigS13A_plot)
  if (PRINT_PROFILE_PAGES_IN_RSTUDIO) for (plot in FigS13B_plots) print(plot)
  message("RStudio: use the Plots pane back/forward arrows to browse panels A and B.\n",
          "To redisplay: print(FigS13A_plot); print(FigS13B_plots[[1]]).\n",
          "LOESS diagnostic: print(FigS13_LOESS_diagnostic_plot).")
}
