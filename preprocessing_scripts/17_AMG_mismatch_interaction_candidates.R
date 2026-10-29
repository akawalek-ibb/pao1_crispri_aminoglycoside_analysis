# Identify exploratory interactions supported by mismatch-guide response profiles.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/17_AMG_mismatch_interaction_candidates.R")
# Inputs:
# - preprocessing_scripts/out_03_AMG_relevant_results/b*_AG_guides.tsv.gz,
#   AG_gene_results.tsv.gz, AG_gene_recurrence.tsv and Comparison_manifest.tsv.
# Upstream / execution order:
# - preprocessing_scripts/03_AMG_relevant_results.R with guide exports enabled; its fitted guide
#   inputs come from preprocessing_scripts/02_AMG_metacounter.R.
# Outputs:
# - preprocessing_scripts/out_17_AMG_mismatch_interaction_candidates/ —
#   Mismatch_based_evidence.tsv, candidate_support.tsv, candidate_bins.tsv.gz,
#   candidate_fit_curves.tsv.gz and settings.
# Interpretation / operation:
# - Exclude conventional same-direction calls and genes recurrent in either direction.
# - Apply the configured bin-crossing and LOESS excursion rules; fit only selected candidates.
# - Requires data.table and drc. This candidate analysis is separate from script 11
#   vulnerability estimates.

INPUT_DIR <- "preprocessing_scripts/out_03_AMG_relevant_results"

OUTPUT_DIR <- "preprocessing_scripts/out_17_AMG_mismatch_interaction_candidates"

Y_FLOOR <- 0.05

Y_CAP <- 1

BIN_WIDTH <- 0.05

CUTOFF_SENS <- -1

RECOVER_THR_SENS <- -1

K_SENS <- 2L

R_SENS <- 2L

MIN_LASTBIN_N <- 8L

STRONG_LASTBIN <- -1

CUTOFF_PROT <- 1

K_PROT <- 2L

MIN_BIN_N <- 3L

MIN_TOTAL_GUIDES <- 10L

MIN_QUALIFIED_BINS <- 4L

PARTIAL_REPRESSION_MAX <- 0.8

LOESS_SPAN <- 0.75

LOESS_DEGREE <- 2L

LOESS_FAMILY <- "gaussian"

LOESS_GRID_N <- 401L

LOESS_MIN_BINS <- 5L

LOESS_AUC_X_MIN <- 0.1

LOESS_AUC_X_MAX <- PARTIAL_REPRESSION_MAX

LOESS_AUC_MIN_SENS <- 0.05

LOESS_AUC_MIN_PROT <- 0.05

LOESS_MIN_WIDTH_SENS <- 0.2

LOESS_MIN_WIDTH_PROT <- 0.2

# Anchor residual multiplier: 3 gives nine times one guide's influence in the
# least-squares objective. The installed drc weights argument is checked below.
ANCHOR_RESIDUAL_MULTIPLIER <- 3

EMAX_THRESHOLD <- 1

DOSE_GRID <- seq(0.1, 1, by = 0.01)

OVERWRITE <- TRUE

# Reproducibility: drc 3.0-1 multiplies continuous-response residuals by weights,
# whereas drc 4.0-0 divides them by weights before squaring. Thus argument 3
# upweights the anchor in the older release but downweights it in the newer one.
# Check the installed objective directly: a unit residual with argument 3 gives
# 9 under multiplication or 1/9 under division. Preserve the intended anchor
# strength across both conventions; all real guide weights remain 1.
anchor_weight_probe <- getFromNamespace("drmEMls", "drc")(
  dose = 1, resp = 1, multCurves = function(dose, parm) rep(0, length(dose)),
  startVec = 0, robustFct = function(residual) residual^2,
  weights = 3, rmNA = FALSE
)$opfct(0)

if (isTRUE(all.equal(as.numeric(anchor_weight_probe), 9, tolerance = 1e-10))) {
  ANCHOR_W <- ANCHOR_RESIDUAL_MULTIPLIER
} else if (isTRUE(all.equal(as.numeric(anchor_weight_probe), 1/9, tolerance = 1e-10))) {
  ANCHOR_W <- 1 / ANCHOR_RESIDUAL_MULTIPLIER
} else {
  stop("Unrecognised drc weighting convention; check the installed package before fitting.")
}
message("drc ", utils::packageVersion("drc"), ": anchor weights argument = ",
        signif(ANCHOR_W, 6), "; relative least-squares weight = ", ANCHOR_RESIDUAL_MULTIPLIER^2)

# Measure negative/positive threshold crossings from bin medians for one gene.
# Require enough total guides and enough bins each supported by MIN_BIN_N guides.
# For eligible profiles, smooth the qualified bins with guide-count-weighted LOESS;
# integrate the area beyond each cutoff and the total width of that excursion
# within the partial-repression range. Width need not be one continuous segment.
s8_profile_metrics <- function(bins) {
  data.table::setorder(bins, bin_i)
  si <- find_sustained_cross_idx_le(bins$bin_i, bins$med_LFC_bin, bins$n_bin, CUTOFF_SENS, RECOVER_THR_SENS, 
    K_SENS, R_SENS, MIN_BIN_N, MIN_LASTBIN_N, STRONG_LASTBIN)
  pi <- find_run_start_adj_ge(bins$bin_i, bins$med_LFC_bin, bins$n_bin, CUTOFF_PROT, K_PROT, 
    MIN_BIN_N)
  st <- interpolate_threshold(bins$y_bin_mid, bins$med_LFC_bin, si, CUTOFF_SENS)
  pt <- interpolate_threshold(bins$y_bin_mid, bins$med_LFC_bin, pi, CUTOFF_PROT)
  qualified <- bins[n_bin >= MIN_BIN_N]
  partial <- qualified[y_bin_mid <= PARTIAL_REPRESSION_MAX]
  qc <- sum(bins$n_bin) >= MIN_TOTAL_GUIDES && nrow(qualified) >= MIN_QUALIFIED_BINS
  out <- data.table::data.table(sensitizing_threshold = st, protective_threshold = pt, n_guides_in_bins = sum(bins$n_bin), 
    n_qualified_bins = nrow(qualified), threshold_qc = qc, minimum_partial_bin_median = if (nrow(partial)) 
      min(partial$med_LFC_bin)
    else NA_real_, maximum_partial_bin_median = if (nrow(partial)) 
      max(partial$med_LFC_bin)
    else NA_real_, loess_auc_sens = NA_real_, loess_auc_prot = NA_real_, loess_width_sens = NA_real_, 
    loess_width_prot = NA_real_, loess_status = "not_eligible")
  if (!qc || (!is.finite(st) && !is.finite(pt))) 
    return(out)
  if (nrow(qualified) < LOESS_MIN_BINS || data.table::uniqueN(qualified$y_bin_mid) < LOESS_MIN_BINS) {
    out$loess_status <- "insufficient_bins"
    return(out)
  }
  fit <- tryCatch(stats::loess(med_LFC_bin ~ y_bin_mid, data = as.data.frame(qualified), weights = n_bin, 
    span = LOESS_SPAN, degree = LOESS_DEGREE, family = LOESS_FAMILY, na.action = stats::na.exclude, 
    control = stats::loess.control(surface = "direct", statistics = "approximate")), error = function(e) NULL)
  if (is.null(fit)) {
    out$loess_status <- "fit_failed"
    return(out)
  }
  lo <- max(Y_FLOOR, min(qualified$y_bin_mid))
  hi <- min(Y_CAP, max(qualified$y_bin_mid))
  if (hi <= lo) {
    out$loess_status <- "insufficient_range"
    return(out)
  }
  grid <- seq(lo, hi, length.out = LOESS_GRID_N)
  prediction <- suppressWarnings(stats::predict(fit, newdata = data.frame(y_bin_mid = grid)))
  keep <- is.finite(prediction) & grid >= LOESS_AUC_X_MIN & grid <= LOESS_AUC_X_MAX
  grid <- grid[keep]
  prediction <- as.numeric(prediction[keep])
  if (length(grid) < 2L) {
    out$loess_status <- "no_auc_domain"
    return(out)
  }
  out[, `:=`(loess_auc_sens = trapz_numeric(grid, pmax(0, CUTOFF_SENS - prediction)), loess_auc_prot = trapz_numeric(grid, 
    pmax(0, prediction - CUTOFF_PROT)), loess_width_sens = trapz_numeric(grid, as.numeric(prediction <= 
    CUTOFF_SENS)), loess_width_prot = trapz_numeric(grid, as.numeric(prediction >= CUTOFF_PROT)), 
    loess_status = "estimated")]
  out
}

# Local functions for this analysis.
s8_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) 
    stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

s8_required <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) 
    stop(label, " lacks: ", paste(missing, collapse = ", "), call. = FALSE)
}

s8_read <- function(path, columns = NULL) {
  if (!file.exists(path)) 
    stop("Missing input: ", path, call. = FALSE)
  read_once <- function(nrows = -1L, classes = NA) {
    con <- if (endsWith(path, ".gz")) 
      gzfile(path, "rt")
    else file(path, "rt")
    on.exit(close(con))
    utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = "", 
      na.strings = c("", "NA"), nrows = nrows, colClasses = classes)
  }
  if (is.null(columns)) 
    return(data.table::as.data.table(read_once()))
  header <- names(read_once(0L))
  s8_required(as.data.frame(stats::setNames(rep(list(character()), length(header)), header)), 
    columns, path)
  data.table::as.data.table(read_once(classes = ifelse(header %in% columns, "character", "NULL")))
}

s8_unique <- function(x, keys, label) {
  if (anyNA(x[, keys, with = FALSE]) || anyDuplicated(x[, keys, with = FALSE])) 
    stop(label, ": missing or duplicate key: ", paste(keys, collapse = ", "), call. = FALSE)
  invisible(x)
}

s8_settings <- function(settings, inputs = character()) {
  data.table::data.table(setting = c(names(settings), paste0("input:", names(inputs)), paste0("md5:", 
    names(inputs)), "R_version"), value = c(as.character(settings), unname(inputs), unname(tools::md5sum(inputs)), 
    R.version.string))
}

s8_setting <- function(directory, name, numeric = FALSE) {
  x <- s8_read(file.path(directory, "Analysis_settings.tsv"))
  value <- x$value[x$setting == name]
  if (length(value) != 1L || is.na(value)) 
    stop("Missing or ambiguous setting: ", name)
  if (numeric) 
    as.numeric(value)
  else as.character(value)
}

s8_write <- function(tables, directory, overwrite = TRUE, gzip = character()) {
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

# Fit a log-logistic response to individual guides plus the low-efficacy zero
# anchor. These fits describe selected candidates after the binned/LOESS screen;
# they are not used to decide which candidates pass that screen.
generate_predicted_profile <- function(temp) {
  if (is.null(temp) || nrow(temp) == 0L) {
    return(tibble::tibble())
  }
  temp <- temp %>% dplyr::mutate(w = 1)
  anchor_point <- tibble::tibble(y_pred = 0.1, LFC.adj = 0, w = ANCHOR_W)
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

# Find the first fitted half-effect crossing relative to zero, using the
# response at efficacy 1 as Emax. Return NA when the curve does not support it.
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

main <- function() {
  s8_require(c("data.table", "dplyr", "tibble", "drc"))
  if (s8_setting(INPUT_DIR, "EXPORT_GUIDE_TABLES") != "TRUE") 
    stop("Run 03 with EXPORT_GUIDE_TABLES <- TRUE before mismatch analysis.")
  suppressPackageStartupMessages(library(dplyr))
  if (Y_FLOOR <= 0 || Y_CAP <= Y_FLOOR || BIN_WIDTH <= 0 || PARTIAL_REPRESSION_MAX > Y_CAP || 
    LOESS_AUC_X_MAX <= LOESS_AUC_X_MIN || MIN_BIN_N < 1 || MIN_TOTAL_GUIDES < 1 || MIN_QUALIFIED_BINS < 
    1) 
    stop("Invalid mismatch-profile settings.")
  # 1. Read the five independent T1 interaction comparisons and primary gene
  # calls. Exclude genes already recurrent in either direction from the exploratory
  # candidate set; conventional same-direction calls are also excluded below.
  manifest <- s8_read(file.path(INPUT_DIR, "Comparison_manifest.tsv"))[count_for_recurrence == 
    TRUE]
  genes <- s8_read(file.path(INPUT_DIR, "AG_gene_results.tsv.gz"))[count_for_recurrence == TRUE]
  recurrent <- s8_read(file.path(INPUT_DIR, "AG_gene_recurrence.tsv"))
  excluded <- recurrent$PAO1_ID[recurrent$recurrent_sensitizing | recurrent$recurrent_protective]
  paths <- file.path(INPUT_DIR, paste0(unique(manifest$batch), "_AG_guides.tsv.gz"))
  guides <- data.table::rbindlist(lapply(paths, s8_read), use.names = TRUE)
  s8_required(guides, c("PAO1_ID", "comparison_id", "spacer", "type", "y_pred", "LFC", "LFC.adj"), 
    "Extracted AG guides")
  # 2. Retain full-match/mismatch guide effects with valid scores and effects.
  # This candidate workflow excludes raw LFC==0 guides before all profile analyses.
  # Pool parent-guide families by gene/comparison, with full-match efficacy set to 1.
  guides <- guides[comparison_id %in% manifest$comparison_id & !is.na(PAO1_ID) & type %in% c("perfect", 
    "mismatch") & is.finite(y_pred) & is.finite(LFC.adj) & is.finite(LFC) & LFC != 0]
  s8_unique(guides, c("comparison_id", "spacer"), "Mismatch input")
  guides[type == "perfect", `:=`(y_pred, 1)]
  guides[, `:=`(y_pred, pmin(pmax(y_pred, Y_FLOOR), Y_CAP))]
  breaks <- unique(c(seq(Y_FLOOR, Y_CAP, by = BIN_WIDTH), Y_CAP))
  mids <- (head(breaks, -1L) + tail(breaks, -1L))/2
  guides[, `:=`(bin_i, as.integer(cut(y_pred, breaks = breaks, include.lowest = TRUE, right = TRUE)))]
  # 3. Calculate median interaction effects in efficacy bins, then threshold and
  # LOESS-support metrics. Negative effects denote sensitization relative to
  # matched untreated cells; positive effects denote protection.
  bins <- guides[!is.na(bin_i), .(med_LFC_bin = stats::median(LFC.adj), n_bin = .N), by = .(PAO1_ID, 
    comparison_id, bin_i)]
  bins[, `:=`(y_bin_mid, mids[bin_i])]
  data.table::setorder(bins, PAO1_ID, comparison_id, bin_i)
  metrics <- bins[, s8_profile_metrics(data.table::copy(.SD)), by = .(PAO1_ID, comparison_id)]
  metrics <- merge(genes[, .(PAO1_ID, comparison_id, gene_name, product_name, medLFC, p_stouffer, 
    stouffer_call)], metrics, by = c("PAO1_ID", "comparison_id"), all.x = TRUE, sort = FALSE)
  metrics[, `:=`(recurrent_any, PAO1_ID %in% excluded)]
  # 4. Select candidates with a crossing within the partial-repression range
  # and sufficient LOESS excursion area/width. Apply the conventional-call and
  # recurrence exclusions as well. These are exploratory profile criteria,
  # not a replacement for the primary gene-level significance calls.
  metrics[, `:=`(new_sensitizing_candidate, !is.na(threshold_qc) & threshold_qc & is.finite(sensitizing_threshold) & 
    sensitizing_threshold <= PARTIAL_REPRESSION_MAX & stouffer_call != "Sensitizing" & !recurrent_any & 
    is.finite(loess_auc_sens) & loess_auc_sens >= LOESS_AUC_MIN_SENS & is.finite(loess_width_sens) & 
    loess_width_sens >= LOESS_MIN_WIDTH_SENS)]
  metrics[, `:=`(new_protective_candidate, !is.na(threshold_qc) & threshold_qc & is.finite(protective_threshold) & 
    protective_threshold <= PARTIAL_REPRESSION_MAX & stouffer_call != "Protective" & !recurrent_any & 
    is.finite(loess_auc_prot) & loess_auc_prot >= LOESS_AUC_MIN_PROT & is.finite(loess_width_prot) & 
    loess_width_prot >= LOESS_MIN_WIDTH_PROT)]
  candidates <- data.table::rbindlist(lapply(c("Sensitizing", "Protective"), function(direction) {
    if (direction == "Sensitizing") 
      x <- metrics[new_sensitizing_candidate == TRUE]
    else x <- metrics[new_protective_candidate == TRUE]
    z <- x[, .(comparison_id, PAO1_ID, gene_name, product_name, interaction_medLFC = medLFC, 
      conventional_gene_level_p_stouffer = p_stouffer, `interaction_threshold (sensitivity)` = if (direction == 
        "Sensitizing") 
        sensitizing_threshold
      else protective_threshold, partial_bin_extreme = if (direction == "Sensitizing") 
        minimum_partial_bin_median
      else maximum_partial_bin_median)]
    z[, `:=`(direction, direction)]
    z
  }))
  keys <- unique(candidates[, .(PAO1_ID, comparison_id)])
  fit_rows <- curves <- vector("list", nrow(keys))
  # 5. Fit only the selected gene/comparison profiles for descriptive Emax/EC50.
  # Keep warning/failure information; reject half-effect estimates below the
  # lowest observed mismatch efficacy rather than extrapolating beyond support.
  for (i in seq_len(nrow(keys))) {
    dat <- guides[PAO1_ID == keys$PAO1_ID[i] & comparison_id == keys$comparison_id[i]]
    fit_warnings <- character()
    predicted <- withCallingHandlers(generate_predicted_profile(as.data.frame(dat)), warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    estimated <- calc_ec50_from_pred(predicted, EMAX_THRESHOLD)
    minimum_mismatch <- if (any(dat$type == "mismatch")) 
      min(dat$y_pred[dat$type == "mismatch"])
    else NA_real_
    fit_status <- if (!nrow(predicted)) 
      "fit_or_prediction_failed"
    else if (is.finite(estimated$EC50)) 
      "estimated"
    else "no_valid_half_max_crossing"
    if (is.finite(estimated$EC50) && is.finite(minimum_mismatch) && estimated$EC50 < minimum_mismatch) {
      estimated$EC50 <- NA_real_
      fit_status <- "below_mismatch_support"
    }
    fit_rows[[i]] <- data.table::data.table(PAO1_ID = keys$PAO1_ID[i], comparison_id = keys$comparison_id[i], 
      EC50 = estimated$EC50, Emax = estimated$Emax, EC50_fit_status = fit_status, fit_warnings = if (length(fit_warnings)) 
        paste(unique(fit_warnings), collapse = "; ")
      else NA_character_)
    if (nrow(predicted)) 
      curves[[i]] <- data.table::data.table(PAO1_ID = keys$PAO1_ID[i], comparison_id = keys$comparison_id[i], 
        predicted)
  }
  fits <- if (length(fit_rows)) 
    data.table::rbindlist(fit_rows)
  else data.table::data.table(PAO1_ID = character(), comparison_id = character(), EC50 = numeric(), 
    Emax = numeric(), EC50_fit_status = character(), fit_warnings = character())
  metrics <- merge(metrics, fits, by = c("PAO1_ID", "comparison_id"), all.x = TRUE, sort = FALSE)
  candidates <- merge(candidates, fits, by = c("PAO1_ID", "comparison_id"), all.x = TRUE, sort = FALSE)
  columns <- c("comparison_id", "PAO1_ID", "gene_name", "product_name", "direction", "interaction_medLFC", 
    "conventional_gene_level_p_stouffer", "interaction_threshold (sensitivity)", "partial_bin_extreme", 
    "EC50", "Emax")
  candidates <- candidates[, ..columns]
  candidates[, `:=`(panel_order, manifest$plot_order[match(comparison_id, manifest$comparison_id)])]
  data.table::setorder(candidates, panel_order, direction, PAO1_ID)
  candidates[, `:=`(panel_order, NULL)]
  # 6. Export the candidate table, complete screening-support metrics and the
  # selected candidates' bins/fit curves. Task 11 remains the separate source of
  # the shared genome-wide vulnerability measures.
  candidate_bins <- merge(bins, keys, by = c("PAO1_ID", "comparison_id"), all = FALSE, sort = FALSE)
  fit_curves <- if (any(vapply(curves, function(x) !is.null(x), logical(1)))) 
    data.table::rbindlist(curves)
  else data.table::data.table(PAO1_ID = character(), comparison_id = character(), dose = numeric(), 
    p = numeric(), pmin = numeric(), pmax = numeric())
  settings <- s8_settings(c(Y_FLOOR = Y_FLOOR, Y_CAP = Y_CAP, BIN_WIDTH = BIN_WIDTH, CUTOFF_SENS = CUTOFF_SENS, 
    RECOVER_THR_SENS = RECOVER_THR_SENS, K_SENS = K_SENS, R_SENS = R_SENS, MIN_LASTBIN_N = MIN_LASTBIN_N, 
    STRONG_LASTBIN = STRONG_LASTBIN, CUTOFF_PROT = CUTOFF_PROT, K_PROT = K_PROT, MIN_BIN_N = MIN_BIN_N, 
    MIN_TOTAL_GUIDES = MIN_TOTAL_GUIDES, MIN_QUALIFIED_BINS = MIN_QUALIFIED_BINS, PARTIAL_REPRESSION_MAX = PARTIAL_REPRESSION_MAX, 
    LOESS_MIN_BINS = LOESS_MIN_BINS, LOESS_SPAN = LOESS_SPAN, LOESS_DEGREE = LOESS_DEGREE, 
    LOESS_FAMILY = LOESS_FAMILY, LOESS_GRID_N = LOESS_GRID_N, LOESS_AUC_X_MIN = LOESS_AUC_X_MIN, 
    LOESS_AUC_X_MAX = LOESS_AUC_X_MAX, LOESS_AUC_MIN_SENS = LOESS_AUC_MIN_SENS, LOESS_AUC_MIN_PROT = LOESS_AUC_MIN_PROT, 
    LOESS_MIN_WIDTH_SENS = LOESS_MIN_WIDTH_SENS, LOESS_MIN_WIDTH_PROT = LOESS_MIN_WIDTH_PROT, 
    recurrence_exclusion = "recurrent in either direction", LFC_zero_filter = "exclude original LFC == 0", 
    ANCHOR_RESIDUAL_MULTIPLIER = ANCHOR_RESIDUAL_MULTIPLIER, ANCHOR_W = ANCHOR_W, EMAX_THRESHOLD = EMAX_THRESHOLD, DOSE_GRID = paste(DOSE_GRID, collapse = ";"), 
    drc_version = as.character(utils::packageVersion("drc"))), c(gene_results = file.path(INPUT_DIR, 
    "AG_gene_results.tsv.gz"), recurrence = file.path(INPUT_DIR, "AG_gene_recurrence.tsv"), 
    manifest = file.path(INPUT_DIR, "Comparison_manifest.tsv"), stats::setNames(paths, unique(manifest$batch))))
  s8_write(list(Mismatch_based_evidence = candidates, candidate_support = metrics, candidate_bins = candidate_bins, 
    candidate_fit_curves = fit_curves, Analysis_settings = settings), OUTPUT_DIR, OVERWRITE, 
    gzip = c("candidate_support", "candidate_bins", "candidate_fit_curves"))
  message("S8 mismatch evidence: ", nrow(candidates), " candidate gene/comparison/direction records.")
  invisible(candidates)
}

main()

