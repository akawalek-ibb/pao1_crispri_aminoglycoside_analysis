# Figure 6E: NPN fluorescence relative to each batch's control.
# Read the NPN/OD600 measurements, normalise, test and plot in this one script.
# Run from the RStudio project's working directory:
# source("Figure_scripts/Fig6/AMG_Fig6E_NPN_permeability.R")
# Requires data.table, ggplot2 and svglite; openxlsx only for optional XLSX input.
# No preprocessing script or prepared NPN result files are required.
#
# Input values must already be NPN fluorescence / OD600, not raw fluorescence
# and not percentages of a control. Wide input: Guide, Batch, control_for_batch,
# Rep1, Rep2, ... . Mark controls with 1/TRUE and other guides with 0/FALSE.
# Long input may instead provide Rep and NPN_over_OD600. Blank cells are missing.
# TSV (including .gz), CSV and the first sheet of XLSX are accepted.
# Each replicate is treated as an independent observation; matching Rep labels
# across guides do not imply paired measurements.
#
# Flow:
# 1. Set the input, statistical choices and appearance below.
# 2. Read and check the measurements; reshape wide replicate columns to rows.
# 3. Express each observation as a percentage of its batch's control mean.
# 4. Calculate group means/SDs and independent two-sided Welch tests on log2
#    NPN/OD600 values; apply Holm correction separately within each batch.
# 5. Optionally select guides, then plot means, SDs, observations and test labels.
# 6. Save PNG/PDF/SVG figures and their underlying TSV data beside the script.
# Tests/correction use all input guides before any display-only selection.
# Missing observations are excluded. Nonpositive measured values stay in the
# display but cannot enter the log2 tests; the input audit records them.

# 1. INPUT AND USER SETTINGS ---------------------------------------------------
# This path matches the measurement file in your Fig6 input folder.
INPUT_FILE <- "Figure_scripts/Fig6/input_AMG_Fig6/Exp2_NPN_measurements_for_manuscript.tsv"

OUTPUT_DIR <- "Figure_scripts/Fig6/out_AMG_Fig6E_NPN_permeability"

# Replicate-column names in wide input. Long input uses its Rep column instead.
REP_COLUMN_PATTERN <- "^Rep"

# Stop on repeated Batch/Guide/Rep IDs. Use "mean" only for technical readings
# that should first be averaged into one independent observation.
DUPLICATE_REPLICATES <- "error"

# Minimum positive observations in BOTH guide and control groups for a test.
MIN_MEASUREMENTS <- 2L

# Correct across all noncontrol guides separately within each batch.
P_ADJUST_METHOD <- "holm"

# Leave NULL to show all groups, or provide guide IDs/labels to show a subset.
# Selecting fewer guides changes only the figure, not normalisation or test correction.
GUIDES_TO_PLOT <- NULL

# Appearance and output size.
FONT_FAMILY <- "Arial"

BAR_FILL <- "#E1E5E8"

OUTLINE_COLOUR <- "#30343B"

POINT_COLOUR <- "#247C83"

# Show or hide the saved nonsignificant test labels.
SHOW_NS <- TRUE

SHOW_BATCH_LABELS <- TRUE

WIDTH_IN <- 7.5

HEIGHT_IN <- 5

DPI <- 300

SEED <- 123L

PRINT_IN_RSTUDIO <- TRUE

OVERWRITE <- TRUE

# 2. LOCAL FUNCTIONS AND INPUT CHECKS ------------------------------------------
# Check that the required R packages are installed before starting.
amg_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) 
    stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

# Read a tab-separated table, including compressed tables ending in .gz.
amg_read <- function(path) {
  if (!file.exists(path)) 
    stop("Missing input: ", path, call. = FALSE)
  con <- if (endsWith(path, ".gz")) 
    gzfile(path, "rt")
  else file(path, "rt")
  on.exit(close(con))
  data.table::as.data.table(utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE, 
    quote = "\"", comment.char = "", na.strings = c("", "NA")))
}

# Stop if a batch/guide group or replicate is missing its identifier or appears twice.
amg_unique <- function(x, columns, label) {
  if (anyNA(x[, columns, with = FALSE]) || anyDuplicated(x[, columns, with = FALSE])) 
    stop(label, ": missing or duplicate key: ", paste(columns, collapse = ", "), call. = FALSE)
  invisible(x)
}

# Save the data behind the figure as tab-separated tables.
amg_write <- function(tables, directory, overwrite = TRUE, gzip = character()) {
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

# Record the settings, input filenames and file checksums so the run can be traced.
amg_settings <- function(settings, inputs) {
  data.table::data.table(setting = c(names(settings), paste0("input:", names(inputs)), paste0("md5:", 
    names(inputs)), "R_version"), value = c(as.character(settings), unname(inputs), unname(tools::md5sum(inputs)), 
    R.version.string))
}

# First build all figure files in a temporary folder, then copy the finished files
# into the output folder. Save PNG, PDF and SVG, plus any accompanying data tables.
rf_export <- function(plot, tables, output, stem, width, height, dpi, overwrite = TRUE, grob = FALSE) {
  paths <- c(paste0(stem, c(".png", ".pdf", ".svg")), if (length(tables)) paste0(names(tables), 
    ".tsv"))
  if (!overwrite && any(file.exists(file.path(output, paths)))) 
    stop("An output exists; set OVERWRITE <- TRUE.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(".rnaseq_figure_", tmpdir = output)
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE))
  if (!grob) {
    ggplot2::ggsave(file.path(stage, paths[1L]), plot, width = width, height = height, dpi = dpi, 
      bg = "white")
    ggplot2::ggsave(file.path(stage, paths[2L]), plot, device = grDevices::cairo_pdf, width = width, 
      height = height, bg = "white")
    ggplot2::ggsave(file.path(stage, paths[3L]), plot, device = svglite::svglite, width = width, 
      height = height, bg = "white")
  } else {
    grDevices::png(file.path(stage, paths[1L]), width = width, height = height, units = "in", 
      res = dpi, type = "cairo")
    tryCatch(grid::grid.draw(plot), finally = grDevices::dev.off())
    grDevices::cairo_pdf(file.path(stage, paths[2L]), width = width, height = height)
    tryCatch(grid::grid.draw(plot), finally = grDevices::dev.off())
    svglite::svglite(file.path(stage, paths[3L]), width = width, height = height)
    tryCatch(grid::grid.draw(plot), finally = grDevices::dev.off())
  }
  if (length(tables)) 
    amg_write(tables, stage, TRUE)
  for (name in paths) if (!file.copy(file.path(stage, name), file.path(output, name), overwrite = overwrite)) 
    stop("Could not save ", name)
  invisible(plot)
}

# Format guide labels so the gene name is italic and the guide number stays upright.
rf_guide_plotmath <- function(x) {
  out <- vapply(x, function(label) {
    m <- regmatches(label, regexec("^g_(.+?)([0-9]+)(-m[0-9]+)?$", label, perl = TRUE))[[1L]]
    if (length(m)) 
      paste0("paste('g_',italic('", m[2L], "'), '", m[3L], m[4L], "')")
    else paste0("'", gsub("'", "", label, fixed = TRUE), "'")
  }, character(1))
  parse(text = out)
}

# Require the named input columns before reshaping the measurements.
amg_required <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) 
    stop(label, " lacks: ", paste(missing, collapse = ", "), call. = FALSE)
}


# Read your measurement table. Non-TSV formats are optional conveniences;
# the supplied manuscript measurements are read directly as TSV.
npn_read <- function(path) {
  if (!file.exists(path)) 
    stop("Missing input: ", path)
  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    amg_require("openxlsx")
    return(data.table::as.data.table(openxlsx::read.xlsx(path, sheet = 1, check.names = FALSE)))
  }
  if (grepl("\\.csv$", path, ignore.case = TRUE)) 
    return(data.table::as.data.table(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, 
      na.strings = c("", "NA"))))
  amg_read(path)
}


# Interpret the control flag explicitly; empty or unexpected flags are errors.
rf_logical <- function(x, label) {
  x <- tolower(trimws(as.character(x)))
  if (anyNA(x) || any(!x %in% c("true", "false", "1", "0"))) 
    stop("Invalid logical values: ", label)
  x %in% c("true", "1")
}


# Stars use adjusted p values: * <0.05, ** <0.01, *** <0.001.
rf_stars <- function(p) {
  ifelse(!is.finite(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, 
    "*", ""))))
}


# Convert original names such as accC_g1m3 to the display label g_accC1-m3.
rf_guide_name <- function(x) {
  x <- sub("^([A-Za-z0-9]+)_g([0-9]+)m([0-9]+)$", "g_\\1\\2-m\\3", x)
  sub("^([A-Za-z0-9]+)_g([0-9]+)$", "g_\\1\\2", x)
}


amg_require(c("data.table", "ggplot2", "svglite"))

if (!DUPLICATE_REPLICATES %in% c("error", "mean") || MIN_MEASUREMENTS < 2L || !P_ADJUST_METHOD %in% 
  stats::p.adjust.methods) 
  stop("Invalid NPN settings.")
# Read Guide/Batch IDs and check that every group has a consistent control
# flag. source_row retains its original location for checking the input later.
raw <- npn_read(INPUT_FILE)
amg_required(raw, c("Guide", "Batch", "control_for_batch"), INPUT_FILE)
if (!nrow(raw)) 
  stop("NPN input is empty. Copy your final NPN/OD600 measurements into the input folder.")
raw[, `:=`(Guide = trimws(as.character(Guide)), Batch = trimws(as.character(Batch)), source_row = seq_len(.N))]
if (anyNA(raw$Guide) || anyNA(raw$Batch) || any(!nzchar(raw$Guide)) || any(!nzchar(raw$Batch))) 
  stop("Missing Guide/Batch in NPN input.")
raw[, `:=`(is_control, rf_logical(control_for_batch, "control_for_batch"))]
if (any(raw[, data.table::uniqueN(is_control), by = .(Guide, Batch)]$V1 != 1L)) 
  stop("Conflicting control flags for the same Guide/Batch.")
# Keep long input as one observation per row, or turn Rep1, Rep2, ... columns
# into that same layout. A missing cell is not converted into a measured zero.
if (all(c("Rep", "NPN_over_OD600") %in% names(raw))) {
  cells <- raw[, .(Guide, Batch, is_control, Rep = as.character(Rep), NPN_over_OD600, source_row)]
} else {
  columns <- grep(REP_COLUMN_PATTERN, names(raw), value = TRUE)
  if (!length(columns)) 
    stop("No Rep columns or long-format NPN_over_OD600 field.")
  cells <- data.table::melt(raw, id.vars = c("Guide", "Batch", "is_control", "source_row"), 
    measure.vars = columns, variable.name = "Rep", value.name = "NPN_over_OD600", variable.factor = FALSE)
}
# Check numeric measurements and count missing/nonpositive cells per group.
# Preserve observed nonpositive ratios for display; log2 testing below uses
# only positive ratios. Nonnumeric or infinite nonmissing values are errors.
value_text <- trimws(as.character(cells$NPN_over_OD600))
missing <- is.na(value_text) | value_text %in% c("", "NA", "NaN", "NULL")
value <- suppressWarnings(as.numeric(value_text))
if (any(!missing & !is.finite(value))) 
  stop("Invalid nonnumeric/nonfinite NPN measurement.")
cells[, `:=`(NPN_over_OD600, value)]
audit <- cells[, .(n_input_cells = .N, n_missing_cells = sum(!is.finite(NPN_over_OD600)), n_nonpositive_cells = sum(is.finite(NPN_over_OD600) & 
  NPN_over_OD600 <= 0)), by = .(Batch, Guide, is_control)]
data <- cells[is.finite(NPN_over_OD600)]
if (!nrow(data) || anyNA(data$Rep) || any(!nzchar(data$Rep))) 
  stop("No measured NPN ratios or missing replicate IDs.")
# A Batch/Guide/Rep key should identify one independent observation.
# Technical averaging occurs only when explicitly requested in the settings.
duplicates <- anyDuplicated(data[, .(Batch, Guide, Rep)])
if (duplicates && DUPLICATE_REPLICATES == "error") 
  stop("Duplicate Batch/Guide/Rep. Resolve biological IDs or explicitly choose technical averaging.")
if (duplicates) {
  data <- data[, .(NPN_over_OD600 = mean(NPN_over_OD600), n_source_values = .N, source_row = min(source_row), 
    is_control = is_control[1L]), by = .(Batch, Guide, Rep)]
} else {
  data[, `:=`(n_source_values, 1L)]
}
# 3. NORMALISE WITHIN EACH BATCH
# Calculate the mean NPN/OD600 of the marked control observations in each
# batch, then divide each observation by that mean and multiply by 100.
# Control observations keep their own variation; only their mean equals 100%.
# Keep the original input order for the eventual batch and guide arrangement.
controls <- data[is_control == TRUE, .(batch_control_mean = mean(NPN_over_OD600), n_control_observations = .N), 
  by = Batch]
if (!setequal(controls$Batch, data$Batch) || any(!is.finite(controls$batch_control_mean) | 
  controls$batch_control_mean <= 0)) 
  stop("Each batch needs a positive mean control ratio.")
data <- merge(data, controls, by = "Batch", all.x = TRUE, sort = FALSE)
data[, `:=`(percent_control = 100 * NPN_over_OD600/batch_control_mean, stat_eligible = NPN_over_OD600 > 
  0, guide_label = rf_guide_name(Guide), batch_order = match(Batch, unique(raw$Batch)))]
data[, `:=`(guide_order, match(Guide, unique(raw$Guide)))]
# 4. SUMMARISE AND TEST
# Bars show the arithmetic mean percentage, with SD across observations.
# Tests compare each noncontrol guide with controls from its own batch using
# the ORIGINAL log2(NPN/OD600) values, not the plotted percentages.
# The Welch test is independent, two-sided and does not assume equal variance.
summary <- data[, .(guide_label = guide_label[1L], is_control = is_control[1L], batch_order = batch_order[1L], 
  guide_order = guide_order[1L], n = .N, mean_percent_control = mean(percent_control), sd_percent_control = stats::sd(percent_control), 
  mean_NPN_over_OD600 = mean(NPN_over_OD600)), by = .(Batch, Guide)]
tests <- summary[is_control == FALSE, .(Batch, Guide)]
statistics <- data.table::rbindlist(lapply(seq_len(nrow(tests)), function(i) {
  batch_i <- tests$Batch[i]
  guide_i <- tests$Guide[i]
  treatment <- data[Batch == batch_i & Guide == guide_i & !is_control & stat_eligible, log2(NPN_over_OD600)]
  control <- data[Batch == batch_i & is_control & stat_eligible, log2(NPN_over_OD600)]
  out <- data.table::data.table(Batch = batch_i, Guide = guide_i, n_treatment = length(treatment), 
    n_control = length(control), log2_difference = NA_real_, statistic = NA_real_, df = NA_real_, 
    p_value = NA_real_, status = "Insufficient positive measurements")
  if (length(treatment) >= MIN_MEASUREMENTS && length(control) >= MIN_MEASUREMENTS) {
    test <- tryCatch(stats::t.test(treatment, control, paired = FALSE, var.equal = FALSE, 
      alternative = "two.sided"), error = identity)
    if (inherits(test, "error")) 
      out[, `:=`(status, conditionMessage(test))]
    else out[, `:=`(log2_difference = mean(treatment) - mean(control), statistic = unname(test$statistic), 
      df = unname(test$parameter), p_value = test$p.value, status = "Tested")]
  }
  out
}), fill = TRUE)
# Adjust the noncontrol-guide p values separately within each batch. This is
# done before selecting guides for display, so hiding a guide cannot change
# another guide's adjusted p value or significance label.
if (nrow(statistics)) {
  statistics[, `:=`(p_adjusted, stats::p.adjust(p_value, method = P_ADJUST_METHOD)), by = Batch]
  statistics[, `:=`(significance, ifelse(!is.finite(p_adjusted), "Not tested", ifelse(p_adjusted < 
    0.05, rf_stars(p_adjusted), "ns")))]
} else statistics <- data.table::data.table(Batch = character(), Guide = character(), n_treatment = integer(), 
  n_control = integer(), log2_difference = numeric(), statistic = numeric(), df = numeric(), 
  p_value = numeric(), p_adjusted = numeric(), significance = character(), status = character())
data.table::setorder(data, batch_order, guide_order, Rep)
data.table::setorder(summary, batch_order, guide_order)

# Confirm unique summary/test groups and observation IDs before plotting.
amg_unique(summary, c("Batch", "Guide"), "NPN summary")
amg_unique(data, c("Batch", "Guide", "Rep"), "NPN measurements")
amg_unique(statistics, c("Batch", "Guide"), "NPN statistics")

# 5. SELECT GUIDES FOR DISPLAY AND DRAW THE FIGURE ------------------------------
# All normalisation and statistical correction is already complete. Leave
# GUIDES_TO_PLOT as NULL to display every batch/guide group in the input.
if (!is.null(GUIDES_TO_PLOT)) {
  if (!length(GUIDES_TO_PLOT) || any(!GUIDES_TO_PLOT %in% c(summary$Guide, summary$guide_label))) 
    stop("Unknown/empty NPN guide selection.")
  summary <- summary[Guide %in% GUIDES_TO_PLOT | guide_label %in% GUIDES_TO_PLOT]
  selected_keys <- paste(summary$Batch, summary$Guide, sep = "\r")
  data <- data[paste(Batch, Guide, sep = "\r") %in% selected_keys]
  statistics <- statistics[paste(Batch, Guide, sep = "\r") %in% selected_keys]
}

# Keep the batch and guide order recorded from the input.
data.table::setorder(summary, batch_order, guide_order)

summary[, `:=`(x, seq_len(.N))]

# Give each group an x-axis position, then place its replicate points at that position.
keys <- paste(summary$Batch, summary$Guide, sep = "\r")

data[, `:=`(x, summary$x[match(paste(Batch, Guide, sep = "\r"), keys)])]

# Check that measurements and summaries describe the same groups and replicate counts.
if (anyNA(data$x)) stop("NPN observations and summary groups disagree.")

n <- data[, .N, by = .(Batch, Guide)]

if (any(n$N != summary$n[match(paste(n$Batch, n$Guide, sep = "\r"), keys)])) stop("NPN group counts disagree.")

# Attach the just-calculated within-batch adjusted test labels to matching bars.
stars <- merge(statistics, summary[, .(Batch, Guide, x, mean_percent_control, sd_percent_control)], 
  by = c("Batch", "Guide"), all.x = TRUE, sort = FALSE)

if (!SHOW_NS) stars <- stars[significance != "ns"]

# Leave enough vertical space above the SD bars for the test and batch labels.
span <- range(c(0, data$percent_control, summary$mean_percent_control + summary$sd_percent_control), 
  finite = TRUE)

padding <- max(5, diff(span) * 0.055)

if (nrow(stars)) stars[, `:=`(y, mean_percent_control + ifelse(is.finite(sd_percent_control), sd_percent_control, 0) + padding)]

# Find the boundaries of each batch so different experiments can be separated visually.
blocks <- summary[, .(xmin = min(x), xmax = max(x), center = mean(range(x))), by = Batch]

# Draw the mean as a bar, SD as an error bar and each independent observation
# as a point. The dotted line at 100% marks the batch-control mean.
p <- ggplot2::ggplot(summary, ggplot2::aes(x = x, y = mean_percent_control)) +
  ggplot2::geom_col(width = 0.64, 
  fill = BAR_FILL, colour = OUTLINE_COLOUR, linewidth = 0.45) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_percent_control - 
  sd_percent_control, ymax = mean_percent_control + sd_percent_control), width = 0.22, colour = OUTLINE_COLOUR, 
  linewidth = 0.45, na.rm = TRUE) +
  ggplot2::geom_point(data = data, ggplot2::aes(x = x, y = percent_control), 
  inherit.aes = FALSE, position = ggplot2::position_jitter(width = 0.13, height = 0, seed = SEED), 
  colour = POINT_COLOUR, size = 1.8, alpha = 0.9) +
  ggplot2::geom_hline(yintercept = 100, linetype = "dotted", 
  colour = "grey60", linewidth = 0.3) +
  ggplot2::scale_x_continuous(breaks = summary$x, labels = rf_guide_plotmath(summary$guide_label), 
  expand = ggplot2::expansion(add = 0.5)) +
  ggplot2::theme_minimal(base_size = 10, base_family = FONT_FAMILY) + 
  ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"), axis.title.y = ggplot2::element_text(face = "bold")) + 
  ggplot2::labs(x = NULL, y = "NPN fluorescence / OD600\n(% of batch control)")

# Add the significance labels calculated above, without repeating the tests.
if (nrow(stars)) p <- p +
  ggplot2::geom_text(data = stars, ggplot2::aes(x = x, y = y, label = significance), 
  inherit.aes = FALSE, size = 3, family = FONT_FAMILY)

# Separate different batches with vertical lines; optionally show their labels.
if (nrow(blocks) > 1L) p <- p +
  ggplot2::geom_vline(xintercept = head(blocks$xmax, -1L) + 0.5, 
  colour = "grey65", linetype = "dashed", linewidth = 0.3)

if (SHOW_BATCH_LABELS && nrow(blocks) > 1L) p <- p +
  ggplot2::geom_text(data = blocks, ggplot2::aes(x = center, 
  y = span[2L] + 3 * padding, label = paste("Batch", Batch)), inherit.aes = FALSE, family = FONT_FAMILY, 
  size = 3)

# 6. SAVE FIGURES AND THE DATA BEHIND THEM --------------------------------------
# Keep the original figure/source-data names. The audit describes all input
# groups, including any guides hidden by a display selection. No separate
# preprocessing output folder or duplicate intermediate tables are produced.
settings <- amg_settings(c(
  test = "Independent two-sided Welch t-test on log2(NPN/OD600)",
  P_ADJUST_METHOD = P_ADJUST_METHOD,
  correction_family = "All noncontrol guides separately within each batch; before plot selection",
  normalization = "100 * ratio / mean biological batch-control ratio",
  missing_cells = "Excluded as missing; observed zero retained for display",
  nonpositive_ratios = "Retained in display; ineligible for log2 statistics",
  DUPLICATE_REPLICATES = DUPLICATE_REPLICATES,
  MIN_MEASUREMENTS = MIN_MEASUREMENTS,
  REP_COLUMN_PATTERN = REP_COLUMN_PATTERN,
  GUIDES_TO_PLOT = if (is.null(GUIDES_TO_PLOT)) "All" else paste(GUIDES_TO_PLOT, collapse = "; "),
  SHOW_NS = SHOW_NS,
  SEED = SEED), c(measurements = INPUT_FILE))

rf_export(p, list(Fig6E_plot_summary = summary, Fig6E_plot_statistics = stars,
  Fig6E_plot_measurements = data, NPN_input_audit = audit, Analysis_settings = settings),
  OUTPUT_DIR, "Fig6E_NPN_permeability", WIDTH_IN, HEIGHT_IN, DPI, OVERWRITE)

if (PRINT_IN_RSTUDIO) print(p)

message("Saved Fig. 6E for ", nrow(summary), " guide/batch groups; ", nrow(data), " plotted observations.")

# Print the locations of the saved figure data at completion.
figure_data_files <- file.path(OUTPUT_DIR, c("Fig6E_plot_summary.tsv",
  "Fig6E_plot_statistics.tsv", "Fig6E_plot_measurements.tsv", "NPN_input_audit.tsv",
  "Analysis_settings.tsv"))
message("Underlying figure data:\n", paste0("  ", normalizePath(figure_data_files, winslash = "/"), collapse = "\n"))
