#!/usr/bin/env Rscript
# Supplementary Figure S14: aminoglycoside-associated repression thresholds.
# Run from the RStudio project root. No other R script is sourced.
#
# Compare sensitivity measured independently in untreated/T0 and treated/T0
# profiles. Colour comes from the separate treated/untreated interaction test
# used in Table S8. These are different contrasts and answer different questions.
#
# Sensitivity estimates/statuses are reused from preprocessing 11; interaction
# calls and T1 recurrence flags are reused from preprocessing 03. No curves,
# interaction cutoffs or recurrence calls are refitted or recalculated here.
# A compact file of absolute gene-level tests against T0 is included alongside
# this script because the vulnerability master does not contain their p values.
#
# Inputs: preprocessing_scripts/out_11_AMG_gene_vulnerability/
# gene_vulnerability_long.tsv.gz; preprocessing_scripts/out_03_AMG_relevant_results/
# Stouffer_calls.tsv; Figure_scripts/FigS14/input_AMG_FigS14/absolute_gene_effects.tsv.gz.
# Run: source("Figure_scripts/FigS14/AMG_FigS14_aminoglycoside_threshold_shifts.R")
# Outputs: combined/individual PNG/PDF/SVG and plotted data, exclusions and settings
# TSVs in Figure_scripts/FigS14/out_AMG_FigS14_repression_thresholds/.
# REBUILD_GENE_EFFECTS is FALSE by default; the optional refresh reads the existing
# metacounter gene results and replaces the compact absolute-effects input only.
#
# R packages used/checked: ggplot2, ggrepel, svglite.

# 1. User settings -----------------------------------------------------------
VULNERABILITY_FILE <- "preprocessing_scripts/out_11_AMG_gene_vulnerability/gene_vulnerability_long.tsv.gz"
INTERACTION_FILE <- "preprocessing_scripts/out_03_AMG_relevant_results/Stouffer_calls.tsv"
GENE_EFFECTS_FILE <- "Figure_scripts/FigS14/input_AMG_FigS14/absolute_gene_effects.tsv.gz"
OUTPUT_DIR <- "Figure_scripts/FigS14/out_AMG_FigS14_repression_thresholds"

# The supplied compact input is ready to use. Set TRUE only when deliberately
# refreshing the gene-level tests from existing metacounter median outputs.
# This optional step reads b1/b2/b4_median_melted_results.tsv.gz; it does not
# rerun metacounter, edgeR or guide fitting. Keep FALSE for routine plotting.
REBUILD_GENE_EFFECTS <- FALSE
METACOUNTER_DIR <- "edgeR_data"

FITNESS_P_STOUFFER_CUT <- 0.05
NEGATIVE_LFC_CUT <- -1
SHIFT_BAND <- 0.15                         # descriptive lines, not a significance cutoff
Y_FLOOR <- 0.05
Y_CAP <- 1
AMK_MAX_LABELS <- 20L                     # largest displayed absolute shifts, including boundary hits

PRINT_IN_RSTUDIO <- TRUE
SAVE_INDIVIDUAL_PANELS <- TRUE
OVERWRITE <- TRUE
FONT_FAMILY <- "Arial"
DPI <- 300
COMBINED_WIDTH_IN <- 16
COMBINED_HEIGHT_IN <- 9.5
PANEL_WIDTH_IN <- 7.2
PANEL_HEIGHT_IN <- 6.0

DIRECTION_COLOURS <- c("Sensitizing interaction" = "#B2182B", "Protective interaction" = "#2166AC")
SCOPE_SHAPES <- c("Recurrent interaction" = 17, "Condition-specific interaction" = 15,
  "Other non-recurrent interaction" = 18, "Other" = 16)
BACKGROUND_COLOUR <- "#BDBDBD"
BACKGROUND_POINT_SIZE <- 1.45
HIGHLIGHT_POINT_SIZE <- 2.55

LABEL_REPEL_SEED <- 1L
LABEL_REPEL_BOX_PADDING <- 0.35
LABEL_REPEL_POINT_PADDING <- 0.18
LABEL_REPEL_FORCE <- 1.5
LABEL_REPEL_FORCE_PULL <- 0.15
LABEL_REPEL_MAX_TIME <- 5
LABEL_REPEL_MAX_ITER <- 20000L

# Each row is a matched pair from the same batch and time point. Both absolute
# profiles share T0. T2 is shown but does not contribute to T1 recurrence.
COMPARISONS <- utils::read.table(text = "analysis_id|batch|panel_label|control_condition|treated_condition|call_column
B4_TOB_0.025_T1|b4|TOB 0.025 (B4)|None_0_T1 - None_0_T0|TOB_0.025_T1 - None_0_T0|TOB_0.025_call
B4_TOB_0.1_T1|b4|TOB 0.1 (B4)|None_0_T1 - None_0_T0|TOB_0.1_T1 - None_0_T0|TOB_0.1_call
B1_TOB_0.4_T1|b1|TOB 0.4 (B1)|None_0_T1 - None_0_T0|TOB_0.4_T1 - None_0_T0|TOB_0.4_call
B4_TOB_1_T1|b4|TOB 1 (B4)|None_0_T1 - None_0_T0|TOB_1_T1 - None_0_T0|TOB_1_call
B1_TOB_0.4_T2|b1|TOB 0.4 T2 (B1)|None_0_T2 - None_0_T0|TOB_0.4_T2 - None_0_T0|TOB_T2_call
B2_AMK_2_T1|b2|AMK 2 (B2)|None_0_T1 - None_0_T0|AMK_2_T1 - None_0_T0|AMK_call",
  sep = "|", header = TRUE, quote = "", stringsAsFactors = FALSE)

# 2. Read the existing numerical results -------------------------------------
packages <- c("ggplot2", "ggrepel", "svglite")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required packages: ", paste(missing_packages, collapse = ", "))

# Local reader for the plain or compressed TSV files used by this script.
read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  con <- if (endsWith(path, ".gz")) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE,
    quote = "\"", comment.char = "", na.strings = c("", "NA"))
}

if (REBUILD_GENE_EFFECTS) {
  gene_tables <- list()
  for (batch_i in unique(COMPARISONS$batch)) {
    path_i <- file.path(METACOUNTER_DIR, paste0(batch_i, "_median_melted_results.tsv.gz"))
    raw_i <- read_tsv(path_i)
    needed <- c("PAO1_ID", "condition", "medLFC", "p_stouffer")
    if (!all(needed %in% names(raw_i))) stop("Metacounter gene table lacks: ", paste(setdiff(needed, names(raw_i)), collapse = ", "))
    selected_conditions <- unique(c(COMPARISONS$control_condition[COMPARISONS$batch == batch_i],
      COMPARISONS$treated_condition[COMPARISONS$batch == batch_i]))
    raw_i <- raw_i[raw_i$condition %in% selected_conditions, needed]
    raw_i$batch <- batch_i
    gene_tables[[batch_i]] <- raw_i[c("PAO1_ID", "batch", "condition", "medLFC", "p_stouffer")]
  }
  gene_effects <- do.call(rbind, gene_tables)
  if (anyDuplicated(gene_effects[c("PAO1_ID", "batch", "condition")])) stop("Duplicate metacounter gene/contrast records.")
  dir.create(dirname(GENE_EFFECTS_FILE), recursive = TRUE, showWarnings = FALSE)
  if (!OVERWRITE && file.exists(GENE_EFFECTS_FILE)) stop("Compact gene-effects file already exists; set OVERWRITE <- TRUE.")
  con <- if (endsWith(GENE_EFFECTS_FILE, ".gz")) gzfile(GENE_EFFECTS_FILE, "wt") else file(GENE_EFFECTS_FILE, "wt")
  tryCatch(utils::write.table(gene_effects, con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA"),
    finally = close(con))
}

gene_effects <- read_tsv(GENE_EFFECTS_FILE)
vulnerability <- read_tsv(VULNERABILITY_FILE)
interactions <- read_tsv(INTERACTION_FILE)
required_columns <- list(
  gene_effects = c("PAO1_ID", "batch", "condition", "medLFC", "p_stouffer"),
  vulnerability = c("PAO1_ID", "gene_name", "batch", "condition", "severity", "sensitivity_y_pred", "sensitivity_status"),
  interactions = c("PAO1_ID", "gene_name", "recurrent_sensitizing", "recurrent_protective", COMPARISONS$call_column))
for (name_i in names(required_columns)) {
  missing_i <- setdiff(required_columns[[name_i]], names(get(name_i)))
  if (length(missing_i)) stop(name_i, " lacks: ", paste(missing_i, collapse = ", "))
}
if (anyNA(gene_effects$PAO1_ID) || anyDuplicated(gene_effects[c("PAO1_ID", "batch", "condition")])) stop("Invalid gene-effects keys.")
if (anyNA(vulnerability$PAO1_ID) || anyDuplicated(vulnerability[c("PAO1_ID", "batch", "condition")])) stop("Invalid vulnerability keys.")
if (anyNA(interactions$PAO1_ID) || anyDuplicated(interactions$PAO1_ID)) stop("Invalid interaction-table gene keys.")

# Only Table S8's sensitizing/protective calls are highlighted. Do not replace
# these with tests against T0 or add a new multiple-testing filter here.
for (column_i in COMPARISONS$call_column) {
  if (any(!is.na(interactions[[column_i]]) & !interactions[[column_i]] %in% c("Sensitizing", "Protective", "Not significant"))) {
    stop("Unrecognised Table S8 call in ", column_i)
  }
}
for (flag_i in c("recurrent_sensitizing", "recurrent_protective")) {
  flag <- toupper(as.character(interactions[[flag_i]]))
  if (any(!is.na(flag) & !flag %in% c("TRUE", "FALSE", "1", "0"))) stop("Invalid recurrence flag: ", flag_i)
  interactions[[flag_i]] <- flag %in% c("TRUE", "1")
}
interactions$n_significant_comparisons <- rowSums(as.data.frame(lapply(interactions[COMPARISONS$call_column],
  function(x) x %in% c("Sensitizing", "Protective"))))

# Check the imported flags against the five T1 calls, leaving TOB T2 out.
# This verifies the source; the imported flags remain the recurrence authority.
t1_columns <- COMPARISONS$call_column[COMPARISONS$analysis_id != "B1_TOB_0.4_T2"]
for (direction_i in c("Sensitizing", "Protective")) {
  n_i <- rowSums(as.data.frame(lapply(interactions[t1_columns], function(x) x %in% direction_i)))
  if (any(interactions[[paste0("recurrent_", tolower(direction_i))]] != (n_i >= 2))) {
    stop("Table S8 recurrence flags disagree with the five T1 calls.")
  }
}

# 3. Build the six matched pairs and select evaluable genes --------------------
# Sensitivity comes directly from the shared, genome-wide calculation:
# >=10 guides and >=4 populated efficacy bins, then the sustained -1 crossing
# rule. No full-match/mismatch series is recalculated in this figure script.
all_panels <- list()
for (i in seq_len(nrow(COMPARISONS))) {
  cfg <- COMPARISONS[i, ]
  control <- gene_effects[gene_effects$batch == cfg$batch & gene_effects$condition == cfg$control_condition, ]
  treated <- gene_effects[gene_effects$batch == cfg$batch & gene_effects$condition == cfg$treated_condition, ]
  if (!nrow(control) || !nrow(treated)) stop("Missing absolute gene effects for ", cfg$panel_label)

  # Keep the union of gene records so that missing/non-evaluable profiles are
  # visible in the audit export, even though they cannot enter the plot.
  panel <- data.frame(PAO1_ID = union(control$PAO1_ID, treated$PAO1_ID), stringsAsFactors = FALSE)
  panel$analysis_id <- cfg$analysis_id
  panel$panel_label <- cfg$panel_label
  panel$batch <- cfg$batch
  for (role_i in c("control", "treated")) {
    effect_i <- if (role_i == "control") control else treated
    condition_i <- if (role_i == "control") cfg$control_condition else cfg$treated_condition
    profile_i <- vulnerability[vulnerability$batch == cfg$batch & vulnerability$condition == condition_i, ]
    if (!nrow(profile_i)) stop("Missing master vulnerability profile for ", cfg$panel_label, ": ", condition_i)
    rows <- match(panel$PAO1_ID, effect_i$PAO1_ID)
    fits <- match(panel$PAO1_ID, profile_i$PAO1_ID)
    panel[[paste0(role_i, "_condition")]] <- condition_i
    panel[[paste0(role_i, "_medLFC")]] <- effect_i$medLFC[rows]
    panel[[paste0(role_i, "_p_stouffer")]] <- effect_i$p_stouffer[rows]
    panel[[paste0(role_i, "_sensitivity")]] <- profile_i$sensitivity_y_pred[fits]
    status_i <- profile_i$sensitivity_status[fits]
    status_i[is.na(status_i)] <- "no_guide_data"
    panel[[paste0(role_i, "_sensitivity_status")]] <- status_i
    if (role_i == "control") panel$gene_name <- profile_i$gene_name[fits]

    # Catch an outdated/mismatched compact gene table before plotting it with
    # a different master analysis. Compare all available full-match effects.
    effect <- panel[[paste0(role_i, "_medLFC")]]
    severity <- profile_i$severity[fits]
    comparable <- !is.na(rows)
    if (any(is.finite(effect[comparable]) != is.finite(severity[comparable])) ||
        any(abs(effect[comparable] - severity[comparable]) > 1e-8, na.rm = TRUE)) {
      stop("Gene effects and vulnerability master disagree for ", cfg$panel_label,
        ". Refresh the compact input from the same metacounter run (REBUILD_GENE_EFFECTS <- TRUE).")
    }
    value <- panel[[paste0(role_i, "_sensitivity")]]
    if (any(is.finite(value) & !(status_i %in% "estimated")) ||
        any(status_i %in% "estimated" & !is.finite(value))) stop("Inconsistent sensitivity/status in ", cfg$panel_label)
    if (any(is.finite(value) & (value < Y_FLOOR | value > Y_CAP))) stop("Sensitivity outside the expected efficacy range.")

    # Only an evaluable profile with no sustained crossing is placed at 1.
    # Insufficient guides/bins are missing profiles, not boundary estimates.
    panel[[paste0(role_i, "_sensitivity_plot")]] <- ifelse(status_i %in% "estimated", value,
      ifelse(status_i %in% "no_sustained_crossing", Y_CAP, NA_real_))
    p_i <- panel[[paste0(role_i, "_p_stouffer")]]
    panel[[paste0(role_i, "_significant_negative")]] <- is.finite(p_i) & p_i < FITNESS_P_STOUFFER_CUT &
      is.finite(effect) & effect <= NEGATIVE_LFC_CUT
  }

  panel$retained_by_fitness_filter <- panel$control_significant_negative | panel$treated_significant_negative
  panel$plotted <- panel$retained_by_fitness_filter & is.finite(panel$control_sensitivity_plot) &
    is.finite(panel$treated_sensitivity_plot)
  panel$exclusion_reason <- ifelse(panel$plotted, "plotted",
    ifelse(!panel$retained_by_fitness_filter, "no significant negative fitness in either condition",
      "non-evaluable sensitivity in at least one condition"))

  # Directions are condition-specific. A gene that is recurrent elsewhere but
  # not significant here remains grey. Recurrence must match THIS hit's direction.
  hit_rows <- match(panel$PAO1_ID, interactions$PAO1_ID)
  panel$interaction_call <- interactions[[cfg$call_column]][hit_rows]
  panel$interaction_call[is.na(panel$interaction_call)] <- "Not significant"
  panel$interaction_hit <- panel$interaction_call %in% c("Sensitizing", "Protective")
  panel$interaction_direction <- ifelse(panel$interaction_call == "Sensitizing", "Sensitizing interaction",
    ifelse(panel$interaction_call == "Protective", "Protective interaction", "Not a significant interaction"))
  recurrent_sensitizing <- interactions$recurrent_sensitizing[hit_rows] %in% TRUE
  recurrent_protective <- interactions$recurrent_protective[hit_rows] %in% TRUE
  panel$recurrent_same_direction <- (panel$interaction_call == "Sensitizing" & recurrent_sensitizing) |
    (panel$interaction_call == "Protective" & recurrent_protective)
  panel$n_significant_comparisons <- interactions$n_significant_comparisons[hit_rows]
  panel$n_significant_comparisons[is.na(panel$n_significant_comparisons)] <- 0L
  panel$interaction_plot_scope <- "Other"
  panel$interaction_plot_scope[panel$interaction_hit & panel$n_significant_comparisons == 1] <- "Condition-specific interaction"
  panel$interaction_plot_scope[panel$interaction_hit & panel$n_significant_comparisons >= 2] <- "Other non-recurrent interaction"
  panel$interaction_plot_scope[panel$recurrent_same_direction] <- "Recurrent interaction"

  names_i <- interactions$gene_name[hit_rows]
  use_name <- !is.na(names_i) & nzchar(names_i)
  panel$gene_name[use_name] <- names_i[use_name]
  panel$label_txt <- ifelse(is.na(panel$gene_name) | !nzchar(panel$gene_name), panel$PAO1_ID,
    gsub("'", "", panel$gene_name, fixed = TRUE))
  panel$sensitivity_delta_display <- panel$treated_sensitivity_plot - panel$control_sensitivity_plot
  both_numeric <- panel$control_sensitivity_status == "estimated" & panel$treated_sensitivity_status == "estimated"
  panel$sensitivity_delta_numeric <- ifelse(both_numeric, panel$treated_sensitivity - panel$control_sensitivity, NA_real_)
  both_no_crossing <- panel$control_sensitivity_status == "no_sustained_crossing" &
    panel$treated_sensitivity_status == "no_sustained_crossing"
  panel$sensitivity_shift_class <- ifelse(!panel$plotted, NA_character_,
    ifelse(both_no_crossing, "Both no sustained crossing",
      ifelse(panel$sensitivity_delta_display <= -SHIFT_BAND, "Lower threshold under antibiotic",
        ifelse(panel$sensitivity_delta_display >= SHIFT_BAND, "Higher threshold under antibiotic", "Threshold-concordant"))))

  # Label every plotted interaction hit in the five TOB panels. For AMK the
  # caption's largest-shift rule is applied across ALL hits, including boundary
  # points. A boundary value alone does not guarantee an AMK label.
  panel$selected_for_label <- panel$plotted & panel$interaction_hit
  if (cfg$analysis_id == "B2_AMK_2_T1") {
    candidates <- which(panel$selected_for_label)
    ranked <- candidates[order(-abs(panel$sensitivity_delta_display[candidates]), panel$PAO1_ID[candidates])]
    panel$selected_for_label <- seq_len(nrow(panel)) %in% head(ranked, AMK_MAX_LABELS)
  }
  all_panels[[i]] <- panel
}

all_genes <- do.call(rbind, all_panels)
all_genes$panel_label <- factor(all_genes$panel_label, levels = COMPARISONS$panel_label)
plot_genes <- all_genes[all_genes$plotted, ]
if (!nrow(plot_genes)) stop("No evaluable matched profiles passed the fitness filter.")

# 4. Arrange labels and summarise the displayed data ---------------------------
# A point displayed at (1,1) has two unresolved thresholds. Do not count its
# zero display difference as evidence that its true thresholds are equal.
labels <- plot_genes[plot_genes$selected_for_label, ]
labels$label_lane <- ifelse(labels$control_sensitivity_status == "no_sustained_crossing", "right",
  ifelse(labels$treated_sensitivity_status == "no_sustained_crossing", "top", "inside"))
labels$label_x <- labels$control_sensitivity_plot
labels$label_y <- labels$treated_sensitivity_plot
for (panel_i in COMPARISONS$panel_label) {
  top <- which(labels$panel_label == panel_i & labels$label_lane == "top")
  top <- top[order(labels$control_sensitivity_plot[top], labels$label_txt[top])]
  if (length(top)) {
    labels$label_x[top] <- if (length(top) == 1L) labels$control_sensitivity_plot[top] else seq(0.08, 0.98, length.out = length(top))
    labels$label_y[top] <- 1.035
  }
  right <- which(labels$panel_label == panel_i & labels$label_lane == "right")
  right <- right[order(labels$treated_sensitivity_plot[right], labels$label_txt[right])]
  if (length(right)) {
    labels$label_x[right] <- 1.035
    labels$label_y[right] <- if (length(right) == 1L) labels$treated_sensitivity_plot[right] else seq(0.10, 0.94, length.out = length(right))
  }
}

summary_rows <- list()
for (i in seq_len(nrow(COMPARISONS))) {
  a <- all_genes[all_genes$analysis_id == COMPARISONS$analysis_id[i], ]
  z <- a[a$plotted, ]
  numeric_pairs <- z[z$control_sensitivity_status == "estimated" & z$treated_sensitivity_status == "estimated", ]
  summary_rows[[i]] <- data.frame(analysis_id = COMPARISONS$analysis_id[i], panel_label = COMPARISONS$panel_label[i],
    n_control_significant_negative = sum(a$control_significant_negative),
    n_treated_significant_negative = sum(a$treated_significant_negative),
    n_retained_fitness_either = sum(a$retained_by_fitness_filter), n_evaluable_both = nrow(z),
    n_estimated_both = nrow(numeric_pairs), n_current_interaction_hits = sum(z$interaction_hit),
    n_current_recurrent_same_direction = sum(z$recurrent_same_direction),
    n_condition_specific = sum(z$interaction_plot_scope == "Condition-specific interaction"),
    n_other_nonrecurrent = sum(z$interaction_plot_scope == "Other non-recurrent interaction"),
    n_labels = sum(z$selected_for_label),
    n_both_no_crossing = sum(z$sensitivity_shift_class == "Both no sustained crossing"),
    median_control_sensitivity = if (nrow(numeric_pairs)) stats::median(numeric_pairs$control_sensitivity) else NA_real_,
    median_treated_sensitivity = if (nrow(numeric_pairs)) stats::median(numeric_pairs$treated_sensitivity) else NA_real_)
}
statistics <- do.call(rbind, summary_rows)

# 5. Draw the combined figure and individual panels ----------------------------
# The original colours, point sizes, fonts, dimensions, band and label lanes
# are retained. The shape scale has explicit limits so absent categories do
# not break legend construction in an individual panel.
neutral_band <- data.frame(x = seq(Y_FLOOR, Y_CAP, length.out = 300L))
neutral_band$lower <- pmax(Y_FLOOR, neutral_band$x - SHIFT_BAND)
neutral_band$upper <- pmin(Y_CAP, neutral_band$x + SHIFT_BAND)
# Draw reference segments only inside the tested efficacy range. An unbounded
# abline can extend into the label margins when clipping is disabled for labels.
equality_line <- data.frame(x = Y_FLOOR, y = Y_FLOOR, xend = Y_CAP, yend = Y_CAP)
shift_lines <- data.frame(x = c(Y_FLOOR, Y_FLOOR + SHIFT_BAND), y = c(Y_FLOOR + SHIFT_BAND, Y_FLOOR),
  xend = c(Y_CAP - SHIFT_BAND, Y_CAP), yend = c(Y_CAP, Y_CAP - SHIFT_BAND))

# One plotting function keeps the combined and single-panel layers identical.
# All biological selection and labelling decisions have already been made above.
draw_threshold_plot <- function(points, label_data, combined = FALSE) {
  inside <- label_data[label_data$label_lane == "inside", ]
  top <- label_data[label_data$label_lane == "top", ]
  right <- label_data[label_data$label_lane == "right", ]
  text_size <- if (combined) 2.25 else 2.4
  p <- ggplot2::ggplot(points, ggplot2::aes(control_sensitivity_plot, treated_sensitivity_plot)) +
    ggplot2::geom_ribbon(data = neutral_band, ggplot2::aes(x = x, ymin = lower, ymax = upper),
      inherit.aes = FALSE, fill = "grey94", colour = NA) +
    ggplot2::geom_segment(data = equality_line, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE, linewidth = 0.55, colour = "grey40") +
    ggplot2::geom_segment(data = shift_lines, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE, linetype = "dashed", linewidth = 0.40, colour = "grey65") +
    ggplot2::geom_point(ggplot2::aes(shape = "Other"), colour = BACKGROUND_COLOUR,
      size = BACKGROUND_POINT_SIZE, alpha = 0.55) +
    ggplot2::geom_point(data = points[points$interaction_hit, ],
      ggplot2::aes(colour = interaction_direction, shape = interaction_plot_scope),
      size = HIGHLIGHT_POINT_SIZE, alpha = 0.95, stroke = 0.75) +
    ggplot2::scale_colour_manual(values = DIRECTION_COLOURS, limits = names(DIRECTION_COLOURS),
      breaks = names(DIRECTION_COLOURS), labels = c("Sensitizing", "Protective"),
      drop = FALSE, name = "Interaction direction") +
    ggplot2::scale_shape_manual(values = SCOPE_SHAPES, limits = names(SCOPE_SHAPES),
      breaks = names(SCOPE_SHAPES), labels = c("Recurrent", "Condition-specific", "Other non-recurrent", "No current interaction"),
      drop = FALSE, name = "Interaction recurrence",
      guide = ggplot2::guide_legend(override.aes = list(
        colour = c("grey30", "grey30", "grey30", BACKGROUND_COLOUR),
        alpha = c(0.95, 0.95, 0.95, 0.55),
        size = c(HIGHLIGHT_POINT_SIZE, HIGHLIGHT_POINT_SIZE, HIGHLIGHT_POINT_SIZE, BACKGROUND_POINT_SIZE)))) +
    ggrepel::geom_text_repel(data = inside,
      ggplot2::aes(control_sensitivity_plot, treated_sensitivity_plot, label = label_txt),
      size = text_size, fontface = "italic", colour = "black", inherit.aes = FALSE,
      box.padding = LABEL_REPEL_BOX_PADDING, point.padding = LABEL_REPEL_POINT_PADDING,
      force = LABEL_REPEL_FORCE, force_pull = LABEL_REPEL_FORCE_PULL,
      max.time = LABEL_REPEL_MAX_TIME, max.iter = LABEL_REPEL_MAX_ITER, max.overlaps = Inf,
      min.segment.length = 0, segment.color = "grey55", segment.alpha = 0.75,
      segment.size = 0.30, seed = LABEL_REPEL_SEED, direction = "both",
      xlim = c(Y_FLOOR, 0.995), ylim = c(Y_FLOOR, 0.995)) +
    ggplot2::geom_segment(data = top,
      ggplot2::aes(x = control_sensitivity_plot, y = 1, xend = label_x, yend = 1.015),
      inherit.aes = FALSE, colour = "grey55", alpha = 0.70, linewidth = 0.35) +
    ggplot2::geom_text(data = top, ggplot2::aes(label_x, label_y, label = label_txt),
      inherit.aes = FALSE, angle = 90, hjust = 0, vjust = 0.5, size = text_size,
      fontface = "italic", colour = "black") +
    ggplot2::geom_segment(data = right,
      ggplot2::aes(x = 1, y = treated_sensitivity_plot, xend = 1.015, yend = label_y),
      inherit.aes = FALSE, colour = "grey55", alpha = 0.70, linewidth = 0.35) +
    ggplot2::geom_text(data = right, ggplot2::aes(label_x, label_y, label = label_txt),
      inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = text_size, fontface = "italic", colour = "black") +
    ggplot2::scale_x_continuous(breaks = seq(0.1, 1, by = 0.2)) +
    ggplot2::scale_y_continuous(breaks = seq(0.1, 1, by = 0.2)) +
    ggplot2::coord_equal(xlim = c(Y_FLOOR, Y_CAP), ylim = c(Y_FLOOR, Y_CAP), clip = "off") +
    ggplot2::theme_classic(base_family = FONT_FAMILY, base_size = 9) +
    ggplot2::theme(legend.position = "right", legend.box.spacing = grid::unit(50, "pt"),
      axis.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold", margin = ggplot2::margin(t = 3, b = 45)),
      strip.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", margin = ggplot2::margin(b = 40))) +
    ggplot2::labs(x = "Silencing-efficiency threshold (untreated)", y = "Silencing-efficiency threshold (antibiotic)")
  if (combined) {
    p <- p + ggplot2::facet_wrap(~panel_label, ncol = 3, nrow = 2) +
      ggplot2::theme(panel.spacing = grid::unit(3.5, "lines"), plot.margin = ggplot2::margin(120, 105, 25, 10))
  } else {
    p <- p + ggplot2::labs(title = as.character(unique(points$panel_label))) +
      ggplot2::theme(plot.margin = ggplot2::margin(140, 125, 25, 10))
  }
  p
}

FigS14_plot <- draw_threshold_plot(plot_genes, labels, combined = TRUE)
FigS14_panels <- list()
if (SAVE_INDIVIDUAL_PANELS) {
  for (i in seq_len(nrow(COMPARISONS))) {
    id <- COMPARISONS$analysis_id[i]
    FigS14_panels[[id]] <- draw_threshold_plot(plot_genes[plot_genes$analysis_id == id, ], labels[labels$analysis_id == id, ])
  }
}

# 6. Save figures and machine-readable underlying data -------------------------
# No workbook or RDS duplicate is created. The full table records exclusions;
# the plotted table retains true estimates separately from display coordinates.
input_paths <- c(gene_effects = GENE_EFFECTS_FILE, vulnerability = VULNERABILITY_FILE, interactions = INTERACTION_FILE)
settings <- data.frame(setting = c(paste0("input:", names(input_paths)), paste0("md5:", names(input_paths)),
    "FITNESS_P_STOUFFER_CUT", "NEGATIVE_LFC_CUT", "SHIFT_BAND", "AMK_MAX_LABELS", "AMK_label_rule",
    "recurrence_rule", "condition_specific_rule", "R_version", "ggplot2_version", "ggrepel_version"),
  value = c(unname(input_paths), unname(tools::md5sum(input_paths)), FITNESS_P_STOUFFER_CUT, NEGATIVE_LFC_CUT,
    SHIFT_BAND, AMK_MAX_LABELS, "largest absolute displayed shifts across all AMK hits; ties by PAO1_ID",
    "Table S8 same-direction recurrence among five T1 comparisons; T2 excluded",
    "significant in only one of the six plotted comparisons", R.version.string,
    as.character(utils::packageVersion("ggplot2")), as.character(utils::packageVersion("ggrepel"))))
tables <- list(all_gene_statuses = all_genes, plotted_genes = plot_genes, labels = labels,
  statistics = statistics, comparison_manifest = COMPARISONS, Analysis_settings = settings)
table_paths <- file.path(OUTPUT_DIR, paste0("FigS14_", names(tables),
  ifelse(names(tables) == "all_gene_statuses", ".tsv.gz", ".tsv")))
plots <- c(list(FigS14_aminoglycoside_threshold_shifts = FigS14_plot), FigS14_panels)
figure_stems <- "FigS14_aminoglycoside_threshold_shifts"
if (length(FigS14_panels)) figure_stems <- c(figure_stems, paste0("FigS14_", names(FigS14_panels)))
figure_paths <- as.vector(outer(file.path(OUTPUT_DIR, figure_stems), c(".png", ".pdf", ".svg"), paste0))
if (!OVERWRITE && any(file.exists(c(table_paths, figure_paths)))) stop("Outputs already exist; set OVERWRITE <- TRUE.")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (i in seq_along(tables)) {
  con <- if (endsWith(table_paths[i], ".gz")) gzfile(table_paths[i], "wt") else file(table_paths[i], "wt")
  tryCatch(utils::write.table(tables[[i]], con, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA"),
    finally = close(con))
}
for (i in seq_along(plots)) {
  width_i <- if (i == 1L) COMBINED_WIDTH_IN else PANEL_WIDTH_IN
  height_i <- if (i == 1L) COMBINED_HEIGHT_IN else PANEL_HEIGHT_IN
  stem_i <- file.path(OUTPUT_DIR, figure_stems[i])
  ggplot2::ggsave(paste0(stem_i, ".png"), plots[[i]], width = width_i, height = height_i, units = "in", dpi = DPI, bg = "white")
  ggplot2::ggsave(paste0(stem_i, ".pdf"), plots[[i]], width = width_i, height = height_i, units = "in",
    device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(paste0(stem_i, ".svg"), plots[[i]], width = width_i, height = height_i, units = "in",
    device = svglite::svglite, bg = "white")
}
if (isTRUE(PRINT_IN_RSTUDIO)) {
  for (p_i in FigS14_panels) print(p_i)
  print(FigS14_plot)                       # leave the complete figure visible last
}
print(statistics, row.names = FALSE)
message("Saved Fig. S14 figures:\n", paste(figure_paths, collapse = "\n"))
message("Underlying figure data:\n", paste(table_paths, collapse = "\n"))
