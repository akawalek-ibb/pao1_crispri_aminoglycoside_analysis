#!/usr/bin/env Rscript
# Show inferred PRECISE411 iModulon activity changes in Figure 6C.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/Fig6/AMG_Fig6C_iModulon_activity.R")
# Inputs:
# - RNAseq_data/iModulon_analysis/DiMA_primary_all.tsv: all 116 modules × five contrasts with
#   q-values and coverage flags.
# - annotation_input/precise411-im_table.csv: matching module descriptions.
# Upstream / execution order:
# - Use the supplied DiMA table directly; script 21 is not required for this standalone figure.
# - To rebuild: RNAseq_analysis/AMG_RNAseq_featureCounts_edgeR.R →
#   RNAseq_data/iModulon_analysis/preprocesing/imodulon_activity_pipeline.py; reconcile the
#   configured output/input locations.
# Outputs:
# - Figure_scripts/Fig6/out_AMG_Fig6C_iModulon_activity/ — PNG/PDF and plotted
#   activity/annotation TSVs.
# Interpretation / operation:
# - Select coverage-qualified modules with q < 0.05 and |delta_activity| > 5 in any comparison;
#   display all five contrasts for selected modules.
#
# R packages used/checked: dplyr, ggplot2, readr, stringr, tibble, tidyr.

needed <- c("readr", "dplyr", "tidyr", "ggplot2", "stringr")

missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) stop("Install packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# User settings: edit these relative paths if the input files are stored elsewhere.
INPUT_FILE <- "RNAseq_data/iModulon_analysis/DiMA_primary_all.tsv"

OUTPUT_DIR <- "Figure_scripts/Fig6/out_AMG_Fig6C_iModulon_activity"

ATLAS_FILE <- "annotation_input/precise411-im_table.csv"

ATLAS_URL <- "https://legacy.imodulondb.org/dataset.html?organism=p_aeruginosa&dataset=precise411"

# A main response must pass this adjusted-p cutoff AND the activity-change cutoff below.
Q_CUTOFF <- 0.05

MIN_ABS_DELTA <- 5

PNG_DPI <- 450

BASE_SIZE <- 10

FONT_FAMILY <- "Arial"

PRINT_IN_RSTUDIO <- TRUE

# Show the reference function beside each module name; retain all qualifying modules.
SHOW_FUNCTION_IN_MAIN <- TRUE

# Wrap long labels to fit the plot without shortening the saved descriptions.
CATEGORY_WRAP <- 22L

FUNCTION_WRAP <- 44L

if (!file.exists(INPUT_FILE)) stop("Edit INPUT_FILE or pass the full DiMA_primary_all.tsv path: ", 
  INPUT_FILE)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 1. Read all saved activity estimates and adjusted p-values. No model is refitted.
d <- readr::read_tsv(INPUT_FILE, show_col_types = FALSE)

req <- c("contrast_id", "module_id", "delta_activity", "q_primary", "q_global", "projection_coverage_ok")

if (!all(req %in% names(d))) stop("Required columns missing: ", paste(setdiff(req, names(d)), collapse = ", "))

# Read the saved flag indicating whether enough reference information was available.
truth <- function(x) tolower(as.character(x)) %in% c("true", "t", "1")

# Keep the four knockdowns together, followed by the TOB-alone control.
# All five comparisons use the induced non-targeting control as their reference.
CONTRASTS <- tibble::tribble(~contrast_id, ~condition, ~x, ~axis_label, "atpIg1m2_vs_control1", "atpI", 
  1, "g_atpI1-m2", "accCg1m3_vs_control1", "accC", 2, "g_accC1-m3", "secAg1m3_vs_control1", "secA", 3, 
  "g_secA1-m3", "ftsZg1m4_vs_control1", "ftsZ", 4, "g_ftsZ1-m4", "control1_TOB_vs_control1", "TOB", 5.5, 
  "TOB\n0.4 µg/mL")

# Recognize the supported aliases for the same biological comparisons.
ALIASES <- c(atpIg1m2_vs_control1 = "atpIg1m2_vs_control1", atpI_vs_control1 = "atpIg1m2_vs_control1", accC_vs_control1 = "accCg1m3_vs_control1", 
  secA_vs_control1 = "secAg1m3_vs_control1", ftsZ_vs_control1 = "ftsZg1m4_vs_control1")

d$contrast_id_input <- d$contrast_id

d$contrast_id <- ifelse(d$contrast_id %in% names(ALIASES), unname(ALIASES[d$contrast_id]), d$contrast_id)

if (anyDuplicated(d[c("contrast_id", "module_id")])) stop("Duplicate module/contrast rows.")

if (!all(CONTRASTS$contrast_id %in% d$contrast_id)) stop("Contrasts not found: ", paste(setdiff(CONTRASTS$contrast_id, 
  d$contrast_id), collapse = ", "))

# Mark changes passing the adjusted-p cutoff and the reference-coverage check.
# Mark those also exceeding the activity-change cutoff as the main responses.
d <- d %>% inner_join(CONTRASTS, by = "contrast_id") %>% mutate(coverage_ok = truth(projection_coverage_ok), 
  significant = is.finite(q_primary) & q_primary < Q_CUTOFF & coverage_ok, primary_DiMA = significant & 
    abs(delta_activity) > MIN_ABS_DELTA, statistical_small_effect = significant & !primary_DiMA)

# Require the full 116-by-five result table, not a table of significant results only.
if (n_distinct(d$module_id) != 116L || nrow(d) != 116L * nrow(CONTRASTS)) stop("Expected the full PRECISE411 results: 116 modules by five contrasts, not a significance-filtered table.")

# Check that the supplied adjusted p-values use the correction described for this figure.
if (!isTRUE(all.equal(d$q_primary, d$q_global, tolerance = 1e-12))) stop("This figure caption assumes primary global BH correction. Match the analysis FDR scope before plotting.")

if (any(!is.finite(d$delta_activity))) stop("Missing/nonfinite activity estimates; resolve before plotting.")

# 2. Match module IDs to the original reference categories and descriptions.
# Allow harmless differences in column spelling, such as spaces or underscores.
normalize_column <- function(x) gsub("[^a-z0-9]", "", tolower(x))

pick_column <- function(dat, aliases, required = TRUE) {
  idx <- match(normalize_column(aliases), normalize_column(names(dat)))
  idx <- idx[!is.na(idx)]
  if (length(idx)) 
    return(names(dat)[idx[[1L]]])
  if (required) 
    stop("Official metadata column missing: ", paste(aliases, collapse = " / "))
  NULL
}

# The current script stops if this annotation file is missing.
# The alternative branch below is not reached when the file is absent.
if (!file.exists(ATLAS_FILE)) stop("Explicit atlas metadata file not found: ", ATLAS_FILE)

if (file.exists(ATLAS_FILE)) {
  delim <- if (grepl("\\.tsv(\\.gz)?$", ATLAS_FILE, ignore.case = TRUE)) 
    "\t"
  else ","
  atlas_raw <- readr::read_delim(ATLAS_FILE, delim = delim, col_types = cols(.default = col_character()), 
    show_col_types = FALSE)
  # Find the reference column containing the module IDs used in the result table.
  id_candidates <- names(atlas_raw)[vapply(atlas_raw, function(v) all(unique(d$module_id) %in% trimws(as.character(v))), logical(1))]
  if (!length(id_candidates)) 
    stop("No atlas column matches all result module IDs. Use the PRECISE411 metadata corresponding to the projected M matrix.")
  id_col <- id_candidates[[1L]]
  cat_col <- pick_column(atlas_raw, c("category", "functional_category"))
  fun_col <- pick_column(atlas_raw, c("function_description", "function", "description"))
  reg_col <- pick_column(atlas_raw, c("regulator", "regulators", "regulator_readable"), required = FALSE)
  atlas <- tibble(module_id = trimws(as.character(atlas_raw[[id_col]])), official_category = as.character(atlas_raw[[cat_col]]), 
    official_function = as.character(atlas_raw[[fun_col]]), official_regulator = if (is.null(reg_col)) 
      NA_character_
    else as.character(atlas_raw[[reg_col]]))
  ANNOTATION_SOURCE <- normalizePath(ATLAS_FILE, mustWork = TRUE)
} else {
  if (!all(c("category", "function_description") %in% names(d))) 
    stop("Provide precise411-im_table.csv as argument 3; the DiMA table lacks embedded atlas annotations.")
  atlas <- tibble(module_id = as.character(d$module_id), official_category = as.character(d$category), 
    official_function = as.character(d$function_description), official_regulator = if ("regulator" %in% 
      names(d)) 
      as.character(d$regulator)
    else NA_character_)
  ANNOTATION_SOURCE <- "Original PRECISE411 annotations embedded in DiMA_primary_all.tsv"
  message("Default atlas file not found; using the reference annotation columns embedded in the DiMA table.")
}

# Keep one unambiguous reference description for every module in the analysis.
atlas <- atlas %>% distinct() %>% filter(module_id %in% d$module_id)

if (anyDuplicated(atlas$module_id)) stop("Duplicate or conflicting official annotations for the same module ID.")

if (!all(unique(d$module_id) %in% atlas$module_id)) stop("Official annotation missing for result modules.")

if (any(is.na(atlas$official_category) | !nzchar(trimws(atlas$official_category)))) stop("Official category missing; provide complete PRECISE411 metadata. No substitute categories will be invented.")

# 3. Select every module with a sufficiently large, significant change in at least
# one comparison. No modules are added or removed by biological preference.
module_summary <- d %>% group_by(module_id) %>% summarise(n_FDR_significant_contrasts = sum(significant), 
  n_large_significant_contrasts = sum(primary_DiMA), max_abs_delta = max(abs(delta_activity)), 
  minimum_q = if (any(is.finite(q_primary))) min(q_primary[is.finite(q_primary)]) else NA_real_, 
  selected_for_main = any(primary_DiMA), .groups = "drop")

SELECTED_MODULES <- module_summary %>% filter(selected_for_main) %>% pull(module_id)

d <- d %>% left_join(module_summary %>% dplyr::select(module_id, selected_for_main), by = "module_id")

atlas <- atlas %>% mutate(annotation_source = ANNOTATION_SOURCE, atlas_url = ATLAS_URL)

# Attach the original reference descriptions without changing their wording.
d <- d %>% left_join(atlas, by = "module_id")

# Show categories alphabetically. This order does not imply biological importance.
category_order <- sort(unique(atlas$official_category))

module_summary <- module_summary %>% left_join(atlas, by = "module_id")

# Within each category, show modules with the largest activity changes first.
# Wrapping long descriptions changes their display only.
MODULE_LAYOUT <- module_summary %>% filter(selected_for_main) %>% mutate(category_order_index = match(official_category, 
  category_order)) %>% arrange(category_order_index, desc(max_abs_delta), module_id) %>% mutate(selection_order = row_number()) %>% 
  mutate(theme = official_category, label = ifelse(SHOW_FUNCTION_IN_MAIN & !is.na(official_function) & 
    nzchar(trimws(official_function)), paste0(module_id, " | ", stringr::str_wrap(official_function, 
    width = FUNCTION_WRAP)), module_id))

# Keep all five comparisons for each selected module, including nonsignificant changes.
source_main <- d %>% inner_join(MODULE_LAYOUT %>% dplyr::select(module_id, theme, label, selection_order), 
  by = "module_id")

if (nrow(source_main) != nrow(MODULE_LAYOUT) * nrow(CONTRASTS)) stop("Incomplete selected module/contrast grid.")

source_main <- source_main %>% mutate(label = factor(label, levels = rev(MODULE_LAYOUT$label)), 
  theme = factor(theme, levels = unique(MODULE_LAYOUT$theme)))

# 4. Save the plotted estimates, module order, reference descriptions and settings.
write_plot_table <- function(x, path) utils::write.table(as.data.frame(x), path, sep = "\t", quote = TRUE, 
  qmethod = "double", row.names = FALSE, na = "NA")

source_export <- source_main %>% dplyr::select(module_id, contrast_label, condition, x, delta_activity, 
  standard_error, ci95_low, ci95_high, p_value, q_primary, coverage_ok, primary_DiMA, statistical_small_effect, 
  official_category, official_function, selection_order)

write_plot_table(MODULE_LAYOUT, file.path(OUTPUT_DIR, "Fig6C_module_layout.tsv"))

write_plot_table(atlas, file.path(OUTPUT_DIR, "Fig6C_official_annotations_used.tsv"))

write_plot_table(source_export, file.path(OUTPUT_DIR, "Fig6C_source.tsv"))

parameters <- tibble(parameter = c("input", "q_cutoff", "abs_delta_cutoff", "display_selection", 
  "activity_scale", "annotation_source", "atlas_url", "category_order", "within_category_order", 
  "selected_module_count", "large_significant_change_count", "tested_comparison_count"), value = as.character(c(INPUT_FILE, 
  Q_CUTOFF, MIN_ABS_DELTA, "all modules with q < cutoff AND abs(delta) > cutoff AND coverage OK in at least one comparison; all five comparisons displayed", 
  "native reference M scale; no row standardization", ANNOTATION_SOURCE, ATLAS_URL, "alphabetical; categories taken verbatim from atlas", 
  "descending maximum absolute delta, then module ID", length(SELECTED_MODULES), sum(d$primary_DiMA), 
  nrow(d))))

write_plot_table(parameters, file.path(OUTPUT_DIR, "Fig6C_parameters.tsv"))

# 5. Draw activity changes on the original reference scale.
# Filled circles mark large significant changes; open circles mark smaller significant
# changes; crosses flag insufficient reference coverage. Colours keep the reference
# sign: a negative activity change need not mean that every module gene is downregulated.
make_plot <- function(dat, title, limit, faceted = TRUE) {
  p <- ggplot(dat, aes(x = x, y = label, fill = delta_activity)) + geom_tile(width = 0.94, height = 0.94, 
    colour = "white", linewidth = 0.3) + geom_point(data = dat %>% filter(primary_DiMA), shape = 16, 
    size = 1.9, colour = "black") + geom_point(data = dat %>% filter(statistical_small_effect), 
    shape = 1, size = 2.1, stroke = 0.55, colour = "black") + geom_point(data = dat %>% filter(!coverage_ok), 
    shape = 4, size = 2.1, colour = "black") + scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", 
    high = "#B2182B", midpoint = 0, limits = c(-limit, limit), name = expression(Delta * A)) + 
    scale_x_continuous(breaks = CONTRASTS$x, labels = CONTRASTS$axis_label, expand = expansion(add = 0.5)) + 
    labs(title = title, x = NULL, y = NULL, caption = paste0("Relative to induced control1. Filled circle: q < ", 
      Q_CUTOFF, " and |delta A| > ", MIN_ABS_DELTA, "; open circle: q < ", Q_CUTOFF, " with smaller effect. Cross: projection coverage flagged.\nNative atlas activity scale; module direction follows reference gene weights. No KD + TOB transcriptomes.")) + 
    theme_minimal(base_size = BASE_SIZE, base_family = FONT_FAMILY) + theme(panel.grid = element_blank(), 
    axis.text.y = element_text(colour = "black", size = BASE_SIZE), axis.text.x = element_text(colour = "black", 
      size = BASE_SIZE - 0.5), plot.title = element_text(face = "bold", size = 12), plot.caption = element_text(hjust = 0, 
      size = 8), strip.text.y = element_text(angle = 0, face = "bold", hjust = 0, size = 9), 
    strip.background = element_rect(fill = "grey95", colour = NA), strip.text.y.left = element_text(angle = 0, 
      face = "bold", hjust = 0, size = 8), legend.position = "right")
  if (faceted) 
    p <- p + facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y", labeller = labeller(theme = label_wrap_gen(width = CATEGORY_WRAP))) + 
      theme(strip.placement = "outside")
  p
}

# Use the same symmetric colour range across the complete main panel.
limit <- if (nrow(source_main)) max(1, ceiling(max(abs(source_main$delta_activity)))) else 1

p <- if (nrow(source_main)) {
  make_plot(source_main, "C  Regulatory-module responses by PRECISE411 functional category", 
    limit)
} else {
  ggplot() + annotate("text", x = 0, y = 0, label = "No modules pass the FDR, effect-size and coverage criteria.") + 
    theme_void()
}

# Give long category and module descriptions enough room; height follows module count.
MAIN_WIDTH <- 12

MAIN_HEIGHT <- max(7.5, nrow(MODULE_LAYOUT) * 0.32 + n_distinct(MODULE_LAYOUT$theme) * 0.22 + 1.8)

# 6. Save PNG and PDF. Use Cairo for the PDF when available to support special characters.
PDF_DEVICE <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

PDF_ARGS <- if (capabilities("cairo")) list() else list(useDingbats = FALSE)

grDevices::png(file.path(OUTPUT_DIR, "Fig6C_iModulon_biology.png"), width = MAIN_WIDTH, height = MAIN_HEIGHT, 
  units = "in", res = PNG_DPI, type = "cairo-png")

png_device <- grDevices::dev.cur()

tryCatch(print(p), finally = grDevices::dev.off(png_device))

do.call(ggsave, c(list(filename = file.path(OUTPUT_DIR, "Fig6C_iModulon_biology.pdf"), plot = p, 
  width = MAIN_WIDTH, height = MAIN_HEIGHT, units = "in", device = PDF_DEVICE), PDF_ARGS))

Fig6C_plot <- p

if (PRINT_IN_RSTUDIO) print(p)


# Print the locations of the saved figure-data tables.
figure_data_files <- file.path(OUTPUT_DIR,c("Fig6C_source.tsv","Fig6C_module_layout.tsv","Fig6C_official_annotations_used.tsv","Fig6C_parameters.tsv"))
message("Underlying figure data:\n",paste0("  ",normalizePath(figure_data_files,winslash="/"),collapse="\n"))
