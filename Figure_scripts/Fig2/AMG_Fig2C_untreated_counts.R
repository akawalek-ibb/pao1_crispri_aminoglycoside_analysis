# Count untreated negative/positive fitness hits across five comparisons.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/Fig2/AMG_Fig2C_untreated_counts.R")
# Inputs:
# - preprocessing_scripts/out_12_AMG_MH_fitness/MH_fitness.tsv: one row per PAO1 locus, category
#   and B1 T1/T2 plus B2–B4 T1 effects/p-values.
# Upstream / execution order:
# - preprocessing_scripts/11_AMG_gene_vulnerability.R →
#   preprocessing_scripts/12_AMG_MH_fitness.R; use existing edgeR_data fits or rebuild them with
#   preprocessing_scripts/02_AMG_metacounter.R.
# Outputs:
# - Figure_scripts/Fig2/out_AMG_Fig2C_untreated_counts/ — PNG/PDF plus category counts, totals
#   and diagnostics TSVs.
# Interpretation / operation:
# - Calls use raw Stouffer p < 0.05 and median log2FC ≤ −1 or ≥ 1.
#
# R packages used/checked: R.utils, data.table, dplyr, ggplot2, openxlsx, tibble, tidyr.

OUTPUT_DIR <- "Figure_scripts/Fig2/out_AMG_Fig2C_untreated_counts"

DEFAULT_FIG4_CONDITIONS <- c("TOB 0.025", "TOB 0.1", "TOB 0.4", "TOB 1", "AMK")

MH_FITNESS_FILE <- "preprocessing_scripts/out_12_AMG_MH_fitness/MH_fitness.tsv"

P_STOUFFER_CUT <- 0.05

LFC_CUT <- 1

FONT_FAMILY <- "Arial"

DPI <- 300

PRINT_IN_RSTUDIO <- TRUE

PLOT_WIDTH <- 9

PLOT_HEIGHT <- 5.8

CATEGORY_LEVELS <- c("E_core", "E_core_conditional", "E", "rest")

CATEGORY_COLORS <- c(E_core = "#8A0000", E_core_conditional = "#0070B0", E = "#FF6600", rest = "#A8A8A8")

CATEGORY_LABELS <- c(E_core = "core", E_core_conditional = "conditional", E = "effect on fitness\n(other studies)", 
  rest = "rest")

DIRECTION_LEVELS <- c("Negative", "Positive")

# HELPERS, INPUTS AND ANALYSIS
# Read, prepare, plot and export in order.
packages <- c("R.utils", "data.table", "dplyr", "ggplot2", "openxlsx", "tibble", "tidyr")

missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(tidyr)
})

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

combined_gene_file <- MH_FITNESS_FILE

# Convert the prepared wide gene table into one row per gene/comparison.
# B1 T1 and T2 are separate time points, not independent batches; each panel
# counts its own gene effects rather than the union across experiments.
read_relevant_gene_results <- function() {
  fitness <- data.table::fread(MH_FITNESS_FILE, data.table = FALSE)
  prefixes <- c("B1_T1_vs_T0", "B1_T2_vs_T0", "B2_T1_vs_T0", "B3_T1_vs_T0", "B4_T1_vs_T0")
  comparisons <- vector("list", 5)
  for (i in 1:5) comparisons[[i]] <- data.frame(PAO1_ID = fitness$PAO1_ID, category = fitness$category, 
    batch = c("b1", "b1", "b2", "b3", "b4")[i], condition = c("None_0_T1 - None_0_T0", "None_0_T2 - None_0_T0", 
      rep("None_0_T1 - None_0_T0", 3))[i], medLFC = fitness[[paste0(prefixes[i], "_medLFC")]], 
    p_stouffer = fitness[[paste0(prefixes[i], "_p_stouffer")]])
  dplyr::bind_rows(comparisons)
}

# 2. Read and reshape the five untreated comparisons
gene_results <- read_relevant_gene_results()

required_columns <- c("batch", "condition", "category", "medLFC", "p_stouffer")

missing_columns <- setdiff(required_columns, names(gene_results))

if (length(missing_columns) > 0L) {
  stop("Gene-level metacounter table is missing columns: ", paste(missing_columns, collapse = ", "))
}

gene_results <- gene_results %>% dplyr::mutate(batch = tolower(trimws(as.character(batch))), condition = trimws(as.character(condition)), 
  category = trimws(as.character(category)), medLFC = suppressWarnings(as.numeric(medLFC)), p_stouffer = suppressWarnings(as.numeric(p_stouffer)))

comparison_manifest <- tibble::tribble(~batch, ~condition, ~plot_label, ~plot_order, "b1", "None_0_T1 - None_0_T0", 
  "None B1: T1 vs T0", 1L, "b1", "None_0_T2 - None_0_T0", "None B1: T2 vs T0", 2L, "b2", "None_0_T1 - None_0_T0", 
  "None B2: T1 vs T0", 3L, "b3", "None_0_T1 - None_0_T0", "None B3: T1 vs T0", 4L, "b4", "None_0_T1 - None_0_T0", 
  "None B4: T1 vs T0", 5L)

available_comparisons <- gene_results %>% dplyr::distinct(batch, condition)

missing_comparisons <- comparison_manifest %>% dplyr::anti_join(available_comparisons, by = c("batch", 
  "condition"))

if (nrow(missing_comparisons) > 0L) {
  stop("The following requested None comparisons are absent:\n", paste(paste0(missing_comparisons$batch, 
    ": ", missing_comparisons$condition), collapse = "\n"))
}

recode_category <- function(x) {
  dplyr::case_when(x %in% c("E_AMR", "AMR") ~ "E", x == "E_AMR_core" ~ "E_core", x == "E_AMR_core_conditional" ~ 
    "E_core_conditional", is.na(x) | x == "" | x == "NA" ~ "rest", TRUE ~ x)
}

# 3. Count effects by comparison and category
# Apply the configured raw-p-value and absolute-effect criteria to each
# comparison. Assign direction from medLFC and group the called genes by their
# existing essentiality annotation. No guide p-values are combined here.
observed_counts <- lapply(seq_len(nrow(comparison_manifest)), function(i) {
  cfg <- comparison_manifest[i, ]
  gene_results %>% dplyr::filter(batch == cfg$batch, condition == cfg$condition, is.finite(medLFC), 
    is.finite(p_stouffer), p_stouffer < P_STOUFFER_CUT, abs(medLFC) >= LFC_CUT) %>% dplyr::mutate(category = recode_category(category), 
    change = dplyr::if_else(medLFC > 0, "Positive", "Negative")) %>% dplyr::count(category, 
    change, name = "count") %>% dplyr::mutate(batch = cfg$batch, condition = cfg$condition, 
    comparison = cfg$plot_label, plot_order = cfg$plot_order)
}) %>% dplyr::bind_rows()

# Include category/direction combinations with zero hits so their absence
# does not change the bar layout. Negative signed counts only place depletion
# bars on the left; the underlying gene counts remain nonnegative.
expected_grid <- comparison_manifest %>% dplyr::transmute(batch, condition, comparison = plot_label, 
  plot_order) %>% tidyr::crossing(category = CATEGORY_LEVELS, change = DIRECTION_LEVELS)

all_counts <- expected_grid %>% dplyr::left_join(observed_counts, by = c("batch", "condition", 
  "comparison", "plot_order", "category", "change")) %>% dplyr::mutate(count = tidyr::replace_na(count, 
  0L), category = factor(category, levels = CATEGORY_LEVELS), change = factor(change, levels = DIRECTION_LEVELS), 
  signed_count = dplyr::if_else(change == "Negative", -as.numeric(count), as.numeric(count)))

direction_totals <- all_counts %>% dplyr::group_by(comparison, batch, condition, plot_order, change) %>% 
  dplyr::summarise(n_genes = sum(count, na.rm = TRUE), .groups = "drop")

totals <- all_counts %>% dplyr::group_by(comparison, plot_order) %>% dplyr::summarise(total = sum(count, 
  na.rm = TRUE), .groups = "drop") %>% dplyr::arrange(plot_order)

# Record how many genes were tested, passed the p-value criterion, and passed
# both criteria. These diagnostics explain why plotted counts can differ
# without changing the input gene-level estimates.
count_diagnostics <- lapply(seq_len(nrow(comparison_manifest)), function(i) {
  cfg <- comparison_manifest[i, ]
  x <- gene_results %>% dplyr::filter(batch == cfg$batch, condition == cfg$condition)
  tibble::tibble(batch = cfg$batch, condition = cfg$condition, plot_label = cfg$plot_label, plot_order = cfg$plot_order, 
    n_genes_tested = sum(is.finite(x$p_stouffer) & is.finite(x$medLFC)), n_p_below_cut = sum(is.finite(x$p_stouffer) & 
      x$p_stouffer < P_STOUFFER_CUT, na.rm = TRUE), n_p_and_absLFC_cut = sum(is.finite(x$p_stouffer) & 
      x$p_stouffer < P_STOUFFER_CUT & is.finite(x$medLFC) & abs(x$medLFC) >= LFC_CUT, na.rm = TRUE))
}) %>% dplyr::bind_rows()

comparison_order <- comparison_manifest %>% dplyr::arrange(plot_order) %>% dplyr::pull(plot_label)

all_counts <- all_counts %>% dplyr::mutate(comparison = factor(comparison, levels = rev(comparison_order)))

totals <- totals %>% dplyr::mutate(comparison = factor(comparison, levels = rev(comparison_order)), 
  comparison_character = as.character(comparison), axis_label = paste0(comparison_character, 
    " (n=", total, ")"))

axis_label_lookup <- stats::setNames(totals$axis_label, totals$comparison_character)

# 4. Draw the diverging barplot
p_counts <- ggplot2::ggplot(all_counts, ggplot2::aes(x = comparison, y = signed_count, fill = category)) + 
  ggplot2::geom_col(position = "stack") +
  ggplot2::coord_flip() +
  ggplot2::geom_hline(yintercept = 0, 
  color = "black", linetype = "dashed", linewidth = 1) +
  ggplot2::theme_minimal(base_size = 12, 
  base_family = FONT_FAMILY) +
  ggplot2::labs(title = "CRISPRi fitness effects in untreated control populations", 
  subtitle = "Raw p_stouffer < 0.05 and |medLFC| >= 1", x = "Comparison", y = "Number of genes", 
  fill = "Essential category") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10), 
  axis.text.y = ggplot2::element_text(size = 10), legend.position = "right") +
  ggplot2::scale_y_continuous(labels = abs) + 
  ggplot2::scale_x_discrete(labels = axis_label_lookup) +
  ggplot2::scale_fill_manual(values = CATEGORY_COLORS, 
  breaks = CATEGORY_LEVELS, labels = CATEGORY_LABELS, drop = FALSE)

# Export the exact category counts, totals and diagnostic denominators used
# in the plot. A gene can contribute in more than one comparison, so sums
# across bars are not the number of unique fitness-associated genes.
data.table::fwrite(data.table::as.data.table(all_counts), file.path(OUTPUT_DIR, "None_significant_hit_counts_by_category.tsv"), 
  sep = "\t")

data.table::fwrite(data.table::as.data.table(totals), file.path(OUTPUT_DIR, "None_significant_hit_totals.tsv"), 
  sep = "\t")

data.table::fwrite(data.table::as.data.table(direction_totals), file.path(OUTPUT_DIR, "None_significant_hit_direction_totals.tsv"), 
  sep = "\t")

data.table::fwrite(data.table::as.data.table(count_diagnostics), file.path(OUTPUT_DIR, "None_significant_hit_diagnostics.tsv"), 
  sep = "\t")

settings <- tibble::tibble(setting = c("input", "P_STOUFFER_CUT", "LFC_CUT", "comparisons"), value = c(combined_gene_file, 
  as.character(P_STOUFFER_CUT), as.character(LFC_CUT), paste(comparison_manifest$plot_label, 
    collapse = "; ")))

ggplot2::ggsave(filename = file.path(OUTPUT_DIR, "None_significant_hit_counts.png"), plot = p_counts, 
  width = PLOT_WIDTH, height = PLOT_HEIGHT, units = "in", dpi = DPI, bg = "white")

ggplot2::ggsave(filename = file.path(OUTPUT_DIR, "None_significant_hit_counts.pdf"), plot = p_counts, 
  width = PLOT_WIDTH, height = PLOT_HEIGHT, units = "in", device = grDevices::cairo_pdf, bg = "white")

if (isTRUE(PRINT_IN_RSTUDIO)) print(p_counts)

message("Completed: AMG_Fig2C_untreated_counts")

data_paths <- normalizePath(file.path(OUTPUT_DIR, c("None_significant_hit_counts_by_category.tsv", 
  "None_significant_hit_totals.tsv", "None_significant_hit_direction_totals.tsv", "None_significant_hit_diagnostics.tsv")), 
  winslash = "/", mustWork = TRUE)

message("Underlying figure data (Fig. 2C):\n", paste0("  ", data_paths, collapse = "\n"))

Fig2C_plot <- p_counts

