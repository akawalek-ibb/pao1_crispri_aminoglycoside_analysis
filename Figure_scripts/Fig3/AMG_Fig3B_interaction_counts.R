# Plot gene–aminoglycoside interaction counts by essentiality category.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/Fig3/AMG_Fig3B_interaction_counts.R")
# Inputs:
# - preprocessing_scripts/out_03_AMG_relevant_results/AG_gene_results.tsv.gz,
#   AG_gene_recurrence.tsv, AG_thresholds.tsv, Comparison_manifest.tsv and Analysis_settings.tsv.
# Upstream / execution order:
# - preprocessing_scripts/03_AMG_relevant_results.R; existing fitted screen results can be reused.
# Outputs:
# - Figure_scripts/Fig3/out_AMG_Fig3B_interaction_counts/ — PNG/PDF, category counts, totals,
#   thresholds and settings TSVs.
# Interpretation / operation:
# - Reuse raw-Stouffer calls and comparison-specific thresholds; T2 is displayed separately from
#   T1 recurrence.
#
# R packages used/checked: data.table, ggplot2.

INPUT_DIR <- "preprocessing_scripts/out_03_AMG_relevant_results"

OUTPUT_DIR <- "Figure_scripts/Fig3/out_AMG_Fig3B_interaction_counts"

CONDITIONS_TO_PLOT <- NULL

CATEGORY_COLOURS <- c(E_core = "#8A0000", E_core_conditional = "#0070B0", E = "#FF6600", rest = "#A8A8A8")

CATEGORY_LABELS <- c(E_core = "core", E_core_conditional = "conditional", E = "effect on fitness\n(other studies)", 
  rest = "rest")

FONT_FAMILY <- "Arial"

WIDTH_IN <- 10

HEIGHT_IN <- 6.5

DPI <- 300

BAR_END_LABEL_SIZE <- 3.4

BAR_END_LABEL_OFFSET_FRACTION <- 0.025

SHOW_TITLE <- TRUE

SHOW_METHOD_SUBTITLE <- FALSE

PRINT_IN_RSTUDIO <- TRUE

OVERWRITE <- TRUE

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

s8_numeric <- function(x, columns, label) {
  for (name in columns) {
    original <- x[[name]]
    value <- suppressWarnings(as.numeric(original))
    if (any(!is.na(original) & !is.finite(value))) 
      stop(label, ": invalid numeric value in ", name)
    data.table::set(x, j = name, value = value)
  }
  invisible(x)
}

s8_call <- function(lfc, p, p_cut, lower, upper) {
  data.table::fcase(is.finite(p) & p < p_cut & is.finite(lfc) & lfc <= lower, "Sensitizing", 
    is.finite(p) & p < p_cut & is.finite(lfc) & lfc >= upper, "Protective", default = "Not significant")
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

s8_stage <- function(tables, directory, gzip = character()) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  for (name in names(tables)) {
    path <- file.path(directory, paste0(name, if (name %in% gzip) 
      ".tsv.gz"
    else ".tsv"))
    con <- if (endsWith(path, ".gz")) 
      gzfile(path, "wt")
    else file(path, "wt")
    tryCatch(utils::write.table(as.data.frame(tables[[name]]), con, sep = "\t", quote = TRUE, 
      qmethod = "double", row.names = FALSE, na = "NA"), finally = close(con))
  }
}

fig3_logical <- function(x, label) {
  value <- tolower(trimws(as.character(x)))
  if (anyNA(value) || any(!value %in% c("true", "false", "1", "0"))) 
    stop("Invalid logical values in ", label)
  value %in% c("true", "1")
}

fig3_category <- function(x) {
  x <- trimws(as.character(x))
  out <- rep("rest", length(x))
  out[x %in% c("E_core", "E_AMR_core")] <- "E_core"
  out[x %in% c("E_core_conditional", "E_AMR_core_conditional")] <- "E_core_conditional"
  out[x %in% c("E", "E_AMR", "AMR")] <- "E"
  known <- c("E_core", "E_AMR_core", "E_core_conditional", "E_AMR_core_conditional", "E", "E_AMR", 
    "AMR", "rest", "NA", "")
  if (any(!is.na(x) & !x %in% known)) 
    stop("Unsupported essentiality category: ", paste(unique(x[!is.na(x) & !x %in% known]), 
      collapse = ", "))
  out
}

# Read the task-03 effects, calls, recurrence, manifest and cutoffs together.
# Reconstruct decisions only to check consistency with the stored calls/counts;
# do not estimate new cutoffs or replace the primary screening decisions.
fig3_inputs <- function(directory) {
  s8_require("data.table")
  paths <- file.path(directory, c("AG_gene_results.tsv.gz", "AG_gene_recurrence.tsv", "Comparison_manifest.tsv", 
    "AG_thresholds.tsv", "Analysis_settings.tsv"))
  if (any(!file.exists(paths))) 
    stop("Run 03_AMG_relevant_results.R. Missing:\n", paste(paths[!file.exists(paths)], collapse = "\n"))
  genes <- s8_read(paths[1L])
  recurrent <- s8_read(paths[2L])
  manifest <- s8_read(paths[3L])
  thresholds <- s8_read(paths[4L])
  s8_required(genes, c("PAO1_ID", "comparison_id", "gene_name", "product_name", "essential_category_code", 
    "medLFC", "p_stouffer", "stouffer_call", "stouffer_evaluated"), paths[1L])
  s8_required(recurrent, c("PAO1_ID", "gene_name", "product_name", "essential_category_code", 
    "n_independent_T1_sensitizing", "n_independent_T1_protective", "recurrent_sensitizing", 
    "recurrent_protective"), paths[2L])
  s8_required(manifest, c("comparison_id", "display_label", "plot_order", "count_for_recurrence"), 
    paths[3L])
  s8_required(thresholds, c("comparison_id", "p_stouffer_cutoff", "sensitizing_LFC_cutoff_used", 
    "protective_LFC_cutoff_used"), paths[4L])
  s8_unique(genes, c("PAO1_ID", "comparison_id"), paths[1L])
  s8_unique(recurrent, "PAO1_ID", paths[2L])
  s8_unique(manifest, "comparison_id", paths[3L])
  s8_unique(thresholds, "comparison_id", paths[4L])
  if (anyDuplicated(manifest$display_label) || anyDuplicated(manifest$plot_order)) 
    stop("Ambiguous comparison labels/order in task 03.")
  manifest[, `:=`(count_for_recurrence, fig3_logical(count_for_recurrence, "count_for_recurrence"))]
  for (column in c("recurrent_sensitizing", "recurrent_protective")) data.table::set(recurrent, 
    j = column, value = fig3_logical(recurrent[[column]], column))
  genes[, `:=`(stouffer_evaluated, fig3_logical(stouffer_evaluated, "stouffer_evaluated"))]
  s8_numeric(genes, c("medLFC", "p_stouffer"), paths[1L])
  if (any(!genes$stouffer_call %in% c("Sensitizing", "Protective", "Not significant"))) 
    stop("Unknown task-03 interaction call.")
  if (length(setdiff(manifest$comparison_id, genes$comparison_id)) || length(setdiff(genes$comparison_id, 
    manifest$comparison_id)) || !setequal(manifest$comparison_id, thresholds$comparison_id)) 
    stop("Task-03 comparison products disagree.")
  if (length(setdiff(genes$PAO1_ID, recurrent$PAO1_ID))) 
    stop("Task-03 gene and recurrence universes disagree.")
  p_cut <- s8_setting(directory, "P_STOUFFER_CUT", numeric = TRUE)
  min_recurrence <- s8_setting(directory, "RECURRENCE_MIN", numeric = TRUE)
  if (!is.finite(p_cut) || !is.finite(min_recurrence)) 
    stop("Invalid task-03 settings.")
  ti <- match(genes$comparison_id, thresholds$comparison_id)
  expected <- s8_call(genes$medLFC, genes$p_stouffer, p_cut, thresholds$sensitizing_LFC_cutoff_used[ti], 
    thresholds$protective_LFC_cutoff_used[ti])
  if (any(expected != genes$stouffer_call) || any(thresholds$p_stouffer_cutoff != p_cut)) 
    stop("Stored calls/thresholds disagree. Regenerate task 03 as one run.")
  if (any(genes$stouffer_evaluated != (is.finite(genes$medLFC) & is.finite(genes$p_stouffer)))) 
    stop("Stored evaluation flags disagree with finite effect/p values.")
  primary <- genes[comparison_id %in% manifest$comparison_id[manifest$count_for_recurrence]]
  check <- primary[, .(n_sens = sum(stouffer_call == "Sensitizing"), n_prot = sum(stouffer_call == "Protective")), by = PAO1_ID]
  ix <- match(recurrent$PAO1_ID, check$PAO1_ID)
  ns <- check$n_sens[ix]
  np <- check$n_prot[ix]
  ns[is.na(ns)] <- 0L
  np[is.na(np)] <- 0L
  if (any(recurrent$n_independent_T1_sensitizing != ns) || any(recurrent$n_independent_T1_protective != 
    np) || any(recurrent$recurrent_sensitizing != (ns >= min_recurrence)) || any(recurrent$recurrent_protective != 
    (np >= min_recurrence))) 
    stop("Stored recurrence is inconsistent with the task-03 T1 calls.")
  data.table::setorder(manifest, plot_order)
  genes[, `:=`(category, fig3_category(essential_category_code))]
  recurrent[, `:=`(category, fig3_category(essential_category_code))]
  list(genes = genes, recurrent = recurrent, manifest = manifest, thresholds = thresholds, source_paths = paths, 
    p_cut = p_cut, recurrence_min = min_recurrence)
}

fig3_conditions <- function(manifest, labels = NULL, primary_only = FALSE) {
  if (is.null(labels)) 
    return(if (primary_only) manifest[count_for_recurrence == TRUE] else data.table::copy(manifest))
  if (!length(labels) || anyDuplicated(labels) || length(setdiff(labels, manifest$display_label))) 
    stop("Unknown, empty or duplicate comparison selection.")
  manifest[match(labels, display_label)]
}

# Read, prepare, plot and export in order.
s8_require(c("data.table", "ggplot2"))

if (!file.exists("AMG_cleaned_Git.Rproj")) stop("Open AMG_cleaned_Git.Rproj and run from its working directory.")

# 2. Read final calls without recomputing significance
input <- fig3_inputs(INPUT_DIR)

manifest <- fig3_conditions(input$manifest, CONDITIONS_TO_PLOT)

if (!identical(names(CATEGORY_LABELS), names(CATEGORY_COLOURS))) stop("Category labels and colours must have the same named order.")

z <- input$genes[comparison_id %in% manifest$comparison_id]

# 3. Count calls by direction and essentiality
# Count the stored interaction calls within each comparison and essentiality
# category. T2 can be displayed, but it is not added to T1 recurrence. These
# bars count gene-comparison calls, not unique recurrent genes across all bars.
observed <- z[stouffer_call != "Not significant", .(n_genes = .N), by = .(comparison_id, category, direction = stouffer_call)]

# Complete the category/direction grid with zeros for empty combinations.
# Use negative plotting coordinates for sensitizing counts and positive ones
# for protective counts; this sign convention positions the bars only.
counts <- data.table::CJ(comparison_id = manifest$comparison_id, category = names(CATEGORY_COLOURS), 
  direction = c("Sensitizing", "Protective"), unique = TRUE)

counts <- merge(counts, observed, by = c("comparison_id", "category", "direction"), all.x = TRUE, 
  sort = FALSE)

counts[is.na(n_genes), `:=`(n_genes, 0L)]

counts[, `:=`(display_label = manifest$display_label[match(comparison_id, manifest$comparison_id)], 
  plot_order = match(comparison_id, manifest$comparison_id), count_for_recurrence = manifest$count_for_recurrence[match(comparison_id, manifest$comparison_id)], signed_count = ifelse(direction == "Sensitizing", -n_genes, n_genes))]

data.table::setorder(counts, plot_order, direction, category)

totals <- counts[, .(n_genes = sum(n_genes)), by = .(comparison_id, display_label, plot_order, count_for_recurrence, direction)]

totals[, `:=`(signed_count, ifelse(direction == "Sensitizing", -n_genes, n_genes))]

gap <- max(1, BAR_END_LABEL_OFFSET_FRACTION * max(totals$n_genes))

totals[, `:=`(label_position = signed_count + ifelse(direction == "Sensitizing", -gap, gap), label = paste0("n=", n_genes))]

plot_data <- data.table::copy(counts)

plot_totals <- data.table::copy(totals)

plot_data[, `:=`(display_label = factor(display_label, levels = rev(manifest$display_label)), category = factor(category, levels = names(CATEGORY_COLOURS)))]

plot_totals[, `:=`(display_label, factor(display_label, levels = rev(manifest$display_label)))]

subtitle <- if (SHOW_METHOD_SUBTITLE) paste0("Task-03 calls; raw Stouffer p < ", input$p_cut, "; comparison-specific effect cutoffs") else NULL

# 4. Draw the stacked counts
plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = display_label, y = signed_count, fill = category)) + 
  ggplot2::geom_col(position = "stack") +
  ggplot2::geom_text(data = plot_totals[direction == "Sensitizing"], ggplot2::aes(x = display_label, y = label_position, label = label), inherit.aes = FALSE, 
  hjust = 1, size = BAR_END_LABEL_SIZE, fontface = "bold", colour = "black", family = FONT_FAMILY) + 
  ggplot2::geom_text(data = plot_totals[direction == "Protective"], ggplot2::aes(x = display_label, 
    y = label_position, label = label), inherit.aes = FALSE, hjust = 0, size = BAR_END_LABEL_SIZE, 
    fontface = "bold", colour = "black", family = FONT_FAMILY) +
  ggplot2::geom_hline(yintercept = 0, 
  colour = "black", linetype = "dashed", linewidth = 1) +
  ggplot2::coord_flip(clip = "off") + 
  ggplot2::scale_y_continuous(labels = abs, expand = ggplot2::expansion(mult = c(0.1, 0.1))) + 
  ggplot2::scale_fill_manual(values = CATEGORY_COLOURS, breaks = names(CATEGORY_COLOURS), labels = CATEGORY_LABELS, 
    drop = FALSE) +
  ggplot2::labs(title = if (SHOW_TITLE) "Significant CRISPRi aminoglycoside interactions" else NULL, 
  subtitle = subtitle, x = "Comparison", y = "Number of genes", fill = "Essential category") + 
  ggplot2::theme_minimal(base_size = 12, base_family = FONT_FAMILY) +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 10), 
  legend.position = "right", plot.margin = ggplot2::margin(10, 24, 10, 24))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

stage <- OUTPUT_DIR

ggplot2::ggsave(file.path(stage, "Fig3B_interaction_counts.png"), plot, width = WIDTH_IN, height = HEIGHT_IN, 
  dpi = DPI, bg = "white")

ggplot2::ggsave(file.path(stage, "Fig3B_interaction_counts.pdf"), plot, device = grDevices::cairo_pdf, 
  width = WIDTH_IN, height = HEIGHT_IN, bg = "white")

settings <- s8_settings(c(P_STOUFFER_CUT = input$p_cut, MAD_MULTIPLIER = s8_setting(INPUT_DIR, 
  "MAD_MULTIPLIER"), MIN_ABS_MEDLFC = s8_setting(INPUT_DIR, "MIN_ABS_MEDLFC"), displayed_comparisons = paste(manifest$display_label, 
  collapse = "; "), calls = "Stored task-03 Stouffer decisions; no BH filter"), stats::setNames(input$source_paths, 
  basename(input$source_paths)))

# Save the plotted counts and totals beside the inherited effect thresholds
# and settings, allowing the displayed bars to be checked against task 03.
s8_stage(list(Fig3B_counts_by_category = counts, Fig3B_direction_totals = totals, Fig3B_thresholds = input$thresholds[comparison_id %in% manifest$comparison_id], Analysis_settings = settings), stage)

if (PRINT_IN_RSTUDIO) print(plot)

message("Saved Fig. 3B: ", nrow(manifest), " comparisons. ", OUTPUT_DIR)

invisible(list(plot = plot, counts = counts, totals = totals))


figure_data_files <- file.path(OUTPUT_DIR,c("Fig3B_counts_by_category.tsv","Fig3B_direction_totals.tsv","Fig3B_thresholds.tsv","Analysis_settings.tsv"))
message("Underlying figure data:\n",paste0("  ",normalizePath(figure_data_files,winslash="/"),collapse="\n"))
