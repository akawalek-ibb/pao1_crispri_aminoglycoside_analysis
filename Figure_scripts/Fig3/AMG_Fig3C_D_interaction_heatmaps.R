# Plot interaction-severity heatmaps for recurrent sensitizing/protective genes.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/Fig3/AMG_Fig3C_D_interaction_heatmaps.R")
# Inputs:
# - preprocessing_scripts/out_03_AMG_relevant_results/ — effects, calls, recurrence and
#   settings; preprocessing_scripts/out_12_AMG_MH_fitness/MH_fitness.tsv.
# - Figure_scripts/Fig3/input_AMG_Fig3/manual_functional_annotations.tsv: PAO1_ID, gene_name,
#   annotation_category; matching palette TSV.
# - annotation_input/pathoadaptive_mutations_science.adi0908_table_s5.csv; optional KEGG
#   annotation is configured below.
# Upstream / execution order:
# - preprocessing_scripts/03_AMG_relevant_results.R;
#   preprocessing_scripts/11_AMG_gene_vulnerability.R →
#   preprocessing_scripts/12_AMG_MH_fitness.R.
# Outputs:
# - Figure_scripts/Fig3/out_AMG_Fig3C_D_interaction_heatmaps/ — PNG/PDF heatmaps and colour
#   legend plus plotted data/settings TSVs.
# Interpretation / operation:
# - Manual is the default annotation mode; KEGG and None are alternatives.
# - Cluster categorical calls using average-linkage Hamming distance; significant cells display
#   capped median adjusted log2FC.
# - No sensitivity analysis, enrichment test or guide-response fit is run.
#
# R packages used/checked: data.table, pheatmap.

INPUT_DIR <- "preprocessing_scripts/out_03_AMG_relevant_results"

OUR_FITNESS_FILE <- "preprocessing_scripts/out_12_AMG_MH_fitness/MH_fitness.tsv"

OUTPUT_DIR <- "Figure_scripts/Fig3/out_AMG_Fig3C_D_interaction_heatmaps"

CONDITIONS_TO_PLOT <- NULL

DIRECTIONS_TO_PLOT <- c("Sensitizing", "Protective")

LFC_DISPLAY_CAP <- 3

CLUSTER_METHOD <- "average"

SHOW_GENE_TREE <- TRUE

SHOW_PUBLISHED_ANNOTATION <- TRUE

SHOW_OUR_FITNESS_ANNOTATION <- TRUE

SHOW_PATHOADAPTIVE_ANNOTATION <- TRUE

PATHOADAPTIVE_FILE <- "annotation_input/pathoadaptive_mutations_science.adi0908_table_s5.csv"

PATHOADAPTIVE_FIRST_N <- 224L

FUNCTIONAL_ANNOTATION_MODE <- "Manual"

MANUAL_ANNOTATION_FILE <- "Figure_scripts/Fig3/input_AMG_Fig3/manual_functional_annotations.tsv"

MANUAL_MISSING_CATEGORY <- "error"

USE_MANUAL_GENE_NAMES <- TRUE

KEGG_FILE <- "annotation_input/PAO1_pathways_all_withCOG.tsv"

KEGG_PRODUCT_FALLBACK <- TRUE

FUNCTIONAL_PALETTE_FILE <- "Figure_scripts/Fig3/input_AMG_Fig3/functional_annotation_palette.tsv"

CUSTOM_FUNCTIONAL_COLOURS <- character()

FONT_FAMILY <- "Arial"

WIDTH_IN <- 18

HEIGHT_IN <- 4.5

CELL_HEIGHT_PT <- 20

TREE_HEIGHT_PT <- 20

FONT_SIZE_GENE <- 7

FONT_SIZE_CONDITION <- 14

FONT_SIZE_ANNOTATION <- 14

FONT_SIZE_ANNOTATION_DETAIL <- 9.5

ANNOTATION_TRACK_HEIGHT_PT <- 20

ANNOTATION_TRACK_GAP_PT <- 8

ANNOTATION_TOP_GAP_PT <- 6

ANNOTATION_LABELS <- list(`Published essentiality` = c("Essential", "(other)"), `Our CRISPRi fitness` = c("Essential", 
  "(this study)"), `Pathoadaptive status` = c("Pathoadaptive", "(Weimann et al., 2024)"), `Functional class` = "Functional class")

DPI <- 300

SHOW_TITLES <- FALSE

SHOW_ANNOTATION_LEGENDS <- FALSE

SHOW_VALUE_LEGEND <- TRUE

EXPORT_COLOUR_LEGEND <- TRUE

LEGEND_WIDTH_IN <- 10

LEGEND_FONT_SIZE <- 9.5

PRINT_IN_RSTUDIO <- TRUE

OVERWRITE <- TRUE

# DISPLAY AND EXPORT HELPERS
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

fig3_label <- function(name, id) {
  name <- trimws(gsub("'", "", as.character(name), fixed = TRUE))
  name[is.na(name) | !nzchar(name)] <- id[is.na(name) | !nzchar(name)]
  name
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

# Measure how often two genes have different categorical interaction calls
# across conditions. Cluster these call patterns, not the numerical effect
# sizes used for colour; the configured linkage then determines gene order.
fig3_hamming <- function(mat) {
  n <- nrow(mat)
  distance <- matrix(0, n, n, dimnames = list(rownames(mat), rownames(mat)))
  if (n > 1L) 
    for (i in seq_len(n - 1L)) for (j in seq.int(i + 1L, n)) {
      comparable <- !is.na(mat[i, ]) & !is.na(mat[j, ])
      distance[i, j] <- distance[j, i] <- if (any(comparable)) 
        mean(mat[i, comparable] != mat[j, comparable])
      else 1
    }
  stats::as.dist(distance)
}

# Select genes recurrent in the requested direction, then keep their stored
# effects/calls across the displayed comparisons. A gene recurrently sensitizing
# can still have a protective call in another comparison; that cell is retained.
# Keep numerical effects, categorical calls and evaluation flags separately.
fig3_heatmap_data <- function(input, direction, manifest, cap = 3, cluster_method = "average") {
  if (!direction %in% c("Sensitizing", "Protective")) 
    stop("Unknown interaction direction.")
  flag <- paste0("recurrent_", tolower(direction))
  count <- paste0("n_independent_T1_", tolower(direction))
  meta <- data.table::copy(input$recurrent[input$recurrent[[flag]]])
  if (!nrow(meta)) 
    stop("No recurrent ", tolower(direction), " genes in task 03.")
  z <- input$genes[PAO1_ID %in% meta$PAO1_ID & comparison_id %in% manifest$comparison_id]
  strength <- z[, .(max_abs_medLFC = if (any(is.finite(medLFC))) 
    max(abs(medLFC[is.finite(medLFC)]))
  else -Inf), by = PAO1_ID]
  meta[, `:=`(n_significant_conditions = get(count), max_abs_medLFC = strength$max_abs_medLFC[match(PAO1_ID, strength$PAO1_ID)])]
  meta[is.na(max_abs_medLFC), `:=`(max_abs_medLFC, -Inf)]
  data.table::setorder(meta, -n_significant_conditions, -max_abs_medLFC, PAO1_ID)
  meta[, `:=`(priority_order, seq_len(.N))]
  effects <- matrix(NA_real_, nrow(meta), nrow(manifest), dimnames = list(meta$PAO1_ID, manifest$display_label))
  calls <- matrix("Not significant", nrow(meta), nrow(manifest), dimnames = dimnames(effects))
  evaluated <- matrix(FALSE, nrow(meta), nrow(manifest), dimnames = dimnames(effects))
  ij <- cbind(match(z$PAO1_ID, meta$PAO1_ID), match(z$comparison_id, manifest$comparison_id))
  effects[ij] <- z$medLFC
  calls[ij] <- z$stouffer_call
  evaluated[ij] <- z$stouffer_evaluated
  # Mask unevaluated/nonsignificant cells rather than assigning them zero effect.
  # Cap displayed colours only; raw medLFC remains available in the source table.
  # For clustering, an absent comparison retains the default Not significant state;
  # the exported evaluated flag distinguishes absence from a tested no-call.
  significant <- evaluated & calls != "Not significant"
  shown <- pmax(pmin(effects, cap), -cap)
  shown[!significant] <- NA_real_
  tree <- if (nrow(meta) > 1L) 
    stats::hclust(fig3_hamming(calls), method = cluster_method)
  else NULL
  order <- if (is.null(tree)) 
    seq_len(nrow(meta))
  else tree$order
  meta[, `:=`(display_order, match(seq_len(.N), order))]
  list(metadata = meta, effects = effects, calls = calls, evaluated = evaluated, significant = significant, 
    shown = shown, tree = tree, order = order)
}

fig3_manual <- function(path) {
  x <- s8_read(path)
  aliases <- c("annotation category", "category", "functional_class")
  if (!"annotation_category" %in% names(x)) {
    found <- intersect(aliases, names(x))
    if (length(found) == 1L) 
      data.table::setnames(x, found, "annotation_category")
  }
  s8_required(x, c("PAO1_ID", "gene_name", "annotation_category"), path)
  x[, `:=`(PAO1_ID = trimws(as.character(PAO1_ID)), gene_name = trimws(as.character(gene_name)), annotation_category = trimws(as.character(annotation_category)))]
  s8_unique(x, "PAO1_ID", path)
  if (any(!nzchar(x$PAO1_ID))) 
    stop("Empty PAO1_ID in manual annotation.")
  x[, .(PAO1_ID, gene_name, annotation_category)]
}

# Prepare the functional annotation track without altering the selected genes.
# Manual reads user categories by PAO1_ID; KEGG assigns broad classes through
# the text-matching rules below; None omits this functional track. Other enabled
# annotation tracks remain independent. This is annotation, not enrichment.
fig3_functional <- function(meta, mode, manual_file, kegg_file, missing = "error", product_fallback = TRUE) {
  mode <- tolower(mode)
  if (!mode %in% c("manual", "kegg", "none") || !missing %in% c("error", "unclassified")) 
    stop("Invalid functional annotation mode/missing policy.")
  out <- data.table::data.table(PAO1_ID = meta$PAO1_ID, manual_gene_name = NA_character_, annotation_category = NA_character_, 
    annotation_present = FALSE)
  if (mode == "none") 
    return(out)
  if (mode == "manual") {
    x <- fig3_manual(manual_file)
    ix <- match(meta$PAO1_ID, x$PAO1_ID)
    out[, `:=`(manual_gene_name = x$gene_name[ix], annotation_category = x$annotation_category[ix])]
    out[, `:=`(annotation_present, !is.na(annotation_category) & nzchar(annotation_category))]
    absent <- out$PAO1_ID[!out$annotation_present]
    if (length(absent) && missing == "error") 
      stop("Manual annotation is missing categories for: ", paste(absent, collapse = ", "))
    if (length(absent)) {
      warning("Missing manual assignments use Other / unclassified: ", paste(absent, collapse = ", "))
      out[annotation_present == FALSE, `:=`(annotation_category, "Other / unclassified")]
    }
    return(out)
  }
  raw <- s8_read(kegg_file)
  normalized <- tolower(gsub("[^A-Za-z0-9]+", "_", names(raw)))
  data.table::setnames(raw, names(raw), normalized)
  s8_required(raw, c("locus_tag", "pathway_name", "pathway_database"), kegg_file)
  raw <- raw[toupper(trimws(pathway_database)) == "KEGG" & locus_tag %in% meta$PAO1_ID]
  terms <- raw[, .(terms = paste(sort(unique(pathway_name[!is.na(pathway_name)])), collapse = "; ")), 
    by = locus_tag]
  text <- terms$terms[match(meta$PAO1_ID, terms$locus_tag)]
  text[is.na(text)] <- ""
  if (product_fallback) 
    text <- paste(text, ifelse(is.na(meta$product_name), "", meta$product_name), sep = "; ")
  classify <- function(value) {
    value <- tolower(value)
    rules <- c(Lipopolysaccharide = "lipopolysaccharide|lipid a|kdo|octulosonate|lipopolysaccharide transport", 
      `Protein export / quality control` = "protein export|secretion system|secretion protein|preprotein|protein translocase|chaperon|protein folding|protease|protein turnover", 
      `Fatty acid & phospholipid metabolism` = "fatty acid|phospholipid|glycerolipid|acetyl.coa carboxylase|acyl carrier protein", 
      `Cell division / envelope maintenance` = "peptidoglycan|cell wall|cell division|septum|septal|murein|muram|glucosamine|penicillin.binding", 
      `DNA replication, repair, chromosome segregation` = "dna replication|dna repair|homologous recombination|mismatch repair|excision repair|chromosome segregation|chromosome partition|dna gyrase|topoisomerase", 
      `Translation & RNA metabolism` = "ribosom|translation|rna polymerase|aminoacyl.trna|trna synthetase|rna degradation|rna helicase|rna processing", 
      `Transport / efflux` = "transporter|transport system|efflux|permease|abc transport|porin", 
      `Regulation / stress response` = "regulator|two.component|signal transduction|stress response|sigma factor|quorum sensing", 
      `Energy metabolism` = "oxidative phosphorylation|atp synthase|glycolysis|gluconeogenesis|citrate cycle|tca cycle|pyruvate|glyoxylate|pentose phosphate|electron transport|respiration|carbon metabolism", 
      `Other metabolism` = "metabolism|biosynthesis|synthetase|synthase|dehydrogenase|hydrolase|isomerase|kinase")
    hit <- which(vapply(rules, function(pattern) grepl(pattern, value, perl = TRUE), logical(1)))
    if (length(hit)) 
      names(rules)[hit[1L]]
    else "Other / unclassified"
  }
  out[, `:=`(annotation_category = vapply(text, classify, character(1)), annotation_present = meta$PAO1_ID %in% terms$locus_tag)]
  out
}

fig3_palette <- function(categories, palette_file = NULL, overrides = character()) {
  palette <- c(`Energy metabolism` = "#00A64F", `DNA replication, repair, chromosome segregation` = "#7D7D7D", 
    Lipopolysaccharide = "#8A4F07", `Protein export / quality control` = "#F36B00", `Fatty acid & phospholipid metabolism` = "#A6A100", 
    `Cell division / envelope maintenance` = "#D72E26", `Translation & RNA metabolism` = "#0057D9", 
    `Transport / efflux` = "#00CFE3", `Regulation / stress response` = "#E7298A", `Other metabolism` = "#333333", 
    `Other / unclassified` = "#D9D9D9")
  if (!is.null(palette_file)) {
    x <- s8_read(palette_file)
    s8_required(x, c("annotation_category", "colour"), palette_file)
    s8_unique(x, "annotation_category", palette_file)
    palette[x$annotation_category] <- x$colour
  }
  if (length(overrides)) {
    if (is.null(names(overrides)) || any(!nzchar(names(overrides))) || anyDuplicated(names(overrides))) 
      stop("Custom colours must have unique category names.")
    palette[names(overrides)] <- overrides
  }
  added <- sort(setdiff(unique(categories[!is.na(categories)]), names(palette)))
  if (length(added)) 
    palette[added] <- grDevices::hcl.colors(max(3L, length(added)), "Dark 3")[seq_along(added)]
  if (anyNA(palette) || inherits(try(grDevices::col2rgb(palette), silent = TRUE), "try-error")) 
    stop("Invalid functional annotation colour.")
  palette
}

fig3_pathoadaptive <- function(path, first_n = 224L) {
  raw <- data.table::fread(path, sep = "\t", header = TRUE, fill = TRUE, data.table = FALSE, 
    showProgress = FALSE)
  if (!ncol(raw)) 
    stop("Empty pathoadaptive table.")
  if (!is.null(first_n)) {
    if (length(first_n) != 1L || !is.finite(first_n) || first_n < 1L || first_n%%1 != 0 || 
      nrow(raw) < first_n) 
      stop("Invalid pathoadaptive first-block size.")
    raw <- raw[seq_len(first_n), , drop = FALSE]
  }
  ids <- trimws(as.character(raw[[1L]]))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) || any(!grepl("^PA[0-9]{4}(\\.[0-9]+)?$", 
    ids))) 
    stop("The selected pathoadaptive rows need unique PAO1 locus IDs.")
  ids
}

fig3_style_labels <- function(plot, font_family) {
  ix <- which(plot$gtable$layout$name == "col_names")
  if (length(ix) != 1L) 
    stop("Cannot locate gene labels in the heatmap.")
  g <- plot$gtable$grobs[[ix]]
  old <- as.list(g$gp)
  old$font <- NULL
  old$fontface <- "italic"
  old$fontfamily <- font_family
  old$col <- "black"
  g$gp <- do.call(grid::gpar, old)
  plot$gtable$grobs[[ix]] <- g
  plot
}

fig3_font_grob <- function(grob, family) {
  if (inherits(grob, "text")) {
    gp <- as.list(grob$gp)
    gp$fontfamily <- family
    grob$gp <- do.call(grid::gpar, gp)
  }
  if (inherits(grob, "gtable")) 
    grob$grobs <- lapply(grob$grobs, fig3_font_grob, family = family)
  if (inherits(grob, "gTree") && length(grob$children)) 
    grob$children <- do.call(grid::gList, lapply(grob$children, fig3_font_grob, family = family))
  grob
}

fig3cd_draw_heatmap <- function(plot, flexible = TRUE) {
  gt <- plot$gtable
  if (isTRUE(flexible)) {
    ix <- which(gt$layout$name == "matrix")
    if (length(ix) != 1L || gt$layout$t[ix] != gt$layout$b[ix]) 
      stop("Cannot locate the heatmap matrix row.")
    gt$heights[gt$layout$t[ix]] <- grid::unit(1, "null")
    gt$respect <- FALSE
  }
  grid::grid.newpage()
  grid::grid.draw(gt)
  invisible(plot)
}

print.amg_fig3cd <- function(x, ...) fig3cd_draw_heatmap(x, flexible = TRUE)

fig3cd_close_device <- function(device) {
  active <- grDevices::dev.list()
  if (!is.null(active) && device %in% unname(active)) 
    grDevices::dev.off(which = device)
  invisible(NULL)
}

fig3cd_bottom_annotations <- function(plot, annotation, colours, gene_order) {
  gt <- plot$gtable
  if (length(gt$heights) != 5L) 
    stop("Unsupported pheatmap layout: expected five layout rows.")
  row_order <- c(1L, 2L, 4L, 5L, 3L)
  old_to_new <- match(seq_len(5L), row_order)
  old_layout <- gt$layout
  gt$heights <- gt$heights[row_order]
  for (i in seq_len(nrow(old_layout))) {
    rows <- old_to_new[seq.int(old_layout$t[i], old_layout$b[i])]
    gt$layout$t[i] <- min(rows)
    gt$layout$b[i] <- max(rows)
  }
  if (!ncol(annotation)) {
    plot$gtable <- gt
    return(plot)
  }
  keys <- names(annotation)
  if (length(setdiff(keys, names(ANNOTATION_LABELS)))) 
    stop("Missing display labels for an annotation track.")
  n <- nrow(annotation)
  k <- length(keys)
  block_height <- ANNOTATION_TOP_GAP_PT + k * ANNOTATION_TRACK_HEIGHT_PT + (k - 1L) * ANNOTATION_TRACK_GAP_PT + 
    2
  centres <- block_height - ANNOTATION_TOP_GAP_PT - ANNOTATION_TRACK_HEIGHT_PT/2 - (seq_len(k) - 
    1L) * (ANNOTATION_TRACK_HEIGHT_PT + ANNOTATION_TRACK_GAP_PT)
  fill <- matrix(NA_character_, n, k, dimnames = list(rownames(annotation)[gene_order], keys))
  for (i in seq_along(keys)) {
    fill[, i] <- unname(colours[[keys[i]]][as.character(annotation[[keys[i]]][gene_order])])
  }
  if (anyNA(fill)) 
    stop("An annotation category has no display colour.")
  ix <- which(gt$layout$name == "col_annotation")
  nx <- which(gt$layout$name == "col_annotation_names")
  if (length(ix) != 1L || length(nx) != 1L) 
    stop("Cannot locate the annotation block in pheatmap.")
  gt$heights[5L] <- grid::unit(block_height, "bigpts")
  gt$grobs[[ix]] <- grid::rectGrob(x = grid::unit(rep((seq_len(n) - 0.5)/n, k), "npc"), y = grid::unit(rep(centres, 
    each = n), "bigpts"), width = grid::unit(1/n, "npc"), height = grid::unit(ANNOTATION_TRACK_HEIGHT_PT, 
    "bigpts"), gp = grid::gpar(fill = fill, col = NA))
  label_grobs <- list()
  label_widths <- numeric()
  for (i in seq_along(keys)) {
    text <- ANNOTATION_LABELS[[keys[i]]]
    if (!length(text) || length(text) > 2L || anyNA(text)) 
      stop("Annotation labels need one or two nonmissing lines.")
    face <- if (keys[i] == "Pathoadaptive status") 
      "bold"
    else "plain"
    y <- if (length(text) == 1L) 
      centres[i]
    else centres[i] + FONT_SIZE_ANNOTATION_DETAIL/2 + 1
    g <- grid::textGrob(text[1L], x = grid::unit(4, "bigpts"), y = grid::unit(y, "bigpts"), 
      hjust = 0, gp = grid::gpar(fontsize = FONT_SIZE_ANNOTATION, fontfamily = FONT_FAMILY, 
        fontface = face, col = "#4D4D4D"))
    label_grobs[[length(label_grobs) + 1L]] <- g
    label_widths <- c(label_widths, grid::convertWidth(grid::grobWidth(g), "bigpts", valueOnly = TRUE))
    if (length(text) == 2L) {
      g <- grid::textGrob(text[2L], x = grid::unit(4, "bigpts"), y = grid::unit(centres[i] - 
        FONT_SIZE_ANNOTATION/2 - 1, "bigpts"), hjust = 0, gp = grid::gpar(fontsize = FONT_SIZE_ANNOTATION_DETAIL, 
        fontfamily = FONT_FAMILY, fontface = face, col = "#4D4D4D"))
      label_grobs[[length(label_grobs) + 1L]] <- g
      label_widths <- c(label_widths, grid::convertWidth(grid::grobWidth(g), "bigpts", valueOnly = TRUE))
    }
  }
  gt$grobs[[nx]] <- grid::gTree(children = do.call(grid::gList, label_grobs))
  ri <- which(gt$layout$name == "row_names")
  if (length(ri) == 1L) {
    g <- gt$grobs[[ri]]
    g$gp$col <- "#4D4D4D"
    gt$grobs[[ri]] <- g
    label_widths <- c(label_widths, grid::convertWidth(grid::grobWidth(g), "bigpts", valueOnly = TRUE))
  }
  old_right <- gt$widths[4L]
  new_right <- grid::unit(max(label_widths) + 12, "bigpts")
  gt$widths[4L] <- new_right
  gt$widths[3L] <- gt$widths[3L] + old_right - new_right
  plot$gtable <- gt
  plot
}

fig3cd_build_heatmap <- function(result, labels, annotation, colours, title, tree) {
  measure_file <- tempfile("fig3cd_measure_", fileext = ".pdf")
  grDevices::cairo_pdf(measure_file, width = WIDTH_IN, height = HEIGHT_IN, family = FONT_FAMILY)
  device <- grDevices::dev.cur()
  on.exit({
    fig3cd_close_device(device)
    unlink(measure_file)
  }, add = TRUE)
  p <- pheatmap::pheatmap(t(result$shown), cluster_rows = FALSE, cluster_cols = tree, scale = "none", 
    color = (grDevices::colorRampPalette(c("#053061", "white", "#67001F")))(100), breaks = seq(-LFC_DISPLAY_CAP, 
      LFC_DISPLAY_CAP, length.out = 101), na_col = "grey88", labels_col = labels, annotation_col = if (ncol(annotation)) 
      annotation[, rev(names(annotation)), drop = FALSE]
    else NA, annotation_colors = if (ncol(annotation)) 
      colours
    else NA, annotation_names_col = TRUE, annotation_legend = SHOW_ANNOTATION_LEGENDS, legend = SHOW_VALUE_LEGEND, 
    main = title, fontsize_row = FONT_SIZE_CONDITION, fontsize_col = FONT_SIZE_GENE, cellheight = CELL_HEIGHT_PT, 
    angle_col = 90, treeheight_col = if (SHOW_GENE_TREE && !is.null(result$tree)) 
      TREE_HEIGHT_PT
    else 0, border_color = NA, fontfamily = "sans", silent = TRUE)
  grDevices::dev.set(device)
  p$gtable <- fig3_font_grob(p$gtable, FONT_FAMILY)
  p <- fig3_style_labels(p, FONT_FAMILY)
  gene_order <- if (is.null(result$tree)) 
    seq_len(nrow(result$metadata))
  else p$tree_col$order
  p <- fig3cd_bottom_annotations(p, annotation, colours, gene_order)
  class(p) <- c("amg_fig3cd", class(p))
  p
}

fig3cd_save_heatmap <- function(plot, file, type) {
  if (type == "pdf") {
    grDevices::cairo_pdf(file, width = WIDTH_IN, height = HEIGHT_IN, family = FONT_FAMILY, 
      onefile = TRUE)
  } else if (type == "png") {
    grDevices::png(file, width = WIDTH_IN, height = HEIGHT_IN, units = "in", res = DPI, type = "cairo")
  } else stop("Unsupported heatmap export format.")
  device <- grDevices::dev.cur()
  on.exit(fig3cd_close_device(device), add = TRUE)
  grDevices::dev.set(device)
  fig3cd_draw_heatmap(plot, flexible = FALSE)
  invisible(file)
}

show_Fig3CD <- function(direction = c("Sensitizing", "Protective"), plots = Fig3CD_plots) {
  direction <- match.arg(direction)
  if (is.null(plots[[direction]])) 
    stop("The requested direction was not generated.")
  fig3cd_draw_heatmap(plots[[direction]]$plot)
}

# Shared annotation key: combine the mappings used by the requested panels.
fig3cd_colour_key <- function(plots) {
  parts <- list()
  for (panel in plots) for (track in names(panel$colours)) {
    map <- panel$colours[[track]]
    parts[[length(parts) + 1L]] <- data.table::data.table(legend_group = track, category = names(map), 
      colour = unname(map), kind = "Annotation")
  }
  key <- if (length(parts)) 
    data.table::rbindlist(parts)
  else data.table::data.table(legend_group = character(), category = character(), colour = character(), 
    kind = character())
  if (nrow(key) && any(key[, .(n_colours = data.table::uniqueN(colour)), by = .(legend_group, category)]$n_colours != 1L)) 
    stop("Conflicting legend colours between panels.")
  key <- unique(key, by = c("legend_group", "category"))
  data.table::rbindlist(list(key, data.table::data.table(legend_group = "Interaction severity", 
    category = "Not significant or not evaluated", colour = "grey88", kind = "Masked cell")))
}

fig3cd_legend_layout <- function(key) {
  keys <- unique(key$legend_group[key$kind == "Annotation"])
  right <- if ("Functional class" %in% keys) 
    "Functional class"
  else if (length(keys) > 1L) 
    tail(keys, 2L)
  else character()
  left <- setdiff(keys, right)
  line_height <- LEGEND_FONT_SIZE + 3
  wrapped <- function(track) lapply(key$category[key$legend_group == track], function(x) strwrap(x, 
    width = 48L))
  block_height <- function(track) 22 + sum(vapply(wrapped(track), function(x) max(1L, length(x)) * 
    line_height, numeric(1))) + 12
  sum_height <- function(tracks) sum(vapply(tracks, block_height, numeric(1)))
  height_pt <- max(88 + sum_height(left), sum_height(right)) + 36
  list(left = left, right = right, height_in = height_pt/72, line_height = line_height, wrap = wrapped)
}

fig3cd_draw_colour_legend <- function(key, layout) {
  grid::grid.newpage()
  width <- LEGEND_WIDTH_IN * 72
  height <- layout$height_in * 72
  margin <- 18
  gap <- 24
  column_width <- (width - 2 * margin - gap)/2
  x_left <- margin
  x_right <- margin + column_width + gap
  top <- height - margin
  text <- function(label, x, y, size = LEGEND_FONT_SIZE, face = "plain", just = c("left", "top")) {
    grid::grid.text(label, x = grid::unit(x, "bigpts"), y = grid::unit(y, "bigpts"), just = just, 
      gp = grid::gpar(fontfamily = FONT_FAMILY, fontsize = size, fontface = face, col = "black", 
        lineheight = 1.2))
  }
  box <- function(colour, x, y) grid::grid.rect(x = grid::unit(x + 5, "bigpts"), y = grid::unit(y - 
    5, "bigpts"), width = grid::unit(10, "bigpts"), height = grid::unit(10, "bigpts"), gp = grid::gpar(fill = colour, 
    col = "grey50", lwd = 0.4))
  text("Interaction severity", x_left, top, LEGEND_FONT_SIZE + 1.5, "bold")
  text("Median adjusted log2FC", x_left, top - 16)
  bar_width <- column_width - 16
  gradient <- (grDevices::colorRampPalette(c("#053061", "white", "#67001F")))(100)
  grid::grid.rect(x = grid::unit(x_left + (seq_len(100L) - 0.5) * bar_width/100, "bigpts"), y = grid::unit(top - 
    38, "bigpts"), width = grid::unit(bar_width/100 + 0.05, "bigpts"), height = grid::unit(12, 
    "bigpts"), gp = grid::gpar(fill = gradient, col = NA))
  ticks <- seq(-LFC_DISPLAY_CAP, LFC_DISPLAY_CAP, length.out = 7L)
  for (i in seq_along(ticks)) text(format(signif(ticks[i], 3), trim = TRUE), x_left + (ticks[i] + 
    LFC_DISPLAY_CAP)/(2 * LFC_DISPLAY_CAP) * bar_width, top - 48, just = c("centre", "top"))
  box("grey88", x_left, top - 66)
  text("Not significant or not evaluated", x_left + 17, top - 66)
  draw_blocks <- function(tracks, x, y) {
    for (track in tracks) {
      title <- if (track %in% names(ANNOTATION_LABELS)) 
        paste(ANNOTATION_LABELS[[track]], collapse = " ")
      else track
      text(title, x, y, LEGEND_FONT_SIZE + 1.5, "bold")
      y <- y - 22
      rows <- key[key$legend_group == track]
      lines <- layout$wrap(track)
      for (i in seq_len(nrow(rows))) {
        box(rows$colour[i], x, y)
        text(paste(lines[[i]], collapse = "\n"), x + 17, y)
        y <- y - max(1L, length(lines[[i]])) * layout$line_height
      }
      y <- y - 12
    }
  }
  draw_blocks(layout$left, x_left, top - 88)
  draw_blocks(layout$right, x_right, top)
  invisible(key)
}

fig3cd_save_colour_legend <- function(key, directory) {
  layout <- fig3cd_legend_layout(key)
  for (type in c("pdf", "png")) {
    file <- file.path(directory, paste0("Fig3C_D_colour_legend.", type))
    if (type == "pdf") 
      grDevices::cairo_pdf(file, width = LEGEND_WIDTH_IN, height = layout$height_in, family = FONT_FAMILY)
    else grDevices::png(file, width = LEGEND_WIDTH_IN, height = layout$height_in, units = "in", 
      res = DPI, type = "cairo")
    device <- grDevices::dev.cur()
    tryCatch({
      grDevices::dev.set(device)
      fig3cd_draw_colour_legend(key, layout)
    }, finally = fig3cd_close_device(device))
  }
  invisible(layout)
}

# Read, prepare, plot and export in order.
s8_require(c("data.table", "pheatmap"))

if (!file.exists("AMG_cleaned_Git.Rproj")) stop("Open AMG_cleaned_Git.Rproj and run from its working directory.")

if (!is.finite(LFC_DISPLAY_CAP) || LFC_DISPLAY_CAP <= 0 || !length(DIRECTIONS_TO_PLOT) || anyDuplicated(DIRECTIONS_TO_PLOT) || 
  any(!DIRECTIONS_TO_PLOT %in% c("Sensitizing", "Protective"))) stop("Invalid heatmap settings.")

display_settings <- c(WIDTH_IN, HEIGHT_IN, CELL_HEIGHT_PT, TREE_HEIGHT_PT, FONT_SIZE_GENE, FONT_SIZE_CONDITION, 
  FONT_SIZE_ANNOTATION, FONT_SIZE_ANNOTATION_DETAIL, ANNOTATION_TRACK_HEIGHT_PT, LEGEND_WIDTH_IN, 
  LEGEND_FONT_SIZE)

if (any(!is.finite(display_settings)) || any(display_settings <= 0) || !is.finite(ANNOTATION_TRACK_GAP_PT) || 
  ANNOTATION_TRACK_GAP_PT < 0 || !is.finite(ANNOTATION_TOP_GAP_PT) || ANNOTATION_TOP_GAP_PT < 
  0) stop("Invalid plot dimensions, fonts or spacing.")

mode <- tolower(FUNCTIONAL_ANNOTATION_MODE)

# 2. Read final effects, recurrence and annotation
input <- fig3_inputs(INPUT_DIR)

manifest <- fig3_conditions(input$manifest, CONDITIONS_TO_PLOT, primary_only = TRUE)

selected <- input$recurrent[("Sensitizing" %in% DIRECTIONS_TO_PLOT & recurrent_sensitizing) | ("Protective" %in% DIRECTIONS_TO_PLOT & recurrent_protective)]

# Match functional categories to the union of displayed recurrent genes so
# both panels use one consistent colour mapping. Optional manual gene names
# change display labels only; PAO1_ID remains the key for all joins.
functional <- fig3_functional(selected, mode, MANUAL_ANNOTATION_FILE, KEGG_FILE, MANUAL_MISSING_CATEGORY, 
  KEGG_PRODUCT_FALLBACK)

palette <- if (mode == "none") character() else fig3_palette(functional$annotation_category, FUNCTIONAL_PALETTE_FILE, 
  CUSTOM_FUNCTIONAL_COLOURS)

palette <- palette[names(palette) %in% functional$annotation_category]

unresolved <- functional$PAO1_ID[!is.na(functional$annotation_category) & grepl("^Unresolved", functional$annotation_category)]

if (length(unresolved)) warning("Preserved unresolved manual category: ", paste(unresolved, collapse = ", "), 
  ". See manual annotation README.")

own <- NULL

if (SHOW_OUR_FITNESS_ANNOTATION) {
  own <- s8_read(OUR_FITNESS_FILE)
  s8_required(own, c("PAO1_ID", "CRISPRi_negative_fitness_union"), OUR_FITNESS_FILE)
  s8_unique(own, "PAO1_ID", OUR_FITNESS_FILE)
  own[, `:=`(CRISPRi_negative_fitness_union, fig3_logical(CRISPRi_negative_fitness_union, "CRISPRi_negative_fitness_union"))]
}

pathoadaptive <- if (SHOW_PATHOADAPTIVE_ANNOTATION) fig3_pathoadaptive(PATHOADAPTIVE_FILE, PATHOADAPTIVE_FIRST_N) else character()

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

stage <- OUTPUT_DIR

all_values <- all_annotations <- all_trees <- list()

plots <- list()

# 3. Draw each directional map with categorical clustering
# Build sensitizing and protective maps independently with the same workflow.
# Annotation tracks show published essentiality, this study's fitness calls,
# pathoadaptive membership and/or functional class according to the settings.
for (direction in DIRECTIONS_TO_PLOT) {
  panel <- if (direction == "Sensitizing") 
    "Fig3C"
  else "Fig3D"
  result <- fig3_heatmap_data(input, direction, manifest, LFC_DISPLAY_CAP, CLUSTER_METHOD)
  meta <- result$metadata
  fi <- match(meta$PAO1_ID, functional$PAO1_ID)
  labels <- fig3_label(meta$gene_name, meta$PAO1_ID)
  if (mode == "manual" && USE_MANUAL_GENE_NAMES) {
    custom <- functional$manual_gene_name[fi]
    keep <- !is.na(custom) & nzchar(custom)
    labels[keep] <- fig3_label(custom[keep], meta$PAO1_ID[keep])
  }
  annotation <- data.frame(row.names = meta$PAO1_ID)
  colours <- list()
  published <- c(E_core = "Core essential", E_core_conditional = "Core conditionally essential", 
    E = "Effect on fitness (other studies)", rest = "Other")[meta$category]
  our_flag <- our_evaluated <- rep(FALSE, nrow(meta))
  if (SHOW_PUBLISHED_ANNOTATION) {
    annotation[["Published essentiality"]] <- factor(published, levels = c("Core essential", 
      "Core conditionally essential", "Effect on fitness (other studies)", "Other"))
    colours[["Published essentiality"]] <- c(`Core essential` = "#D73027", `Core conditionally essential` = "#FC8D59", 
      `Effect on fitness (other studies)` = "#0070B0", Other = "#D9D9D9")
  }
  if (SHOW_OUR_FITNESS_ANNOTATION) {
    oi <- match(meta$PAO1_ID, own$PAO1_ID)
    our_evaluated <- !is.na(oi)
    our_flag <- own$CRISPRi_negative_fitness_union[oi]
    our_flag[is.na(our_flag)] <- FALSE
    annotation[["Our CRISPRi fitness"]] <- factor(ifelse(our_flag, "Identified", "Not identified"), 
      levels = c("Identified", "Not identified"))
    colours[["Our CRISPRi fitness"]] <- c(Identified = "#542788", `Not identified` = "#E6E6E6")
  }
  if (SHOW_PATHOADAPTIVE_ANNOTATION) {
    annotation[["Pathoadaptive status"]] <- factor(ifelse(meta$PAO1_ID %in% pathoadaptive, 
      "Pathoadaptive", "Not pathoadaptive"), levels = c("Pathoadaptive", "Not pathoadaptive"))
    colours[["Pathoadaptive status"]] <- c(Pathoadaptive = "#1B9E77", `Not pathoadaptive` = "#E6E6E6")
  }
  if (mode != "none") {
    panel_palette <- palette[names(palette) %in% functional$annotation_category[fi]]
    annotation[["Functional class"]] <- factor(functional$annotation_category[fi], levels = names(panel_palette))
    colours[["Functional class"]] <- panel_palette
  }
  title <- if (SHOW_TITLES) 
    paste0(direction, " recurrent genes\nSignificant interaction severity: median adjusted log2FC")
  else NA_character_
  tree <- if (is.null(result$tree)) 
    FALSE
  else result$tree
  stem <- paste0(panel, "_", tolower(direction), "_interaction_severity")
  plot <- fig3cd_build_heatmap(result, labels, annotation, colours, title, tree)
  # Check that the gene order in the rendered heatmap matches the recorded
  # cluster order before exporting its cells and annotations.
  actual_order <- if (is.null(result$tree)) 
    seq_len(nrow(meta))
  else plot$tree_col$order
  if (!identical(actual_order, result$order)) 
    stop("Heatmap and exported gene orders disagree.")
  fig3cd_save_heatmap(plot, file.path(stage, paste0(stem, ".pdf")), "pdf")
  fig3cd_save_heatmap(plot, file.path(stage, paste0(stem, ".png")), "png")
  ix <- cbind(rep(seq_len(nrow(meta)), each = nrow(manifest)), rep(seq_len(nrow(manifest)), times = nrow(meta)))
  # Retain each original effect, call, evaluation flag and displayed capped value.
  # These records distinguish masked cells from small or truly zero effects.
  values <- data.table::data.table(panel = panel, direction = direction, PAO1_ID = meta$PAO1_ID[ix[, 1L]], gene_name = labels[ix[, 1L]], display_order = meta$display_order[ix[, 1L]], comparison_id = manifest$comparison_id[ix[, 2L]], display_label = manifest$display_label[ix[, 2L]], medLFC = result$effects[ix], stouffer_call = result$calls[ix], 
    evaluated = result$evaluated[ix], significant = result$significant[ix], displayed_value = result$shown[ix])
  annotations <- data.table::data.table(panel = panel, direction = direction, PAO1_ID = meta$PAO1_ID, 
    gene_name = labels, priority_order = meta$priority_order, display_order = meta$display_order, 
    n_significant_conditions = meta$n_significant_conditions, published_essentiality = unname(published), 
    our_CRISPRi_negative_fitness = our_flag, our_fitness_evaluated = our_evaluated, pathoadaptive = meta$PAO1_ID %in% 
      pathoadaptive, annotation_category = functional$annotation_category[fi], functional_colour = if (mode == 
      "none") 
      NA_character_
    else unname(palette[functional$annotation_category[fi]]), manual_assignment_present = if (mode == 
      "manual") 
      functional$annotation_present[fi]
    else NA)
  data.table::setorder(values, display_order, comparison_id)
  data.table::setorder(annotations, display_order)
  all_values[[direction]] <- values
  all_annotations[[direction]] <- annotations
  if (!is.null(result$tree)) 
    all_trees[[direction]] <- data.table::data.table(panel = panel, direction = direction, 
      merge_step = seq_len(nrow(result$tree$merge)), left = result$tree$merge[, 1L], right = result$tree$merge[, 2L], height = result$tree$height)
  plots[[direction]] <- list(plot = plot, data = result, annotation = annotation, colours = colours, 
    labels = labels)
  if (isTRUE(PRINT_IN_RSTUDIO)) 
    fig3cd_draw_heatmap(plot)
  message(panel, ": ", nrow(meta), " recurrent ", tolower(direction), " genes; functional mode ", 
    FUNCTIONAL_ANNOTATION_MODE, ".")
}

# 4. Export the shared legend and exact plot data
# Export one shared legend for the effect scale and enabled annotation tracks,
# plus the palette and clustering tables. The RStudio drawing helper uses a
# flexible layout; export dimensions remain controlled by the settings above.
colour_key <- fig3cd_colour_key(plots)

if (isTRUE(EXPORT_COLOUR_LEGEND)) fig3cd_save_colour_legend(colour_key, stage)

source_paths <- input$source_paths

if (SHOW_OUR_FITNESS_ANNOTATION) source_paths <- c(source_paths, OUR_FITNESS_FILE)

if (SHOW_PATHOADAPTIVE_ANNOTATION) source_paths <- c(source_paths, PATHOADAPTIVE_FILE)

if (mode == "manual") source_paths <- c(source_paths, MANUAL_ANNOTATION_FILE)

if (mode == "kegg") source_paths <- c(source_paths, KEGG_FILE)

if (mode != "none" && !is.null(FUNCTIONAL_PALETTE_FILE)) source_paths <- c(source_paths, FUNCTIONAL_PALETTE_FILE)

settings <- s8_settings(c(functional_mode = FUNCTIONAL_ANNOTATION_MODE, KEGG_product_fallback = KEGG_PRODUCT_FALLBACK, 
  manual_missing_policy = MANUAL_MISSING_CATEGORY, manual_gene_names = USE_MANUAL_GENE_NAMES, 
  P_STOUFFER_CUT = input$p_cut, RECURRENCE_MIN = input$recurrence_min, displayed_conditions = paste(manifest$display_label, 
    collapse = "; "), recurrence_conditions = paste(input$manifest$display_label[input$manifest$count_for_recurrence], 
    collapse = "; "), LFC_display_cap = LFC_DISPLAY_CAP, clustering = paste("Hamming distance on stored categorical calls;", 
    CLUSTER_METHOD), missing_comparison_state_for_clustering = "Not significant; evaluated flag retained separately", 
  width_in = WIDTH_IN, height_in = HEIGHT_IN, cell_height_pt = CELL_HEIGHT_PT, annotation_position = "Below gene labels", 
  annotation_track_height_pt = ANNOTATION_TRACK_HEIGHT_PT, annotation_track_gap_pt = ANNOTATION_TRACK_GAP_PT, 
  annotation_legends = SHOW_ANNOTATION_LEGENDS, annotation_strip_exports = "Disabled", separate_colour_legend = EXPORT_COLOUR_LEGEND, 
  RStudio_layout = "Flexible matrix height and independent pane width; no fixed aspect ratio", 
  pheatmap_version = as.character(utils::packageVersion("pheatmap"))), stats::setNames(source_paths, 
  basename(source_paths)))

tables <- list(Fig3C_D_plotted_values = data.table::rbindlist(all_values), Fig3C_D_annotations = data.table::rbindlist(all_annotations), 
  Fig3C_D_clustering = if (length(all_trees)) data.table::rbindlist(all_trees) else data.table::data.table(panel = character(), 
    direction = character(), merge_step = integer(), left = integer(), right = integer(), height = numeric()), 
  Fig3_functional_palette_used = data.table::data.table(annotation_category = names(palette), 
    colour = unname(palette)), Analysis_settings = settings)

tables$Fig3C_D_colour_key <- colour_key

s8_stage(tables, stage)

invisible(OUTPUT_DIR)

message("Saved Fig. 3C/D and plain TSV data: ", OUTPUT_DIR)

if (isTRUE(EXPORT_COLOUR_LEGEND)) message("Colour legend files:\n", paste0("  ", normalizePath(file.path(OUTPUT_DIR, 
  paste0("Fig3C_D_colour_legend.", c("pdf", "png"))), winslash = "/", mustWork = TRUE), collapse = "\n"))

message("Underlying figure data:\n", paste0("  ", normalizePath(file.path(OUTPUT_DIR, paste0(names(tables), 
  ".tsv")), winslash = "/", mustWork = TRUE), collapse = "\n"))

Fig3CD_plots <- plots

