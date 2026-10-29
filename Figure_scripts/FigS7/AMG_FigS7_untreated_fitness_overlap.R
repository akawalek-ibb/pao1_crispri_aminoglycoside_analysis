# Supplementary Figure S7: overlap of untreated negative-fitness gene sets.
# Run from the AMG RStudio project working directory:
# source("Figure_scripts/FigS7/AMG_FigS7_untreated_fitness_overlap.R")
#
# One prepared input: MH_fitness.tsv from 12_AMG_MH_fitness.R, also used for
# Table S5. If that file is already present, this script runs directly.
# No other R script, raw guide table or external essential-gene dataset is read.
#
# Flow:
# 1. Read the five untreated gene-level comparisons from the prepared table.
# 2. Select genes with raw Stouffer p < 0.05 AND median log2FC <= -1.
# 3. Record each gene's membership and count exclusive/pairwise overlaps.
# 4. Draw the five-set Venn, print it in RStudio and export the source data.
# B1 T1 and T2 remain separate sets; they are time points of the same batch.
# No Stouffer combination, new multiple-testing adjustment or curve fit is run.
#
# R packages used/checked: ggVennDiagram, ggplot2.

# 1. USER SETTINGS -------------------------------------------------------------
# This is the same prepared fitness table used to assemble Supplementary S5.
MH_FITNESS_FILE <- "preprocessing_scripts/out_12_AMG_MH_fitness/MH_fitness.tsv"

OUTPUT_DIR <- "Figure_scripts/FigS7/out_AMG_FigS7_untreated_fitness_overlap"

# Raw gene-level p values; the depletion cutoff is inclusive.
P_STOUFFER_CUT <- 0.05
LFC_CUT <- 1

# Reproduce the manuscript union. Set NULL when deliberately using other data
# or thresholds; this check does not determine which genes enter the sets.
EXPECTED_UNION_N <- 743L

# Keep the five sets in their original order. prefix identifies columns in
# MH_fitness.tsv; set_label is the text drawn beside each Venn set.
comparison_manifest <- data.frame(
  set_id = c("B1_T1", "B1_T2", "B2_T1", "B3_T1", "B4_T1"),
  set_label = c("B1 control T1 vs T0", "B1 control T2 vs T0",
                "B2 control T1 vs T0", "B3 control T1 vs T0",
                "B4 control T1 vs T0"),
  prefix = c("B1_T1_vs_T0", "B1_T2_vs_T0", "B2_T1_vs_T0",
             "B3_T1_vs_T0", "B4_T1_vs_T0"),
  batch = c("B1", "B1", "B2", "B3", "B4"),
  condition = c("None_0_T1 - None_0_T0", "None_0_T2 - None_0_T0",
                "None_0_T1 - None_0_T0", "None_0_T1 - None_0_T0",
                "None_0_T1 - None_0_T0"),
  stringsAsFactors = FALSE
)

# Original Venn colour gradient, font and image dimensions.
VENN_LOW_COLOR <- "#FDE725FF"
VENN_HIGH_COLOR <- "#440154FF"
FONT_FAMILY <- "sans"
PLOT_TITLE <- "Negative fitness effects across untreated CRISPRi comparisons"
WIDTH_IN <- 12
HEIGHT_IN <- 10
DPI <- 300

# Space for long set labels outside the Venn panel.
SIDE_MARGIN_PT <- 72

PRINT_IN_RSTUDIO <- TRUE
OVERWRITE <- TRUE

# 2. READ THE PREPARED GENE RESULTS --------------------------------------------
packages <- c("ggplot2", "ggVennDiagram")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
if (!file.exists(MH_FITNESS_FILE)) {
  stop("Missing prepared fitness table: ", MH_FITNESS_FILE,
       "\nUse the existing task-12 output, or run 12_AMG_MH_fitness.R to create it.")
}

fitness <- utils::read.delim(MH_FITNESS_FILE, check.names = FALSE,
                           stringsAsFactors = FALSE, na.strings = c("", "NA"))
required_columns <- c("PAO1_ID", paste0(comparison_manifest$prefix, "_medLFC"),
                      paste0(comparison_manifest$prefix, "_p_stouffer"))
missing_columns <- setdiff(required_columns, names(fitness))
if (length(missing_columns)) {
  stop("Fitness table lacks: ", paste(missing_columns, collapse = ", "))
}
fitness$PAO1_ID <- trimws(as.character(fitness$PAO1_ID))
if (anyNA(fitness$PAO1_ID) || any(!nzchar(fitness$PAO1_ID)) || anyDuplicated(fitness$PAO1_ID)) {
  stop("MH_fitness.tsv must contain one nonblank PAO1_ID per gene.")
}

# 3. MAKE THE FIVE NEGATIVE-FITNESS SETS ----------------------------------------
# Read the stored effects/p values rather than relying on prewritten call
# labels. A gene needs finite values for BOTH fields to be evaluated.
# Keep locus tags as supplied, including small-gene IDs such as PA1552.1.
control_set_list <- vector("list", nrow(comparison_manifest))
names(control_set_list) <- comparison_manifest$set_label
statistics_list <- vector("list", nrow(comparison_manifest))
set_sizes <- comparison_manifest
set_sizes$n_genes_tested <- 0L
set_sizes$n_negative_fitness <- 0L

for (i in seq_len(nrow(comparison_manifest))) {
  cfg <- comparison_manifest[i, ]
  medLFC <- as.numeric(fitness[[paste0(cfg$prefix, "_medLFC")]])
  p_stouffer <- as.numeric(fitness[[paste0(cfg$prefix, "_p_stouffer")]])
  tested <- is.finite(medLFC) & is.finite(p_stouffer)
  negative_fitness_call <- tested & p_stouffer < P_STOUFFER_CUT & medLFC <= -LFC_CUT

  control_set_list[[i]] <- sort(fitness$PAO1_ID[negative_fitness_call])
  set_sizes$n_genes_tested[i] <- sum(tested)
  set_sizes$n_negative_fitness[i] <- sum(negative_fitness_call)

  # One source-data row per gene/comparison, including nonsignificant and
  # unevaluated genes so the threshold decisions remain inspectable.
  statistics_list[[i]] <- data.frame(
    PAO1_ID = fitness$PAO1_ID,
    set_id = cfg$set_id,
    set_label = cfg$set_label,
    batch = cfg$batch,
    condition = cfg$condition,
    medLFC = medLFC,
    p_stouffer = p_stouffer,
    tested = tested,
    negative_fitness_call = negative_fitness_call,
    stringsAsFactors = FALSE
  )
}
gene_statistics <- do.call(rbind, statistics_list)
all_negative_genes <- sort(unique(unlist(control_set_list, use.names = FALSE)))
if (!length(all_negative_genes)) stop("No genes pass the negative-fitness criteria.")
if (!is.null(EXPECTED_UNION_N) && length(all_negative_genes) != EXPECTED_UNION_N) {
  stop("Negative-fitness union contains ", length(all_negative_genes),
       " genes; expected ", EXPECTED_UNION_N, ". Check the fitness input and thresholds.")
}

# 4. RECORD MEMBERSHIP AND OVERLAP COUNTS --------------------------------------
# TRUE/FALSE columns describe which sets contain each gene in the union.
# RegionPattern follows comparison_manifest order: e.g. 11000 means B1 T1
# and B1 T2 only. RegionNumber counts genes in that exact exclusive region.
membership_matrix <- vapply(control_set_list, function(x) all_negative_genes %in% x,
                            logical(length(all_negative_genes)))
colnames(membership_matrix) <- comparison_manifest$set_id
membership <- data.frame(PAO1_ID = all_negative_genes,
                         membership_matrix, check.names = FALSE)

# Reuse any available gene annotations from the SAME input; no extra file.
annotation_columns <- intersect(c("PA14_ID", "gene_name", "product_name"), names(fitness))
if (length(annotation_columns)) {
  annotations <- fitness[match(all_negative_genes, fitness$PAO1_ID), annotation_columns, drop = FALSE]
  membership <- cbind(membership["PAO1_ID"], annotations,
                      membership[comparison_manifest$set_id])
}
membership$RegionPattern <- apply(membership_matrix, 1, function(x) paste0(as.integer(x), collapse = ""))
membership$n_sets <- rowSums(membership_matrix)
region_frequency <- table(membership$RegionPattern)
membership$RegionNumber <- as.integer(region_frequency[membership$RegionPattern])
membership <- membership[order(-membership$n_sets, -membership$RegionNumber, membership$PAO1_ID), ]
rownames(membership) <- NULL

# List every possible nonempty region, including regions with zero genes.
# Thus this file accounts for all 31 regions of the five-set Venn.
combinations <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), nrow(comparison_manifest))))
combinations <- combinations[rowSums(combinations) > 0, , drop = FALSE]
patterns <- apply(combinations, 1, function(x) paste0(as.integer(x), collapse = ""))
region_counts <- data.frame(
  RegionPattern = patterns,
  n_sets = rowSums(combinations),
  set_ids = apply(combinations, 1, function(x) paste(comparison_manifest$set_id[x], collapse = "; ")),
  n_genes = as.integer(table(factor(membership$RegionPattern, levels = patterns))),
  stringsAsFactors = FALSE
)
region_counts <- region_counts[order(region_counts$RegionPattern), ]

# Pairwise intersections/Jaccard indices summarise each pair of comparisons.
# Jaccard = shared genes / genes present in either set. Diagonal rows compare
# a set with itself; both directions are kept, as in the original script.
pairwise_rows <- list()
for (i in seq_along(control_set_list)) {
  for (j in seq_along(control_set_list)) {
    a <- control_set_list[[i]]
    b <- control_set_list[[j]]
    union_n <- length(union(a, b))
    pairwise_rows[[length(pairwise_rows) + 1L]] <- data.frame(
      set_A = comparison_manifest$set_id[i],
      set_B = comparison_manifest$set_id[j],
      n_A = length(a), n_B = length(b),
      intersection_n = length(intersect(a, b)), union_n = union_n,
      jaccard = if (union_n > 0) length(intersect(a, b)) / union_n else NA_real_,
      A_only_n = length(setdiff(a, b)), B_only_n = length(setdiff(b, a)),
      stringsAsFactors = FALSE
    )
  }
}
pairwise_overlap <- do.call(rbind, pairwise_rows)

# 5. DRAW THE VENN -------------------------------------------------------------
# Keep the original ggVennDiagram call and display choices. Counts are
# exclusive-region counts; the fill gradient reflects those same counts.
FigS7_plot <- ggVennDiagram::ggVennDiagram(
  control_set_list, label_alpha = 0, label = "count"
) +
  ggplot2::scale_fill_gradient(low = VENN_LOW_COLOR, high = VENN_HIGH_COLOR) +
  ggplot2::theme_minimal(base_size = 8, base_family = FONT_FAMILY) +
  ggplot2::labs(title = PLOT_TITLE, x = "", y = "") +
  ggplot2::coord_fixed(clip = "off") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                 plot.margin = ggplot2::margin(15, SIDE_MARGIN_PT, 15, SIDE_MARGIN_PT))

# 6. EXPORT FIGURE AND UNDERLYING DATA -----------------------------------------
# Six source tables: all gene statistics; union-gene memberships; all 31
# region counts; set sizes/comparison definitions; pairwise overlaps; settings.
# Compress only the larger full gene-statistics table. No workbook is created.
settings <- data.frame(
  setting = c("input_file", "input_md5", "p_value_metric", "P_STOUFFER_CUT",
              "negative_LFC_threshold", "n_negative_union", "R_version",
              "ggplot2_version", "ggVennDiagram_version"),
  value = c(MH_FITNESS_FILE, unname(tools::md5sum(MH_FITNESS_FILE)),
            "raw Stouffer p; no BH/FDR filter", as.character(P_STOUFFER_CUT),
            paste0("medLFC <= -", LFC_CUT), as.character(length(all_negative_genes)),
            R.version.string, as.character(utils::packageVersion("ggplot2")),
            as.character(utils::packageVersion("ggVennDiagram"))),
  stringsAsFactors = FALSE
)
source_tables <- list(
  FigS7_gene_statistics.tsv.gz = gene_statistics,
  FigS7_negative_fitness_membership.tsv = membership,
  FigS7_region_counts.tsv = region_counts,
  FigS7_set_sizes.tsv = set_sizes,
  FigS7_pairwise_overlap.tsv = pairwise_overlap,
  Analysis_settings.tsv = settings
)
figure_paths <- file.path(OUTPUT_DIR, paste0("FigS7_untreated_fitness_overlap", c(".png", ".tiff", ".pdf")))
source_paths <- file.path(OUTPUT_DIR, names(source_tables))
if (!OVERWRITE && any(file.exists(c(figure_paths, source_paths)))) {
  stop("An output exists; set OVERWRITE <- TRUE to replace it.")
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(source_tables)) {
  connection <- if (endsWith(source_paths[i], ".gz")) gzfile(source_paths[i], "wt") else file(source_paths[i], "wt")
  tryCatch(utils::write.table(source_tables[[i]], connection, sep = "\t", quote = TRUE,
                             row.names = FALSE, na = "NA"), finally = close(connection))
}
ggplot2::ggsave(figure_paths[1L], FigS7_plot, width = WIDTH_IN, height = HEIGHT_IN,
                units = "in", dpi = DPI, bg = "white")
ggplot2::ggsave(figure_paths[2L], FigS7_plot, width = WIDTH_IN, height = HEIGHT_IN,
                units = "in", dpi = DPI, compression = "lzw", bg = "white")
ggplot2::ggsave(figure_paths[3L], FigS7_plot, width = WIDTH_IN, height = HEIGHT_IN,
                units = "in", device = grDevices::cairo_pdf, bg = "white")

if (isTRUE(PRINT_IN_RSTUDIO)) print(FigS7_plot)
message("Saved Supplementary Fig. S7: ", length(all_negative_genes), " genes in the union.")
message("Figure files:\n", paste0("  ", normalizePath(figure_paths, winslash = "/"), collapse = "\n"))
message("Underlying figure data:\n", paste0("  ", normalizePath(source_paths, winslash = "/"), collapse = "\n"))
