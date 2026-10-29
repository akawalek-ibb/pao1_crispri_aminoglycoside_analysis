# RNA-seq analysis from featureCounts: CRISPRi knockdowns and TOB control
# Run from the manuscript project root, or Source in its RStudio project:
#   source("RNAseq_analysis/AMG_RNAseq_featureCounts_edgeR.R")
#
# Inputs: 18 per-sample count files in RNAseq_analysis/input/featureCounts/,
# annotation_input/PAO1_curated.tsv, and RNAseq_analysis/input/NC_002516.2_fixed2.gtf when counts
# lack Length. Count tables use Geneid and a final count column, or two columns.
# Upstream: supplied count inputs, or preprocessing/01_download_SRA_reads.sh ->
# 02_rename_SRA_reads.sh -> 03_preprocess_RNAseq_reads.sh and the manual count-file
# hand-off described in README_RNAseq_analysis.md. No upstream rerun is needed
# when using the supplied counts.
#
# All new results stay in RNAseq_analysis/out/. A rerun replaces results ONLY
# there. It does not copy to RNAseq_data/DE/, the frozen iModulon input folder,
# or any preprocessing/figure output folder. See README_RNAseq_analysis.md.
#

# Load libraries
suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(openxlsx)
  library(data.table)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(ggplot2)
  library(forcats)
  library(tibble)
})

# =========================================================
# RNA-seq analysis: CRISPRi silencing and tobramycin response
# =========================================================
# Purpose:
#   Differential-expression analysis for the CRISPRi-silencing RNA-seq
#   experiment. The analysis asks:
#     1. how tobramycin affects the non-targeting/control CRISPRi strain
#        control1, and
#     2. how silencing selected genes changes the transcriptome relative to
#        the same control1 reference condition.
#
#   Analysis flow:
#     1. load featureCounts tables,
#     2. remove core rRNA genes,
#     3. calculate FPKM-like normalized expression values,
#     4. generate PCA and expression-distribution plots,
#     5. perform edgeR differential-expression analysis,
#     6. generate volcano plots,
#     7. generate heatmaps and cluster-trend plots for significant genes.
#
# Input:
#   FeatureCount files are expected in featurecount_dir. File names are derived
#   from the trimmed FASTQ names by removing _R1.fastq.gz / _R2.fastq.gz and
#   adding .txt, for example:
#     10_control1.txt
#     13_ftsZg1m4.txt
#     25_control1_TOB.txt
#
#   Expected groups:
#     control1       control CRISPRi strain/reference
#     control1_TOB   control1 plus tobramycin
#     ftsZg1m4      ftsZ-silencing condition
#     atpIg1m2     atpI-silencing condition
#     accCg1m3      accC-silencing condition
#     secAg1m3      secA-silencing condition
#
# Output:
#   All tables and plots are written to out_dir.
# =========================================================

# =========================================================
# settings
# =========================================================

# Paths are relative to the manuscript project root.
analysis_dir <- "RNAseq_analysis"
featurecount_dir <- file.path(analysis_dir, "input", "featureCounts")
out_dir <- file.path(analysis_dir, "out")

# Reuse the curated annotation already present in the main bundle.
annotation_file <- "annotation_input/PAO1_curated.tsv"

# Needed for TPM only if the count files do not retain featureCounts Length.
# Supply the EXACT GTF used for featureCounts -a, not a replacement annotation.
counting_gtf_file <- file.path(analysis_dir, "input", "NC_002516.2_fixed2.gtf")

if (!dir.exists(featurecount_dir)) {
  stop("FeatureCount directory not found: ", featurecount_dir)
}
if (!file.exists(annotation_file)) {
  stop("Annotation file not found: ", annotation_file)
}

# Keep the project working directory unchanged. Every output/readback below
# uses an explicit path within this analysis folder.
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# Plot display and export settings
# =========================================================

# Print every ggplot and every pheatmap object in the RStudio
# Plots pane while the script runs.
PRINT_PLOTS_IN_RSTUDIO <- TRUE

# Retain the existing PNG and SVG outputs and add a robust PDF
# export for every ggplot. Heatmaps are exported as PNG and PDF.
SAVE_PNG <- TRUE
SAVE_PDF <- TRUE
SAVE_SVG <- TRUE

# Use cairo_pdf explicitly. This is more reliable for text and
# vector output than relying on the default interactive device.
save_ggplot_outputs <- function(
    plot_object,
    basename,
    width,
    height
) {
  if (isTRUE(SAVE_PNG)) {
    ggplot2::ggsave(
      filename = file.path(
        out_dir,
        paste0(
          basename,
          ".png"
        )
      ),
      plot = plot_object,
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  }

  if (isTRUE(SAVE_PDF)) {
    ggplot2::ggsave(
      filename = file.path(
        out_dir,
        paste0(
          basename,
          ".pdf"
        )
      ),
      plot = plot_object,
      width = width,
      height = height,
      units = "in",
      device = grDevices::cairo_pdf,
      bg = "white"
    )
  }

  if (isTRUE(SAVE_SVG)) {
    ggplot2::ggsave(
      filename = file.path(
        out_dir,
        paste0(
          basename,
          ".svg"
        )
      ),
      plot = plot_object,
      width = width,
      height = height,
      units = "in",
      bg = "white"
    )
  }
}

print_ggplot_in_rstudio <- function(
    plot_object
) {
  if (isTRUE(PRINT_PLOTS_IN_RSTUDIO) && interactive()) {
    print(
      plot_object
    )
  }
}

# pheatmap returns a grid gtable rather than a ggplot object.
# Build it with silent = TRUE, then explicitly draw the gtable
# both in RStudio and on each file device.
print_pheatmap_in_rstudio <- function(
    heatmap_object
) {
  if (isTRUE(PRINT_PLOTS_IN_RSTUDIO) && interactive()) {
    grid::grid.newpage()
    grid::grid.draw(
      heatmap_object$gtable
    )
  }
}

save_pheatmap_outputs <- function(
    heatmap_object,
    basename,
    width,
    height
) {
  if (isTRUE(SAVE_PNG)) {
    grDevices::png(
      filename = file.path(
        out_dir,
        paste0(
          basename,
          ".png"
        )
      ),
      width = width,
      height = height,
      units = "in",
      res = 300,
      bg = "white"
    )

    tryCatch(
      {
        grid::grid.newpage()
        grid::grid.draw(
          heatmap_object$gtable
        )
      },
      finally = {
        grDevices::dev.off()
      }
    )
  }

  if (isTRUE(SAVE_PDF)) {
    grDevices::cairo_pdf(
      filename = file.path(
        out_dir,
        paste0(
          basename,
          ".pdf"
        )
      ),
      width = width,
      height = height,
      bg = "white"
    )

    tryCatch(
      {
        grid::grid.newpage()
        grid::grid.draw(
          heatmap_object$gtable
        )
      },
      finally = {
        grDevices::dev.off()
      }
    )
  }
}

# Load annotation file
PAO1 <- read_tsv(annotation_file, show_col_types = FALSE)

# Define core rRNA genes
core_rRNA_tags <- c(
  "PA0668.1", "PA0668.4", "PA0668.5",
  "PA4280.1", "PA4280.2", "PA4280.5",
  "PA4690.1", "PA4690.2", "PA4690.5",
  "PA5369.1", "PA5369.2", "PA5369.5"
)

# =========================================================
# Load count files
# =========================================================

sample_info <- tibble::tribble(
  ~file,                 ~sample,          ~strain,
  "10_control1.txt",        "10_control1",       "control1",
  "11_control1.txt",        "11_control1",       "control1",
  "12_control1.txt",        "12_control1",       "control1",
  "13_ftsZg1m4.txt",       "13_ftsZg1m4",      "ftsZg1m4",
  "14_ftsZg1m4.txt",       "14_ftsZg1m4",      "ftsZg1m4",
  "15_ftsZg1m4.txt",       "15_ftsZg1m4",      "ftsZg1m4",
  "16_atpIg1m2.txt",      "16_atpIg1m2",     "atpIg1m2",
  "17_atpIg1m2.txt",      "17_atpIg1m2",     "atpIg1m2",
  "18_atpIg1m2.txt",      "18_atpIg1m2",     "atpIg1m2",
  "19_accCg1m3.txt",       "19_accCg1m3",      "accCg1m3",
  "20_accCg1m3.txt",       "20_accCg1m3",      "accCg1m3",
  "21_accCg1m3.txt",       "21_accCg1m3",      "accCg1m3",
  "22_secAg1m3.txt",       "22_secAg1m3",      "secAg1m3",
  "23_secAg1m3.txt",       "23_secAg1m3",      "secAg1m3",
  "24_secAg1m3.txt",       "24_secAg1m3",      "secAg1m3",
  "25_control1_TOB.txt",    "25_control1_TOB",   "control1_TOB",
  "26_control1_TOB.txt",    "26_control1_TOB",   "control1_TOB",
  "27_control1_TOB.txt",    "27_control1_TOB",   "control1_TOB"
) %>%
  mutate(
    path = file.path(featurecount_dir, file),
    strain = factor(
      strain,
      levels = c("control1", "control1_TOB", "ftsZg1m4", "atpIg1m2", "accCg1m3", "secAg1m3")
    )
  )

missing_files <- sample_info$path[!file.exists(sample_info$path)]
if (length(missing_files) > 0) {
  stop("Missing count files:\n", paste(missing_files, collapse = "\n"))
}

# Reader for the original two-column files or standard per-sample featureCounts
# tables. fread() automatically detects the table after the featureCounts
# comment preamble; count values, including fractional values, are unchanged.
read_count_file <- function(path, sample_name) {
  dt <- data.table::fread(
    path,
    data.table = FALSE,
    showProgress = FALSE
  )

  if (ncol(dt) < 2) {
    stop("Count file has fewer than two columns: ", path)
  }

  # Standard featureCounts output usually has Geneid in the first column and counts in the last column.
  if ("Geneid" %in% names(dt)) {
    locus <- as.character(dt[["Geneid"]])
    count <- dt[[ncol(dt)]]
  } else {
    locus <- as.character(dt[[1]])
    count <- dt[[ncol(dt)]]
  }

  count <- suppressWarnings(as.numeric(count))
  out <- data.frame(
    locus_tag = locus,
    count = count,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$locus_tag) & out$locus_tag != "", ]
  out <- out[!is.na(out$count), ]
  out <- out[!duplicated(out$locus_tag), ]

  names(out)[2] <- sample_name
  out
}

count_tables <- lapply(seq_len(nrow(sample_info)), function(i) {
  read_count_file(sample_info$path[i], sample_info$sample[i])
})

alldata_df <- Reduce(function(x, y) full_join(x, y, by = "locus_tag"), count_tables)
alldata_df[is.na(alldata_df)] <- 0

# =========================================================
# Export merged raw featureCounts matrix
# =========================================================
# This is the unnormalized gene-level count matrix directly assembled
# from the individual featureCounts files, before rRNA removal and
# before edgeR filtering/normalization.
merged_featurecounts <- alldata_df %>%
  dplyr::rename(gene_id = locus_tag) %>%
  dplyr::select(gene_id, dplyr::all_of(sample_info$sample))

readr::write_tsv(
  merged_featurecounts,
  file.path(out_dir, "CRISPRi_RNAseq_merged_featureCounts.tsv")
)

cat(
  "Merged raw featureCounts matrix saved:",
  file.path(out_dir, "CRISPRi_RNAseq_merged_featureCounts.tsv"),
  "\n"
)
cat(
  "Genes:", nrow(merged_featurecounts),
  "| Samples:", ncol(merged_featurecounts) - 1,
  "\n"
)

alldata <- as.data.frame(alldata_df)
rownames(alldata) <- alldata$locus_tag
alldata$locus_tag <- NULL

# Make sure columns are ordered according to sample_info.
alldata <- alldata[, sample_info$sample, drop = FALSE]

# Filter out rRNA genes
alldata <- alldata[!rownames(alldata) %in% core_rRNA_tags, ]

write.csv(
  data.frame(PAO1_ID = rownames(alldata), alldata, check.names = FALSE),
  file = file.path(out_dir, "CRISPRi_silencing_RNAseq_raw_counts_after_rRNA_filter.csv"),
  row.names = FALSE
)

# =========================================================
# FPKM normalization
# =========================================================

fpm <- sweep(as.matrix(alldata), 2, colSums(alldata) / 1e6, "/") %>%
  as.data.frame()

# Join gene length info from annotation
temp <- data.frame(PAO1_ID = rownames(alldata), stringsAsFactors = FALSE)

required_annot_cols <- c("PAO1_ID", "start", "end")
if (!all(required_annot_cols %in% names(PAO1))) {
  stop("Annotation file must contain columns: ", paste(required_annot_cols, collapse = ", "))
}

PAO1_subset <- PAO1[, c("PAO1_ID", "start", "end")]
temp <- left_join(temp, PAO1_subset, by = "PAO1_ID")

# Calculate gene length and handle strand orientation
temp <- mutate(temp, length = abs(end - start) + 1) %>%
  mutate(length = length / 1000)

# Normalize by gene length
fpm <- sweep(as.matrix(fpm), 1, temp$length, "/") %>%
  as.data.frame()

# Remove any NA rows
keep_complete <- complete.cases(fpm)
fpm <- fpm[keep_complete, ]

write.csv(
  data.frame(PAO1_ID = rownames(fpm), fpm, check.names = FALSE),
  file = file.path(out_dir, "CRISPRi_silencing_RNAseq_FPKM_like_values.csv"),
  row.names = FALSE
)

# =========================================================
# iModulon-compatible expression exports
# =========================================================
# Derive TPM independently from counts BEFORE expression filtering.
# Do not reuse fpm: its complete.cases() step can drop genes after the
# annotation join. Prefer the Length column from the same featureCounts
# files that supplied the counts. For trimmed two-column files, obtain CDS
# lengths from the EXACT GTF used for counting, rather than the annotation table.
# TPM = (count / length_bp) / sum(count / length_bp per sample) * 1e6.
# All new objects live in a local scope; failures warn without altering the
# existing FPKM, edgeR, PCA or downstream objects.
imodulon_export_ok <- FALSE
tryCatch(local({
  # Adjust only this export-specific path if the counting GTF is elsewhere.
  imodulon_gtf_file <- counting_gtf_file
  imodulon_counts <- as.matrix(alldata)
  storage.mode(imodulon_counts) <- "double"
  imodulon_gene_ids <- rownames(imodulon_counts)
  if (nrow(imodulon_counts) == 0L || ncol(imodulon_counts) == 0L ||
      any(!is.finite(imodulon_counts)) || any(imodulon_counts < 0) ||
      is.null(imodulon_gene_ids) || anyNA(imodulon_gene_ids) ||
      any(!nzchar(imodulon_gene_ids)) || anyDuplicated(imodulon_gene_ids) ||
      !identical(colnames(imodulon_counts), as.character(sample_info$sample))) {
    stop("Invalid count matrix or sample-column order for the iModulon export.")
  }

  imodulon_lengths <- rep(NA_real_, length(imodulon_gene_ids))
  imodulon_length_source <- rep("missing", length(imodulon_gene_ids))
  for (imodulon_i in seq_len(nrow(sample_info))) {
    imodulon_fc <- read.delim(
      sample_info$path[imodulon_i], comment.char = "#",
      check.names = FALSE, stringsAsFactors = FALSE
    )
    if (!all(c("Geneid", "Length") %in% names(imodulon_fc))) next
    if (anyDuplicated(imodulon_fc$Geneid)) {
      stop("Duplicate Geneid entries in ", sample_info$path[imodulon_i])
    }
    imodulon_match <- match(imodulon_gene_ids, as.character(imodulon_fc$Geneid))
    imodulon_fc_length <- suppressWarnings(as.numeric(imodulon_fc$Length))[imodulon_match]
    imodulon_fc_valid <- is.finite(imodulon_fc_length) & imodulon_fc_length > 0

    # Check the source counts too: the export must come from this same run.
    imodulon_fc_counts <- suppressWarnings(as.numeric(imodulon_fc[[ncol(imodulon_fc)]]))[imodulon_match]
    imodulon_present <- !is.na(imodulon_match)
    if (any(!is.finite(imodulon_fc_counts[imodulon_present])) ||
        any(abs(imodulon_fc_counts[imodulon_present] -
                imodulon_counts[imodulon_present, imodulon_i]) > 1e-8)) {
      stop("Source counts disagree with the merged matrix for ", sample_info$sample[imodulon_i])
    }
    imodulon_conflict <- imodulon_fc_valid & is.finite(imodulon_lengths)
    if (any(abs(imodulon_lengths[imodulon_conflict] -
                imodulon_fc_length[imodulon_conflict]) > 1e-8)) {
      stop("Feature lengths differ between samples; use one consistent counting annotation.")
    }
    imodulon_lengths[imodulon_fc_valid] <- imodulon_fc_length[imodulon_fc_valid]
    imodulon_length_source[imodulon_fc_valid] <- "featureCounts_Length"
  }

  # The user's shell command retains only Geneid and counts (no Length).
  # Match featureCounts -t CDS with the default gene_id grouping.
  # Merge overlapping CDS intervals within each gene/chromosome so bases are
  # counted once if a gene has multiple CDS records.
  if (any(!is.finite(imodulon_lengths))) {
    if (!file.exists(imodulon_gtf_file)) {
      stop("Counting GTF not found: ", imodulon_gtf_file,
           ". Set imodulon_gtf_file to the GTF used with featureCounts -a.")
    }
    imodulon_gtf <- read.delim(
      imodulon_gtf_file, header = FALSE, sep = "\t", comment.char = "#",
      quote = "", check.names = FALSE, stringsAsFactors = FALSE,
      col.names = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
    )
    imodulon_gtf <- imodulon_gtf[imodulon_gtf$feature == "CDS", , drop = FALSE]
    if (nrow(imodulon_gtf) == 0L) stop("Counting GTF has no CDS records.")
    imodulon_attr_matches <- regmatches(
      imodulon_gtf$attributes,
      regexec('(^|;)[[:space:]]*gene_id[[:space:]]+"([^"]+)"', imodulon_gtf$attributes)
    )
    imodulon_gtf$gene_id <- vapply(imodulon_attr_matches, function(x) {
      if (length(x) >= 3L) x[[3]] else NA_character_
    }, character(1))
    imodulon_gtf$start <- suppressWarnings(as.numeric(imodulon_gtf$start))
    imodulon_gtf$end <- suppressWarnings(as.numeric(imodulon_gtf$end))
    if (anyNA(imodulon_gtf$gene_id) || any(!nzchar(imodulon_gtf$gene_id)) ||
        any(!is.finite(imodulon_gtf$start)) || any(!is.finite(imodulon_gtf$end)) ||
        any(imodulon_gtf$start < 1) || any(imodulon_gtf$end < imodulon_gtf$start)) {
      stop("Invalid CDS gene_id or coordinates in the counting GTF.")
    }
    imodulon_interval_length <- function(rows) {
      chromosome_rows <- split(rows, imodulon_gtf$seqname[rows])
      sum(vapply(chromosome_rows, function(idx) {
        idx <- idx[order(imodulon_gtf$start[idx], imodulon_gtf$end[idx])]
        left <- imodulon_gtf$start[idx[1L]]
        right <- imodulon_gtf$end[idx[1L]]
        total <- 0
        for (k in idx[-1L]) {
          if (imodulon_gtf$start[k] <= right + 1) {
            right <- max(right, imodulon_gtf$end[k])
          } else {
            total <- total + right - left + 1
            left <- imodulon_gtf$start[k]
            right <- imodulon_gtf$end[k]
          }
        }
        total + right - left + 1
      }, numeric(1)))
    }
    imodulon_gtf_lengths <- vapply(
      split(seq_len(nrow(imodulon_gtf)), imodulon_gtf$gene_id),
      imodulon_interval_length, numeric(1)
    )
    imodulon_matched_gtf_lengths <- unname(imodulon_gtf_lengths[imodulon_gene_ids])
    imodulon_use_gtf <- !is.finite(imodulon_lengths) &
      is.finite(imodulon_matched_gtf_lengths) & imodulon_matched_gtf_lengths > 0
    imodulon_lengths[imodulon_use_gtf] <- imodulon_matched_gtf_lengths[imodulon_use_gtf]
    imodulon_length_source[imodulon_use_gtf] <- "counting_GTF_CDS_union"
  }
  imodulon_keep <- is.finite(imodulon_lengths) & imodulon_lengths > 0
  imodulon_length_audit <- data.frame(
    gene_id = imodulon_gene_ids,
    length_bp = imodulon_lengths,
    length_source = imodulon_length_source,
    included_in_TPM = imodulon_keep,
    present_in_existing_FPKM = imodulon_gene_ids %in% rownames(fpm),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  readr::write_tsv(imodulon_length_audit,
                  file.path(out_dir, "CRISPRi_silencing_RNAseq_iModulon_length_audit.tsv"))
  if (any(!imodulon_keep)) {
    warning(sum(!imodulon_keep), " genes lack valid lengths and are omitted from TPM. ",
            "See CRISPRi_silencing_RNAseq_iModulon_length_audit.tsv. ",
            "Do not substitute counts or invent lengths for these genes.")
  }
  if (!any(imodulon_keep)) stop("No valid gene lengths; no expression file generated.")

  imodulon_rates <- sweep(imodulon_counts[imodulon_keep, , drop = FALSE],
                          1, imodulon_lengths[imodulon_keep], "/")
  imodulon_totals <- colSums(imodulon_rates)
  if (any(!is.finite(imodulon_rates)) || any(!is.finite(imodulon_totals)) ||
      any(imodulon_totals <= 0)) stop("Invalid TPM rates or zero sample totals.")
  imodulon_tpm <- sweep(imodulon_rates, 2, imodulon_totals, "/") * 1e6
  imodulon_log2_tpm <- log2(imodulon_tpm + 1)
  imodulon_expression_table <- data.frame(
    gene_id = imodulon_gene_ids[imodulon_keep], imodulon_log2_tpm,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  imodulon_expression_file <- file.path(out_dir, "CRISPRi_silencing_RNAseq_log2_TPM.tsv")
  readr::write_tsv(imodulon_expression_table, imodulon_expression_file)

  # Use the same GEO-style group labels in the metadata and Python config.
  # The expression columns and metadata sample_id values refer to the same samples.
  imodulon_group_map <- c(
    "control1" = "control1",
    "control1_TOB" = "control1_TOB",
    "ftsZg1m4" = "ftsZg1m4",
    "atpIg1m2" = "atpIg1m2",
    "accCg1m3" = "accCg1m3",
    "secAg1m3" = "secAg1m3"
  )
  imodulon_groups <- unname(imodulon_group_map[as.character(sample_info$strain)])

  if (anyNA(imodulon_groups)) {
    warning(
      "iModulon expression was exported, but sample metadata was skipped: ",
      "a strain is missing from imodulon_group_map."
    )
  } else {
    imodulon_sample_table <- data.frame(
      sample_id = as.character(sample_info$sample),
      group = imodulon_groups,
      source_strain = as.character(sample_info$strain),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    imodulon_sample_file <- file.path(
      out_dir,
      "CRISPRi_silencing_RNAseq_iModulon_samples.tsv"
    )
    readr::write_tsv(imodulon_sample_table, imodulon_sample_file)
    cat("iModulon sample metadata saved:", imodulon_sample_file, "\n")
  }

  imodulon_export_ok <<- TRUE
  cat("iModulon log2(TPM + 1) matrix saved:", imodulon_expression_file, "\n")
  cat("Genes:", nrow(imodulon_expression_table),
      "| Samples:", ncol(imodulon_expression_table) - 1L, "\n")
  cat("Length sources: featureCounts =", sum(imodulon_length_source == "featureCounts_Length"),
      "| counting GTF =", sum(imodulon_length_source == "counting_GTF_CDS_union"),
      "| missing =", sum(!imodulon_keep), "\n")
  cat("Genes recovered relative to the existing FPKM matrix:",
      sum(imodulon_keep & !imodulon_length_audit$present_in_existing_FPKM), "\n")
}), error = function(e) {
  warning("iModulon export failed: ", conditionMessage(e),
          " The existing analysis objects have not been modified. ",
          "Do not use an older TPM file from this output directory.")
})


# =========================================================
# PCA plot
# =========================================================

groups <- sample_info$strain
names(groups) <- sample_info$sample

# Perform PCA on normalized expression values
t_fpm <- data.table::transpose(fpm)
colnames(t_fpm) <- rownames(fpm)
rownames(t_fpm) <- colnames(fpm)
t_fpm <- as.data.frame(t_fpm)
t_fpm <- t_fpm[, apply(t_fpm, 2, var) != 0, drop = FALSE]

pca_obj <- prcomp(t_fpm, center = TRUE, scale. = TRUE)

# Extract PCA data for plotting manually
pca_df <- as.data.frame(pca_obj$x[, 1:2])
pca_df$Sample <- rownames(pca_df)
pca_df$Group <- factor(groups[pca_df$Sample], levels = levels(sample_info$strain))

# Plot with non-overlapping labels
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = Sample), size = 4, max.overlaps = 100) +
  stat_ellipse(aes(group = Group), linetype = "dashed", na.rm = TRUE) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of CRISPRi-silencing RNA-seq expression profiles",
    x = "PC1",
    y = "PC2",
    color = "Condition"
  ) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

print_ggplot_in_rstudio(
  pca_plot
)

save_ggplot_outputs(
  plot_object = pca_plot,
  basename = "PCA_plot_CRISPRi_silencing_RNAseq",
  width = 8,
  height = 6
)

# =========================================================
# Boxplot FPKM
# =========================================================

fpm_plot_df <- fpm
fpm_plot_df$group <- rownames(fpm_plot_df)
fpm_m <- reshape2::melt(fpm_plot_df, id.vars = "group")

fpkm_plot <- ggplot(fpm_m, aes(variable, value)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title = "", x = "", y = "FPKM-like expression") +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.ticks.x = element_blank()
  )

print_ggplot_in_rstudio(
  fpkm_plot
)

save_ggplot_outputs(
  plot_object = fpkm_plot,
  basename = "FPKM_boxplot_CRISPRi_silencing_RNAseq",
  width = 10,
  height = 6
)

# =========================================================
# Differential expression
# =========================================================

# Define sample names in the same order as in alldata
sample_names <- sample_info$sample

# View the sample assignment
print(sample_info[, c("sample", "strain")])

# Model design with strain-specific coefficients
design <- model.matrix(~ 0 + strain, data = sample_info)
colnames(design) <- levels(sample_info$strain)

print(design)

# Create DGE object and preprocess
y <- DGEList(counts = alldata)
y$samples$group <- sample_info$strain

# edgeR expression filtering
keep <- filterByExpr(y, design = design, group = sample_info$strain)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

write.csv(
  data.frame(PAO1_ID = rownames(y$counts), y$counts, check.names = FALSE),
  file = file.path(out_dir, "CRISPRi_silencing_RNAseq_edgeR_counts_after_filterByExpr.csv"),
  row.names = FALSE
)

# Estimate dispersions and fit
y <- estimateDisp(y, design)
fit <- glmFit(y, design)

# Define contrasts
# Only planned comparisons are run:
#   1. tobramycin response in the control1 control background,
#   2. each silenced-gene condition versus the same control1 reference.
contrast_matrix <- makeContrasts(
  control1_TOB_vs_control1 = control1_TOB - control1,
  ftsZg1m4_vs_control1    = ftsZg1m4 - control1,
  atpIg1m2_vs_control1   = atpIg1m2 - control1,
  accCg1m3_vs_control1    = accCg1m3 - control1,
  secAg1m3_vs_control1    = secAg1m3 - control1,
  levels = design
)

# Extract expression matrices
logcpm <- as.data.frame(cpm(y, log = TRUE))
rawcounts <- as.data.frame(y$counts)

# Add locus_tag column for merging
logcpm$locus_tag <- rownames(logcpm)
rawcounts$locus_tag <- rownames(rawcounts)

# Rename contrasts for output files
contrast_renames <- list(
  "control1_TOB_vs_control1" = "control1_TOB_vs_control1",
  "ftsZg1m4_vs_control1"    = "ftsZg1m4_vs_control1",
  "atpIg1m2_vs_control1"   = "atpIg1m2_vs_control1",
  "accCg1m3_vs_control1"    = "accCg1m3_vs_control1",
  "secAg1m3_vs_control1"    = "secAg1m3_vs_control1"
)

# Run DE for contrasts
de_results_list <- list()

for (name in colnames(contrast_matrix)) {
  cat("Running:", name, "\n")

  lrt <- glmLRT(fit, contrast = contrast_matrix[, name])
  results <- lrt$table

  results$FDR <- p.adjust(results$PValue, method = "BH")
  results$PAO1_ID <- rownames(results)

  results <- results %>%
    mutate(FoldChange = 2^logFC) %>%
    left_join(PAO1, by = "PAO1_ID") %>%
    left_join(logcpm, by = c("PAO1_ID" = "locus_tag"), suffix = c("", "_logCPM")) %>%
    left_join(rawcounts, by = c("PAO1_ID" = "locus_tag"), suffix = c("", "_rawCount")) %>%
    dplyr::select(
      PAO1_ID,
      any_of(c("PA14_ID", "gene_name", "product_name", "strand", "start", "end")),
      logFC, logCPM, LR, PValue, FDR, FoldChange,
      matches("_logCPM$"), matches("_rawCount$"), everything()
    ) %>%
    distinct(PAO1_ID, .keep_all = TRUE)

  outname <- contrast_renames[[name]]
  de_results_list[[outname]] <- results

  openxlsx::write.xlsx(results, file.path(out_dir, paste0("PAO1_", outname, "_All.xlsx")))
  write.csv(results, file.path(out_dir, paste0("PAO1_", outname, "_All.csv")), row.names = FALSE)

  fc_2 <- filter(results, (FoldChange > 2 | FoldChange < 0.5) & FDR < 0.05) %>%
    arrange(desc(logFC))

  openxlsx::write.xlsx(fc_2, file.path(out_dir, paste0("PAO1_", outname, "_FC2.xlsx")))
  write.csv(fc_2, file.path(out_dir, paste0("PAO1_", outname, "_FC2.csv")), row.names = FALSE)
}

# =========================================================
# Volcano plots
# =========================================================

for (orig_name in names(contrast_renames)) {
  outname <- contrast_renames[[orig_name]]
  cat("Generating volcano plot for:", outname, "\n")

  # Load DE results
  res <- openxlsx::read.xlsx(file.path(out_dir, paste0("PAO1_", outname, "_All.xlsx")), sheet = 1)

  # Define gene label fallback
  res <- res %>%
    mutate(
      Label = ifelse(!is.na(gene_name) & gene_name != "", gene_name, PAO1_ID),
      negLogFDR = -log10(pmax(FDR, .Machine$double.xmin)),
      SigType = case_when(
        FDR < 0.05 & logFC > 1  ~ "Up",
        FDR < 0.05 & logFC < -1 ~ "Down",
        TRUE                    ~ "NS"
      )
    )

  sig_colors <- c("Up" = "darkred", "Down" = "navyblue", "NS" = "gray80")

  volcano <- ggplot(res, aes(x = logFC, y = negLogFDR)) +
    geom_point(aes(color = SigType), size = 2.2, alpha = 0.8) +
    scale_color_manual(values = sig_colors) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_text_repel(
      data = filter(res, SigType %in% c("Up", "Down")),
      aes(label = Label),
      size = 3,
      max.overlaps = 66,
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "gray60"
    ) +
    labs(
      title = paste("Volcano Plot:", outname),
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10]~FDR),
      color = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  print_ggplot_in_rstudio(
    volcano
  )

  save_ggplot_outputs(
    plot_object = volcano,
    basename = paste0(
      "Volcano_",
      outname
    ),
    width = 8,
    height = 6
  )
}

# =========================================================
# Heatmaps and clustered expression trends for each contrast
# =========================================================

fc_threshold <- 2
fdr_threshold <- 0.05
k_clusters <- 6
palette_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(25)

all_files <- file.path(
  out_dir,
  paste0(
    "PAO1_",
    unname(
      unlist(
        contrast_renames
      )
    ),
    "_All.xlsx"
  )
)

all_files <- all_files[
  file.exists(
    all_files
  )
]

if (length(all_files) == 0L) {
  stop(
    "No planned current-run _All.xlsx contrast files were found."
  )
}

for (all_file in all_files) {
  outname <- sub("^PAO1_(.*)_All\\.xlsx$", "\\1", basename(all_file))
  cat("Processing:", outname, "\n")

  # Load DE results
  de_results <- openxlsx::read.xlsx(all_file, sheet = 1)

  # Filter DE genes by PAO1_ID
  de_genes <- de_results %>%
    filter((FoldChange > fc_threshold | FoldChange < 1 / fc_threshold) & FDR < fdr_threshold) %>%
    pull(PAO1_ID) %>%
    unique()

  if (length(de_genes) < 5) {
    cat("Skipping", outname, "- too few DE genes.\n")
    next
  }

  # Normalize counts
  raw_counts <- as.data.frame(y$counts)
  lib_size <- colSums(y$counts)
  norm_counts <- sweep(raw_counts, 2, lib_size / 1e6, `/`)
  norm_counts$PAO1_ID <- rownames(norm_counts)

  # Filter for DE genes
  heatmap_data <- norm_counts %>%
    filter(PAO1_ID %in% de_genes)

  rownames(heatmap_data) <- heatmap_data$PAO1_ID
  heatmap_data <- heatmap_data[, setdiff(colnames(heatmap_data), "PAO1_ID"), drop = FALSE]

  # Z-score scaling
  scaled_counts <- t(scale(t(heatmap_data)))
  scaled_counts <- scaled_counts[complete.cases(scaled_counts), , drop = FALSE]

  if (nrow(scaled_counts) < 5) {
    cat("Skipping", outname, "- too few complete rows after scaling.\n")
    next
  }

  # Heatmap with clustering.
  # Build once, then explicitly draw the same gtable in RStudio,
  # PNG and PDF. This avoids blank PDF devices in loops.
  hm <- pheatmap(
    mat = scaled_counts,
    color = palette_colors,
    border_color = NA,
    scale = "none",
    clustering_method = "ward.D2",
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    treeheight_row = 30,
    treeheight_col = 30,
    cutree_rows = k_clusters,
    fontsize_row = 4,
    fontsize_col = 10,
    cellwidth = 20,
    cellheight = 2,
    main = paste(
      "Heatmap of DE Genes:",
      outname
    ),
    silent = TRUE
  )

  print_pheatmap_in_rstudio(
    hm
  )

  save_pheatmap_outputs(
    heatmap_object = hm,
    basename = paste0(
      "Heatmap_PAO1ids_",
      outname
    ),
    width = 12,
    height = 20
  )

  # Cluster assignments by PAO1_ID
  cluster_assignments <- cutree(hm$tree_row, k = k_clusters)
  clustered_ids <- names(cluster_assignments)

  scaled_df <- scaled_counts[clustered_ids, , drop = FALSE] %>%
    as.data.frame()
  scaled_df$PAO1_ID <- rownames(scaled_df)
  scaled_df$cluster <- cluster_assignments[scaled_df$PAO1_ID]

  # Long format for plotting
  long_data <- scaled_df %>%
    pivot_longer(cols = -c(PAO1_ID, cluster), names_to = "sample", values_to = "zscore") %>%
    mutate(sample = factor(sample, levels = colnames(scaled_counts)))

  cluster_sizes <- long_data %>%
    distinct(PAO1_ID, cluster) %>%
    count(cluster, name = "n_genes")

  long_data <- long_data %>%
    left_join(cluster_sizes, by = "cluster") %>%
    mutate(
      cluster_label = paste0("Cluster ", cluster, " (n = ", n_genes, ")"),
      cluster_label = factor(cluster_label, levels = unique(cluster_label[order(cluster)]))
    )

  # Plot cluster trends
  g <- ggplot(long_data, aes(x = sample, y = zscore, group = PAO1_ID)) +
    geom_line(alpha = 0.2, color = "gray50", linewidth = 0.4) +
    stat_summary(aes(group = 1), fun = mean, geom = "line", color = "darkred", linewidth = 1.2) +
    facet_wrap(~ cluster_label, scales = "free_y") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold")
    ) +
    labs(
      title = paste("Expression Patterns Across Clusters:", outname),
      y = "Z-score (scaled expression)",
      x = "Sample"
    )

  print_ggplot_in_rstudio(
    g
  )

  save_ggplot_outputs(
    plot_object = g,
    basename = paste0(
      "Cluster_LinePlot_",
      outname
    ),
    width = 10,
    height = 6
  )

  # Save table with cluster annotations
  clusters_df <- data.frame(
    PAO1_ID = names(cluster_assignments),
    cluster = cluster_assignments,
    stringsAsFactors = FALSE
  )

  de_results_with_clusters <- de_results %>%
    left_join(clusters_df, by = "PAO1_ID") %>%
    filter(!is.na(cluster))

  openxlsx::write.xlsx(
    de_results_with_clusters,
    file = file.path(out_dir, paste0("PAO1_", outname, "_FC2_withClusters.xlsx")),
    rowNames = FALSE
  )

  write.csv(
    de_results_with_clusters,
    file = file.path(out_dir, paste0("PAO1_", outname, "_FC2_withClusters.csv")),
    row.names = FALSE
  )
}

# =========================================================
# Merged selected DE gene set across planned CRISPRi contrasts
# =========================================================

merged_name <- "CRISPRi_silencing_selected_contrasts"

planned_result_files <- file.path(
  out_dir,
  paste0("PAO1_", unname(unlist(contrast_renames)), "_All.xlsx")
)

planned_result_files <- planned_result_files[file.exists(planned_result_files)]

if (length(planned_result_files) > 0) {

  de_tables <- lapply(planned_result_files, function(f) {
    nm <- sub("^PAO1_(.*)_All\\.xlsx$", "\\1", basename(f))
    dt <- openxlsx::read.xlsx(f, sheet = 1)
    dt$contrast_name <- nm
    dt
  })

  names(de_tables) <- sub("^PAO1_(.*)_All\\.xlsx$", "\\1", basename(planned_result_files))

  sig_gene_sets <- lapply(de_tables, function(dt) {
    dt %>%
      filter((FoldChange > fc_threshold | FoldChange < 1 / fc_threshold) & FDR < fdr_threshold) %>%
      pull(PAO1_ID) %>%
      unique()
  })

  de_genes <- unique(unlist(sig_gene_sets))

  sig_summary <- tibble(
    contrast = names(sig_gene_sets),
    n_sig_genes = vapply(sig_gene_sets, length, integer(1))
  )

  write.csv(
    sig_summary,
    file = file.path(out_dir, paste0("PAO1_", merged_name, "_significant_gene_counts.csv")),
    row.names = FALSE
  )

  if (length(de_genes) >= 5) {

    # Normalize counts
    raw_counts <- as.data.frame(y$counts)
    lib_size <- colSums(y$counts)
    norm_counts <- sweep(raw_counts, 2, lib_size / 1e6, `/`)
    norm_counts$PAO1_ID <- rownames(norm_counts)

    # Subset to DE genes and z-score
    heatmap_data <- norm_counts %>%
      filter(PAO1_ID %in% de_genes)

    rownames(heatmap_data) <- heatmap_data$PAO1_ID
    heatmap_data <- heatmap_data[, setdiff(colnames(heatmap_data), "PAO1_ID"), drop = FALSE]
    scaled_counts <- t(scale(t(heatmap_data)))
    scaled_counts <- scaled_counts[complete.cases(scaled_counts), , drop = FALSE]

    if (nrow(scaled_counts) >= 5) {

      # Heatmap and clustering.
      # Build once and explicitly render the gtable everywhere.
      hm <- pheatmap(
        mat = scaled_counts,
        color = palette_colors,
        border_color = NA,
        scale = "none",
        clustering_method = "ward.D2",
        cluster_cols = FALSE,
        cluster_rows = TRUE,
        cutree_rows = k_clusters,
        fontsize_row = 4,
        fontsize_col = 10,
        cellwidth = 20,
        cellheight = 2,
        main = paste(
          "Heatmap of DE Genes:",
          merged_name
        ),
        silent = TRUE
      )

      print_pheatmap_in_rstudio(
        hm
      )

      save_pheatmap_outputs(
        heatmap_object = hm,
        basename = paste0(
          "Heatmap_",
          merged_name
        ),
        width = 12,
        height = 20
      )

      # Clustering output
      cluster_assignments <- cutree(hm$tree_row, k = k_clusters)
      clustered_ids <- names(cluster_assignments)

      scaled_df <- scaled_counts[clustered_ids, , drop = FALSE] %>%
        as.data.frame()
      scaled_df$PAO1_ID <- rownames(scaled_df)
      scaled_df$cluster <- cluster_assignments[scaled_df$PAO1_ID]

      # Long format and line plot
      long_data <- scaled_df %>%
        pivot_longer(cols = -c(PAO1_ID, cluster), names_to = "sample", values_to = "zscore") %>%
        mutate(sample = factor(sample, levels = colnames(scaled_counts)))

      cluster_sizes <- long_data %>%
        distinct(PAO1_ID, cluster) %>%
        count(cluster, name = "n_genes")

      long_data <- long_data %>%
        left_join(cluster_sizes, by = "cluster") %>%
        mutate(
          cluster_label = paste0("Cluster ", cluster, " (n = ", n_genes, ")"),
          cluster_label = factor(cluster_label, levels = unique(cluster_label[order(cluster)]))
        )

      g <- ggplot(long_data, aes(x = sample, y = zscore, group = PAO1_ID)) +
        geom_line(alpha = 0.2, color = "gray50", linewidth = 0.4) +
        stat_summary(aes(group = 1), fun = mean, geom = "line", color = "darkred", linewidth = 1.2) +
        facet_wrap(~ cluster_label, scales = "free_y") +
        theme_minimal(base_size = 14) +
        theme(
          panel.grid.major.y = element_line(color = "gray90"),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          strip.text = element_text(face = "bold", size = 12),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 14, face = "bold")
        ) +
        labs(
          title = paste("Expression Patterns Across Clusters:", merged_name),
          y = "Z-score (scaled expression)",
          x = "Sample"
        )

      print_ggplot_in_rstudio(
        g
      )

      save_ggplot_outputs(
        plot_object = g,
        basename = paste0(
          "Cluster_LinePlot_",
          merged_name
        ),
        width = 10,
        height = 6
      )

      # Save cluster table
      clusters_df <- data.frame(
        PAO1_ID = names(cluster_assignments),
        cluster = cluster_assignments,
        stringsAsFactors = FALSE
      )

      # Merge DE metrics from all planned contrasts for de_genes
      merged_de <- tibble(PAO1_ID = de_genes)

      for (nm in names(de_tables)) {
        add_tbl <- de_tables[[nm]] %>%
          filter(PAO1_ID %in% de_genes) %>%
          dplyr::select(PAO1_ID, logFC, FoldChange, FDR, logCPM, PValue, LR) %>%
          rename_with(~ paste0(., "_", nm), -PAO1_ID)

        merged_de <- full_join(merged_de, add_tbl, by = "PAO1_ID")
      }

      merged_de <- merged_de %>%
        left_join(clusters_df, by = "PAO1_ID") %>%
        arrange(cluster, PAO1_ID)

      # Add annotations from PAO1 reference
      annotation_cols <- c("PAO1_ID", "PA14_ID", "gene_name", "product_name", "start", "end", "strand")
      annotation_cols <- annotation_cols[annotation_cols %in% colnames(PAO1)]

      merged_de <- merged_de %>%
        left_join(PAO1[, annotation_cols], by = "PAO1_ID") %>%
        relocate(PAO1_ID, any_of(c("PA14_ID", "gene_name", "product_name", "start", "end", "strand")), .before = everything())

      openxlsx::write.xlsx(
        merged_de,
        file = file.path(out_dir, paste0("PAO1_", merged_name, "_CombinedDE_withClusters.xlsx")),
        rowNames = FALSE
      )

      write.csv(
        merged_de,
        file = file.path(out_dir, paste0("PAO1_", merged_name, "_CombinedDE_withClusters.csv")),
        row.names = FALSE
      )
    }
  } else {
    cat("Too few DE genes in selected CRISPRi contrast set.\n")
  }
}

cat("\nAnalysis complete.\n")
cat("All ggplots and heatmaps were printed in RStudio:", PRINT_PLOTS_IN_RSTUDIO, "\n")
cat("PNG export:", SAVE_PNG, "\n")
cat("PDF export using cairo_pdf:", SAVE_PDF, "\n")
cat("SVG export:", SAVE_SVG, "\n")
cat("Output directory:", out_dir, "\n")

# Record this run's installed versions; using another edgeR release can change
# numerical results. This records the environment without changing any tests.
writeLines(capture.output(utils::sessionInfo()),
           file.path(out_dir, "RNAseq_sessionInfo.txt"))

cat("\nUnderlying data generated by this run:\n")
cat("Merged counts: ", file.path(out_dir, "CRISPRi_RNAseq_merged_featureCounts.tsv"), "\n", sep = "")
cat("edgeR input after expression filtering: ",
    file.path(out_dir, "CRISPRi_silencing_RNAseq_edgeR_counts_after_filterByExpr.csv"), "\n", sep = "")
current_tpm_file <- file.path(out_dir, "CRISPRi_silencing_RNAseq_log2_TPM.tsv")
if (isTRUE(imodulon_export_ok) && file.exists(current_tpm_file)) {
  cat("log2(TPM + 1): ", current_tpm_file, "\n", sep = "")
} else {
  cat("TPM export was not generated; review the export warning and counting-GTF input.\n")
}
cat("Complete DE results (five comparisons):\n")
cat(paste0("  ", file.path(out_dir, paste0("PAO1_", unname(unlist(contrast_renames)), "_All.xlsx")),
           collapse = "\n"), "\n")
cat("The existing frozen manuscript inputs have not been replaced.\n")
