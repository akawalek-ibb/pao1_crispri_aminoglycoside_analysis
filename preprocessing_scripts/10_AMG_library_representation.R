# Purpose: Calculate shared inputs for Supplementary Table S4 and Figure S6.
# Inputs: The generated count matrix and Sample_key.tsv, plus PAO1_curated.tsv.
# The raw-count exporter creates both required count/sample files. Headers define
# sample names and order. All guides are retained, including zero-count guides.
# Method: Detection is >=1 read; cloning is mean untreated-T0 count >=10
# in at least one batch. Counts are not normalized or filtered by fitness.
# Missing count entries are treated as zero for calculations, following the
# original QC analysis. The deposited values, including NA, are kept separately.
# Packages: Base R only.
# First run preprocessing_scripts/09_AMG_export_GEO_raw_counts.R.
# Then run from the open AMG_cleaned_Git.Rproj:
# source("preprocessing_scripts/10_AMG_library_representation.R")
# Output: Plain TSV tables beside this script in out_10_AMG_library_representation/.
# USER SETTINGS
COUNTS_FILE <- "preprocessing_scripts/out_09_AMG_export_GEO_raw_counts/CRISPRi_guide_raw_counts.tsv.gz"

ANNOTATION_FILE <- "annotation_input/PAO1_curated.tsv"

SAMPLE_KEY_FILE <- "preprocessing_scripts/out_09_AMG_export_GEO_raw_counts/Sample_key.tsv"

OUTPUT_DIR <- "preprocessing_scripts/out_10_AMG_library_representation"

DETECTION_MIN_READS <- 1L

CLONED_MEAN_T0_MIN_READS <- 10

DESIGN_CALL_ORDER <- c("main_essential", "other_essential", "non_essential")

read_tsv <- function(path) {
  if (!file.exists(path)) 
    stop("Missing input: ", path, call. = FALSE)
  connection <- if (endsWith(path, ".gz")) 
    gzfile(path, "rt")
  else file(path, "rt")
  on.exit(close(connection))
  utils::read.delim(connection, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "", 
    na.strings = "NA")
}

require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) 
    stop(label, " lacks: ", paste(missing, collapse = ", "), call. = FALSE)
}

# Local functions for this analysis.
main <- function() {
  if (!is.finite(DETECTION_MIN_READS) || DETECTION_MIN_READS < 1 || DETECTION_MIN_READS%%1 != 
    0) 
    stop("Detection cutoff must be a positive integer.")
  if (!is.finite(CLONED_MEAN_T0_MIN_READS) || CLONED_MEAN_T0_MIN_READS <= 0) 
    stop("Cloning cutoff must be positive.")
  # 1. Read the unfiltered guide-by-sample matrix from task 09.
  # Check its guide IDs, classes and GEO-format sample headers before counting
  # library representation. All designed key guides remain in the calculations.
  raw_counts <- read_tsv(COUNTS_FILE)
  require_columns(raw_counts, c("spacer", "type", "locus_tag"), "GEO count file")
  if (!identical(names(raw_counts)[1:3], c("spacer", "type", "locus_tag"))) 
    stop("The first three GEO columns must be spacer, type and locus_tag.")
  if (!nrow(raw_counts) || anyNA(raw_counts$spacer) || anyDuplicated(raw_counts$spacer)) 
    stop("GEO spacers must be present and unique.")
  if (any(!grepl("^[ACGT]{20}$", raw_counts$spacer))) 
    stop("Invalid 20-nt spacer sequence.")
  if (anyNA(raw_counts$type) || any(!raw_counts$type %in% c("perfect", "mismatch", "control"))) 
    stop("Unexpected guide type.")
  targeted <- raw_counts$type != "control"
  if (anyNA(raw_counts$locus_tag[targeted]) || any(!nzchar(raw_counts$locus_tag[targeted]))) 
    stop("Targeted guides lack primary gene assignments.")
  sample_names <- names(raw_counts)[-(1:3)]
  pattern <- "^(B[0-9]+)_([0-9]+)_(T[0-9]+)_(MH|FBS|TOB[0-9]+(?:\\.[0-9]+)?|AMK[0-9]+(?:\\.[0-9]+)?)$"
  if (!length(sample_names) || anyDuplicated(sample_names) || any(!grepl(pattern, sample_names, 
    perl = TRUE))) 
    stop("Count headers must use the canonical GEO sample names.")
  if (any(!vapply(raw_counts[sample_names], is.numeric, logical(1)))) 
    stop("Sample counts must be numeric or NA.")
  original_matrix <- as.matrix(raw_counts[sample_names])
  observed <- original_matrix[!is.na(original_matrix)]
  if (any(!is.finite(observed)) || any(observed < 0) || any(observed%%1 != 0)) 
    stop("Raw counts must be finite, nonnegative integers or NA.")
  # 2. Treat missing entries as zero for coverage calculations, while retaining
  # the original missing-value matrix for reporting and for unavailable T0 batches.
  counts <- original_matrix
  counts[is.na(counts)] <- 0
  fields <- regmatches(sample_names, regexec(pattern, sample_names, perl = TRUE))
  batch <- vapply(fields, `[`, character(1), 2)
  timing <- vapply(fields, `[`, character(1), 4)
  treatment <- vapply(fields, `[`, character(1), 5)
  sample_key <- read_tsv(SAMPLE_KEY_FILE)
  require_columns(sample_key, c("label", "batch", "drug or treatment", "dose", "timing", "rep", 
    "original_sample_id"), "Generated sample key")
  if (anyNA(sample_key) || anyDuplicated(sample_key$label) || anyDuplicated(sample_key$original_sample_id) || 
    !identical(sample_names, sample_key$label)) 
    stop("Generated Sample_key.tsv does not match the count headers; rerun the raw-count exporter.")
  if (!identical(batch, sample_key$batch) || !identical(timing, sample_key$timing)) 
    stop("Sample metadata disagree with count headers.")
  batches <- unique(batch)
  mean_columns <- paste0(batches, "_mean_no_drug_T0")
  batch_means <- matrix(NA_real_, nrow(counts), length(batches), dimnames = list(NULL, mean_columns))
  t0_indices <- list()
  # 3. For each batch, calculate each guide's mean count across untreated T0
  # replicates. A guide with no recorded T0 entries in that batch retains NA.
  for (i in seq_along(batches)) {
    columns <- which(batch == batches[i] & timing == "T0" & treatment == "MH")
    if (!length(columns)) 
      stop("No untreated T0 samples for ", batches[i])
    t0_indices[[i]] <- columns
    batch_means[, i] <- rowMeans(counts[, columns, drop = FALSE])
    batch_means[rowSums(!is.na(original_matrix[, columns, drop = FALSE])) == 0, i] <- NA_real_
  }
  batch_means[is.nan(batch_means)] <- NA_real_
  t0_columns <- unlist(t0_indices, use.names = FALSE)
  # 4. Distinguish detection from the stronger library-representation criterion.
  # Detected means at least the minimum count in any T0 sample; cloned means the
  # mean T0 count reaches the configured cutoff in at least one batch. These
  # descriptive labels do not remove guides or call a fitness phenotype.
  detected <- rowSums(counts[, t0_columns, drop = FALSE] >= DETECTION_MIN_READS) > 0
  cloned <- rowSums(batch_means >= CLONED_MEAN_T0_MIN_READS, na.rm = TRUE) > 0
  guide_coverage <- data.frame(spacer = raw_counts$spacer, guide_category = raw_counts$type, 
    detected_any_T0 = detected, cloned_any_batch_ge10 = cloned, check.names = FALSE)
  guide_coverage <- cbind(guide_coverage, as.data.frame(round(batch_means, 1), check.names = FALSE))
  annotation <- read_tsv(ANNOTATION_FILE)
  require_columns(annotation, c("PAO1_ID", "gene_name", "design_call"), "Gene annotation")
  if (anyNA(annotation$PAO1_ID) || anyDuplicated(annotation$PAO1_ID)) 
    stop("Gene annotation IDs must be present and unique.")
  genes <- sort(unique(raw_counts$locus_tag[targeted]))
  annotation_index <- match(genes, annotation$PAO1_ID)
  if (anyNA(annotation_index)) 
    stop("Some targeted genes are absent from PAO1_curated.tsv.")
  design <- annotation$design_call[annotation_index]
  if (anyNA(design) || any(!design %in% DESIGN_CALL_ORDER)) 
    stop("Missing or unsupported design_call for a targeted gene.")
  gene_index <- match(raw_counts$locus_tag, genes)
  n_by_gene <- function(keep) tabulate(gene_index[keep], nbins = length(genes))
  perfect <- raw_counts$type == "perfect"
  mismatch <- raw_counts$type == "mismatch"
  # 5. Count full-match and mismatch guides per target gene: designed, detected
  # and passing the T0-mean criterion. Summarize these by library-design category
  # and guide class for the supplementary table.
  gene_coverage <- data.frame(locus_tag = genes, gene_name = annotation$gene_name[annotation_index], 
    design_call = design, n_perfect_designed = n_by_gene(perfect), n_mismatch_designed = n_by_gene(mismatch), 
    n_perfect_cloned_10 = n_by_gene(perfect & cloned), n_mismatch_cloned_10 = n_by_gene(mismatch & 
      cloned), n_perfect_detected_any_T0 = n_by_gene(perfect & detected), n_mismatch_detected_any_T0 = n_by_gene(mismatch & 
      detected))
  design_summary <- do.call(rbind, lapply(DESIGN_CALL_ORDER, function(group) {
    x <- gene_coverage[gene_coverage$design_call == group & gene_coverage$n_perfect_designed > 
      0, , drop = FALSE]
    if (!nrow(x)) 
      stop("No full-match-targeted genes for design group ", group)
    data.frame(design_call = group, n_genes_with_designed_perfect = nrow(x), n_genes_with_perfect_cloned_10 = sum(x$n_perfect_cloned_10 > 
      0), n_genes_without_perfect_cloned_10 = sum(x$n_perfect_cloned_10 == 0), total_perfect_designed = sum(x$n_perfect_designed), 
      total_perfect_cloned_10 = sum(x$n_perfect_cloned_10), fraction_perfect_cloned_10 = sum(x$n_perfect_cloned_10)/sum(x$n_perfect_designed), 
      mean_perfect_designed_per_gene = mean(x$n_perfect_designed), median_perfect_designed_per_gene = stats::median(x$n_perfect_designed), 
      mean_perfect_cloned_10_per_gene = mean(x$n_perfect_cloned_10), median_perfect_cloned_10_per_gene = stats::median(x$n_perfect_cloned_10), 
      total_mismatch_designed = sum(x$n_mismatch_designed), total_mismatch_cloned_10 = sum(x$n_mismatch_cloned_10))
  }))
  library_summary <- do.call(rbind, lapply(c("perfect", "control", "mismatch", "ALL"), function(kind) {
    keep <- if (kind == "ALL") 
      rep(TRUE, nrow(raw_counts))
    else raw_counts$type == kind
    data.frame(guide_type = kind, n_guides = sum(keep), n_detected_any_T0 = sum(detected[keep]), 
      pct_detected_any_T0 = 100 * mean(detected[keep]), n_cloned_meanT0_ge10 = sum(cloned[keep]), 
      pct_cloned_meanT0_ge10 = 100 * mean(cloned[keep]))
  }))
  # 6. Prepare the underlying Figure S6 data rather than drawing the plots here.
  # A: total assigned counts; B: read-count bins by guide class; C: detected
  # full-match guides per gene/sample; D: represented full-match guides by design
  # category. These views describe library composition, not gene essentiality.
  panel_A <- data.frame(label = sample_names, batch = batch, total_reads = colSums(counts))
  read_bins <- c("0 reads", "1-10", "11-100", ">100")
  panel_B <- do.call(rbind, lapply(seq_along(sample_names), function(i) {
    do.call(rbind, lapply(c("perfect", "mismatch", "control"), function(kind) {
      x <- counts[raw_counts$type == kind, i]
      n <- tabulate(1L + (x > 0) + (x > 10) + (x > 100), nbins = 4L)
      data.frame(label = sample_names[i], sgRNA_class = kind, count_bin = read_bins, n_guides = n, 
        fraction = n/sum(n))
    }))
  }))
  per_gene_detected <- vapply(seq_along(sample_names), function(i) n_by_gene(perfect & counts[, 
    i] >= DETECTION_MIN_READS), integer(length(genes)))
  dimnames(per_gene_detected) <- list(genes, sample_names)
  panel_C <- do.call(rbind, lapply(seq_along(sample_names), function(i) {
    x <- per_gene_detected[gene_coverage$n_perfect_designed > 0, i]
    n <- tabulate(1L + (x >= 1) + (x >= 2) + (x >= 3) + (x > 10), nbins = 5L)
    data.frame(label = sample_names[i], guides_bin = c("0", "1", "2", "3-10", ">10"), n_genes = n, 
      fraction = n/sum(n))
  }))
  panel_D <- stats::aggregate(list(n_genes = rep(1L, sum(gene_coverage$n_perfect_designed > 0))), 
    by = list(design_call = design[gene_coverage$n_perfect_designed > 0], n_perfect_cloned_10 = gene_coverage$n_perfect_cloned_10[gene_coverage$n_perfect_designed > 
      0]), FUN = sum)
  panel_D <- panel_D[order(match(panel_D$design_call, DESIGN_CALL_ORDER), panel_D$n_perfect_cloned_10), 
    , drop = FALSE]
  settings <- data.frame(parameter = c("detection_min_reads", "cloned_mean_T0_min_reads", "count_source", 
    "count_source_md5", "missing_count_entries"), value = c(DETECTION_MIN_READS, CLONED_MEAN_T0_MIN_READS, 
    COUNTS_FILE, unname(tools::md5sum(COUNTS_FILE)), sum(is.na(original_matrix))))
  # 7. Write shared TSV inputs for Table S4 and Figure S6, with the sample key
  # and count-source checksum so both outputs use the same representation data.
  tables <- list(Sample_key = sample_key, Guide_coverage = guide_coverage, Gene_coverage = gene_coverage, 
    Design_summary = design_summary, Library_summary = library_summary, Figure_S6A = panel_A, 
    Figure_S6B = panel_B, Figure_S6C = panel_C, Figure_S6D = panel_D, Analysis_settings = settings)
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  for (name in names(tables)) utils::write.table(tables[[name]], file.path(OUTPUT_DIR, paste0(name, 
    ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  message("Saved ", length(tables), " TSV tables: ", OUTPUT_DIR, "\n", nrow(raw_counts), " guides; ", 
    length(sample_names), " samples; ", sum(detected), " detected; ", sum(cloned), " passing the cloning criterion.")
  invisible(tables)
}

main()

