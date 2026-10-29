# Compare recurrent screen hits with Murray et al. aminoglycoside results.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("preprocessing_scripts/16_AMG_Murray2015_comparison.R")
# Inputs:
# - preprocessing_scripts/out_03_AMG_relevant_results/AG_gene_recurrence.tsv.
# - annotation_input/Murray_2015_mbo005152522sd2.xlsx: Supplementary Data Set 2, header row 3.
# - annotation_input/PA14_PA01_orthologs.csv: PA14-to-PAO1 mappings.
# Upstream / execution order:
# - preprocessing_scripts/03_AMG_relevant_results.R; Murray results and orthology are supplied
#   reference inputs.
# Outputs:
# - preprocessing_scripts/out_16_AMG_Murray2015_comparison/ — Murray_comparison.tsv,
#   Murray_mapped_results.tsv and Analysis_settings.tsv.
# Interpretation / operation:
# - GEN, NEO and TOB calls require FDR < 0.05 and |logFC| ≥ 1.
# - Multiple mapped matches select lowest FDR then largest |logFC|; all mapped IDs remain reported.

GENE_DIR <- "preprocessing_scripts/out_03_AMG_relevant_results"

MURRAY_FILE <- "annotation_input/Murray_2015_mbo005152522sd2.xlsx"

ORTHOLOG_FILE <- "annotation_input/PA14_PA01_orthologs.csv"

MURRAY_SHEET <- "Supplementary Data Set 2"

MURRAY_HEADER_ROW <- 3L

MURRAY_ANTIBIOTICS <- c("GEN", "NEO", "TOB")

MURRAY_FDR_CUT <- 0.05

MURRAY_ABS_LFC_CUT <- 1

OUTPUT_DIR <- "preprocessing_scripts/out_16_AMG_Murray2015_comparison"

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

s8_clean <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | toupper(x) %in% c("", "NA", "NAN", "NONE", "NULL")] <- NA_character_
  x
}

s8_collapse <- function(x, sep = "; ") {
  x <- sort(unique(s8_clean(x)))
  x <- x[!is.na(x)]
  if (length(x)) 
    paste(x, collapse = sep)
  else NA_character_
}

s8_settings <- function(settings, inputs = character()) {
  data.table::data.table(setting = c(names(settings), paste0("input:", names(inputs)), paste0("md5:", 
    names(inputs)), "R_version"), value = c(as.character(settings), unname(inputs), unname(tools::md5sum(inputs)), 
    R.version.string))
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

main <- function() {
  s8_require(c("data.table", "openxlsx"))
  if (!file.exists(MURRAY_FILE) || !file.exists(ORTHOLOG_FILE)) 
    stop("Required Murray source/mapping input is missing:\n", paste(c(MURRAY_FILE, ORTHOLOG_FILE)[!file.exists(c(MURRAY_FILE, 
      ORTHOLOG_FILE))], collapse = "\n"))
  # 1. Read this study's recurrent genes and the published aminoglycoside data.
  # Retain the configured GEN, NEO and TOB conditions for external comparison.
  genes <- s8_read(file.path(GENE_DIR, "AG_gene_recurrence.tsv"))
  raw <- data.table::as.data.table(openxlsx::read.xlsx(MURRAY_FILE, sheet = MURRAY_SHEET, startRow = MURRAY_HEADER_ROW, 
    check.names = FALSE))
  s8_required(raw, c("Condition", "Locus", "Gene", "logFC", "FDR"), MURRAY_FILE)
  external <- raw[trimws(as.character(Condition)) %in% MURRAY_ANTIBIOTICS, .(antibiotic = trimws(as.character(Condition)), 
    PA14_ID = s8_clean(Locus), Murray_gene = s8_clean(Gene), Murray_logFC = suppressWarnings(as.numeric(logFC)), 
    Murray_FDR = suppressWarnings(as.numeric(FDR)))]
  external <- external[!is.na(PA14_ID)]
  if (any(external$Murray_FDR < 0 | external$Murray_FDR > 1, na.rm = TRUE)) 
    stop("Murray FDR outside [0,1].")
  # 2. Map published PA14 loci to PAO1 using the supplied orthology table.
  # Only mapped records can be compared; missing mappings are not negative results.
  mapping <- data.table::as.data.table(utils::read.csv(ORTHOLOG_FILE, check.names = FALSE, stringsAsFactors = FALSE, 
    na.strings = c("", "NA")))
  s8_required(mapping, c("PA14_ID", "PAO1_ID"), ORTHOLOG_FILE)
  mapping <- unique(mapping[, .(PA14_ID = s8_clean(PA14_ID), PAO1_ID = s8_clean(PAO1_ID))])
  mapping <- mapping[!is.na(PA14_ID) & !is.na(PAO1_ID)]
  mapped <- merge(external, mapping, by = "PA14_ID", all = FALSE, allow.cartesian = TRUE, sort = FALSE)
  if (!nrow(mapped)) 
    stop("No Murray records map to PAO1; check the orthology input.")
  mapped[, `:=`(abs_LFC, abs(Murray_logFC))]
  # 3. When multiple published records map to one PAO1 gene/antibiotic, select
  # the smallest FDR, then the largest absolute effect. Retain all mapped PA14 IDs
  # as annotation, but do not average their effects or combine their p-values.
  data.table::setorderv(mapped, c("antibiotic", "PAO1_ID", "Murray_FDR", "abs_LFC"), c(1L, 1L, 
    1L, -1L), na.last = TRUE)
  best <- mapped[, .(Murray_PA14_IDs = s8_collapse(PA14_ID, sep = ";"), Murray_gene = Murray_gene[1L], 
    Murray_logFC = Murray_logFC[1L], Murray_FDR = Murray_FDR[1L]), by = .(antibiotic, PAO1_ID)]
  best[, `:=`(Murray_tested, is.finite(Murray_logFC) & is.finite(Murray_FDR))]
  best[, `:=`(Murray_direction, data.table::fcase(Murray_tested & Murray_FDR < MURRAY_FDR_CUT & 
    Murray_logFC <= -MURRAY_ABS_LFC_CUT, "Sensitizing", Murray_tested & Murray_FDR < MURRAY_FDR_CUT & 
    Murray_logFC >= MURRAY_ABS_LFC_CUT, "Protective", Murray_tested, "Not significant", default = "Not tested"))]
  # 4. Compare each recurrent gene/direction with the published directional calls.
  # Same-direction support in at least two external antibiotics is summarized
  # separately from any opposite-direction evidence.
  queries <- data.table::rbindlist(lapply(c("Sensitizing", "Protective"), function(direction) {
    x <- genes[get(paste0("recurrent_", tolower(direction))) == TRUE, .(PAO1_ID, gene_name)]
    x[, `:=`(recurrent_direction, direction)]
    x
  }))
  collapse_external_conditions <- function(x) s8_collapse(x, sep = ";")
  comparison <- lapply(seq_len(nrow(queries)), function(i) {
    id <- queries$PAO1_ID[i]
    direction <- queries$recurrent_direction[i]
    x <- best[PAO1_ID == id]
    same <- x$Murray_direction == direction & x$Murray_tested
    opposite <- x$Murray_direction %in% c("Sensitizing", "Protective") & x$Murray_direction != 
      direction
    data.table::data.table(PAO1_ID = id, gene_name = queries$gene_name[i], recurrent_direction = direction, 
      Murray_n_conditions_tested = sum(x$Murray_tested), Murray_same_direction_conditions = collapse_external_conditions(x$antibiotic[same]), 
      Murray_consensus_same_direction_min2 = sum(same) >= 2L, Murray_any_opposite_direction = any(opposite), 
      Murray_opposite_direction_conditions = collapse_external_conditions(x$antibiotic[opposite]))
  })
  out <- if (length(comparison)) 
    data.table::rbindlist(comparison)
  else data.table::data.table(PAO1_ID = character(), gene_name = character(), recurrent_direction = character(), 
    Murray_n_conditions_tested = integer(), Murray_same_direction_conditions = character(), 
    Murray_consensus_same_direction_min2 = logical(), Murray_any_opposite_direction = logical(), 
    Murray_opposite_direction_conditions = character())
  # 5. Add the individual published effects, FDRs and calls for each antibiotic,
  # then export the recurrent-gene comparison and the complete mapped reference.
  for (drug in MURRAY_ANTIBIOTICS) {
    x <- best[antibiotic == drug]
    idx <- match(out$PAO1_ID, x$PAO1_ID)
    if (drug == "GEN") 
      out[, `:=`(Murray_GEN_PA14_IDs = x$Murray_PA14_IDs[idx], Murray_GEN_gene = x$Murray_gene[idx])]
    out[, `:=`((paste0("Murray_", drug, c("_logFC", "_FDR", "_direction"))), list(x$Murray_logFC[idx], 
      x$Murray_FDR[idx], x$Murray_direction[idx]))]
  }
  columns <- c("PAO1_ID", "gene_name", "recurrent_direction", "Murray_n_conditions_tested", "Murray_same_direction_conditions", 
    "Murray_consensus_same_direction_min2", "Murray_any_opposite_direction", "Murray_opposite_direction_conditions", 
    "Murray_GEN_PA14_IDs", "Murray_GEN_gene", paste0("Murray_", MURRAY_ANTIBIOTICS, "_logFC"), 
    paste0("Murray_", MURRAY_ANTIBIOTICS, "_FDR"), paste0("Murray_", MURRAY_ANTIBIOTICS, "_direction"))
  out <- out[, ..columns]
  out[, `:=`(order, match(recurrent_direction, c("Sensitizing", "Protective")))]
  data.table::setorder(out, order, PAO1_ID)
  out[, `:=`(order, NULL)]
  settings <- s8_settings(c(MURRAY_SHEET = MURRAY_SHEET, MURRAY_HEADER_ROW = MURRAY_HEADER_ROW, 
    MURRAY_FDR_CUT = MURRAY_FDR_CUT, MURRAY_ABS_LFC_CUT = MURRAY_ABS_LFC_CUT, MURRAY_ANTIBIOTICS = paste(MURRAY_ANTIBIOTICS, 
      collapse = "; ")), c(recurrence = file.path(GENE_DIR, "AG_gene_recurrence.tsv"), Murray = MURRAY_FILE, 
    orthologs = ORTHOLOG_FILE))
  s8_write(list(Murray_comparison = out, Murray_mapped_results = best, Analysis_settings = settings), 
    OUTPUT_DIR, OVERWRITE)
  message("S8 Murray comparison: ", nrow(out), " recurrent gene/direction records.")
  invisible(out)
}

main()

