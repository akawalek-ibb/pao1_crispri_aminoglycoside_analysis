# Compile essentiality and antimicrobial-fitness evidence for PAO1 library design.
# Prerequisites: data.table, readxl; open AMG_cleaned_Git.Rproj in RStudio.
# Run from the project root. Source-study calls and mapping multiplicity are retained.
# Outputs: library_design/out_AMG_collect_library_design_evidence/ contains
# PAO1_library_design_evidence.tsv, PAO1_library_source_evidence.tsv.gz and
# PAO1_library_evidence_summary.tsv. This script does not generate guide designs
# or assemble a formatted Supplementary Table S2 workbook.
# Run: source("library_design/AMG_collect_library_design_evidence.R")
#
# SOURCE FILES
# Source tables and annotations are in library_design/input/; orthology is read
# from annotation_input/PA14_PA01_orthologs.csv. No upstream script is required.
# FILES below lists each required filename. Published shortlists come from the
# named studies; PAO1161 results are an author-provided unpublished input.
# Prefer the existing shared annotation_input/PAO1_curated.tsv when available.
# Large guide-prediction files are not required by this evidence-collection step.
# USER SETTINGS
INPUT_DIR <- file.path("library_design", "input")

OUTPUT_DIR <- file.path("library_design", "out_AMG_collect_library_design_evidence")

CURATED_ANNOTATION_FILE <- file.path("annotation_input", "PAO1_curated.tsv")


FILES <- c(annotation = "Pseudomonas_aeruginosa_PAO1_107.csv", annotation_extra = "PAO1.csv", orthologues = "PA14_PA01_orthologs.csv", 
  poulsen_categories = "pnas.1900570116.sd05.xlsx", poulsen_comparison = "pnas.1900570116.sd06.xlsx", 
  turner = "Turner_2015_S1_short.xlsx", skurnik = "Skurnik_2013_table_short.xlsx", pao1161 = "All_T1_gumbel.xlsx", 
  gallagher = "Gallagher_2011_short.xlsx", murray = "Murray_2015mbo005152522sd2.xlsx", vitale = "Table_sum.XLSX")

PAO1161_DESIGN_CALLS <- c("E", "EB", "U")

PRIMARY_MIN_ESSENTIALITY_SOURCES <- 4L

MURRAY_CONDITION <- "TOB"

MURRAY_FDR_CUTOFF <- 0.05

MURRAY_MIN_ABS_LOGFC <- 1

run_evidence_collection <- function() {
  packages <- c("data.table", "readxl")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) 
    stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!file.exists("AMG_cleaned_Git.Rproj")) 
    stop("Open AMG_cleaned_Git.Rproj and run from its working directory.", call. = FALSE)
  paths <- setNames(file.path(INPUT_DIR, FILES), names(FILES))
  paths[["orthologues"]] <- "annotation_input/PA14_PA01_orthologs.csv"
  absent <- paths[!file.exists(paths)]
  if (length(absent)) 
    stop("Required inputs are missing:\n", paste(absent, collapse = "\n"), call. = FALSE)
  if (!length(PAO1161_DESIGN_CALLS) || any(!PAO1161_DESIGN_CALLS %in% c("E", "EB", "U", "NE", 
    "S"))) 
    stop("Invalid PAO1161_DESIGN_CALLS.")
  if (PRIMARY_MIN_ESSENTIALITY_SOURCES < 1L || PRIMARY_MIN_ESSENTIALITY_SOURCES > 6L) 
    stop("Invalid primary-group threshold.")
  clean <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x) | x %in% c("", "NA", "-")] <- NA_character_
    x
  }
  require_columns <- function(x, required, label) {
    absent <- setdiff(required, names(x))
    if (length(absent)) 
      stop(label, " lacks columns: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  read_text <- function(path) as.data.frame(data.table::fread(file = path, colClasses = "character", 
    na.strings = c("", "NA"), check.names = FALSE))
  read_excel <- function(path, skip = 0L) as.data.frame(suppressMessages(readxl::read_excel(path, 
    skip = skip, col_types = "text", .name_repair = "unique")))
  write_table <- function(x, name) {
    data.table::fwrite(x, file.path(OUTPUT_DIR, name), sep = "\t", na = "NA", quote = "auto")
  }
  annotation <- read_text(paths[["annotation"]])
  require_columns(annotation, c("Locus Tag", "Feature Type", "Gene Name", "Product Name"), "PAO1 annotation")
  genes <- data.frame(PAO1_ID = clean(annotation[["Locus Tag"]]), feature_type = clean(annotation[["Feature Type"]]), 
    gene_name = clean(annotation[["Gene Name"]]), product_name = clean(annotation[["Product Name"]]), 
    stringsAsFactors = FALSE)
  if (anyNA(genes$PAO1_ID) || anyDuplicated(genes$PAO1_ID)) 
    stop("PAO1 annotation must have one nonempty record per locus.")
  genes$gene_name[is.na(genes$gene_name)] <- genes$PAO1_ID[is.na(genes$gene_name)]
  annotation_extra_path <- CURATED_ANNOTATION_FILE
  extra <- read_text(annotation_extra_path)
  require_columns(extra, c("PAO1_ID", "PA14_ID", "TU"), "Additional annotation")
  extra$PAO1_ID <- clean(extra$PAO1_ID)
  extra$PA14_ID <- clean(extra$PA14_ID)
  extra$TU <- clean(extra$TU)
  first_nonmissing <- function(x) {
    x <- unique(clean(x))
    x <- x[!is.na(x)]
    if (length(x)) 
      x[[1L]]
    else NA_character_
  }
  extra_groups <- split(seq_len(nrow(extra)), extra$PAO1_ID)
  annotation_value <- function(id, column) {
    ix <- extra_groups[[id]]
    if (is.null(ix)) 
      return(NA_character_)
    values <- unique(clean(extra[[column]][ix]))
    values <- values[!is.na(values)]
    if (column == "TU" && length(values) > 1L) 
      stop("Conflicting TU annotation for ", id, "; supply the current PAO1_curated.tsv.")
    first_nonmissing(values)
  }
  genes$PA14_ID <- vapply(genes$PAO1_ID, annotation_value, character(1), column = "PA14_ID")
  genes$TU <- vapply(genes$PAO1_ID, annotation_value, character(1), column = "TU")
  orthologues <- read_text(paths[["orthologues"]])
  require_columns(orthologues, c("PA14_ID", "PAO1_ID"), "Orthologue table")
  orthologues <- unique(data.frame(PA14_ID = clean(orthologues$PA14_ID), PAO1_ID = clean(orthologues$PAO1_ID)))
  orthologues <- orthologues[!is.na(orthologues$PA14_ID) & !is.na(orthologues$PAO1_ID), , drop = FALSE]
  chunks <- list()
  add_source <- function(source, source_gene_id, source_call, supports_design, kind, organism = "PAO1", 
    source_logFC = NA_real_, source_FDR = NA_real_, mapped_gene_id = source_gene_id) {
    n <- length(source_gene_id)
    if (!n) 
      return(invisible(NULL))
    x <- data.frame(source = source, source_gene_id = clean(source_gene_id), source_call = clean(source_call), 
      supports_design = as.logical(supports_design), evidence_kind = kind, source_logFC = source_logFC, 
      source_FDR = source_FDR, stringsAsFactors = FALSE)
    if (organism != "PA14") 
      x$PAO1_ID <- clean(mapped_gene_id)
    x <- x[!is.na(x$source_gene_id), , drop = FALSE]
    x$supports_design[is.na(x$supports_design)] <- FALSE
    if (organism == "PA14") {
      x <- merge(x, orthologues, by.x = "source_gene_id", by.y = "PA14_ID", all.x = TRUE, 
        sort = FALSE)
    }
    x$in_annotation <- !is.na(x$PAO1_ID) & x$PAO1_ID %in% genes$PAO1_ID
    chunks[[length(chunks) + 1L]] <<- unique(x)
    invisible(NULL)
  }
  poulsen <- read_excel(paths[["poulsen_categories"]])
  require_columns(poulsen, c("PA14_ID", "Essential_Category"), "Poulsen categories")
  add_source("Poulsen_2019", poulsen$PA14_ID, poulsen$Essential_Category, TRUE, "essentiality", 
    "PA14")
  poulsen_mapped <- merge(poulsen, orthologues, by = "PA14_ID", all.x = TRUE, sort = FALSE)
  core <- unique(poulsen_mapped$PAO1_ID[poulsen_mapped$Essential_Category == "Core"])
  conditional <- unique(poulsen_mapped$PAO1_ID[poulsen_mapped$Essential_Category == "Conditional"])
  published <- read_excel(paths[["poulsen_comparison"]])
  study_columns <- c(Lee_2015 = "Lee S.A. et al 2015", Liberati_2006 = "Liberati N.T. et al 2006")
  require_columns(published, c("PA14_ID", unname(study_columns)), "Poulsen published-study comparison")
  for (study in names(study_columns)) {
    call <- clean(published[[study_columns[[study]]]])
    add_source(study, published$PA14_ID, call, call == "E", "essentiality", "PA14")
  }
  turner <- read_excel(paths[["turner"]])
  require_columns(turner, "PAO1_ID", "Turner shortlist")
  add_source("Turner_2015", turner$PAO1_ID, "E", TRUE, "essentiality")
  skurnik <- read_excel(paths[["skurnik"]])
  require_columns(skurnik, c("PA14_ID", "Skurnik_2013"), "Skurnik shortlist")
  add_source("Skurnik_2013", skurnik$PA14_ID, skurnik$Skurnik_2013, skurnik$Skurnik_2013 == "E", 
    "essentiality", "PA14")
  pao1161 <- read_excel(paths[["pao1161"]], skip = 20L)
  require_columns(pao1161, c("#Orf", "Name", "Call"), "PAO1161 Gumbel results")
  pao1161$Name <- clean(pao1161$Name)
  add_source("PAO1161_TnSeq_unpublished", pao1161[["#Orf"]], pao1161$Call, pao1161$Call %in% 
    PAO1161_DESIGN_CALLS, "essentiality", mapped_gene_id = pao1161$Name)
  gallagher <- read_excel(paths[["gallagher"]])
  require_columns(gallagher, c("PAO1_ID", "Gallagher_2011"), "Gallagher shortlist")
  add_source("Gallagher_2011", gallagher$PAO1_ID, gallagher$Gallagher_2011, TRUE, "antimicrobial_fitness")
  murray <- read_excel(paths[["murray"]])
  require_columns(murray, c("Condition", "PA14_ID", "logFC", "FDR"), "Murray results")
  murray <- murray[murray$Condition == MURRAY_CONDITION & !is.na(murray$Condition), , drop = FALSE]
  murray_logFC <- suppressWarnings(as.numeric(murray$logFC))
  murray_FDR <- suppressWarnings(as.numeric(murray$FDR))
  selected <- is.finite(murray_logFC) & is.finite(murray_FDR) & murray_FDR < MURRAY_FDR_CUTOFF & 
    abs(murray_logFC) >= MURRAY_MIN_ABS_LOGFC
  add_source("Murray_2015", murray$PA14_ID, ifelse(selected, "AMR", "not_selected"), selected, 
    "antimicrobial_fitness", "PA14", murray_logFC, murray_FDR)
  vitale <- read_excel(paths[["vitale"]])
  require_columns(vitale, c("PAO1_ID", "Vitale_2020"), "Vitale shortlist")
  add_source("Vitale_2020", vitale$PAO1_ID, vitale$Vitale_2020, TRUE, "antimicrobial_fitness")
  evidence <- as.data.frame(data.table::rbindlist(chunks, use.names = TRUE, fill = TRUE))
  evidence <- evidence[order(evidence$source, evidence$PAO1_ID, evidence$source_gene_id, na.last = TRUE), 
    , drop = FALSE]
  essential_sources <- c("Poulsen_2019", "Lee_2015", "Liberati_2006", "Skurnik_2013", "Turner_2015", 
    "PAO1161_TnSeq_unpublished")
  antimicrobial_sources <- c("Gallagher_2011", "Murray_2015", "Vitale_2020")
  for (study in c(essential_sources, antimicrobial_sources)) {
    support <- unique(evidence$PAO1_ID[evidence$source == study & evidence$supports_design & 
      evidence$in_annotation])
    genes[[study]] <- ifelse(genes$PAO1_ID %in% support, if (study %in% essential_sources) 
      "E"
    else "AMR", NA_character_)
  }
  genes$n_essentiality_evidence <- rowSums(!is.na(genes[, essential_sources, drop = FALSE]))
  genes$n_amr_fitness_evidence <- rowSums(!is.na(genes[, antimicrobial_sources, drop = FALSE]))
  category <- ifelse(genes$n_essentiality_evidence > 0L, ifelse(genes$n_amr_fitness_evidence > 
    0L, "E_AMR", "E"), ifelse(genes$n_amr_fitness_evidence > 0L, "AMR", "rest"))
  category[genes$PAO1_ID %in% conditional] <- paste0(category[genes$PAO1_ID %in% conditional], 
    "_core_conditional")
  category[genes$PAO1_ID %in% core] <- paste0(sub("_core_conditional$", "", category[genes$PAO1_ID %in% 
    core]), "_core")
  genes$evidence_category <- category
  genes$essential_category <- ifelse(genes$PAO1_ID %in% core, "E_core", ifelse(genes$PAO1_ID %in% 
    conditional, "E_core_conditional", ifelse(category == "rest", "rest", "E")))
  genes$library_group_planned <- ifelse(genes$PAO1_ID %in% core | genes$n_essentiality_evidence >= 
    PRIMARY_MIN_ESSENTIALITY_SOURCES, "primary_essential", ifelse(genes$n_essentiality_evidence + 
    genes$n_amr_fitness_evidence > 0L, "fitness_associated", "non_essential"))
  summary <- do.call(rbind, lapply(c(essential_sources, antimicrobial_sources), function(study) {
    x <- evidence[evidence$source == study, , drop = FALSE]
    data.frame(source = study, evidence_kind = x$evidence_kind[[1L]], n_PAO1_supported = length(unique(x$PAO1_ID[x$supports_design & 
      x$in_annotation])), n_supported_unmapped_records = sum(x$supports_design & !x$in_annotation))
  }))
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  write_table(genes, "PAO1_library_design_evidence.tsv")
  write_table(evidence, "PAO1_library_source_evidence.tsv.gz")
  write_table(summary, "PAO1_library_evidence_summary.tsv")
  message("Evidence compiled for ", nrow(genes), " PAO1 features from ", nrow(summary), " sources. Output: ", 
    OUTPUT_DIR)
  invisible(genes)
}

run_evidence_collection()

