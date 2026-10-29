#!/usr/bin/env Rscript
############################################################
## Goal: append Pasteur EcoWG1 activity predictions to a guide catalogue.
## Run from the manuscript R project root:
##   source("guide_efficacy_prediction/pasteur/AMG_score_guides_with_Pasteur_model.R")
## Inputs: a guide CSV/TSV and its matching GenBank genome sequence.
## Required columns: locus_tag, sgrna, pam, downstream, strand, position_offset.
## Recommended columns: record_id, feature_start, feature_end.
## The supplied example contains all 46 LESB58 efp guides from the source catalogue.
## For all genes, set INPUT_FILE to the full catalogue; no gene filter is applied.
## Upstream: supply a guide catalogue with the coordinate conventions illustrated
## in the example. This script scores guides; it does not design them.
## Dependency: R package data.table. Model coefficients are embedded below;
## no Python environment or separate model file is needed.
## Outputs: scored table and optional status counts in out_Pasteur_scores/.
## Input rows/columns are retained; failed guides carry an NA score and status.
## Model context: spacer positions 15-20, PAM base 1 and 16 downstream bases.
## A 15-base downstream sequence is extended using the matching genome.
## Coefficient source: Pasteur EcoWG1 reg_coef.pkl,
## MD5 541e9ff987e1ea069c08e54e8faca744.
############################################################

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# USER SETTINGS
# ============================================================

# Small LESB58 example. Replace this path with a full catalogue to score all genes.
INPUT_FILE <- file.path(
  "guide_efficacy_prediction", "pasteur", "input", "LESB58_efp_guides.csv"
)

# Matching reference genome. NULL uses the input `genbank` column to locate it.
GENBANK_FILE <- file.path(
  "guide_efficacy_prediction", "pasteur", "input",
  "Pseudomonas_aeruginosa_LESB58_125.gbk"
)

# Results are saved beside this script in a dedicated output folder.
OUTPUT_FILE <- file.path(
  "guide_efficacy_prediction", "pasteur", "out_Pasteur_scores",
  "LESB58_efp_with_Pasteur_scores.tsv"
)

# Stop when any guide cannot be scored.
# FALSE is recommended because failed rows are retained with a status code.
STRICT_MODE <- FALSE

# Show an RStudio-compatible progress bar.
SHOW_PROGRESS <- TRUE

# Number of guides processed per scoring block.
# Larger values are faster but use more memory.
CHUNK_SIZE <- 50000L

# Use a sequence-search fallback when coordinate lookup fails.
# FALSE is safest for large catalogues: failed coordinate lookups are
# retained with an NA score and a diagnostic status instead of initiating
# a potentially slow whole-genome sequence search.
ALLOW_SEQUENCE_SEARCH_FALLBACK <- FALSE

# Print the exact input row, locus and spacer before each fallback search.
# This is useful only when ALLOW_SEQUENCE_SEARCH_FALLBACK = TRUE.
REPORT_FALLBACK_GUIDES <- TRUE

# Optional explicit exclusions.
#
# Skip by original input row number:
SKIP_INPUT_ROWS <- integer()
# Example:
# SKIP_INPUT_ROWS <- c(154321L, 154322L)
#
# Skip exact locus_tag + sgrna combinations:
SKIP_GUIDES <- data.table(
  locus_tag = character(),
  sgrna = character()
)
# Example:
# SKIP_GUIDES <- data.table(
#   locus_tag = c("PA0001", "PALES_05242"),
#   sgrna = c(
#     "GCACTGCTGCCAAAGTTCCA",
#     "AAAAAAAAAAAAAAAAAAAA"
#   )
# )

# Save a separate status-count table.
WRITE_QC_SUMMARY <- TRUE

# ============================================================
# EMBEDDED PASTEUR MODEL
# ============================================================

PASTEUR_MODEL_NAME <- "Pasteur EcoWG1 reg_coef.pkl"
PASTEUR_MODEL_MD5 <- "541e9ff987e1ea069c08e54e8faca744"
BASE_ORDER <- c("A", "T", "G", "C")

# Exact flattened coefficient order used by the original Python code:
#   all 23 A coefficients,
#   all 23 T coefficients,
#   all 23 G coefficients,
#   all 23 C coefficients.
PASTEUR_COEFFICIENTS <- c(
  -0, 0.12276915331274969, 0.014196762707291635, 0.18069189016342443, 0.059271194019993187, -0.25341302793996195,
  -0.35545337438599711, -0.10888686703080286, -0.084238854396293236, -0.27446400292497086, -0.084398550894202423, -0.00062282263032629223,
  -0.00042593857874643099, -0, 0.10796816854567665, 0.0153438578550551, 0.29680365264124675, 0,
  -0.00035596035126217549, -0, 0, 0.19051425110374787, 0.12775335573103072, 0,
  0, 0, -0.00041902637991167262, 0, 0.3817037559596001, -0.096833843620341586,
  0.069825934675871962, 0.024440278574431713, 0.10551623923402585, 0.059281404307832714, -0.073718314544998234, 0.11317070936371149,
  0.020163150756442805, 0.1061097472690767, 0, 0.043290269525785997, 0, 0.061346473858938648,
  -0.040887360090005764, -0.19332056286569915, -0.09799025814575843, -0.071277442854183554, 0, -0,
  -0, -0.15445413947067085, -0.29853819839469281, -0.5073586937010538, 0.10279400279059485, -0.68678034316579173,
  -0.011155058645955731, -0, -0.0017338821172699253, 0.040655149824689797, 0.056645845103197479, 0,
  -0.096950625330359544, -0, -0.00063684631107936457, 0, 0.048015798916964676, 0,
  0.079895064969003948, 0.14499518518758675, 0.08657428759945146, -0.22638795583239227, -0.042210701019074932, -0.13324981400048325,
  0.03963148074210246, -0, 0.41925114854710671, 0.31320171169911776, 0.24560391988965935, 0.22529551619922544,
  0, 0.050941720477355297, 0.034443813545656246, -0.049282201038877149, -0.10155781146955048, -0.00088148800725202021,
  -0.31602431723684982, -0.17791445341240131, -0.24437144423998292, -0.074764489167238088, 0.072212988911235163, -0,
  -0.00030181175674919616, -0.00020575680798407697
)

if (length(PASTEUR_COEFFICIENTS) != 92L) {
  stop("Internal Pasteur coefficient vector must contain exactly 92 values.")
}

PASTEUR_COEF_MATRIX <- matrix(
  PASTEUR_COEFFICIENTS,
  nrow = 4L,
  ncol = 23L,
  byrow = TRUE,
  dimnames = list(BASE_ORDER, paste0("position_", seq_len(23L)))
)

# ============================================================
# GENERAL HELPERS
# ============================================================

normalize_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  toupper(gsub("\\s+", "", x))
}

valid_dna <- function(x, expected_length = NULL) {
  x <- normalize_text(x)
  valid <- nzchar(x) & grepl("^[ACGT]+$", x)
  
  if (!is.null(expected_length)) {
    valid <- valid & nchar(x) == expected_length
  }
  
  valid
}

reverse_complement <- function(x) {
  x <- normalize_text(x)
  
  if (length(x) == 0L) {
    return(character())
  }
  
  lengths <- nchar(x)
  max_length <- max(lengths, na.rm = TRUE)
  result <- rep("", length(x))
  
  for (position in seq.int(max_length, 1L)) {
    present <- lengths >= position
    if (any(present)) {
      result[present] <- paste0(
        result[present],
        chartr(
          "ACGTN",
          "TGCAN",
          substr(x[present], position, position)
        )
      )
    }
  }
  
  result
}

strip_table_extensions <- function(path) {
  name <- basename(path)
  name <- sub("\\.gz$", "", name, ignore.case = TRUE)
  sub("\\.(csv|tsv|txt)$", "", name, ignore.case = TRUE)
}

make_qc_output <- function(output_file) {
  file.path(
    dirname(output_file),
    paste0(strip_table_extensions(output_file), "_status_summary.tsv")
  )
}

output_separator <- function(path) {
  if (grepl("\\.csv(\\.gz)?$", tolower(path))) "," else "\t"
}

# ============================================================
# LIGHTWEIGHT GENBANK SEQUENCE PARSER
# ============================================================

read_text_lines <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    connection <- gzfile(path, open = "rt")
    on.exit(close(connection), add = TRUE)
    readLines(connection, warn = FALSE)
  } else {
    readLines(path, warn = FALSE)
  }
}

extract_genbank_token <- function(block, prefix) {
  line <- grep(
    paste0("^", prefix, "\\s+"),
    block,
    value = TRUE
  )
  
  if (length(line) == 0L) {
    return(NA_character_)
  }
  
  fields <- strsplit(trimws(line[[1L]]), "\\s+")[[1L]]
  
  if (length(fields) < 2L) {
    NA_character_
  } else {
    fields[[2L]]
  }
}

read_genbank_sequences <- function(genbank_file) {
  lines <- read_text_lines(genbank_file)
  record_ends <- which(trimws(lines) == "//")
  
  if (length(record_ends) == 0L) {
    stop("No complete GenBank records ending with // were found in: ",
         genbank_file)
  }
  
  record_starts <- c(1L, head(record_ends, -1L) + 1L)
  records <- vector("list", length(record_ends))
  
  for (i in seq_along(record_ends)) {
    block <- lines[record_starts[[i]]:record_ends[[i]]]
    
    origin_line <- grep("^ORIGIN\\b", block)
    
    if (length(origin_line) == 0L) {
      stop("GenBank record ", i, " has no ORIGIN section.")
    }
    
    sequence_lines <- block[
      seq.int(origin_line[[1L]] + 1L, length(block) - 1L)
    ]
    
    sequence <- toupper(
      gsub(
        "[^A-Za-z]",
        "",
        paste(sequence_lines, collapse = "")
      )
    )
    sequence <- chartr("U", "T", sequence)
    sequence <- gsub("[^ACGTN]", "N", sequence)
    
    if (!nzchar(sequence)) {
      stop("No nucleotide sequence parsed from GenBank record ", i, ".")
    }
    
    version <- extract_genbank_token(block, "VERSION")
    accession <- extract_genbank_token(block, "ACCESSION")
    locus <- extract_genbank_token(block, "LOCUS")
    
    aliases <- unique(na.omit(c(
      toupper(version),
      toupper(accession),
      toupper(locus),
      toupper(sub("\\.[0-9]+$", "", version))
    )))
    
    canonical_id <- if (!is.na(version)) {
      version
    } else if (!is.na(accession)) {
      accession
    } else {
      locus
    }
    
    records[[i]] <- list(
      canonical_id = canonical_id,
      aliases = aliases,
      sequence = sequence,
      length = nchar(sequence)
    )
  }
  
  sequence_environment <- new.env(parent = emptyenv())
  
  for (record in records) {
    for (alias in record$aliases) {
      if (exists(alias, envir = sequence_environment, inherits = FALSE)) {
        existing <- get(alias, envir = sequence_environment, inherits = FALSE)
        if (!identical(existing, record$sequence)) {
          stop("GenBank alias maps to different sequences: ", alias)
        }
      } else {
        assign(alias, record$sequence, envir = sequence_environment)
      }
    }
  }
  
  list(
    records = records,
    sequence_environment = sequence_environment
  )
}

resolve_record_sequence <- function(record_id, reference) {
  key <- toupper(normalize_text(record_id))
  
  if (nzchar(key) &&
      exists(key, envir = reference$sequence_environment, inherits = FALSE)) {
    return(get(
      key,
      envir = reference$sequence_environment,
      inherits = FALSE
    ))
  }
  
  key_without_version <- sub("\\.[0-9]+$", "", key)
  
  if (nzchar(key_without_version) &&
      exists(
        key_without_version,
        envir = reference$sequence_environment,
        inherits = FALSE
      )) {
    return(get(
      key_without_version,
      envir = reference$sequence_environment,
      inherits = FALSE
    ))
  }
  
  if (length(reference$records) == 1L) {
    return(reference$records[[1L]]$sequence)
  }
  
  NULL
}

resolve_genbank_file <- function(input_dt) {
  if (!is.null(GENBANK_FILE) && nzchar(GENBANK_FILE)) {
    return(normalizePath(GENBANK_FILE, mustWork = TRUE))
  }
  
  if (!"genbank" %in% names(input_dt)) {
    stop(
      "GENBANK_FILE is NULL and the input has no `genbank` column. ",
      "Set GENBANK_FILE to the full GenBank path."
    )
  }
  
  listed <- unique(as.character(input_dt$genbank))
  listed <- listed[!is.na(listed) & nzchar(listed)]
  
  candidates <- unique(c(
    listed[file.exists(listed)],
    file.path(dirname(INPUT_FILE), listed),
    file.path(dirname(INPUT_FILE), basename(listed))
  ))
  candidates <- candidates[file.exists(candidates)]
  
  if (length(candidates) == 0L) {
    stop(
      "Could not locate the GenBank file automatically.\n",
      "The input reports: ", paste(listed, collapse = ", "), "\n",
      "Set GENBANK_FILE explicitly in USER SETTINGS."
    )
  }
  
  if (length(candidates) > 1L) {
    stop(
      "More than one GenBank candidate was found:\n",
      paste(candidates, collapse = "\n"),
      "\nSet GENBANK_FILE explicitly."
    )
  }
  
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

# ============================================================
# 16TH DOWNSTREAM-BASE RECOVERY
# ============================================================

recover_coordinate_context <- function(
    dt,
    reference
) {
  n <- nrow(dt)
  
  result <- data.table(
    extra_base = rep(NA_character_, n),
    source = rep("coordinate_metadata_unavailable", n),
    failure = rep(NA_character_, n)
  )
  
  positions <- suppressWarnings(as.integer(as.numeric(dt$position_offset)))
  strands <- suppressWarnings(as.integer(as.numeric(dt$strand)))
  record_ids <- if ("record_id" %in% names(dt)) {
    normalize_text(dt$record_id)
  } else {
    rep("", n)
  }
  
  valid_metadata <- is.finite(positions) & strands %in% c(-1L, 1L)
  result[!is.finite(positions), failure := "missing_position_offset"]
  result[is.finite(positions) & !strands %in% c(-1L, 1L),
         failure := "invalid_feature_strand"]
  
  keys <- unique(record_ids[valid_metadata])
  if (length(keys) == 0L && any(valid_metadata)) {
    keys <- ""
  }
  
  for (key in keys) {
    idx <- which(valid_metadata & record_ids == key)
    
    if (length(idx) == 0L && key == "") {
      idx <- which(valid_metadata)
    }
    
    genome <- resolve_record_sequence(key, reference)
    
    if (is.null(genome)) {
      result[idx, `:=`(
        failure = "reference_record_ambiguous",
        source = "coordinate_lookup_failed"
      )]
      next
    }
    
    genome_length <- nchar(genome)
    
    plus_idx <- idx[strands[idx] == 1L]
    if (length(plus_idx) > 0L) {
      reverse_genome <- positions[plus_idx] - 4L
      first <- reverse_genome - 16L + 1L
      last <- reverse_genome
      in_bounds <- first >= 1L & last <= genome_length
      
      if (any(in_bounds)) {
        selected <- plus_idx[in_bounds]
        genomic16 <- substring(
          genome,
          first = first[in_bounds],
          last = last[in_bounds]
        )
        expected16 <- reverse_complement(genomic16)
        matches <- (
          nchar(expected16) == 16L &
            valid_dna(expected16, 16L) &
            substr(expected16, 1L, 15L) == dt$downstream_norm[selected]
        )
        
        if (any(matches)) {
          accepted <- selected[matches]
          result[accepted, `:=`(
            extra_base = substr(expected16[matches], 16L, 16L),
            source = "genbank_coordinate_lookup",
            failure = NA_character_
          )]
        }
        
        if (any(!matches)) {
          rejected <- selected[!matches]
          result[rejected, `:=`(
            source = "coordinate_lookup_failed",
            failure = "coordinate_context_mismatch"
          )]
        }
      }
      
      if (any(!in_bounds)) {
        result[plus_idx[!in_bounds], `:=`(
          source = "coordinate_lookup_failed",
          failure = "coordinate_out_of_bounds"
        )]
      }
    }
    
    minus_idx <- idx[strands[idx] == -1L]
    if (length(minus_idx) > 0L) {
      reverse_genome <- positions[minus_idx] + 3L
      first <- reverse_genome + 1L
      last <- reverse_genome + 16L
      in_bounds <- first >= 1L & last <= genome_length
      
      if (any(in_bounds)) {
        selected <- minus_idx[in_bounds]
        expected16 <- substring(
          genome,
          first = first[in_bounds],
          last = last[in_bounds]
        )
        matches <- (
          nchar(expected16) == 16L &
            valid_dna(expected16, 16L) &
            substr(expected16, 1L, 15L) == dt$downstream_norm[selected]
        )
        
        if (any(matches)) {
          accepted <- selected[matches]
          result[accepted, `:=`(
            extra_base = substr(expected16[matches], 16L, 16L),
            source = "genbank_coordinate_lookup",
            failure = NA_character_
          )]
        }
        
        if (any(!matches)) {
          rejected <- selected[!matches]
          result[rejected, `:=`(
            source = "coordinate_lookup_failed",
            failure = "coordinate_context_mismatch"
          )]
        }
      }
      
      if (any(!in_bounds)) {
        result[minus_idx[!in_bounds], `:=`(
          source = "coordinate_lookup_failed",
          failure = "coordinate_out_of_bounds"
        )]
      }
    }
  }
  
  result
}

all_fixed_occurrences <- function(text, pattern) {
  positions <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
  if (length(positions) == 1L && positions[[1L]] == -1L) {
    integer()
  } else {
    as.integer(positions)
  }
}

recover_by_sequence_search <- function(row, reference) {
  prefix <- paste0(
    row$sgrna_norm,
    row$pam_norm,
    row$downstream_norm
  )
  
  if (nchar(prefix) != 38L || !valid_dna(prefix, 38L)) {
    return(list(
      extra_base = NA_character_,
      source = "invalid_search_context",
      match_count = 0L
    ))
  }
  
  preferred_sequence <- resolve_record_sequence(
    if ("record_id" %in% names(row)) row$record_id else "",
    reference
  )
  
  sequences <- if (!is.null(preferred_sequence)) {
    list(preferred_sequence)
  } else {
    lapply(reference$records, `[[`, "sequence")
  }
  
  candidates <- list()
  candidate_n <- 0L
  
  feature_start <- if ("feature_start" %in% names(row)) {
    suppressWarnings(as.integer(as.numeric(row$feature_start)))
  } else {
    NA_integer_
  }
  feature_end <- if ("feature_end" %in% names(row)) {
    suppressWarnings(as.integer(as.numeric(row$feature_end)))
  } else {
    NA_integer_
  }
  
  add_candidates <- function(sequence, orientation) {
    oriented <- if (orientation == "forward") {
      sequence
    } else {
      reverse_complement(sequence)
    }
    
    positions <- all_fixed_occurrences(oriented, prefix)
    
    if (length(positions) == 0L) {
      return(invisible(NULL))
    }
    
    sequence_length <- nchar(sequence)
    prefix_length <- nchar(prefix)
    
    for (start_1based in positions) {
      extra_position <- start_1based + prefix_length
      
      if (extra_position > nchar(oriented)) {
        next
      }
      
      if (orientation == "forward") {
        forward_start_0 <- start_1based - 1L
        forward_end_0 <- forward_start_0 + prefix_length
      } else {
        oriented_start_0 <- start_1based - 1L
        forward_start_0 <- sequence_length -
          (oriented_start_0 + prefix_length)
        forward_end_0 <- sequence_length - oriented_start_0
      }
      
      overlaps_feature <- FALSE
      if (is.finite(feature_start) && is.finite(feature_end)) {
        overlaps_feature <- (
          max(forward_start_0, feature_start) <
            min(forward_end_0, feature_end)
        )
      }
      
      candidate_n <<- candidate_n + 1L
      candidates[[candidate_n]] <<- list(
        extra_base = substr(oriented, extra_position, extra_position),
        overlaps_feature = overlaps_feature
      )
    }
    
    invisible(NULL)
  }
  
  for (sequence in sequences) {
    add_candidates(sequence, "forward")
    add_candidates(sequence, "reverse")
  }
  
  if (length(candidates) == 0L) {
    return(list(
      extra_base = NA_character_,
      source = "context_not_found_in_genbank",
      match_count = 0L
    ))
  }
  
  overlapping <- candidates[
    vapply(candidates, `[[`, logical(1), "overlaps_feature")
  ]
  selected <- if (length(overlapping) > 0L) overlapping else candidates
  bases <- unique(vapply(selected, `[[`, character(1), "extra_base"))
  
  if (length(bases) == 1L) {
    return(list(
      extra_base = bases[[1L]],
      source = if (length(overlapping) > 0L) {
        "genbank_sequence_search_feature_match"
      } else if (length(selected) == 1L) {
        "genbank_sequence_search_unique_match"
      } else {
        "genbank_sequence_search_same_extra_base"
      },
      match_count = length(selected)
    ))
  }
  
  list(
    extra_base = NA_character_,
    source = "ambiguous_genbank_context",
    match_count = length(selected)
  )
}

# ============================================================
# MODEL SCORING
# ============================================================

score_model_contexts <- function(context23) {
  if (length(context23) == 0L) {
    return(numeric())
  }
  
  if (!all(valid_dna(context23, 23L))) {
    stop("score_model_contexts received invalid 23-nt context.")
  }
  
  score <- numeric(length(context23))
  
  for (position in seq_len(23L)) {
    base_index <- match(
      substr(context23, position, position),
      BASE_ORDER
    )
    
    score <- score + PASTEUR_COEF_MATRIX[
      cbind(base_index, rep(position, length(base_index)))
    ]
  }
  
  score
}

score_chunk <- function(chunk, reference) {
  n <- nrow(chunk)
  
  chunk[, `:=`(
    locus_norm = normalize_text(locus_tag),
    sgrna_norm = normalize_text(sgrna),
    pam_norm = normalize_text(pam),
    downstream_norm = normalize_text(downstream)
  )]
  
  status <- rep("scorable", n)
  
  assign_status <- function(mask, label) {
    eligible <- status == "scorable" & mask
    status[eligible] <<- label
  }
  
  if ("skip_pasteur__" %in% names(chunk)) {
    assign_status(
      !is.na(chunk$skip_pasteur__) & chunk$skip_pasteur__,
      "skipped_by_user"
    )
  }
  
  assign_status(!nzchar(chunk$locus_norm), "missing_locus_tag")
  assign_status(!valid_dna(chunk$sgrna_norm, 20L), "invalid_sgrna")
  assign_status(
    !valid_dna(chunk$pam_norm, 3L) |
      substr(chunk$pam_norm, 2L, 3L) != "GG",
    "invalid_pam"
  )
  assign_status(
    !nzchar(chunk$downstream_norm) |
      !grepl("^[ACGT]+$", chunk$downstream_norm),
    "invalid_downstream"
  )
  assign_status(
    nchar(chunk$downstream_norm) < 15L,
    "downstream_shorter_than_15"
  )
  
  downstream16 <- ifelse(
    nchar(chunk$downstream_norm) >= 16L,
    substr(chunk$downstream_norm, 1L, 16L),
    ""
  )
  context_source <- ifelse(
    nchar(chunk$downstream_norm) >= 16L,
    "input_downstream",
    ""
  )
  context_match_count <- integer(n)
  
  lookup_idx <- which(
    status == "scorable" &
      nchar(chunk$downstream_norm) == 15L
  )
  
  if (length(lookup_idx) > 0L) {
    coordinate <- recover_coordinate_context(
      chunk[lookup_idx],
      reference
    )
    
    coordinate_ok <- !is.na(coordinate$extra_base)
    
    if (any(coordinate_ok)) {
      accepted <- lookup_idx[coordinate_ok]
      downstream16[accepted] <- paste0(
        chunk$downstream_norm[accepted],
        coordinate$extra_base[coordinate_ok]
      )
      context_source[accepted] <- coordinate$source[coordinate_ok]
    }
    
    failed_local <- which(!coordinate_ok)
    
    if (length(failed_local) > 0L &&
        isTRUE(ALLOW_SEQUENCE_SEARCH_FALLBACK)) {
      for (j in failed_local) {
        global_index <- lookup_idx[[j]]
        
        if (isTRUE(REPORT_FALLBACK_GUIDES)) {
          message(
            "\nFallback sequence search for input row ",
            chunk$pasteur_input_row__[[global_index]],
            ": locus_tag=",
            as.character(chunk$locus_tag[[global_index]]),
            ", sgrna=",
            as.character(chunk$sgrna[[global_index]])
          )
        }
        
        fallback <- recover_by_sequence_search(
          as.list(chunk[global_index]),
          reference
        )
        
        context_match_count[[global_index]] <- fallback$match_count
        
        if (!is.na(fallback$extra_base) &&
            fallback$extra_base %in% BASE_ORDER) {
          downstream16[[global_index]] <- paste0(
            chunk$downstream_norm[[global_index]],
            fallback$extra_base
          )
          context_source[[global_index]] <- fallback$source
        } else {
          status[[global_index]] <- fallback$source
          context_source[[global_index]] <- fallback$source
        }
      }
    } else if (length(failed_local) > 0L) {
      failed_global <- lookup_idx[failed_local]
      status[failed_global] <- coordinate$failure[failed_local]
      context_source[failed_global] <- coordinate$source[failed_local]
    }
  }
  
  assign_status(
    nchar(downstream16) != 16L,
    "invalid_downstream_length"
  )
  
  context23 <- rep("", n)
  score <- rep(NA_real_, n)
  scorable <- status == "scorable"
  
  if (any(scorable)) {
    context23[scorable] <- paste0(
      substr(chunk$sgrna_norm[scorable], 15L, 20L),
      substr(chunk$pam_norm[scorable], 1L, 1L),
      downstream16[scorable]
    )
    
    invalid_context <- !valid_dna(context23[scorable], 23L)
    
    if (any(invalid_context)) {
      bad_indices <- which(scorable)[invalid_context]
      status[bad_indices] <- "invalid_23nt_model_context"
      context23[bad_indices] <- ""
    }
    
    final_scorable <- status == "scorable"
    score[final_scorable] <- score_model_contexts(
      context23[final_scorable]
    )
  }
  
  data.table(
    pasteur_activity_score = score,
    pasteur_model_context_23nt = context23,
    pasteur_downstream_16nt = downstream16,
    pasteur_score_status = status,
    pasteur_context_source = context_source,
    pasteur_context_match_count = context_match_count
  )
}

assign_global_quartiles <- function(scores) {
  result <- rep(NA_integer_, length(scores))
  finite <- is.finite(scores)
  
  if (!any(finite)) {
    return(result)
  }
  
  cutoffs <- as.numeric(
    quantile(
      scores[finite],
      probs = c(0.25, 0.50, 0.75),
      names = FALSE,
      na.rm = TRUE,
      type = 7
    )
  )
  
  result[finite & scores >= cutoffs[[3L]]] <- 1L
  result[
    finite &
      scores < cutoffs[[3L]] &
      scores >= cutoffs[[2L]]
  ] <- 2L
  result[
    finite &
      scores < cutoffs[[2L]] &
      scores >= cutoffs[[1L]]
  ] <- 3L
  result[finite & scores < cutoffs[[1L]]] <- 4L
  
  result
}

# ============================================================
# MAIN WORKFLOW
# ============================================================

run_pasteur_scoring <- function() {
  input_file <- normalizePath(INPUT_FILE, mustWork = TRUE)
  output_file <- normalizePath(OUTPUT_FILE, mustWork = FALSE)
  
  if (identical(input_file, output_file)) {
    stop("OUTPUT_FILE must not overwrite INPUT_FILE.")
  }
  
  message("Reading guide table:")
  message("  ", input_file)
  
  guides <- fread(input_file, header = TRUE)
  
  if (nrow(guides) == 0L) {
    stop("The input file contains no guide rows.")
  }
  
  required_columns <- c(
    "locus_tag",
    "sgrna",
    "pam",
    "downstream",
    "strand",
    "position_offset"
  )
  
  missing_columns <- setdiff(required_columns, names(guides))
  
  if (length(missing_columns) > 0L) {
    stop(
      "Input file is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      "\nAvailable columns: ",
      paste(names(guides), collapse = ", ")
    )
  }
  
  # Preserve the original input row for diagnostics and optional skipping.
  guides[, pasteur_input_row__ := .I]
  guides[, skip_pasteur__ := FALSE]
  
  if (length(SKIP_INPUT_ROWS) > 0L) {
    skip_rows <- unique(
      suppressWarnings(as.integer(SKIP_INPUT_ROWS))
    )
    skip_rows <- skip_rows[
      is.finite(skip_rows) &
        skip_rows >= 1L &
        skip_rows <= nrow(guides)
    ]
    
    if (length(skip_rows) > 0L) {
      guides[
        pasteur_input_row__ %in% skip_rows,
        skip_pasteur__ := TRUE
      ]
    }
  }
  
  if (nrow(SKIP_GUIDES) > 0L) {
    required_skip_columns <- c("locus_tag", "sgrna")
    missing_skip_columns <- setdiff(
      required_skip_columns,
      names(SKIP_GUIDES)
    )
    
    if (length(missing_skip_columns) > 0L) {
      stop(
        "SKIP_GUIDES is missing: ",
        paste(missing_skip_columns, collapse = ", ")
      )
    }
    
    guide_keys <- paste(
      normalize_text(guides$locus_tag),
      normalize_text(guides$sgrna),
      sep = "|"
    )
    
    skip_keys <- unique(
      paste(
        normalize_text(SKIP_GUIDES$locus_tag),
        normalize_text(SKIP_GUIDES$sgrna),
        sep = "|"
      )
    )
    
    guides[guide_keys %chin% skip_keys, skip_pasteur__ := TRUE]
  }
  
  skipped_n <- guides[, sum(skip_pasteur__)]
  
  if (skipped_n > 0L) {
    message(
      "Explicitly skipping ",
      format(skipped_n, big.mark = ","),
      " guide row(s). They will remain in the output with ",
      "pasteur_score_status = skipped_by_user."
    )
  }
  
  genbank_file <- resolve_genbank_file(guides)
  
  message("Reading GenBank reference:")
  message("  ", genbank_file)
  reference <- read_genbank_sequences(genbank_file)
  
  message("Guides:  ", format(nrow(guides), big.mark = ","))
  message("Records: ", length(reference$records))
  message("Scoring with embedded Pasteur coefficients...")
  
  guide_count <- nrow(guides)
  chunk_starts <- seq.int(1L, guide_count, by = CHUNK_SIZE)
  result_blocks <- vector("list", length(chunk_starts))
  
  progress <- NULL
  if (isTRUE(SHOW_PROGRESS)) {
    progress <- txtProgressBar(
      min = 0L,
      max = guide_count,
      initial = 0L,
      style = 3L
    )
    on.exit(close(progress), add = TRUE)
  }
  
  start_time <- proc.time()[["elapsed"]]
  
  for (block_i in seq_along(chunk_starts)) {
    first <- chunk_starts[[block_i]]
    last <- min(first + CHUNK_SIZE - 1L, guide_count)
    
    result_blocks[[block_i]] <- score_chunk(
      copy(guides[first:last]),
      reference
    )
    
    if (!is.null(progress)) {
      setTxtProgressBar(progress, last)
    }
  }
  
  score_results <- rbindlist(result_blocks, use.names = TRUE)
  scored <- cbind(copy(guides), score_results)
  
  setnames(
    scored,
    "pasteur_input_row__",
    "pasteur_input_row"
  )
  scored[, skip_pasteur__ := NULL]
  
  scored[
    ,
    pasteur_score_quartile_global :=
      assign_global_quartiles(pasteur_activity_score)
  ]
  
  scored[
    ,
    pasteur_score_centered_within_gene := {
      values <- pasteur_activity_score
      finite <- is.finite(values)
      
      if (!any(finite)) {
        rep(NA_real_, .N)
      } else {
        values - median(values[finite])
      }
    },
    by = locus_tag
  ]
  
  scored[
    ,
    pasteur_activity_rank_within_gene := frank(
      -pasteur_activity_score,
      ties.method = "min",
      na.last = "keep"
    ),
    by = locus_tag
  ]
  
  scored[
    ,
    `:=`(
      pasteur_model_name = PASTEUR_MODEL_NAME,
      pasteur_model_md5 = PASTEUR_MODEL_MD5
    )
  ]
  
  appended_columns <- c(
    "pasteur_input_row",
    "pasteur_activity_score",
    "pasteur_score_centered_within_gene",
    "pasteur_score_quartile_global",
    "pasteur_activity_rank_within_gene",
    "pasteur_model_context_23nt",
    "pasteur_downstream_16nt",
    "pasteur_score_status",
    "pasteur_context_source",
    "pasteur_context_match_count",
    "pasteur_model_name",
    "pasteur_model_md5"
  )
  
  original_columns <- setdiff(names(scored), appended_columns)
  setcolorder(scored, c(original_columns, appended_columns))
  
  failed <- scored$pasteur_score_status != "scorable"
  
  if (isTRUE(STRICT_MODE) && any(failed)) {
    failure_summary <- scored[
      failed,
      .N,
      by = pasteur_score_status
    ][order(-N)]
    
    print(failure_summary)
    stop(
      "STRICT_MODE is TRUE and ",
      sum(failed),
      " guide(s) could not be scored."
    )
  }
  
  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  separator <- output_separator(output_file)
  
  fwrite(
    scored,
    file = output_file,
    sep = separator,
    quote = separator == ",",
    na = "NA"
  )
  
  status_summary <- scored[
    ,
    .N,
    by = pasteur_score_status
  ][order(-N)]
  
  if (isTRUE(WRITE_QC_SUMMARY)) {
    fwrite(
      status_summary,
      file = make_qc_output(output_file),
      sep = "\t",
      na = "NA"
    )
  }
  
  elapsed <- proc.time()[["elapsed"]] - start_time
  scored_n <- sum(is.finite(scored$pasteur_activity_score))
  
  cat("\n\nPasteur scoring summary\n")
  cat("-----------------------\n")
  cat("Input rows:    ", format(nrow(scored), big.mark = ","), "\n", sep = "")
  cat("Scored rows:   ", format(scored_n, big.mark = ","), "\n", sep = "")
  cat(
    "Unscored rows: ",
    format(nrow(scored) - scored_n, big.mark = ","),
    "\n",
    sep = ""
  )
  cat(
    "Loci:          ",
    format(uniqueN(scored$locus_tag), big.mark = ","),
    "\n",
    sep = ""
  )
  cat("Elapsed time:  ", sprintf("%.2f seconds", elapsed), "\n", sep = "")
  cat(
    "Rate:          ",
    sprintf("%.1f guides/second", nrow(scored) / max(elapsed, 1e-9)),
    "\n\n",
    sep = ""
  )
  
  print(status_summary)
  
  cat("\nSaved scored table:\n", output_file, "\n", sep = "")
  
  if (isTRUE(WRITE_QC_SUMMARY)) {
    cat(
      "Saved status summary:\n",
      make_qc_output(output_file),
      "\n",
      sep = ""
    )
  }
  
  cat(
    "\nHigher pasteur_activity_score indicates higher predicted ",
    "on-target activity.\n",
    sep = ""
  )
  
  invisible(scored)
}

# Run when this script is sourced in RStudio.
pasteur_scored_guides <- run_pasteur_scoring()
