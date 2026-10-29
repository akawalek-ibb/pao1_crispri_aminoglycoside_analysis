#!/usr/bin/env Rscript
# accC PAO1 growth analysis with tobramycin: one-sample t-tests.
# Run from the RStudio project working directory:
# source("Figure_scripts/Fig4/input_growth_curves/accC/accC_growth_analysis_pvalue_plus_effectsize.R")
#
# Inputs, all under Figure_scripts/Fig4/input_growth_curves/accC/:
# - 2026_03_13_accC_PAO1_TOB_data.xlsx
# - 2026_03_13_accC_PAO1_TOB_samplesheet.xlsx
# - accC_guide_efficiency.xlsx
#
# Follow the flow:
# 1. Read guide sequences/efficacies and match plate wells to the sample sheet.
# 2. Inspect the raw plate, apply the listed exclusions and subtract blanks.
# 3. Draw corrected growth curves and calculate one set of metrics per well.
# 4. Test combined-treatment metrics against the multiplicative expected mean.
# 5. Draw panels A-F and export their underlying data and individual images.
#
# Expected metric = mean(TOB only) * mean(rhamnose only) / mean(untreated).
# The one-sample t-test compares combined-treatment wells against that fixed
# expected mean. T_TEST_ALTERNATIVE = "less" asks whether combined growth is
# lower. A call additionally requires the configured effect cutoff:
# p < 0.05 AND (observed mean - expected mean) / untreated mean <= -0.10.
# DISPLAY_BH_ADJUSTED_P selects the p value used for labels AND calls; the
# default FALSE uses raw p values. BH-adjusted values are also exported.
# The one-sample test treats the estimated expectation as fixed; variation in
# the other three means is not propagated. This script does not bootstrap.
#
# Growth rates are dOD600/dt in OD600 per hour, not specific log-OD growth rates.
# Rate metrics use the selected time/OD window. Capped AUC integrates its
# separately specified time interval. Time points are not test replicates.
# Each well is treated as an independent observation.
#
# Visible labels: accC_g1 for the full-match guide and accC_g1m1,
# accC_g1m2, ... for mismatch variants. Original input keys remain intact.
# Panels: A guides; B six growth curves; C maximum rate; D median rate;
# E capped AUC; F the selected guide's growth and expected OD curve.
# The condition-specific efficacy fits are visual summaries; t-tests use
# observed well metrics, not fitted predictions.
#
# Outputs: combined/individual PNG and PDF figures, plus
# <OUTPUT_BASENAME>_data_and_statistics.xlsx, under
# Figure_scripts/Fig4/input_growth_curves/accC/out_accC_growth_analysis_pvalue_plus_effectsize/.
# These same outputs supply Figure 4 and Supplementary FigS16 source data;
# see README_FigS16_underlying_data.txt in the parent input_growth_curves folder.

# Input format: plate headers follow nine introductory rows; time is in seconds.
# Sample-sheet columns are well, strain and conditions, including blank wells.
# Upstream: experimental workbooks supplied externally; no preprocessing script is required.
#
# R packages used/checked: conflicted, dplyr, drc, ggplot2, patchwork, purrr, readxl,
# scales, tidyr, writexl.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(writexl)
  library(scales)
  library(drc)
  library(conflicted)
})
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflicts_prefer(dplyr::filter)
# ============================================================
# SCRIPT AND INPUT PATHS
# ============================================================
#
# Edit only these four entries when analysing another plate.
# The script does not scan the script directory or try
# alternative filenames.


# Paths are relative to the RStudio project root. source() does not change
# the working directory to the directory containing this script.
INPUT_FOLDER <- "Figure_scripts/Fig4/input_growth_curves/accC"

RAW_DATA_FILENAME <-
  "2026_03_13_accC_PAO1_TOB_data.xlsx"

SAMPLE_SHEET_FILENAME <-
  "2026_03_13_accC_PAO1_TOB_samplesheet.xlsx"

PREDICTED_EFFICIENCY_FILENAME <-
  "accC_guide_efficiency.xlsx"

# Outputs are saved in this
# subfolder inside INPUT_FOLDER.
OUTPUT_SUBFOLDER <-
  "out_accC_growth_analysis_pvalue_plus_effectsize"

RAW_DATA_FILE <- file.path(
  INPUT_FOLDER,
  RAW_DATA_FILENAME
)

SAMPLE_SHEET_FILE <- file.path(
  INPUT_FOLDER,
  SAMPLE_SHEET_FILENAME
)

PREDICTED_EFFICIENCY_FILE <- file.path(
  INPUT_FOLDER,
  PREDICTED_EFFICIENCY_FILENAME
)

required_input_files <- c(
  RAW_DATA_FILE,
  SAMPLE_SHEET_FILE,
  PREDICTED_EFFICIENCY_FILE
)

missing_input_files <- required_input_files[!file.exists(required_input_files)]

if (length(missing_input_files) > 0L) {
  stop(
    "The following explicitly specified input files were not found:\n",
    paste(
      missing_input_files,
      collapse = "\n"
    ),
    "\n\nEdit INPUT_FOLDER and the three filename settings ",
    "at the beginning of the script."
  )
}

OUTPUT_DIR <- file.path(
  INPUT_FOLDER,
  OUTPUT_SUBFOLDER
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# ANALYSIS SETTINGS
# ============================================================

PLOT_TITLE <-
  "PAO1 accC FM + 5 mismatches, effect of TOB in MH"

# Historical guide keys in display order. They match mis_name in the guide
# workbook and strain in the sample sheet after removing PAO1_. Visible
# manuscript labels are applied separately by pretty_guide_name().
ORDER_ON_PLOT <- c(
  "accCm1",
  "accCm2",
  "accCm3",
  "accCm4",
  "accCm5",
  "accCFM"
)

# Preserve this order: untreated, TOB only, rhamnose only, then both.
# The test helper identifies combined-treatment observations from entry 4.
ORDER_CONDITIONS <- c(
  "MH_0",
  "MH_TOB0.4",
  "MH_Rha0.1",
  "MH_Rha0.1_TOB0.4"
)

CONDITION_LABELS <- c(
  "MH_0" = "MH",
  "MH_TOB0.4" = "MH + TOB",
  "MH_Rha0.1" = "MH + rhamnose",
  "MH_Rha0.1_TOB0.4" =
    "MH + rhamnose + TOB"
)

# Growth-curve colours.
CURVE_COLORS <- c(
  "MH_0" = "#000066ff",
  "MH_TOB0.4" = "#E69F00",
  "MH_Rha0.1" = "#ba0000ff",
  "MH_Rha0.1_TOB0.4" = "#7B3294"
)

# Apply the identical condition mapping in panels C-E.
BAR_COLORS <- CURVE_COLORS

# Print the uncorrected plate overview immediately after the
# input data are loaded. Use this first RStudio plot to identify
# wells that should be entered in WELLS_TO_DISCARD.
PRINT_INITIAL_PLATE_QC_IN_RSTUDIO <- TRUE

# Wells can be listed here after quality-control inspection.
# Example:
# WELLS_TO_DISCARD <- c("B04", "F09")
WELLS_TO_DISCARD <- character()

# Display range for growth curves.
PLOT_TIME_MIN <- 0
PLOT_TIME_MAX <- 16
PLOT_OD_MIN <- 0
PLOT_OD_MAX <- 1.25

# Dashed analysis rectangle and growth-rate window. Both time and OD bounds
# are applied before the rate spline is fitted; effective time coverage can
# differ among wells. The existing bounds are retained.
ANALYSIS_TIME_LEFT <- 4
ANALYSIS_TIME_RIGHT <- 10
ANALYSIS_OD_LOWER <- 0
ANALYSIS_OD_UPPER <- 0.8

# AUC uses the complete curve from the experiment start to the
# selected right-hand limit. OD values above the cap are truncated
# rather than removed. This ensures that every condition is
# integrated over the same time interval and limits contributions
# from stationary-phase fluctuations.
AUC_TIME_LEFT <- 0
AUC_TIME_RIGHT <- ANALYSIS_TIME_RIGHT
AUC_OD_LOWER <- 0
AUC_OD_CAP <- ANALYSIS_OD_UPPER

# Blank correction.
BLANK_SMOOTHING_SPAR <- 0.85
ENFORCE_NONDECREASING_BLANK <- TRUE

# Smoothing used only to calculate growth-rate derivatives.
# Plotted biological curves remain unsmoothed.
GROWTH_RATE_SPLINE_SPAR <- 0.70

# The expected curve is not calculated when the untreated mean
# is below this early-growth threshold.
EXPECTED_BASELINE_THRESHOLD <- 0.05

# Panel F matches the configured input guide key. Its visible label follows
# the accC_g1mN convention.
SELECTED_GUIDE_FOR_PANEL_F <- "accC_m3"

# Statistical display.
#
# Compare observed combined-treatment well metrics with the predicted
# multiplicative mean using a one-sample t-test. The predicted mean is treated
# as fixed; uncertainty in that estimated expectation is not propagated.
#
# For the prespecified hypothesis that the observed combined
# response is LOWER than the multiplicative expectation, use:
#   "less"
#
# Use "two.sided" when no direction was prespecified.
T_TEST_ALTERNATIVE <- "less"

# Biological relevance criterion for sensitizing interactions.
# Normalized interaction effect =
#   (observed combined - multiplicative expectation) /
#   untreated value.
# A value <= -0.10 corresponds to an additional reduction
# of at least 10% of the untreated value.
NORMALIZED_INTERACTION_CUTOFF <- -0.10

# FALSE prints raw p-values.
# BH-adjusted values are still calculated and exported.
DISPLAY_BH_ADJUSTED_P <- FALSE

# Predicted-efficacy curve settings.
EFFICACY_PREDICTION_STEP <- 0.01

# Output.
PRINT_ALL_SUBPANELS_IN_RSTUDIO <- TRUE
PRINT_COMBINED_FIGURE_IN_RSTUDIO <- TRUE
OUTPUT_BASENAME <- "Figure_accC_growth_analysis_pvalue_plus_effectsize"

# ============================================================
# GENERAL HELPERS
# ============================================================

normalize_guide_name <- function(x) {
  tolower(
    gsub(
      "[^A-Za-z0-9]",
      "",
      as.character(
        x
      )
    )
  )
}

# Display only: accCFM -> accC_g1; accCm3 -> accC_g1m3.
# Internal keys and joins keep their original values so input records match.
pretty_guide_name <- function(x) {
  x <- as.character(
    x
  )
  
  x <- sub(
    "FM$",
    "_g1",
    x
  )
  
  x <- sub(
    "m([0-9]+)$",
    "_g1m\\1",
    x
  )
  
  x
}

# These helpers format existing p values and stars; they do not run tests
# or change which raw/adjusted p value is used by the display rule.
format_single_p_value <- function(p_value) {
  if (
    length(p_value) != 1L ||
    !is.finite(
      p_value
    )
  ) {
    return(
      "NA"
    )
  }
  
  if (p_value < 0.001) {
    return(
      format(
        p_value,
        scientific = TRUE,
        digits = 2,
        trim = TRUE
      )
    )
  }
  
  formatC(
    p_value,
    format = "f",
    digits = 3
  )
}

format_p_value <- function(p_value) {
  unname(
    vapply(
      as.numeric(
        p_value
      ),
      format_single_p_value,
      character(1)
    )
  )
}

single_p_to_stars <- function(p_value) {
  if (
    length(p_value) != 1L ||
    !is.finite(
      p_value
    )
  ) {
    return(
      ""
    )
  }
  
  if (p_value < 0.001) {
    return(
      "***"
    )
  }
  
  if (p_value < 0.01) {
    return(
      "**"
    )
  }
  
  if (p_value < 0.05) {
    return(
      "*"
    )
  }
  
  "ns"
}

p_to_stars <- function(p_value) {
  unname(
    vapply(
      as.numeric(
        p_value
      ),
      single_p_to_stars,
      character(1)
    )
  )
}

save_ggplot_all_formats <- function(
    plot_object,
    basename,
    width,
    height
) {
  ggsave(
    filename = file.path(
      OUTPUT_DIR,
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
  
  ggsave(
    filename = file.path(
      OUTPUT_DIR,
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

# ============================================================
# LOAD AND STANDARDIZE GUIDE ANNOTATION
# ============================================================

# Read supplied sequences and efficacy predictions. Derive mismatch
# positions directly from the full-match/variant sequences for panel A.
# Warnings flag supplied-position disagreements; input sequences and
# efficacy predictions are retained rather than corrected or reassigned.
guide_annotation_raw <- readxl::read_xlsx(
  PREDICTED_EFFICIENCY_FILE
) %>%
  as_tibble()

required_guide_columns <- c(
  "mis_name",
  "sgRNA",
  "pred_efficacy"
)

missing_guide_columns <- setdiff(
  required_guide_columns,
  names(
    guide_annotation_raw
  )
)

if (length(missing_guide_columns) > 0L) {
  stop(
    "Guide-efficiency table is missing columns: ",
    paste(
      missing_guide_columns,
      collapse = ", "
    )
  )
}

guide_annotation <- guide_annotation_raw %>%
  transmute(
    guide = trimws(
      as.character(
        mis_name
      )
    ),
    
    sequence = toupper(
      trimws(
        as.character(
          sgRNA
        )
      )
    ),
    
    supplied_position = if (
      "position" %in%
      names(
        guide_annotation_raw
      )
    ) {
      suppressWarnings(
        as.integer(
          position
        )
      )
    } else {
      NA_integer_
    },
    
    original_base = if (
      "original" %in%
      names(
        guide_annotation_raw
      )
    ) {
      toupper(
        trimws(
          as.character(
            original
          )
        )
      )
    } else {
      NA_character_
    },
    
    substituted_base = if (
      "subst" %in%
      names(
        guide_annotation_raw
      )
    ) {
      toupper(
        trimws(
          as.character(
            subst
          )
        )
      )
    } else {
      NA_character_
    },
    
    pred_efficacy = suppressWarnings(
      as.numeric(
        pred_efficacy
      )
    )
  ) %>%
  filter(
    !is.na(
      guide
    ),
    guide != "",
    !is.na(
      sequence
    ),
    sequence != ""
  )

full_match_candidates <- guide_annotation %>%
  filter(
    normalize_guide_name(
      guide
    ) ==
      normalize_guide_name(
        "accCFM"
      ) |
      toupper(
        original_base
      ) ==
      "FULL"
  )

if (nrow(full_match_candidates) != 1L) {
  stop(
    "Expected exactly one full-match guide in the ",
    "guide-efficiency table, but found ",
    nrow(
      full_match_candidates
    ),
    "."
  )
}

FULL_MATCH_SEQUENCE <-
  full_match_candidates$sequence[[1L]]

find_mismatch_position <- function(
    sequence,
    full_match_sequence
) {
  sequence_bases <- strsplit(
    sequence,
    ""
  )[[1L]]
  
  full_match_bases <- strsplit(
    full_match_sequence,
    ""
  )[[1L]]
  
  if (
    length(sequence_bases) !=
    length(full_match_bases)
  ) {
    return(
      NA_integer_
    )
  }
  
  mismatch_positions <- which(
    sequence_bases !=
      full_match_bases
  )
  
  if (length(mismatch_positions) == 0L) {
    return(
      NA_integer_
    )
  }
  
  if (length(mismatch_positions) > 1L) {
    warning(
      "Guide ",
      sequence,
      " differs from the full-match sequence at more than ",
      "one position. The first mismatch is displayed."
    )
  }
  
  mismatch_positions[[1L]]
}

guide_annotation <- guide_annotation %>%
  mutate(
    is_full_match =
      sequence ==
      FULL_MATCH_SEQUENCE,
    
    mismatch_position = vapply(
      sequence,
      find_mismatch_position,
      integer(1),
      full_match_sequence =
        FULL_MATCH_SEQUENCE
    ),
    
    full_match_base = vapply(
      mismatch_position,
      function(position_i) {
        if (!is.finite(position_i)) {
          return(NA_character_)
        }
        
        substr(
          FULL_MATCH_SEQUENCE,
          position_i,
          position_i
        )
      },
      character(1)
    ),
    
    mismatch_base = vapply(
      seq_along(sequence),
      function(row_i) {
        position_i <- mismatch_position[[row_i]]
        
        if (!is.finite(position_i)) {
          return(NA_character_)
        }
        
        substr(
          sequence[[row_i]],
          position_i,
          position_i
        )
      },
      character(1)
    ),
    
    # The full-match efficacy is fixed at 1 by definition.
    pred_efficacy = if_else(
      is_full_match,
      1,
      pred_efficacy
    ),
    
    guide = factor(
      guide,
      levels =
        ORDER_ON_PLOT
    ),
    
    guide_display = pretty_guide_name(
      guide
    ),
    
    mismatch_label = case_when(
      is_full_match ~
        "full match",
      
      is.finite(
        mismatch_position
      ) ~
        paste0(
          mismatch_position,
          " ",
          full_match_base,
          "\u2192",
          mismatch_base
        ),
      
      TRUE ~
        "not resolved"
    ),
    
    efficacy_label = if_else(
      is_full_match,
      "1",
      formatC(
        pred_efficacy,
        format = "f",
        digits = 2
      )
    )
  ) %>%
  arrange(
    guide
  )

if (any(is.na(guide_annotation$guide))) {
  stop(
    "The guide-efficiency table contains guide names outside ",
    "ORDER_ON_PLOT."
  )
}

# Warn if the supplied mismatch position disagrees with the
# position calculated directly from the full-match sequence.
position_disagreement <- guide_annotation %>%
  filter(
    !is_full_match,
    is.finite(
      supplied_position
    ),
    is.finite(
      mismatch_position
    ),
    supplied_position !=
      mismatch_position
  )

if (nrow(position_disagreement) > 0L) {
  warning(
    "Supplied and sequence-derived mismatch positions disagree ",
    "for: ",
    paste(
      as.character(
        position_disagreement$guide
      ),
      collapse = ", "
    ),
    ". Sequence-derived positions are displayed."
  )
}

# Make one row per nucleotide for the guide-table drawing. Sequence
# highlighting and labels are visual annotations, not inputs to growth tests.
sequence_long <- purrr::map_dfr(
  seq_len(
    nrow(
      guide_annotation
    )
  ),
  function(row_i) {
    sequence_i <-
      guide_annotation$sequence[[row_i]]
    
    bases_i <- strsplit(
      sequence_i,
      ""
    )[[1L]]
    
    tibble(
      guide =
        guide_annotation$guide[[row_i]],
      guide_display =
        guide_annotation$guide_display[[row_i]],
      y_position =
        nrow(
          guide_annotation
        ) -
        row_i +
        1L,
      base_position =
        seq_along(
          bases_i
        ),
      base =
        bases_i,
      mismatch_position =
        guide_annotation$mismatch_position[[row_i]],
      is_mismatch =
        seq_along(
          bases_i
        ) ==
        guide_annotation$mismatch_position[[row_i]]
    )
  }
)

guide_display_rows <- guide_annotation %>%
  mutate(
    y_position =
      rev(
        seq_len(
          n()
        )
      )
  )

# ============================================================
# PANEL A: GUIDE TABLE
# ============================================================

panel_A_guides <- ggplot() +
  geom_tile(
    data = sequence_long %>%
      filter(
        is_mismatch
      ),
    aes(
      x = base_position,
      y = y_position
    ),
    width = 0.88,
    height = 0.72,
    fill = "#FDE0DD",
    color = "#D73027",
    linewidth = 0.45
  ) +
  geom_text(
    data = sequence_long,
    aes(
      x = base_position,
      y = y_position,
      label = base,
      color = is_mismatch
    ),
    family = "mono",
    size = 3.25
  ) +
  scale_color_manual(
    values = c(
      "FALSE" = "grey35",
      "TRUE" = "#2166AC"
    ),
    guide = "none"
  ) +
  geom_text(
    data = guide_display_rows,
    aes(
      x = -1.0,
      y = y_position,
      label = guide_display
    ),
    hjust = 1,
    family = "mono",
    size = 3.25,
    color = "grey25"
  ) +
  geom_text(
    data = guide_display_rows,
    aes(
      x = 23.2,
      y = y_position,
      label = mismatch_label
    ),
    hjust = 0.5,
    family = "mono",
    size = 3.0,
    color = "grey30"
  ) +
  geom_text(
    data = guide_display_rows,
    aes(
      x = 28.0,
      y = y_position,
      label = efficacy_label
    ),
    hjust = 0.5,
    family = "mono",
    size = 3.25,
    color = "grey25"
  ) +
  annotate(
    "text",
    x = -1.0,
    y = nrow(
      guide_annotation
    ) + 0.9,
    label = "guide",
    hjust = 1,
    family = "mono",
    size = 3.3,
    color = "grey25"
  ) +
  annotate(
    "text",
    x = 10.5,
    y = nrow(
      guide_annotation
    ) + 0.9,
    label = "sequence",
    hjust = 0.5,
    family = "mono",
    size = 3.3,
    color = "grey25"
  ) +
  annotate(
    "text",
    x = 23.2,
    y = nrow(
      guide_annotation
    ) + 0.9,
    label = "mismatch\nposition",
    hjust = 0.5,
    family = "mono",
    size = 3.0,
    lineheight = 0.9,
    color = "grey25"
  ) +
  annotate(
    "text",
    x = 28.0,
    y = nrow(
      guide_annotation
    ) + 0.9,
    label = "predicted\nefficacy",
    hjust = 0.5,
    family = "mono",
    size = 3.0,
    lineheight = 0.9,
    color = "grey25"
  ) +
  coord_cartesian(
    xlim = c(
      -6,
      30
    ),
    ylim = c(
      0.45,
      nrow(
        guide_annotation
      ) + 1.25
    ),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(
      t = 8,
      r = 6,
      b = 8,
      l = 8
    )
  )

# ============================================================
# LOAD AND CLEAN PLATE DATA
# ============================================================

# Skip nine introductory rows to reach the supplied plate-export headers.
# Convert seconds to hours and match strain/condition metadata by well ID.
# The sample sheet must also identify the blank wells.
raw_plate <- readxl::read_xlsx(
  RAW_DATA_FILE,
  skip = 9,
  col_names = TRUE
) %>%
  as_tibble()

sample_sheet <- readxl::read_xlsx(
  SAMPLE_SHEET_FILE,
  col_names = TRUE
) %>%
  as_tibble() %>%
  transmute(
    well = trimws(
      as.character(
        well
      )
    ),
    
    strain = trimws(
      as.character(
        strain
      )
    ),
    
    conditions = trimws(
      as.character(
        conditions
      )
    )
  )

time_columns <- names(raw_plate)[grepl("time", names(raw_plate), ignore.case = TRUE)]

if (length(time_columns) == 0L) {
  stop(
    "No time column was found in the raw plate-data table."
  )
}

TIME_COLUMN <- time_columns[[1L]]

well_columns <- names(raw_plate)[grepl("\\([A-H][0-9]{2}\\)", names(raw_plate))]

if (length(well_columns) == 0L) {
  stop(
    "No well columns were found in the raw plate-data table."
  )
}

plate_long_raw <- raw_plate %>%
  rename(
    time_seconds =
      all_of(
        TIME_COLUMN
      )
  ) %>%
  dplyr::select(
    time_seconds,
    all_of(
      well_columns
    )
  ) %>%
  pivot_longer(
    cols = all_of(
      well_columns
    ),
    names_to = "measurement_name",
    values_to = "OD_raw"
  ) %>%
  mutate(
    time_seconds = suppressWarnings(
      as.numeric(
        time_seconds
      )
    ),
    
    time =
      time_seconds /
      3600,
    
    OD_raw = suppressWarnings(
      as.numeric(
        OD_raw
      )
    ),
    
    well = sub(
      ".*\\(([A-H][0-9]{2})\\).*",
      "\\1",
      measurement_name
    )
  ) %>%
  dplyr::select(
    well,
    time_seconds,
    time,
    OD_raw
  ) %>%
  filter(
    is.finite(
      time
    ),
    is.finite(
      OD_raw
    )
  ) %>%
  left_join(
    sample_sheet,
    by = "well"
  )

missing_metadata <- plate_long_raw %>%
  filter(
    is.na(
      strain
    ) |
      is.na(
        conditions
      )
  ) %>%
  distinct(
    well
  )

if (nrow(missing_metadata) > 0L) {
  stop(
    "Sample-sheet metadata are missing for wells: ",
    paste(
      missing_metadata$well,
      collapse = ", "
    )
  )
}

# ============================================================
# INITIAL RAW-PLATE QC PLOT
# ============================================================
#
# This is intentionally printed before WELLS_TO_DISCARD is
# applied and before blank correction. Each facet represents one
# physical well in its plate position. The plot is not included
# in the final A-F figure. It is saved separately with the individual panels
# at the end of the script.

raw_plate_qc_data <- plate_long_raw %>%
  mutate(
    plate_row = substr(
      well,
      1,
      1
    ),
    
    plate_column = suppressWarnings(
      as.integer(
        substr(
          well,
          2,
          3
        )
      )
    ),
    
    condition_for_fill = if_else(
      is.na(
        conditions
      ) |
        conditions == "",
      "unassigned",
      as.character(
        conditions
      )
    )
  )

raw_plate_qc_background <- raw_plate_qc_data %>%
  distinct(
    well,
    plate_row,
    plate_column,
    condition_for_fill
  )

raw_plate_qc_fill_values <- c(
  CURVE_COLORS,
  "blank" = "#D9D9D9",
  "unassigned" = "white"
)

initial_raw_plate_qc_plot <- ggplot(
  raw_plate_qc_data,
  aes(
    x = time,
    y = OD_raw
  )
) +
  geom_rect(
    data =
      raw_plate_qc_background,
    aes(
      fill =
        condition_for_fill
    ),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.07,
    color = NA,
    inherit.aes = FALSE
  ) +
  geom_line(
    color = "black",
    linewidth = 0.35
  ) +
  scale_fill_manual(
    values =
      raw_plate_qc_fill_values,
    breaks = c(
      ORDER_CONDITIONS,
      "blank",
      "unassigned"
    ),
    labels = c(
      CONDITION_LABELS,
      "Blank",
      "Unassigned"
    ),
    drop = FALSE,
    name = "Condition",
    guide = guide_legend(
      override.aes = list(
        alpha = 0.50
      )
    )
  ) +
  facet_grid(
    plate_row ~
      plate_column,
    drop = FALSE,
    labeller = label_both
  ) +
  coord_cartesian(
    xlim = c(
      PLOT_TIME_MIN,
      PLOT_TIME_MAX
    ),
    ylim = c(
      0,
      PLOT_OD_MAX
    )
  ) +
  labs(
    title =
      "Raw plate overview for well-quality inspection",
    subtitle = paste0(
      "Uncorrected OD values before blank subtraction and before ",
      "WELLS_TO_DISCARD is applied"
    ),
    x =
      "Time [h]",
    y =
      "Raw OD 600nm"
  ) +
  theme_linedraw(
    base_size = 8
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 5
    ),
    
    axis.text.y = element_text(
      size = 5
    ),
    
    strip.text = element_text(
      size = 6
    ),
    
    legend.position = "top",
    
    legend.background =
      element_blank(),
    
    plot.title = element_text(
      face = "bold"
    )
  )

if (
  isTRUE(
    PRINT_INITIAL_PLATE_QC_IN_RSTUDIO
  )
) {
  print(
    initial_raw_plate_qc_plot
  )
}

# Apply exactly the exclusions in WELLS_TO_DISCARD. The preceding QC plot
# still shows those wells; it does not automatically exclude extra curves.
wells_to_discard <- unique(
  trimws(
    as.character(
      WELLS_TO_DISCARD
    )
  )
)

wells_to_discard <- wells_to_discard[!is.na(wells_to_discard) & wells_to_discard != ""]

plate_long_raw <- plate_long_raw %>%
  filter(
    !well %in%
      wells_to_discard
  )

# Average the blank wells at each time and smooth that background. With
# the existing nondecreasing option, prevent the correction from falling.
# Subtract the time-matched correction from each retained nonblank well.
blank_long <- plate_long_raw %>%
  filter(
    conditions ==
      "blank"
  )

if (nrow(blank_long) == 0L) {
  stop(
    "No wells with conditions == 'blank' were found."
  )
}

blank_curve <- blank_long %>%
  group_by(
    time_seconds,
    time
  ) %>%
  summarise(
    blank_mean_raw = mean(
      OD_raw,
      na.rm = TRUE
    ),
    
    blank_sd_raw = sd(
      OD_raw,
      na.rm = TRUE
    ),
    
    n_blank = sum(
      is.finite(
        OD_raw
      )
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    time
  )

if (nrow(blank_curve) < 4L) {
  stop(
    "At least four blank time points are required."
  )
}

blank_fit <- smooth.spline(
  x = blank_curve$time,
  y = blank_curve$blank_mean_raw,
  spar = BLANK_SMOOTHING_SPAR
)

blank_curve$blank_correction <- as.numeric(
  predict(
    blank_fit,
    x = blank_curve$time
  )$y
)

if (
  isTRUE(
    ENFORCE_NONDECREASING_BLANK
  )
) {
  blank_curve$blank_correction <- cummax(
    blank_curve$blank_correction
  )
}

growth_long <- plate_long_raw %>%
  filter(
    conditions !=
      "blank"
  ) %>%
  left_join(
    blank_curve %>%
      dplyr::select(
        time_seconds,
        blank_mean_raw,
        blank_correction
      ),
    by = "time_seconds"
  ) %>%
  mutate(
    OD =
      OD_raw -
      blank_correction,
    
    strain = gsub(
      "^PAO1_",
      "",
      strain
    ),
    
    strain = gsub(
      "^PA14_",
      "",
      strain
    ),
    
    strain = factor(
      strain,
      levels =
        ORDER_ON_PLOT
    ),
    
    conditions = factor(
      conditions,
      levels =
        ORDER_CONDITIONS
    )
  )

if (any(is.na(growth_long$strain))) {
  stop(
    "One or more strain names are not present in ORDER_ON_PLOT."
  )
}

if (any(is.na(growth_long$conditions))) {
  stop(
    "One or more condition names are not present in ",
    "ORDER_CONDITIONS."
  )
}

strain_display_labels <- setNames(
  pretty_guide_name(
    ORDER_ON_PLOT
  ),
  ORDER_ON_PLOT
)

# Summarise corrected curves by guide/condition/time for mean lines and SD
# ribbons. The later metric calculations retain individual wells.
growth_summary <- growth_long %>%
  group_by(
    strain,
    conditions,
    time
  ) %>%
  summarise(
    mean = mean(
      OD,
      na.rm = TRUE
    ),
    
    sd = sd(
      OD,
      na.rm = TRUE
    ),
    
    n = sum(
      is.finite(
        OD
      )
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    high =
      mean +
      sd,
    
    low =
      mean -
      sd
  )

# ============================================================
# EXPECTED MULTIPLICATIVE CURVE
# ============================================================

# Form the expected OD curve from TOB-only, rhamnose-only and untreated
# mean traces. Keep the supplied plotting convention: undefined early-time
# values become zero, a zero floor is applied, then the curve is capped
# at untreated mean OD. Metric-level expected values below are separate
# and do not use this OD-curve cap.
expected_curve_all <- growth_summary %>%
  dplyr::select(
    strain,
    time,
    conditions,
    mean
  ) %>%
  pivot_wider(
    names_from = conditions,
    values_from = mean
  ) %>%
  mutate(
    expected_multiplicative = if_else(
      is.finite(
        MH_0
      ) &
        MH_0 >
        EXPECTED_BASELINE_THRESHOLD &
        is.finite(
          MH_TOB0.4
        ) &
        is.finite(
          MH_Rha0.1
        ),
      
      MH_TOB0.4 *
        MH_Rha0.1 /
        MH_0,
      
      0
    ),
    
    expected_multiplicative = pmax(
      expected_multiplicative,
      0
    ),
    
    expected_multiplicative = pmin(
      expected_multiplicative,
      MH_0
    )
  )

# ============================================================
# PANEL B: SIX GROWTH CURVES
# 
# ============================================================

panel_B_growth_curves <- ggplot(
  growth_summary,
  aes(
    x = time,
    y = mean,
    color = conditions,
    fill = conditions
  )
) +
  geom_rect(
    aes(
      xmin = ANALYSIS_TIME_LEFT,
      xmax = ANALYSIS_TIME_RIGHT,
      ymin = ANALYSIS_OD_LOWER,
      ymax = ANALYSIS_OD_UPPER
    ),
    fill = NA,
    color = "grey50",
    linetype = "dashed",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  stat_summary(
    aes(
      group = conditions
    ),
    fun = mean,
    geom = "line"
  ) +
  scale_colour_manual(
    values =
      CURVE_COLORS,
    breaks =
      ORDER_CONDITIONS,
    labels =
      CONDITION_LABELS
  ) +
  scale_fill_manual(
    values =
      CURVE_COLORS,
    breaks =
      ORDER_CONDITIONS,
    labels =
      CONDITION_LABELS
  ) +
  geom_ribbon(
    aes(
      x = time,
      ymax = high,
      ymin = low,
      group = conditions
    ),
    alpha = 0.35,
    color = NA
  ) +
  geom_line(
    data = expected_curve_all,
    aes(
      x = time,
      y = expected_multiplicative,
      linetype =
        "Multiplicative expectation",
      group = strain
    ),
    color = "black",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  scale_linetype_manual(
    values = c(
      "Multiplicative expectation" =
        "dashed"
    ),
    name = "Legend"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = "solid"
      )
    ),
    linetype = guide_legend(
      override.aes = list(
        color = "black"
      )
    )
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    ),
    legend.position = "right",
    legend.background =
      element_blank()
  ) +
  coord_cartesian(
    xlim = c(
      PLOT_TIME_MIN,
      PLOT_TIME_MAX
    ),
    ylim = c(
      PLOT_OD_MIN,
      PLOT_OD_MAX
    )
  ) +
  xlab(
    "Time [h]"
  ) +
  ylab(
    "OD 600nm"
  ) +
  facet_wrap(
    ~strain,
    ncol = 3,
    labeller = as_labeller(
      strain_display_labels
    )
  ) +
  theme_linedraw() +
  theme(
    panel.grid.minor = element_line(
      color = "grey80",
      linewidth = 0.25
    ),
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    ),
    
    legend.position = "right",
    
    legend.background = element_blank()
  ) +
  ggtitle(
    "PAO1"
  )

# ============================================================
# CURVE METRICS IN THE SELECTED ANALYSIS WINDOW
# ============================================================

calculate_curve_metrics <- function(data) {
  clean_data <- data %>%
    filter(
      is.finite(
        time
      ),
      is.finite(
        OD
      )
    ) %>%
    arrange(
      time
    )
  
  # ----------------------------------------------------------
  # Maximum and median growth rates
  # ----------------------------------------------------------
  #
  # Preserve the selected time and OD window for derivative-based
  # metrics. AUC uses its separately specified integration interval and cap.
  growth_rate_data <- clean_data %>%
    filter(
      time >=
        ANALYSIS_TIME_LEFT,
      time <=
        ANALYSIS_TIME_RIGHT,
      OD >=
        ANALYSIS_OD_LOWER,
      OD <=
        ANALYSIS_OD_UPPER
    )
  
  if (
    nrow(growth_rate_data) >=
    4L &&
    dplyr::n_distinct(
      growth_rate_data$time
    ) >=
    4L
  ) {
    spline_fit <- smooth.spline(
      x =
        growth_rate_data$time,
      y =
        growth_rate_data$OD,
      spar =
        GROWTH_RATE_SPLINE_SPAR
    )
    
    derivative <- as.numeric(
      predict(
        spline_fit,
        x =
          growth_rate_data$time,
        deriv = 1
      )$y
    )
    
    maximum_growth_rate <- max(
      derivative,
      na.rm = TRUE
    )
    
    median_growth_rate <- median(
      derivative,
      na.rm = TRUE
    )
  } else {
    maximum_growth_rate <- NA_real_
    median_growth_rate <- NA_real_
  }
  
  # ----------------------------------------------------------
  # Area under the curve
  # ----------------------------------------------------------
  #
  # Integrate the complete selected AUC time interval. Clip low corrected
  # OD at AUC_OD_LOWER and cap high values at AUC_OD_CAP without dropping
  # their time points. The time interval is the same for every condition.
  auc_data <- clean_data %>%
    filter(
      time >=
        AUC_TIME_LEFT,
      time <=
        AUC_TIME_RIGHT
    ) %>%
    mutate(
      OD_auc = pmin(
        pmax(
          OD,
          AUC_OD_LOWER
        ),
        AUC_OD_CAP
      )
    )
  
  if (
    nrow(auc_data) >=
    2L &&
    dplyr::n_distinct(
      auc_data$time
    ) >=
    2L
  ) {
    auc_value <- sum(
      diff(
        auc_data$time
      ) *
        (
          head(
            auc_data$OD_auc,
            -1L
          ) +
            tail(
              auc_data$OD_auc,
              -1L
            )
        ) /
        2,
      na.rm = TRUE
    )
  } else {
    auc_value <- NA_real_
  }
  
  tibble(
    maximum_growth_rate =
      maximum_growth_rate,
    median_growth_rate =
      median_growth_rate,
    auc =
      auc_value,
    n_growth_rate_points =
      nrow(
        growth_rate_data
      ),
    n_auc_points =
      nrow(
        auc_data
      )
  )
}

# Calculate maximum/median dOD600/dt and capped AUC for each whole-well
# curve. Summarising each well once keeps the well as the test replicate.
curve_metrics_wide <- growth_long %>%
  group_by(
    well,
    strain,
    conditions
  ) %>%
  group_modify(
    ~calculate_curve_metrics(
      .x
    )
  ) %>%
  ungroup()

curve_metrics_long <- curve_metrics_wide %>%
  pivot_longer(
    cols = c(
      maximum_growth_rate,
      median_growth_rate,
      auc
    ),
    names_to = "metric",
    values_to = "metric_value"
  ) %>%
  filter(
    is.finite(
      metric_value
    )
  )

metric_summary <- curve_metrics_long %>%
  group_by(
    metric,
    strain,
    conditions
  ) %>%
  summarise(
    mean_metric = mean(
      metric_value,
      na.rm = TRUE
    ),
    
    sd_metric = sd(
      metric_value,
      na.rm = TRUE
    ),
    
    n_replicates = sum(
      is.finite(
        metric_value
      )
    ),
    
    .groups = "drop"
  )

# ============================================================
# GUIDE EFFICACY MAPPING
# ============================================================

# Attach a guide's supplied efficacy to its well metrics for plotting.
# This does not rescale observed metrics or contribute to the t-test.
guide_efficacy_map <- guide_annotation %>%
  transmute(
    strain = factor(
      as.character(
        guide
      ),
      levels =
        ORDER_ON_PLOT
    ),
    
    pred_efficacy =
      pred_efficacy
  )

metric_efficacy_raw <- curve_metrics_long %>%
  left_join(
    guide_efficacy_map,
    by = "strain"
  ) %>%
  filter(
    is.finite(
      pred_efficacy
    )
  )

metric_efficacy_summary <- metric_summary %>%
  left_join(
    guide_efficacy_map,
    by = "strain"
  ) %>%
  filter(
    is.finite(
      pred_efficacy
    )
  )

# ============================================================
# MULTIPLICATIVE EXPECTATION FOR EACH METRIC
# ============================================================

# For each guide/metric, calculate the fixed reference from the three
# untreated/single-treatment means. The observed combined group does not
# enter this expectation, and fitted efficacy curves are not used here.
metric_expected_effects <- metric_summary %>%
  dplyr::select(
    metric,
    strain,
    conditions,
    mean_metric
  ) %>%
  pivot_wider(
    names_from = conditions,
    values_from = mean_metric
  ) %>%
  mutate(
    predicted_multiplicative =
      MH_TOB0.4 *
      MH_Rha0.1 /
      MH_0
  ) %>%
  dplyr::select(
    metric,
    strain,
    predicted_multiplicative
  )

# ============================================================
# ONE-SAMPLE T-TEST AGAINST THE FIXED EXPECTED MEAN
# ============================================================

# Take the observed combined-treatment well metrics and test them against
# the fixed expected mean. Use the configured direction and require at least two
# observations; failed tests return the configured missing result.
one_sample_multiplicative_test <- function(
    data,
    predicted_value
) {
  observed_combined <- data %>%
    filter(
      as.character(
        conditions
      ) ==
        ORDER_CONDITIONS[[4L]],

      is.finite(
        metric_value
      )
    ) %>%
    pull(
      metric_value
    )

  p_value <- if (
    length(
      observed_combined
    ) >
    1L &&
    is.finite(
      predicted_value
    )
  ) {
    tryCatch(
      stats::t.test(
        observed_combined,
        mu =
          predicted_value,
        alternative =
          T_TEST_ALTERNATIVE
      )$p.value,
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }

  observed_mean <- if (
    length(
      observed_combined
    ) >
    0L
  ) {
    mean(
      observed_combined,
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

  tibble(
    n_combined =
      length(
        observed_combined
      ),

    mean_observed_combined =
      observed_mean,

    sd_observed_combined = if (
      length(
        observed_combined
      ) >
      1L
    ) {
      sd(
        observed_combined,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },

    predicted_multiplicative =
      predicted_value,

    observed_minus_predicted =
      observed_mean -
      predicted_value,

    p_value_one_sample_t =
      p_value
  )
}

# Run the test separately for each guide/metric and normalise the
# observed-minus-expected difference by untreated mean. BH adjustment
# is across guides separately within each metric. Display/call decisions
# use raw p values by default, as set above. A displayed interaction must
# satisfy BOTH p < 0.05 and the existing effect <= -0.10 criterion.
# 'ns' can therefore indicate failure of either requirement.
metric_statistics <- curve_metrics_long %>%
  group_by(
    metric,
    strain
  ) %>%
  group_modify(
    function(data_i, key_i) {
      metric_i <-
        as.character(
          key_i$metric[[1L]]
        )

      strain_i <-
        as.character(
          key_i$strain[[1L]]
        )

      predicted_value <-
        metric_expected_effects %>%
        filter(
          metric ==
            metric_i,

          as.character(
            strain
          ) ==
            strain_i
        ) %>%
        pull(
          predicted_multiplicative
        )

      predicted_value <- if (
        length(
          predicted_value
        ) ==
        1L
      ) {
        predicted_value[[1L]]
      } else {
        NA_real_
      }

      untreated_value <-
        metric_summary %>%
        filter(
          metric ==
            metric_i,

          as.character(
            strain
          ) ==
            strain_i,

          as.character(
            conditions
          ) ==
            "MH_0"
        ) %>%
        pull(
          mean_metric
        )

      untreated_value <- if (
        length(
          untreated_value
        ) ==
        1L
      ) {
        untreated_value[[1L]]
      } else {
        NA_real_
      }

      one_sample_multiplicative_test(
        data =
          data_i,
        predicted_value =
          predicted_value
      ) %>%
        mutate(
          untreated_metric =
            untreated_value,

          normalized_interaction_effect = if (
            is.finite(
              untreated_value
            ) &&
            untreated_value != 0
          ) {
            observed_minus_predicted /
              untreated_value
          } else {
            NA_real_
          }
        )
    }
  ) %>%
  ungroup() %>%
  group_by(
    metric
  ) %>%
  mutate(
    p_value_BH = p.adjust(
      p_value_one_sample_t,
      method = "BH"
    ),

    p_value_display = if (
      isTRUE(
        DISPLAY_BH_ADJUSTED_P
      )
    ) {
      p_value_BH
    } else {
      p_value_one_sample_t
    },

    passes_p_value =
      is.finite(
        p_value_display
      ) &
      p_value_display <
        0.05,

    passes_effect_size =
      is.finite(
        normalized_interaction_effect
      ) &
      normalized_interaction_effect <=
        NORMALIZED_INTERACTION_CUTOFF,

    interaction_significant =
      passes_p_value &
      passes_effect_size,

    significance_label = if_else(
      interaction_significant,
      paste0(
        p_to_stars(
          p_value_display
        ),
        "
(",
        if (
          isTRUE(
            DISPLAY_BH_ADJUSTED_P
          )
        ) {
          "q = "
        } else {
          "p = "
        },
        format_p_value(
          p_value_display
        ),
        ")"
      ),
      "ns"
    )
  ) %>%
  ungroup()

# ============================================================
# CONDITION-SPECIFIC EFFICACY-RESPONSE CURVES
# ============================================================

efficacy_min <- min(
  guide_efficacy_map$pred_efficacy,
  na.rm = TRUE
)

efficacy_max <- max(
  guide_efficacy_map$pred_efficacy,
  na.rm = TRUE
)

EFFICACY_PREDICTION_GRID <- seq(
  from =
    efficacy_min,
  to =
    efficacy_max,
  by =
    EFFICACY_PREDICTION_STEP
)

# Fit a separate four-parameter log-logistic line for each condition and
# metric to visualise the observed guide series. No artificial anchor or
# nonunit weights are supplied. These fits do not determine t-test p values.
fit_metric_efficacy_curve <- function(
    data,
    metric_name,
    condition_name
) {
  data <- data %>%
    filter(
      is.finite(
        pred_efficacy
      ),
      is.finite(
        metric_value
      )
    )
  
  if (
    nrow(
      data
    ) <
    6L ||
    n_distinct(
      data$pred_efficacy
    ) <
    4L
  ) {
    warning(
      "Insufficient observations for ",
      metric_name,
      " / ",
      condition_name,
      "."
    )
    
    return(
      tibble()
    )
  }
  
  fitted_model <- tryCatch(
    drc::drm(
      metric_value ~
        pred_efficacy,
      data = data,
      fct =
        drc::LL2.4()
    ),
    error = function(e) {
      warning(
        "LL2.4 fit failed for ",
        metric_name,
        " / ",
        condition_name,
        ": ",
        conditionMessage(
          e
        )
      )
      
      NULL
    }
  )
  
  if (is.null(fitted_model)) {
    return(
      tibble()
    )
  }
  
  predicted_values <- tryCatch(
    as.numeric(
      predict(
        fitted_model,
        newdata = data.frame(
          pred_efficacy =
            EFFICACY_PREDICTION_GRID
        )
      )
    ),
    error = function(e) {
      warning(
        "Prediction failed for ",
        metric_name,
        " / ",
        condition_name,
        ": ",
        conditionMessage(
          e
        )
      )
      
      rep(
        NA_real_,
        length(
          EFFICACY_PREDICTION_GRID
        )
      )
    }
  )
  
  tibble(
    pred_efficacy =
      EFFICACY_PREDICTION_GRID,
    predicted_metric =
      predicted_values
  ) %>%
    filter(
      is.finite(
        predicted_metric
      )
    )
}

metric_model_predictions <- metric_efficacy_raw %>%
  group_by(
    metric,
    conditions
  ) %>%
  group_modify(
    function(data_i, key_i) {
      fit_metric_efficacy_curve(
        data =
          data_i,
        metric_name =
          as.character(
            key_i$metric[[1L]]
          ),
        condition_name =
          as.character(
            key_i$conditions[[1L]]
          )
      )
    }
  ) %>%
  ungroup()

metric_multiplicative_curves <-
  metric_model_predictions %>%
  dplyr::select(
    metric,
    conditions,
    pred_efficacy,
    predicted_metric
  ) %>%
  pivot_wider(
    names_from =
      conditions,
    values_from =
      predicted_metric
  ) %>%
  mutate(
    predicted_multiplicative =
      MH_TOB0.4 *
      MH_Rha0.1 /
      MH_0
  ) %>%
  dplyr::select(
    metric,
    pred_efficacy,
    predicted_multiplicative
  ) %>%
  filter(
    is.finite(
      predicted_multiplicative
    )
  )

# ============================================================
# SIGNIFICANCE-LABEL POSITIONS
# ============================================================

metric_plot_ranges <- metric_efficacy_summary %>%
  group_by(
    metric
  ) %>%
  summarise(
    metric_min = min(
      mean_metric -
        replace_na(
          sd_metric,
          0
        ),
      na.rm = TRUE
    ),
    
    metric_max = max(
      mean_metric +
        replace_na(
          sd_metric,
          0
        ),
      na.rm = TRUE
    ),
    
    metric_range = max(
      metric_max -
        metric_min,
      0.01
    ),
    
    .groups = "drop"
  )

metric_label_data <- metric_statistics %>%
  left_join(
    guide_efficacy_map,
    by = "strain"
  ) %>%
  left_join(
    metric_plot_ranges,
    by = "metric"
  ) %>%
  mutate(
    y_position =
      mean_observed_combined -
      0.05 *
      metric_range
  )

# ============================================================
# EFFICACY-RELATION PLOT FACTORY
# ============================================================

make_efficacy_metric_plot <- function(
    metric_name,
    panel_title,
    y_axis_label,
    show_legend = TRUE
) {
  summary_i <- metric_efficacy_summary %>%
    filter(
      metric ==
        metric_name
    )
  
  predictions_i <- metric_model_predictions %>%
    filter(
      metric ==
        metric_name
    )
  
  multiplicative_i <- metric_multiplicative_curves %>%
    filter(
      metric ==
        metric_name
    )
  
  labels_i <- metric_label_data %>%
    filter(
      metric ==
        metric_name
    )
  
  ggplot(
    summary_i,
    aes(
      x =
        pred_efficacy,
      y =
        mean_metric,
      color =
        conditions
    )
  ) +
    geom_point() +
    geom_errorbar(
      aes(
        ymin =
          mean_metric -
          sd_metric,
        ymax =
          mean_metric +
          sd_metric
      ),
      width = 0.02
    ) +
    geom_line(
      data =
        predictions_i,
      aes(
        x =
          pred_efficacy,
        y =
          predicted_metric,
        group =
          conditions,
        color =
          conditions
      ),
      inherit.aes = FALSE
    ) +
    geom_line(
      data =
        multiplicative_i,
      aes(
        x =
          pred_efficacy,
        y =
          predicted_multiplicative,
        linetype =
          "Multiplicative expectation"
      ),
      color = "black",
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    geom_text(
      data =
        labels_i,
      aes(
        x =
          pred_efficacy,
        y =
          y_position,
        label =
          significance_label
      ),
      size = 2.6,
      hjust = 0.5,
      angle = 0,
      color = "black",
      inherit.aes = FALSE
    ) +
    labs(
      x =
        "Predicted efficacy",
      y =
        y_axis_label,
      title =
        panel_title,
      color =
        "Condition",
      linetype =
        "Effect type"
    ) +
    scale_x_continuous(
      limits = c(
        efficacy_min - 0.02,
        efficacy_max + 0.02
      ),
      breaks = seq(
        0.2,
        1,
        by = 0.2
      )
    ) +
    scale_colour_manual(
      values =
        BAR_COLORS,
      breaks =
        ORDER_CONDITIONS,
      labels =
        CONDITION_LABELS
    ) +
    scale_linetype_manual(
      values = c(
        "Multiplicative expectation" =
          "solid"
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1
      ),
      
      legend.position = if (
        isTRUE(
          show_legend
        )
      ) {
        "right"
      } else {
        "none"
      },
      
      plot.title = element_text(
        face = "bold"
      )
    ) +
    coord_cartesian(
      clip = "off"
    )
}

# ============================================================
# PANELS C, D AND E
# ============================================================

panel_C_maximum_growth_rate <-
  make_efficacy_metric_plot(
    metric_name =
      "maximum_growth_rate",
    panel_title =
      "Maximum growth rate versus predicted silencing efficacy",
    y_axis_label =
      "Maximum growth rate [OD h\u207b\u00b9]",
    show_legend = TRUE
  )

panel_D_median_growth_rate <-
  make_efficacy_metric_plot(
    metric_name =
      "median_growth_rate",
    panel_title =
      "Median growth rate versus predicted silencing efficacy",
    y_axis_label =
      "Median growth rate [OD h\u207b\u00b9]",
    show_legend = FALSE
  )

panel_E_auc <-
  make_efficacy_metric_plot(
    metric_name =
      "auc",
    panel_title =
      "Capped whole-curve AUC versus predicted silencing efficacy",
    y_axis_label =
      "Capped AUC [OD\u00b7h]",
    show_legend = FALSE
  )

# ============================================================
# PANEL F: SELECTED GUIDE
# ============================================================

# Match the selected input guide key, then use its prepared
# mean growth and expected OD curves. No additional test is run for panel F.
selected_guide_index <- match(
  normalize_guide_name(
    SELECTED_GUIDE_FOR_PANEL_F
  ),
  normalize_guide_name(
    ORDER_ON_PLOT
  )
)

if (is.na(selected_guide_index)) {
  stop(
    "SELECTED_GUIDE_FOR_PANEL_F = '",
    SELECTED_GUIDE_FOR_PANEL_F,
    "' does not match any guide in ORDER_ON_PLOT."
  )
}

selected_guide <- ORDER_ON_PLOT[[selected_guide_index]]

selected_growth_summary <- growth_summary %>%
  filter(
    as.character(
      strain
    ) ==
      selected_guide
  )

selected_expected_curve <- expected_curve_all %>%
  filter(
    as.character(
      strain
    ) ==
      selected_guide
  )

panel_F_selected_guide <- ggplot(
  selected_growth_summary,
  aes(
    x = time,
    y = mean,
    color = conditions,
    fill = conditions
  )
) +
  geom_rect(
    aes(
      xmin = ANALYSIS_TIME_LEFT,
      xmax = ANALYSIS_TIME_RIGHT,
      ymin = ANALYSIS_OD_LOWER,
      ymax = ANALYSIS_OD_UPPER
    ),
    fill = NA,
    color = "grey50",
    linetype = "dashed",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = low,
      ymax = high,
      group = conditions
    ),
    alpha = 0.35,
    color = NA
  ) +
  geom_line(
    aes(
      group = conditions
    ),
    linewidth = 1.2
  ) +
  geom_line(
    data =
      selected_expected_curve,
    aes(
      x = time,
      y =
        expected_multiplicative,
      linetype =
        "Multiplicative expectation"
    ),
    color = "black",
    linewidth = 0.95,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values =
      CURVE_COLORS,
    breaks =
      ORDER_CONDITIONS,
    labels =
      CONDITION_LABELS
  ) +
  scale_fill_manual(
    values =
      CURVE_COLORS,
    breaks =
      ORDER_CONDITIONS,
    labels =
      CONDITION_LABELS
  ) +
  scale_linetype_manual(
    values = c(
      "Multiplicative expectation" =
        "dashed"
    ),
    name = NULL
  ) +
  coord_cartesian(
    xlim = c(
      PLOT_TIME_MIN,
      PLOT_TIME_MAX
    ),
    ylim = c(
      PLOT_OD_MIN,
      PLOT_OD_MAX
    )
  ) +
  labs(
    title = paste0(
      pretty_guide_name(
        selected_guide
      ),
      ": observed and multiplicative expected growth"
    ),
    x = "Time [h]",
    y = "Growth [OD 600 nm]",
    color = "Condition",
    fill = "Condition"
  ) +
  theme_linedraw(
    base_size = 9
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    ),
    
    legend.position = "top",
    
    legend.background =
      element_blank(),
    
    plot.title = element_text(
      face = "bold",
      size = 9.5
    ),
    
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )

# ============================================================
# DISPLAY EACH SUBPANEL IN RSTUDIO
# ============================================================

if (
  isTRUE(
    PRINT_ALL_SUBPANELS_IN_RSTUDIO
  )
) {
  print(
    panel_A_guides
  )
  
  print(
    panel_B_growth_curves
  )
  
  print(
    panel_C_maximum_growth_rate
  )
  
  print(
    panel_D_median_growth_rate
  )
  
  print(
    panel_E_auc
  )
  
  print(
    panel_F_selected_guide
  )
}

# ============================================================
# COMBINED FIGURE A-F
# ============================================================

figure_design <- "
AAABBBBB
AAABBBBB
CCCCDDDD
EEEEFFFF
"

combined_figure <- panel_A_guides +
  panel_B_growth_curves +
  panel_C_maximum_growth_rate +
  panel_D_median_growth_rate +
  panel_E_auc +
  panel_F_selected_guide +
  plot_layout(
    design =
      figure_design
  ) +
  plot_annotation(
    title =
      PLOT_TITLE,
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 14
      ),
      
      plot.tag = element_text(
        face = "bold",
        size = 15
      )
    )
  )

if (
  isTRUE(
    PRINT_COMBINED_FIGURE_IN_RSTUDIO
  )
) {
  print(
    combined_figure
  )
}

save_ggplot_all_formats(
  plot_object =
    combined_figure,
  basename =
    OUTPUT_BASENAME,
  width = 16,
  height = 13.5
)

# ============================================================
# EXCEL EXPORT
# ============================================================

# Main underlying-data output: <OUTPUT_BASENAME>_data_and_statistics.xlsx.
# Its nine sheets retain guides, blank correction, mean growth curves,
# selected expected curve, per-well metrics, summaries/statistics, and
# fitted efficacy/expectation lines. guide_display gives manuscript labels;
# guide/strain retain original matching IDs. Input workbooks are unchanged.
guide_export <- guide_annotation %>%
  transmute(
    guide =
      as.character(
        guide
      ),
    guide_display,
    sequence,
    mismatch_position,
    full_match_base,
    mismatch_base,
    mismatch_label,
    predicted_efficacy =
      pred_efficacy,
    is_full_match
  )

metric_summary_export <- metric_summary %>%
  left_join(
    metric_statistics,
    by = c(
      "metric",
      "strain"
    )
  ) %>%
  arrange(
    metric,
    strain,
    conditions
  )

write_xlsx(
  list(
    guide_annotation =
      as.data.frame(
        guide_export
      ),
    
    blank_correction =
      as.data.frame(
        blank_curve
      ),
    
    growth_curve_summary =
      as.data.frame(
        growth_summary
      ),
    
    selected_guide_expected_curve =
      as.data.frame(
        selected_expected_curve
      ),
    
    per_well_metrics =
      as.data.frame(
        curve_metrics_wide
      ),
    
    metric_summary_and_statistics =
      as.data.frame(
        metric_summary_export
      ),
    
    multiplicative_statistics =
      as.data.frame(
        metric_statistics
      ),
    
    efficacy_model_predictions =
      as.data.frame(
        metric_model_predictions
      ),
    
    fitted_multiplicative_curves =
      as.data.frame(
        metric_multiplicative_curves
      )
  ),
  path = file.path(
    OUTPUT_DIR,
    paste0(
      OUTPUT_BASENAME,
      "_data_and_statistics.xlsx"
    )
  ),
  format_headers = TRUE
)

# ============================================================
# GLOBAL OBJECTS AND STATUS MESSAGE
# ============================================================

assign(
  "accC_initial_raw_plate_qc",
  initial_raw_plate_qc_plot,
  envir = .GlobalEnv
)

assign(
  "accC_panel_A_guides",
  panel_A_guides,
  envir = .GlobalEnv
)

assign(
  "accC_panel_B_growth_curves",
  panel_B_growth_curves,
  envir = .GlobalEnv
)

assign(
  "accC_panel_C_maximum_growth_rate",
  panel_C_maximum_growth_rate,
  envir = .GlobalEnv
)

assign(
  "accC_panel_D_median_growth_rate",
  panel_D_median_growth_rate,
  envir = .GlobalEnv
)

assign(
  "accC_panel_E_auc",
  panel_E_auc,
  envir = .GlobalEnv
)

assign(
  "accC_panel_F_selected_guide",
  panel_F_selected_guide,
  envir = .GlobalEnv
)

assign(
  "accC_combined_figure_A_to_F",
  combined_figure,
  envir = .GlobalEnv
)

message(
  "Completed accC growth analysis.",
  "\nInitial raw-plate QC plot printed first: ",
  PRINT_INITIAL_PLATE_QC_IN_RSTUDIO,
  ".",
  "\nPanels generated: A-F only.",
  "\nAll six subpanels were printed separately in RStudio: ",
  PRINT_ALL_SUBPANELS_IN_RSTUDIO,
  ".",
  "\nPanel B uses the original growth-curve style.",
  "\nPanels C-E show each metric versus predicted efficacy.",
  "\nPanel F guide: ",
  pretty_guide_name(
    selected_guide
  ),
  ".",
  "\nGrowth-rate window: ",
  ANALYSIS_TIME_LEFT,
  "-",
  ANALYSIS_TIME_RIGHT,
  " h and OD ",
  ANALYSIS_OD_LOWER,
  "-",
  ANALYSIS_OD_UPPER,
  ".",
  "\nAUC window: ",
  AUC_TIME_LEFT,
  "-",
  AUC_TIME_RIGHT,
  " h, with OD capped at ",
  AUC_OD_CAP,
  ".",
  "\nDisplayed statistics: one-sample t-test, alternative = '",
  T_TEST_ALTERNATIVE,
  "'.",
  "\nInteraction call additionally requires normalized effect <= ",
  NORMALIZED_INTERACTION_CUTOFF,
  ".",
  "\nAdditive expected-effect curves are not calculated or plotted.",
  "\nFigure formats: PNG and PDF only.",
  "\nOutput directory:\n  ",
  OUTPUT_DIR
)

# ============================================================
# SAVE INDIVIDUAL PLOTS
# PNG AND PDF
# ============================================================

# Export the QC plot and panels A-F at their existing sizes. These images
# and the same underlying-data workbook support Figure 4 and FigS16.
INDIVIDUAL_PLOT_DIR <- file.path(
  OUTPUT_DIR,
  "individual_panels"
)

dir.create(
  INDIVIDUAL_PLOT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

save_individual_plot <- function(
    plot_object,
    filename,
    width,
    height,
    dpi = 300
) {
  ggplot2::ggsave(
    filename = file.path(
      INDIVIDUAL_PLOT_DIR,
      paste0(
        filename,
        ".png"
      )
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  
  ggplot2::ggsave(
    filename = file.path(
      INDIVIDUAL_PLOT_DIR,
      paste0(
        filename,
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

# ------------------------------------------------------------
# Define the size of each exported panel in inches
# ------------------------------------------------------------

save_individual_plot(
  plot_object = initial_raw_plate_qc_plot,
  filename = "QC_raw_plate_overview",
  width = 14,
  height = 9
)

save_individual_plot(
  plot_object = panel_A_guides,
  filename = "Panel_A_guide_sequences",
  width = 6,
  height = 4
)

save_individual_plot(
  plot_object = panel_B_growth_curves,
  filename = "Panel_B_growth_curves",
  width = 7,
  height = 3.7
)

save_individual_plot(
  plot_object = panel_C_maximum_growth_rate,
  filename = "Panel_C_maximum_growth_rate",
  width = 4.5,
  height = 2.25
)

save_individual_plot(
  plot_object = panel_D_median_growth_rate,
  filename = "Panel_D_median_growth_rate",
  width = 2.5,
  height = 2.25
)

save_individual_plot(
  plot_object = panel_E_auc,
  filename = "Panel_E_AUC",
  width = 2.5,
  height = 2.25
)

save_individual_plot(
  plot_object = panel_F_selected_guide,
  filename = paste0(
    "Panel_F_",
    selected_guide,
    "_growth_curve"
  ),
  width = 7,
  height = 5
)

message(
  "Individual plots saved to:\n  ",
  INDIVIDUAL_PLOT_DIR
)
