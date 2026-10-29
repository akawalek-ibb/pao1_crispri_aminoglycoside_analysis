#!/usr/bin/env Rscript
# Figure 4: tobramycin killing assay for ftsZ_g1m4 in PAO1.
# Working directory: manuscript project root containing AMG_cleaned_Git.Rproj.
# Run in R/RStudio: source("Figure_scripts/Fig4/input_killing_assays/ftsZ_g1m4_TOB_PAO1_killing_assay_plotting_and_endpoint_statistics.R")
# Input: Figure_scripts/Fig4/input_killing_assays/ftsZm4_PAO1_TOB.xlsx
# First-sheet columns: culture, time (hours), treatment, CFU (CFU/mL).
# Upstream: supplied experimental CFU workbook; no preprocessing script is required.
# Packages: readxl, dplyr, tidyr, purrr, stringr, ggplot2, patchwork,
# scales, writexl, tibble and readr. No external helper script is sourced.
#
# Flow: read counts -> identify strains/replicates -> summarize time courses
# -> compare endpoint reductions -> draw figures -> export source data.
# Display labels are CTXcA and ftsZ_g1m4. Input culture names are retained;
# CONTROL_REGEX/TEST_REGEX recognize their accepted aliases. The final culture
# letter identifies the biological replicate used for pairing.
#
# Defaults: endpoint 3 h; the main two-sided paired test compares target and
# control replicate-specific log10(CFU without TOB / CFU with TOB) reductions.
# This measures the extra treatment-associated reduction in the target strain,
# not a comparison with starting CFU. Optional within-strain tests are separate.
# Zero CFU observations use half the 500 CFU/mL detection limit for log analysis;
# original CFU observations are also retained in exported data. Copied time-zero
# points are for plotting only and do not enter the endpoint tests.
#
# Outputs: PNG/PDF figures, source TSVs and an XLSX workbook in out/ beside
# INPUT_FILE. Output prefixes come from the input workbook stem.
# Source the complete file; gene/strain labels do not change the test definition.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(writexl)
  library(tibble)
})


# 1. User settings: input location, target label and analysis choices.
# Run from the project root, not from the input_killing_assays folder.
INPUT_FILE <- "Figure_scripts/Fig4/input_killing_assays/ftsZm4_PAO1_TOB.xlsx"

# Save results in the out folder beside the input workbook.
# Each workbook supplies its own output filename prefix, keeping assay outputs separate.
OUTPUT_DIR <- file.path(
  dirname(
    INPUT_FILE
  ),
  "out"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Names used on plots and in output tables. Source culture names remain unchanged.
CONTROL_NAME <- "CTXcA"
TEST_NAME <- "ftsZ_g1m4"

# Recognize the original workbook names and the standardized guide names.
# The final letter of each culture name identifies replicate a, b or c.
CONTROL_REGEX <- "^CTX"
TEST_REGEX <- "^(ftsZm4|ftsZ_g1m4)"

REPLICATE_REGEX <- "([A-Za-z])$"

# Identify tobramycin-treated and untreated cultures from the treatment column.
ANTIMICROBIAL_NAME <- "TOB"
ANTIMICROBIAL_PATTERN <- "TOB"

NO_ANTIMICROBIAL_PATTERN <- paste(
  c(
    "control",
    "no[_ ]?drug",
    "no[_ ]?antibiotic",
    "untreated"
  ),
  collapse = "|"
)

# Leave NA to read the TOB concentration from labels such as TOB_1.5.
ANTIMICROBIAL_CONCENTRATION_OVERRIDE <- NA_character_

ANTIMICROBIAL_CONCENTRATION_UNIT <- "\u00b5g mL\u207b\u00b9"

# Choose the measured time point used for the killing comparison and fold labels.
ENDPOINT_TIME_H <- 3

# Keep zero observations in the input. For logarithmic analysis, replace them
# using the detection-limit settings below; measured CFU is also exported separately.
ZERO_CFU_HANDLING <- "replace"

# Detection limit in CFU/mL. With the fraction below, zero is replaced by 250 CFU/mL.
CFU_DETECTION_LIMIT <- 500

ZERO_REPLACEMENT_FRACTION_OF_LOD <- 0.5

# Use the shared starting count to draw both curves from time zero.
# These copied points are for plotting only and do not enter the endpoint tests.
DUPLICATE_TIME_ZERO_FOR_TREATED_CURVE <- TRUE

# TRUE compares matched experimental runs: control a with target a, and so on.
# FALSE uses an independent Welch test. Keep the setting matching the experiment.
PAIRED_CONTROL_AND_TEST <- TRUE

RUN_STATISTICS <- TRUE

# Optional secondary tests compare TOB with no TOB within each strain.
# The main sensitization test compares the treatment-associated reductions between strains.
CALCULATE_WITHIN_STRAIN_TESTS <- TRUE
SHOW_WITHIN_STRAIN_P_VALUES_ON_BRACKETS <- FALSE

SHOW_BETWEEN_STRAIN_P_VALUE <- TRUE

# Display the measured replicates as well as the group means.
SHOW_INDIVIDUAL_POINTS <- TRUE

SHOW_INDIVIDUAL_PAIR_LINES <- TRUE

PRINT_PLOTS_IN_RSTUDIO <- TRUE

# Optional plot titles. NULL keeps the automatically generated descriptions.
PLOT_TITLE <- NULL
PLOT_SUBTITLE <- NULL

# Plot appearance. These choices do not change the numerical analysis.
BASE_FAMILY <- "sans"
BASE_SIZE <- 10

CONTROL_COLOR <- "#4D4D4D"
TEST_COLOR <- "#0072B2"

NO_DRUG_LINE_TYPE <- "solid"
DRUG_LINE_TYPE <- "22"

NO_DRUG_SHAPE <- 16
DRUG_SHAPE <- 17

INDIVIDUAL_POINT_SIZE <- 1.45
INDIVIDUAL_POINT_ALPHA <- 0.45
MEAN_POINT_SIZE <- 2.5

RIGHT_X_EXPANSION <- 1.15
CONTROL_BRACKET_X_OFFSET <- 0.12
TEST_BRACKET_X_OFFSET <- 0.52
BRACKET_TICK_LENGTH <- 0.07


# Small formatting functions keep p-values and fold labels readable.
format_p_value <- function(
    p_value
) {
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

format_fold <- function(
    fold_value
) {
  if (
    length(fold_value) != 1L ||
    !is.finite(
      fold_value
    )
  ) {
    return(
      "NA"
    )
  }
  
  if (fold_value >= 1000) {
    return(
      scales::comma(
        round(
          fold_value
        ),
        accuracy = 1
      )
    )
  }
  
  if (fold_value >= 100) {
    return(
      formatC(
        fold_value,
        format = "f",
        digits = 0
      )
    )
  }
  
  if (fold_value >= 10) {
    return(
      formatC(
        fold_value,
        format = "f",
        digits = 1
      )
    )
  }
  
  formatC(
    fold_value,
    format = "f",
    digits = 2
  )
}

# Place the CFU axis ticks at powers of ten spanning the observed values.
log10_breaks_for_data <- function(
    x
) {
  x <- x[
    is.finite(
      x
    ) &
      x > 0
  ]
  
  if (length(x) == 0L) {
    return(
      numeric()
    )
  }
  
  lower_power <- floor(
    log10(
      min(
        x
      )
    )
  )
  
  upper_power <- ceiling(
    log10(
      max(
        x
      )
    )
  )
  
  10^seq(
    lower_power,
    upper_power,
    by = 1
  )
}

log10_axis_labels <- function(
    x
) {
  powers <- round(
    log10(
      x
    )
  )
  
  parse(
    text = paste0(
      "10^",
      powers
    )
  )
}

# Read the concentration from the first TOB treatment label, unless overridden above.
extract_antimicrobial_concentration <- function(
    treatment_values
) {
  if (
    !is.na(
      ANTIMICROBIAL_CONCENTRATION_OVERRIDE
    ) &&
    ANTIMICROBIAL_CONCENTRATION_OVERRIDE != ""
  ) {
    return(
      ANTIMICROBIAL_CONCENTRATION_OVERRIDE
    )
  }
  
  treated_values <- unique(
    as.character(
      treatment_values[
        stringr::str_detect(
          as.character(
            treatment_values
          ),
          stringr::regex(
            ANTIMICROBIAL_PATTERN,
            ignore_case = TRUE
          )
        )
      ]
    )
  )
  
  if (length(treated_values) == 0L) {
    return(
      NA_character_
    )
  }
  
  stringr::str_extract(
    treated_values[[1L]],
    "[0-9]+(?:\\.[0-9]+)?"
  )
}

# Run a two-sided t-test after removing unusable values. Return no test
# when fewer than two observations are available or the test cannot be calculated.
safe_t_test <- function(
    x,
    y,
    paired
) {
  x <- as.numeric(
    x
  )
  
  y <- as.numeric(
    y
  )
  
  if (isTRUE(paired)) {
    keep <- stats::complete.cases(
      x,
      y
    )
    
    x <- x[
      keep
    ]
    
    y <- y[
      keep
    ]
  } else {
    x <- x[
      is.finite(
        x
      )
    ]
    
    y <- y[
      is.finite(
        y
      )
    ]
  }
  
  if (
    length(x) < 2L ||
    length(y) < 2L
  ) {
    return(
      NULL
    )
  }
  
  tryCatch(
    stats::t.test(
      x,
      y,
      paired =
        paired,
      alternative =
        "two.sided"
    ),
    error = function(e) {
      NULL
    }
  )
}

# Save the finished plots in the requested image formats and dimensions.
save_plot_outputs <- function(
    plot_object,
    basename,
    width,
    height,
    save_png,
    save_pdf,
    save_tiff,
    png_dpi,
    tiff_dpi
) {
  if (isTRUE(save_png)) {
    ggplot2::ggsave(
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
      dpi = png_dpi,
      bg = "white"
    )
  }
  
  if (isTRUE(save_pdf)) {
    ggplot2::ggsave(
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
  
  if (isTRUE(save_tiff)) {
    ggplot2::ggsave(
      filename = file.path(
        OUTPUT_DIR,
        paste0(
          basename,
          ".tiff"
        )
      ),
      plot = plot_object,
      width = width,
      height = height,
      units = "in",
      dpi = tiff_dpi,
      compression = "lzw",
      bg = "white"
    )
  }
}

# Show plots in RStudio when requested.
print_plot_in_rstudio <- function(
    plot_object
) {
  if (isTRUE(PRINT_PLOTS_IN_RSTUDIO)) {
    print(
      plot_object
    )
  }
}


# 2. Read and check the workbook. Its first sheet must contain
# culture, time, treatment and CFU columns.
if (!file.exists(INPUT_FILE)) {
  stop(
    "Input workbook was not found:\n",
    INPUT_FILE
  )
}

raw_data <- readxl::read_xlsx(
  INPUT_FILE,
  sheet = 1
) %>%
  as_tibble()

required_columns <- c(
  "culture",
  "time",
  "treatment",
  "CFU"
)

missing_columns <- setdiff(
  required_columns,
  names(
    raw_data
  )
)

if (length(missing_columns) > 0L) {
  stop(
    "Input workbook is missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

# Keep complete records, read time and CFU as numbers, and retain the original
# CFU value before applying the zero-count rule.
cleaned_data <- raw_data %>%
  dplyr::select(
    culture,
    time,
    treatment,
    CFU
  ) %>%
  filter(
    !is.na(
      culture
    ),
    !is.na(
      time
    ),
    !is.na(
      treatment
    ),
    !is.na(
      CFU
    )
  ) %>%
  mutate(
    culture = trimws(
      as.character(
        culture
      )
    ),
    
    time = suppressWarnings(
      as.numeric(
        time
      )
    ),
    
    treatment_raw = trimws(
      as.character(
        treatment
      )
    ),
    
    CFU = suppressWarnings(
      as.numeric(
        CFU
      )
    ),
    
    CFU_raw = CFU,
    
    CFU_below_detection =
      CFU == 0
  )

if (
  any(
    !is.finite(
      cleaned_data$time
    )
  )
) {
  stop(
    "Non-numeric time values were detected."
  )
}

if (
  any(
    !is.finite(
      cleaned_data$CFU
    ) |
    cleaned_data$CFU < 0
  )
) {
  stop(
    "CFU values must be finite and cannot be negative."
  )
}

# Count zero observations and apply the chosen replacement only where needed.
n_zero_cfu <- sum(
  cleaned_data$CFU == 0,
  na.rm = TRUE
)

zero_replacement_value <- NA_real_

if (n_zero_cfu > 0L) {
  if (ZERO_CFU_HANDLING == "error") {
    stop(
      n_zero_cfu,
      " zero-CFU observation(s) were detected. ",
      "Set ZERO_CFU_HANDLING <- \"replace\" and provide ",
      "CFU_DETECTION_LIMIT before running the log-scale analysis."
    )
  }
  
  if (ZERO_CFU_HANDLING != "replace") {
    stop(
      "ZERO_CFU_HANDLING must be either \"error\" or ",
      "\"replace\"."
    )
  }
  
  if (
    length(CFU_DETECTION_LIMIT) != 1L ||
    !is.finite(
      CFU_DETECTION_LIMIT
    ) ||
    CFU_DETECTION_LIMIT <= 0
  ) {
    stop(
      "Zero CFU values are present, but CFU_DETECTION_LIMIT ",
      "has not been set to a positive finite value."
    )
  }
  
  if (
    length(ZERO_REPLACEMENT_FRACTION_OF_LOD) != 1L ||
    !is.finite(
      ZERO_REPLACEMENT_FRACTION_OF_LOD
    ) ||
    ZERO_REPLACEMENT_FRACTION_OF_LOD <= 0 ||
    ZERO_REPLACEMENT_FRACTION_OF_LOD > 1
  ) {
    stop(
      "ZERO_REPLACEMENT_FRACTION_OF_LOD must be greater than ",
      "0 and no greater than 1."
    )
  }
  
  zero_replacement_value <-
    CFU_DETECTION_LIMIT *
    ZERO_REPLACEMENT_FRACTION_OF_LOD
  
  cleaned_data <- cleaned_data %>%
    mutate(
      CFU = if_else(
        CFU_below_detection,
        zero_replacement_value,
        CFU
      )
    )
  
  message(
    "Replaced ",
    n_zero_cfu,
    " zero-CFU observation(s) with ",
    format(
      zero_replacement_value,
      scientific = FALSE,
      trim = TRUE
    ),
    " CFU/mL (",
    ZERO_REPLACEMENT_FRACTION_OF_LOD,
    " × detection limit)."
  )
}

# Stop if the same culture, time and treatment has been entered more than once.
duplicate_rows <- cleaned_data %>%
  count(
    culture,
    time,
    treatment_raw,
    name = "n_rows"
  ) %>%
  filter(
    n_rows > 1L
  )

if (nrow(duplicate_rows) > 0L) {
  stop(
    "Duplicate culture/time/treatment rows were detected."
  )
}

# 3. Assign each record to control/target, biological replicate and TOB/no TOB.
DRUG_TREATMENT_LABEL <- paste0(
  "+",
  ANTIMICROBIAL_NAME
)

NO_DRUG_LABEL <- paste0(
  "No ",
  ANTIMICROBIAL_NAME
)

antimicrobial_concentration <- extract_antimicrobial_concentration(
  cleaned_data$treatment_raw
)

if (is.na(antimicrobial_concentration)) {
  stop(
    "The antimicrobial concentration could not be extracted. ",
    "Set ANTIMICROBIAL_CONCENTRATION_OVERRIDE manually."
  )
}

# Extract replicate IDs from culture names so endpoint measurements can be paired.
replicate_match <- stringr::str_match(
  cleaned_data$culture,
  REPLICATE_REGEX
)

if (ncol(replicate_match) < 2L) {
  stop(
    "REPLICATE_REGEX must contain one capture group."
  )
}

cleaned_data <- cleaned_data %>%
  mutate(
    strain = case_when(
      stringr::str_detect(
        culture,
        CONTROL_REGEX
      ) ~
        CONTROL_NAME,
      
      stringr::str_detect(
        culture,
        TEST_REGEX
      ) ~
        TEST_NAME,
      
      TRUE ~
        NA_character_
    ),
    
    replicate =
      replicate_match[, 2],
    
    treatment_key = case_when(
      stringr::str_detect(
        treatment_raw,
        stringr::regex(
          ANTIMICROBIAL_PATTERN,
          ignore_case = TRUE
        )
      ) ~
        "drug",
      
      stringr::str_detect(
        treatment_raw,
        stringr::regex(
          NO_ANTIMICROBIAL_PATTERN,
          ignore_case = TRUE
        )
      ) ~
        "no_drug",
      
      TRUE ~
        NA_character_
    )
  )

# Stop rather than silently dropping cultures or treatments that are not recognized.
unclassified_rows <- cleaned_data %>%
  filter(
    is.na(
      strain
    ) |
      is.na(
        replicate
      ) |
      replicate == "" |
      is.na(
        treatment_key
      )
  )

if (nrow(unclassified_rows) > 0L) {
  stop(
    "Some rows could not be assigned to strain, replicate, ",
    "or treatment. Check CONTROL_REGEX, TEST_REGEX, ",
    "REPLICATE_REGEX and treatment patterns."
  )
}

cleaned_data <- cleaned_data %>%
  mutate(
    strain = factor(
      strain,
      levels = c(
        CONTROL_NAME,
        TEST_NAME
      )
    ),
    
    treatment_key = factor(
      treatment_key,
      levels = c(
        "no_drug",
        "drug"
      )
    ),
    
    treatment_label = factor(
      if_else(
        treatment_key ==
          "no_drug",
        NO_DRUG_LABEL,
        DRUG_TREATMENT_LABEL
      ),
      levels = c(
        NO_DRUG_LABEL,
        DRUG_TREATMENT_LABEL
      )
    )
  )


# 4. Prepare the time-course display. Keep a separate plotting copy so adding
# the shared time-zero baseline cannot alter the measured endpoint data.
plot_data <- cleaned_data

if (isTRUE(DUPLICATE_TIME_ZERO_FOR_TREATED_CURVE)) {
  minimum_time <- min(
    cleaned_data$time,
    na.rm = TRUE
  )
  
  has_treated_time_zero <- cleaned_data %>%
    filter(
      time ==
        minimum_time,
      
      treatment_key ==
        "drug"
    ) %>%
    nrow() >
    0L
  
  if (!has_treated_time_zero) {
    duplicated_time_zero <- cleaned_data %>%
      filter(
        time ==
          minimum_time,
        
        treatment_key ==
          "no_drug"
      ) %>%
      mutate(
        treatment_key = factor(
          "drug",
          levels = c(
            "no_drug",
            "drug"
          )
        ),
        
        treatment_label = factor(
          DRUG_TREATMENT_LABEL,
          levels = c(
            NO_DRUG_LABEL,
            DRUG_TREATMENT_LABEL
          )
        ),
        
        treatment_raw = paste0(
          ANTIMICROBIAL_NAME,
          "_",
          antimicrobial_concentration,
          "_time0_plot_copy"
        )
      )
    
    plot_data <- bind_rows(
      cleaned_data,
      duplicated_time_zero
    )
  }
}

# Calculate arithmetic mean CFU and SD for each strain, treatment and time.
# A lower error bar that cannot be shown on a log axis is kept above zero for display.
timecourse_summary <- plot_data %>%
  group_by(
    strain,
    treatment_key,
    treatment_label,
    time
  ) %>%
  summarise(
    n = n(),
    
    mean_CFU = mean(
      CFU
    ),
    
    sd_CFU = stats::sd(
      CFU
    ),
    
    minimum_observed_CFU = min(
      CFU
    ),
    
    lower_CFU = pmax(
      mean_CFU -
        sd_CFU,
      minimum_observed_CFU /
        10
    ),
    
    upper_CFU =
      mean_CFU +
      sd_CFU,
    
    .groups = "drop"
  )


# 5. Select the endpoint from the measured data, not the plotting copy.
endpoint_long <- cleaned_data %>%
  filter(
    time ==
      ENDPOINT_TIME_H
  ) %>%
  dplyr::select(
    strain,
    replicate,
    treatment_key,
    treatment_label,
    CFU
  )

if (nrow(endpoint_long) == 0L) {
  stop(
    "No measurements were found at ENDPOINT_TIME_H = ",
    ENDPOINT_TIME_H,
    "."
  )
}

endpoint_duplicates <- endpoint_long %>%
  count(
    strain,
    replicate,
    treatment_key,
    name = "n_rows"
  ) %>%
  filter(
    n_rows > 1L
  )

if (nrow(endpoint_duplicates) > 0L) {
  stop(
    "Duplicate strain/replicate/treatment observations were ",
    "detected at the selected endpoint."
  )
}

# Put each replicate's treated and untreated CFU side by side.
# Each replicate must have both measurements at the same endpoint.
endpoint_reductions <- endpoint_long %>%
  dplyr::select(
    strain,
    replicate,
    treatment_key,
    CFU
  ) %>%
  pivot_wider(
    names_from =
      treatment_key,
    values_from =
      CFU,
    names_prefix =
      "CFU_"
  )

if (
  !"CFU_no_drug" %in%
  names(
    endpoint_reductions
  ) ||
  !"CFU_drug" %in%
  names(
    endpoint_reductions
  )
) {
  stop(
    "Both no-antimicrobial and antimicrobial-treated endpoint ",
    "measurements are required."
  )
}

incomplete_endpoint_pairs <- endpoint_reductions %>%
  filter(
    is.na(
      CFU_no_drug
    ) |
      is.na(
        CFU_drug
      )
  )

if (nrow(incomplete_endpoint_pairs) > 0L) {
  stop(
    "Some strain/replicate combinations do not contain both ",
    "endpoint treatments."
  )
}

# For each replicate, calculate log10(CFU without TOB) minus log10(CFU with TOB).
# The corresponding fold decrease is CFU without TOB divided by CFU with TOB.
# This compares treatments at the same time; it is not normalized to starting CFU.
endpoint_reductions <- endpoint_reductions %>%
  mutate(
    log10_CFU_no_drug = log10(
      CFU_no_drug
    ),
    
    log10_CFU_drug = log10(
      CFU_drug
    ),
    
    delta_log10_CFU =
      log10_CFU_no_drug -
      log10_CFU_drug,
    
    fold_decrease =
      10^delta_log10_CFU
  )

# Average the replicate log reductions, then convert back to a geometric mean
# fold decrease. These are the fold labels displayed beside the curves.
fold_summary <- endpoint_reductions %>%
  group_by(
    strain
  ) %>%
  summarise(
    n = n(),
    
    mean_delta_log10_CFU = mean(
      delta_log10_CFU
    ),
    
    sd_delta_log10_CFU = stats::sd(
      delta_log10_CFU
    ),
    
    geometric_mean_fold_decrease =
      10^mean_delta_log10_CFU,
    
    .groups = "drop"
  )


# 6. Calculate the optional within-strain paired tests on endpoint log10 CFU.
within_strain_tests <- purrr::map_dfr(
  levels(
    cleaned_data$strain
  ),
  function(strain_i) {
    data_i <- endpoint_reductions %>%
      filter(
        strain ==
          strain_i
      )
    
    test_i <- NULL
    
    if (isTRUE(CALCULATE_WITHIN_STRAIN_TESTS)) {
      test_i <- safe_t_test(
        data_i$log10_CFU_no_drug,
        data_i$log10_CFU_drug,
        paired = TRUE
      )
    }
    
    tibble(
      strain =
        strain_i,
      
      n =
        nrow(
          data_i
        ),
      
      test =
        "paired two-sided t-test on log10 CFU",
      
      statistic = if (
        is.null(
          test_i
        )
      ) {
        NA_real_
      } else {
        unname(
          test_i$statistic
        )
      },
      
      degrees_freedom = if (
        is.null(
          test_i
        )
      ) {
        NA_real_
      } else {
        unname(
          test_i$parameter
        )
      },
      
      p_value = if (
        is.null(
          test_i
        )
      ) {
        NA_real_
      } else {
        test_i$p.value
      }
    )
  }
)

# Prepare the main comparison of killing in the target and control strains.
between_strain_test <- tibble(
  control =
    CONTROL_NAME,
  
  test_strain =
    TEST_NAME,
  
  endpoint_time_h =
    ENDPOINT_TIME_H,
  
  paired =
    PAIRED_CONTROL_AND_TEST,
  
  n_control =
    NA_integer_,
  
  n_test =
    NA_integer_,
  
  mean_delta_control =
    NA_real_,
  
  mean_delta_test =
    NA_real_,
  
  difference_test_minus_control =
    NA_real_,
  
  statistic =
    NA_real_,
  
  degrees_freedom =
    NA_real_,
  
  p_value =
    NA_real_
)

# Compare the replicate-specific log reductions between target and control.
# A larger reduction for the target indicates greater TOB activity under these conditions.
if (isTRUE(RUN_STATISTICS)) {
  control_reductions <- endpoint_reductions %>%
    filter(
      strain ==
        CONTROL_NAME
    ) %>%
    arrange(
      replicate
    )
  
  test_reductions <- endpoint_reductions %>%
    filter(
      strain ==
        TEST_NAME
    ) %>%
    arrange(
      replicate
    )
  
  # For the paired test, match target and control reductions by replicate ID.
  if (isTRUE(PAIRED_CONTROL_AND_TEST)) {
    matched_reductions <- control_reductions %>%
      dplyr::select(
        replicate,
        delta_control =
          delta_log10_CFU
      ) %>%
      inner_join(
        test_reductions %>%
          dplyr::select(
            replicate,
            delta_test =
              delta_log10_CFU
          ),
        by =
          "replicate"
      )
    
    control_delta <-
      matched_reductions$delta_control
    
    test_delta <-
      matched_reductions$delta_test
  } else {
    control_delta <-
      control_reductions$delta_log10_CFU
    
    test_delta <-
      test_reductions$delta_log10_CFU
  }
  
  between_test_result <- safe_t_test(
    test_delta,
    control_delta,
    paired =
      PAIRED_CONTROL_AND_TEST
  )
  
  between_strain_test <- between_strain_test %>%
    mutate(
      n_control =
        length(
          control_delta
        ),
      
      n_test =
        length(
          test_delta
        ),
      
      mean_delta_control = mean(
        control_delta,
        na.rm = TRUE
      ),
      
      mean_delta_test = mean(
        test_delta,
        na.rm = TRUE
      ),
      
      difference_test_minus_control =
        mean_delta_test -
        mean_delta_control,
      
      statistic = if (
        is.null(
          between_test_result
        )
      ) {
        NA_real_
      } else {
        unname(
          between_test_result$statistic
        )
      },
      
      degrees_freedom = if (
        is.null(
          between_test_result
        )
      ) {
        NA_real_
      } else {
        unname(
          between_test_result$parameter
        )
      },
      
      p_value = if (
        is.null(
          between_test_result
        )
      ) {
        NA_real_
      } else {
        between_test_result$p.value
      }
    )
}


# 7. Position the endpoint brackets using the arithmetic mean CFU values.
# Their labels come from the geometric mean fold decreases calculated above.
endpoint_mean_wide <- timecourse_summary %>%
  filter(
    time ==
      ENDPOINT_TIME_H
  ) %>%
  dplyr::select(
    strain,
    treatment_key,
    mean_CFU
  ) %>%
  pivot_wider(
    names_from =
      treatment_key,
    values_from =
      mean_CFU,
    names_prefix =
      "mean_"
  )

if (
  !"mean_no_drug" %in%
  names(
    endpoint_mean_wide
  ) ||
  !"mean_drug" %in%
  names(
    endpoint_mean_wide
  )
) {
  stop(
    "Could not construct endpoint bracket coordinates."
  )
}

endpoint_brackets <- endpoint_mean_wide %>%
  left_join(
    fold_summary,
    by =
      "strain"
  ) %>%
  left_join(
    within_strain_tests %>%
      dplyr::select(
        strain,
        within_p_value =
          p_value
      ),
    by =
      "strain"
  ) %>%
  mutate(
    bracket_x = if_else(
      strain ==
        CONTROL_NAME,
      ENDPOINT_TIME_H +
        CONTROL_BRACKET_X_OFFSET,
      ENDPOINT_TIME_H +
        TEST_BRACKET_X_OFFSET
    ),
    
    y_top = pmax(
      mean_no_drug,
      mean_drug
    ),
    
    y_bottom = pmin(
      mean_no_drug,
      mean_drug
    ),
    
    y_label = sqrt(
      y_top *
        y_bottom
    ),
    
    fold_text = paste0(
      vapply(
        geometric_mean_fold_decrease,
        format_fold,
        character(1)
      ),
      "\u00d7"
    ),
    
    bracket_label = if (
      isTRUE(
        SHOW_WITHIN_STRAIN_P_VALUES_ON_BRACKETS
      )
    ) {
      paste0(
        fold_text,
        "\nwithin-strain p = ",
        vapply(
          within_p_value,
          format_p_value,
          character(1)
        )
      )
    } else {
      fold_text
    }
  )

between_test_label <- paste0(
  TEST_NAME,
  " vs ",
  CONTROL_NAME,
  " \u0394log\u2081\u2080 CFU",
  "\n",
  "p = ",
  format_p_value(
    between_strain_test$p_value[[1L]]
  )
)


# 8. Draw the main time course with mean CFU, SD, individual measurements,
# endpoint fold labels and the main between-strain p-value.
strain_colors <- c(
  CONTROL_COLOR,
  TEST_COLOR
)

names(strain_colors) <- c(
  CONTROL_NAME,
  TEST_NAME
)

treatment_line_types <- c(
  NO_DRUG_LINE_TYPE,
  DRUG_LINE_TYPE
)

names(treatment_line_types) <- c(
  NO_DRUG_LABEL,
  DRUG_TREATMENT_LABEL
)

treatment_shapes <- c(
  NO_DRUG_SHAPE,
  DRUG_SHAPE
)

names(treatment_shapes) <- c(
  NO_DRUG_LABEL,
  DRUG_TREATMENT_LABEL
)

automatic_title <- paste0(
  ANTIMICROBIAL_NAME,
  " killing assay"
)

automatic_subtitle <- paste0(
  CONTROL_NAME,
  " and ",
  TEST_NAME,
  "; ",
  ANTIMICROBIAL_NAME,
  " ",
  antimicrobial_concentration,
  " ",
  ANTIMICROBIAL_CONCENTRATION_UNIT
)

title_to_use <- if (
  is.null(
    PLOT_TITLE
  )
) {
  automatic_title
} else {
  PLOT_TITLE
}

subtitle_to_use <- if (
  is.null(
    PLOT_SUBTITLE
  )
) {
  automatic_subtitle
} else {
  PLOT_SUBTITLE
}

timecourse_y_breaks <- log10_breaks_for_data(
  c(
    timecourse_summary$lower_CFU,
    timecourse_summary$upper_CFU
  )
)

maximum_time <- max(
  timecourse_summary$time,
  na.rm = TRUE
)

minimum_time <- min(
  timecourse_summary$time,
  na.rm = TRUE
)

final_killing_plot <- ggplot(
  timecourse_summary,
  aes(
    x =
      time,
    y =
      mean_CFU,
    color =
      strain,
    linetype =
      treatment_label,
    shape =
      treatment_label,
    group =
      interaction(
        strain,
        treatment_label
      )
  )
) +
  geom_errorbar(
    aes(
      ymin =
        lower_CFU,
      ymax =
        upper_CFU
    ),
    width = 0.07,
    linewidth = 0.45
  ) +
  geom_line(
    linewidth = 0.8
  )

if (isTRUE(SHOW_INDIVIDUAL_POINTS)) {
  final_killing_plot <- final_killing_plot +
    geom_point(
      data =
        cleaned_data,
      aes(
        x =
          time,
        y =
          CFU,
        color =
          strain,
        shape =
          treatment_label
      ),
      position = position_jitter(
        width = 0.025,
        height = 0,
        seed = 1
      ),
      size =
        INDIVIDUAL_POINT_SIZE,
      alpha =
        INDIVIDUAL_POINT_ALPHA,
      inherit.aes = FALSE
    )
}

final_killing_plot <- final_killing_plot +
  geom_point(
    size =
      MEAN_POINT_SIZE,
    stroke = 0.35
  ) +
  geom_segment(
    data =
      endpoint_brackets,
    aes(
      x =
        bracket_x,
      xend =
        bracket_x,
      y =
        y_bottom,
      yend =
        y_top,
      color =
        strain
    ),
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data =
      endpoint_brackets,
    aes(
      x =
        bracket_x -
        BRACKET_TICK_LENGTH,
      xend =
        bracket_x,
      y =
        y_top,
      yend =
        y_top,
      color =
        strain
    ),
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data =
      endpoint_brackets,
    aes(
      x =
        bracket_x -
        BRACKET_TICK_LENGTH,
      xend =
        bracket_x,
      y =
        y_bottom,
      yend =
        y_bottom,
      color =
        strain
    ),
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  geom_text(
    data =
      endpoint_brackets,
    aes(
      x =
        bracket_x +
        0.035,
      y =
        y_label,
      label =
        bracket_label,
      color =
        strain
    ),
    hjust = 0,
    vjust = 0.5,
    size = 3.0,
    lineheight = 1.0,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  scale_y_log10(
    breaks =
      timecourse_y_breaks,
    labels =
      log10_axis_labels,
    expand = expansion(
      mult = c(
        0.07,
        0.18
      )
    )
  ) +
  scale_x_continuous(
    breaks = sort(
      unique(
        timecourse_summary$time
      )
    ),
    limits = c(
      minimum_time,
      maximum_time +
        RIGHT_X_EXPANSION
    ),
    expand = expansion(
      mult = c(
        0.025,
        0
      )
    )
  ) +
  scale_color_manual(
    values =
      strain_colors,
    drop = FALSE
  ) +
  scale_linetype_manual(
    values =
      treatment_line_types,
    drop = FALSE
  ) +
  scale_shape_manual(
    values =
      treatment_shapes,
    drop = FALSE
  ) +
  labs(
    title =
      title_to_use,
    
    subtitle =
      subtitle_to_use,
    
    x =
      "Time [h]",
    
    y =
      expression(
        "Viable cells [CFU mL"^{-1}*"]"
      ),
    
    color =
      "Strain",
    
    linetype =
      "Treatment",
    
    shape =
      "Treatment"
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  theme_bw(
    base_size =
      BASE_SIZE,
    base_family =
      BASE_FAMILY
  ) +
  theme(
    panel.grid.major = element_line(
      color = "grey88",
      linewidth = 0.3
    ),
    
    panel.grid.minor = element_blank(),
    
    plot.title = element_text(
      face = "bold"
    ),
    
    legend.position =
      "top",
    
    legend.box =
      "vertical",
    
    legend.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    ),
    
    plot.margin = margin(
      t = 8,
      r = 26,
      b = 8,
      l = 8
    )
  )

if (
  isTRUE(
    SHOW_BETWEEN_STRAIN_P_VALUE
  )
) {
  final_killing_plot <- final_killing_plot +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label =
        between_test_label,
      hjust = -0.03,
      vjust = 1.25,
      size = 3.0
    )
}


# 9. Draw a second view of the endpoint measurements, with lines linking each
# replicate's TOB and no-TOB observations.
endpoint_plot_data <- endpoint_long %>%
  mutate(
    x_group = case_when(
      strain ==
        CONTROL_NAME &
        treatment_key ==
        "no_drug" ~
        1,
      
      strain ==
        CONTROL_NAME &
        treatment_key ==
        "drug" ~
        2,
      
      strain ==
        TEST_NAME &
        treatment_key ==
        "no_drug" ~
        4,
      
      strain ==
        TEST_NAME &
        treatment_key ==
        "drug" ~
        5
    )
  )

endpoint_summary <- endpoint_plot_data %>%
  group_by(
    strain,
    treatment_key,
    treatment_label,
    x_group
  ) %>%
  summarise(
    n = n(),
    
    mean_CFU = mean(
      CFU
    ),
    
    sd_CFU = stats::sd(
      CFU
    ),
    
    minimum_observed_CFU = min(
      CFU
    ),
    
    lower_CFU = pmax(
      mean_CFU -
        sd_CFU,
      minimum_observed_CFU /
        10
    ),
    
    upper_CFU =
      mean_CFU +
      sd_CFU,
    
    .groups = "drop"
  )

endpoint_plot_brackets <- endpoint_brackets %>%
  mutate(
    x_start = if_else(
      strain ==
        CONTROL_NAME,
      1,
      4
    ),
    
    x_end = if_else(
      strain ==
        CONTROL_NAME,
      2,
      5
    ),
    
    endpoint_bracket_y =
      y_top *
      2.5,
    
    endpoint_tick_bottom =
      endpoint_bracket_y /
      1.3,
    
    endpoint_label_y =
      endpoint_bracket_y *
      1.25
  )

endpoint_y_breaks <- log10_breaks_for_data(
  c(
    endpoint_plot_data$CFU,
    endpoint_plot_brackets$endpoint_label_y
  )
)

endpoint_diagnostic_plot <- ggplot()

if (
  isTRUE(
    SHOW_INDIVIDUAL_POINTS
  ) &&
  isTRUE(
    SHOW_INDIVIDUAL_PAIR_LINES
  )
) {
  endpoint_diagnostic_plot <- endpoint_diagnostic_plot +
    geom_line(
      data =
        endpoint_plot_data,
      aes(
        x =
          x_group,
        y =
          CFU,
        group =
          interaction(
            strain,
            replicate
          )
      ),
      color = "grey72",
      linewidth = 0.45,
      alpha = 0.75
    )
}

if (isTRUE(SHOW_INDIVIDUAL_POINTS)) {
  endpoint_diagnostic_plot <- endpoint_diagnostic_plot +
    geom_point(
      data =
        endpoint_plot_data,
      aes(
        x =
          x_group,
        y =
          CFU,
        color =
          strain,
        shape =
          treatment_label
      ),
      size = 2.0,
      alpha = 0.85
    )
}

endpoint_diagnostic_plot <- endpoint_diagnostic_plot +
  geom_errorbar(
    data =
      endpoint_summary,
    aes(
      x =
        x_group,
      ymin =
        lower_CFU,
      ymax =
        upper_CFU,
      color =
        strain
    ),
    width = 0.12,
    linewidth = 0.65
  ) +
  geom_point(
    data =
      endpoint_summary,
    aes(
      x =
        x_group,
      y =
        mean_CFU,
      color =
        strain
    ),
    shape = 21,
    fill = "white",
    stroke = 0.9,
    size = 3.0
  ) +
  geom_segment(
    data =
      endpoint_plot_brackets,
    aes(
      x =
        x_start,
      xend =
        x_end,
      y =
        endpoint_bracket_y,
      yend =
        endpoint_bracket_y
    ),
    linewidth = 0.45
  ) +
  geom_segment(
    data =
      endpoint_plot_brackets,
    aes(
      x =
        x_start,
      xend =
        x_start,
      y =
        endpoint_tick_bottom,
      yend =
        endpoint_bracket_y
    ),
    linewidth = 0.45
  ) +
  geom_segment(
    data =
      endpoint_plot_brackets,
    aes(
      x =
        x_end,
      xend =
        x_end,
      y =
        endpoint_tick_bottom,
      yend =
        endpoint_bracket_y
    ),
    linewidth = 0.45
  ) +
  geom_text(
    data =
      endpoint_plot_brackets,
    aes(
      x =
        (
          x_start +
            x_end
        ) /
        2,
      y =
        endpoint_label_y,
      label =
        bracket_label
    ),
    size = 3.0,
    lineheight = 1.0
  ) +
  scale_y_log10(
    breaks =
      endpoint_y_breaks,
    labels =
      log10_axis_labels,
    expand = expansion(
      mult = c(
        0.07,
        0.18
      )
    )
  ) +
  scale_x_continuous(
    breaks = c(
      1,
      2,
      4,
      5
    ),
    labels = c(
      paste0(
        CONTROL_NAME,
        "\n",
        NO_DRUG_LABEL
      ),
      paste0(
        CONTROL_NAME,
        "\n",
        DRUG_TREATMENT_LABEL
      ),
      paste0(
        TEST_NAME,
        "\n",
        NO_DRUG_LABEL
      ),
      paste0(
        TEST_NAME,
        "\n",
        DRUG_TREATMENT_LABEL
      )
    )
  ) +
  scale_color_manual(
    values =
      strain_colors,
    drop = FALSE
  ) +
  scale_shape_manual(
    values =
      treatment_shapes,
    drop = FALSE
  ) +
  labs(
    title = paste0(
      ENDPOINT_TIME_H,
      "-h endpoint"
    ),
    
    subtitle = if (
      isTRUE(
        SHOW_BETWEEN_STRAIN_P_VALUE
      )
    ) {
      between_test_label
    } else {
      NULL
    },
    
    x =
      NULL,
    
    y =
      expression(
        "Viable cells [CFU mL"^{-1}*"]"
      ),
    
    color =
      "Strain",
    
    shape =
      "Treatment"
  ) +
  theme_bw(
    base_size =
      BASE_SIZE,
    base_family =
      BASE_FAMILY
  ) +
  theme(
    panel.grid.major = element_line(
      color = "grey88",
      linewidth = 0.3
    ),
    
    panel.grid.minor = element_blank(),
    
    plot.title = element_text(
      face = "bold",
      size = 10
    ),
    
    plot.subtitle = element_text(
      size = 8.5
    ),
    
    axis.text.x = element_text(
      size = 8.2
    ),
    
    legend.position =
      "none"
  )


# Combine the time course and the endpoint detail into a second figure.
dual_killing_plot <- final_killing_plot +
  endpoint_diagnostic_plot +
  patchwork::plot_layout(
    widths = c(
      1.55,
      1
    )
  )


# Display the plots when requested. The plots are also saved below.
print_plot_in_rstudio(
  final_killing_plot
)

print_plot_in_rstudio(
  endpoint_diagnostic_plot
)

print_plot_in_rstudio(
  dual_killing_plot
)


# 10. Save the numerical results. Keep output prefixes based on the original
# workbook names so existing figure-assembly scripts can still locate them.
output_prefix <- tools::file_path_sans_ext(
  basename(
    INPUT_FILE
  )
)

# Record the settings, including target label, endpoint and zero-count handling.
analysis_settings <- tibble(
  parameter = c(
    "Input file",
    "Control strain",
    "Test strain",
    "Antimicrobial",
    "Antimicrobial concentration",
    "Concentration unit",
    "Endpoint time [h]",
    "Zero-CFU handling",
    "CFU detection limit [CFU/mL]",
    "Zero replacement fraction of LOD",
    "Zero replacement value [CFU/mL]",
    "Number of zero-CFU observations",
    "Paired target-versus-control test",
    "Show individual points",
    "Duplicate time-zero for treated curve"
  ),
  
  value = c(
    INPUT_FILE,
    CONTROL_NAME,
    TEST_NAME,
    ANTIMICROBIAL_NAME,
    antimicrobial_concentration,
    ANTIMICROBIAL_CONCENTRATION_UNIT,
    as.character(
      ENDPOINT_TIME_H
    ),
    ZERO_CFU_HANDLING,
    as.character(
      CFU_DETECTION_LIMIT
    ),
    as.character(
      ZERO_REPLACEMENT_FRACTION_OF_LOD
    ),
    as.character(
      zero_replacement_value
    ),
    as.character(
      n_zero_cfu
    ),
    as.character(
      PAIRED_CONTROL_AND_TEST
    ),
    as.character(
      SHOW_INDIVIDUAL_POINTS
    ),
    as.character(
      DUPLICATE_TIME_ZERO_FOR_TREATED_CURVE
    )
  )
)

# Save one workbook containing the tests, fold summaries, replicate-level
# reductions, time-course summaries, cleaned source data and analysis settings.
writexl::write_xlsx(
  list(
    between_strain_test =
      as.data.frame(
        between_strain_test
      ),
    
    within_strain_tests =
      as.data.frame(
        within_strain_tests
      ),
    
    fold_summary =
      as.data.frame(
        fold_summary
      ),
    
    endpoint_reductions =
      as.data.frame(
        endpoint_reductions
      ),
    
    timecourse_summary =
      as.data.frame(
        timecourse_summary
      ),
    
    cleaned_data =
      as.data.frame(
        cleaned_data
      ),
    
    analysis_settings =
      as.data.frame(
        analysis_settings
      )
  ),
  path = file.path(
    OUTPUT_DIR,
    paste0(
      output_prefix,
      "_killing_analysis.xlsx"
    )
  ),
  format_headers = TRUE
)

# Also save the main test and replicate reductions as plain-text tables.
readr::write_tsv(
  between_strain_test,
  file.path(
    OUTPUT_DIR,
    paste0(
      output_prefix,
      "_between_strain_test.tsv"
    )
  )
)

readr::write_tsv(
  endpoint_reductions,
  file.path(
    OUTPUT_DIR,
    paste0(
      output_prefix,
      "_endpoint_reductions.tsv"
    )
  )
)


# 11. Figure export settings. Change these dimensions to adjust output size.
FINAL_PLOT_WIDTH_IN <- 3.2
FINAL_PLOT_HEIGHT_IN <- 4

DUAL_PLOT_WIDTH_IN <- 7
DUAL_PLOT_HEIGHT_IN <- 4

SAVE_PNG <- TRUE
SAVE_PDF <- TRUE
SAVE_TIFF <- FALSE

PNG_DPI <- 300
TIFF_DPI <- 300

# Save the main time course and the combined figure using the chosen formats.
save_plot_outputs(
  plot_object =
    final_killing_plot,
  basename = paste0(
    output_prefix,
    "_killing_final"
  ),
  width =
    FINAL_PLOT_WIDTH_IN,
  height =
    FINAL_PLOT_HEIGHT_IN,
  save_png =
    SAVE_PNG,
  save_pdf =
    SAVE_PDF,
  save_tiff =
    SAVE_TIFF,
  png_dpi =
    PNG_DPI,
  tiff_dpi =
    TIFF_DPI
)

save_plot_outputs(
  plot_object =
    dual_killing_plot,
  basename = paste0(
    output_prefix,
    "_killing_dual"
  ),
  width =
    DUAL_PLOT_WIDTH_IN,
  height =
    DUAL_PLOT_HEIGHT_IN,
  save_png =
    SAVE_PNG,
  save_pdf =
    SAVE_PDF,
  save_tiff =
    SAVE_TIFF,
  png_dpi =
    PNG_DPI,
  tiff_dpi =
    TIFF_DPI
)


# Keep the finished plots and key result tables available in the R session.
assign(
  "killing_final_plot",
  final_killing_plot,
  envir = .GlobalEnv
)

assign(
  "killing_endpoint_plot",
  endpoint_diagnostic_plot,
  envir = .GlobalEnv
)

assign(
  "killing_dual_plot",
  dual_killing_plot,
  envir = .GlobalEnv
)

assign(
  "killing_endpoint_reductions",
  endpoint_reductions,
  envir = .GlobalEnv
)

assign(
  "killing_between_strain_test",
  between_strain_test,
  envir = .GlobalEnv
)

# Print a short completion message with the input-specific settings and output folder.
message(
  "Completed killing-assay analysis.",
  "\nControl: ",
  CONTROL_NAME,
  ".",
  "\nTest strain: ",
  TEST_NAME,
  ".",
  "\nAntimicrobial: ",
  ANTIMICROBIAL_NAME,
  " ",
  antimicrobial_concentration,
  " ",
  ANTIMICROBIAL_CONCENTRATION_UNIT,
  ".",
  "\nEndpoint: ",
  ENDPOINT_TIME_H,
  " h.",
  "\nIndividual points shown: ",
  SHOW_INDIVIDUAL_POINTS,
  ".",
  "\nZero-CFU observations: ",
  n_zero_cfu,
  ".",
  if (n_zero_cfu > 0L) {
    paste0(
      "\nZero replacement used: ",
      zero_replacement_value,
      " CFU/mL."
    )
  } else {
    ""
  },
  "\nOutput directory:\n  ",
  OUTPUT_DIR
)
