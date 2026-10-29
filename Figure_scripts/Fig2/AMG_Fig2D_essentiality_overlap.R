# Draw the four-set essentiality overlap for Figure 2D.
# Working directory: manuscript project root (the directory containing AMG_cleaned_Git.Rproj).
# Run in R/RStudio: source("Figure_scripts/Fig2/AMG_Fig2D_essentiality_overlap.R")
# Inputs:
# - preprocessing_scripts/out_12_AMG_MH_fitness/Essentiality_membership.tsv: unique PAO1_ID and
#   four logical membership columns.
# Upstream / execution order:
# - preprocessing_scripts/11_AMG_gene_vulnerability.R → preprocessing_scripts/12_AMG_MH_fitness.R.
# Outputs:
# - Figure_scripts/Fig2/out_AMG_Fig2D_essentiality_overlap/ — PNG/PDF/SVG, membership/overlap
#   counts and missed-core TSVs.
# Interpretation / operation:
# - Reference sets are the supplied Poulsen core, mapped Zhang and Kaczmarczyk shared299
#   memberships.
#
# R packages used/checked: data.table, ggVennDiagram, ggplot2, svglite.

INPUT_FILE <- "preprocessing_scripts/out_12_AMG_MH_fitness/Essentiality_membership.tsv"

OUTPUT_DIR <- "Figure_scripts/Fig2/out_AMG_Fig2D_essentiality_overlap"

SET_LABELS <- c(Poulsen_core = "Poulsen Core", Zhang_2025 = "Zhang 2025", This_study = "Our study", 
  Kaczmarczyk_2025 = "Kaczmarczyk 2025 (all media)")

CHECK_MANUSCRIPT_SET_SIZES <- TRUE

EXPECTED_SET_SIZES <- c(Poulsen_core = 320L, Zhang_2025 = 580L, This_study = 743L, Kaczmarczyk_2025 = 299L)

# Outline/text colours and overlap tones from the supplied Fig. 2D reference.
SET_COLOURS <- c(Poulsen_core = "#8A0000", Zhang_2025 = "#504416",
  This_study = "#000000", Kaczmarczyk_2025 = "#009E73")
REGION_FILLS <- c(
  "1" = "#DEB3B3", "2" = "#CDC7B9", "3" = "#D2D2D2", "4" = "#B3E3D5",
  "1/2" = "#BA8E84", "1/3" = "#BB9C9C", "1/4" = "#9BAFA0",
  "2/3" = "#AFACA0", "2/4" = "#91BDA4", "3/4" = "#94C2B7",
  "1/2/3" = "#9F8079", "1/2/4" = "#82917D", "1/3/4" = "#829D91",
  "2/3/4" = "#79A693", "1/2/3/4" = "#708977"
)
SET_LABEL_SIZE_PT <- 14
COUNT_SIZE_PT <- 16
OUTLINE_WIDTH <- 0.8

FONT_FAMILY <- "Arial"

WIDTH_IN <- 4

HEIGHT_IN <- 2.9

DPI <- 300

PRINT_IN_RSTUDIO <- TRUE

OVERWRITE <- TRUE

# Local functions for this analysis.
amg_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) 
    stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

amg_required <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) 
    stop(label, " lacks: ", paste(missing, collapse = ", "), call. = FALSE)
}

amg_read <- function(path) {
  if (!file.exists(path)) 
    stop("Missing input: ", path, call. = FALSE)
  con <- if (endsWith(path, ".gz")) 
    gzfile(path, "rt")
  else file(path, "rt")
  on.exit(close(con))
  data.table::as.data.table(utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE, 
    quote = "\"", comment.char = "", na.strings = c("", "NA")))
}

amg_write <- function(tables, directory, overwrite = TRUE, gzip = character()) {
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

amg_settings <- function(settings, inputs) {
  data.table::data.table(setting = c(names(settings), paste0("input:", names(inputs)), paste0("md5:", 
    names(inputs)), "R_version"), value = c(as.character(settings), unname(inputs), unname(tools::md5sum(inputs)), 
    R.version.string))
}

fig2_logical <- function(x, label) {
  x <- tolower(as.character(x))
  if (anyNA(x) || any(!x %in% c("true", "false", "1", "0"))) 
    stop("Invalid logical values: ", label)
  x %in% c("true", "1")
}

fig2_export <- function(plot, tables, output, stem, width, height, dpi, overwrite = TRUE) {
  names_out <- c(paste0(stem, c(".png", ".pdf", ".svg")), paste0(names(tables), ".tsv"))
  if (!overwrite && any(file.exists(file.path(output, names_out)))) 
    stop("An output exists; set OVERWRITE <- TRUE to replace it.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(".fig2_", tmpdir = output)
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE))
  ggplot2::ggsave(file.path(stage, names_out[1L]), plot, width = width, height = height, dpi = dpi, 
    bg = "white")
  ggplot2::ggsave(file.path(stage, names_out[2L]), plot, device = grDevices::cairo_pdf, width = width, 
    height = height, bg = "white")
  ggplot2::ggsave(file.path(stage, names_out[3L]), plot, device = svglite::svglite, width = width, 
    height = height, bg = "white")
  amg_write(tables, stage, TRUE)
  for (name in names_out) if (!file.copy(file.path(stage, name), file.path(output, name), overwrite = overwrite)) 
    stop("Could not save ", name)
  invisible(plot)
}

# Use the stored TRUE/FALSE memberships as four distinct sets; no new
# essentiality call or guide-coverage filter is applied. Each four-digit binary
# pattern identifies an exclusive overlap region, counted once per gene.
prepare_overlap <- function(input) {
  members <- data.table::copy(input$membership)
  amg_required(members, names(SET_LABELS), "Essentiality_membership.tsv")
  sets <- lapply(names(SET_LABELS), function(field) members$PAO1_ID[members[[field]]])
  names(sets) <- unname(SET_LABELS)
  sizes <- stats::setNames(lengths(sets), names(SET_LABELS))
  if (CHECK_MANUSCRIPT_SET_SIZES && (!identical(names(sizes), names(EXPECTED_SET_SIZES)) || any(sizes != 
    EXPECTED_SET_SIZES))) 
    stop("Set sizes differ from the manuscript: ", paste(names(sizes), sizes, collapse = "; "), 
      ". Check S5 inputs before changing the expected sizes.")
  members[, `:=`(region_pattern, apply(as.matrix(.SD) * 1L, 1L, paste0, collapse = "")), .SDcols = names(SET_LABELS)]
  regions <- members[, .(n_genes = .N), by = region_pattern]
  data.table::setorder(regions, region_pattern)
  # Identify the Poulsen core genes absent from all three CRISPRi-derived sets.
  # This is a set-membership comparison, not a test of why a gene was missed.
  missed <- members[Poulsen_core & !Zhang_2025 & !This_study & !Kaczmarczyk_2025]
  list(sets = sets, members = members, regions = regions, sizes = data.table::data.table(set = names(sizes), 
    label = unname(SET_LABELS), n_genes = unname(sizes)), missed = missed)
}

# Read, prepare, plot and export in order.
amg_require(c("data.table", "ggplot2", "ggVennDiagram", "svglite"))

# 2. Read the fixed gene-set memberships
membership <- amg_read(INPUT_FILE)

for (field in names(SET_LABELS)) membership[[field]] <- fig2_logical(membership[[field]], field)

input <- list(membership = membership, paths = INPUT_FILE)

# 3. Count intersections and construct the Venn diagram
result <- prepare_overlap(input)

# Use the four-set geometry and counts from ggVennDiagram, with the reference
# fills, outlines and text instead of a count-dependent colour gradient.
venn_data <- ggVennDiagram::process_data(ggVennDiagram::Venn(result$sets))
region_shapes <- as.data.frame(ggVennDiagram::venn_regionedge(venn_data))
set_edges <- as.data.frame(ggVennDiagram::venn_setedge(venn_data))
count_labels <- as.data.frame(ggVennDiagram::venn_regionlabel(venn_data))
region_shapes$fill <- unname(REGION_FILLS[region_shapes$id])
set_edges$colour <- unname(SET_COLOURS[as.integer(set_edges$id)])
if (anyNA(region_shapes$fill) || anyNA(set_edges$colour)) stop("Missing reference Venn colour.")

sizes <- result$sizes$n_genes
set_labels <- data.frame(
  X = c(0, 0.33, 0.83, 1.025),
  Y = c(0.57, 0.963, 0.963, 0.57),
  label = c(
    paste0("Poulsen et al., 2019\ncore essential (", sizes[1], ")"),
    paste0("Zhang et al., 2025\nessential (", sizes[2], ")"),
    paste0("This study\n(", sizes[3], ")"),
    paste0("Kaczmarczyk et al., 2025\nshared fitness (", sizes[4], ")")
  ),
  colour = unname(SET_COLOURS),
  angle = c(90, 0, 0, 270)
)
pt_per_mm <- 72.27 / 25.4
p <- ggplot2::ggplot() +
  ggplot2::geom_polygon(data = region_shapes,
    ggplot2::aes(X, Y, group = id, fill = fill), colour = NA) +
  ggplot2::geom_path(data = set_edges,
    ggplot2::aes(X, Y, group = id, colour = colour),
    linewidth = OUTLINE_WIDTH, lineend = "round") +
  ggplot2::geom_text(data = count_labels,
    ggplot2::aes(X, Y, label = count), family = FONT_FAMILY,
    fontface = "plain", colour = "black", size = COUNT_SIZE_PT / pt_per_mm) +
  ggplot2::geom_text(data = set_labels,
    ggplot2::aes(X, Y, label = label, colour = colour, angle = angle),
    family = FONT_FAMILY, fontface = "bold", lineheight = 0.9,
    size = SET_LABEL_SIZE_PT / pt_per_mm, hjust = 0.5, vjust = 0.5) +
  ggplot2::scale_fill_identity() +
  ggplot2::scale_colour_identity() +
  ggplot2::coord_fixed(xlim = c(-0.125, 1.125), ylim = c(0.14, 1.04625),
    expand = FALSE, clip = "off") +
  ggplot2::theme_void(base_family = FONT_FAMILY) +
  ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(0, 0, 0, 0))

settings <- amg_settings(c(set_order = paste(names(SET_LABELS), collapse = "; "), Kaczmarczyk_reference = "Frozen shared299 intersection across LB, glucose, succinate", 
  calls = "Stored S5 membership; no additional filter",
  font_family = FONT_FAMILY, set_label_size_pt = SET_LABEL_SIZE_PT,
  count_size_pt = COUNT_SIZE_PT, outline_width = OUTLINE_WIDTH,
  set_colours = paste(names(SET_COLOURS), SET_COLOURS, collapse = "; "),
  region_fills = paste(names(REGION_FILLS), REGION_FILLS, collapse = "; ")), input$paths)

# Save both the memberships and region counts so the diagram is reproducible
# without inferring gene identities from the displayed overlap numbers.
tables <- list(Fig2D_membership = result$members, Fig2D_region_counts = result$regions, Fig2D_set_sizes = result$sizes, 
  Fig2D_core_genes_missed_by_three_screens = result$missed, Analysis_settings = settings)

fig2_export(p, tables, OUTPUT_DIR, "Fig2D_essentiality_overlap", WIDTH_IN, HEIGHT_IN, DPI, OVERWRITE)

if (isTRUE(PRINT_IN_RSTUDIO)) print(p)

message("Saved Fig. 2D; ", nrow(result$missed), " Poulsen core genes missed by all three CRISPRi sets.")

data_paths <- normalizePath(file.path(OUTPUT_DIR, paste0(names(tables), ".tsv")), winslash = "/", 
  mustWork = TRUE)

message("Underlying figure data (Fig. 2D):\n", paste0("  ", data_paths, collapse = "\n"))

Fig2D_result <- c(list(plot = p), result)

Fig2D_plot <- Fig2D_result$plot

