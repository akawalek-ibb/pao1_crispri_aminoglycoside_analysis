# Figure S14: aminoglycoside-associated repression thresholds

Run from the project root:

```r
source("Figure_scripts/FigS14/AMG_FigS14_aminoglycoside_threshold_shifts.R")
```

Requires ggplot2, ggrepel and svglite. Individual panels and the combined figure
are displayed in RStudio and remain as `FigS14_panels` and `FigS14_plot`.

## Inputs and preparation

- `preprocessing_scripts/out_11_AMG_gene_vulnerability/gene_vulnerability_long.tsv.gz`:
  sensitivity estimates and evaluation statuses from independently calculated
  untreated/T0 and treated/T0 profiles.
- `preprocessing_scripts/out_03_AMG_relevant_results/Stouffer_calls.tsv`:
  direct treated/untreated interaction calls and T1 recurrence flags.
- `Figure_scripts/FigS14/input_AMG_FigS14/absolute_gene_effects.tsv.gz`:
  included compact gene-level effects and p-values against T0, with PAO1_ID,
  batch, condition, medLFC and p_stouffer columns.

The first two prepared outputs are not included. Generate them using
`preprocessing_scripts/11_AMG_gene_vulnerability.R` and
`preprocessing_scripts/03_AMG_relevant_results.R`, respectively, reusing the
supplied fitted screen exports. The figure script does not refit sensitivities
or recalculate direct interaction calls.

The compact input contains the ten absolute contrasts needed for the six
paired panels. Absolute-fitness p-values cannot be replaced with direct
interaction p-values. To rebuild this input from the supplied b1/b2/b4 median
metacounter exports, set `REBUILD_GENE_EFFECTS <- TRUE`; this overwrites the
compact file when `OVERWRITE` is TRUE. Keep it FALSE for ordinary plotting.
The figure checks agreement between its absolute effects and the vulnerability
master's severity values before combining them.

## Selection and display

Include genes with negative fitness in at least one matched absolute profile
(raw Stouffer p < 0.05, median log2FC <= -1) and evaluable sensitivities in both.
Only `no_sustained_crossing` is displayed at 1; insufficient-information profiles
are excluded. The underlying missing estimates remain NA.

Colours represent the current direct interaction: sensitizing red, protective
blue, otherwise grey. Shapes distinguish same-direction recurrent interactions,
condition-specific interactions and other non-recurrent interactions. Recurrence
requires at least two of the five T1 calls; T2 does not contribute. Here,
condition-specific means significant in only one of the six displayed comparisons.

The equality diagonal and +/-0.15 band are descriptive, not significance tests.
Label current interaction hits in TOB panels; for AMK, label the 20 largest
absolute displayed threshold shifts, including eligible boundary points.

## Outputs

`Figure_scripts/FigS14/out_AMG_FigS14_repression_thresholds/` contains PNG/PDF/SVG
figures and the all-gene status, plotted-gene, label, statistics, comparison
manifest and settings tables. The all-gene table is gzip-compressed. True
numeric differences are retained separately from display-boundary differences.

Set `SAVE_INDIVIDUAL_PANELS <- FALSE` to export/display only the combined figure.
The script prints the underlying-data file locations on completion.
