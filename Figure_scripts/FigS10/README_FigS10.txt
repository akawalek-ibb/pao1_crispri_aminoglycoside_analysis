# Figure S10: comparison of MH and FBS vulnerability

Run each panel script from the project root:

```r
source("Figure_scripts/FigS10/AMG_FigS10A_MH_FBS_severity.R")
source("Figure_scripts/FigS10/AMG_FigS10B_MH_FBS_sensitivity.R")
source("Figure_scripts/FigS10/AMG_FigS10C_MH_FBS_EK50.R")
```

## Inputs and preparation

All panels read:

- `preprocessing_scripts/out_13_AMG_MH_FBS_gene_effects/MH_FBS_gene_effects.tsv`
- `annotation_input/PAO1_curated.tsv`

The shared out_13 table is not included in this bundle. Generate it by running
`preprocessing_scripts/11_AMG_gene_vulnerability.R`, then
`preprocessing_scripts/13_AMG_MH_FBS_gene_effects.R`. These use the supplied fitted
screen results; metacounter does not need to be rerun for those existing inputs.
If the prepared table is already available, no upstream calculation is required.

Packages: ggplot2, ggrepel and svglite. Each panel runs independently, displays
its plot in RStudio and prints source-data paths. Guide-response curves are not
refitted here. No per-generation normalization is applied.

## Selection and interpretation

- **A, severity:** compare B3 MH T1/T0 with FBS T1/T0. Labels and point emphasis
  use direct FBS T1 versus MH T1 tests: raw Stouffer p < 0.05 and |median log2FC|
  >= 1. The solid diagonal shows equality; dashed lines show +/-1 differences.
- **B, sensitivity:** require significant negative fitness in either medium,
  with an evaluable profile in both. Only `no_sustained_crossing` is displayed
  at 1; its underlying estimate remains NA. Insufficient-guide/bin profiles are
  excluded. Correlations and the paired Wilcoxon test use numeric pairs only.
- **C, EK50:** require significant negative fitness and valid numeric EK50 in
  both media. No missing value is substituted with 1. Correlations and the
  paired Wilcoxon test use the plotted numeric pairs.

Negative fitness uses raw Stouffer p < 0.05 and median log2FC <= -1.
Essentiality categories come from the supplied curated annotation. Display
palettes and dimensions differ among panels according to their own settings.
Gene counts and statistics are derived from the selected input rather than
forced to a figure-caption count.

## Outputs

The three scripts write to their respective `out_AMG_FigS10A_severity/`,
`out_AMG_FigS10B_sensitivity/` and `out_AMG_FigS10C_EK50/` folders beside the
scripts. Each exports PNG/PDF/SVG plus all-gene, plotted-gene, statistics and
settings TSVs. Panel B also exports the numeric pairs used for statistics.
RStudio objects are `FigS10A_plot`, `FigS10B_plot` and `FigS10C_plot`.
