# Figure S9: fitness-response descriptors and transcriptional-unit content

Run from the project root using prepared task-12 tables:

```r
source("Figure_scripts/FigS9/AMG_FigS9A_severity_EK50.R")
source("Figure_scripts/FigS9/AMG_FigS9B_C_D_TU_content.R")
```

## Inputs

Both scripts read
`preprocessing_scripts/out_12_AMG_MH_fitness/Fig2G_input.tsv` and
`annotation_input/PAO1_curated.tsv`. Panels B–D also use
`preprocessing_scripts/out_12_AMG_MH_fitness/MH_fitness.tsv`.
If rebuilding these tables, run task 11 followed by task 12.

## Panels and selection

- **A:** severity versus EK50 for genes with finite values of both descriptors.
  Point size reflects TU gene content; the curve is a GAM fit to the plotted
  association. Spearman correlation is also reported.
- **B:** severity versus TU gene content.
- **C:** sensitivity versus TU gene content.
- **D:** EK50 versus TU gene content.

Consensus values are gene-wise medians across B1–B4 untreated T1/T0 comparisons;
B1 T2 is excluded. Severity and sensitivity require at least one finite
estimate, while EK50 requires at least two valid estimates.

Panels B–D include genes with raw Stouffer p < 0.05 and medLFC ≤ −1 in at
least one of the four T1 controls. TU counts are displayed as 0, 1, 2, 3, 4,
5 and 6+; Spearman tests use uncapped counts. Gene-level tests include each
selected gene, whereas TU-level tests use one median per TU to avoid counting
the same TU repeatedly. Missing TUs are excluded from B–D and assigned zero
for point sizing in A.

## TU count definition

Both scripts default to `TU_COUNT_DEFINITION <- "design_call"`, which counts
genes in the `main_essential` and `other_essential` library-design groups.
These include fitness-associated genes, so the resulting count must be
described as **essential/fitness-associated genes**, not exclusively genes
classified as essential in published Tn-seq studies.

For the latter definition, set `TU_COUNT_DEFINITION <- "Nr_E"` in both scripts.
This counts genes with `Nr_E >= 1` and changes point sizes, group membership
and potentially association statistics. Use the same definition in related
figures and their legends.

## Outputs and settings

Outputs are saved beside the scripts in
`out_AMG_FigS9A_severity_EK50/` and `out_AMG_FigS9B_C_D_TU_content/`.
They include SVG/PDF/PNG figures, plotted values, statistical summaries and
run settings. Panel A also exports its fitted smooth curve.

Edit input paths, selection thresholds, TU definition, labels, dimensions and
font near the top of each script. `PRINT_IN_RSTUDIO` controls plot display;
`OVERWRITE` controls replacement of existing outputs. These scripts reuse the
prepared gene-response estimates; they do not refit guide-response curves.
