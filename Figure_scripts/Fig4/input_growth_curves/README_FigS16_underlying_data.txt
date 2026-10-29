Figure 4 and Supplementary Figure S16: growth-curve source data

The same gene-specific growth-analysis scripts generate the individual
plots and underlying data used for Figure 4 and Supplementary Figure S16.
Run them from the RStudio project working directory:

  source("Figure_scripts/Fig4/input_growth_curves/accC/accC_growth_analysis_pvalue_plus_effectsize.R")
  source("Figure_scripts/Fig4/input_growth_curves/atpI/atpI_growth_analysis_pvalue_plus_effectsize.R")
  source("Figure_scripts/Fig4/input_growth_curves/ftsZ/ftsZ_growth_analysis_pvalue_plus_effectsize.R")
  source("Figure_scripts/Fig4/input_growth_curves/secA/secA_growth_analysis_pvalue_plus_effectsize.R")

All four scripts use a one-sample t-test comparing observed combined-treatment
well metrics with the fixed multiplicative expected mean. The default
alternative is "less". Raw p values are used for labels by default.
All four scripts require BOTH selected p < 0.05 and normalized interaction
effect <= -0.10 for sensitizing labels. The normalized effect is
(mean combined - multiplicative expected mean) / mean untreated.
Both requirements must pass; otherwise the displayed label is "ns".


