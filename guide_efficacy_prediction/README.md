# CRISPRi guide-efficacy prediction

Three tools are provided in separate folders. Run commands from the manuscript
project root; each folder has its own inputs, instructions and dependencies.

| Folder | Purpose | Example |
| --- | --- | --- |
| [this_study/](this_study/README.md) | Apply the PAO1 model developed in this study to full-match guides | Ten PAO1 *efp* guides |
| [pasteur/](pasteur/README.md) | Score an existing guide catalogue with the published Pasteur EcoWG1 model | 46 LESB58 *efp* guides |
| [ciao/](ciao/README.md) | Design and score guides with the published CIAO/MERF model, using a wrapper for parallel processing of multiple genes | LESB58 *efp*, *secA* and *accC* |

For this study's model, source:

```r
source("guide_efficacy_prediction/this_study/AMG_score_any_sgRNA_set_with_model.R")
```

For CIAO, run the numbered scripts in order: `01` verifies or downloads the
model, `02` runs predictions and `03` combines the results. The optional
`00` script makes a gene inventory from a GFF. The core predictor is called
automatically by `02`.

The fitted PAO1 model, its training-guide table and library-guide scores are
in `this_study/`. Scripts outside this folder that use those files must use
that path. The supplied CIAO model files and completed example outputs are
retained; the new three-gene example writes to a separate output folder.

Scores have different scales: higher scores indicate stronger predicted
activity for this study's and Pasteur models; lower CIAO `Activity score`
indicates stronger predicted depletion. None is a measured transcript-repression
fraction.

## Published models

- **Pasteur:** Calvo-Villamañán et al. (2020), *Nucleic Acids Research* 48:e64.
  [On-target activity predictions enable improved CRISPR–dCas9 screens in bacteria](https://doi.org/10.1093/nar/gkaa294).
- **CIAO/MERF:** Yu et al. (2024), *Genome Biology* 25:13.
  [Improved prediction of bacterial CRISPRi guide efficiency from depletion screens through mixed-effect machine learning and data integration](https://doi.org/10.1186/s13059-023-03153-y).
