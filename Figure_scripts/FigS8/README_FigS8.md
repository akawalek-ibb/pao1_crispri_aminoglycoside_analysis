# Figure S8: guide features, model development and evaluation

Run the scripts from the project root. The four figure scripts are in
`Figure_scripts/FigS8/`.

| Script | Purpose |
| --- | --- |
| `AMG_FigS8A_B_guide_position_offtargets.R` | Target-position and predicted off-target analyses |
| `AMG_FigS8C_F_train_evaluate_model.R` | Lasso model development and repeated evaluation holding out whole genes |
| `AMG_FigS8G_selected_gene_examples.R` | Guide predictions for efp, lptF and parE, refitting with the target gene excluded |
| `AMG_FigS8H_three_model_performance.R` | Comparison of guide rankings from the study model, CRISPRi-Pasteur and MERF/CIAO |

Candidate-guide scoring with the saved model is a separate workflow described
in the [guide-scoring README](../../guide_efficacy_prediction/this_study/README.md).

## Panels A–F

Supply the eight processed screen files in `edgeR_data/`:

```text
b1_melted_results.tsv.gz             b1_median_melted_results.tsv.gz
b2_melted_results.tsv.gz             b2_median_melted_results.tsv.gz
b3_melted_results.tsv.gz             b3_median_melted_results.tsv.gz
b4_melted_results.tsv.gz             b4_median_melted_results.tsv.gz
```

These fitted exports are included and are produced by
`preprocessing_scripts/02_AMG_metacounter.R` when rebuilding. Reuse them for
figure work. S8A/B requires the separate count/context layout from script 25:

```r
source("preprocessing_scripts/25_AMG_guide_efficacy_inputs.R")
source("Figure_scripts/FigS8/AMG_FigS8A_B_guide_position_offtargets.R")
source("Figure_scripts/FigS8/AMG_FigS8C_F_train_evaluate_model.R")
```

Task 25 reads `Barcode_counting/count_tables/b*_overall_stats.tsv`, batch
designs and `annotation_input/PAO1_key.tsv`. It writes annotated count matrices
and design copies to `preprocessing_scripts/out_25_AMG_guide_efficacy_inputs/`.
The guide key supplies spacer, PAM, downstream context, offset and off-target
annotations. Filtering and model fitting take place in the figure scripts.

S8C–F prefers the supplied frozen `edgeR_data/raw_counts/b*_pao_join_raw_counts.tsv.gz`
files and accepts the original `Barcode_counting/count_tables/` designs. S8A/B
does not search that frozen-count folder, so run script 25 for A/B. If changing
the counts, note that regenerating 25 does not overwrite C–F's preferred frozen
files. Select matching count/context files and fitted results explicitly.

The optional genome-wide scoring block in C–F requires the externally supplied
`Figure_scripts/FigS8/input_AMG_FigS8/PAO1_107_relaxed_unfiletered.csv`.
When that file is absent, the block is skipped with a warning.

A–B outputs go to `out_AMG_FigS8A_B_guide_properties/`; C–F outputs go to
`out_AMG_FigS8C_F_model/`, both under the Figure S8 folder. C–F exports a fitted
model and its training table there.

## Panels G–H

Set `MODEL_DIR` to a directory containing this matched pair:

```text
CRISPRi_guide_efficacy_model.rds
model_training_guides_with_scores.tsv.gz
```

The default directory is `guide_efficacy_prediction/this_study/`. The study training
table contains 6,868 eligible guides from 460 genes, including their sequence
context and centered response `y`. Preserve row order when reproducing seeded
splits. `Lib1_all_guides_with_model_scores.tsv.gz` contains deployment scores
and cannot substitute for the training-response table.

The bundle includes the two external-model prediction tables at:

```text
Figure_scripts/FigS8/input_AMG_FigS8/PAO1_sgRNA_with_Pasteur_scores.tsv.gz
Figure_scripts/FigS8/input_AMG_FigS8/PAO1_all_sgRNA_predictions_CIAO.csv
```

Then run:

```r
source("Figure_scripts/FigS8/AMG_FigS8G_selected_gene_examples.R")
source("Figure_scripts/FigS8/AMG_FigS8H_three_model_performance.R")
```

To evaluate a newly trained model, point `MODEL_DIR` in both scripts to the
C–F output folder, which contains the model and training table. G–H outputs
are saved in `out_AMG_FigS8G_selected_genes/` and
`out_AMG_FigS8H_model_comparison/` under the Figure S8 folder. Console messages
identify the underlying figure-data tables.

## Interpretation and adaptation

- The modeled response is within-gene-centered guide-associated fitness,
  aggregated within and across batches. It is not a measured fraction of
  transcript repression.
- Lambda tuning holds out whole genes. Feature-set selection uses repeated
  outer evaluations, so performance of the selected feature set is not an
  independent final test. D–F use the split closest to the median outer
  Spearman correlation.
- Panel B uses the explicit `non_essential` design category when available;
  this is not equivalent to every gene outside the negative-fitness set.
- Panel G defaults to efp, lptF and parE, in that order. Panel H derives its
  gene set from common guide predictions and valid held-out estimates; the
  evaluated count is reported in its output rather than forced to 460.
- CIAO is the software label used for MERF predictions. Its Activity score is
  negated so that higher plotting scores indicate stronger predicted activity.

Edit input paths, filters, seeds, model features and plotting settings in the
relevant script. Use matching R/package versions when reproducing fits.
Required packages are listed or loaded in each script.
