# Guide-efficacy model from this study

Run the included ten-guide PAO1 *efp* (PA2851) example from the project root:

```r
source("guide_efficacy_prediction/this_study/AMG_score_any_sgRNA_set_with_model.R")
```

Required R packages are `data.table` and `glmnet`. The script reads
`guide_efficacy_prediction/this_study/CRISPRi_guide_efficacy_model.rds` and
`guide_efficacy_prediction/this_study/input/example_PA2851_efp_guides.csv`.
Results are saved in `guide_efficacy_prediction/this_study/out_user_scores/` and returned
as `kawalek_scored_guides` in RStudio. This example runs independently of the
screen and figure analyses.

## Input format

To score another table, change `INPUT_FILE` and `OUTPUT_FILE` at the top of
the script. CSV, TSV and gzip-compressed versions are supported.
Setting `OUTPUT_FILE <- NULL` writes the result beside the input file.

Use these Hawkins-style column names:

| Column | Contents |
| --- | --- |
| `locus_tag` | Locus identifier used to group candidates for ranking |
| `sgrna` | 20-nt spacer in the design-output orientation, without PAM or cloning adapters |
| `pam` | NGG PAM in the same orientation as the spacer |
| `downstream` | At least 15 nt immediately after the PAM in that orientation; the first 15 are used |
| `offset` | Distance from the annotated gene start, using the Hawkins-style design convention |
| `off_target_n` | Predicted off-target count; the column is required but can be NA for the supplied model, which does not use this predictor |

Additional columns, such as `gene_name`, are retained. Rename columns explicitly
when using another designer and ensure that sequence orientation and offsets
follow the expected convention. Supply full-match candidates; the model is
not calibrated for mismatch-titration variants or non-targeting controls.

## Interpreting scores

Higher `kawalek_efficacy_score` indicates stronger predicted relative guide
performance. Scores are unconstrained regression predictions, not
probabilities or measured transcript-repression fractions.

Centered scores and ranks are calculated among the candidates supplied for
each locus. Changing that candidate set can change its centered scores and
ranks; raw model scores remain unchanged.

The supplied model uses sequence and target offset and was trained with
6,868 guides from 460 genes. Applying it to the included library-guide example
demonstrates scoring, not independent model validation. See the
[Figure S8 guide](../../Figure_scripts/FigS8/README_FigS8.md) for model development
and evaluation.

The fitted model and the supplied `Lib1_all_guides_with_model_scores.tsv.gz`
and `model_training_guides_with_scores.tsv.gz` tables are retained in this
folder. Figure scripts using these files should point to
`guide_efficacy_prediction/this_study/`.
