# iModulon analysis for Figure 6C and Supplementary Table S9

`imodulon_activity_pipeline.py` estimates changes in the activity of PRECISE411
reference gene-expression modules (iModulons) from sample-level RNA-seq data.
Its script, configuration and requirements are in
`RNAseq_data/iModulon_analysis/preprocesing/`.

## Inputs

Paths below are relative to the project root.

| File | Contents |
| --- | --- |
| `RNAseq_data/iModulon_analysis/input/CRISPRi_silencing_RNAseq_log2_TPM.tsv` | Study expression matrix |
| `RNAseq_data/iModulon_analysis/input/CRISPRi_silencing_RNAseq_iModulon_samples.tsv` | Matching sample identifiers and condition groups |
| `annotation_input/precise411-M.csv` | Reference gene weights for the 116 modules |
| `annotation_input/precise411-im_table.csv` | Module names and annotations |
| `annotation_input/precise411-M_thresholds.csv` | Gene-membership thresholds used for module diagnostics |

The expression matrix must have `gene_id` as its first column and sample IDs
as the remaining headers. Values must be log2(TPM + 1), before control centering
and without differential-expression filtering. The sample sheet must contain
`sample_id` and `group`, with identifiers matching the expression columns.
The study groups are `control1`, `control1_TOB`, `ftsZg1m4`, `atpIg1m2`,
`accCg1m3` and `secAg1m3`; these match `config.json` and the upstream R export.

The expression matrix and sample sheet are exported by
`RNAseq_analysis/AMG_RNAseq_featureCounts_edgeR.R`. Copy both files from
`RNAseq_analysis/out/` into the input folder above, or update their paths in
`config.json`.

## Reference weights

Use the PRECISE411 weights with the matching annotation and threshold files.
The reference is available from [iModulonDB](https://imodulondb.org/p_aeruginosa/precise411)
and its [dataset download page](https://legacy.imodulondb.org/dataset.html?organism=p_aeruginosa&dataset=precise411)
([Rajput et al., 2022](https://doi.org/10.1093/nar/gkac743)).

If `annotation_input/precise411-M.csv` is missing, activate the Python environment
described below and run this from the project root. It downloads the weights
and replaces numbered component headers with names from the annotation table.

```bash
python - <<'PY'
import pandas as pd

url = "https://legacy.imodulondb.org/organisms/p_aeruginosa/precise411/data_files/M.csv"
weights = pd.read_csv(url, index_col=0)
atlas = pd.read_csv("annotation_input/precise411-im_table.csv")
names = dict(zip(atlas["k"].astype(str), atlas["imodulon_name"]))
if set(weights.columns) != set(names):
    raise ValueError("The weights and module annotation do not match.")
weights = weights.rename(columns=names).sort_index()
weights.index.name = "gene_id"
weights.to_csv("annotation_input/precise411-M.csv")
PY
```

Keep the reference weights in their native sign and scale.

## Run

From a terminal in the project root, create and activate the Python 3.9
environment, then install the pinned requirements:

```bash
python3.9 -m venv .venv-imodulon
source .venv-imodulon/bin/activate
python -m pip install -r RNAseq_data/iModulon_analysis/preprocesing/requirements.txt
```

After preparing the inputs, run:

```bash
python RNAseq_data/iModulon_analysis/preprocesing/imodulon_activity_pipeline.py \
  --config RNAseq_data/iModulon_analysis/preprocesing/config.json
```

File paths inside `config.json` are resolved relative to that file. The supplied
paths use `../input/` for study data and `../../../annotation_input/` for
reference files.

## Analysis and outputs

For each gene, the pipeline subtracts mean log2(TPM + 1) in the control group
and jointly estimates all 116 module activities using the pseudoinverse of the
reference weight matrix. Group differences are tested using sample-level
ordinary least squares and two-sided t tests.

The study configuration applies Benjamini–Hochberg correction across all
580 module–comparison tests. A differential iModulon activity (DiMA) call
requires `q_primary < 0.05`, `|delta_activity| > 5` and adequate reference-gene
coverage. The study results contain 18 samples, five comparisons and 29 DiMA
calls. Activity changes describe inferred expression signatures; they do not
directly measure regulator activity.

The supplied configuration sets `"output_dir": "../recomputed"`, so new results
go to `RNAseq_data/iModulon_analysis/recomputed/`. Reruns can replace files
there. The retained `DiMA_primary_all.tsv` in the parent folder is not replaced.
Setting `output_dir` to `..` instead would write into that parent folder and
can overwrite its retained results.

| Output | Contents |
| --- | --- |
| `DiMA_primary_all.tsv` | All module–comparison results, including activity differences, uncertainty, p-values, adjusted p-values and DiMA calls |
| `DiMA_primary_significant.tsv` | Results passing the DiMA criteria |
| `A_inferred_sample_activities.tsv` | Inferred activity of each module in each sample |
| `module_projection_QC.tsv`, `sample_projection_QC.tsv`, `contrast_projection_QC.tsv` | Reference coverage and projection diagnostics |
| `figures/` | Diagnostic plots |
| `run_provenance.json` and `run.log` | Configuration, input hashes, package versions and run messages |

The Figure 6C and Table S9 scripts use the complete `DiMA_primary_all.tsv`,
including nonsignificant results. The Figure 6C script reads
`RNAseq_data/iModulon_analysis/DiMA_primary_all.tsv`; adjust its input path if
using results from a separate output folder.

## Adapting the analysis

- Update `expression` and `samples` together. Keep `input_units` as `log2_TPM`
  and `input_pseudocount` as `1` for the upstream log2(TPM + 1) export.
- Set `control_group` and the group names in each contrast's `weights` to match
  the sample sheet. For a treatment-minus-control comparison, use weights
  `1` for the treatment and `-1` for the control.
- Keep `M`, `imodulon_table` and `M_thresholds` consistent. Review
  `min_gene_overlap_fraction` and `min_module_weight_fraction` when checking
  reference coverage.
- Set `fdr`, `fdr_scope` and `min_abs_delta` for the intended calling criteria.
  The activity cutoff depends on the reference weight scale.
- Choose `output_dir` before running. Set `plots` to control diagnostic exports
  and `heatmap_max_modules` to limit the diagnostic heatmap.
