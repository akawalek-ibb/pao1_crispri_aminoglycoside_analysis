# CIAO/MERF predictions for multiple genes

This adaptation wraps the published CIAO/MERF guide-design predictor of
Yu et al. (2024; reference below) for parallel assessment of multiple genes.
Each gene runs in a separate process with isolated RNA-folding files. Guide
selection, feature calculations and the released prediction model are retained;
the wrapper does not train a model or change the scores.

Upstream: [BarquistLab/CRISPRi_guide_efficiency_bacteria](https://github.com/BarquistLab/CRISPRi_guide_efficiency_bacteria),
`5_gRNA_design`, commit `2dad9153ea61db49001e2cb48095fbb2eac46003`.

MIT License : Copyright (c) 2023 Yanying Yu, Lars Barquist

## Execution order

| Script | Action |
| --- | --- |
| `00_list_genes_from_gff.sh` | Optional: list locus tags from a GFF |
| `01_fetch_CIAO_model.py` | Verify existing model files or download missing files |
| `02_run_CIAO_batch.py` | Run the predictor for selected genes, optionally in parallel |
| `03_concatenate_CIAO_gRNA_predictions.py` | Combine completed per-gene predictions |

`CRISPRi_design_MERF_progress.py` is the core predictor called automatically by
`02`; it is not an additional manual step.

## Setup and model files

Use the working CIAO environment or create it from `environment.yml`:

```bash
conda env create -f guide_efficacy_prediction/ciao/environment.yml
conda activate ciao_design
python guide_efficacy_prediction/ciao/01_fetch_CIAO_model.py
```

The downloader checks SHA256 hashes against the pinned upstream release.
Existing matching files are reused. Only `CRISPRi_headers.sav` and
`CRISPRi_model.sav` are needed; the larger `Merf_model.sav` is not used.
Model files are retained in this bundle under `saved_model/`, which is
excluded from Git and can be restored with script `01`.

Python package versions for the completed *efp* example are recorded in
`reference_outputs/LESB58_efp/run_config.json`. The supplied environment uses
Python 3.8.12 and scikit-learn 0.24.2, matching the released model.

These executables must also be available on `PATH`:

| Executable | Version used by this adaptation |
| --- | --- |
| `RNAduplex2.1.9h` | ViennaRNA 2.1.9h, RNA–DNA hybridization |
| `RNAduplex2.4.14` | ViennaRNA 2.4.14, RNA–RNA hybridization |
| `RNAfold` | ViennaRNA 2.1.9h |

For installation, see the [upstream instructions](https://github.com/BarquistLab/CRISPRi_guide_efficiency_bacteria#requirements).
The unversioned `RNAfold` must resolve to 2.1.9h; script `02` checks this.

## Run the three-gene example

From the manuscript project root, run:

```bash
python guide_efficacy_prediction/ciao/02_run_CIAO_batch.py \
  --fasta guide_efficacy_prediction/ciao/input/LESB58_reference.fna \
  --gff guide_efficacy_prediction/ciao/input/LESB58_genomic.gff \
  --genes guide_efficacy_prediction/ciao/input/LESB58_example_genes.txt \
  --jobs 3 \
  --output guide_efficacy_prediction/ciao/out_LESB58_three_genes
```

The example selects *efp* (`PLES_RS11455`), *secA* (`PLES_RS24605`) and
*accC* (`PLES_RS26955`). `--jobs 3` permits three concurrent gene processes.
Use fewer jobs if needed. Add `--check-only` to check inputs, dependencies and
model loading without starting predictions.

Each completed gene has `<output>/<locus_tag>/gRNAs.csv` and `run.log`.
Despite its extension, `gRNAs.csv` is tab-separated. Combine the outputs:

```bash
python guide_efficacy_prediction/ciao/03_concatenate_CIAO_gRNA_predictions.py \
  --root guide_efficacy_prediction/ciao/out_LESB58_three_genes \
  --id-column locus_tag --strict-files \
  --output LESB58_three_genes_CIAO.csv
```

The combined CSV is written inside the `--root` folder. Original scores and
guide columns are preserved. **Lower, more negative `Activity score` predicts
stronger depletion.**

The supplied, completed single-gene *efp* run is retained unchanged in
`reference_outputs/LESB58_efp/`; it is a reference example, not the output of
the three-gene command above.

## Other genes or strains

Edit the gene list, or omit `--genes` to process all compatible genes. Choose
a new output folder for each new run. Add `--resume` to reuse completed outputs
from an interrupted run with matching inputs, predictor, environment and settings.
No-guide outcomes and failures are reported separately.

The reference must be a single-record genome FASTA with matching NCBI GFF.
The supplied LESB58 gene inventory contains 6,195 loci; 6,191 have supported
biotypes. To generate an inventory for another reference:

```bash
bash guide_efficacy_prediction/ciao/00_list_genes_from_gff.sh genome.gff > genes.txt
```

Defaults remain 20-nt guides, NGG PAM, 30–85% GC and operon-first flag 0.
`--if-promoter` sets that flag uniformly; it is not inferred from annotations.

For PAO1, retain the concatenator's default `PAO1_ID` column and use
`--output PAO1_all_sgRNA_predictions_CIAO.csv` for the Figure S8 comparison.

## Reference

Yu Y, Gawlitt S, de Andrade e Sousa LB, Merdivan E, Piraud M, Beisel CL,
Barquist L. 2024. Improved prediction of bacterial CRISPRi guide efficiency
from depletion screens through mixed-effect machine learning and data integration.
*Genome Biology* 25:13.
[doi:10.1186/s13059-023-03153-y](https://doi.org/10.1186/s13059-023-03153-y).
