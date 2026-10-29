# CRISPRi analysis in *Pseudomonas aeruginosa*

Analysis scripts and supporting data for pooled PAO1 CRISPRi fitness and
aminoglycoside-response screens, guide-performance modeling, RNA-seq and
individual validation assays.

**Due to space constraints, frozen supplementary tables and figure PDFs are
provided as reference outputs. Underlying numerical figure data and
associated outputs that can be regenerated using the supplied scripts.**

## Getting started

Open `AMG_cleaned_Git.Rproj` in RStudio and use the project root as the working
directory. Start with the figure or table of interest. Prepared inputs allow
all outputs to be reproduced without rerunning upstream analyses.
Script naming follows the table and figure order in the manuscript.  

Each script defines its input paths, settings and output location. Script
headers describe required inputs and upstream dependencies. See
[REBUILDING.md](REBUILDING.md) when changing inputs or rebuilding analyses.

## Repository contents

| Folder | Contents |
| --- | --- |
| `Barcode_counting/` | Read retrieval, barcode counting and count tables |
| `annotation_input/` | Gene, transcriptional-unit, pathway and reference annotations |
| `preprocessing_scripts/` | Screen summaries, enrichment analyses and prepared figure/table inputs |
| `edgeR_data/` | Processed guide-level and gene-level screen results and annotated raw counts |
| `Figure_scripts/` | Main and supplementary figure scripts |
| `Table_scripts/` | Supplementary-table assemblers and legends |
| `RNAseq_analysis/` | Differential-expression analysis |
| `RNAseq_data/` | RNA-seq inputs, coverage and iModulon analysis |
| `guide_efficacy_prediction/` | Saved guide model and candidate-scoring example |
| `library_design/` | Essentiality/fitness evidence and library records |
| `MAGeCK_data/` | MAGeCK preparation, execution and comparison with Stouffer calls |

Folder-specific READMEs provide input formats and run instructions.
Supplementary-table assemblers cover S4, S5, S7, S8 and S9.

## Software

The supplied R environment records describe R 4.5.2 and Bioconductor 3.22 on
Ubuntu 22.04.5 LTS (x86_64). Key packages include edgeR 4.8.2, limma 3.66.0,
poolr 1.2-0, glmnet 5.1, mgcv 1.9-4, clusterProfiler 4.18.4 and fgsea 1.36.2.
Full records are in `R_sessionInfo.txt` and `R_package_versions.tsv`; the
installed-package list is a snapshot, not a list of required installations.

[SOFTWARE.md](SOFTWARE.md) separates requirements by workflow. iModulon
inference uses a Python environment with packages pinned in
`RNAseq_data/iModulon_analysis/preprocesing/requirements.txt`; the recorded
reference used Python 3.9.13. MAGeCK and external guide-prediction tools use
separate environments. They are needed only when regenerating their outputs.

## Data and rebuilding

Pooled-screen CRISPRi sequencing data are available under
[GSE349013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349013).
Read-download and sample-renaming instructions are in `Barcode_counting/seq_data/`.

RNA-seq data are available under
[GSE349001]. Read-download, sample-renaming and preprocessing instructions are in RNAseq_analysis/preprocessing/.

Analyses using processed screen results expect the batch-specific exports in
`edgeR_data/`. If these are not included in a downloaded bundle, regenerate
them from the count tables and experimental designs with
`preprocessing_scripts/02_AMG_metacounter.R`, following [REBUILDING.md](REBUILDING.md).

