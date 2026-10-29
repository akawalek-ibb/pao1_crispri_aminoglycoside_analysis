# Software and environments

## R analyses

The supplied records describe R 4.5.2, Bioconductor 3.22 and Ubuntu 22.04.5
LTS (x86_64). Key analysis packages are edgeR 4.8.2, limma 3.66.0,
poolr 1.2-0, glmnet 5.1, mgcv 1.9-4, clusterProfiler 4.18.4 and fgsea 1.36.2.
Plotting includes ggplot2 4.0.3 and ComplexHeatmap 2.26.1.

Script headers identify direct requirements. The supplied
session loaded Matrix 1.7-6. BLAS/LAPACK details are retained in the session file.

These records describe the supplied environment, rather than certifying the
environment of every historical analysis.

## Other workflows

| Workflow | Requirements and version record |
| --- | --- |
| iModulon inference | Python 3.9.13 is the recorded reference; NumPy 1.26.4, pandas 2.3.3, SciPy 1.13.1, matplotlib 3.9.4 and seaborn 0.13.2 are pinned in `RNAseq_data/iModulon_analysis/preprocesing/requirements.txt`|
| MAGeCK RRA | Best is to use a separate environment. MAGeCK version 0.5.9.5` |
| Barcode counting | Bash, cutadapt 3.5, BBMap 38.95 and Java OpenJDK 11.0.13, plus standard shell utilities. Java is required here independently of R. |
| SRA read retrieval | SRA Toolkit (`prefetch` and `fasterq-dump`) and gzip. Version 2.11.3 was used.|
| RNA-seq read processing | featureCounts/Subread v2.0.3. The R analysis can start from featureCounts tables|
| CRISPRi-Pasteur and CIAO/MERF predictions | Prepared predictions suffice for downstream comparisons.

