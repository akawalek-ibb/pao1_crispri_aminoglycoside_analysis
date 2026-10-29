# RNA-seq analysis from featureCounts

`AMG_RNAseq_featureCounts_edgeR.R` analyses per-sample gene counts with edgeR,
produces expression plots and exports log2(TPM + 1) values for iModulon analysis.
Place `RNAseq_analysis/` in the project root and run the script from that root.

## Rebuild from SRA reads

The three scripts in `preprocessing/` download the 18 SRA runs, assign the
GEO-style filenames, and run fastp, Bowtie2, SAMtools, featureCounts and
strand-specific BigWig export. See [preprocessing/README.md](preprocessing/README.md)
for commands, dependencies and the manual hand-off of newly generated files.
The supplied count inputs below remain usable without downloading raw reads.

## Inputs

Place the 18 count files in `RNAseq_analysis/input/featureCounts/`.
The experiment has six conditions, each with three biological replicates.

| Count files | Group in `sample_info` | Condition |
| --- | --- | --- |
| `10_control1.txt`, `11_control1.txt`, `12_control1.txt` | `control1` | Induced non-targeting control |
| `13_ftsZg1m4.txt`, `14_ftsZg1m4.txt`, `15_ftsZg1m4.txt` | `ftsZg1m4` | ftsZ_g1m4 |
| `16_atpIg1m2.txt`, `17_atpIg1m2.txt`, `18_atpIg1m2.txt` | `atpIg1m2` | atpI_g1m2 |
| `19_accCg1m3.txt`, `20_accCg1m3.txt`, `21_accCg1m3.txt` | `accCg1m3` | accC_g1m3 |
| `22_secAg1m3.txt`, `23_secAg1m3.txt`, `24_secAg1m3.txt` | `secAg1m3` | secA_g1m3 |
| `25_control1_TOB.txt`, `26_control1_TOB.txt`, `27_control1_TOB.txt` | `control1_TOB` | Non-targeting control plus tobramycin |

Use one count column per file. Standard featureCounts tables must contain
`Geneid` and have the sample counts in the final column. Two-column tables
with gene identifiers followed by counts are also accepted. Use unique gene
identifiers consistent with the annotation. Fractional counts are retained.

The annotation file is `annotation_input/PAO1_curated.tsv`. It must contain
`PAO1_ID`, `start` and `end`; gene names, products and other available annotations
are included in the result tables.

TPM calculation uses the `Length` column in the featureCounts files. For
two-column inputs without lengths, supply the GTF used for counting at
`RNAseq_analysis/input/NC_002516.2_fixed2.gtf`. The script reconstructs lengths
from the union of CDS intervals for each `gene_id`, matching counting with
`-t CDS`. A different annotation may give different lengths.

## Run

With the project root as the working directory, run in R or RStudio:

```r
source("RNAseq_analysis/AMG_RNAseq_featureCounts_edgeR.R")
```

Required R packages: edgeR, limma, dplyr, readr, tidyr, openxlsx, data.table,
pheatmap, RColorBrewer, ggrepel, ggplot2, forcats, tibble, reshape2 and svglite.
The script records package versions in `RNAseq_sessionInfo.txt`.

## Analysis and outputs

The script merges counts, removes the 12 rRNA loci listed in `core_rRNA_tags`
and calculates FPKM-like values for PCA and expression-distribution plots.
Differential expression is tested from counts using `filterByExpr`, TMM
normalization, dispersion estimation and edgeR likelihood-ratio tests. Each
knockdown and the tobramycin-treated control are compared with `control1`.
Benjamini–Hochberg correction is applied separately within each comparison.

`_All` tables contain all tested genes. `_FC2` tables select genes with
FDR < 0.05 and FoldChange > 2 or < 0.5, equivalent to |log2FC| > 1.
Volcano classifications and heatmap gene selection use the same strict cutoff.
Heatmaps and cluster plots are generated when enough genes remain for plotting.

All outputs are saved in `RNAseq_analysis/out/`. Rerunning overwrites files
with the same names there.

| Output | Contents |
| --- | --- |
| `CRISPRi_RNAseq_merged_featureCounts.tsv` | Merged counts before rRNA removal |
| `CRISPRi_silencing_RNAseq_raw_counts_after_rRNA_filter.csv` | Counts after rRNA removal |
| `CRISPRi_silencing_RNAseq_edgeR_counts_after_filterByExpr.csv` | Counts retained for differential-expression testing |
| `CRISPRi_silencing_RNAseq_FPKM_like_values.csv` | Length-normalized expression used for diagnostic plots |
| `PAO1_*_All.xlsx` and `.csv` | Complete results for each of the five comparisons |
| `PAO1_*_FC2.xlsx` and `.csv` | Results passing the FDR and fold-change cutoffs |
| `CRISPRi_silencing_RNAseq_log2_TPM.tsv` | log2(TPM + 1) matrix for iModulon analysis |
| `CRISPRi_silencing_RNAseq_iModulon_samples.tsv` | Matching sample identifiers and condition groups |
| `CRISPRi_silencing_RNAseq_iModulon_length_audit.tsv` | Gene lengths, their sources and inclusion in the TPM matrix |
| Plot files and `RNAseq_sessionInfo.txt` | PCA, expression distributions, volcanoes, heatmaps, cluster plots and package versions |

TPM is calculated after rRNA removal and before edgeR expression filtering.
Genes without valid lengths are omitted from TPM and recorded in the length
audit. The edgeR analysis can continue if TPM export fails; use the TPM matrix
only when the current run reports a successful export.

## Adapting the analysis

- **Files and samples:** edit `featurecount_dir` and `sample_info`. Each row
  assigns a filename, sample identifier and condition (`strain`). Update the
  factor levels when changing conditions.
- **Annotation and exclusions:** edit `annotation_file`, `counting_gtf_file`
  and `core_rRNA_tags` to match the organism and counting annotation.
- **Comparisons:** update `contrast_matrix` and `contrast_renames` together.
  These determine the tested comparisons and their output filenames.
- **Cutoffs:** update the `_FC2` export filter, volcano classification and
  heatmap settings together. The heatmap block also defines `k_clusters`.
- **Plots:** set `PRINT_PLOTS_IN_RSTUDIO`, `SAVE_PNG`, `SAVE_PDF` and `SAVE_SVG`
  near the top of the script.
- **iModulon analysis:** update `imodulon_group_map` when changing condition
  names. Its group labels must match the Python pipeline configuration.

The count filenames, edgeR groups/contrasts and iModulon sample sheet use the
labels `control1`, `control1_TOB`, `ftsZg1m4`, `atpIg1m2`, `accCg1m3` and `secAg1m3`.
These are the GEO-style library names; manuscript guide labels may include
an underscore, for example atpI_g1m2.

To use the TPM exports in the iModulon pipeline, copy the matrix and sample
sheet to `RNAseq_data/iModulon_analysis/input/`, or point the Python configuration
to their locations in `RNAseq_analysis/out/`. For downstream edgeR analyses,
set input paths and comparison names to match the exported `_All` tables.
