# Rebuilding analyses

Run from the project root containing `AMG_cleaned_Git.Rproj`. 

The supplied bundle includes raw barcode counts, batch designs, fitted screen
exports in `edgeR_data/`, RNA-seq DE workbooks, BigWigs, DiMA results and the saved
guide model. Figure/table scripts that read shared `preprocessing_scripts/out_*/`
tables can reuse compatible prepared outputs already present. Tables and intermediate
files can be generated using the preparation steps first.

## Pooled CRISPRi screens

The supplied fitted exports are `edgeR_data/b*_melted_results.tsv.gz` and
`edgeR_data/b*_median_melted_results.tsv.gz`, for batches b1–b4. They can be
reused. Run
`preprocessing_scripts/02_AMG_metacounter.R` only when rebuilding the fits from
raw counts/designs; it replaces fitted results in `edgeR_data/`. Script
`01_AMG_sample_QC_PCA.R` reviews samples but does not remove them or change designs.

The following filenames are under `preprocessing_scripts/`:

| Task/script | Inputs and downstream use |
| --- | --- |
| `03_AMG_relevant_results.R` | Complete fitted screen exports → drug calls, recurrence, thresholds and guide subsets. Needed for Figure 3, S11–S14, Table S8 analyses and Figure 6A/B |
| `04_AMG_compare_MAGeCK_Stouffer.R` | Script 03 plus the complete supplied MAGeCK catalogue → method-comparison tables for S8 and Figure S12 |
| `05_AMG_TU_interactions.R` | Script 03 with guide export enabled → pooled-guide TU interaction tests for S8 |
| `06_AMG_TU_functional_annotation.R` | Supplied curated/pathway annotations → gene/TU memberships for scripts 07–08 and S8 |
| `07_AMG_recurrent_gene_TU_collapse.R` | Scripts 03 and 06 → grouping recurrent gene hits by TU for S8 |
| `08_AMG_recurrent_KEGG_enrichment.R` | Scripts 03 and 06 → recurrent-hit KEGG overrepresentation for S8 |
| `09_AMG_export_GEO_raw_counts.R` | Unfiltered barcode counts, designs and guide key → count matrix, sample key and reference comparison |
| `10_AMG_library_representation.R` | Script 09 plus annotation → shared Table S4/Figure S6 summaries |
| `11_AMG_gene_vulnerability.R` | Fitted guide effects and annotation → per-comparison vulnerability, consensus, bins and curves; fits may be slow |
| `12_AMG_MH_fitness.R` | Fitted gene effects, script 11 and published membership reference → S5 and Figures 2C/D/G, S7/S9 |
| `13_AMG_MH_FBS_gene_effects.R` | B3 fitted gene effects and script 11 → S7 and Figure S10 |
| `14_AMG_MH_FBS_KEGG.R` | Script 13 plus KEGG annotation → S7 and Figure 2H |
| `15_AMG_MH_FBS_MAGeCK.R` | Complete supplied MAGeCK catalogue → S7 FBS comparison; independent of scripts 13–14 |
| `16_AMG_Murray2015_comparison.R` | Script 03, Murray workbook and orthology → published-data comparison for S8 |
| `17_AMG_mismatch_interaction_candidates.R` | Script 03 guide exports → exploratory candidates for S8. Figure S13 separately calculates its sensitizing-only analysis |
| `18_AMG_RNAseq_results.R` | Five complete DE workbooks → standardized RNA-seq results for Figure 6 and S9 |
| `19_AMG_RNAseq_GO_ORA.R` | Script 18 and GO reference → directional GO overrepresentation |
| `20_AMG_RNAseq_GO_GSEA.R` | Script 18 and GO reference → ranked GO enrichment; script 19 is not required |
| `21_AMG_RNAseq_iModulons.R` | Complete DiMA table and atlas → S9 iModulon tables; not needed by standalone Figure 6C |
| `22_AMG_RNAseq_selected_GO.R` | Scripts 18–19 and selected GO input → Figure 6D and S9 display summaries |
| `25_AMG_guide_efficacy_inputs.R` | Raw counts/designs and guide key → annotated count/context matrices and design copies for S8A–F |

To run all the preparation tasks, use a terminal at the project root:
```bash
bash preprocessing_scripts/run_03-25.sh
```

The runner lists tasks 03–22 and 25; The runner sets the project root as the working
directory but does not load `.Rproj` settings or interactive RStudio objects.
QC (01), metacounter (02), MAGeCK execution and upstream RNA-seq/iModulon analyses and MAGecK 
are separate steps. Reruns will most often replace outputs according to each script's settings.

### Raw reads and barcode-count

Raw FASTQs are external inputs. Download/conversion and renaming are described
by `Barcode_counting/seq_data/sequencing_reads_howtoget_info.txt` and
`rename_sra_fastq.sh`.

`Barcode_counting/count_Lib1_barcodes.sh` writes new count tables under
`Barcode_counting/out_count_Lib1_barcodes/b*/`. Downstream R scripts use the
retained historical count tables in `Barcode_counting/count_tables/`. 
Batch-prefixed sample IDs, such as `B1_1-T0_MH`, are used in both the read basenames and the design `condition` fields.

To rebuild analyses from newly counted reads, compare the completed tables with
the retained inputs, confirm that their third-column sample IDs match the batch
designs, and manually copy the required `b*_overall_stats.tsv` files into
`Barcode_counting/count_tables/` before running downstream R scripts.

### MAGeCK

See [README_MAGeCK.md](MAGeCK_data/README_MAGeCK.md). Run
`MAGeCK_data/01_MAGeCK_preprocessing_and_run.R`, activate the MAGeCK environment
in a terminal, execute the generated master shell script, then run
`MAGeCK_data/02_MAGeCK_postprocessing.R`.

The complete supplied hand-off catalogue can be used directly by scripts 04/15:
`MAGeCK_data/mageck_rra_postprocessing_AG/MAGeCK_drug_contrasts_all_genes.tsv.gz`

## RNA-seq rebuild and frozen inputs

[RNAseq_analysis](RNAseq_analysis/README_RNAseq_analysis.md) starts from the 18
included per-sample featureCounts files and writes only to `RNAseq_analysis/out/`.
The supplied two-column counts require the included counting GTF for TPM lengths.
The R script explicitly selects 18 GEO-style sample filenames; each sample
should have one selected count file. The
[RNA-seq README](RNAseq_analysis/README_RNAseq_analysis.md) lists the sample
names, input formats and output files.

For rebuilding from raw reads, follow the
[RNA-seq preprocessing instructions](RNAseq_analysis/preprocessing/README.md):
`01_download_SRA_reads.sh` → `02_rename_SRA_reads.sh` →
`03_preprocess_RNAseq_reads.sh`, all in `RNAseq_analysis/preprocessing/`.
These steps download paired reads, create named links, trim and align reads,
process BAMs, run featureCounts and optionally export strand-specific RPKM
BigWigs. New results stay in `RNAseq_analysis/out_preprocessing/`.

After reviewing a new run, copy its count tables deliberately into
`RNAseq_analysis/input/featureCounts/` before running the DEG analysis. DE
exports use the GEO-style comparison filenames expected by script 18 and
Figure 5, including `PAO1_atpIg1m2_vs_control1_All.xlsx`. To adopt them downstream,
copy reviewed complete DE workbooks from `RNAseq_analysis/out/` to
`RNAseq_data/DE/`, and BigWigs from
`RNAseq_analysis/out_preprocessing/coverage/` to `RNAseq_data/coverage/`.
No automatic replacement of retained manuscript inputs occurs.

The [iModulon pipeline](RNAseq_data/iModulon_analysis/preprocesing/README_iModulon.md)
uses log2(TPM + 1), matching samples and PRECISE411 references. To use a new
RNA-seq run, copy both `CRISPRi_silencing_RNAseq_log2_TPM.tsv` and
`CRISPRi_silencing_RNAseq_iModulon_samples.tsv` from `RNAseq_analysis/out/` into
`RNAseq_data/iModulon_analysis/input/`, or explicitly update their paths in
`config.json`. Figure 6C and script 21 still read the retained `DiMA_primary_all.tsv` one directory above.

## Supplementary-table assembly

| Table | Preparation required before its `Table_scripts/AMG_TableS*_*.R` assembler |
| --- | --- |
| S4 | 09 → 10 |
| S5 | 11 → 12 |
| S7 | 11 → 13 → 14, plus 15 using supplied MAGeCK results |
| S8 | 03 → 04/05; 06 → 07/08 with 03 also available; 16/17 after 03 |
| S9 | 18 → 19/20; 22 after 18/19; independent 21 using supplied DiMA |

Each assembler reads its editable legend TSV from `Table_scripts/`. Preserve
those files, the reference annotations and all-gene inputs. Output paths and
overwrite behavior are set within each script. The R package inventory and
session information in the project root describe the supplied software environment.
