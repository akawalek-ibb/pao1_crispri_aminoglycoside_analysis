# External sources and credits

## External code, models and analytical methods

| Repository file or workflow | Source and credit | What is included here |
| --- | --- | --- |
| `guide_efficacy_prediction/ciao/CRISPRi_design_MERF_progress.py` | Yu et al. (2024), *Improved prediction of bacterial CRISPRi guide efficiency from depletion screens through mixed-effect machine learning and data integration*. [doi:10.1186/s13059-023-03153-y](https://doi.org/10.1186/s13059-023-03153-y); [BarquistLab source](https://github.com/BarquistLab/CRISPRi_guide_efficiency_bacteria/tree/2dad9153ea61db49001e2cb48095fbb2eac46003/5_gRNA_design). | Adaptation of upstream `CRISPRi_design_MERF.py`, with local execution/progress changes.CIAO models downloaded by `guide_efficacy_prediction/ciao/01_fetch_CIAO_model.py` | Same Yu et al. source and pinned commit; upstream `5_gRNA_design/saved_model/`. | External `CRISPRi_headers.sav` and `CRISPRi_model.sav`, fetched separately. |
| `guide_efficacy_prediction/pasteur/AMG_score_guides_with_Pasteur_model.R` | Calvo-Villamañán et al. (2020), *On-target activity predictions enable improved CRISPR–dCas9 screens in bacteria*. [doi:10.1093/nar/gkaa294](https://doi.org/10.1093/nar/gkaa294); [EcoWG1 repository](https://gitlab.pasteur.fr/dbikard/ecowg1). | Study R implementation with embedded coefficients from upstream [`reg_coef.pkl`](https://gitlab.pasteur.fr/dbikard/ecowg1/-/raw/master/reg_coef.pkl). |
| `preprocessing_scripts/02_AMG_metacounter.R` and the fitted gene-vulnerability workflow | Ward et al. (2024), *Essential gene knockdowns reveal genetic vulnerabilities and antibiotic sensitivities in Acinetobacter baumannii*. [doi:10.1128/mbio.02051-23](https://doi.org/10.1128/mbio.02051-23); [publication-linked code](https://github.com/jasonpeterslab/acinetobacter_baumannii_CRISPRi_seq). | Adapted analysis framework, with study-specific filtering, aggregation and vulnerability rules. It is not an identical reproduction of the Ward implementation. |
| Mismatch-guide design and relative-efficacy annotations in `annotation_input/PAO1_key.tsv` and `library_design/historical_outputs/` | Hawkins et al. (2020), *Mismatch-CRISPRi reveals the co-varying expression–fitness relationships of essential genes in Escherichia coli and Bacillus subtilis*. [doi:10.1016/j.cels.2020.09.009](https://doi.org/10.1016/j.cels.2020.09.009); [mismatch-CRISPRi code](https://github.com/traeki/mismatch_crispri); [guide designer](https://github.com/traeki/sgrna_design). | Study guide designs and relative-efficacy predictions using the published framework. |
| `RNAseq_data/iModulon_analysis/preprocesing/imodulon_activity_pipeline.py` | [PyModulon](https://github.com/SBRG/pymodulon), particularly `pymodulon.util.infer_activities`; Sastry et al. (2024), *iModulonMiner and PyModulon: Software for unsupervised mining of gene expression compendia*. [doi:10.1371/journal.pcbi.1012546](https://doi.org/10.1371/journal.pcbi.1012546). | Study pipeline implementing the projection with NumPy and testing activity changes using biological replicates. |


`Figure_scripts/FigS8/input_AMG_FigS8/PAO1_all_sgRNA_predictions_CIAO.csv`
and `PAO1_sgRNA_with_Pasteur_scores.tsv.gz` were generated for this study
using those external models.

## Published essentiality and antimicrobial-fitness evidence

The following files are under `library_design/input/` unless another path is given. 

| File | Published source | Local representation or limitation |
| --- | --- | --- |
| `pnas.1900570116.sd05.xlsx` | Poulsen et al. (2019), *Defining the core essential genome of Pseudomonas aeruginosa*, Dataset S5. [doi:10.1073/pnas.1900570116](https://doi.org/10.1073/pnas.1900570116). | Essentiality catalogue with local additional annotation. |
| `pnas.1900570116.sd06.xlsx` | Poulsen et al. (2019), Dataset S6, at the DOI above. Also credits the original Lee and Liberati studies | Immediate source for the Lee/Liberati evidence read by the collector. |
| Lee calls extracted from that Dataset S6 | Lee et al. (2015), *General and condition-specific essential functions of Pseudomonas aeruginosa*. [doi:10.1073/pnas.1422186112](https://doi.org/10.1073/pnas.1422186112). |
| Liberati calls extracted from that Dataset S6 | Liberati et al. (2006), *An ordered, nonredundant library of Pseudomonas aeruginosa strain PA14 transposon insertion mutants*. [doi:10.1073/pnas.0511100103](https://doi.org/10.1073/pnas.0511100103). |
| `Turner_2015_S1_short.xlsx` | Turner et al. (2015), *Essential genome of Pseudomonas aeruginosa in cystic fibrosis sputum*, Table S1. [doi:10.1073/pnas.1419677112](https://doi.org/10.1073/pnas.1419677112). | Local extract retaining the initial, MOPS and synthetic-sputum evidence columns. |
| `Skurnik_2013_table_short.xlsx` | Skurnik et al. (2013), *A comprehensive analysis of in vitro and in vivo genetic fitness of Pseudomonas aeruginosa using high-throughput sequencing of transposon libraries*. [doi:10.1371/journal.ppat.1003582](https://doi.org/10.1371/journal.ppat.1003582). | Local PA14 gene-list extract, mapped for the PAO1 analysis. |
| `Gallagher_2011_short.xlsx` | Gallagher, Shendure and Manoil (2011), *Genome-scale identification of resistance functions in Pseudomonas aeruginosa using Tn-seq*. [doi:10.1128/mBio.00315-10](https://doi.org/10.1128/mBio.00315-10). | Local tobramycin-resistance evidence extract.|
| `Murray_2015mbo005152522sd2.xlsx`; `annotation_input/Murray_2015_mbo005152522sd2.xlsx` | Murray et al. (2015), *Intrinsic Antimicrobial Resistance Determinants in the Superbug Pseudomonas aeruginosa*, Data Set S2. [doi:10.1128/mBio.01603-15](https://doi.org/10.1128/mBio.01603-15). |
| `Table_sum.XLSX` | Vitale et al. (2020), *Identification of Genes Required for Resistance to Peptidomimetic Antibiotics by Transposon Sequencing*. [doi:10.3389/fmicb.2020.01681](https://doi.org/10.3389/fmicb.2020.01681). | Local PAO1-identifier shortlist from a study performed in PA14.|

### Reference memberships for Figure 2D

`annotation_input/AMG_published_fitness_membership.tsv` is a study-assembled
membership table with 320 Poulsen, 580 Zhang and 299 Kaczmarczyk members.
Its columns refer to:

- Poulsen et al. (2019), core essential catalogue, DOI above.
- Zhang et al. (2025), *CRISPRi screen identifies FprB as a synergistic target
  for gallium therapy in Pseudomonas aeruginosa*.
  [doi:10.1038/s41467-025-61208-z](https://doi.org/10.1038/s41467-025-61208-z).
- Kaczmarczyk et al. (2025), *Genome-wide high-density CRISPR interference
  screens reveal condition-specific metabolic vulnerabilities in Pseudomonas
  aeruginosa PAO1*. The cited source is the bioRxiv preprint
  [doi:10.1101/2025.08.12.669819](https://doi.org/10.1101/2025.08.12.669819).
  The retained 299-gene set is the shared intersection across LB, glucose
  and succinate comparisons, as specified by the Figure 2D script.

The original Kaczmarczyk comparison tables and the extraction script for this
membership file are not included. 

## Shared annotations and other published data

| Repository file | Source and credit | Version, transformation or information still needed |
| --- | --- | --- |
| `library_design/input/Pseudomonas_aeruginosa_PAO1_107.csv` | [Pseudomonas Genome Database](https://www.pseudomonas.com/); Winsor et al. (2016), *Enhanced annotations and features for comparing thousands of Pseudomonas genomes in the Pseudomonas genome database*. [doi:10.1093/nar/gkv1227](https://doi.org/10.1093/nar/gkv1227). | The file explicitly records PGD version 22.1, assembly GCF_000006765.1 and chromosome NC_002516.2. |
| `annotation_input/PAO1_curated.tsv`; `library_design/input/PAO1.csv` | Study-curated annotation tables incorporating PGD annotations and other references. | These merge coordinates/products, orthology, pathways, operons/TUs and study-specific categories. They are not single unmodified PGD downloads. The operon/TU source, release and curation steps require confirmation. |
| `annotation_input/PAO1_GO_unified_propagated.tsv.gz` | [Gene Ontology Consortium](https://geneontology.org/docs/go-citation-policy/) and the providers of the original gene annotations. 
| `annotation_input/PAO1_pathways_all_withCOG.tsv` | Mixed-source mapping containing [KEGG](https://www.kegg.jp/kegg/), [COG](https://www.ncbi.nlm.nih.gov/research/cog/), eggNOG, PseudoCAP, PseudoCyc and MetaCyc annotations, as identified in the file's `Pathway_Database` column. | 
| `annotation_input/pathoadaptive_mutations_science.adi0908_table_s5.csv` | Weimann et al. (2024), *Evolution and host-specific adaptation of Pseudomonas aeruginosa*, Supplementary Table S5. [doi:10.1126/science.adi0908](https://doi.org/10.1126/science.adi0908). | Local CSV representation.|
| `annotation_input/precise411-M.csv`, `precise411-im_table.csv`, `precise411-M_thresholds.csv` | Rajput et al. (2022), *Advanced transcriptomic analysis reveals the role of efflux pumps and media composition in antibiotic responses of Pseudomonas aeruginosa*. [doi:10.1093/nar/gkac743](https://doi.org/10.1093/nar/gkac743); [PRECISE411 on iModulonDB](https://imodulondb.org/p_aeruginosa/precise411); [source repository](https://github.com/akanksha-r/modulome_paeru2.0). | Reference for 116 modules. The retained M matrix uses named module columns; the annotation table includes local convenience fields.|


## Reference genomes and counting annotations

| Repository file | Reference source | Local status |
| --- | --- | --- |
| `RNAseq_analysis/input/NC_002516.2.fa` | PAO1 chromosome [NC_002516.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_002516.2), RefSeq assembly [GCF_000006765.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000006765.1/). Genome publication: Stover et al. (2000), [doi:10.1038/35023079](https://doi.org/10.1038/35023079). | 
| `RNAseq_analysis/input/NC_002516.2_fixed2.gtf`; `annotation_input/PAO1_orfs.bed` | PAO1 annotation matching NC_002516.2. | Locally prepared counting/coordinate annotations. |
| `guide_efficacy_prediction/ciao/input/LESB58_reference.fna`; `LESB58_genomic.gff` | LESB58 chromosome [NC_011770.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_011770.1), assembly [GCF_000026645.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000026645.1/). Genome publication: Winstanley et al. (2009), [doi:10.1101/gr.086082.108](https://doi.org/10.1101/gr.086082.108). | The GFF header records RefSeq annotation `GCF_000026645.1-RS_2025_08_25`, dated 2025-08-25. |
| `guide_efficacy_prediction/pasteur/input/Pseudomonas_aeruginosa_LESB58_125.gbk` | Same LESB58 assembly; PGD GenBank-format export. Credit the genome publication and PGD above. |
