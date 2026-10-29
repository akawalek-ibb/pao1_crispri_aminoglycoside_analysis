# Pasteur guide-activity predictions

This script scores an existing bacterial CRISPRi guide catalogue with the
embedded Pasteur EcoWG1 coefficients from Calvo-Villamañán et al. (2020;
reference below). Upstream inputs are the guide-design
table and its matching GenBank genome. The supplied example contains all
46 LESB58 *efp* guides (`PLES_RS11455`), in their original order.

## Run the example

Copy `guide_efficacy_prediction/` into the manuscript project root, preserving
the subfolders. Open the R project and run:

```r
source("guide_efficacy_prediction/pasteur/AMG_score_guides_with_Pasteur_model.R")
```

The only additional R package is `data.table`. Coefficients are embedded in the
script; no Python environment or separate coefficient file is required.
The supplied genome is `input/Pseudomonas_aeruginosa_LESB58_125.gbk`.

## Score a full catalogue

In **USER SETTINGS**, change `INPUT_FILE` to the complete guide catalogue,
`GENBANK_FILE` to that strain's matching reference and `OUTPUT_FILE` to the
desired result filename. Source the script again. It processes every input
row, across all genes, without an additional gene-selection step. The example
CSV is a subset; supply a full catalogue to score all genes.

Required columns:

| Column | Content |
| --- | --- |
| `locus_tag` | Gene identifier used for within-gene summaries |
| `sgrna` | 20-nt guide sequence |
| `pam` | Three-base NGG PAM |
| `downstream` | At least 15 downstream bases, in guide/PAM orientation |
| `strand` | Feature strand, `1` or `-1`, using the supplied catalogue's convention |
| `position_offset` | Genomic coordinate using the supplied catalogue's convention; not the gene-relative `offset` |

Retain `record_id`, `feature_start` and `feature_end` when available. The model
uses spacer bases 15–20, the first PAM base and 16 downstream bases. If only
15 downstream bases are supplied, the script recovers the next base from the
genome and checks the adjacent sequence. Coordinates must match that genome.

## Outputs

Results appear in `out_Pasteur_scores/` as
`LESB58_efp_with_Pasteur_scores.tsv` and a status-count TSV. All input rows
and columns are retained. Check `pasteur_score_status`; failed rows have an
`NA` score. With the default settings, coordinate failures do not trigger a
whole-genome sequence search.

Higher `pasteur_activity_score` predicts stronger activity; it is not a measured
repression fraction. Within-gene rank 1 and global quartile 1 are highest.
Global quartiles depend on the input catalogue, so example-only quartiles will
differ from those obtained when scoring the full catalogue.

For Figure S8 G–H, score the full PAO1 catalogue and supply that output at the
Pasteur prediction path configured in those scripts. The LESB58 subset is a
standalone usage example.

## Reference

Calvo-Villamañán A, Ng JW, Planel R, Ménager H, Chen A, Cui L, Bikard D.
2020. On-target activity predictions enable improved CRISPR–dCas9 screens in
bacteria. *Nucleic Acids Research* 48:e64.
[doi:10.1093/nar/gkaa294](https://doi.org/10.1093/nar/gkaa294).

Original tool: [CRISPR Pasteur](https://crispr-browser.pasteur.cloud/guide-rna-design).
