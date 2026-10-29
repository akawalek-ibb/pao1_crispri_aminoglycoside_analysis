# Library-design evidence and guide records

From the project root, run:

```r
source("library_design/AMG_collect_library_design_evidence.R")
```

The script collects nine essentiality/antimicrobial-fitness sources from
`library_design/input/` and the orthology mapping from
`annotation_input/PA14_PA01_orthologs.csv`. It uses
`annotation_input/PAO1_curated.tsv` for additional annotation when available.
Required R packages are `data.table` and `readxl`. Source calls and mapping information are written as
TSV/TSV.gz files in `library_design/out_AMG_collect_library_design_evidence/`.
This evidence-collection step does not redesign guides or create a formatted
Table S2 workbook. No downstream script automatically replaces the curated
guide key or published reference memberships with these outputs.
Edit the source paths and settings in the script to adapt the collection.

The library guide records are in `library_design/historical_outputs/`:

| File | Records | Unique spacers |
| --- | ---: | ---: |
| `Lib1_FLguides.tsv.gz` | 26,656 locus assignments | 26,602 |
| `Lib1_Mismatched_guides.tsv.gz` | 54,800 locus assignments | 54,783 |
| `Lib1_Control_guides.tsv.gz` | 1,693 | 1,693 |
| `Lib1_all_guides.tsv.gz` | 83,078 adapter-containing 78-nt oligos | 83,078 |

The all-guides file contains 78-nt oligos: a 29-nt prefix, the 20-nt spacer
and a 29-nt suffix. Its central spacers define the same 83,078-guide set as
the analysis key. It is not a spacer-only barcode reference.

A spacer assigned to multiple loci can appear on separate rows in the
assignment tables. `annotation_input/PAO1_key.tsv` defines the primary target
assignment used in screen summaries. For activity prediction, use the
[guide-scoring workflow](../guide_efficacy_prediction/this_study/README.md).
