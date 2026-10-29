# Editable Figure 3 functional annotations

`manual_functional_annotations.tsv` contains three columns:

| Column | Contents |
| --- | --- |
| `PAO1_ID` | Unique locus identifier used to match the screen results |
| `gene_name` | Displayed gene label |
| `annotation_category` | Functional class used for the annotation colour |

Edit gene labels and categories in this table. Row order does not control
heatmap ordering or clustering. Another TSV can be used if it has the same
columns and unique locus identifiers.

`functional_annotation_palette.tsv` maps category labels to colours. For a new
category, add its colour there or allow the plotter to assign one. Category
spelling must match between the two tables.

The supplied table labels PA2172 as `Unresolved category (PA2172)`, with peach
colour `#FA8C57`. Update the annotation and palette together when assigning its
biological category.
