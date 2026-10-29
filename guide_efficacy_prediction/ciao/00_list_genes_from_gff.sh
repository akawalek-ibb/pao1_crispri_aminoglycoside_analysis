#!/usr/bin/env bash
# List unique locus tags from gene/CDS/pseudogene features in an NCBI GFF.
# Input: GFF file argument. Output: one locus tag per line on stdout.
# Run from the project root; redirect stdout to a gene-list file.
set -euo pipefail

gff="${1:?usage: $0 genome.gff3 > locus_tags.txt}"

awk -F'\t' '
  $0 ~ /^#/ { next }
  NF < 9 { next }

  # Restrict to gene-like features.
  # If your GFF contains locus_tag only in CDS entries,
  # change this to ($3=="CDS").
  ($3 == "gene" || $3 == "CDS" || $3 == "pseudogene") {

    if (match($9, /(^|;)locus_tag=[^;]+/)) {
      tag = substr($9, RSTART, RLENGTH)
      sub(/(^|;)locus_tag=/, "", tag)
      print tag
    }
  }
' "$gff" \
| sed '/^$/d' \
| sort -u
