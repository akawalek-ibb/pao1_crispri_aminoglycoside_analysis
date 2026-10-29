#!/usr/bin/env bash
# Assign GEO-style sample filenames without moving or duplicating raw reads.
# Input: script 01's SRR*_1.fastq.gz / SRR*_2.fastq.gz and RNAseq_samples.tsv.
# Output: named symbolic links in RNAseq_analysis/out_preprocessing/fastq_named/.
# Preview: bash RNAseq_analysis/preprocessing/02_rename_SRA_reads.sh --dry-run
# Apply:   bash RNAseq_analysis/preprocessing/02_rename_SRA_reads.sh
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
SAMPLE_TABLE="RNAseq_analysis/preprocessing/RNAseq_samples.tsv"
SOURCE_DIR="RNAseq_analysis/out_preprocessing/download/fastq_sra"
OUTPUT_DIR="RNAseq_analysis/out_preprocessing/fastq_named"
DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then DRY_RUN=true
elif [[ $# -gt 0 ]]; then echo "Usage: $0 [--dry-run]" >&2; exit 1; fi
# Validate the complete map and every available destination before making links.
while IFS=$'\t' read -r run gsm experiment geo_name sample group rep; do
  [[ "$run" == "run" ]] && continue
  [[ -n "$run" ]] || continue
  for mate in 1 2; do
    source="$SOURCE_DIR/${run}_${mate}.fastq.gz"
    target="$OUTPUT_DIR/${sample}_R${mate}.fastq.gz"
    [[ -s "$source" ]] || { echo "Missing input: $source" >&2; exit 1; }
    if [[ -e "$target" || -L "$target" ]]; then
      [[ -L "$target" && "$(readlink -f "$target")" == "$(realpath "$source")" ]] || {
        echo "Destination already contains another file: $target" >&2; exit 1;
      }
    fi
    printf '%s -> %s\n' "$source" "$target"
  done
done < "$SAMPLE_TABLE"
if "$DRY_RUN"; then echo "Preview only: no links created."; exit 0; fi
mkdir -p "$OUTPUT_DIR"
while IFS=$'\t' read -r run gsm experiment geo_name sample group rep; do
  [[ "$run" == "run" ]] && continue
  [[ -n "$run" ]] || continue
  for mate in 1 2; do
    target="$OUTPUT_DIR/${sample}_R${mate}.fastq.gz"
    [[ -L "$target" ]] || ln -s "$(realpath "$SOURCE_DIR/${run}_${mate}.fastq.gz")" "$target"
  done
done < "$SAMPLE_TABLE"
echo "Next: bash RNAseq_analysis/preprocessing/03_preprocess_RNAseq_reads.sh"
