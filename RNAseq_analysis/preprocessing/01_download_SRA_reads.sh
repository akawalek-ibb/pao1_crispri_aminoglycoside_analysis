#!/usr/bin/env bash
# Download the 18 runs in RNAseq_samples.tsv and convert them to paired FASTQs.
# Requires SRA Toolkit (prefetch, fasterq-dump) and gzip. Run from the project root:
# bash RNAseq_analysis/preprocessing/01_download_SRA_reads.sh
# Outputs stay in RNAseq_analysis/out_preprocessing/download/. Original SRA names
# are retained here; script 02 creates the named FASTQ links used by script 03.
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
THREADS=8
SAMPLE_TABLE="RNAseq_analysis/preprocessing/RNAseq_samples.tsv"
SRA_DIR="RNAseq_analysis/out_preprocessing/download/sra"
FASTQ_DIR="RNAseq_analysis/out_preprocessing/download/fastq_sra"
TMP_DIR="RNAseq_analysis/out_preprocessing/download/tmp"
for tool in prefetch fasterq-dump gzip; do
  command -v "$tool" >/dev/null || { echo "Missing executable: $tool" >&2; exit 1; }
done
mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$TMP_DIR"
while IFS=$'\t' read -r run gsm experiment geo_name sample group rep; do
  [[ "$run" == "run" ]] && continue
  [[ -n "$run" ]] || continue
  [[ "$run" =~ ^SRR[0-9]+$ ]] || { echo "Invalid run: $run" >&2; exit 1; }
  r1="$FASTQ_DIR/${run}_1.fastq.gz"
  r2="$FASTQ_DIR/${run}_2.fastq.gz"
  if [[ -s "$r1" && -s "$r2" ]]; then
    echo "Already converted: $run ($sample)"
    continue
  fi
  # Do not combine a partial conversion with a different set of mate files.
  if [[ -e "$r1" || -e "$r2" || -e "$FASTQ_DIR/${run}_1.fastq" || -e "$FASTQ_DIR/${run}_2.fastq" ]]; then
    echo "Partial FASTQ output for $run; inspect/remove its incomplete FASTQs before retrying." >&2
    exit 1
  fi
  echo "Downloading $run -> $sample"
  prefetch "$run" --output-directory "$SRA_DIR"
  fasterq-dump "$SRA_DIR/$run" --split-files --threads "$THREADS" \
    --temp "$TMP_DIR" --outdir "$FASTQ_DIR"
  [[ -s "$FASTQ_DIR/${run}_1.fastq" && -s "$FASTQ_DIR/${run}_2.fastq" ]] || {
    echo "Missing paired FASTQs for $run" >&2; exit 1;
  }
  gzip "$FASTQ_DIR/${run}_1.fastq" "$FASTQ_DIR/${run}_2.fastq"
done < "$SAMPLE_TABLE"
echo "Next: bash RNAseq_analysis/preprocessing/02_rename_SRA_reads.sh"
