#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ============================================================
# Lib1 CRISPRi barcode counting
#
# Purpose: Filter paired reads for the guide cassette, count exact
# barcode matches with SEAL, and export barcode/count/condition rows.
# Inputs: Lib1_barcodes.fasta, Lib1_barcodes.tsv, paired FASTQ.gz files.
# Software: Bash, cutadapt, BBMap (seal.sh), Java, awk and sort.
# Run: BATCH=b1 bash Barcode_counting/count_Lib1_barcodes.sh
# Reads: Barcode_counting/seq_data/b1/<sample>_R1/R2.fastq.gz
# Outputs: out_count_Lib1_barcodes/b1/ beside this script.
# Existing count_tables/ files are never modified.
# Upstream: supply FASTQ pairs externally; seq_data/sequencing_reads_howtoget_info.txt
# and seq_data/rename_sra_fastq.sh describe the download/renaming route.
# The output/input directory separation is intentional: downstream R scripts
# use the retained historical count tables in count_tables/b*_overall_stats.tsv.
# Recounting writes to out_count_Lib1_barcodes/b*/ and does not replace those inputs.
# Both FASTQ basenames and design condition fields use batch-prefixed sample IDs
# (e.g. B1_1-T0_MH). To use a recount downstream, compare it with the retained
# table and manually copy the completed b*_overall_stats.tsv into count_tables/
# after confirming that its third-column sample IDs match the batch design.
# ============================================================

# USER SETTINGS

# Paths are anchored to the script, independently of the working directory.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

# Select one batch per run: b1, b2, b3 or b4.
BATCH="${BATCH:-b1}"

# Override READS_DIR if FASTQs are stored elsewhere.
READS_DIR="${READS_DIR:-${SCRIPT_DIR}/seq_data/${BATCH}}"
BARCODES_FASTA="${SCRIPT_DIR}/Lib1_barcodes.fasta"
BARCODES_TSV="${SCRIPT_DIR}/Lib1_barcodes.tsv"
OUTPUT_DIR="${SCRIPT_DIR}/out_count_Lib1_barcodes/${BATCH}"

# Computational settings.
CUTADAPT_THREADS="${CUTADAPT_THREADS:-24}"
SEAL_THREADS="${SEAL_THREADS:-4}"
SEAL="${SEAL:-}"  # Optional executable path; otherwise discover seal.sh.

# Retain filtered FASTQs by default; false removes them after a successful run.
KEEP_FILTERED_READS="${KEEP_FILTERED_READS:-true}"

# true validates files and software without running analysis or making outputs.
CHECK_ONLY="${CHECK_ONLY:-false}"

# INTERNAL CONSTANTS AND HELPERS

# SEAL statistics are divided by two when exported. The count/pair convention
# must match the SEAL output used; it is not checked against FASTQ contents here.
COUNT_DIVISOR=2

fail() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
trap 'printf "ERROR: Failed at line %s. Inspect the working directory printed above.\n" "$LINENO" >&2' ERR

# VALIDATE INPUTS AND SOFTWARE

case "${BATCH}" in b1|b2|b3|b4) ;; *) fail "BATCH must be b1, b2, b3 or b4." ;; esac
for flag in "${KEEP_FILTERED_READS}" "${CHECK_ONLY}"; do
  case "${flag}" in true|false) ;; *) fail "Boolean settings must be true or false." ;; esac
done
for threads in "${CUTADAPT_THREADS}" "${SEAL_THREADS}"; do
  [[ "${threads}" =~ ^[1-9][0-9]*$ ]] || fail "Thread counts must be positive integers."
done

[[ -s "${BARCODES_FASTA}" ]] || fail "Missing or empty barcode FASTA: ${BARCODES_FASTA}"
# The TSV reference is required by this input check; counting itself uses the FASTA.
[[ -s "${BARCODES_TSV}" ]] || fail "Missing or empty barcode TSV: ${BARCODES_TSV}"
[[ -d "${READS_DIR}" ]] || fail "Reads directory not found: ${READS_DIR}. Supply FASTQs or set READS_DIR."

R1_FILES=("${READS_DIR}"/*_R1.fastq.gz)
R2_FILES=("${READS_DIR}"/*_R2.fastq.gz)
(( ${#R1_FILES[@]} > 0 )) || fail "No *_R1.fastq.gz files found in ${READS_DIR}."

SAMPLES=()
for r1 in "${R1_FILES[@]}"; do
  r2="${r1%_R1.fastq.gz}_R2.fastq.gz"
  [[ -s "${r1}" ]] || fail "Empty R1 file: ${r1}"
  [[ -s "${r2}" ]] || fail "Missing or empty R2 file: ${r2}"
  sample="$(basename -- "${r1}" _R1.fastq.gz)"
  [[ -n "${sample}" ]] || fail "FASTQ filenames must include a sample identifier."
  SAMPLES+=("${sample}")
done
for r2 in "${R2_FILES[@]}"; do
  [[ -s "${r2%_R2.fastq.gz}_R1.fastq.gz" ]] || fail "R2 has no nonempty R1 partner: ${r2}"
done

for executable in cutadapt java awk sort mktemp; do
  command -v "${executable}" >/dev/null 2>&1 || fail "Required software not found: ${executable}"
done
if [[ -n "${SEAL}" ]]; then
  SEAL="$(command -v "${SEAL}")" || fail "SEAL executable not found."
elif command -v seal.sh >/dev/null 2>&1; then
  SEAL="$(command -v seal.sh)"
elif [[ -x /usr/share/bbmap/seal.sh ]]; then
  SEAL=/usr/share/bbmap/seal.sh
else
  fail "seal.sh not found. Add BBMap to PATH or set SEAL."
fi
[[ -x "${SEAL}" ]] || fail "SEAL is not executable: ${SEAL}"

printf 'Batch %s: %s paired samples; SEAL: %s\n' "${BATCH}" "${#SAMPLES[@]}" "${SEAL}"
if [[ "${CHECK_ONLY}" == true ]]; then
  printf 'Input paths, pair names and software checks passed. FASTQ contents were not checked.\n'
  exit 0
fi

# Work in a new directory so failed runs cannot masquerade as completed outputs.
# Each run regenerates results; existing outputs are not used as cached inputs.
mkdir -p -- "${OUTPUT_DIR}"
RUN_DIR="$(mktemp -d "${OUTPUT_DIR}/.run_XXXXXX")"
printf 'Working directory: %s\n' "${RUN_DIR}"

# FILTER AND COUNT EACH SAMPLE

for sample in "${SAMPLES[@]}"; do
  printf 'Processing %s\n' "${sample}"
  r1="${READS_DIR}/${sample}_R1.fastq.gz"
  r2="${READS_DIR}/${sample}_R2.fastq.gz"
  out1="${RUN_DIR}/${sample}_cut_R1.fastq.gz"
  out2="${RUN_DIR}/${sample}_cut_R2.fastq.gz"

  # Select exact cassette matches with the configured pair-filter rule.
  cutadapt \
    -g TAGTNNNNNNNNNNNNNNNNNNNNGTTT \
    -G AAACNNNNNNNNNNNNNNNNNNNNACTA \
    --overlap 28 \
    --discard-untrimmed \
    --pair-filter both \
    --cores="${CUTADAPT_THREADS}" \
    --no-trim \
    --no-indels \
    -e 0 \
    -o "${out1}" \
    -p "${out2}" \
    "${r1}" "${r2}" \
    > "${RUN_DIR}/${sample}_cutadapt.txt" 2>&1

  "${SEAL}" \
    overwrite=t \
    in1="${out1}" \
    in2="${out2}" \
    k=20 mm=f hdist=0 ambiguous=toss \
    t="${SEAL_THREADS}" \
    trd=t rename=t fbm=t int=f \
    ref="${BARCODES_FASTA}" \
    stats="${RUN_DIR}/${sample}.seal.tsv" \
    > "${RUN_DIR}/${sample}.seal.log" 2>&1

  [[ -s "${RUN_DIR}/${sample}.seal.tsv" ]] || fail "No SEAL statistics produced for ${sample}; see ${RUN_DIR}."
done

# EXPORT COUNT TABLES

# Headerless three-column format matches the original overall_stats.tsv.
# Iterate only over this run's samples to exclude unrelated statistics files.
OVERALL="${RUN_DIR}/${BATCH}_overall_stats.tsv"
{
  for sample in "${SAMPLES[@]}"; do
    awk -v sample="${sample}" -v divisor="${COUNT_DIVISOR}" '
      BEGIN { FS="\t"; OFS="\t" }
      !/^#/ {
        if ($1 != "#Matched" && $1 != "#Name" && $1 != "Name") {
          print $1, $2 / divisor, sample
        }
      }
    ' "${RUN_DIR}/${sample}.seal.tsv" | sort -k1,1
  done
} > "${OVERALL}"
[[ -s "${OVERALL}" ]] || fail "Count table is empty; inspect ${RUN_DIR}."

# Keep the summary header first and sort the sample rows.
SUMMARY="${RUN_DIR}/${BATCH}_sample_count_summary.tsv"
{
  printf 'condition\ttotal_barcode_counts\n'
  awk 'BEGIN { FS="\t"; OFS="\t" } { counts[$3] += $2 }
       END { for (s in counts) print s, counts[s] }' "${OVERALL}" | sort -k1,1
} > "${SUMMARY}"

if [[ "${KEEP_FILTERED_READS}" == false ]]; then
  for sample in "${SAMPLES[@]}"; do
    rm -- "${RUN_DIR}/${sample}_cut_R1.fastq.gz" "${RUN_DIR}/${sample}_cut_R2.fastq.gz"
  done
fi

# Promote outputs only after every sample and table completes successfully.
for result in "${RUN_DIR}"/*; do
  mv -- "${result}" "${OUTPUT_DIR}/"
done
rmdir -- "${RUN_DIR}"

# Remove retained filtered reads for these samples from earlier successful runs.
if [[ "${KEEP_FILTERED_READS}" == false ]]; then
  for sample in "${SAMPLES[@]}"; do
    rm -f -- "${OUTPUT_DIR}/${sample}_cut_R1.fastq.gz" "${OUTPUT_DIR}/${sample}_cut_R2.fastq.gz"
  done
fi

printf 'Completed %s.\nCounts: %s/%s_overall_stats.tsv\nSummary: %s/%s_sample_count_summary.tsv\n' \
  "${BATCH}" "${OUTPUT_DIR}" "${BATCH}" "${OUTPUT_DIR}" "${BATCH}"
