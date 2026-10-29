#!/usr/bin/env bash
# Move downloaded SRR read pairs to the filenames used by barcode counting.
# Run from Barcode_counting/seq_data/ after conversion and gzip compression.
# Preview: bash rename_sra_fastq.sh --dry-run
# Rename:  bash rename_sra_fastq.sh
# Optional final argument: path to sequencing_reads_howtoget_info.txt.
# Requires Bash 4+, awk and mv. No sequencing files are downloaded here.

set -euo pipefail

DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    shift
fi
(( $# <= 1 )) || { printf '%s\n' 'Usage: bash rename_sra_fastq.sh [--dry-run] [mapping-file.txt]' >&2; exit 1; }

MAP_FILE="${1:-$(dirname -- "${BASH_SOURCE[0]}")/sequencing_reads_howtoget_info.txt}"
FASTQ_DIR="$PWD/downloaded_fastq"
TARGET_DIR="$PWD"

die() { printf '%s\n' "$*" >&2; exit 1; }
[[ -r "$MAP_FILE" ]] || die "Cannot read the sample map: $MAP_FILE"
[[ -d "$FASTQ_DIR" ]] || die "FASTQ directory not found: $FASTQ_DIR. Run from seq_data/."

# Comment-table columns: # GSM sample_header original_sample_ID SRR.
mapfile -t SAMPLE_MAP < <(awk '$1 == "#" && $2 ~ /^GSM[0-9]+$/ {print $5 "\t" $4}' "$MAP_FILE")
(( ${#SAMPLE_MAP[@]} > 0 )) || die "No GSM-to-SRR sample rows found in $MAP_FILE."

declare -A seen_runs=() seen_samples=()
sources=()
targets=()
already_named=0

# Check every available pair and destination before moving any files.
for entry in "${SAMPLE_MAP[@]}"; do
    read -r run sample extra <<< "$entry"
    [[ "$run" =~ ^SRR[0-9]+$ && "$sample" =~ ^B[1-4]_[A-Za-z0-9_.-]+$ && -z "${extra:-}" ]] ||
        die "Invalid sample-map row: $entry"
    [[ -z "${seen_runs[$run]:-}" && -z "${seen_samples[$sample]:-}" ]] ||
        die "Duplicate run or sample in the map: $entry"
    seen_runs[$run]=1
    seen_samples[$sample]=1

    batch="${sample%%_*}"
    batch="${batch,,}"
    dest1="$TARGET_DIR/$batch/${sample}_R1.fastq.gz"
    dest2="$TARGET_DIR/$batch/${sample}_R2.fastq.gz"

    source1=""
    source2=""
    for stem in "$run" "$run.sra" "$run.sralite"; do
        [[ ! -e "$FASTQ_DIR/${stem}_1.fastq" && ! -e "$FASTQ_DIR/${stem}_2.fastq" ]] ||
            die "Uncompressed reads found for $run; finish compression before renaming."
        r1="$FASTQ_DIR/${stem}_1.fastq.gz"
        r2="$FASTQ_DIR/${stem}_2.fastq.gz"
        if [[ -e "$r1" || -e "$r2" ]]; then
            [[ -z "$source1" ]] || die "Multiple FASTQ sets found for $run; keep one matched pair."
            [[ -s "$r1" && -s "$r2" ]] || die "Missing or empty R1/R2 mate for $run. No files were moved."
            source1="$r1"
            source2="$r2"
        fi
    done

    if [[ -z "$source1" ]]; then
        # Runs not downloaded are skipped. Completed renames can be left in place.
        if [[ -s "$dest1" && -s "$dest2" ]]; then
            already_named=$((already_named + 1))
        elif [[ -e "$dest1" || -L "$dest1" || -e "$dest2" || -L "$dest2" ]]; then
            die "Incomplete destination pair for $sample. No files were moved."
        fi
        continue
    fi

    [[ ! -e "$dest1" && ! -L "$dest1" && ! -e "$dest2" && ! -L "$dest2" ]] ||
        die "Destination exists for $sample. No files were moved."
    sources+=("$source1" "$source2")
    targets+=("$dest1" "$dest2")
done

for i in "${!sources[@]}"; do
    printf '%s -> %s\n' "${sources[$i]}" "${targets[$i]}"
    if ! "$DRY_RUN"; then
        mkdir -p -- "$(dirname -- "${targets[$i]}")"
        mv -n -- "${sources[$i]}" "${targets[$i]}"
        [[ ! -e "${sources[$i]}" ]] || die "Move did not complete; destination may have appeared: ${targets[$i]}"
    fi
done

if "$DRY_RUN"; then
    printf 'Preview: %d pairs would be renamed; %d named pairs already present.\n' "$(( ${#sources[@]} / 2 ))" "$already_named"
else
    printf 'Renamed %d pairs; %d named pairs already present.\n' "$(( ${#sources[@]} / 2 ))" "$already_named"
fi
