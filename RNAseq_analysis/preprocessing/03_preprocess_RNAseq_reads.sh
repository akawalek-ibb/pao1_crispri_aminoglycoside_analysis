#!/usr/bin/env bash
# Trim, align, mark duplicates, count genes and generate strand-specific coverage.
# First run scripts 01 and 02, or supply their named paired FASTQs.
# Reference inputs: RNAseq_analysis/input/NC_002516.2.fa and NC_002516.2_fixed2.gtf.
# Requires fastp, Bowtie2, SAMtools, featureCounts/Subread and optionally deepTools.
# All generated files stay in out_preprocessing/; retained R/Figure 5 inputs are
# not replaced. See README.md for the deliberate count/BigWig hand-off.
# Run: bash RNAseq_analysis/preprocessing/03_preprocess_RNAseq_reads.sh
# Check the 18 sample/file assignments without processing: append --dry-run.
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
# User settings: analytical options below retain the supplied shell commands.
THREADS=22
COUNT_THREADS=3
MAKE_BIGWIGS=true
SAMPLE_TABLE="RNAseq_analysis/preprocessing/RNAseq_samples.tsv"
FASTA="RNAseq_analysis/input/NC_002516.2.fa"
GTF="RNAseq_analysis/input/NC_002516.2_fixed2.gtf"
FASTQ_DIR="RNAseq_analysis/out_preprocessing/fastq_named"
OUTPUT_DIR="RNAseq_analysis/out_preprocessing"
INDEX="$OUTPUT_DIR/index/PAO1"
DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then DRY_RUN=true
elif [[ $# -gt 0 ]]; then echo "Usage: $0 [--dry-run]" >&2; exit 1; fi
[[ -s "$SAMPLE_TABLE" && -s "$FASTA" && -s "$GTF" ]] || { echo "Missing sample map, FASTA or GTF." >&2; exit 1; }
# Check the full input list before starting work. Existing BAM/count/coverage
# results are not silently overwritten; use a fresh OUTPUT_DIR for a new run.
while IFS=$'\t' read -r run gsm experiment geo_name sample group rep; do
  [[ "$run" == "run" ]] && continue
  [[ -n "$run" ]] || continue
  printf '%s: %s_R1/R2.fastq.gz -> %s.txt and strand tracks\n' "$run" "$sample" "$sample"
  if ! "$DRY_RUN"; then
    [[ -s "$FASTQ_DIR/${sample}_R1.fastq.gz" && -s "$FASTQ_DIR/${sample}_R2.fastq.gz" ]] || {
      echo "Missing named FASTQs for $sample; run script 02." >&2; exit 1;
    }
    for target in "$OUTPUT_DIR/bam/${sample}.markdup.bam" "$OUTPUT_DIR/featureCounts/${sample}.txt" \
      "$OUTPUT_DIR/coverage/${sample}_fw.bigwig" "$OUTPUT_DIR/coverage/${sample}_re.bigwig"; do
      [[ ! -e "$target" ]] || { echo "Output already exists: $target. Use a fresh OUTPUT_DIR." >&2; exit 1; }
    done
  fi
done < "$SAMPLE_TABLE"
if "$DRY_RUN"; then echo "Dry run: no reads processed or outputs created."; exit 0; fi
for tool in fastp bowtie2 bowtie2-build samtools featureCounts; do
  command -v "$tool" >/dev/null || { echo "Missing executable: $tool" >&2; exit 1; }
done
if "$MAKE_BIGWIGS"; then command -v bamCoverage >/dev/null || { echo "Missing executable: bamCoverage" >&2; exit 1; }; fi
mkdir -p "$OUTPUT_DIR"/{index,trimmed,bam,featureCounts,coverage,logs}
# Record actual executable versions; the README lists the reported reference versions.
{
  fastp --version 2>&1
  bowtie2 --version
  samtools --version
  featureCounts -v 2>&1
  if "$MAKE_BIGWIGS"; then bamCoverage --version; fi
} > "$OUTPUT_DIR/software_versions.txt"
# 1. Build the index from the exact supplied PAO1 reference.
bowtie2-build "$FASTA" "$INDEX" > "$OUTPUT_DIR/logs/bowtie2-build.log" 2>&1
while IFS=$'\t' read -r run gsm experiment geo_name sample group rep; do
  [[ "$run" == "run" ]] && continue
  [[ -n "$run" ]] || continue
  echo "Processing $sample ($run)"
  r1="$OUTPUT_DIR/trimmed/${sample}_R1.filtered.fastq.gz"
  r2="$OUTPUT_DIR/trimmed/${sample}_R2.filtered.fastq.gz"
  # 2. Apply the original fastp defaults to paired reads. Keep raw FASTQs intact.
  fastp --thread "$THREADS" \
    -h "$OUTPUT_DIR/logs/${sample}_fastp.html" -j "$OUTPUT_DIR/logs/${sample}_fastp.json" \
    -i "$FASTQ_DIR/${sample}_R1.fastq.gz" -I "$FASTQ_DIR/${sample}_R2.fastq.gz" \
    -o "$r1" -O "$r2" 2> "$OUTPUT_DIR/logs/${sample}_fastp.log"
  # 3. Bowtie2 default reporting: no -k/-a multimapping expansion. Pipe SAM into
  # name sorting to avoid retaining a large intermediate SAM file.
  bowtie2 -x "$INDEX" -1 "$r1" -2 "$r2" --verbose --threads "$THREADS" \
    2> "$OUTPUT_DIR/logs/${sample}_bowtie2.log" \
    | samtools sort -n --threads "$THREADS" -o "$OUTPUT_DIR/bam/${sample}.namesort.bam" -
  # 4. Fix mate fields, coordinate-sort and MARK duplicates. Neither markdup -r
  # nor featureCounts --ignoreDup is used: duplicate-marked reads are retained.
  samtools fixmate -m --threads "$THREADS" "$OUTPUT_DIR/bam/${sample}.namesort.bam" "$OUTPUT_DIR/bam/${sample}.fixmate.bam"
  samtools sort --threads "$THREADS" -o "$OUTPUT_DIR/bam/${sample}.positionsort.bam" "$OUTPUT_DIR/bam/${sample}.fixmate.bam"
  bam="$OUTPUT_DIR/bam/${sample}.markdup.bam"
  samtools markdup --threads "$THREADS" "$OUTPUT_DIR/bam/${sample}.positionsort.bam" "$bam"
  samtools index "$bam"
  rm -- "$OUTPUT_DIR/bam/${sample}.namesort.bam" "$OUTPUT_DIR/bam/${sample}.fixmate.bam" "$OUTPUT_DIR/bam/${sample}.positionsort.bam"
  # 5. Preserve the supplied counting flags, including -p without --countReadPairs.
  # In Subread >=2.0.2 this specifies paired-end INPUT, but counts reads rather
  # than read pairs. Do not add --countReadPairs when reproducing this command.
  # -Q 10 is retained; -M --fraction does not mean ambiguous reads pass that gate.
  # Keep the full table (including Length) and .summary for R/TPM and QC.
  featureCounts -t CDS -s2 -p -B -P -C -T "$COUNT_THREADS" -Q 10 -D 1000 \
    --verbose -M --fraction -a "$GTF" -o "$OUTPUT_DIR/featureCounts/${sample}.txt" "$bam" \
    > "$OUTPUT_DIR/logs/${sample}_featureCounts.log" 2>&1
  # 6. Use the original 50-bp RPKM coverage settings. Coverage does not inherit
  # featureCounts' MAPQ/length filters; no new filtering is added here.
  if "$MAKE_BIGWIGS"; then
    bamCoverage --bam "$bam" --outFileName "$OUTPUT_DIR/coverage/${sample}_fw.bigwig" \
      --binSize 50 --filterRNAstrand forward --normalizeUsing RPKM \
      --numberOfProcessors "$THREADS" --verbose 2> "$OUTPUT_DIR/logs/${sample}_coverage_fw.log"
    bamCoverage --bam "$bam" --outFileName "$OUTPUT_DIR/coverage/${sample}_re.bigwig" \
      --binSize 50 --filterRNAstrand reverse --normalizeUsing RPKM \
      --numberOfProcessors "$THREADS" --verbose 2> "$OUTPUT_DIR/logs/${sample}_coverage_re.log"
  fi
done < "$SAMPLE_TABLE"
echo "Saved new counts: $OUTPUT_DIR/featureCounts/"
echo "Saved new BigWigs: $OUTPUT_DIR/coverage/"
echo "Review these outputs before replacing the retained inputs; see preprocessing/README.md."
