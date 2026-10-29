# RNA-seq reads to gene counts and BigWigs

Use this folder inside `RNAseq_analysis/` in the manuscript project root.
The workflow processes 18 paired-end libraries: three biological replicates
of control1 (the non-targeting control), control1 + TOB, and the four knockdowns. It uses the supplied PAO1
FASTA and counting GTF, not a newly downloaded annotation.

## Sample names and public data

RNA-seq data: [GSE349001](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349001).
`RNAseq_samples.tsv` is the complete run-to-sample map. SRR accessions come
from the supplied SRA Run Selector table; GSM titles, experiment accessions
and original library names were checked against the public GEO sample records.
The original GEO names remain in the map for traceability. Generated filenames
use the GEO-style labels already used by the manuscript figure scripts.

| Sample IDs (three replicates) | SRR accessions in the same order | Condition |
| --- | --- | --- |
| `10_control1`, `11_control1`, `12_control1` | SRR40905244, SRR40905243, SRR40905242 | Induced non-targeting control |
| `13_ftsZg1m4`, `14_ftsZg1m4`, `15_ftsZg1m4` | SRR40905241, SRR40905240, SRR40905239 | ftsZ_g1m4 |
| `16_atpIg1m2`, `17_atpIg1m2`, `18_atpIg1m2` | SRR40905238, SRR40905237, SRR40905236 | atpI_g1m2 |
| `19_accCg1m3`, `20_accCg1m3`, `21_accCg1m3` | SRR40905235, SRR40905234, SRR40905233 | accC_g1m3 |
| `22_secAg1m3`, `23_secAg1m3`, `24_secAg1m3` | SRR40905232, SRR40905231, SRR40905230 | secA_g1m3 |
| `25_control1_TOB`, `26_control1_TOB`, `27_control1_TOB` | SRR40905229, SRR40905228, SRR40905227 | Control + 0.4 µg/mL tobramycin |

Sample IDs and group names match the GEO library names. Manuscript guide
labels such as atpI_g1m2 describe the same construct as the file/group label
`atpIg1m2`. Keep sample-number prefixes when naming replicate files.

## Software

Run in Linux/Bash with the tools on PATH. The reported processing versions
are fastp 0.20.1, Bowtie2 2.4.4, SAMtools 1.13 and featureCounts/Subread 2.0.3.
Read retrieval requires SRA Toolkit (`prefetch`, `fasterq-dump`) and gzip.
BigWig export requires deepTools (`bamCoverage`); its historical version was
not specified in the supplied shell script. Use the recorded environment
when available. Coreutils supplies `realpath`, `readlink` and `ln`.

Script 03 records actual processing-tool versions in
`RNAseq_analysis/out_preprocessing/software_versions.txt`.
Activate your existing environment before running the scripts; they do not
activate Conda environments automatically. For R dependencies, see
[../README_RNAseq_analysis.md](../README_RNAseq_analysis.md).

## Execution order

The steps below rebuild counts and coverage from raw reads. To analyse the
retained count inputs directly, follow
[../README_RNAseq_analysis.md](../README_RNAseq_analysis.md); no read-processing
step is required for those supplied inputs.

Run these commands from the project root. No raw sequencing files are included.
Allow ample disk space for the SRA cache, uncompressed FASTQs, temporary files,
trimmed FASTQs and BAMs. The downloader processes one run at a time.

```bash
bash RNAseq_analysis/preprocessing/01_download_SRA_reads.sh
bash RNAseq_analysis/preprocessing/02_rename_SRA_reads.sh --dry-run
bash RNAseq_analysis/preprocessing/02_rename_SRA_reads.sh
bash RNAseq_analysis/preprocessing/03_preprocess_RNAseq_reads.sh --dry-run
bash RNAseq_analysis/preprocessing/03_preprocess_RNAseq_reads.sh
```

1. **Download:** `prefetch` caches each accession; `fasterq-dump` splits mates
   and gzip compresses them. Complete gzip pairs are reused. Interrupted FASTQ
   conversions require inspection before retrying. SRA archives are retained.
2. **Assign names:** script 02 creates named symbolic links to the downloaded
   FASTQs. It does not delete, move or duplicate the raw read files. If the
   SRR-named gzip pairs are already downloaded, put them in
   `out_preprocessing/download/fastq_sra/` and start at script 02.
3. **Process:** fastp trimming, Bowtie2 alignment, SAMtools mate correction,
   coordinate sorting, duplicate marking/indexing, gene counting and optional
   strand-specific RPKM BigWigs. Script 03 refuses to overwrite existing final
   BAM/count/coverage files. For a separate run, edit its `OUTPUT_DIR`.

Set `MAKE_BIGWIGS=false` in script 03 if only counts are needed. Thread counts
are also settings at the top of that script. The named input FASTQ directory
is independent of its processing-output directory.

## Counting and coverage choices

The supplied counting flags are retained:

```bash
featureCounts -t CDS -s2 -p -B -P -C -T 3 -Q 10 -D 1000 \
  --verbose -M --fraction -a NC_002516.2_fixed2.gtf -o sample.txt sample.markdup.bam
```

**In Subread 2.0.3, `-p` declares paired-end input; without
`--countReadPairs`, reads rather than read pairs are counted.** The change
was introduced in [Subread 2.0.2](https://subread.sourceforge.net/).
The cleaned script preserves the command, including the absence of
`--countReadPairs`; it does not silently change the counting unit. Use the
reported version when reproducing it. Pair/fragment-related switches are
retained as supplied and should not be described as proof of fragment counting.

Bowtie2 uses default alignment reporting, without `-k` or `-a`. MAPQ filtering
remains `-Q 10`; retaining `-M --fraction` does not establish that ambiguous
alignments survive that filter. Duplicates are marked and retained, not removed.
No new rRNA read-removal step is introduced. The R script retains its existing
rRNA-gene exclusions.

Full featureCounts tables and their `.summary` files are kept. The R reader
accepts these tables directly, and `Length` is available for TPM calculation.
The retained two-column count inputs also work, using the same GTF for lengths.
No second stripped count export is needed.

BigWigs retain the supplied settings: 50-bp bins, RPKM normalization and
`--filterRNAstrand forward/reverse`. The library is reverse-stranded; deepTools
interprets these options as transcript strands for dUTP-type libraries.
Coverage does not inherit the featureCounts MAPQ/length filters. The filenames
end in `_fw.bigwig` and `_re.bigwig`, matching Figure 5.

## Outputs and deliberate hand-off

All new preprocessing results stay under `RNAseq_analysis/out_preprocessing/`:

| Folder | Contents |
| --- | --- |
| `download/` | SRA cache, SRR-named FASTQs and conversion temporary directory |
| `fastq_named/` | GEO-style raw-read links |
| `trimmed/` | Filtered paired FASTQs |
| `index/` | Bowtie2 PAO1 index |
| `bam/` | Coordinate-sorted duplicate-marked BAMs and indexes |
| `featureCounts/` | Eighteen full count tables and assignment summaries |
| `coverage/` | Forward/reverse RPKM BigWigs |
| `logs/` | Per-sample fastp reports, alignment/counting logs and coverage logs |

The included R script and retained count inputs have matching GEO-style
filenames. Their existing count values are unchanged. You can run that analysis
without processing reads again:

```r
source("RNAseq_analysis/AMG_RNAseq_featureCounts_edgeR.R")
```

To adopt a NEW raw-read processing run, first compare counts and QC with the
retained reference inputs. Preserve a local copy of the retained inputs, then
copy the new count tables deliberately:

```bash
cp RNAseq_analysis/out_preprocessing/featureCounts/*.txt RNAseq_analysis/input/featureCounts/
```

The R script then writes to `RNAseq_analysis/out/`; files with the same names
there may be overwritten. It does not replace manuscript inputs automatically.
To adopt newly generated complete DE workbooks or BigWigs downstream, copy them
only after reviewing the run:

```bash
cp RNAseq_analysis/out/PAO1_*_All.xlsx RNAseq_data/DE/
cp RNAseq_analysis/out_preprocessing/coverage/*.bigwig RNAseq_data/coverage/
```

R outputs use the existing comparison filenames, including
`PAO1_atpIg1m2_vs_control1_All.xlsx`. Figure 5 and preprocessing 18 can therefore
read an adopted run without renaming those workbooks. For iModulon input,
copy the matching expression matrix and sample sheet together as described
in the RNA-seq analysis README. A raw-read rebuild has not been asserted to
produce byte-identical historical outputs.

Official tool instructions:
[SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump),
[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml),
[featureCounts](https://subread.sourceforge.net/featureCounts.html),
[bamCoverage](https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html).
