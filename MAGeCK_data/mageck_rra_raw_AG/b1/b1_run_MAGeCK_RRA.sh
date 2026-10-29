#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

'mageck' 'test' '-k' 'b1_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B1_05_T1_TOB0.4,B1_06_T1_TOB0.4' '-c' 'B1_03_T1_MH,B1_04_T1_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b1_TOB_0.4_T1_vs_None_0_T1' '--normcounts-to-file'
'mageck' 'test' '-k' 'b1_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B1_09_T2_TOB0.4,B1_10_T2_TOB0.4' '-c' 'B1_07_T2_MH,B1_08_T2_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b1_TOB_0.4_T2_vs_None_0_T2' '--normcounts-to-file'
'mageck' 'test' '-k' 'b1_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B1_03_T1_MH,B1_04_T1_MH' '-c' 'B1_01_T0_MH,B1_02_T0_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b1_None_0_T1_vs_None_0_T0' '--normcounts-to-file'
'mageck' 'test' '-k' 'b1_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B1_07_T2_MH,B1_08_T2_MH' '-c' 'B1_01_T0_MH,B1_02_T0_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b1_None_0_T2_vs_None_0_T0' '--normcounts-to-file'
'mageck' 'test' '-k' 'b1_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B1_05_T1_TOB0.4,B1_06_T1_TOB0.4' '-c' 'B1_01_T0_MH,B1_02_T0_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b1_TOB_0.4_T1_vs_None_0_T0' '--normcounts-to-file'
'mageck' 'test' '-k' 'b1_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B1_09_T2_TOB0.4,B1_10_T2_TOB0.4' '-c' 'B1_01_T0_MH,B1_02_T0_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b1_TOB_0.4_T2_vs_None_0_T0' '--normcounts-to-file'
