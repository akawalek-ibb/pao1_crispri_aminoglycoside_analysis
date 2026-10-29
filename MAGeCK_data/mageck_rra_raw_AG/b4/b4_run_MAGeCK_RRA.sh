#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

'mageck' 'test' '-k' 'b4_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B4_16_T1_TOB0.025,B4_17_T1_TOB0.025' '-c' 'B4_01_T1_MH,B4_02_T1_MH,B4_03_T1_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b4_TOB_0.025_T1_vs_None_0_T1' '--normcounts-to-file'
'mageck' 'test' '-k' 'b4_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B4_18_T1_TOB0.1,B4_19_T1_TOB0.1' '-c' 'B4_01_T1_MH,B4_02_T1_MH,B4_03_T1_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b4_TOB_0.1_T1_vs_None_0_T1' '--normcounts-to-file'
'mageck' 'test' '-k' 'b4_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B4_20_T1_TOB1,B4_21_T1_TOB1' '-c' 'B4_01_T1_MH,B4_02_T1_MH,B4_03_T1_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b4_TOB_1_T1_vs_None_0_T1' '--normcounts-to-file'
'mageck' 'test' '-k' 'b4_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B4_01_T1_MH,B4_02_T1_MH,B4_03_T1_MH' '-c' 'B4_28_T0_MH,B4_29_T0_MH,B4_30_T0_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b4_None_0_T1_vs_None_0_T0' '--normcounts-to-file'
'mageck' 'test' '-k' 'b4_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B4_16_T1_TOB0.025,B4_17_T1_TOB0.025' '-c' 'B4_28_T0_MH,B4_29_T0_MH,B4_30_T0_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b4_TOB_0.025_T1_vs_None_0_T0' '--normcounts-to-file'
'mageck' 'test' '-k' 'b4_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B4_18_T1_TOB0.1,B4_19_T1_TOB0.1' '-c' 'B4_28_T0_MH,B4_29_T0_MH,B4_30_T0_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b4_TOB_0.1_T1_vs_None_0_T0' '--normcounts-to-file'
'mageck' 'test' '-k' 'b4_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B4_20_T1_TOB1,B4_21_T1_TOB1' '-c' 'B4_28_T0_MH,B4_29_T0_MH,B4_30_T0_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b4_TOB_1_T1_vs_None_0_T0' '--normcounts-to-file'
