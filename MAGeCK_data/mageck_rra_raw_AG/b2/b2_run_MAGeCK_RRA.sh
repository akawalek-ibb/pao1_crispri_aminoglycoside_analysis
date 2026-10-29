#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

'mageck' 'test' '-k' 'b2_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B2_16_T1_AMK2,B2_17_T1_AMK2,B2_18_T1_AMK2' '-c' 'B2_01_T1_MH,B2_02_T1_MH,B2_03_T1_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b2_AMK_2_T1_vs_None_0_T1' '--normcounts-to-file'
'mageck' 'test' '-k' 'b2_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B2_01_T1_MH,B2_02_T1_MH,B2_03_T1_MH' '-c' 'B2_25_T0_MH,B2_26_T0_MH,B2_27_T0_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b2_None_0_T1_vs_None_0_T0' '--normcounts-to-file'
'mageck' 'test' '-k' 'b2_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B2_16_T1_AMK2,B2_17_T1_AMK2,B2_18_T1_AMK2' '-c' 'B2_25_T0_MH,B2_26_T0_MH,B2_27_T0_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b2_AMK_2_T1_vs_None_0_T0' '--normcounts-to-file'
