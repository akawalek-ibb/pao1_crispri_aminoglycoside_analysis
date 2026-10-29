#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

'mageck' 'test' '-k' 'b3_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B3_22_T1_FBS,B3_23_T1_FBS,B3_24_T1_FBS' '-c' 'B3_01_T1_MH,B3_02_T1_MH,B3_03_T1_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b3_FBS_60_T1_vs_None_0_T1' '--normcounts-to-file'
'mageck' 'test' '-k' 'b3_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B3_01_T1_MH,B3_02_T1_MH,B3_03_T1_MH' '-c' 'B3_25_T0_MH,B3_26_T0_MH,B3_27_T0_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b3_None_0_T1_vs_None_0_T0' '--normcounts-to-file'
'mageck' 'test' '-k' 'b3_MAGeCK_fullmatch_raw_counts.tsv' '-t' 'B3_22_T1_FBS,B3_23_T1_FBS,B3_24_T1_FBS' '-c' 'B3_25_T0_MH,B3_26_T0_MH,B3_27_T0_MH' '--norm-method' 'median' '--gene-lfc-method' 'median' '-n' 'results/b3_FBS_60_T1_vs_None_0_T0' '--normcounts-to-file'
