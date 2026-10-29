#!/usr/bin/env bash
# Run shared preprocessing tasks 03–22 and 25, each in a fresh R session.
# QC (01) and metacounter (02) must be run separately when required.
# Requires Rscript on PATH and the R packages/inputs documented in each script.
# Existing outputs may be overwritten according to each script's settings.
# From the project root: bash preprocessing_scripts/run_03-25.sh
# Comment out any run line below to skip that task. Retain its prepared outputs
# if later tasks require them. MAGeCK results and RNA-seq/iModulon inputs must
# already exist; this runner does not execute those upstream workflows.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
cd -- "${SCRIPT_DIR}/.."
command -v Rscript >/dev/null 2>&1 || { echo "Rscript was not found on PATH." >&2; exit 1; }

run() {
  printf '\nRunning %s\n' "$1"
  Rscript --vanilla "preprocessing_scripts/$1"
  printf 'Finished %s\n' "$1"
}

# Aminoglycoside interactions and functional annotation.
run "03_AMG_relevant_results.R"
run "04_AMG_compare_MAGeCK_Stouffer.R"
run "05_AMG_TU_interactions.R"
run "06_AMG_TU_functional_annotation.R"
run "07_AMG_recurrent_gene_TU_collapse.R"
run "08_AMG_recurrent_KEGG_enrichment.R"

# Raw-count export and shared library representation.
run "09_AMG_export_GEO_raw_counts.R"
run "10_AMG_library_representation.R"

# Vulnerability, MH/FBS effects and published-data comparisons.
run "11_AMG_gene_vulnerability.R"
run "12_AMG_MH_fitness.R"
run "13_AMG_MH_FBS_gene_effects.R"
run "14_AMG_MH_FBS_KEGG.R"
run "15_AMG_MH_FBS_MAGeCK.R"
run "16_AMG_Murray2015_comparison.R"
run "17_AMG_mismatch_interaction_candidates.R"

# RNA-seq, GO and iModulon summaries.
run "18_AMG_RNAseq_results.R"
run "19_AMG_RNAseq_GO_ORA.R"
run "20_AMG_RNAseq_GO_GSEA.R"
run "21_AMG_RNAseq_iModulons.R"
run "22_AMG_RNAseq_selected_GO.R"

# Guide-efficacy model input preparation.
run "25_AMG_guide_efficacy_inputs.R"

printf "\nAll enabled preprocessing scripts completed successfully.\n"
