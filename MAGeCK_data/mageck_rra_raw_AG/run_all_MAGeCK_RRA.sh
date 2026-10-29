#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

'mageck' '--version' 2>&1 | tee "$SCRIPT_DIR/MAGeCK_version.txt"

bash "$SCRIPT_DIR/b1/b1_run_MAGeCK_RRA.sh"
bash "$SCRIPT_DIR/b2/b2_run_MAGeCK_RRA.sh"
bash "$SCRIPT_DIR/b3/b3_run_MAGeCK_RRA.sh"
bash "$SCRIPT_DIR/b4/b4_run_MAGeCK_RRA.sh"
