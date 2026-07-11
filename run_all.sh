#!/usr/bin/env bash
# Regenerate every figure and confirmation log from scratch (a few minutes total).
set -euo pipefail
cd "$(dirname "$0")"

scripts=(
  task1_11_closedform
  task2_asymmetric
  task3_initial_state
  task4_comparative_statics
  task5_paper_figures
  task6_summary
  task6_confirmation
  task6_diagnose_signs
  verify_dx_dV11
  task7_transmission_scans
  task8_c_lambda_regen
  task9_br_clean_fig
)

for s in "${scripts[@]}"; do
  echo "======================================================================"
  echo "  running numerics/${s}.py"
  echo "======================================================================"
  python3 "numerics/${s}.py"
  echo
done

echo "======================================================================"
echo "  running verify_scripts/verify_step8_transmission.py"
echo "======================================================================"
python3 "verify_scripts/verify_step8_transmission.py"

echo "All scripts complete. Figures in figs/, logs and data in numerics/."
