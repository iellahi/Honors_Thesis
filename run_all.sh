#!/usr/bin/env bash
# Regenerate every figure and confirmation log from scratch (~40s total).
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
)

for s in "${scripts[@]}"; do
  echo "======================================================================"
  echo "  running numerics/${s}.py"
  echo "======================================================================"
  python3 "numerics/${s}.py"
  echo
done

echo "All scripts complete. Figures in figs/, logs and data in numerics/."
