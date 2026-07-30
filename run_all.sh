#!/usr/bin/env bash
# Regenerate every figure, log, and confirmation table from scratch.
# Runtime is about 90 seconds; task7_transmission_scans.py is the slowest single
# step at roughly 13 seconds.
#
# Requires Python 3.10+ and the packages in requirements.txt.
set -euo pipefail
cd "$(dirname "$0")"

run () {
  echo "======================================================================"
  echo "  $1"
  echo "======================================================================"
  python3 "$1"
  echo
}

# 1. State by state: closed forms against an independent Bellman + FOC solve,
#    and the numerical claims that have no proof.
run numerics/task1_11_closedform.py        # (1,1) closed form, corner regime
run numerics/task2_asymmetric.py           # asymmetric fixed point, leader dominance
run numerics/task3_initial_state.py        # (0,0) quadratic, Delta0, interiority
run numerics/task4_comparative_statics.py  # comparative statics, det J, BR regions

# 2. Leader dominance for mu <= lam (Proposition 3): the symbolic proof chain.
run numerics/prop2_dag_identity.py         # equilibrium identity at G_L = G_F = 0
run numerics/prop2_extension_symbolic.py   # affine in mu, both endpoints positive
run numerics/prop2_extension_check.py      # random search for a counterexample

# 3. Figures used in the thesis.
run numerics/task5_paper_figures.py        # beta11, beta00, c/lambda sweeps
run numerics/task8_c_lambda_regen.py       # c/lambda figure variants
run numerics/task9_br_clean_fig.py         # best-response figure for Section 4.2

# 4. Transmission: the Jacobian, the Cramer numerators, and how far the
#    sqrt(2) condition reaches.
run numerics/task7_transmission_scans.py         # sqrt(2) frontier, dz*/dbeta11 map
run numerics/task7_margin_addendum.py            # is sqrt(2) sharp?
run verify_scripts/verify_step8_transmission.py  # symbolic identities, exact zero

# 5. Consolidated confirmation on two independent seeds, each followed by a
#    50-digit mpmath pass over every finite-difference flag.
run numerics/task6_summary.py              # seed 202
run numerics/task6_diagnose_signs.py       # adjudicate the seed-202 flags
run numerics/task6_seed707.py              # seed 707
run numerics/task6_diagnose_seed707.py     # adjudicate the seed-707 flags
run numerics/task6_confirmation.py         # seed 11, machine-readable table

# 6. Symbolic verification of the appendices, in the order they appear in the
#    thesis. Each prints PASS/FAIL per identity; all should be PASS.
run verify_scripts/verify_step2.py         # App. A.1  Bellman -> FOC -> quadratic
run verify_scripts/verify_step2b.py        # App. A.3  root selection, interiority
run verify_scripts/verify_step2c.py        # App. A.3  sufficiency of (A2)
run verify_scripts/verify_step3.py         # App. B.2  effort comparative statics
run verify_scripts/verify_step3c.py        # App. B.2  exact polynomial division
run verify_scripts/verify_step3d.py        # App. B.4  dV11/dbeta11 > 0 (Proposition 2)
run verify_scripts/verify_step4.py         # App. D.1  the (0,0) quadratic
run verify_scripts/verify_step4b.py        # App. D.2, D.3  master quadratic, root selection
run verify_scripts/verify_step4c.py        # App. D.4  small-mu limit
run verify_scripts/verify_step5.py         # App. C.1  best-response quadratics
run verify_scripts/verify_step5d.py        # App. C.4  laggard BR slope
run verify_scripts/verify_step6.py         # App. D    (0,0) comparative statics

echo "Done. Figures in figs/, logs and tables in numerics/."
