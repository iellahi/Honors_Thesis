# Numerical Validation — Summary Memo

**Model.** Discrete-time, two-firm, two-stage innovation race; simultaneous breakthroughs, coin-flip
tiebreaker, direct `(0,0)→(1,1)` jump. Maintained assumptions: **(A1)** `x∈[0,1]`, `λ,μ<½`;
**(A2)** `λW≤c`; **(A2′)** `μW≤c`.
**Scope.** Tasks 1–9, reusing the `verify_scripts/` solvers (`verify_step5b.py` reference).
Scripts + logs in `/numerics`, figures in `/figs`. No `.tex` files modified.
Tasks 1–6: confirmation of every theorem and numerical claim. Tasks 7–9 + Step 8: the
transmission package (Section 5.5) and figure updates — see §4 below.

**Bottom line.** Within the maintained assumptions, **every theorem passes and every numerical claim
holds with zero violations.** The two apparent sign "failures" in the first pass were finite-difference
artifacts, confirmed clean by 50-digit mpmath re-checks (0 real violations). Outside the assumptions,
each claim fails exactly where the theory predicts it can. An **independent consolidated re-scan**
(`task6_confirmation.py`, seed 11, N=60,000) reproduces every result with **0 hard violations**
(min Δ₀=0.0213; provocation 99.92%; `dx*/dβ₁₁` 100% positive), corroborating the seed-202 tables below
under a different seed. A post-edit re-run against the 7/11 draft (`task6_rerun_7-11.txt`) reproduces
everything again: 0 hard violations.

**Status upgrades since the first pass (now in the 7/11 draft, Appendix C.5).** Three claims that
were numerical-only in Tasks 1–6 are now theorems: `det J > 0`, uniqueness of the asymmetric fixed
point, and `dx*/dV₁₁ > 0` conditional on `x* < √2·z*` (which implies `dx*/dβ₁₁ > 0` at `μ=λ`).
The numerics' job for these moved from sign-checking to verifying the `√2` condition region-wide (§4).

---

## 1. Confirmation table — THEOREMS (numerically cross-checked, in-region)

Scan: seed 202, N = 60,000 draws in (A1)/(A2)/(A2′). "Real viol." = violations surviving a
50-digit mpmath re-check of every flagged case.

| # | Result | Check | Float flags | Real viol. | Status |
|---|---|---|---|---|---|
| Lem. δ | `Δ=(u−c)²+8c²>0`; form `=9c²−2cu+u²` | closed form | 0 | 0 | ✅ |
| Lem. root | subtracted root `Mx<1`; additive `Mx>1` | closed form | 0 | 0 | ✅ |
| Lem. key | `√Δ>3c−u` (`Δ−(3c−u)²=4cu`) | closed form | 0 | 0 | ✅ |
| Prop 1 | `x₁₁*` signs `(W,λ,c,β₁₁)=(+,+,−,−)` | finite diff | 0 | 0 | ✅ |
| Prop 1′ | `V₁₁` signs `(β₁₁,W,c,λ)=(+,+,+,−)` | FD + **mpmath** | 221 | **0** | ✅ |
| Lem. ord | `0≤V₀₁≤V₁₁<W/2≤V₁₀≤W` | closed form | 0 | 0 | ✅ |
| Prop 2 | `z*>x*` at `μ=λ` | fixed point | 0 | 0 | ✅ |
| §3 BR | laggard BR slope `sign(X′)=sign(N)` | FD vs identity | 0 | 0 | ✅ |
| Prop 3 | `x₀₀*` signs `(V₁₀,V₁₁,V₀₁)=(+,+,−)` | FD + **mpmath** | 18 | **0** | ✅ |
| Lem. master | `F₀₀[terminal subs]−F₁₁ ≡ 0` | **sympy** exact | — | 0 | ✅ |
| §4 BR | additive-root BRs solve `G_L=G_F=0` | residual | 0 | 0 | ✅ |

`min Δ = 0.0229`, `min Δ₀ = 0.0203` over the scan.
**Note on the two "flag" rows.** Prop 1′ flags are all the `dV₁₁/dc` slot where `V₁₁` is near its
`c→∞` saturation (true derivative ~1e-6, below float FD resolution). Prop 3 flags are all the
`dx₀₀*/dV₁₁` slot at very small `μ` (≈0.01), where `F₀₀,V₁₁=2μm(1+β₀₀)x ∝ μ²` is tiny. mpmath at
50 digits recovers the correct signs `(+,+,+,−)` and `(+,+,−)` in every flagged case.

## 2. Confirmation table — NUMERICAL claims (in-region)

| Claim | Result | Prior | Status |
|---|---|---|---|
| `dx*/dβ₁₁ > 0` (laggard catch-up) | **100.0%** (60000/60000) | ~100% | ✅ |
| `dz*/dβ₁₁ < 0` (leader, ambiguous) | **95.56%** negative | ~96% | ✅ |
| `dx*/dλ` total-vs-direct reversals | **1.76%** | rare | ✅ |
| provocation region `N>0` | **99.84–99.92%** | ~99.8% | ✅ |
| `Δ₀>0` | 0 violations | 0 | ✅ |
| `x₀₀*∈(0,1)`, `mx₀₀*<1` | 0 violations | 0 | ✅ |
| `det J > 0` (regularity) | 0 violations | 0 | ✅ |

## 3. Break-region map — where claims fail OUTSIDE the assumptions

Each assumption relaxed one at a time (fresh scans). This is the boundary of validity.

| Relaxation | What breaks | Frequency | Reading |
|---|---|---|---|
| **(A2)** `λW>c` | interior `x₁₁*>1` → corner `x*=1` | **19.8%** of tested band | (A2) is exactly the interior-vs-corner guard for the `(1,1)` state. In-region: 0%. |
| **`μ=λ` → `μ≠λ`** | Prop 2 `z*>x*` fails | **22.1%** | Prop 2's equal-rate hypothesis is binding; failures concentrate at `μ≫λ`. |
| **(A2′)** `μW>c` (A2 kept) | `(0,0)` interiority (`Δ₀>0`, `x₀₀*∈(0,1)`, `mx₀₀*<1`) | **0 violations** | (A2′) is load-bearing for the *analytical* root-selection lemma and the unconditional `dx₀₀*/dV₁₀`, but numerical `(0,0)` interiority is **robust** to its relaxation in the tested band. |

**Takeaways for the writeup.**
1. No parameter region *inside* the maintained assumptions breaks any numerical claim.
2. The assumptions are not merely convenient — (A2) and the `μ=λ` hypothesis each mark a real boundary
   beyond which a specific result fails. (A2), moreover, is *conservative*: the true `(1,1)` corner
   onset is at `λW/c ≈ 1.5–1.7`, not 1 (Task 1).
3. (A2′) buys the analytical `(0,0)` proofs; the numerical `(0,0)` equilibrium survives its relaxation,
   so the assumption is sufficient but not numerically tight for interiority.

---

## 4. Transmission package (Tasks 7–9, Step 8) — Section 5.5 of the draft

**Symbolic identities (`verify_scripts/verify_step8_transmission.py`, log `numerics/step8_results.txt`).**
All verified as exact symbolic zeros in sympy, independently of `verify_dx_dV11.py`:

| ID | Identity | Consequence |
|---|---|---|
| I1 | `G_L,x\|_{G_L=0} = −cλz²/(2x)` | leader BR slope |
| I2 | `G_F,z\|_{G_F=0} = −cxN/(2z(1−λz))` | laggard BR slope (`N` = provocation bracket) |
| R1 | `det J = c²·[core]` + positive-term decomposition | **`det J > 0` proven under (A1)** |
| R2 | Cramer numerator for `dx*/dV₁₁` | **`dx*/dV₁₁ > 0` proven whenever `x* < √2·z*`** |
| R4 | Cramer numerator for `dz*/dV₁₁` | exact sign boundary for `dz*/dβ₁₁` (Lem. leadersign) |

Numeric spot checks (seed 8, N=20,000 solved equilibria): 0 `det J ≤ 0`; 0 `num_x ≤ 0` with
`x* < √2·z*`; Cramer vs finite difference max rel-err `3.6e-8` (x), `6.7e-6` (z); 0 sign mismatches
for the R4 bracket.

**How far the `√2` condition reaches (`task7_transmission_scans.py`, log `task7_results.txt`).**
Seed 41, N=60,000 over the `μ≠λ` region: `x* < √2·z*` holds in **60,000/60,000 draws (100%)**,
including **all 13,176** draws where Prop 2 itself fails (`z* ≤ x*`). On the task2 grid config the
max ratio `x*/(√2·z*)` is **0.813**. Transmission is proven at every equilibrium the scans reach.
Figure: `figs/prop2_frontier_sqrt2_overlay.png`; grid data: `task7_sqrt2_grid.csv` (regenerable).

**Is `√2` tight? (`task7_margin_addendum.txt`).** Yes — exactly right, not slack. Adversarial corner
sweep pushes the ratio to **0.9894**; the limit probe (`λ→0`, `β₁₁→1`, `μ→½`) drives it toward 1
(0.99999 at the most extreme point) without attaining it. So the result is reported as a conditional
theorem plus region-wide numerical verification of the condition, not an unconditional claim.

**Leader sign boundary (Task 7 Part C).** At 40,000 scanned equilibria the sign of the FD total
`dz*/dβ₁₁` matches the sign of the analytical R4 bracket
`B = λz³(1−λz)/2 − μ(1−λz)x³ − λzx²` in **40,000/40,000** cases (0 near-zero flags).
Negative share 95.48%, consistent with the ~95.5% prior. Figure: `figs/dz_dbeta11_map_R4overlay.png`
(exact boundary `B=0` overlaid — a derived curve, not a fit).

**Figure updates (Tasks 8–9).**
- Task 8: two candidate replacements for `c_lambda_V11_sweep.png` that strip all continuous-time
  framing — `c_lambda_V11_sweep_norefline.png` (no benchmark line) and
  `c_lambda_V11_sweep_limitline.png` (benchmark relabeled as the exact rare-breakthrough limit
  `V₁₁ → (1+2β₁₁)W/(3(1+β₁₁))`, which is also the exact `c→∞` saturation value). The original paper
  figure is not overwritten; **the choice between variants is pending**.
- Task 9: `figs/best_response_clean.png` — decluttered best-response figure for Section 4.2
  (smooth BR curves, equilibrium point, 45° line showing `z* > x*`).

**Remaining numerical-only claims after the upgrades.** Convergence of best-response iteration
(stability) — converged in every draw, supports but does not prove stability; the region-wide
coverage of the `√2` condition; the `dz*/dβ₁₁` negative share (~95.5%, no general claim made);
the `dx*/dλ` reversal share (~2%); the provocation-region share (~99.9%).

---

## 5. Reproducibility

| Task | Script | Key output |
|---|---|---|
| 1 | `task1_11_closedform.py` | `figs/corner_regime.png`, `task1_corner_reference.csv` |
| 2 | `task2_asymmetric.py` | `figs/prop2_frontier.png`, `task2_frontier_grid.csv` |
| 3 | `task3_initial_state.py` | `figs/task3_00_diagnostics.png`, `task3_violations.csv` (empty) |
| 4 | `task4_comparative_statics.py` | `figs/dz_dbeta11_map.png`, `dx_dlam_reversals.png`, `laggard_BR.png` |
| 5 | `task5_paper_figures.py` | `figs/beta11_sweep.png`, `beta00_sweep.png`, `c_lambda_V11_sweep.png` |
| 6 | `task6_summary.py`, `task6_diagnose_signs.py`, `task6_confirmation.py` | `task6_confirmation.txt`, `task6_confirmation_table.csv`, this memo |
| 6 (7/11 re-run) | same scripts, post-edit re-confirmation | `task6_rerun_7-11.txt`, `task6_rerun_log.txt` |
| 7 | `task7_transmission_scans.py` | `task7_results.txt`, `task7_margin_addendum.txt`, `figs/prop2_frontier_sqrt2_overlay.png`, `figs/dz_dbeta11_map_R4overlay.png`, `task7_sqrt2_grid.csv` |
| 8 | `task8_c_lambda_regen.py` | `figs/c_lambda_V11_sweep_norefline.png`, `figs/c_lambda_V11_sweep_limitline.png` |
| 9 | `task9_br_clean_fig.py` | `figs/best_response_clean.png` |
| Step 8 | `verify_scripts/verify_step8_transmission.py`, `numerics/verify_dx_dV11.py` | `step8_results.txt` |

Seeds/sizes match `verify_scripts` for parity; residual tol `1e-9`, iteration tol `1e-13`,
sign re-checks at 50-digit mpmath precision. Task 7 uses seed 41 (N=60,000 / 40,000); the step8
spot checks use seed 8 (N=20,000); the margin addendum uses seed 43 (N=60,000).
