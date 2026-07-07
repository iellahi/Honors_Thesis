# Numerical Validation — Summary Memo

**Model.** Discrete-time, two-firm, two-stage innovation race; simultaneous breakthroughs, coin-flip
tiebreaker, direct `(0,0)→(1,1)` jump. Maintained assumptions: **(A1)** `x∈[0,1]`, `λ,μ<½`;
**(A2)** `λW≤c`; **(A2′)** `μW≤c`.
**Scope.** Tasks 1–6, reusing the `verify_scripts/` solvers (`verify_step5b.py` reference).
Scripts + logs in `/numerics`, figures in `/figs`. No `.tex` files modified.

**Bottom line.** Within the maintained assumptions, **every theorem passes and every numerical claim
holds with zero violations.** The two apparent sign "failures" in the first pass were finite-difference
artifacts, confirmed clean by 50-digit mpmath re-checks (0 real violations). Outside the assumptions,
each claim fails exactly where the theory predicts it can. An **independent consolidated re-scan**
(`task6_confirmation.py`, seed 11, N=60,000) reproduces every result with **0 hard violations**
(min Δ₀=0.0213; provocation 99.92%; `dx*/dβ₁₁` 100% positive), corroborating the seed-202 tables below
under a different seed.

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

## 4. Reproducibility

| Task | Script | Key output |
|---|---|---|
| 1 | `task1_11_closedform.py` | `figs/corner_regime.png`, `task1_corner_reference.csv` |
| 2 | `task2_asymmetric.py` | `figs/prop2_frontier.png`, `task2_frontier_grid.csv` |
| 3 | `task3_initial_state.py` | `figs/task3_00_diagnostics.png`, `task3_violations.csv` (empty) |
| 4 | `task4_comparative_statics.py` | `figs/dz_dbeta11_map.png`, `dx_dlam_reversals.png`, `laggard_BR.png` |
| 5 | `task5_paper_figures.py` | `figs/beta11_sweep.png`, `beta00_sweep.png`, `c_lambda_V11_sweep.png` |
| 6 | `task6_summary.py`, `task6_diagnose_signs.py`, `task6_confirmation.py` | `task6_confirmation.txt`, `task6_confirmation_table.csv`, this memo |

Seeds/sizes match `verify_scripts` for parity; residual tol `1e-9`, iteration tol `1e-13`,
sign re-checks at 50-digit mpmath precision.
