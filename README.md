# Numerical Validation — Innovation Race with Knowledge Spillovers

This repository contains the numerical validation and figure generation for the theory in the
thesis. It confirms the closed-form equilibria and comparative-static results of a **discrete-time,
two-firm, two-stage innovation race** with knowledge spillovers, solved as a Markov Perfect
Equilibrium (MPE). Every closed form is cross-checked against an *independent* numerical solve of the
underlying Bellman + first-order-condition system, so the figures rest on the model primitives rather
than on the algebra being validated.

> **Read this file top-to-bottom to see every result — no code needs to be run.** The scripts are
> provided for full reproducibility (`bash run_all.sh` regenerates everything).

**Full writeup:** the complete thesis PDF (model, proofs, discussion) will be added here after the
defense. This repo is the reproducible numerical companion to Section 5.

---

## The model in one picture

![State transitions of the innovation race](figs/state_diagram_preview.png)

Two firms race through two stages. Each period a firm makes a breakthrough with probability `Pᵢ`
(increasing in its effort and in the rival's effort via the spillover `β`). Simultaneous research
breakthroughs jump the game directly from `(0,0)` to `(1,1)`; in the asymmetric states the leader's
development breakthrough takes precedence; a simultaneous breakthrough in `(1,1)` is settled by a fair
coin flip (expected prize `W/2`). Costs are quadratic (`c/2·x²`), the prize is winner-take-all `W`, no
discounting.

**Maintained assumptions.** (A1) `x∈[0,1]`, `λ,μ<½`; (A2) `λW≤c`; (A2′) `μW≤c`. All scans draw
`W∈(0, c/max(λ,μ)]` so (A2) and (A2′) both hold.

---

## Results at a glance

A consolidated scan of **60,000** parameter draws re-checks every theorem and every numerical claim.
Finite-difference sign flags are re-verified at **50-digit precision** (`mpmath`).

**Within the maintained assumptions, no numerical claim breaks and every theorem passes — 0 hard
violations across all checks.** (Independently reproduced under a second random seed.)

| | Result | Outcome |
|---|---|---|
| **Closed forms** | `(1,1)` root vs independent Bellman+FOC solve | max diff **3.5×10⁻¹⁴** |
| **Theorems** | Δ>0, root selection, value ordering, Prop 1/2, master quadratic, BR slope | **0 / 60000** violations |
| **Prop 2** | leader out-invests laggard `z*>x*` at `μ=λ` | **0 / 30000** failures |
| **Catch-up** | `dx*/dβ₁₁>0` (higher spillover → more laggard effort) | **100%** of draws |
| **(0,0) state** | `Δ₀>0`, `x₀₀*∈(0,1)`, `mx₀₀*<1` | **0** violations, `min Δ₀≈0.02` |
| **Regularity** | `det J>0` at the equilibrium | **0** violations |

---

## Key figures

### 1. The headline trade-off: frontier effort vs competitive balance
![beta11 sweep](figs/beta11_sweep.png)

As the development-phase spillover `β₁₁` rises, symmetric **frontier effort `x₁₁*` falls** while the
**laggard's catch-up effort `x*` rises** and the **tied-state value `V₁₁` increases** toward `W/2`.
This is the central mechanism: spillovers trade frontier intensity for competitive balance, and the
laggard-catch-up channel operates through the anticipated value of the tied state.

### 2. Proposition 2 and its frontier
![Prop 2 frontier](figs/prop2_frontier.png)

Leader dominance `z*>x*` holds on and around `μ=λ` (Prop 2, proven). It fails only when the laggard's
research rate `μ` sufficiently exceeds the leader's development rate `λ` — never for `μ≤λ`, rising to
81% of draws for `μ/λ≥3` (22% overall). The equal-rate hypothesis is *binding, not cosmetic*.

### 3. The non-monotone laggard best response
![Laggard BR](figs/laggard_BR.png)

The laggard's best response `X(z)` is non-monotone: `sign X′(z) = sign(μx(1−λz)² − 2λ²z²)`. Below the
threshold, higher leader effort **provokes** the laggard (erosion of `V₀₁` dominates); above it,
classic **discouragement**. The equilibrium sits on the provocation arm in **99.9%** of the region.

### 4. Initial-stage spillovers
![beta00 sweep](figs/beta00_sweep.png)

Higher initial-phase spillover `β₀₀` lowers own initial effort `x₀₀*` and the breakthrough probability,
yet **raises** the initial-state value `V₀₀` — the spillover benefit outweighs the effort reduction.

### 5. A qualitatively new dependence: `V₁₁` on `c` and `λ`
![V11 vs c and lambda](figs/c_lambda_V11_sweep.png)

In discrete time `V₁₁` is **increasing in `c`** and **decreasing in `λ`** — a dependence that is
identically flat in the continuous-time model. The magnitude is modest (a few percent; the figure uses
`β₁₁=0.1` for visibility and annotates the exact change). The contribution is the *existence and sign*
of the dependence, not its size.

<details>
<summary><b>Supporting / diagnostic figures</b> (click to expand)</summary>

- **`figs/corner_regime.png`** — outside (A2), interior `x₁₁*` exceeds 1 and the equilibrium is the
  corner `x*=1`. (A2) is *conservative*: the true corner onset is at `λW/c≈1.5–1.7`, not 1.
- **`figs/task3_00_diagnostics.png`** — `Δ₀>0`, interior `x₀₀*`, and `mx₀₀*<1−β₀₀` everywhere (0/60k).
- **`figs/dz_dbeta11_map.png`** — laggard response to `β₁₁` is positive everywhere; the leader's is
  sign-ambiguous (negative in ~96% of draws).
- **`figs/dx_dlam_reversals.png`** — the indirect `V₁₁` channel flips the sign of `dx*/dλ` in a thin
  high-`λ`/low-`μ` band (~2% of draws).
</details>

---

## Claim → script → figure map

| Script | Validates | Output |
|---|---|---|
| `numerics/task1_11_closedform.py` | `(1,1)` closed form vs independent solve; corner regime | `figs/corner_regime.png`, `task1_corner_reference.csv` |
| `numerics/task2_asymmetric.py` | asymmetric fixed point; Prop 2; `μ≠λ` frontier | `figs/prop2_frontier.png` |
| `numerics/task3_initial_state.py` | `(0,0)`: `Δ₀>0`, interiority, unconditional `dx₀₀*/dV₁₀` | `figs/task3_00_diagnostics.png`, `task3_violations.csv` (empty) |
| `numerics/task4_comparative_statics.py` | `dx*/dβ₁₁`, `dz*/dβ₁₁`, `dx*/dλ` reversals, provocation, `det J` | `figs/dz_dbeta11_map.png`, `dx_dlam_reversals.png`, `laggard_BR.png` |
| `numerics/task5_paper_figures.py` | `β₁₁`, `β₀₀`, `c`, `λ` sweeps | `figs/beta11_sweep.png`, `beta00_sweep.png`, `c_lambda_V11_sweep.png` |
| `numerics/task6_summary.py` | consolidated confirmation + break-region map (seed 202) | `task6_confirmation.txt` |
| `numerics/task6_confirmation.py` | independent consolidated re-scan (seed 11) | `task6_confirmation_table.csv` |
| `numerics/task6_diagnose_signs.py` | 50-digit mpmath re-check of FD sign flags | console |

The full write-up of the confirmation, the break-region map, and honest caveats is in
**[`numerics/SUMMARY_memo.md`](numerics/SUMMARY_memo.md)**. The `verify_scripts/` folder holds the
original **SymPy** scripts that verify the *analytical* derivations (Steps 2–6) symbolically.

---

## Confirmation detail

**Theorems (0 violations, cross-checked over 60k draws).** Δ=(u−c)²+8c²>0; subtracted-root selection;
key inequality √Δ>3c−u; value ordering `0≤V₀₁≤V₁₁<W/2≤V₁₀≤W`; Prop 1 (`x₁₁*` signs); Prop 2; laggard
BR slope identity; additive-root best responses solve their quadratics (residual ≤9.5×10⁻¹⁵); the
master-quadratic identity `F₀₀[terminal subs] − F₁₁ ≡ 0` (verified symbolically).

**Sign-flag caveat, resolved.** Naïve finite differences flag 221 cases for Prop 1′ (`V₁₁` signs) and 18
for Prop 3 (`x₀₀*` signs). All are near-zero-derivative artifacts — `dV₁₁/dc` near its `c→∞`
saturation, and `dx₀₀*/dV₁₁ ∝ μ²` at tiny `μ`. At 50-digit precision **every flagged case recovers the
correct sign: 0 real violations.**

**Break-region map (validity boundary).** Relaxing each assumption one at a time shows they mark real
boundaries, not mere convenience:

| Relaxation | What breaks | Frequency |
|---|---|---|
| (A2) `λW>c` | interior `x₁₁*>1` → corner `x*=1` | 19.8% of tested band (in-region: 0%) |
| `μ=λ → μ≠λ` | Prop 2 `z*>x*` fails | 22.1% |
| (A2′) `μW>c` | `(0,0)` interiority | 0 violations (analytical proofs need it; numerics robust) |

**Not proven (numerical only).** Three results are numerically supported but not proven — the
non-monotone laggard BR blocks an analytical sign: full-equilibrium `dx*/dβ₁₁>0`, `dx*/dV₁₁>0`, and
`det J>0` (uniqueness/stability of the interior MPE).

---

## Reproduce

Requires Python 3.10+ and the packages in `requirements.txt` (`numpy`, `matplotlib`, `sympy`,
`mpmath` — **no SciPy needed**; the independent solvers are hand-rolled bisection / value iteration).

```bash
pip install -r requirements.txt
bash run_all.sh          # regenerates every figure and confirmation log (~40s total)
```

Each script is self-contained and writes to `figs/` and `numerics/`. Scans reuse fixed seeds for
parity (residual tolerance `1e-9`, iteration tolerance `1e-13`).

---

## Repository layout

```
numerical-validation/
├── README.md                 this guide
├── requirements.txt
├── run_all.sh
├── figs/                     all figures (PNG, 300 dpi)
├── numerics/                 validation scripts, logs, small data, SUMMARY_memo.md
└── verify_scripts/           original SymPy scripts (analytical verification, Steps 2–6)
```
