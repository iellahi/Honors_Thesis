# Numerical validation — innovation race with knowledge spillovers

This is the numerical side of the thesis (Section 5). I took every closed-form and comparative-statics
result and checked it against an independent solve of the model, to make sure the algebra actually
holds up.

**Short version: it all holds. 0 violations across 60,000 parameter draws.** A handful of edge cases
tripped the finite-difference check, but those turned out to be numerical artifacts; I re-ran them at
50-digit precision and they're clean.

## The model

![Figure 1 — state transitions of the innovation race](figs/state_diagram_preview.png)

Two firms race through two stages. Each period a firm breaks through with probability `Pᵢ` (rising in
its own effort and the rival's, through the spillover `β`). Simultaneous research breakthroughs jump
straight from `(0,0)` to `(1,1)`; in the asymmetric states the leader wins ties; a tie in `(1,1)` is a
coin flip. Winner takes `W`, quadratic costs, no discounting.

## The main result

![Frontier effort vs competitive balance](figs/beta11_sweep.png)

This is the one that matters. As the development-phase spillover `β₁₁` goes up, frontier effort `x₁₁*`
drops, but the laggard's catch-up effort `x*` rises and the tied-state value `V₁₁` climbs with it.
That's the frontier-effort vs competitive-balance trade-off — the core mechanism of the paper.

## Leader dominance, and where it breaks (Prop 2)

![Prop 2 frontier](figs/prop2_frontier.png)

The leader out-invests the laggard (`z*>x*`) whenever the two move at the same rate — proven, and 0
failures in 30,000 draws. It only breaks once the laggard's rate `μ` gets well above the leader's `λ`
(~22% of the `μ≠λ` cases). So the equal-rate assumption is doing real work, not just tidying up.

## The laggard's best response is non-monotone

![Laggard best response](figs/laggard_BR.png)

Higher leader effort can actually *provoke* the laggard into working harder rather than discouraging
it, up to a threshold. The equilibrium lands in that provocation region ~99.9% of the time.

## What else I checked (all clean)

- the `(1,1)` closed form matches an independent Bellman+FOC solve to ~1e-14
- `(0,0)` state: `Δ₀>0` and `x₀₀*` stays in `(0,1)` everywhere
- laggard catch-up `dx*/dβ₁₁>0`: 100% of draws
- `det J>0` (regularity): 0 violations
- `β₀₀` and `c`/`λ` effects: see `figs/beta00_sweep.png` and `figs/c_lambda_V11_sweep.png`

## If you want more

- the full confirmation tables, break-region map, and caveats → [`numerics/SUMMARY_memo.md`](numerics/SUMMARY_memo.md)
- to re-run it all yourself → `bash run_all.sh` (Python 3.10+, `pip install -r requirements.txt`, ~40s)

Full thesis PDF coming once the writeup's done.
