# Numerical validation — innovation race with knowledge spillovers

This is the numerical side of my thesis (Section 5). I took every closed-form and
comparative-statics result in the paper and checked it against an independent solve of the model,
to make sure the algebra actually holds up.

**Short version: it all holds. 0 violations across 60,000 parameter draws.** A handful of edge
cases tripped the finite-difference check, but those turned out to be numerical artifacts. I re-ran
them at 50-digit precision and they're clean.

Want the deeper version, with every figure and every check explained?
See **[README_detailed.md](README_detailed.md)**. The full technical memo, with confirmation tables
and caveats, is **[`numerics/SUMMARY_memo.md`](numerics/SUMMARY_memo.md)**.

## The model

![Figure 1 — state transitions of the innovation race](figs/state_diagram_preview.png)

Two firms race through two stages: research, then development. Each period a firm breaks through
with some probability. That probability rises with its own effort, and also with the rival's effort
through the spillover `β` (knowledge leaks between firms). Both firms can break through in the same
period, which jumps the game straight from `(0,0)` to `(1,1)`. A tie at the finish is a coin flip.
Winner takes the prize `W`. Effort has quadratic cost.

## The main result

![Frontier effort vs competitive balance](figs/beta11_sweep.png)

This is the one that matters. As the development-phase spillover `β₁₁` goes up, frontier effort
`x₁₁*` drops. But the laggard's catch-up effort `x*` rises, and the value of being tied at the
front, `V₁₁`, climbs with it. That's the trade-off at the core of the paper: spillovers make the
leaders lazier but keep the race alive.

## How far the catch-up theorem reaches

![sqrt(2) frontier overlay](figs/prop2_frontier_sqrt2_overlay.png)

New since the last version. The catch-up result (`dx*/dV₁₁ > 0`) is now a proven theorem whenever
the laggard's effort stays below √2 times the leader's. So the question for the numerics changed:
how much of the parameter space satisfies that condition? Answer: **all of it**. 60,000 out of
60,000 draws, including every single draw where leader dominance itself fails. And the √2 constant
is not slack. Pushing the parameters to their most extreme corner drives the ratio to 0.99, but
never past 1.

## Leader dominance, and where it breaks (Prop 2)

![Prop 2 frontier](figs/prop2_frontier.png)

The leader out-invests the laggard (`z* > x*`) whenever the two firms move at the same rate.
That's proven, and 0 failures in 30,000 draws confirm it. It only breaks once the laggard's rate
`μ` gets well above the leader's `λ` (about 22% of the `μ≠λ` cases). So the equal-rate assumption
is doing real work, not just tidying up.

## The laggard's best response is non-monotone

![Laggard best response](figs/laggard_BR.png)

Higher leader effort can actually *provoke* the laggard into working harder rather than
discouraging it, up to a threshold. The equilibrium lands in that provocation region ~99.9% of
the time. There's also a cleaner version of this picture for the paper:
`figs/best_response_clean.png`.

## What else I checked (all clean)

- the `(1,1)` closed form matches an independent Bellman+FOC solve to ~1e-14
- `(0,0)` state: `Δ₀ > 0` and `x₀₀*` stays in `(0,1)` everywhere
- laggard catch-up `dx*/dβ₁₁ > 0`: 100% of draws
- `det J > 0`: now proven analytically, and 0 violations numerically
- the exact sign boundary for the leader's response `dz*/dβ₁₁` matches the numerics in
  40,000 out of 40,000 draws (`figs/dz_dbeta11_map_R4overlay.png`)
- `β₀₀` and `c`/`λ` effects: see `figs/beta00_sweep.png` and `figs/c_lambda_V11_sweep.png`
  (two candidate updates of the `c`/`λ` figure are included; the final pick is pending)

## If you want more

- the full walkthrough of every figure and check → [README_detailed.md](README_detailed.md)
- confirmation tables, break-region map, and caveats → [`numerics/SUMMARY_memo.md`](numerics/SUMMARY_memo.md)
- to re-run everything yourself → `bash run_all.sh` (Python 3.10+,
  `pip install -r requirements.txt`, a few minutes)

Full thesis PDF coming once the writeup's done.
