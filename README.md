# Numerical validation — knowledge spillovers in a two-stage innovation race

This is the numerical side of my honors thesis (Section 5). Every closed form and
every comparative-statics sign in the paper is checked here against an independent
solve of the model.

## What "checked" means

The independent solver never touches the closed forms. It rebuilds the equilibrium
from the Bellman equations and first-order conditions directly — bisection in the
symmetric states, best-response iteration in the asymmetric ones. If the closed
forms and the solver agree, the algebra in the paper is right.

The checks come in three layers:

- **SymPy.** Every algebraic identity in the appendices is required to simplify to
  exact zero. Those claims rest on no tolerance.
- **Finite differences.** Comparative-statics signs are checked at solved equilibria.
- **mpmath at 50 digits.** A finite difference near zero cannot resolve a sign, so
  every flagged case is re-solved at high precision and decided there.

**Result: nothing breaks.** Closed forms match the independent solver to about
1e-14. Every symbolic identity is exactly zero. Every proven sign holds. No
numerical claim fails anywhere in the maintained parameter region, across three
independent seeds. The finite-difference flags that did appear were all resolution
artifacts and cleared at 50 digits.

## The model

![State transitions of the innovation race](figs/state_diagram_preview.png)

Two firms race through two stages: research, then development. Each period a firm
breaks through with a probability that rises with its own effort and, through the
spillover `β`, with its rival's. Spillovers only work between firms in the same
stage — a firm still in research cannot absorb its rival's late-stage advances.
Both firms can break through in the same period, which sends the game straight from
`(0,0)` to `(1,1)`. A tie at the finish is a coin flip. Winner takes `W`. Effort
costs `c/2 · x²`.

## The trade-off

![Frontier effort against competitive balance](figs/beta11_sweep.png)

This is the result the paper is built on. As the development-stage spillover `β₁₁`
rises, frontier effort `x₁₁*` falls — each firm free-rides on what spills over from
the other. But the value of being tied at the frontier, `V₁₁`, rises, and so does
the laggard's catch-up effort `x*`. Spillovers buy competitive balance and pay for
it in frontier effort.

The laggard gets nothing from the leader directly. Spillovers are switched off when
the race is asymmetric. It works harder only because the state it is racing toward
is worth more.

## How far the catch-up theorem reaches

![Margin to the sqrt(2) frontier](figs/prop2_frontier_sqrt2_overlay.png)

The catch-up result `dx*/dV₁₁ > 0` is proven whenever the laggard's effort stays
below √2 times the leader's. That condition is weaker than leader dominance, so the
numerical question is how much of the parameter space satisfies it.

All of it. 60,000 out of 60,000 draws, including all 13,176 draws where leader
dominance itself fails. Pushing the parameters into their most extreme corner drives
the ratio to 0.99 without crossing 1, so the condition comes close to binding.
Whether the result survives past √2 is untested.

## Leader dominance, and where it breaks (Proposition 3)

![Leader-dominance frontier](figs/prop2_frontier.png)

The leader out-invests the laggard, `z* > x*`, whenever the laggard's research rate
does not exceed the leader's development rate (`μ ≤ λ`). That is proven, and 0
failures across the `μ ≤ λ` draws confirm it.

It breaks once `μ` climbs above `λ`: about 22% of `μ ≠ λ` draws fail, and **every
failure has `μ > λ`**. The failure share climbs with the ratio — 0% at the boundary,
81% once `μ/λ ≥ 3`. The hypothesis is one-sided, and the numerics map exactly where
it bites.

## The laggard's best response is not monotone

![Laggard best response](figs/laggard_BR.png)

Up to a threshold, a harder-working leader pushes the laggard to work harder rather
than discouraging it. Past the threshold the laggard retreats. Equilibria sit on the
first arm 99.8–99.9% of the time, so discouragement in this model shows up as a
lower level of effort, not as a collapsing response.

## Everything else

- `(1,1)` closed form against an independent Bellman + FOC solve: agrees to ~1e-14
- `(0,0)` state: `Δ₀ > 0` and `x₀₀* ∈ (0,1)` everywhere, minimum `Δ₀ ≈ 0.021`
- laggard catch-up `dx*/dβ₁₁ > 0`: 100% of draws
- `det J > 0`: proven analytically, and 0 violations numerically
- the exact sign boundary for the leader's response `dz*/dβ₁₁` matches the numerics
  at all 40,000 scanned equilibria (`figs/dz_dbeta11_map_R4overlay.png`)
- initial-state regularity (`figs/task3_00_diagnostics.png`) and the `β₀₀` sweep
  (`figs/beta00_sweep.png`)

## Reproduce

```bash
pip install -r requirements.txt   # Python 3.10+
bash run_all.sh                   # ~90 seconds
```

Every scan uses a fixed seed and reports its sample size. Figures land in `figs/`,
logs and tables in `numerics/`.

## More

- every figure and check, walked through, with the confirmation tables and caveats
  → [README_detailed.md](README_detailed.md)
- the thesis itself → `Ellahi_2026_Innovation_Race_Spillovers.pdf`

## Acknowledgement

I thank my advisor, Professor Luke Boosey, for his guidance on the model and the
proofs throughout this project.
