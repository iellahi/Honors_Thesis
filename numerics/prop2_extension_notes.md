# Leader dominance for μ ≤ λ — proof record

This is the working record for what is Proposition 3 in the thesis. The result was
originally proved only at `μ = λ`; this note extends it to all `μ ≤ λ`.

**Status: PROVED (analytically, with symbolic + numerical verification).**
Leader dominance `z* > x*` now holds for all `μ ≤ λ`, not just `μ = λ`.

## Why the old proof did not directly generalize
The existing App C.3 proof proves *pointwise* best-response dominance `Z(y) > X(y)`
for every common rival effort `y` at `μ = λ` (the two BRs share a functional form there).
For `μ < λ` this pointwise claim is **false near y = 0**: the small-`y` expansions give
`Z(y) ≈ √(2μy(W−V₁₁)/c)` and `X(y) ≈ √(2λy·V₁₁/c)`, so `Z(y) > X(y)` there requires
`V₁₁ < μW/(λ+μ)`, which is strictly stronger than the available bound `V₁₁ < W/2` once
`μ < λ`. So a new (equilibrium-level, not pointwise-for-all-y) argument is needed.
Steps (ii) Z′>0 and (iii) contradiction of the old proof are reused unchanged; only the
pointwise step is replaced.

## The new argument (equilibrium-level)
Let `B_L = (λ/2)(1−μx)z + μx`, `B_R = (μ/2)(1−λz)x + λz`, and define
- `L(z,x) = z²(1−λz)·B_L`  (call it sec_LHS)
- `R(z,x) = x²·B_R`         (call it sec_RHS)
Both are strictly positive on the feasible region (`λz<1`, `μx<1`).

**(dag) Equilibrium identity.** Eliminating W and V₁₁ from `G_L = G_F = 0` gives exactly
`L·V₁₁ = R·(W − V₁₁)`.  Verified symbolically: direct substitution → 0, and reduction of
`L·V − R·(W−V)` modulo the ideal ⟨G_L, G_F⟩ → 0 (`prop2_dag_identity.py`).

**(§) Key inequality.** For all `x ≥ z > 0`, `0 < μ ≤ λ`, `λz < 1`:  `Φ := R − L ≥ 0`.
Proof: `Φ` is **affine in μ** (`∂²Φ/∂μ² = 0`, verified symbolically), so on `μ ∈ [0,λ]`
it is the convex combination `Φ(μ) = (1−μ/λ)Φ(0) + (μ/λ)Φ(λ)`; it suffices that both
endpoints are ≥ 0 for `x ≥ z`:
- `Φ(0) = λz(x² − ½(1−λz)z²) ≥ λz(z² − ½z²) = ½λz³ > 0`.
- `Φ(λ) = λ·ψ(z,x)` where `ψ` is increasing in x on `x ≥ z`
  (`∂²ψ/∂x² = 3(1−λz)x + 2z > 0`, `∂ψ/∂x|_{x=z} = ½z²(5−λ²z²) > 0`)
  and `ψ(z,z) = ½λz⁴(3−λz) > 0`.  Hence `ψ(z,x) > 0` for `x ≥ z`, so `Φ(λ) > 0`.
All endpoint/derivative identities verified symbolically (`prop2_extension_symbolic.py`).
`μ ≤ λ` is used *only here* (it keeps `μ` inside `[0,λ]`; for `μ > λ`, `Φ` can turn
negative — this is the one-sided necessity).

**Contradiction.** Suppose `x* ≥ z*`. By (§), `R(z*,x*) ≥ L(z*,x*) > 0`. Substituting into
(dag): `L·V₁₁ = R·(W−V₁₁) ≥ L·(W−V₁₁)`, and cancelling `L > 0` gives `V₁₁ ≥ W − V₁₁`,
i.e. `V₁₁ ≥ W/2` — contradicting Lemma 3 (`V₁₁ < W/2`). Hence `z* > x*`.  ∎

## Numerical verification
- `z* > x*`: 0 failures in 60,000 μ≤λ equilibrium draws (seed 11). (`prop2_extension_check.py`)
- (§) `R − L ≥ 0`: 0 violations in 3,000,000 free draws of `(λ, μ≤λ, z, x≥z)`;
  min margin `1.3e−11` (attained only at the x=z boundary, where it is exactly 0).
- (dag) identity: max relative residual `5.0e−10` at solved equilibria (solver tolerance).

## Where this lands in the thesis
- Proposition 3 (leader dominance): hypothesis `μ ≤ λ`, proof in Appendix C.3.
- The necessity remark is one-sided — only `μ > λ` can break the ordering.
- Corollary 1, `dx*/dβ₁₁ > 0`, extends to `μ ≤ λ`: it needs only `z* > x* ⇒ x* < √2 z*`.
- Section 5.4 maps the failure region, which lies entirely in the `μ > λ` half.
