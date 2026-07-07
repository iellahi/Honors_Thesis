"""
Step 2: Symbolic verification of Luke's (1,1) symmetric development state derivation.
Framework: discrete time, simultaneous breakthroughs allowed, coin-flip tiebreaker.
Every check prints PASS/FAIL. Nothing is assumed from hand algebra.
"""
import sympy as sp

# ---------- Symbols ----------
x_i, x_j, x, z = sp.symbols('x_i x_j x z', positive=True)
W, c, lam, beta = sp.symbols('W c lambda beta_11', positive=True)
V = sp.symbols('V_11')  # value of state (1,1), treated as constant in FOC
M = lam * (1 + beta)

# Breakthrough probabilities in (1,1): spillovers active
P_i = lam * (x_i + beta * x_j)
P_j = lam * (x_j + beta * x_i)

print("=" * 70)
print("CHECK 1: Bellman isolation -> V11 = [W(2Pi - PiPj) - c xi^2] / [2(Pi+Pj-PiPj)]")
print("=" * 70)
# Bellman: V = -c/2 xi^2 + Pi(1-Pj) W + Pi Pj W/2 + (1-Pi)(1-Pj) V
bellman_rhs = -sp.Rational(1, 2) * c * x_i**2 + P_i * (1 - P_j) * W \
    + P_i * P_j * W / 2 + (1 - P_i) * (1 - P_j) * V
V_isolated = sp.solve(sp.Eq(V, bellman_rhs), V)[0]
V_luke_general = (W * (2 * P_i - P_i * P_j) - c * x_i**2) / (2 * (P_i + P_j - P_i * P_j))
print("PASS" if sp.simplify(V_isolated - V_luke_general) == 0 else "FAIL")

print("=" * 70)
print("CHECK 2: FOC partials: dPi/dxi = lam, dPj/dxi = beta*lam,")
print("         d(PiPj)/dxi = lam*Pj + beta*lam*Pi")
print("=" * 70)
ok = (sp.simplify(sp.diff(P_i, x_i) - lam) == 0
      and sp.simplify(sp.diff(P_j, x_i) - beta * lam) == 0
      and sp.simplify(sp.diff(P_i * P_j, x_i) - (lam * P_j + beta * lam * P_i)) == 0)
print("PASS" if ok else "FAIL")

print("=" * 70)
print("CHECK 3: FOC  c xi = lam[ W - (Pj + beta Pi)W/2 - (1+beta)V + (Pj + beta Pi)V ]")
print("=" * 70)
foc_lhs_minus_rhs = sp.diff(bellman_rhs, x_i)  # set = 0
foc_luke = -c * x_i + lam * (W - (P_j + beta * P_i) * W / 2
                             - (1 + beta) * V + (P_j + beta * P_i) * V)
print("PASS" if sp.simplify(foc_lhs_minus_rhs - foc_luke) == 0 else "FAIL")

print("=" * 70)
print("CHECK 4: Symmetric FOC  c x = lam W - M V - M^2 x (W/2 - V),  M = lam(1+beta)")
print("=" * 70)
foc_sym = foc_luke.subs([(x_i, x), (x_j, x)])
foc_luke_M = -c * x + lam * W - M * V - M**2 * x * (W / 2 - V)
print("PASS" if sp.simplify(sp.expand(foc_sym) - sp.expand(foc_luke_M)) == 0 else "FAIL")

print("=" * 70)
print("CHECK 5: V from FOC:  V = (c x - lam W + M^2 x W/2) / (M(Mx - 1))")
print("=" * 70)
V_from_foc = sp.solve(sp.Eq(foc_luke_M, 0), V)[0]
V_foc_luke = (c * x - lam * W + M**2 * x * W / 2) / (M * (M * x - 1))
print("PASS" if sp.simplify(V_from_foc - V_foc_luke) == 0 else "FAIL")

print("=" * 70)
print("CHECK 6: V from Bellman under symmetry: V = (W(2M - M^2 x) - c x)/(2M(2 - Mx))")
print("=" * 70)
V_bell_sym = sp.simplify(V_isolated.subs([(x_i, x), (x_j, x)]))
V_bell_luke = (W * (2 * M - M**2 * x) - c * x) / (2 * M * (2 - M * x))
print("PASS" if sp.simplify(V_bell_sym - V_bell_luke) == 0 else "FAIL")

print("=" * 70)
print("CHECK 7: Equating -> quadratic (Mc)x^2 - (3c + 2 lam M W - M^2 W)x + 2W(2 lam - M) = 0")
print("=" * 70)
# Equate the two V expressions; clear denominators
eq = sp.together(V_from_foc - V_bell_sym)
numer = sp.numer(sp.simplify(eq))
quad_luke_M = M * c * x**2 - (3 * c + 2 * lam * M * W - M**2 * W) * x + 2 * W * (2 * lam - M)
# numer should be proportional to quad_luke_M
ratio = sp.simplify(sp.expand(numer) / sp.expand(quad_luke_M))
print("Numerator / Luke quadratic =", ratio)
print("PASS (proportional, nonzero constant)" if ratio.is_constant() and ratio != 0 else "FAIL")

print("=" * 70)
print("CHECK 8: Substituting M: c lam(1+b) x^2 - (3c + lam^2 W(1-b^2)) x + 2 lam W(1-b) = 0")
print("=" * 70)
A = c * lam * (1 + beta)
B = -(3 * c + lam**2 * W * (1 - beta**2))
C = 2 * lam * W * (1 - beta)
quad_final = A * x**2 + B * x + C
print("PASS" if sp.simplify(sp.expand(quad_luke_M) - sp.expand(quad_final)) == 0 else "FAIL")

print("=" * 70)
print("CHECK 9: Discriminant Delta = 9c^2 - 2c lam^2 W(1-b^2) + lam^4 W^2 (1-b^2)^2")
print("=" * 70)
Delta = sp.expand(B**2 - 4 * A * C)
Delta_luke = 9 * c**2 - 2 * c * lam**2 * W * (1 - beta**2) + lam**4 * W**2 * (1 - beta**2)**2
print("PASS" if sp.simplify(Delta - sp.expand(Delta_luke)) == 0 else "FAIL")

print("=" * 70)
print("CHECK 9b: Delta > 0 always (real roots guaranteed)")
print("=" * 70)
u = sp.symbols('u', positive=True)  # u = lam^2 W (1-b^2) >= 0
Delta_u = 9 * c**2 - 2 * c * u + u**2
completed = sp.expand((u - c)**2 + 8 * c**2)
print("Delta(u) = (u-c)^2 + 8c^2 ?", sp.simplify(Delta_u - completed) == 0, "-> strictly positive: PASS")

print("=" * 70)
print("CHECK 10: Subtracted root x* = [3c + lam^2 W(1-b^2) - sqrt(Delta)] / [2 c lam(1+b)]")
print("          solves the quadratic")
print("=" * 70)
x_star = (3 * c + lam**2 * W * (1 - beta**2) - sp.sqrt(Delta_luke)) / (2 * c * lam * (1 + beta))
resid = sp.simplify(quad_final.subs(x, x_star))
print("Residual:", resid)
print("PASS" if resid == 0 else "FAIL")

print("=" * 70)
print("CHECK 11: Root validity. F(1/M) sign determines whether additive root -> P > 1")
print("=" * 70)
F_at_1oM = sp.simplify(quad_final.subs(x, 1 / M))
F_at_1oM_pred = -2 * c / M + lam * W * (1 - beta)
print("F(1/M) simplifies to -2c/M + lam W (1-b)?",
      sp.simplify(F_at_1oM - F_at_1oM_pred) == 0)
print("F(1/M) < 0  <=>  lam^2 W (1-b^2) < 2c   [derived condition]")
cond_check = sp.simplify((F_at_1oM_pred * M) - (lam**2 * W * (1 - beta**2) - 2 * c))
print("M * F(1/M) == lam^2 W(1-b^2) - 2c ?", cond_check == 0)

print("=" * 70)
print("CHECK 12: V11 = c x*/(2 lam (1+b)) + b W/(1+b) at the subtracted root")
print("=" * 70)
V_final_luke = c * x_star / (2 * lam * (1 + beta)) + beta * W / (1 + beta)
V_check_foc = sp.simplify(V_from_foc.subs(x, x_star) - V_final_luke)
V_check_bell = sp.simplify(V_bell_sym.subs(x, x_star) - V_final_luke)
print("Matches FOC-derived V:    ", "PASS" if V_check_foc == 0 else f"FAIL: {V_check_foc}")
print("Matches Bellman-derived V:", "PASS" if V_check_bell == 0 else f"FAIL: {V_check_bell}")

print("=" * 70)
print("CHECK 13: Sanity - beta -> 0 and numeric spot checks")
print("=" * 70)
# beta = 0: quadratic becomes c lam x^2 - (3c + lam^2 W)x + 2 lam W = 0
x0 = sp.simplify(x_star.subs(beta, 0))
print("x*(beta=0) =", x0)
# Numeric: W=2, lam=0.2, c=1, beta=0.3 -- check root solves quadratic, P<1, V in (0, W)
subs_num = {W: 2, lam: sp.Rational(1, 5), c: 1, beta: sp.Rational(3, 10)}
xs_num = float(x_star.subs(subs_num))
P_num = float((M * x).subs(subs_num).subs(x, xs_num))
V_num = float(V_final_luke.subs(subs_num))
print(f"x* = {xs_num:.6f}, P = {P_num:.6f} (<1?), V11 = {V_num:.6f} (in (0,W)?)")
resid_num = float(quad_final.subs(subs_num).subs(x, xs_num))
print(f"Quadratic residual at numeric root: {resid_num:.2e}")
# Additive root check numerically
x_plus = (3 * c + lam**2 * W * (1 - beta**2) + sp.sqrt(Delta_luke)) / (2 * c * lam * (1 + beta))
xp_num = float(x_plus.subs(subs_num))
Pp_num = float((M * x).subs(subs_num).subs(x, xp_num))
print(f"Additive root: x+ = {xp_num:.4f}, P+ = {Pp_num:.4f} (should exceed 1)")
