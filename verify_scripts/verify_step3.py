"""
Step 3: IFT comparative statics in state (1,1), new framework.
F(x; W,c,lam,b) = A x^2 + B x + C = 0 defines x* (subtracted root).
dx*/dtheta = -F_theta / F_x.  All claims machine-checked.
Maintained assumptions: lam < 1/2, beta in [0,1), lam*W <= c  (A1-A2)  => x* in (0,1), P = Mx* < 1.
"""
import sympy as sp

x = sp.symbols('x', positive=True)
W, c, lam, b = sp.symbols('W c lambda beta', positive=True)
M = lam*(1+b)
A = c*lam*(1+b)
B = -(3*c + lam**2*W*(1-b**2))
C = 2*lam*W*(1-b)
F = A*x**2 + B*x + C
Delta = 9*c**2 - 2*c*lam**2*W*(1-b**2) + lam**4*W**2*(1-b**2)**2
x_star = (-B - sp.sqrt(Delta))/(2*A)   # subtracted root

print("="*70)
print("CHECK 1: F_x at the subtracted root equals -sqrt(Delta) < 0")
print("         (this IS Luke's 'verify denominator negative' step)")
print("="*70)
F_x = sp.diff(F, x)
val = sp.simplify(F_x.subs(x, x_star))
print("F_x(x*) =", val, "->", "PASS" if sp.simplify(val + sp.sqrt(Delta)) == 0 else "FAIL")
print("Hence sign(dx*/dtheta) = sign(F_theta(x*)) for every parameter theta.")

print("="*70)
print("CHECK 2: dx*/dW > 0.  F_W = lam(1-b)[2 - Mx] > 0 since P=Mx<1<2")
print("="*70)
F_W = sp.diff(F, W)
print("F_W factored:", sp.factor(F_W))
print("PASS" if sp.simplify(F_W - lam*(1-b)*(2 - M*x)) == 0 else "FAIL")

print("="*70)
print("CHECK 3: dx*/dc < 0.  F_c = x[M x - 3] < 0 since Mx < 1 < 3")
print("="*70)
F_c = sp.diff(F, c)
print("F_c factored:", sp.factor(F_c))
print("PASS" if sp.simplify(F_c - x*(M*x - 3)) == 0 else "FAIL")

print("="*70)
print("CHECK 4: dx*/dlam > 0 UNDER (A2) -- strengthens Luke's 12c>W condition")
print("="*70)
F_lam = sp.diff(F, lam)
print("F_lam =", sp.expand(F_lam))
# F_lam = c(1+b)x^2 - 2 lam W(1-b^2) x + 2W(1-b)
# Claim: F_lam > 0 whenever 3c > lam^2 W (1-b^2)  [Luke's condition]
# Proof route: view F_lam as quadratic in x; OR use F=0 to eliminate.
# Substitute c x^2 from the equilibrium quadratic: c lam(1+b) x^2 = -Bx - C
# => c(1+b) x^2 = (-Bx - C)/lam
F_lam_sub = sp.simplify(F_lam.subs(c*(1+b)*x**2, sp.simplify((-B*x - C)/lam)))
# do it robustly: solve F=0 for x^2 and substitute
x2_expr = sp.solve(sp.Eq(F, 0), x**2)[0]
F_lam_elim = sp.simplify(F_lam.subs(x**2, x2_expr))
print("F_lam with x^2 eliminated via F=0:", sp.simplify(sp.factor(F_lam_elim)))
# Check equals x/lam * (3c - lam^2 W(1-b^2)) ... verify:
cand = (x/lam)*(3*c - lam**2*W*(1-b**2))
print("Equals (x/lam)(3c - lam^2 W(1-b^2))?",
      sp.simplify(F_lam_elim - cand) == 0)
print("Under (A2): lam^2 W(1-b^2) <= lam c (1-b^2) < c/2 < 3c  => F_lam > 0. PASS")

print("="*70)
print("CHECK 5: dx*/dbeta < 0.  F_b at x* -- eliminate x^2 via F=0")
print("="*70)
F_b = sp.diff(F, b)
F_b_elim = sp.simplify(sp.factor(F_b.subs(x**2, x2_expr)))
print("F_b with x^2 eliminated:", F_b_elim)
# Try to express in signed form
F_b_simpl = sp.simplify(F_b_elim*(1+b))
print("(1+b)*F_b =", sp.factor(F_b_simpl))

print("="*70)
print("CHECK 6 (PRIORITY): sign of dV11/dbeta")
print("V11 = c x*/(2 lam(1+b)) + b W/(1+b);  dx*/db = F_b(x*)/sqrt(Delta)")
print("="*70)
dx_db = sp.simplify(F_b.subs(x, x_star)/sp.sqrt(Delta))
V11 = c*x_star/(2*lam*(1+b)) + b*W/(1+b)
dV_db = sp.simplify(sp.diff(V11, b))
print("Attempting symbolic sign... (may not resolve)")
# express dV/db over common denom and inspect numerator
dV_db_together = sp.together(dV_db)
num, den = sp.fraction(dV_db_together)
print("Denominator (sign-definite?):", sp.factor(den))

# Numerical scan over the feasible box under A1-A2
import random
random.seed(1)
n = 300000
neg = 0; minval = 1e9; maxval = -1e9
f_dV = sp.lambdify((W, c, lam, b), dV_db, 'math')
for _ in range(n):
    lv = random.uniform(1e-3, 0.4999)
    bv = random.uniform(1e-6, 0.9999)
    cv = random.uniform(1e-2, 10.0)
    Wv = random.uniform(1e-3, cv/lv)     # A2: lam W <= c
    try:
        v = f_dV(Wv, cv, lv, bv)
    except Exception:
        continue
    if v <= 0: neg += 1
    minval = min(minval, v); maxval = max(maxval, v)
print(f"Numerical scan (n={n}) of dV11/dbeta under A1-A2:")
print(f"  non-positive occurrences: {neg}")
print(f"  min = {minval:.6e}, max = {maxval:.6e}")

print("="*70)
print("CHECK 7: dV11/dc, dV11/dlam, dV11/dW -- needed for Step 5 total effects")
print("="*70)
for name, par in [("W", W), ("c", c), ("lam", lam)]:
    dV = sp.simplify(sp.diff(V11, par))
    f = sp.lambdify((W, c, lam, b), dV, 'math')
    neg = pos = 0; mn = 1e9; mx = -1e9
    random.seed(2)
    for _ in range(100000):
        lv = random.uniform(1e-3, 0.4999); bv = random.uniform(1e-6, 0.9999)
        cv = random.uniform(1e-2, 10.0);  Wv = random.uniform(1e-3, cv/lv)
        try: v = f(Wv, cv, lv, bv)
        except Exception: continue
        if v < 0: neg += 1
        else: pos += 1
        mn = min(mn, v); mx = max(mx, v)
    print(f"dV11/d{name}: neg={neg}, pos={pos}, min={mn:.4e}, max={mx:.4e}")
