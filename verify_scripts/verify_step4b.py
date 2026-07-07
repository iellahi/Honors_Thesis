"""
Step 4b: (i) master symmetric-state quadratic + nesting of (1,1); (ii) root selection
for (0,0) under (A2') mu*W <= c, conditional on the value ordering 0<=V01<=V11<=V10<=W.
"""
import sympy as sp

x = sp.symbols('x', positive=True)
c, mu, b0, lam, b1, W = sp.symbols('c mu beta_00 lambda beta_11 W', positive=True)
V10, V01, V11 = sp.symbols('V_10 V_01 V_11')

# ---------- Master quadratic (claimed cleaned form) ----------
m = mu*(1+b0)
G = V10 + b0*V01 - (1+b0)*V11
D = V10 - V01
F00 = c*m*x**2 - (3*c + 2*mu*m*G)*x + 2*mu*(1-b0)*D

print("="*70)
print("N1: cleaned form matches raw S5 polynomial (up to factor (1+b0))")
print("="*70)
# Rebuild raw numerator exactly as in verify_step4.py
x_i, x_j = sp.symbols('x_i x_j', positive=True)
V00 = sp.symbols('V_00')
P_i = mu*(x_i + b0*x_j); P_j = mu*(x_j + b0*x_i)
rhs = -sp.Rational(1,2)*c*x_i**2 + P_i*(1-P_j)*V10 + P_j*(1-P_i)*V01 \
      + P_i*P_j*V11 + (1-P_i)*(1-P_j)*V00
V00_iso = sp.solve(sp.Eq(V00, rhs), V00)[0].subs([(x_i,x),(x_j,x)])
foc = sp.diff(rhs, x_i).subs([(x_i,x),(x_j,x)])
V00_foc = sp.solve(sp.Eq(foc, 0), V00)[0]
numer = sp.expand(sp.numer(sp.together(sp.simplify(V00_foc - V00_iso))))
ratio = sp.simplify(numer / sp.expand((1+b0)*F00))
print("raw numerator / [(1+b0) * cleaned F00] =", ratio,
      "-> PASS" if ratio.is_constant() and ratio != 0 else "FAIL")

print("="*70)
print("N2: NESTING -- F00 with (V10,V01,V11)->(W,0,W/2), mu->lam, b0->b1")
print("    equals the verified (1,1) quadratic")
print("="*70)
F11_claim = F00.subs([(V10, W), (V01, 0), (V11, W/2), (mu, lam), (b0, b1)])
F11 = c*lam*(1+b1)*x**2 - (3*c + lam**2*W*(1-b1**2))*x + 2*lam*W*(1-b1)
print("PASS" if sp.simplify(sp.expand(F11_claim) - sp.expand(F11)) == 0 else "FAIL")

print("="*70)
print("N3: F00(1/m) reduction:  m*F00(1/m) = -2c - 2 mu m [G - (1-b0)D]")
print("    and  G - (1-b0)D = b0(V10 - V11) - (V11 - V01)")
print("="*70)
lhs = sp.expand(m*F00.subs(x, 1/m))
rhs3 = sp.expand(-2*c - 2*mu*m*(G - (1-b0)*D))
print("Reduction identity:", "PASS" if sp.simplify(lhs - rhs3) == 0 else "FAIL")
bracket = sp.expand(G - (1-b0)*D)
bracket_claim = sp.expand(b0*(V10 - V11) - (V11 - V01))
print("Bracket identity:  ", "PASS" if sp.simplify(bracket - bracket_claim) == 0 else "FAIL")
print("""
Conditional Lemma (root selection in (0,0)). Suppose (A1), (A2') mu*W <= c, and the
value ordering 0 <= V01 <= V11 <= V10 <= W [proof obligation: Step 5]. Then:
  bracket = b0(V10-V11) - (V11-V01) <= b0(V10-V11) <= b0*W <= W,
  so  m*F00(1/m) = -2c - 2 mu m*bracket <= -2c + 2 mu m W <= -2c + 2 m c < 0  [mu W<=c, m<1]
  => F00(1/m) < 0 => roots straddle 1/m => additive root has P=mx>1 (infeasible);
     the subtracted root is the unique valid candidate.""")

print("="*70)
print("N4: constant/leading terms: A0 = cm > 0; C0 = 2 mu (1-b0) D > 0 iff V10 > V01")
print("    => both roots same sign; sum = (3c + 2 mu m G)/(cm) -- positive if 3c+2mu m G>0")
print("    Note G = (V10-V11) - b0(V11-V01): sign ambiguous in general; but even if G<0,")
print("    |2 mu m G| <= 2 mu m W <= 2mc < 3c under (A2') => middle coeff stays negative.")
print("="*70)
# Machine-check that |2 mu m G| < 3c under ordering + (A2'): G in [-(1+b0)W/2 ... W] crude;
# use |G| <= max(V10-V11, b0(V11-V01)) <= W  (each term in [0,W], difference of two nonneg <= W)
# 2 mu m |G| <= 2 mu m W <= 2 m c < 3c. Symbolic witness of the chain endpoints:
chain = sp.simplify(3*c - 2*mu*m*W - (3*c - 2*m*c))  # equals 2m(c - mu W) >= 0 under A2'
print("3c - 2 mu m W - (3c - 2mc) = 2m(c - mu W):",
      "PASS" if sp.simplify(chain - 2*m*(c - mu*W)) == 0 else "FAIL")
print("=> sum of roots > 0, product > 0: both roots strictly positive (given V10>V01).")

print("="*70)
print("N5: Discriminant Delta0 = (3c + 2 mu m G)^2 - 8 c m mu (1-b0) D")
print("    Positivity NOT automatic in general values -- deferred to post-Step-5 closure.")
print("    Numeric probe with plausible orderings to see if it can fail in-region:")
print("="*70)
import random
random.seed(3)
Delta0 = (3*c + 2*mu*m*G)**2 - 8*c*m*mu*(1-b0)*D
f = sp.lambdify((c, mu, b0, V10, V01, V11, W), Delta0, 'math')
fails = 0; n = 200000; mn = 1e18
for _ in range(n):
    muv = random.uniform(1e-3, 0.4999); b0v = random.uniform(0, 0.9999)
    cv = random.uniform(1e-2, 10.0);  Wv = random.uniform(1e-3, cv/muv)  # A2'
    # random ordered values 0 <= V01 <= V11 <= V10 <= W with V11 <= W/2 (feasibility)
    v01 = random.uniform(0, Wv/2); v11 = random.uniform(v01, Wv/2); v10 = random.uniform(v11, Wv)
    val = f(cv, muv, b0v, v10, v01, v11, Wv)
    if val <= 0: fails += 1
    mn = min(mn, val)
print(f"scan n={n}: Delta0 <= 0 occurrences: {fails}, min = {mn:.6e}")
