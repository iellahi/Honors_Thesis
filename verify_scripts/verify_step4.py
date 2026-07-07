"""
Step 4: original derivation of the (0,0) state, new framework.
Transitions: i-only -> V10, j-only -> V01, BOTH -> V11 (direct jump), neither -> stay.
Spillovers active: P_i = mu(x_i + b0 x_j).  m = mu(1+b0).
V10, V01, V11 treated as constants (backward induction).
Deliverables: (S1) Bellman isolation, (S2) FOC, (S3) symmetric FOC in m-form,
(S4) equilibrium quadratic, (S5) discriminant, (S6) sanity limits.
"""
import sympy as sp

x_i, x_j, x = sp.symbols('x_i x_j x', positive=True)
c, mu, b0 = sp.symbols('c mu beta_00', positive=True)
V10, V01, V11, V00 = sp.symbols('V_10 V_01 V_11 V_00')
m = mu*(1 + b0)

P_i = mu*(x_i + b0*x_j)
P_j = mu*(x_j + b0*x_i)

print("="*70)
print("S1: Bellman isolation (general, asymmetric efforts)")
print("="*70)
rhs = -sp.Rational(1,2)*c*x_i**2 + P_i*(1-P_j)*V10 + P_j*(1-P_i)*V01 \
      + P_i*P_j*V11 + (1-P_i)*(1-P_j)*V00
V00_isolated = sp.solve(sp.Eq(V00, rhs), V00)[0]
V00_target = (P_i*(1-P_j)*V10 + P_j*(1-P_i)*V01 + P_i*P_j*V11
              - sp.Rational(1,2)*c*x_i**2) / (P_i + P_j - P_i*P_j)
print("V00 = [Pi(1-Pj)V10 + Pj(1-Pi)V01 + PiPj V11 - (c/2)xi^2]/(Pi+Pj-PiPj):",
      "PASS" if sp.simplify(V00_isolated - V00_target) == 0 else "FAIL")

print("="*70)
print("S2: FOC (differentiate Bellman RHS in x_i, continuation V00 held fixed)")
print("="*70)
foc = sp.diff(rhs, x_i)   # = 0 at optimum; foc = -c x_i + mu*[...]
# Candidate compact form:
foc_target = -c*x_i + mu*( (1-P_j)*V10 - b0*P_i*V10
                           + b0*(1-P_i)*V01 - P_j*V01
                           + (P_j + b0*P_i)*V11
                           - (1-P_j)*V00 - b0*(1-P_i)*V00 )
print("c x_i = mu[ (1-Pj-b0 Pi)V10 + (b0(1-Pi)-Pj)V01 + (Pj+b0 Pi)V11 - ((1-Pj)+b0(1-Pi))V00 ]:",
      "PASS" if sp.simplify(foc - foc_target) == 0 else "FAIL")

print("="*70)
print("S3: Symmetric FOC.  x_i=x_j=x, P = m x.")
print("Claim: c x = mu(1-P)(V10 + b0 V01) - mu P(b0 V10 + V01) + m P V11 - m(1-P)V00")
print("="*70)
foc_sym = sp.simplify(foc.subs([(x_i, x), (x_j, x)]))
P = m*x
foc_sym_target = -c*x + mu*(1-P)*(V10 + b0*V01) - mu*P*(b0*V10 + V01) \
                 + m*P*V11 - m*(1-P)*V00
print("PASS" if sp.simplify(foc_sym - foc_sym_target) == 0 else "FAIL")

print("="*70)
print("S4: Symmetric Bellman value")
print("Claim: V00 = [P(1-P)(V10+V01) + P^2 V11 - (c/2)x^2] / [P(2-P)]")
print("="*70)
V00_sym = sp.simplify(V00_isolated.subs([(x_i, x), (x_j, x)]))
V00_sym_target = (P*(1-P)*(V10+V01) + P**2*V11 - sp.Rational(1,2)*c*x**2) / (P*(2-P))
print("PASS" if sp.simplify(V00_sym - V00_sym_target) == 0 else "FAIL")

print("="*70)
print("S5: Equate FOC-implied V00 with Bellman-implied V00 -> equilibrium quadratic in x")
print("="*70)
V00_from_foc = sp.solve(sp.Eq(foc_sym_target, 0), V00)[0]
eq = sp.together(sp.simplify(V00_from_foc - V00_sym_target))
numer = sp.expand(sp.numer(eq))
# collect as polynomial in x
poly = sp.Poly(numer, x)
coeffs = poly.all_coeffs()
print("Degree in x:", poly.degree())
print("Raw coefficients (highest degree first):")
for k, co in enumerate(coeffs):
    print(f"  x^{poly.degree()-k}:", sp.factor(sp.simplify(co)))
