"""
Step 3c: corrected verification.
(i)  R1/R2 via exact polynomial division:  F_theta = quotient*F + remainder.
(ii) R8 priority result with the CORRECT prefactor:
     dV11/dbeta = c (T1 - sqrt(Dt) T2) / (4 lam^2 (1+b)^3 sqrt(Dt)),
     where q = lam^2 W / c,
       T1 = 18 + 3q(1+b)^2 + q(1-b^2) + q^2(1-b^2)(1+b)^2 + q^2(1-b^2)^2 - 8q(1+b),
       T2 = 6 - 2q(1+b),   Dt = 9 - 2q(1-b^2) + q^2(1-b^2)^2.
     Then: T2 > 4 > 0;  T1 > 10 > 0;  T1^2 - Dt T2^2 = 16q(1+b)[q^2(1-b^2)^2 - 2q(1+b^3) + 9] > 0
     for q in (0,1/2), b in [0,1)  =>  T1 > sqrt(Dt) T2  =>  dV11/dbeta > 0.
     (A1)-(A2) imply q = lam*(lam W/c) < 1/2. QED.
"""
import sympy as sp

x = sp.symbols('x', positive=True)
W, c, lam, b, q = sp.symbols('W c lambda beta q', positive=True)
M = lam*(1+b)
u = lam**2*W*(1-b**2)
A = c*M
B = -(3*c + u)
C = 2*lam*W*(1-b)
F = sp.expand(A*x**2 + B*x + C)
Delta = 9*c**2 - 2*c*u + u**2

print("="*70)
print("R1 (redone): poly division  F_lam = (1/lam) F + x(3c-u)/lam")
print("="*70)
F_lam = sp.expand(sp.diff(A, lam)*x**2 + sp.diff(B, lam)*x + sp.diff(C, lam))
quo, rem = sp.div(F_lam, F, x)
print("quotient:", sp.simplify(quo), "  [expect 1/lam]")
target1 = x*(3*c - u)/lam
print("remainder == x(3c-u)/lam:", "PASS" if sp.simplify(sp.expand(rem) - sp.expand(target1)) == 0 else "FAIL")
print("=> at x*: F(x*)=0 so F_lam(x*) = x*(3c-u)/lam > 0 since u < c/2 < 3c under (A1)-(A2).")
print("=> dx*/dlam = F_lam/sqrt(Delta) > 0. QED (Luke's 12c>W condition not needed)")

print("="*70)
print("R2 (redone): poly division  F_b = (1/(1+b)) F + [x(3c+M^2 W) - 4 lam W]/(1+b)")
print("="*70)
F_b = sp.expand(sp.diff(A, b)*x**2 + sp.diff(B, b)*x + sp.diff(C, b))
quo2, rem2 = sp.div(F_b, F, x)
print("quotient:", sp.simplify(quo2), "  [expect 1/(1+b)]")
target2 = (x*(3*c + M**2*W) - 4*lam*W)/(1+b)
print("remainder == [x(3c+M^2W)-4lamW]/(1+b):",
      "PASS" if sp.simplify(sp.expand(rem2) - sp.expand(target2)) == 0 else "FAIL")

print("="*70)
print("R8 (redone): exact prefactor identity for dV11/dbeta")
print("="*70)
S = sp.sqrt(Delta)
x_star = (3*c + u - S)/(2*c*M)
V11 = c*x_star/(2*lam*(1+b)) + b*W/(1+b)
dV_db = sp.diff(V11, b)

T1 = 18 + 3*q*(1+b)**2 + q*(1-b**2) + q**2*(1-b**2)*(1+b)**2 + q**2*(1-b**2)**2 - 8*q*(1+b)
T2 = 6 - 2*q*(1+b)
Dt = 9 - 2*q*(1-b**2) + q**2*(1-b**2)**2

# Substitute q -> lam^2 W / c  (dimensionless spillover-scaled prize)
qval = lam**2*W/c
cand = c*(T1 - sp.sqrt(Dt)*T2)/(4*lam**2*(1+b)**3*sp.sqrt(Dt))
cand = cand.subs(q, qval)
diff_check = sp.simplify(sp.radsimp(dV_db - cand))
print("dV11/dbeta == c(T1 - sqrt(Dt)T2)/(4 lam^2 (1+b)^3 sqrt(Dt)):",
      "PASS" if diff_check == 0 else f"FAIL: {diff_check}")

print("="*70)
print("R8 positivity chain (pure algebra in (q,b)):")
print("="*70)
# (a) T2 > 4 for q<1/2, b<1:  T2 = 6 - 2q(1+b) > 6 - 2*(1/2)*2 = 4
print("(a) T2 = 6 - 2q(1+b) > 4 > 0 for q<1/2, b<1.  [immediate]")
# (b) T1 > 10:  only negative term is -8q(1+b) > -8; all others >= 0; 18 - 8 = 10
print("(b) T1 > 18 - 8q(1+b) > 18 - 8 = 10 > 0.      [immediate]")
# (c) factorization identity
P = sp.expand(T1**2 - Dt*T2**2)
P_claim = sp.expand(16*q*(1+b)*(q**2*(1-b**2)**2 - 2*q*(1+b**3) + 9))
print("(c) T1^2 - Dt*T2^2 == 16q(1+b)[q^2(1-b^2)^2 - 2q(1+b^3) + 9]:",
      "PASS" if sp.simplify(P - P_claim) == 0 else "FAIL")
# (d) inner factor >= 9 - 2q(1+b^3) > 9 - 2*(1/2)*2 = 7 > 0
print("(d) inner factor >= 9 - 2q(1+b^3) > 9 - 2 = 7 > 0 for q<1/2, b<1.  [immediate]")
print("=> T1^2 > Dt T2^2 with T1,T2>0 => T1 > sqrt(Dt) T2 => dV11/dbeta > 0.  QED")

print("="*70)
print("Consistency: numeric spot checks of the identity and sign")
print("="*70)
import random
random.seed(7)
f_exact = sp.lambdify((W, c, lam, b), dV_db, 'math')
f_cand  = sp.lambdify((W, c, lam, b), cand, 'math')
worst = 0
for _ in range(2000):
    lv = random.uniform(1e-3, 0.4999); bv = random.uniform(1e-6, 0.9999)
    cv = random.uniform(1e-2, 10.0);  Wv = random.uniform(1e-3, cv/lv)
    e1, e2 = f_exact(Wv, cv, lv, bv), f_cand(Wv, cv, lv, bv)
    worst = max(worst, abs(e1 - e2)/max(abs(e1), 1e-12))
    assert e1 > 0
print(f"max relative discrepancy exact-vs-candidate over 2000 draws: {worst:.2e}")
print("all dV11/dbeta values positive: PASS")
