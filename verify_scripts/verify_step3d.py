"""
Appendix B.4: the anticipatory catch-up channel, dV11/dbeta11 > 0 (Proposition 2).

The appendix writes the derivative as

    dV11/dbeta11 = c (T1 - sqrt(Dt) T2) / (4 lam^2 (1+b)^3 sqrt(Dt)),

with, for q = lam^2 W / c < 1/2 and b = beta_11 < 1,

    T1 = 18 + 3q(1+b)^2 + q(1-b^2) + q^2(1-b^2)(1+b)^2 + q^2(1-b^2)^2 - 8q(1+b)
    T2 = 6 - 2q(1+b)
    Dt = 9 - 2q(1-b^2) + q^2(1-b^2)^2

and proves the sign by showing T1 > 0, T2 > 0, and T1^2 - Dt*T2^2 > 0, the last
via the factorization checked below. Everything here is polynomial in (q, b), so
each identity is asserted as an exact symbolic zero rather than to a tolerance.

Check (0) closes the loop back to the model: it rebuilds x11* and V11 from the
primitives (c, lam, W, beta), differentiates in beta, and confirms the result is
the T1/T2 expression above. Without it the remaining checks would verify the
algebra of the bound while taking the bounded expression on faith.

The bounds use only (A1) and (A2), which give q < 1/2.

Thesis: Appendix B.4 (Proposition 2, the tied-state value).
"""
import sympy as sp

q, b = sp.symbols('q beta', positive=True)

T1 = (18 + 3*q*(1+b)**2 + q*(1-b**2) + q**2*(1-b**2)*(1+b)**2
      + q**2*(1-b**2)**2 - 8*q*(1+b))
T2 = 6 - 2*q*(1+b)
Dt = 9 - 2*q*(1-b**2) + q**2*(1-b**2)**2

def check(label, expr_is_zero):
    print(f"{label}: {'PASS (symbolic zero)' if sp.simplify(expr_is_zero) == 0 else 'FAIL'}")

print("=" * 70)
print("B.4: dV11/dbeta11 > 0  (Proposition 2)")
print("=" * 70)

# (0) Closure: the T1/T2 expression really is dV11/dbeta11.
#     Rebuild the equilibrium objects from the model primitives,
#         u    = lam^2 W (1-b^2),      Delta = 9c^2 - 2cu + u^2,
#         x11* = (3c + u - sqrt(Delta)) / (2 c lam (1+b)),
#         V11  = c x11* / (2 lam (1+b)) + b W / (1+b),
#     differentiate V11 in b at fixed (c, lam, W) -- note q = lam^2 W / c does
#     not depend on b -- and subtract the appendix expression.
cM, lamM, WM = sp.symbols('c_model lambda_model W_model', positive=True)
uM   = lamM**2*WM*(1 - b**2)
DeM  = 9*cM**2 - 2*cM*uM + uM**2
x11M = (3*cM + uM - sp.sqrt(DeM))/(2*cM*lamM*(1 + b))
V11M = cM*x11M/(2*lamM*(1 + b)) + b*WM/(1 + b)

appendix_form = (cM*(T1 - sp.sqrt(Dt)*T2)/(4*lamM**2*(1 + b)**3*sp.sqrt(Dt))
                 ).subs(q, lamM**2*WM/cM)
check("(0) dV11/db = c(T1 - sqrt(Dt) T2)/(4 lam^2 (1+b)^3 sqrt(Dt))  [closure]",
      sp.radsimp(sp.together(sp.diff(V11M, b) - appendix_form)))

# (1) The factorization that carries the sign.
target = 16*q*(1+b)*(q**2*(1-b**2)**2 - 2*q*(1+b**3) + 9)
check("(1) T1^2 - Dt*T2^2 = 16q(1+b)[q^2(1-b^2)^2 - 2q(1+b^3) + 9]",
      sp.expand(T1**2 - Dt*T2**2) - sp.expand(target))

# (2) Dt is the discriminant in the substituted variables: Dt = Delta/c^2 with
#     u = q c (1-b^2).  Delta = 9c^2 - 2cu + u^2.
c_, u_ = sp.symbols('c u', positive=True)
Delta = 9*c_**2 - 2*c_*u_ + u_**2
check("(2) Dt = Delta/c^2 under u = q c (1-b^2)",
      sp.expand(Delta.subs(u_, q*c_*(1-b**2))/c_**2) - sp.expand(Dt))

print("-" * 70)
print("Sign bounds under (A1)-(A2), which give q < 1/2 and b < 1:")

# (3) T2 > 0: T2 = 6 - 2q(1+b) > 6 - 2*(1/2)*2 = 4.
worst_T2 = T2.subs({q: sp.Rational(1, 2), b: 1})
print(f"  T2 > {worst_T2} > 0  (worst case q -> 1/2, b -> 1)"
      f"   {'PASS' if worst_T2 > 0 else 'FAIL'}")

# (4) T1 > 0: every term is nonnegative except -8q(1+b) > -8.
print(f"  T1 > 18 - 8q(1+b) > 18 - 8 = 10 > 0"
      f"   {'PASS' if sp.simplify((18 - 8*q*(1+b)).subs({q: sp.Rational(1,2), b: 1})) > 0 else 'FAIL'}")

# (5) The bracket in (1) exceeds 9 - 2q(1+b^3) > 9 - 2 = 7 > 0.
bracket_low = (9 - 2*q*(1+b**3)).subs({q: sp.Rational(1, 2), b: 1})
print(f"  bracket > {bracket_low} > 0  (dropping the q^2 term, worst case)"
      f"   {'PASS' if bracket_low > 0 else 'FAIL'}")

print("-" * 70)
print("T1 > 0, T2 > 0 and T1^2 > Dt*T2^2 give T1 > sqrt(Dt) T2, hence")
print("dV11/dbeta11 > 0 under (A1)-(A2).  QED")
