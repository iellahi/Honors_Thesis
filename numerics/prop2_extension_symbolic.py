"""
SYMBOLIC verification of the analytic proof that z* > x* for mu <= lam (Prop. 3).
Filename note: "prop2" is the pre-renumbering label for leader dominance,
which is Proposition 3 in the final thesis. The filenames are kept because
Appendix C.3 cites them by name.


Chain of lemmas to verify symbolically:
 (I)  Equilibrium identity (dag): G_L = 0 and G_F = 0 imply
        sec_LHS * V11 = sec_RHS * (W - V11),  i.e. RatioFn = (W-V11)/V11.
 (II) Key inequality (sec): Phi := sec_RHS - sec_LHS >= 0 for x>=z, 0<mu<=lam.
      Proof structure to verify:
        (a) Phi is AFFINE in mu  (d2 Phi/dmu2 = 0).
        (b) slope  dPhi/dmu = (1-lam z) x (x^2/2 - z^2(1 - lam z/2)).
        (c) Phi(mu=0) = lam z (x^2 - (1/2)(1-lam z) z^2)  > 0 for x>=z.
        (d) Phi(mu=lam) = lam * psi, with
              psi(z,z) = z^3 (lam z/2)(3 - lam z) > 0,
              dpsi/dx|_{x=z} > 0,  and  d2psi/dx2 = 3(1-lam z)x + 2z > 0.
      => min over mu in [0,lam] is an endpoint; both endpoints >=0 for x>=z => Phi>=0.

Thesis: Appendix C.3, step 2 of the leader-dominance proof.
"""
import sympy as sp

z, x, lam, mu, W, V = sp.symbols('z x lam mu W V', positive=True)

# Best-response quadratics (equal zero at equilibrium)
G_L = sp.Rational(1,2)*lam*(1-mu*x)*z**2 + mu*x*z - lam*mu*x*(W-V)     # /c
G_F = sp.Rational(1,2)*mu*(1-lam*z)*x**2 + lam*z*x - mu*lam*z*(1-lam*z)*V  # /c

# ---- (I) Verify the (dag) identity: eliminate to show RatioFn = (W-V)/V ----
# sec brackets
B_L = sp.Rational(1,2)*lam*(1-mu*x)*z + mu*x
B_R = sp.Rational(1,2)*mu*(1-lam*z)*x + lam*z
sec_LHS = z**2*(1-lam*z)*B_L
sec_RHS = x**2*B_R
# From G_F=0: solve for V ; from G_L=0: solve for (W-V). Then form RatioFn and (W-V)/V.
V_from_F  = sp.solve(G_F, V)[0]
WmV_from_L = sp.solve(G_L, W)[0] - V_from_F   # (W) - V, with V substituted
# RatioFn = sec_LHS/sec_RHS ; identity claims sec_LHS * V = sec_RHS * (W-V)
lhs_id = sp.simplify(sec_LHS * V_from_F)
rhs_id = sp.simplify(sec_RHS * WmV_from_L)
print("(I) dag identity  sec_LHS*V - sec_RHS*(W-V) simplifies to:",
      sp.simplify(lhs_id - rhs_id))

# ---- (II) Key inequality Phi = sec_RHS - sec_LHS ----
Phi = sp.expand(sec_RHS - sec_LHS)
print("(IIa) d^2 Phi/dmu^2 =", sp.simplify(sp.diff(Phi, mu, 2)), " (0 => affine in mu)")
slope = sp.simplify(sp.diff(Phi, mu))
slope_claim = (1-lam*z)*x*(x**2/2 - z**2*(1 - lam*z/2))
print("(IIb) dPhi/dmu - claimed slope =", sp.simplify(slope - slope_claim))

Phi0 = sp.simplify(Phi.subs(mu, 0))
Phi0_claim = lam*z*(x**2 - sp.Rational(1,2)*(1-lam*z)*z**2)
print("(IIc) Phi(mu=0) - claim =", sp.simplify(Phi0 - Phi0_claim))
# Phi0 > 0 for x>=z:  x^2 - (1/2)(1-lam z)z^2 >= z^2 - (1/2)z^2 = z^2/2 > 0
print("      Phi(mu=0) at x=z:", sp.simplify(Phi0.subs(x, z)), " (>0)")

PhiL = sp.simplify(Phi.subs(mu, lam))
psi = sp.simplify(PhiL/lam)
print("(IId) psi(z,z) =", sp.factor(psi.subs(x, z)))
dpsi = sp.diff(psi, x)
print("      dpsi/dx at x=z =", sp.factor(sp.simplify(dpsi.subs(x, z))))
print("      d2psi/dx2 =", sp.simplify(sp.diff(psi, x, 2)), " = 3(1-lam z)x+2z")
print("      d2psi/dx2 - (3(1-lam z)x+2z) =", sp.simplify(sp.diff(psi, x, 2) - (3*(1-lam*z)*x+2*z)))

# Positivity certificates (all should be manifestly positive for 0<lam z<1, 0<z, x>=z)
print()
print("Positivity summary (want all > 0 for x>=z, 0<lam*z<1):")
print("  Phi(mu=0)|x=z  =", sp.factor(Phi0.subs(x, z)))
print("  psi|x=z        =", sp.factor(psi.subs(x, z)))
print("  dpsi/dx|x=z    =", sp.factor(sp.simplify(dpsi.subs(x, z))))
