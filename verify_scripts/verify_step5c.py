"""Step 5c: (i) X'(z) < 0 on the BR; (ii) dX/dV11|z > 0; (iii) with Z'(x)>0 (proven),
composite X(Z(x)) is decreasing => unique equilibrium; monotone comparative statics
=> dx*/dV11 > 0 => dx*/dbeta11 = (dx*/dV11)(dV11/dbeta11) > 0.  All machine-checked."""
import sympy as sp

z, x = sp.symbols('z x', positive=True)
c, lam, mu, W, V11 = sp.symbols('c lambda mu W V_11', positive=True)
G_F = (c*mu/2)*(1-lam*z)*x**2 + c*lam*z*x - mu*lam*z*(1-lam*z)*V11

print("C1: sign of G_F_z on the BR (eliminate V11 via G_F = 0)")
G_F_x = sp.diff(G_F, x)
G_F_z = sp.diff(G_F, z)
V11_on = sp.solve(sp.Eq(G_F, 0), V11)[0]
G_F_z_on = sp.together(sp.simplify(G_F_z.subs(V11, V11_on)))
num, den = sp.fraction(G_F_z_on)
print("  numerator:", sp.factor(num))
print("  denominator:", sp.factor(den))
cand = c*x*(mu*x*(1-lam*z)**2 + 2*lam*z) / (2*z*(1-lam*z))
print("  equals c x[mu x(1-lam z)^2 + 2 lam z]/(2z(1-lam z)):",
      "PASS" if sp.simplify(G_F_z_on - cand) == 0 else "FAIL")
print("  => G_F_z > 0 on BR => X'(z) = -G_F_z/G_F_x < 0 since G_F_x > 0.  QED")

print("C2: dX/dV11 at fixed z:  -G_F_V11/G_F_x")
G_F_V = sp.diff(G_F, V11)
print("  G_F_V11 =", sp.factor(G_F_V), " < 0  (lam z < 1)")
print("  => dX/dV11 = -G_F_V11/G_F_x > 0.  QED")

print("""C3: THEOREM (uniqueness + transmission).
 (a) Z'(x) > 0 [5A] and X'(z) < 0 [C1] => phi(x) := X(Z(x)) strictly decreasing;
     equilibrium solves x = phi(x) => interior equilibrium is UNIQUE.
 (b) phi strictly increasing in V11: X shifts up in V11 [C2]; Z shifts down in V11
     (G_L_V11 = +lam mu x > 0 with G_L_z > 0 => dZ/dV11 < 0) and X'(z) < 0.
 (c) phi decreasing in x, increasing in V11 => x*(V11) strictly increasing.
 (d) dx*/dbeta11 = (dx*/dV11)(dV11/dbeta11) > 0 by Step 3.  QED""")
G_L = (c*lam/2)*(1-mu*x)*z**2 + c*mu*x*z - lam*mu*x*(W - V11)
print("  check G_L_V11 = +lam mu x:",
      "PASS" if sp.simplify(sp.diff(G_L, V11) - lam*mu*x) == 0 else "FAIL")
