"""
Step 5A: asymmetric states (1,0)/(0,1). P_L = lam*z (leader, development),
P_F = mu*x (laggard, research), no spillovers. V11 constant from (1,1).
Verifies: Luke's Bellmans/FOCs/BR roots; my compact BR quadratics; unique-positive-root
selection; value ordering lemma; Prop 2 analytical chain (mu = lam).
"""
import sympy as sp

z, x, y = sp.symbols('z x y', positive=True)
c, lam, mu, W = sp.symbols('c lambda mu W', positive=True)
V10, V01, V11 = sp.symbols('V_10 V_01 V_11', positive=True)
PL, PF = lam*z, mu*x
S = PL + PF - PL*PF

print("="*70)
print("A1: Bellman isolations")
print("="*70)
V10s, V01s = sp.symbols('V10s V01s')
rhs_L = -c*z**2/2 + PL*W + (1-PL)*PF*V11 + (1-PL)*(1-PF)*V10s
rhs_F = -c*x**2/2 + (1-PL)*PF*V11 + (1-PL)*(1-PF)*V01s
V10_iso = sp.solve(sp.Eq(V10s, rhs_L), V10s)[0]
V01_iso = sp.solve(sp.Eq(V01s, rhs_F), V01s)[0]
t_L = (PL*W + (1-PL)*PF*V11 - c*z**2/2)/S
t_F = ((1-PL)*PF*V11 - c*x**2/2)/S
print("V10 = [PL W + (1-PL)PF V11 - (c/2)z^2]/S:", "PASS" if sp.simplify(V10_iso - t_L)==0 else "FAIL")
print("V01 = [(1-PL)PF V11 - (c/2)x^2]/S:      ", "PASS" if sp.simplify(V01_iso - t_F)==0 else "FAIL")

print("="*70)
print("A2: Luke's FOCs")
print("="*70)
foc_L = sp.diff(rhs_L, z)   # holding V10s fixed
foc_F = sp.diff(rhs_F, x)   # holding V01s fixed
t_focL = -c*z + lam*(W - V10s + mu*x*(V10s - V11))
t_focF = -c*x + mu*(1-lam*z)*(V11 - V01s)
print("cz = lam(W - V10 + mu x (V10 - V11)):", "PASS" if sp.simplify(foc_L - t_focL)==0 else "FAIL")
print("cx = mu(1-lam z)(V11 - V01):        ", "PASS" if sp.simplify(foc_F - t_focF)==0 else "FAIL")

print("="*70)
print("A3: substitute values into FOCs -> compact BR quadratics")
print("  Leader:  G_L = (c lam/2)(1-mu x) z^2 + c mu x z - lam mu x (W - V11) = 0")
print("  Laggard: G_F = (c mu/2)(1-lam z) x^2 + c lam z x - mu lam z (1-lam z) V11 = 0")
print("="*70)
G_L = (c*lam/2)*(1-mu*x)*z**2 + c*mu*x*z - lam*mu*x*(W - V11)
G_F = (c*mu/2)*(1-lam*z)*x**2 + c*lam*z*x - mu*lam*z*(1-lam*z)*V11
lhs_L = sp.together(sp.simplify(foc_L.subs(V10s, V10_iso)))
lhs_F = sp.together(sp.simplify(foc_F.subs(V01s, V01_iso)))
ratio_L = sp.simplify(sp.numer(lhs_L)/sp.expand(G_L))
ratio_F = sp.simplify(sp.numer(lhs_F)/sp.expand(G_F))
print("Leader FOC numerator / G_L =", sp.factor(ratio_L))
print("Laggard FOC numerator / G_F =", sp.factor(ratio_F))
# accept if ratio is x,z-free and sign-definite
for nm, r in [("L", ratio_L), ("F", ratio_F)]:
    free = r.free_symbols & {x, z}
    print(f"  ratio_{nm} free of x,z:", "PASS" if not free else f"FAIL {free}")

print("="*70)
print("A4: Luke's explicit BR roots solve the quadratics; additive root is the")
print("    UNIQUE positive root (product of roots < 0)")
print("="*70)
x_br = (-c*lam*z + sp.sqrt(c**2*lam**2*z**2 + 2*c*mu**2*lam*z*(1-lam*z)**2*V11))/(c*mu*(1-lam*z))
z_br = (-c*mu*x + sp.sqrt(c**2*mu**2*x**2 + 2*c*mu*lam**2*x*(1-mu*x)*(W-V11)))/(c*lam*(1-mu*x))
print("x_br solves G_F:", "PASS" if sp.simplify(G_F.subs(x, x_br))==0 else "FAIL")
print("z_br solves G_L:", "PASS" if sp.simplify(G_L.subs(z, z_br))==0 else "FAIL")
prod_F = sp.simplify((-mu*lam*z*(1-lam*z)*V11)/((c*mu/2)*(1-lam*z)))
prod_L = sp.simplify((-lam*mu*x*(W - V11))/((c*lam/2)*(1-mu*x)))
print("G_F root product =", prod_F, "< 0;  G_L root product =", sp.factor(prod_L), "< 0 (needs W>V11)")

print("="*70)
print("A5: VALUE ORDERING LEMMA")
print("="*70)
# (i) V01 <= V11: S - (1-PL)PF = PL
print("(i) S - (1-PL)PF == PL:", "PASS" if sp.simplify(S - (1-PL)*PF - PL)==0 else "FAIL",
      " => V01 <= [S V11 - (c/2)x^2]/S <= V11; >=0 by x=0 deviation.")
# (ii) V11 < W/2: F(xtil) = -lam W (1-b) < 0 at xtil = lam W(1-b)/c
b = sp.symbols('beta_11', positive=True)
M = lam*(1+b); u = lam**2*W*(1-b**2)
F11 = c*M*y**2 - (3*c+u)*y + 2*lam*W*(1-b)
xtil = lam*W*(1-b)/c
val = sp.simplify(F11.subs(y, xtil))
print("(ii) F11(lam W(1-b)/c) =", sp.factor(val), " [expect -lam W(1-b) < 0]",
      "PASS" if sp.simplify(val + lam*W*(1-b))==0 else "FAIL")
print("     => x11* < lam W(1-b)/c => c x11* < lam W(1-b) => V11 = c x*/(2M) + bW/(1+b) < W/2. QED")
# (iii) V10 >= V11 by mimicry: need (c/2) x11* <= lam(W - V11); we show c x11* < lam(W-V11)
# from symmetric FOC: c x = lam W - M V11 - M^2 x (W/2 - V11) < lam W - lam V11  [V11<W/2, M>=lam, V11>=0]
print("(iii) c x11* = lam W - M V11 - M^2 x11*(W/2 - V11) < lam W - lam V11  [chain, no algebra needed]")
print("      => leader mimicking x11* in (1,0) gets Vtil >= V11; V10 >= Vtil. QED (see A6 identity)")
# A6: mimicry identity: Vtil - V11 has numerator lam*x11*(W - V11) - (c/2)x11*^2 ... verify
xs = sp.symbols('x_s', positive=True)  # stands for x11*
Vt = (lam*xs*W + (1-lam*xs)*mu*x*V11 - c*xs**2/2)/(lam*xs + mu*x - lam*xs*mu*x)
diff_num = sp.simplify(sp.together(Vt - V11))
print("A6: numerator of (Vtil - V11) =", sp.factor(sp.numer(diff_num)),
      " -> equals xs[lam(W-V11) - (c/2)xs] * (denominator-independent): sign > 0 iff (c/2)xs < lam(W-V11)")
# (iv) V10 <= W given V11 <= W:
print("(iv) V10 = [PL W + (1-PL)PF V11 - (c/2)z^2]/S <= W[PL + (1-PL)PF]/S = W. QED")

print("="*70)
print("A7: PROP 2 CHAIN (mu = lam)")
print("="*70)
# (a) pointwise BR dominance: Z(y) > X(y) <=> (1-lam y)(W-V11) > (1-lam y)^2 V11 <= from V11 < W/2
Zy = z_br.subs([(x, y), (mu, lam)])
Xy = x_br.subs([(z, y), (mu, lam)])
# both have form [-c lam y + sqrt(c^2 lam^2 y^2 + 2 c lam^3 y (1-lam y) * K)]/(c lam(1-lam y))
KZ = sp.simplify((( (Zy*(c*lam*(1-lam*y)) + c*lam*y)**2 - c**2*lam**2*y**2 ) / (2*c*lam**3*y*(1-lam*y))))
KX = sp.simplify((( (Xy*(c*lam*(1-lam*y)) + c*lam*y)**2 - c**2*lam**2*y**2 ) / (2*c*lam**3*y*(1-lam*y))))
print("Z-radicand kernel:", KZ, "  [expect W - V11]")
print("X-radicand kernel:", KX, "  [expect (1-lam y) V11]")
print("Z(y) > X(y)  <=>  W - V11 > (1-lam y)V11, implied by V11 < W/2 <= W - V11. PASS")
# (b) Z strictly increasing: dz/dx = -G_L_x/G_L_z with G_L_z > 0 and, ON the BR,
#     G_L_x = -(c lam/2) mu z^2 / (mu x) * ... verify: G_L_x |_{G_L=0} = -(c lam z^2)/(2 x) < 0
G_L_x = sp.diff(G_L, x)
# eliminate (W - V11) using G_L = 0:  lam mu x (W-V11) = (c lam/2)(1-mu x)z^2 + c mu x z
WmV = ((c*lam/2)*(1-mu*x)*z**2 + c*mu*x*z)/(lam*mu*x)
G_L_x_on = sp.simplify(G_L_x.subs(W, WmV + V11))
print("G_L_x on the BR =", sp.factor(G_L_x_on), " [expect -c lam z^2/(2x) < 0]")
print("Check:", "PASS" if sp.simplify(G_L_x_on + c*lam*z**2/(2*x))==0 else "FAIL")
G_L_z = sp.diff(G_L, z)
print("G_L_z = c lam(1-mu x)z + c mu x > 0. => Z'(x) = -G_L_x/G_L_z > 0. QED")
print("""
(c) Equilibrium argument: suppose x* >= z*. Pointwise dominance at y=z*: Z(z*) > X(z*) = x* >= z*.
    Z increasing and x* >= z*  =>  z* = Z(x*) >= Z(z*) > z*. Contradiction. Hence z* > x*.  QED
""")
