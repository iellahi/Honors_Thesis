"""Clean symbolic check of the (dag) equilibrium identity.
At equilibrium G_L=0 and G_F=0. Show sec_LHS*V - sec_RHS*(W-V) = 0 there.

Filename note: "prop2" is the pre-renumbering label for leader dominance,
which is Proposition 3 in the final thesis. The filenames are kept because
Appendix C.3 cites them by name.

Thesis: Appendix C.3, step 1 of the leader-dominance proof.
"""
import sympy as sp
z, x, lam, mu, W, V = sp.symbols('z x lam mu W V', positive=True)

G_L = sp.Rational(1,2)*lam*(1-mu*x)*z**2 + mu*x*z - lam*mu*x*(W-V)
G_F = sp.Rational(1,2)*mu*(1-lam*z)*x**2 + lam*z*x - mu*lam*z*(1-lam*z)*V

B_L = sp.Rational(1,2)*lam*(1-mu*x)*z + mu*x
B_R = sp.Rational(1,2)*mu*(1-lam*z)*x + lam*z
sec_LHS = z**2*(1-lam*z)*B_L
sec_RHS = x**2*B_R

E = sec_LHS*V - sec_RHS*(W - V)          # want = 0 at equilibrium

# Substitute W from G_L=0, then V from G_F=0
W_sol = sp.solve(G_L, W)[0]
V_sol = sp.solve(G_F, V)[0]
E_sub = E.subs(W, W_sol).subs(V, V_sol)
print("sec_LHS*V - sec_RHS*(W-V) at equilibrium =", sp.simplify(E_sub))

# Cross-check the OTHER way: express as combination and reduce modulo the ideal
# Using groebner/reduce: E should reduce to 0 modulo {G_L, G_F}
GB = sp.groebner([G_L, G_F], W, V, order='lex')
red = GB.reduce(E)[1]
print("E reduced modulo <G_L, G_F> =", sp.simplify(red))
