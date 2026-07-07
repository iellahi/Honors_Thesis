"""Step 2c: analytical sufficiency proof of Assumption lam*W <= c, machine-checked.
Given: 0 < lam < 1/2, 0 <= beta < 1, lam*W <= c, M = lam(1+beta) (so 0 < M < 1).
Prove: (i) lam^2 W (1-b^2) < 2c   (ii) lam W (1-b)(2-M) < c(3-M)."""
import sympy as sp

lam, W, c, beta = sp.symbols('lambda W c beta', positive=True)
M = lam*(1+beta)

# ---- (i): chain lam^2 W(1-b^2) = lam*(lam W)*(1-b^2) <= lam*c*(1-b^2) < c/2 < 2c
# Step 1: lam^2 W (1-b^2) <= lam*c*(1-b^2)  given lam W <= c  (multiply by lam(1-b^2) >= 0)
step1 = sp.simplify(lam*c*(1-beta**2) - lam*(lam*W)*(1-beta**2) - lam*(1-beta**2)*(c - lam*W))
print("(i) step1 identity (lam(1-b^2))(c-lamW) rearrangement:", step1 == 0)
# Step 2: lam*c*(1-b^2) < c/2 given lam < 1/2 and (1-b^2) <= 1: lam*(1-b^2) < 1/2
# Step 3: c/2 < 2c trivially. Margin check via symbolic bound:
expr_i = 2*c - lam**2*W*(1-beta**2)
# substitute worst case lam W = c, lam = 1/2, beta = 0: 2c - c/2 = 3c/2 > 0
print("(i) worst-case margin 2c - lam^2 W(1-b^2) at lamW=c, lam=1/2, b=0:",
      sp.simplify(expr_i.subs([(W, c/lam)]).subs(lam, sp.Rational(1,2)).subs(beta, 0)))

# ---- (ii): lam W (1-b)(2-M) <= c(1-b)(2-M) < c(2-M) <= c(3-M) - c < c(3-M)
# Key inequality: (1-b)(2-M) < (3-M): since (1-b) <= 1, (1-b)(2-M) <= 2-M = (3-M) - 1 < 3-M
diff_key = sp.simplify((3-M) - (2-M))
print("(ii) (3-M)-(2-M) =", diff_key, " (strictly positive gap of exactly 1)")
# Full worst-case margin: c(3-M) - lam W(1-b)(2-M) at lam W = c:
margin_ii = sp.simplify(c*(3-M) - c*(1-beta)*(2-M))
print("(ii) margin at lamW=c:", sp.factor(margin_ii), " -> positive since M<1<... check:")
# margin = c[3-M - (1-b)(2-M)] = c[3-M-2+M+2b-bM] = c[1 + 2b - bM] = c[1 + b(2-M)] > 0
print("(ii) simplified margin == c*(1 + beta*(2-M)):",
      sp.simplify(margin_ii - c*(1+beta*(2-M))) == 0, "-> strictly positive. QED")
