"""Step 2 supplement: root validity + interiority conditions, stated precisely."""
import sympy as sp
import itertools, random

x = sp.symbols('x', positive=True)
W, c, lam, beta = sp.symbols('W c lambda beta_11', positive=True)
M = lam*(1+beta)
A = c*lam*(1+beta); B = -(3*c + lam**2*W*(1-beta**2)); C = 2*lam*W*(1-beta)
F = A*x**2 + B*x + C
Delta = 9*c**2 - 2*c*lam**2*W*(1-beta**2) + lam**4*W**2*(1-beta**2)**2
x_minus = (3*c + lam**2*W*(1-beta**2) - sp.sqrt(Delta))/(2*c*lam*(1+beta))

print("CLAIM A: F(1) < 0  <=>  lam*W*(1-beta)*(2-M) < c*(3-M)   [x- < 1 < x+]")
F1 = sp.simplify(F.subs(x,1) - (c*(M-3) + lam*W*(1-beta)*(2-M)))
print("F(1) identity:", "PASS" if F1 == 0 else f"FAIL {F1}")

print("\nCLAIM B: lam*W <= c is SUFFICIENT for both (i) x- < 1/M < x+ and (ii) x- < 1")
print("  (i) needs lam^2 W (1-b^2) < 2c:  lam^2 W(1-b^2) <= lam*c*(1-b^2)/1 < c/2*... check numerically")
print("  (ii) needs lam W(1-b)(2-M) < c(3-M)")
# Numerical scan over the feasible box
random.seed(0)
viol_i = viol_ii = viol_x1 = 0; n = 200000
for _ in range(n):
    lv = random.uniform(1e-3, 0.499)
    bv = random.uniform(0, 0.999)
    cv = random.uniform(1e-3, 10)
    Wv = random.uniform(1e-3, cv/lv)   # enforce lam*W <= c
    Mv = lv*(1+bv)
    if not (lv**2*Wv*(1-bv**2) < 2*cv): viol_i += 1
    if not (lv*Wv*(1-bv)*(2-Mv) < cv*(3-Mv)): viol_ii += 1
    # direct check of x- < 1
    Dv = 9*cv**2 - 2*cv*lv**2*Wv*(1-bv**2) + lv**4*Wv**2*(1-bv**2)**2
    xm = (3*cv + lv**2*Wv*(1-bv**2) - Dv**0.5)/(2*cv*lv*(1+bv))
    if not (xm < 1): viol_x1 += 1
print(f"  scan n={n}: violations (i)={viol_i}, (ii)={viol_ii}, direct x-<1: {viol_x1}")

print("\nCLAIM C: WITHOUT lam*W <= c, x- can exceed 1 (corner regime exists)")
# try large W
lv, bv, cv, Wv = 0.49, 0.0, 0.1, 50.0
Dv = 9*cv**2 - 2*cv*lv**2*Wv + lv**4*Wv**2
xm = (3*cv + lv**2*Wv - Dv**0.5)/(2*cv*lv)
print(f"  lam={lv}, beta={bv}, c={cv}, W={Wv}: x- = {xm:.4f}  -> corner if > 1")

print("\nCLAIM D: at beta=1, C=0 so x- = 0 exactly (full free-riding boundary)")
print("  x-(beta=1) =", sp.simplify(x_minus.subs(beta,1)))

print("\nCLAIM E: both roots strictly positive for beta<1 (product C/A>0, sum -B/A>0)")
print("  C/A =", sp.simplify(C/A), " -B/A =", sp.simplify(-B/A))
