# Verify candidate analytical proof: dx*/dV11 > 0 at the asymmetric fixed point.
# System (draft eq. GL, GF):
#   G_L = (c*lam/2)(1-mu*x) z^2 + c*mu*x*z - lam*mu*x*(W - V)      = 0
#   G_F = (c*mu/2)(1-lam*z) x^2 + c*lam*z*x - mu*lam*z*(1-lam*z)*V = 0
# Claims to verify symbolically:
#  (1) G_L,x |_{G_L=0} = -c*lam*z^2/(2x)                    [draft identity]
#  (2) G_F,z |_{G_F=0} = -c*x*N/(2z(1-lam z)), N = mu*x(1-lam z)^2 - 2 lam^2 z^2
#  (3) Numerator of dx/dV (Cramer) = (c*lam*mu/(2z(1-lam z))) *
#        [ 2 lam^2 x^2 z^2 + 2 lam (1-mu x) z^3 (1-lam z)^2 + mu x (1-lam z)^2 (2z^2 - x^2) ]
#  (4) det J = c^2 * [ (lam(1-mu x)z + mu x)(mu(1-lam z)x + lam z)
#                      - lam*mu*x*z*(1-lam z)/4 + lam^3 z^3 / (2(1-lam z)) ]
# Then numeric spot-checks that det J > 0 and (3) > 0 when x < sqrt(2) z.

import sympy as sp

c, lam, mu, W, V, x, z = sp.symbols('c lam mu W V x z', positive=True)

G_L = sp.Rational(1,2)*c*lam*(1-mu*x)*z**2 + c*mu*x*z - lam*mu*x*(W - V)
G_F = sp.Rational(1,2)*c*mu*(1-lam*z)*x**2 + c*lam*z*x - mu*lam*z*(1-lam*z)*V

# eliminate (W-V) via G_L=0 and V via G_F=0
WmV_sol = sp.solve(sp.Eq(G_L, 0), W)[0] - V          # W - V on the leader BR
V_sol   = sp.solve(sp.Eq(G_F, 0), V)[0]              # V on the laggard BR

GLx = sp.diff(G_L, x)
GLz = sp.diff(G_L, z)
GFz = sp.diff(G_F, z)
GFx = sp.diff(G_F, x)
GLV = sp.diff(G_L, V)
GFV = sp.diff(G_F, V)

# (1)
GLx_br = sp.simplify(GLx.subs(W, WmV_sol + V))
chk1 = sp.simplify(GLx_br - (-c*lam*z**2/(2*x)))
print("(1) GLx|BR - (-c lam z^2/(2x)) =", chk1)

# (2)
N = mu*x*(1-lam*z)**2 - 2*lam**2*z**2
GFz_br = sp.simplify(GFz.subs(V, V_sol))
chk2 = sp.simplify(GFz_br - (-c*x*N/(2*z*(1-lam*z))))
print("(2) GFz|BR - (-c x N/(2z(1-lam z))) =", chk2)

# (3) Cramer numerator for dx/dV:  GLV*GFz - GLz*GFV   (on the BRs)
num = sp.simplify(GLV*GFz_br - GLz*GFV)
target = (c*lam*mu/(2*z*(1-lam*z))) * ( 2*lam**2*x**2*z**2
        + 2*lam*(1-mu*x)*z**3*(1-lam*z)**2
        + mu*x*(1-lam*z)**2*(2*z**2 - x**2) )
chk3 = sp.simplify(num - target)
print("(3) numerator - target =", chk3)

# (4) det J on the BRs
detJ = sp.simplify(GLz*GFx - GLx_br*GFz_br)
target4 = c**2*( (lam*(1-mu*x)*z + mu*x)*(mu*(1-lam*z)*x + lam*z)
                 - lam*mu*x*z*(1-lam*z)/4 + lam**3*z**3/(2*(1-lam*z)) )
chk4 = sp.simplify(detJ - target4)
print("(4) detJ - target =", chk4)

# (4b) positivity of detJ: expand product; the lone lam*mu*x*z term dominates the
# single negative term lam*mu*x*z*(1-lam z)/4 since (1-lam z)/4 < 1.
resid = sp.expand( (lam*(1-mu*x)*z + mu*x)*(mu*(1-lam*z)*x + lam*z)
                   - lam*mu*x*z*(1-lam*z)/4 )
print("(4b) expanded detJ/c^2 core (should have all + terms except the -1/4 piece):")
print(sp.collect(resid, [lam*mu*x*z]))

# numeric spot checks across a grid incl. mu != lam
import itertools, random, math
random.seed(0)
viol_det, viol_num, tested = 0, 0, 0
for _ in range(20000):
    lv = random.uniform(0.01, 0.49); mv = random.uniform(0.01, 0.49)
    cv = random.uniform(0.05, 5.0);  xv = random.uniform(1e-4, 0.999)
    zv = random.uniform(1e-4, 0.999)
    subs = {c: cv, lam: lv, mu: mv, x: xv, z: zv}
    dj = float(target4.subs(subs))
    if dj <= 0: viol_det += 1
    if xv < math.sqrt(2)*zv:
        tested += 1
        nm = float(target.subs(subs))
        if nm <= 0: viol_num += 1
print(f"detJ>0 violations: {viol_det}/20000; numerator>0 violations (x<sqrt2 z): {viol_num}/{tested}")
