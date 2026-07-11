"""
Step 8: transmission package for Section 5.5 / Appendix C.
Independent rebuild of the 7/10 working results (R1)-(R4). Nothing is reused
from verify_dx_dV11.py; all identities re-derived from the BR quadratics.

System (draft eqs. GL, GF; leader z, laggard x; V = V11 constant from (1,1)):
  G_L(z;x) = (c lam/2)(1-mu x) z^2 + c mu x z - lam mu x (W - V)       = 0
  G_F(x;z) = (c mu/2)(1-lam z) x^2 + c lam z x - mu lam z (1-lam z) V  = 0

Verified here (symbolic zero required for each):
  (I1) G_L,x |_{G_L=0} = -c lam z^2 / (2x)
  (I2) G_F,z |_{G_F=0} = -c x N / (2 z (1-lam z)),  N = mu x(1-lam z)^2 - 2 lam^2 z^2
  (R1) det J = c^2 [ (lam(1-mu x)z + mu x)(mu(1-lam z)x + lam z)
                     - lam mu x z (1-lam z)/4 + lam^3 z^3 / (2(1-lam z)) ]
       + positivity decomposition: det J / c^2 equals a sum of terms each >= 0
       under (A1), the paired term being lam mu x z (3+lam z)/4.
  (R2) Cramer numerator for dx*/dV:
       num_x = G_L,V G_F,z - G_L,z G_F,V
             = (c lam mu/(2z(1-lam z))) [ 2 lam^2 x^2 z^2
               + 2 lam (1-mu x) z^3 (1-lam z)^2 + mu x (1-lam z)^2 (2z^2 - x^2) ]
       => num_x > 0 whenever x < sqrt(2) z  (sufficient condition).
  (R4) Cramer numerator for dz*/dV:
       num_z = -G_L,V G_F,x + G_L,x G_F,V
             = (c lam mu / x) [ lam z^3 (1-lam z)/2 - mu (1-lam z) x^3 - lam z x^2 ]

Numeric spot checks AT SOLVED EQUILIBRIA (not free (x,z) points):
  seed 8, N = 20000 draws in (A1)-(A2'), solver = verify_step5b fixed-point
  iteration (tol 1e-13, residual tol 1e-9). Checks: det J > 0; Cramer
  dx*/dV11, dz*/dV11 vs central finite differences (rel err); sign(num_x) > 0
  when x* < sqrt(2) z*; sign(dz*/dV11) = sign(R4 bracket).
"""
import sympy as sp

c, lam, mu, W, V, x, z = sp.symbols('c lam mu W V x z', positive=True)

G_L = sp.Rational(1, 2)*c*lam*(1 - mu*x)*z**2 + c*mu*x*z - lam*mu*x*(W - V)
G_F = sp.Rational(1, 2)*c*mu*(1 - lam*z)*x**2 + c*lam*z*x - mu*lam*z*(1 - lam*z)*V

# raw partials
GLz, GLx, GLV = sp.diff(G_L, z), sp.diff(G_L, x), sp.diff(G_L, V)
GFx, GFz, GFV = sp.diff(G_F, x), sp.diff(G_F, z), sp.diff(G_F, V)

# on-BR substitutions: W solved from G_L = 0, V solved from G_F = 0
W_on_L = sp.solve(sp.Eq(G_L, 0), W)[0]
V_on_F = sp.solve(sp.Eq(G_F, 0), V)[0]

fails = 0
def check(label, expr):
    global fails
    r = sp.simplify(expr)
    ok = (r == 0)
    fails += (not ok)
    print(f"{label}: {'PASS (symbolic zero)' if ok else 'FAIL -> ' + str(r)}")

# (I1)
GLx_br = sp.simplify(GLx.subs(W, W_on_L))
check("(I1) G_L,x|BR + c lam z^2/(2x)", GLx_br + c*lam*z**2/(2*x))

# (I2)
N = mu*x*(1 - lam*z)**2 - 2*lam**2*z**2
GFz_br = sp.simplify(GFz.subs(V, V_on_F))
check("(I2) G_F,z|BR + c x N/(2z(1-lam z))", GFz_br + c*x*N/(2*z*(1 - lam*z)))

# (R1) det J at the fixed point (GLz, GFx contain no W or V)
detJ = sp.simplify(GLz*GFx - GLx_br*GFz_br)
core = ((lam*(1 - mu*x)*z + mu*x)*(mu*(1 - lam*z)*x + lam*z)
        - lam*mu*x*z*(1 - lam*z)/4 + lam**3*z**3/(2*(1 - lam*z)))
check("(R1) det J - c^2 * core", detJ - c**2*core)

# (R1) positivity: core == sum of manifestly nonnegative terms under (A1)
pos_sum = (lam*mu*x*z*(1 - mu*x)*(1 - lam*z)      # >= 0
           + lam**2*z**2*(1 - mu*x)               # >= 0
           + mu**2*x**2*(1 - lam*z)               # >= 0
           + lam*mu*x*z*(3 + lam*z)/4             # > 0  (pairs bd with -1/4 term)
           + lam**3*z**3/(2*(1 - lam*z)))         # > 0
check("(R1) core - positive-term decomposition", core - pos_sum)

# (R2) Cramer numerator for dx/dV
num_x = sp.simplify(GLV*GFz_br - GLz*GFV)
target_x = (c*lam*mu/(2*z*(1 - lam*z)))*(2*lam**2*x**2*z**2
            + 2*lam*(1 - mu*x)*z**3*(1 - lam*z)**2
            + mu*x*(1 - lam*z)**2*(2*z**2 - x**2))
check("(R2) num_x - target", num_x - target_x)

# (R4) Cramer numerator for dz/dV
num_z = sp.simplify(-GLV*GFx + GLx_br*GFV)
target_z = (c*lam*mu/x)*(lam*z**3*(1 - lam*z)/2 - mu*(1 - lam*z)*x**3 - lam*z*x**2)
check("(R4) num_z - target", num_z - target_z)

print(f"\nSymbolic block: {'ALL PASS' if fails == 0 else str(fails) + ' FAILURES'}\n")

# ---------------------------------------------------------------- numeric part
import random, math

def V11_and_x11(W_, c_, l_, b1):
    u = l_*l_*W_*(1 - b1*b1)
    Delta = 9*c_*c_ - 2*c_*u + u*u
    x11 = (3*c_ + u - math.sqrt(Delta))/(2*c_*l_*(1 + b1))
    return c_*x11/(2*l_*(1 + b1)) + b1*W_/(1 + b1), x11

def br_x(z_, W_, c_, l_, m_, V_):
    rad = c_*c_*l_*l_*z_*z_ + 2*c_*m_*m_*l_*z_*(1 - l_*z_)**2*V_
    return (-c_*l_*z_ + math.sqrt(rad))/(c_*m_*(1 - l_*z_))

def br_z(x_, W_, c_, l_, m_, V_):
    rad = c_*c_*m_*m_*x_*x_ + 2*c_*m_*l_*l_*x_*(1 - m_*x_)*(W_ - V_)
    return (-c_*m_*x_ + math.sqrt(rad))/(c_*l_*(1 - m_*x_))

def solve_fp(W_, c_, l_, m_, V_, tol=1e-13, itmax=5000):
    z_, x_ = 0.5, 0.5
    for _ in range(itmax):
        zn = br_z(x_, W_, c_, l_, m_, V_)
        xn = br_x(zn, W_, c_, l_, m_, V_)
        if abs(zn - z_) < tol and abs(xn - x_) < tol:
            z_, x_ = zn, xn
            break
        z_, x_ = zn, xn
    G_L_ = (c_*l_/2)*(1 - m_*x_)*z_*z_ + c_*m_*x_*z_ - l_*m_*x_*(W_ - V_)
    G_F_ = (c_*m_/2)*(1 - l_*z_)*x_*x_ + c_*l_*z_*x_ - m_*l_*z_*(1 - l_*z_)*V_
    return z_, x_, max(abs(G_L_), abs(G_F_))

def detJ_num(c_, l_, m_, x_, z_):
    return c_*c_*((l_*(1 - m_*x_)*z_ + m_*x_)*(m_*(1 - l_*z_)*x_ + l_*z_)
                  - l_*m_*x_*z_*(1 - l_*z_)/4 + l_**3*z_**3/(2*(1 - l_*z_)))

def numx_num(c_, l_, m_, x_, z_):
    return (c_*l_*m_/(2*z_*(1 - l_*z_)))*(2*l_*l_*x_*x_*z_*z_
            + 2*l_*(1 - m_*x_)*z_**3*(1 - l_*z_)**2
            + m_*x_*(1 - l_*z_)**2*(2*z_*z_ - x_*x_))

def numz_num(c_, l_, m_, x_, z_):
    return (c_*l_*m_/x_)*(l_*z_**3*(1 - l_*z_)/2 - m_*(1 - l_*z_)*x_**3 - l_*z_*x_*x_)

random.seed(8)
N = 20000
RES_TOL, FD_FLAG = 1e-9, 1e-7
viol_det = viol_numx = viol_numz_sign = 0
n_sqrt2, n_fail_sqrt2 = 0, 0
fd_x_bad = fd_z_bad = flagged = solved = 0
max_relerr_x = max_relerr_z = 0.0

for i in range(N):
    l_ = random.uniform(1e-2, 0.4999)
    m_ = random.uniform(1e-2, 0.4999) if i % 2 else l_
    b1 = random.uniform(1e-4, 0.9995)
    c_ = random.uniform(5e-2, 10.0)
    W_ = random.uniform(1e-3, c_/max(l_, m_))
    V_, _ = V11_and_x11(W_, c_, l_, b1)
    z_, x_, res = solve_fp(W_, c_, l_, m_, V_)
    if res > RES_TOL:
        continue
    solved += 1
    dj = detJ_num(c_, l_, m_, x_, z_)
    nx = numx_num(c_, l_, m_, x_, z_)
    nz = numz_num(c_, l_, m_, x_, z_)
    if dj <= 0: viol_det += 1
    if x_ < math.sqrt(2)*z_:
        n_sqrt2 += 1
        if nx <= 0: viol_numx += 1
    # central FD in V11 vs Cramer
    h = 1e-6*max(V_, 1e-4)
    if V_ - h > 0 and V_ + h < W_:
        zp, xp, rp = solve_fp(W_, c_, l_, m_, V_ + h)
        zm, xm, rm = solve_fp(W_, c_, l_, m_, V_ - h)
        if rp < RES_TOL and rm < RES_TOL:
            fd_x, fd_z = (xp - xm)/(2*h), (zp - zm)/(2*h)
            cr_x, cr_z = nx/dj, nz/dj
            if abs(fd_x) > FD_FLAG:
                rel = abs(fd_x - cr_x)/max(abs(fd_x), 1e-300)
                max_relerr_x = max(max_relerr_x, rel)
                if rel > 1e-3: fd_x_bad += 1
            else:
                flagged += 1
            if abs(fd_z) > FD_FLAG:
                rel = abs(fd_z - cr_z)/max(abs(fd_z), 1e-300)
                max_relerr_z = max(max_relerr_z, rel)
                if rel > 1e-3: fd_z_bad += 1
                if (fd_z > 0) != (nz > 0): viol_numz_sign += 1
            else:
                flagged += 1

print(f"Numeric spot checks at solved equilibria: seed 8, N={N}, solved={solved}")
print(f"  det J <= 0 violations:                    {viol_det}/{solved}")
print(f"  num_x <= 0 with x* < sqrt2 z*:            {viol_numx}/{n_sqrt2}")
print(f"  Cramer vs FD rel-err > 1e-3 (x, z):       {fd_x_bad}, {fd_z_bad} "
      f"(max relerr {max_relerr_x:.2e}, {max_relerr_z:.2e})")
print(f"  sign(dz*/dV11) != sign(R4 bracket):       {viol_numz_sign}")
print(f"  near-zero FD flags skipped:               {flagged}")
