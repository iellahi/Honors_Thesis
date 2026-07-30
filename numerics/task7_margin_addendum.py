"""
Sharpness of the sqrt(2) condition.

Proposition 5 gives dx*/dV11 > 0 whenever x* < sqrt(2) z*. The region-wide scan
in task7_transmission_scans.py shows the condition holds everywhere it looks, so
the question left is whether the constant sqrt(2) can be improved. Three probes,
each pushing the ratio x*/(sqrt2 z*) as high as it will go:

  (a) a random scan on a fresh seed, for a baseline maximum;
  (b) an adversarial sweep over the corners of the parameter box;
  (c) a limit probe driving lam -> 0, b11 -> 1, mu -> 1/2, where the ratio is
      largest, to see whether it approaches 1.

The ratio climbs toward 1 without reaching it, so sqrt(2) is the right constant:
no smaller one would cover the region.

Solvers ported from verify_scripts/verify_step5b.py.
Outputs: numerics/task7_margin_addendum.txt

Thesis: Section 5.5, the sharpness of the sqrt(2) condition (Appendix C.5).
"""
import math, random, os

HERE   = os.path.dirname(os.path.abspath(__file__))
THESIS = os.path.dirname(HERE)

log_lines = []
def log(s=""):
    print(s); log_lines.append(str(s))

# ---------------- solvers (reference: verify_step5b.py) ----------------
def V11_and_x11(W, c, lam, b1):
    u = lam*lam*W*(1-b1*b1); Delta = 9*c*c - 2*c*u + u*u
    x11 = (3*c + u - math.sqrt(Delta))/(2*c*lam*(1+b1))
    return c*x11/(2*lam*(1+b1)) + b1*W/(1+b1), x11

def br_x(z, W, c, lam, mu, V11):
    rad = c*c*lam*lam*z*z + 2*c*mu*mu*lam*z*(1-lam*z)**2*V11
    return (-c*lam*z + math.sqrt(rad))/(c*mu*(1-lam*z))

def br_z(x, W, c, lam, mu, V11):
    rad = c*c*mu*mu*x*x + 2*c*mu*lam*lam*x*(1-mu*x)*(W-V11)
    return (-c*mu*x + math.sqrt(rad))/(c*lam*(1-mu*x))

def solve_fp(W, c, lam, mu, V11, tol=1e-13, itmax=5000):
    z, x = 0.5, 0.5
    for _ in range(itmax):
        zn = br_z(x, W, c, lam, mu, V11)
        xn = br_x(zn, W, c, lam, mu, V11)
        if abs(zn-z) < tol and abs(xn-x) < tol:
            z, x = zn, xn; break
        z, x = zn, xn
    G_L = (c*lam/2)*(1-mu*x)*z*z + c*mu*x*z - lam*mu*x*(W-V11)
    G_F = (c*mu/2)*(1-lam*z)*x*x + c*lam*z*x - mu*lam*z*(1-lam*z)*V11
    return z, x, max(abs(G_L), abs(G_F))

def full_solve(W, c, lam, mu, b1):
    V11, _ = V11_and_x11(W, c, lam, b1)
    z, x, res = solve_fp(W, c, lam, mu, V11)
    return V11, z, x, res

SQ2 = math.sqrt(2)

def ratio(W, c, lam, mu, b1):
    """x*/(sqrt2 z*) at the interior equilibrium, or None if it did not solve."""
    _, z, x, res = full_solve(W, c, lam, mu, b1)
    if res > 1e-9 or z <= 0:
        return None
    return x/(SQ2*z)

# ==========================================================================
# (a) random scan on a fresh seed
# ==========================================================================
random.seed(43)
N = 60000
best, arg = 0.0, None
for _ in range(N):
    lam = random.uniform(1e-2, 0.4999)
    mu  = random.uniform(1e-2, 0.4999)
    b1  = random.uniform(1e-4, 0.9995)
    c   = random.uniform(5e-2, 10.0)
    W   = random.uniform(1e-3, c/max(lam, mu))
    r = ratio(W, c, lam, mu, b1)
    if r is not None and r > best:
        best, arg = r, (lam, mu, b1, c, W)
log(f"random scan (seed 43, N={N}): max x*/(sqrt2 z*) = {best:.6f}")
log("  argmax (lam,mu,b1,c,W) = ({:.4f}, {:.4f}, {:.4f}, {:.4f}, {:.4f})".format(*arg))

# ==========================================================================
# (b) adversarial sweep over the corners of the parameter box
# ==========================================================================
LAMS  = [0.01, 0.05, 0.1, 0.25, 0.4999]
MUS   = [0.01, 0.25, 0.4999]
B1S   = [1e-4, 0.5, 0.9, 0.9995]
CS    = [0.05, 1.0, 10.0]
WFRAC = [0.05, 0.5, 0.95]          # W as a fraction of its upper bound c/max(lam,mu)

best_c, arg_c = 0.0, None
for lam in LAMS:
    for mu in MUS:
        for b1 in B1S:
            for c in CS:
                for wf in WFRAC:
                    W = wf * c/max(lam, mu)
                    r = ratio(W, c, lam, mu, b1)
                    if r is not None and r > best_c:
                        best_c, arg_c = r, (lam, mu, b1, c, wf)
log(f"adversarial corner sweep:      max x*/(sqrt2 z*) = {best_c:.6f}")
log("  argmax (lam,mu,b1,c,W-frac) = ({}, {}, {}, {}, {})".format(*arg_c))

# ==========================================================================
# (c) limit probe: lam -> 0, b11 -> 1, at mu just under 1/2
# ==========================================================================
MU_L, C_L, WFRAC_L = 0.4999, 1.0, 0.05
log(f"limit probe: ratio x*/(sqrt2 z*) as (lam, 1-b1, Wfrac) -> 0, mu={MU_L}, c={C_L}")
for lam in (1e-2, 1e-3, 1e-4, 1e-5):
    for gap in (1e-2, 1e-4, 1e-6):
        b1 = 1.0 - gap
        W  = WFRAC_L * C_L/max(lam, MU_L)
        r  = ratio(W, C_L, lam, MU_L, b1)
        log(f"  lam={lam:.0e} 1-b1={gap:.0e}: ratio={r:.8f}" if r is not None
            else f"  lam={lam:.0e} 1-b1={gap:.0e}: no interior solution")

out = os.path.join(THESIS, "numerics", "task7_margin_addendum.txt")
with open(out, "w") as f:
    f.write("\n".join(log_lines) + "\n")
print("\n[written]", os.path.relpath(out, THESIS))
