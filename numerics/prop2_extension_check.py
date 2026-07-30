"""
Extend Prop. 3 (leader dominance z* > x*) from mu = lam to mu <= lam.
Filename note: "prop2" is the pre-renumbering label for leader dominance,
which is Proposition 3 in the final thesis. The filenames are kept because
Appendix C.3 cites them by name.


Strategy (equilibrium-level, NOT pointwise-for-all-y, which is FALSE for mu<lam):
Eliminate (W - V11) and V11 from the two BR quadratics at equilibrium by taking
the ratio G_L/G_F = 0 form.  This yields the exact equilibrium identity

  (dag)  z^2 (1-lam z) V11 [ (lam/2)(1-mu x) z + mu x ]
         = x^2 (W - V11) [ (mu/2)(1-lam z) x + lam z ]

Assume for contradiction x* >= z*.  Using Lemma 3 ( (W-V11)/V11 > 1 ) it suffices
to show the purely-geometric inequality (independent of W, V11):

  (sec)  x^2 [ (mu/2)(1-lam z) x + lam z ]  >=  z^2 (1-lam z) [ (lam/2)(1-mu x) z + mu x ]
         for all feasible x >= z > 0 and mu <= lam.

If (sec) holds then RHS(dag) > x^2 V11 [...] >= LHS(dag), contradicting equality,
so x* >= z* is impossible => z* > x*.

This script:
  (1) Reproduces z* > x* over mu<=lam draws with the reference solver.
  (2) Verifies (dag) is an exact equilibrium identity (residual ~ 0).
  (3) Verifies the KEY inequality (sec) over a large random grid of x>=z, mu<=lam.
  (4) Verifies the full logical chain at each equilibrium.
  (5) Reports the worst-case margin of (sec).

Thesis: Appendix C.3, numerical support for the leader-dominance proof.
"""
import random, math

# ---- reference solver (copied from verify_step5b.py) ----
def V11_and_x11(W, c, lam, b1):
    u = lam*lam*W*(1-b1*b1)
    Delta = 9*c*c - 2*c*u + u*u
    x11 = (3*c + u - math.sqrt(Delta)) / (2*c*lam*(1+b1))
    V11 = c*x11/(2*lam*(1+b1)) + b1*W/(1+b1)
    return V11, x11

def br_x(z, W, c, lam, mu, V11):
    rad = c*c*lam*lam*z*z + 2*c*mu*mu*lam*z*(1-lam*z)**2*V11
    return (-c*lam*z + math.sqrt(rad)) / (c*mu*(1-lam*z))

def br_z(x, W, c, lam, mu, V11):
    rad = c*c*mu*mu*x*x + 2*c*mu*lam*lam*x*(1-mu*x)*(W-V11)
    return (-c*mu*x + math.sqrt(rad)) / (c*lam*(1-mu*x))

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

def dag_LHS(z, x, lam, mu, V11):
    return z*z*(1-lam*z)*V11*((lam/2)*(1-mu*x)*z + mu*x)
def dag_RHS(z, x, lam, mu, W, V11):
    return x*x*(W-V11)*((mu/2)*(1-lam*z)*x + lam*z)

def sec_RHS(z, x, lam, mu):   # x^2 [ (mu/2)(1-lam z) x + lam z ]
    return x*x*((mu/2)*(1-lam*z)*x + lam*z)
def sec_LHS(z, x, lam, mu):   # z^2 (1-lam z)[ (lam/2)(1-mu x) z + mu x ]
    return z*z*(1-lam*z)*((lam/2)*(1-mu*x)*z + mu*x)

# =========================================================================
# (1) & (2) & (4): equilibrium reproduction + identity + chain, mu <= lam
# =========================================================================
random.seed(11)
N = 60000
fail_zx = 0
max_dag_resid = 0.0
fail_chain = 0
min_dag_gap = 1e18   # RHS(dag)-LHS(dag) at equilibrium; should be ~0 (identity) -> check separately
n_eq = 0
worst_ratio_eq = 0.0  # RatioFn = sec_LHS/sec_RHS at equilibrium; want < (W-V11)/V11
for i in range(N):
    lam_ = random.uniform(1e-2, 0.4999)
    mu_  = random.uniform(1e-2, lam_)          # enforce mu <= lam
    b1_  = random.uniform(1e-4, 0.9995)
    c_   = random.uniform(5e-2, 10.0)
    W_   = random.uniform(1e-3, c_/max(lam_, mu_))
    V11, x11 = V11_and_x11(W_, c_, lam_, b1_)
    z, x, res = solve_fp(W_, c_, lam_, mu_, V11)
    if res > 1e-9: continue
    n_eq += 1
    if not (z > x): fail_zx += 1
    # (2) dag identity residual (relative)
    L, R = dag_LHS(z, x, lam_, mu_, V11), dag_RHS(z, x, lam_, mu_, W_, V11)
    scale = max(abs(L), abs(R), 1e-30)
    max_dag_resid = max(max_dag_resid, abs(R-L)/scale)
    # (4) chain: RatioFn at eq vs (W-V11)/V11
    ratio = sec_LHS(z, x, lam_, mu_)/sec_RHS(z, x, lam_, mu_)
    bound = (W_-V11)/V11
    worst_ratio_eq = max(worst_ratio_eq, ratio - 1.0)  # want ratio < 1 <= bound
    if not (ratio < bound): fail_chain += 1

print("="*70)
print(f"(1)(2)(4) Equilibrium scan, mu<=lam, seed=11, N={N}, solved={n_eq}")
print(f"  z*>x* failures: {fail_zx}")
print(f"  max relative residual of identity (dag): {max_dag_resid:.3e}")
print(f"  chain failures (RatioFn < (W-V11)/V11): {fail_chain}")
print(f"  worst (RatioFn-1) at equilibrium: {worst_ratio_eq:.3e}  (want <0)")

# =========================================================================
# (3) KEY inequality (sec) over free grid x>=z>0, mu<=lam (NOT equilibrium)
# =========================================================================
random.seed(777)
M = 3_000_000
viol_sec = 0
min_margin = 1e18   # min of (sec_RHS - sec_LHS), want >= 0
argmin = None
for _ in range(M):
    lam_ = random.uniform(1e-4, 0.5)
    mu_  = random.uniform(1e-4, lam_)
    z_   = random.uniform(1e-4, min(0.999, 0.999/lam_))
    # x >= z, feasible: mu*x<1 and x<=1 (interior effort). scan x up to 1.
    x_   = random.uniform(z_, 1.0)
    if lam_*z_ >= 1 or mu_*x_ >= 1: continue
    m = sec_RHS(z_, x_, lam_, mu_) - sec_LHS(z_, x_, lam_, mu_)
    if m < min_margin:
        min_margin = m; argmin = (lam_, mu_, z_, x_)
    if m < 0: viol_sec += 1

print("="*70)
print(f"(3) KEY inequality (sec) over free grid x>=z, mu<=lam, N={M}")
print(f"  violations (sec_RHS < sec_LHS): {viol_sec}")
print(f"  min margin (sec_RHS - sec_LHS): {min_margin:.3e}  (want >= 0)")
print(f"  argmin (lam,mu,z,x): {argmin}")

# Also test near the boundary x=z and mu=lam (the tight corner)
print("  spot checks near x=z, mu=lam:")
for lam_ in (0.1, 0.3, 0.49):
    for s in (0.2, 0.6, 0.95):
        z_=x_=s
        if lam_*z_>=1: continue
        m = sec_RHS(z_, x_, lam_, lam_) - sec_LHS(z_, x_, lam_, lam_)
        # at x=z,mu=lam margin should equal sec_RHS*(lam*z) i.e. RatioFn=(1-lam z)
        ratio = sec_LHS(z_,x_,lam_,lam_)/sec_RHS(z_,x_,lam_,lam_)
        print(f"    lam=mu={lam_}, x=z={s}: margin={m:.3e}, RatioFn={ratio:.4f} (=1-lam*z={1-lam_*s:.4f})")
print("done")
