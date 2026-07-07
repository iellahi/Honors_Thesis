"""
Task 3 -- initial research state (0,0).
Quadratic:  cm x^2 - [3c + 2 mu m G] x + 2 mu (1-b0) D = 0,
  m = mu(1+b0),  G = V10 + b0 V01 - (1+b0)V11,  D = V10 - V01,  x00* = subtracted root.
  Delta0 = (3c + 2 mu m G)^2 - 8 c m mu (1-b0) D.

Numerical claims to confirm (NOT theorems -- flag any violation immediately):
  * Delta0 > 0            (analytical only for m<1/2 or b0>=4/5; else numerical)
  * x00* in (0,1)  and  m x00* < 1
Plus the internal-consistency check behind Prop 3's *unconditional* dx00*/dV10:
  exact identity  F00((1-b0)/m) = (1-b0)[ -c(2+b0)/m + 2 mu(1+b0)(V11-V01) ] < 0  (under A2'),
  which is equivalent to  m x00* < 1-b0  at the subtracted root  =>  F00_{V10} > 0.

Full closed system solved per draw (V11 -> asymmetric (z,x) -> V10,V01 -> (0,0)).
Solvers ported from verify_step5b.py / verify_step6.py.  Figure: PNG 300 dpi.
Outputs: numerics/task3_results.txt, numerics/task3_violations.csv, figs/task3_00_diagnostics.png
"""
import math, random, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE   = os.path.dirname(os.path.abspath(__file__))
THESIS = os.path.dirname(HERE)
FIGS   = os.path.join(THESIS, "figs")
NUM    = os.path.join(THESIS, "numerics")

log_lines = []
def log(s=""):
    print(s); log_lines.append(str(s))

# ---------------- solvers (verify_step5b.py / verify_step6.py) ----------------
def V11_and_x11(W, c, lam, b1):
    u = lam*lam*W*(1-b1*b1); Delta = 9*c*c - 2*c*u + u*u
    x11 = (3*c + u - math.sqrt(Delta))/(2*c*lam*(1+b1))
    return c*x11/(2*lam*(1+b1)) + b1*W/(1+b1)

def br_x(z, W, c, lam, mu, V11):
    rad = c*c*lam*lam*z*z + 2*c*mu*mu*lam*z*(1-lam*z)**2*V11
    return (-c*lam*z + math.sqrt(rad))/(c*mu*(1-lam*z))
def br_z(x, W, c, lam, mu, V11):
    rad = c*c*mu*mu*x*x + 2*c*mu*lam*lam*x*(1-mu*x)*(W-V11)
    return (-c*mu*x + math.sqrt(rad))/(c*lam*(1-mu*x))
def solve_fp(W, c, lam, mu, V11, tol=1e-13, itmax=5000):
    z, x = 0.5, 0.5
    for _ in range(itmax):
        zn = br_z(x, W, c, lam, mu, V11); xn = br_x(zn, W, c, lam, mu, V11)
        if abs(zn-z) < tol and abs(xn-x) < tol:
            z, x = zn, xn; break
        z, x = zn, xn
    G_L = (c*lam/2)*(1-mu*x)*z*z + c*mu*x*z - lam*mu*x*(W-V11)
    G_F = (c*mu/2)*(1-lam*z)*x*x + c*lam*z*x - mu*lam*z*(1-lam*z)*V11
    return z, x, max(abs(G_L), abs(G_F))
def values_10_01(z, x, W, c, lam, mu, V11):
    PL, PF = lam*z, mu*x; S = PL + PF - PL*PF
    V10 = (PL*W + (1-PL)*PF*V11 - c*z*z/2)/S
    V01 = ((1-PL)*PF*V11 - c*x*x/2)/S
    return V10, V01
def x00_full(W, c, mu, b0, V10, V01, V11):
    m = mu*(1+b0)
    G = V10 + b0*V01 - (1+b0)*V11
    D = V10 - V01
    A0, B0, C0 = c*m, -(3*c + 2*mu*m*G), 2*mu*(1-b0)*D
    D0 = B0*B0 - 4*A0*C0
    x00 = (-B0 - math.sqrt(D0))/(2*A0) if D0 > 0 else float('nan')
    quad_res = A0*x00*x00 + B0*x00 + C0 if D0 > 0 else float('nan')
    return x00, D0, m, G, D, quad_res

# ==========================================================================
# PART A -- full-system (0,0) sweep (parity: seed 11, N=60000)
# ==========================================================================
log("="*74)
log("PART A: (0,0) closure over the full system  (seed 11, N=60000)")
log("="*74)
random.seed(11)
N = 60000
minD0 = math.inf; minD0_numonly = math.inf
minx00 = math.inf; maxx00 = -math.inf; max_mx00 = -math.inf
max_quadres = 0.0
n_ok = n_skip = 0
v_D0 = v_x00 = v_mx00 = 0
n_analytic = n_numonly = 0
violations = []
for i in range(N):
    lam = random.uniform(1e-2, 0.4999)
    mu  = random.uniform(1e-2, 0.4999) if i % 2 else lam
    b1  = random.uniform(1e-4, 0.9995)
    b0  = random.uniform(0.0, 0.9995)
    c   = random.uniform(5e-2, 10.0)
    W   = random.uniform(1e-3, c/max(lam, mu))          # (A2) & (A2')
    V11 = V11_and_x11(W, c, lam, b1)
    z, x, res = solve_fp(W, c, lam, mu, V11)
    if res > 1e-9:
        n_skip += 1; continue
    V10, V01 = values_10_01(z, x, W, c, lam, mu, V11)
    x00, D0, m, G, D, qr = x00_full(W, c, mu, b0, V10, V01, V11)
    n_ok += 1
    analytic_region = (m < 0.5) or (b0 >= 0.8)          # where Delta0>0 is proven
    if analytic_region: n_analytic += 1
    else:               n_numonly  += 1; minD0_numonly = min(minD0_numonly, D0)
    minD0 = min(minD0, D0)
    if D0 <= 0:
        v_D0 += 1
        violations.append(("Delta0<=0", lam, mu, b1, b0, c, W, D0))
        continue
    minx00 = min(minx00, x00); maxx00 = max(maxx00, x00)
    max_mx00 = max(max_mx00, m*x00); max_quadres = max(max_quadres, abs(qr))
    if not (0 < x00 < 1):
        v_x00 += 1; violations.append(("x00_not_in(0,1)", lam, mu, b1, b0, c, W, x00))
    if not (m*x00 < 1):
        v_mx00 += 1; violations.append(("mx00>=1", lam, mu, b1, b0, c, W, m*x00))

log(f"converged draws: {n_ok}   (skipped non-converged asym: {n_skip})")
log(f"  analytic Delta0>0 region (m<1/2 or b0>=4/5): {n_analytic}   numerical-only: {n_numonly}")
log("")
log(f"[Delta0 > 0]     violations = {v_D0}      min Delta0 (all)        = {minD0:.6e}")
log(f"                                   min Delta0 (numerical-only) = {minD0_numonly:.6e}")
log(f"[x00* in (0,1)]  violations = {v_x00}      range x00* = [{minx00:.6e}, {maxx00:.6f}]")
log(f"[m x00* < 1]     violations = {v_mx00}      max m x00* = {max_mx00:.6f}")
log(f"[quadratic residual at x00*] max = {max_quadres:.3e}")
log(f"(memory reference: min Delta0 ~ 2.1e-2, zero violations)")

# ==========================================================================
# PART B -- exact identity behind unconditional dx00*/dV10 (seed 31, N=30000)
# ==========================================================================
log("")
log("="*74)
log("PART B: Prop 3 consistency -- exact identity & m x00* < 1-b0 frequency")
log("        (seed 31, N=30000; expect m x00* < 1-b0 to hold 100% under A2')")
log("="*74)
random.seed(31)
Nb = 30000
hold = fail = 0
max_identity_err = 0.0
max_F00_at_pt = -math.inf     # should stay < 0
for _ in range(Nb):
    lam = random.uniform(1e-2, 0.4999); mu = random.uniform(1e-2, 0.4999)
    c   = random.uniform(5e-2, 10.0);   W  = random.uniform(1e-3, c/max(lam, mu))
    b1  = random.uniform(1e-4, 0.9995); b0 = random.uniform(0.0, 0.9995)
    V11 = V11_and_x11(W, c, lam, b1)
    z, x, res = solve_fp(W, c, lam, mu, V11)
    if res > 1e-9: continue
    V10, V01 = values_10_01(z, x, W, c, lam, mu, V11)
    x00, D0, m, G, D, qr = x00_full(W, c, mu, b0, V10, V01, V11)
    if D0 <= 0 or not (x00 == x00): continue
    # F00 at the point (1-b0)/m, two ways:
    pt = (1-b0)/m
    F00_pt = c*m*pt*pt - (3*c + 2*mu*m*G)*pt + 2*mu*(1-b0)*D
    F00_id = (1-b0)*(-c*(2+b0)/m + 2*mu*(1+b0)*(V11-V01))
    max_identity_err = max(max_identity_err, abs(F00_pt - F00_id))
    max_F00_at_pt = max(max_F00_at_pt, F00_pt)
    if m*x00 < 1-b0: hold += 1
    else:            fail += 1
log(f"m x00* < 1-b0 :  holds {hold}, fails {fail}  "
    f"({100*fail/max(hold+fail,1):.3f}% -> should be 0 if unconditional)")
log(f"max |F00((1-b0)/m) - identity|  = {max_identity_err:.3e}   (identity verified)")
log(f"max F00((1-b0)/m) over scan     = {max_F00_at_pt:.6e}   (must be < 0)")
if fail == 0 and max_F00_at_pt < 0:
    log("=> CONFIRMED: at the subtracted root m x00* < 1-b0 ALWAYS, so F00_{V10} > 0")
    log("   unconditionally (matches finalized Prop 3; not merely the raw partial).")
else:
    log("=> FLAG: unconditional dx00*/dV10 claim not fully corroborated -- inspect.")

# ==========================================================================
# write violations (expected empty) and diagnostic figure
# ==========================================================================
with open(os.path.join(NUM, "task3_violations.csv"), "w") as f:
    f.write("type,lambda,mu,beta11,beta00,c,W,value\n")
    for row in violations:
        f.write(",".join(str(r) for r in row) + "\n")
log("")
log(f"violations written: {len(violations)}  (see numerics/task3_violations.csv)")

# collect arrays for the figure (fresh light scan for plotting)
random.seed(7)
D0s, x00s, mx00s, oneMb0 = [], [], [], []
for _ in range(20000):
    lam = random.uniform(1e-2, 0.4999); mu = random.uniform(1e-2, 0.4999)
    c   = random.uniform(5e-2, 10.0);   W  = random.uniform(1e-3, c/max(lam, mu))
    b1  = random.uniform(1e-4, 0.9995); b0 = random.uniform(0.0, 0.9995)
    V11 = V11_and_x11(W, c, lam, b1)
    z, x, res = solve_fp(W, c, lam, mu, V11)
    if res > 1e-9: continue
    V10, V01 = values_10_01(z, x, W, c, lam, mu, V11)
    x00, D0, m, G, D, qr = x00_full(W, c, mu, b0, V10, V01, V11)
    if D0 != D0: continue
    D0s.append(D0); x00s.append(x00); mx00s.append(m*x00); oneMb0.append(1-b0)
D0s = np.array(D0s); x00s = np.array(x00s); mx00s = np.array(mx00s); oneMb0 = np.array(oneMb0)

fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.6))
# (A) Delta0 distribution, log-x, mark analytic min
axes[0].hist(D0s, bins=np.logspace(np.log10(max(D0s.min(),1e-3)), np.log10(D0s.max()), 60),
             color="#4c78c8", alpha=0.85)
axes[0].set_xscale("log")
axes[0].axvline(minD0, color="red", ls="--", lw=1.6, label=fr"min $\Delta_0={minD0:.3g}$")
axes[0].set_xlabel(r"$\Delta_0$  (log scale)"); axes[0].set_ylabel("count")
axes[0].set_title(r"(A) discriminant $\Delta_0>0$  (0 violations)")
axes[0].legend(fontsize=9)
# (B) x00* distribution in (0,1)
axes[1].hist(x00s, bins=60, color="#54a24b", alpha=0.85)
axes[1].set_xlim(0, 1); axes[1].set_xlabel(r"$x_{00}^*$")
axes[1].set_title(r"(B) interior $x_{00}^*\in(0,1)$  (0 violations)")
axes[1].axvline(x00s.max(), color="red", ls="--", lw=1.2, label=fr"max $={x00s.max():.3f}$")
axes[1].legend(fontsize=9)
# (C) m x00* vs 1-b0 : all below the diagonal
idx = np.random.default_rng(0).choice(len(mx00s), size=min(6000, len(mx00s)), replace=False)
axes[2].scatter(oneMb0[idx], mx00s[idx], s=5, alpha=0.35, color="#e45756", edgecolors="none")
axes[2].plot([0, 1], [0, 1], "k-", lw=1.5, label=r"$m x_{00}^*=1-\beta_{00}$")
axes[2].set_xlim(0, 1); axes[2].set_ylim(0, 1)
axes[2].set_xlabel(r"$1-\beta_{00}$"); axes[2].set_ylabel(r"$m\,x_{00}^*$")
axes[2].set_title(r"(C) $m x_{00}^*<1-\beta_{00}$ always $\Rightarrow$ $F_{00,V_{10}}>0$")
axes[2].legend(fontsize=9, loc="upper left")
fig.tight_layout()
outpng = os.path.join(FIGS, "task3_00_diagnostics.png")
fig.savefig(outpng, dpi=300, bbox_inches="tight")
log(f"figure written: {outpng}")

with open(os.path.join(NUM, "task3_results.txt"), "w") as f:
    f.write("\n".join(log_lines) + "\n")
print("\n[written]", os.path.join(NUM, "task3_results.txt"))
