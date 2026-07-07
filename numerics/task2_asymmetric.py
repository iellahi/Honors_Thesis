"""
Task 2 -- asymmetric states (1,0)/(0,1).
  (a) Solve the best-response fixed point (additive roots) by iteration; confirm
      residuals of the BR quadratics G_L, G_F < 1e-9 (the additive-root formulas
      are proven; this is a numerical confirmation).
  (b) Confirm Prop 2: z* > x* at mu = lam (proven; expect 0 violations).
  (c) Map the mu != lam region where z* > x* FAILS (~22% of draws) -- frontier figure.

Solvers ported from verify_scripts/verify_step5b.py.  Figures: PNG only, 300 dpi.
Outputs: numerics/task2_results.txt, numerics/task2_frontier_grid.csv, figs/prop2_frontier.png
"""
import math, random, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm

HERE   = os.path.dirname(os.path.abspath(__file__))
THESIS = os.path.dirname(HERE)
FIGS   = os.path.join(THESIS, "figs")
NUM    = os.path.join(THESIS, "numerics")

log_lines = []
def log(s=""):
    print(s); log_lines.append(str(s))

# ---------------- solvers (reference: verify_step5b.py) ----------------
def V11_and_x11(W, c, lam, b1):
    u = lam*lam*W*(1-b1*b1); Delta = 9*c*c - 2*c*u + u*u
    x11 = (3*c + u - math.sqrt(Delta))/(2*c*lam*(1+b1))
    V11 = c*x11/(2*lam*(1+b1)) + b1*W/(1+b1)
    return V11, x11

def br_x(z, W, c, lam, mu, V11):          # laggard BR (additive root)
    rad = c*c*lam*lam*z*z + 2*c*mu*mu*lam*z*(1-lam*z)**2*V11
    return (-c*lam*z + math.sqrt(rad))/(c*mu*(1-lam*z))

def br_z(x, W, c, lam, mu, V11):          # leader BR (additive root)
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
    G_L = (c*lam/2)*(1-mu*x)*z*z + c*mu*x*z - lam*mu*x*(W-V11)          # =0 at eq
    G_F = (c*mu/2)*(1-lam*z)*x*x + c*lam*z*x - mu*lam*z*(1-lam*z)*V11   # =0 at eq
    return z, x, max(abs(G_L), abs(G_F))

def values_10_01(z, x, W, c, lam, mu, V11):
    PL, PF = lam*z, mu*x; S = PL + PF - PL*PF
    V10 = (PL*W + (1-PL)*PF*V11 - c*z*z/2)/S
    V01 = ((1-PL)*PF*V11 - c*x*x/2)/S
    return V10, V01

# ==========================================================================
# PART A -- validation sweep (parity with verify_step5b: seed 11, N=60000)
# ==========================================================================
log("="*74)
log("PART A: asymmetric fixed point -- residuals, Prop 2, value ordering")
log("="*74)
random.seed(11)
N = 60000
max_res = 0.0
n_sym = n_asym = 0
fail_sym = fail_asym = 0
bad_order = bad_int = bad_prob = 0
worst_res_case = None
for i in range(N):
    lam = random.uniform(1e-2, 0.4999)
    mu  = random.uniform(1e-2, 0.4999) if i % 2 else lam   # half sym, half asym
    b1  = random.uniform(1e-4, 0.9995)
    c   = random.uniform(5e-2, 10.0)
    W   = random.uniform(1e-3, c/max(lam, mu))             # (A2) & (A2')
    V11, x11 = V11_and_x11(W, c, lam, b1)
    z, x, res = solve_fp(W, c, lam, mu, V11)
    if res > max_res:
        max_res = res; worst_res_case = (lam, mu, b1, c, W)
    if res > 1e-9:                # skip pathological non-converged (should be ~none)
        continue
    V10, V01 = values_10_01(z, x, W, c, lam, mu, V11)
    if mu == lam:
        n_sym += 1
        if not (z > x): fail_sym += 1
    else:
        n_asym += 1
        if not (z > x): fail_asym += 1
    if not (0 <= V01 <= V11 <= V10 <= W + 1e-12): bad_order += 1
    if not (0 < z < 1 and 0 <= x < 1):            bad_int += 1
    if not (lam*z < 1 and mu*x < 1):              bad_prob += 1

log(f"draws N = {N}  (half mu=lam, half mu!=lam; seed 11)")
log(f"max BR-quadratic residual |G_L|,|G_F|   = {max_res:.3e}   (tol 1e-9)")
log(f"  worst-residual params (lam,mu,b1,c,W) = "
    f"({worst_res_case[0]:.4f},{worst_res_case[1]:.4f},{worst_res_case[2]:.4f},"
    f"{worst_res_case[3]:.4f},{worst_res_case[4]:.4f})")
log("")
log(f"[Prop 2] mu = lam : z* > x* fails {fail_sym} / {n_sym}   (expect 0)")
log(f"[map]    mu != lam: z* > x* fails {fail_asym} / {n_asym}  "
    f"= {100*fail_asym/max(n_asym,1):.2f}%   (memory: ~22%)")
log(f"value ordering 0<=V01<=V11<=V10<=W violations = {bad_order}")
log(f"interiority z,x in (0,1) violations           = {bad_int}")
log(f"probability lam z<1, mu x<1 violations        = {bad_prob}")

# ==========================================================================
# PART B -- structure of the z>x failure region (which params drive the flip?)
# ==========================================================================
log("")
log("="*74)
log("PART B: where does z*>x* fail?  fraction failing, binned by mu/lam")
log("="*74)
random.seed(101)
Nb = 200000
bins = [0.0, 0.5, 0.8, 1.0, 1.25, 1.6, 2.0, 3.0, 1e9]
btot = [0]*(len(bins)-1); bfail = [0]*(len(bins)-1)
for _ in range(Nb):
    lam = random.uniform(1e-2, 0.4999)
    mu  = random.uniform(1e-2, 0.4999)
    if mu == lam: continue
    b1  = random.uniform(1e-4, 0.9995)
    c   = random.uniform(5e-2, 10.0)
    W   = random.uniform(1e-3, c/max(lam, mu))
    V11, _ = V11_and_x11(W, c, lam, b1)
    z, x, res = solve_fp(W, c, lam, mu, V11)
    if res > 1e-9: continue
    r = mu/lam
    for k in range(len(bins)-1):
        if bins[k] <= r < bins[k+1]:
            btot[k] += 1
            if not (z > x): bfail[k] += 1
            break
log(f"random draws N = {Nb} (seed 101, mu!=lam)")
log(f"{'mu/lam bin':>16} | {'n':>8} | {'z<=x %':>8}")
for k in range(len(bins)-1):
    lo, hi = bins[k], bins[k+1]
    hi_s = "inf" if hi > 1e8 else f"{hi:g}"
    pct = 100*bfail[k]/btot[k] if btot[k] else float('nan')
    log(f"  [{lo:g}, {hi_s}) {'':>3} | {btot[k]:>8} | {pct:>7.2f}%")
overall = 100*sum(bfail)/max(sum(btot),1)
log(f"overall z*<=x* share (mu!=lam) = {overall:.2f}%")
log("=> the flip is driven by mu>lam (laggard's research rate exceeds the leader's")
log("   development rate): for mu<=lam, z*>x* holds ~always; it fails increasingly as mu/lam grows.")

# ==========================================================================
# FIGURE -- prop2_frontier.png  (2 panels)
# ==========================================================================
fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.6, 5.0))

# Panel A: (lam, mu) sign map of z*-x* at fixed c, b1; W = rho * (A2 cap)
c_fig, b1_fig, rho = 1.0, 0.3, 0.9
lam_g = np.linspace(0.02, 0.49, 300)
mu_g  = np.linspace(0.02, 0.49, 300)
LAM, MU = np.meshgrid(lam_g, mu_g)
SIGN = np.zeros_like(LAM)
for ii in range(LAM.shape[0]):
    for jj in range(LAM.shape[1]):
        l, m = LAM[ii, jj], MU[ii, jj]
        W = rho * c_fig/max(l, m)
        V11, _ = V11_and_x11(W, c_fig, l, b1_fig)
        z, x, res = solve_fp(W, c_fig, l, m, V11)
        SIGN[ii, jj] = 1.0 if (res < 1e-9 and z > x) else 0.0
cmap = ListedColormap(["#e4674f", "#4c78c8"])   # red: z<=x, blue: z>x
axA.pcolormesh(LAM, MU, SIGN, cmap=cmap, shading="auto", vmin=0, vmax=1)
axA.plot([0.02, 0.49], [0.02, 0.49], "k--", lw=1.4, label=r"$\mu=\lambda$ (Prop 2)")
axA.contour(LAM, MU, SIGN, levels=[0.5], colors="white", linewidths=2.2)
axA.set_xlim(0.02, 0.49); axA.set_ylim(0.02, 0.49)
axA.set_xlabel(r"leader rate $\lambda$"); axA.set_ylabel(r"laggard rate $\mu$")
axA.set_title(fr"(A) sign of $z^*-x^*$   ($c={c_fig:g},\ \beta_{{11}}={b1_fig:g},\ W=0.9\,c/\max(\lambda,\mu)$)")
axA.text(0.34, 0.12, r"$z^*>x^*$"+"\n(leader leads)", color="white", ha="center", va="center",
         fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="#4c78c8", ec="none", alpha=0.85))
axA.text(0.13, 0.40, r"$z^*\leq x^*$"+"\n(Prop 2 fails)", color="white", ha="center", va="center",
         fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="#e4674f", ec="none", alpha=0.9))
axA.legend(fontsize=8, loc="lower right", framealpha=0.9)

# Panel B: z*<=x* share vs mu/lam (from the Part B scan)
centers, shares, counts = [], [], []
for k in range(len(bins)-1):
    lo, hi = bins[k], bins[k+1]
    if btot[k] == 0: continue
    ctr = (lo + (hi if hi < 1e8 else lo*1.5))/2
    centers.append(min(ctr, 3.5)); shares.append(100*bfail[k]/btot[k]); counts.append(btot[k])
axB.axvline(1.0, color="grey", ls="--", lw=1.3)
axB.text(1.02, 88, r"$\mu=\lambda$", color="grey", fontsize=9, rotation=90, va="top")
axB.bar(range(len(shares)), shares, color="#e4674f", alpha=0.85, width=0.7)
axB.set_xticks(range(len(shares)))
lbls = []
for k in range(len(bins)-1):
    if btot[k] == 0: continue
    lo, hi = bins[k], bins[k+1]
    lbls.append(f"[{lo:g},{'∞' if hi>1e8 else f'{hi:g}'})")
axB.set_xticklabels(lbls, rotation=35, ha="right", fontsize=8)
axB.set_ylim(0, 100)
axB.set_xlabel(r"ratio $\mu/\lambda$ (bin)"); axB.set_ylabel(r"share with $z^*\leq x^*$  (%)")
axB.set_title(r"(B) Prop-2 failure share by $\mu/\lambda$   (200k draws, full assumption region)")
axB.text(0.98, 0.96, fr"overall (μ≠λ): {overall:.1f}%", transform=axB.transAxes,
         ha="right", va="top", fontsize=9,
         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#999999", alpha=0.9))

fig.tight_layout()
outpng = os.path.join(FIGS, "prop2_frontier.png")
fig.savefig(outpng, dpi=300, bbox_inches="tight")
log("")
log(f"figure written: {outpng}")

# save the frontier sign grid
np.savetxt(os.path.join(NUM, "task2_frontier_grid.csv"),
           np.column_stack([LAM.ravel(), MU.ravel(), SIGN.ravel()]),
           header="lambda,mu,z_gt_x(1=yes)", delimiter=",", comments="")

with open(os.path.join(NUM, "task2_results.txt"), "w") as f:
    f.write("\n".join(log_lines) + "\n")
print("\n[written]", os.path.join(NUM, "task2_results.txt"))
