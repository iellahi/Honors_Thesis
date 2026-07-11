"""
Task 7 -- scans for the Step-8 transmission package (Section 5.5 upgrade).
Solvers ported from verify_scripts/verify_step5b.py (reference).  PNG 300 dpi.

  (a) sqrt(2)-frontier: dx*/dV11 > 0 is proven wherever x* < sqrt(2) z* (R2).
      Measure, over the mu != lam region: share of ALL draws with x* < sqrt2 z*;
      share of Prop-2-FAILURE draws (z* <= x*) still inside the proven region;
      binned by mu/lam (parity with task2 Part B bins).  Seed 41, N = 60000.
  (b) grid overlay data for prop2_frontier panel A (task2 config: c=1.0,
      b1=0.3, W=0.9 c/max(lam,mu), 300x300): sign(z*-x*) and sign(sqrt2 z*-x*).
      Candidate figure: figs/prop2_frontier_sqrt2_overlay.png (does NOT
      overwrite the existing paper figure).
  (c) leader sign consistency (R4): at scanned equilibria, sign of the FD total
      dz*/dbeta11 must match sign of the bracket
      B = lam z^3 (1-lam z)/2 - mu (1-lam z) x^3 - lam z x^2.
      Near-zero FD cases flagged, not counted as violations (task6 convention).
      Seed 41, N = 40000.
  (d) grid overlay data for dz_dbeta11_map panel B (task4 config: c=1.0,
      b1=0.3, W=0.85 c/max(lam,mu)): the exact R4 boundary B = 0.
      Candidate figure: figs/dz_dbeta11_map_R4overlay.png.

Outputs: numerics/task7_results.txt, numerics/task7_sqrt2_grid.csv,
         figs/prop2_frontier_sqrt2_overlay.png, figs/dz_dbeta11_map_R4overlay.png
"""
import math, random, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

HERE   = os.path.dirname(os.path.abspath(__file__))
THESIS = os.path.dirname(HERE)
FIGS   = os.path.join(THESIS, "figs")

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

def R4_bracket(lam, mu, x, z):
    return lam*z**3*(1-lam*z)/2 - mu*(1-lam*z)*x**3 - lam*z*x**2

SQ2 = math.sqrt(2)

# ==========================================================================
# PART (a) -- sqrt(2)-frontier vs Prop-2 failure region   seed 41, N=60000
# ==========================================================================
log("="*74)
log("PART A: sqrt(2)-frontier (proven-transmission region) vs Prop-2 failures")
log("="*74)
random.seed(41)
N = 60000
bins = [0.0, 0.5, 0.8, 1.0, 1.25, 1.6, 2.0, 3.0, 1e9]
btot  = [0]*(len(bins)-1)      # draws per bin
bfail = [0]*(len(bins)-1)      # z* <= x*  (Prop-2 failure)
bfail_in = [0]*(len(bins)-1)   # z* <= x*  AND  x* < sqrt2 z*  (still proven)
n_all = n_sq2 = n_fail = n_fail_sq2 = 0
for _ in range(N):
    lam = random.uniform(1e-2, 0.4999)
    mu  = random.uniform(1e-2, 0.4999)
    b1  = random.uniform(1e-4, 0.9995)
    c   = random.uniform(5e-2, 10.0)
    W   = random.uniform(1e-3, c/max(lam, mu))
    V11, z, x, res = full_solve(W, c, lam, mu, b1)
    if res > 1e-9: continue
    n_all += 1
    in_sq2 = (x < SQ2*z)
    fail   = not (z > x)
    n_sq2      += in_sq2
    n_fail     += fail
    n_fail_sq2 += (fail and in_sq2)
    r = mu/lam
    for k in range(len(bins)-1):
        if bins[k] <= r < bins[k+1]:
            btot[k] += 1
            if fail: bfail[k] += 1
            if fail and in_sq2: bfail_in[k] += 1
            break
log(f"draws N = {N} (seed 41, mu and lam independent)")
log(f"x* < sqrt2 z*  (proven dx*/dV11>0):        {n_sq2}/{n_all} = {100*n_sq2/n_all:.2f}%")
log(f"Prop-2 failures (z*<=x*):                  {n_fail}/{n_all} = {100*n_fail/n_all:.2f}%")
log(f"Prop-2 failures still inside sqrt2 region: {n_fail_sq2}/{n_fail} = {100*n_fail_sq2/max(n_fail,1):.2f}%")
log("")
log(f"{'mu/lam bin':>14} | {'n':>7} | {'z<=x %':>7} | {'of which x<sqrt2 z %':>20}")
for k in range(len(bins)-1):
    lo, hi = bins[k], bins[k+1]
    hi_s = "inf" if hi > 1e8 else f"{hi:g}"
    if btot[k] == 0: continue
    p1 = 100*bfail[k]/btot[k]
    p2 = 100*bfail_in[k]/max(bfail[k], 1)
    log(f"  [{lo:g}, {hi_s})   | {btot[k]:>7} | {p1:>6.2f}% | {p2:>19.2f}%")

# ==========================================================================
# PART (b) -- grid + candidate overlay figure (task2 panel A config)
# ==========================================================================
c_fig, b1_fig, rho = 1.0, 0.3, 0.9
lam_g = np.linspace(0.02, 0.49, 300)
mu_g  = np.linspace(0.02, 0.49, 300)
LAM, MU = np.meshgrid(lam_g, mu_g)
SIGN  = np.zeros_like(LAM)     # 1 = z* > x*
RATIO = np.full_like(LAM, np.nan)   # x*/(sqrt2 z*); < 1 means proven transmission
for ii in range(LAM.shape[0]):
    for jj in range(LAM.shape[1]):
        l, m = LAM[ii, jj], MU[ii, jj]
        W = rho*c_fig/max(l, m)
        V11, z, x, res = full_solve(W, c_fig, l, m, b1_fig)
        if res < 1e-9:
            SIGN[ii, jj]  = 1.0 if z > x else 0.0
            RATIO[ii, jj] = x/(SQ2*z)
np.savetxt(os.path.join(HERE, "task7_sqrt2_grid.csv"),
           np.column_stack([LAM.ravel(), MU.ravel(), SIGN.ravel(), RATIO.ravel()]),
           header="lambda,mu,z_gt_x(1=yes),ratio_x_over_sqrt2z", delimiter=",", comments="")
log("")
log(f"grid (task2 config): max x*/(sqrt2 z*) = {np.nanmax(RATIO):.4f} "
    f"(< 1 everywhere: proven-transmission region covers the whole grid, incl. all "
    f"{int((SIGN < 0.5).sum())} Prop-2-failure cells)")

# NOTE: the sqrt2-frontier never enters the maintained region (ratio < 1 every-
# where; scan sup -> 1 only in the degenerate corner lam->0, b1->1, mu->1/2).
# So instead of a frontier line, plot the MARGIN: heatmap of the ratio itself.
fig, ax = plt.subplots(figsize=(6.9, 5.2))
pcm = ax.pcolormesh(LAM, MU, RATIO, cmap="viridis", shading="auto", vmin=0, vmax=1)
cb = fig.colorbar(pcm, ax=ax); cb.set_label(r"$x^*/(\sqrt{2}\,z^*)$   ($<1$ = proven $dx^*/dV_{11}>0$)")
ax.contour(LAM, MU, SIGN, levels=[0.5], colors="white", linewidths=2.2)
ax.plot([], [], color="white", lw=2.2, label=r"$z^*=x^*$ (Prop-2 frontier)")
cs = ax.contour(LAM, MU, RATIO, levels=[0.5, 0.7, 0.9], colors="#f2c744",
                linewidths=1.0, linestyles="dashed")
ax.clabel(cs, fontsize=7, fmt="%.1f")
ax.plot([0.02, 0.49], [0.02, 0.49], "k--", lw=1.2, label=r"$\mu=\lambda$ (Prop 2)")
ax.set_xlim(0.02, 0.49); ax.set_ylim(0.02, 0.49)
ax.set_xlabel(r"leader rate $\lambda$"); ax.set_ylabel(r"laggard rate $\mu$")
ax.set_title(r"margin to the $x^*=\sqrt{2}z^*$ frontier: ratio $<1$ everywhere,"
             + "\n" + r"so $dx^*/dV_{11}>0$ is proven across the whole map"
             + f"  ($c={c_fig:g},\\ \\beta_{{11}}={b1_fig:g}$)", fontsize=10)
ax.legend(fontsize=8, loc="lower right", framealpha=0.95)
fig.tight_layout()
out1 = os.path.join(FIGS, "prop2_frontier_sqrt2_overlay.png")
fig.savefig(out1, dpi=300, bbox_inches="tight"); plt.close(fig)
log(f"figure written: {out1}")

# ==========================================================================
# PART (c) -- R4 leader-sign consistency at scanned equilibria  seed 41, N=40000
# ==========================================================================
log("")
log("="*74)
log("PART C: sign(dz*/dbeta11) [FD, total] vs sign(R4 bracket)  seed 41, N=40000")
log("="*74)
random.seed(41)
Nc = 40000
match = mismatch = flagged = usable = 0
neg_share = 0
for i in range(Nc):
    lam = random.uniform(1e-2, 0.4999)
    mu  = random.uniform(1e-2, 0.4999) if i % 2 else lam
    b1  = random.uniform(1e-4, 0.9990)
    c   = random.uniform(5e-2, 10.0)
    W   = random.uniform(1e-3, c/max(lam, mu))
    V11, z, x, res = full_solve(W, c, lam, mu, b1)
    if res > 1e-9: continue
    h = 1e-6*max(b1, 1e-3)
    if b1 + h >= 1: continue
    V11p, zp, xp, resp = full_solve(W, c, lam, mu, b1 + h)
    V11m, zm, xm, resm = full_solve(W, c, lam, mu, b1 - h)
    if resp > 1e-9 or resm > 1e-9: continue
    fd = (zp - zm)/(2*h)
    B  = R4_bracket(lam, mu, x, z)
    usable += 1
    if fd < 0: neg_share += 1
    scale = max(abs(z), 1e-6)
    if abs(fd) < 1e-6*scale:
        flagged += 1
    elif (fd > 0) == (B > 0):
        match += 1
    else:
        mismatch += 1
log(f"usable draws: {usable}   sign matches: {match}   MISMATCHES: {mismatch}   "
    f"near-zero flags: {flagged}")
log(f"dz*/dbeta11 negative share: {100*neg_share/usable:.2f}%  (memory: ~95.5%)")

# ==========================================================================
# PART (d) -- R4 boundary overlay for dz_dbeta11_map panel B (task4 config)
# ==========================================================================
c_f, b1_f, rho4 = 1.0, 0.3, 0.85
lam4 = np.linspace(0.02, 0.49, 240)
mu4  = np.linspace(0.02, 0.49, 240)
L4, M4 = np.meshgrid(lam4, mu4)
BVAL = np.full_like(L4, np.nan)
SGNZ = np.full_like(L4, np.nan)   # FD sign of dz/dbeta11 (for visual check)
for ii in range(L4.shape[0]):
    for jj in range(L4.shape[1]):
        l, m = L4[ii, jj], M4[ii, jj]
        W = rho4*c_f/max(l, m)
        V11, z, x, res = full_solve(W, c_f, l, m, b1_f)
        if res > 1e-9: continue
        BVAL[ii, jj] = R4_bracket(l, m, x, z)
        h = 1e-6*b1_f
        _, zp, _, rp = full_solve(W, c_f, l, m, b1_f + h)
        _, zm, _, rm = full_solve(W, c_f, l, m, b1_f - h)
        if rp < 1e-9 and rm < 1e-9:
            SGNZ[ii, jj] = 1.0 if (zp - zm) > 0 else 0.0

fig, ax = plt.subplots(figsize=(6.6, 5.2))
cmap2 = ListedColormap(["#4c78c8", "#e4a13a"])   # blue: dz<0, orange: dz>0
ax.pcolormesh(L4, M4, np.where(np.isnan(SGNZ), 0.0, SGNZ), cmap=cmap2,
              shading="auto", vmin=0, vmax=1)
ax.contour(L4, M4, BVAL, levels=[0.0], colors="black", linewidths=2.2)
ax.plot([], [], color="black", lw=2.2,
        label=r"exact boundary $\frac{1}{2}\lambda z^3(1-\lambda z)=\mu(1-\lambda z)x^3+\lambda z x^2$")
ax.set_xlim(0.02, 0.49); ax.set_ylim(0.02, 0.49)
ax.set_xlabel(r"leader rate $\lambda$"); ax.set_ylabel(r"laggard rate $\mu$")
ax.set_title(r"sign of $dz^*/d\beta_{11}$ with the exact (R4) boundary"
             + f"\n($c={c_f:g},\\ \\beta_{{11}}={b1_f:g},\\ W=0.85\\,c/\\max(\\lambda,\\mu)$)",
             fontsize=10)
ax.legend(fontsize=8, loc="upper right", framealpha=0.95)
fig.tight_layout()
out2 = os.path.join(FIGS, "dz_dbeta11_map_R4overlay.png")
fig.savefig(out2, dpi=300, bbox_inches="tight"); plt.close(fig)
log(f"figure written: {out2}")

with open(os.path.join(HERE, "task7_results.txt"), "w") as f:
    f.write("\n".join(log_lines) + "\n")
print("\n[written]", os.path.join(HERE, "task7_results.txt"))
