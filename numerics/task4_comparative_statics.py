"""
Task 4 -- comparative statics within the assumption region.
  (a) dx*/dbeta11 > 0  (laggard catch-up; direct channel analytical, full sign numerical) -> share +.
  (b) dz*/dbeta11 sign map (ambiguous, ~96% negative) -> (lam,mu) sign map + overall share.
  (c) dx*/dlam total vs direct (V11 frozen): rare sign reversals -> map + share (Luke's caveat).
  (d) non-monotone laggard BR X(z): provocation [mu x(1-lam z)^2 > 2 lam^2 z^2] vs discouragement;
      ~99.8% of equilibria in the provocation region.  Also det J > 0 share (regularity).

Finite differences on the full closed system. Solvers ported from verify_step5b.py / verify_step5d.py.
Figures: PNG 300 dpi. Outputs:
  numerics/task4_results.txt, figs/dz_dbeta11_map.png, figs/dx_dlam_reversals.png, figs/laggard_BR.png
"""
import math, random, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap

HERE   = os.path.dirname(os.path.abspath(__file__))
THESIS = os.path.dirname(HERE)
FIGS   = os.path.join(THESIS, "figs")
NUM    = os.path.join(THESIS, "numerics")
log_lines = []
def log(s=""):
    print(s); log_lines.append(str(s))

# ---------------- solvers (verify_step5b.py) ----------------
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
def full_solve(W, c, lam, mu, b1):
    V11 = V11_and_x11(W, c, lam, b1)
    z, x, res = solve_fp(W, c, lam, mu, V11)
    return V11, z, x, res

# ==========================================================================
# PART A -- transmission wrt beta11 (total, finite difference)   seed 11
# ==========================================================================
log("="*74)
log("PART A: transmission d/dbeta11  (seed 11, N=40000)")
log("="*74)
random.seed(11)
N = 40000
dx_pos = dx_tot = 0
dz_pos = dz_neg = 0
for i in range(N):
    lam = random.uniform(1e-2, 0.4999)
    mu  = random.uniform(1e-2, 0.4999) if i % 2 else lam
    b1  = random.uniform(1e-4, 0.9990)
    c   = random.uniform(5e-2, 10.0)
    W   = random.uniform(1e-3, c/max(lam, mu))
    V11, z, x, res = full_solve(W, c, lam, mu, b1)
    if res > 1e-9: continue
    h = 1e-6*max(b1, 1e-3)
    if b1 + h >= 1: continue
    _, zb, xb, resb = full_solve(W, c, lam, mu, b1 + h)
    if resb > 1e-9: continue
    dx_tot += 1
    if (xb - x)/h > 0: dx_pos += 1
    if (zb - z)/h > 0: dz_pos += 1
    else:              dz_neg += 1
log(f"[a] dx*/dbeta11 > 0 : {dx_pos}/{dx_tot}  = {100*dx_pos/dx_tot:.3f}%  (laggard catch-up; expect ~100%)")
log(f"[b] dz*/dbeta11 sign: positive {dz_pos}, negative {dz_neg}  "
    f"-> negative {100*dz_neg/(dz_pos+dz_neg):.2f}%  (memory ~96%)")

# ==========================================================================
# PART B -- dx*/dlam total vs direct (V11 frozen): rare sign reversals  seed 11
# ==========================================================================
log("")
log("="*74)
log("PART B: total vs direct effects for c and lam (V11 frozen = direct)  (subsample)")
log("="*74)
random.seed(11)
tot = {k: [0, 0] for k in ["dz_dc", "dx_dc", "dz_dlam", "dx_dlam"]}
flip = {k: 0 for k in ["dz_dc", "dx_dc", "dz_dlam", "dx_dlam"]}
nsub = 0
for i in range(N):
    lam = random.uniform(1e-2, 0.4999)
    mu  = random.uniform(1e-2, 0.4999) if i % 2 else lam
    b1  = random.uniform(1e-4, 0.9990)
    c   = random.uniform(5e-2, 10.0)
    W   = random.uniform(1e-3, c/max(lam, mu))
    if i % 5 != 0: continue
    V11, z, x, res = full_solve(W, c, lam, mu, b1)
    if res > 1e-9: continue
    nsub += 1
    for par, kz, kx in [("c", "dz_dc", "dx_dc"), ("lam", "dz_dlam", "dx_dlam")]:
        if par == "c":
            h = 1e-6*c
            _, zp, xp, rp = full_solve(W, c+h, lam, mu, b1)
            zd, xd, rd = solve_fp(W, c+h, lam, mu, V11)          # direct: V11 frozen
        else:
            h = 1e-6*lam
            if lam + h >= 0.5: continue
            _, zp, xp, rp = full_solve(W, c, lam+h, mu, b1)
            zd, xd, rd = solve_fp(W, c, lam+h, mu, V11)
        if rp > 1e-9 or rd > 1e-9: continue
        tz, tx = (zp-z)/h, (xp-x)/h
        dzc, dxc = (zd-z)/h, (xd-x)/h                            # direct
        tot[kz][0 if tz < 0 else 1] += 1
        tot[kx][0 if tx < 0 else 1] += 1
        if (tz < 0) != (dzc < 0): flip[kz] += 1
        if (tx < 0) != (dxc < 0): flip[kx] += 1
log(f"subsample n = {nsub}   ([neg,pos] counts; flips = direct-vs-total sign disagreements)")
for k in ["dz_dc", "dx_dc", "dz_dlam", "dx_dlam"]:
    n = sum(tot[k])
    log(f"  {k:8}: [neg,pos]={tot[k]}   total-negative {100*tot[k][0]/max(n,1):.2f}%   "
        f"sign reversals vs direct: {flip[k]}  ({100*flip[k]/max(n,1):.3f}%)")
log("=> dx*/dlam is the one channel with sign reversals (indirect V11 effect flips the direct sign).")

# ==========================================================================
# PART C -- provocation region & regularity  seed 21   (verify_step5d parity)
# ==========================================================================
log("")
log("="*74)
log("PART C: laggard-BR slope region & det J  (seed 21, N=40000)")
log("="*74)
def jac_det(z, x, W, c, lam, mu, V11, h=1e-7):
    def GL(z_, x_): return (c*lam/2)*(1-mu*x_)*z_*z_ + c*mu*x_*z_ - lam*mu*x_*(W-V11)
    def GF(z_, x_): return (c*mu/2)*(1-lam*z_)*x_*x_ + c*lam*z_*x_ - mu*lam*z_*(1-lam*z_)*V11
    a11 = (GL(z+h,x)-GL(z-h,x))/(2*h); a12 = (GL(z,x+h)-GL(z,x-h))/(2*h)
    a21 = (GF(z+h,x)-GF(z-h,x))/(2*h); a22 = (GF(z,x+h)-GF(z,x-h))/(2*h)
    return a11*a22 - a12*a21
random.seed(21)
Nc = 40000
prov = disc = 0; badJ = 0; margins = []
for _ in range(Nc):
    lam = random.uniform(1e-2, 0.4999); mu = random.uniform(1e-2, 0.4999)
    c   = random.uniform(5e-2, 10.0);   W  = random.uniform(1e-3, c/max(lam, mu))
    b1  = random.uniform(1e-4, 0.9990)
    V11, z, x, res = full_solve(W, c, lam, mu, b1)
    if res > 1e-9: continue
    N_margin = mu*x*(1-lam*z)**2 - 2*lam*lam*z*z    # >0 provocation, <0 discouragement
    margins.append(N_margin)
    if N_margin > 0: prov += 1
    else:            disc += 1
    if jac_det(z, x, W, c, lam, mu, V11) <= 0: badJ += 1
tot_c = prov + disc
log(f"equilibria in PROVOCATION region (N>0): {prov}/{tot_c} = {100*prov/tot_c:.2f}%  (memory ~99.8%)")
log(f"equilibria in DISCOURAGEMENT region (N<=0): {disc}/{tot_c} = {100*disc/tot_c:.2f}%")
log(f"det J <= 0 violations: {badJ}/{tot_c}  (regularity/uniqueness support; expect 0)")

# ==========================================================================
# FIGURE 1 -- dz_dbeta11_map.png : (a) dx/dbeta11>0 vs (b) dz/dbeta11 sign
# ==========================================================================
c_f, b1_f, rho = 1.0, 0.3, 0.85
lam_g = np.linspace(0.02, 0.49, 90); mu_g = np.linspace(0.02, 0.49, 90)
LAM, MU = np.meshgrid(lam_g, mu_g)
SGN_dx = np.zeros_like(LAM); SGN_dz = np.zeros_like(LAM)
for ii in range(LAM.shape[0]):
    for jj in range(LAM.shape[1]):
        l, m = LAM[ii,jj], MU[ii,jj]; W = rho*c_f/max(l, m)
        V11, z, x, res = full_solve(W, c_f, l, m, b1_f)
        h = 1e-6
        _, zb, xb, resb = full_solve(W, c_f, l, m, b1_f+h)
        if res > 1e-9 or resb > 1e-9:
            SGN_dx[ii,jj] = np.nan; SGN_dz[ii,jj] = np.nan; continue
        SGN_dx[ii,jj] = 1.0 if (xb-x) > 0 else 0.0
        SGN_dz[ii,jj] = 1.0 if (zb-z) > 0 else 0.0
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.6, 5.0))
cmap2 = ListedColormap(["#e4674f", "#4c78c8"])
a1.pcolormesh(LAM, MU, SGN_dx, cmap=ListedColormap(["#e4674f", "#54a24b"]), vmin=0, vmax=1, shading="auto")
a1.plot([0.02,0.49],[0.02,0.49],"k--",lw=1.1)
a1.set_title(r"(A) sign of $\mathrm{d}x^*/\mathrm{d}\beta_{11}$  (laggard) — uniformly $>0$")
a1.set_xlabel(r"$\lambda$"); a1.set_ylabel(r"$\mu$")
a1.text(0.28,0.12,r"$\mathrm{d}x^*/\mathrm{d}\beta_{11}>0$"+"\n(catch-up, robust)", color="white",
        ha="center", fontsize=9.5, bbox=dict(boxstyle="round,pad=0.3", fc="#54a24b", ec="none", alpha=0.9))
a2.pcolormesh(LAM, MU, SGN_dz, cmap=cmap2, vmin=0, vmax=1, shading="auto")
a2.contour(LAM, MU, np.nan_to_num(SGN_dz, nan=0.0), levels=[0.5], colors="white", linewidths=2.0)
a2.plot([0.02,0.49],[0.02,0.49],"k--",lw=1.1, label=r"$\mu=\lambda$")
a2.set_title(r"(B) sign of $\mathrm{d}z^*/\mathrm{d}\beta_{11}$  (leader) — mostly $<0$")
a2.set_xlabel(r"$\lambda$"); a2.set_ylabel(r"$\mu$")
a2.legend(fontsize=8, loc="lower right", framealpha=0.9)
prox = [Line2D([0],[0],marker="s",color="none",markerfacecolor="#4c78c8",markersize=11,label=r"$>0$"),
        Line2D([0],[0],marker="s",color="none",markerfacecolor="#e4674f",markersize=11,label=r"$<0$")]
a2.legend(handles=prox+[Line2D([0],[0],ls="--",color="k",label=r"$\mu=\lambda$")],
          fontsize=8, loc="lower right", framealpha=0.9)
fig.suptitle(fr"transmission of $\beta_{{11}}$  ($c={c_f:g},\ \beta_{{11}}={b1_f:g},\ W=0.85\,c/\max(\lambda,\mu)$)",
             fontsize=11)
fig.tight_layout()
f1 = os.path.join(FIGS, "dz_dbeta11_map.png"); fig.savefig(f1, dpi=300, bbox_inches="tight")
log(""); log(f"figure written: {f1}")

# ==========================================================================
# FIGURE 2 -- dx_dlam_reversals.png : where total and direct dx*/dlam disagree
# ==========================================================================
SGN_flip = np.zeros_like(LAM); TX = np.zeros_like(LAM)
for ii in range(LAM.shape[0]):
    for jj in range(LAM.shape[1]):
        l, m = LAM[ii,jj], MU[ii,jj]; W = rho*c_f/max(l, m)
        V11, z, x, res = full_solve(W, c_f, l, m, b1_f)
        hl = 1e-6*l
        if l+hl >= 0.5 or res > 1e-9:
            SGN_flip[ii,jj] = np.nan; TX[ii,jj]=np.nan; continue
        _, zp, xp, rp = full_solve(W, c_f, l+hl, m, b1_f)
        zd, xd, rd = solve_fp(W, c_f, l+hl, m, V11)
        if rp > 1e-9 or rd > 1e-9:
            SGN_flip[ii,jj]=np.nan; TX[ii,jj]=np.nan; continue
        tx = (xp-x)/hl; dxd = (xd-x)/hl
        TX[ii,jj] = tx
        SGN_flip[ii,jj] = 1.0 if ((tx<0)!=(dxd<0)) else 0.0
n_flip = int(np.nansum(SGN_flip)); n_cells = int(np.sum(~np.isnan(SGN_flip)))
fig2, (b1, b2) = plt.subplots(1, 2, figsize=(12.6, 5.0))
pcm = b1.pcolormesh(LAM, MU, TX, cmap="RdBu_r", shading="auto",
                    vmin=-np.nanmax(np.abs(TX)), vmax=np.nanmax(np.abs(TX)))
fig2.colorbar(pcm, ax=b1, label=r"total $\mathrm{d}x^*/\mathrm{d}\lambda$")
b1.plot([0.02,0.49],[0.02,0.49],"k--",lw=1.0)
b1.set_title(r"(A) total $\mathrm{d}x^*/\mathrm{d}\lambda$ (mostly $>0$)")
b1.set_xlabel(r"$\lambda$"); b1.set_ylabel(r"$\mu$")
# reversal cells highlighted
b2.pcolormesh(LAM, MU, np.where(np.isnan(SGN_flip), 0.0, SGN_flip),
              cmap=ListedColormap(["#eeeeee", "#d62728"]), vmin=0, vmax=1, shading="auto")
b2.plot([0.02,0.49],[0.02,0.49],"k--",lw=1.0, label=r"$\mu=\lambda$")
b2.set_title(fr"(B) sign reversals total vs direct  ({n_flip}/{n_cells} cells on this slice)")
b2.set_xlabel(r"$\lambda$"); b2.set_ylabel(r"$\mu$")
b2.legend(handles=[Line2D([0],[0],marker="s",color="none",markerfacecolor="#d62728",markersize=11,
                          label="reversal"), Line2D([0],[0],ls="--",color="k",label=r"$\mu=\lambda$")],
          fontsize=8, loc="lower right", framealpha=0.9)
fig2.suptitle(r"$\mathrm{d}x^*/\mathrm{d}\lambda$: the indirect $V_{11}$ channel flips the direct sign only in a thin band",
              fontsize=10.5)
fig2.tight_layout()
f2 = os.path.join(FIGS, "dx_dlam_reversals.png"); fig2.savefig(f2, dpi=300, bbox_inches="tight")
log(f"figure written: {f2}   (slice reversal cells: {n_flip}/{n_cells})")

# ==========================================================================
# FIGURE 3 -- laggard_BR.png : X(z) provocation vs discouragement + eq; margin hist
# ==========================================================================
lam_r, mu_r, c_r, b1_r = 0.30, 0.30, 1.0, 0.30
W_r = 0.85*c_r/max(lam_r, mu_r)
V11_r = V11_and_x11(W_r, c_r, lam_r, b1_r)
zeq, xeq, _ = solve_fp(W_r, c_r, lam_r, mu_r, V11_r)
zs = np.linspace(1e-3, 0.985/lam_r, 700)
Xz = np.array([br_x(z, W_r, c_r, lam_r, mu_r, V11_r) for z in zs])
Nz = mu_r*Xz*(1-lam_r*zs)**2 - 2*lam_r**2*zs**2          # slope-sign of X at z
fig3, (c1, c2) = plt.subplots(1, 2, figsize=(12.6, 5.0))
prov_mask = Nz > 0
c1.plot(zs[prov_mask], Xz[prov_mask], ".", ms=2.4, color="#54a24b", label="provocation ($X'>0$)")
c1.plot(zs[~prov_mask], Xz[~prov_mask], ".", ms=2.4, color="#e4674f", label="discouragement ($X'<0$)")
# leader BR z=Z(x): plot as (Z(x), x)
xs2 = np.linspace(1e-3, 0.985/mu_r, 700)
Zx = np.array([br_z(x, W_r, c_r, lam_r, mu_r, V11_r) for x in xs2])
c1.plot(Zx, xs2, "-", color="#4c78c8", lw=1.8, label=r"leader BR $z=Z(x)$")
c1.plot([zeq], [xeq], "k*", ms=15, label=fr"equilibrium $(z^*,x^*)$")
c1.axvline(zeq, color="grey", ls=":", lw=1.0)
c1.set_xlim(0, min(1.0, 0.985/lam_r)); c1.set_ylim(0, max(Xz.max(), xeq)*1.1)
c1.set_xlabel(r"leader effort $z$"); c1.set_ylabel(r"laggard effort $x$")
c1.set_title(fr"(A) laggard BR $X(z)$  ($\lambda=\mu={lam_r:g},\ c={c_r:g},\ \beta_{{11}}={b1_r:g}$)")
c1.legend(fontsize=8, loc="upper right", framealpha=0.9)
# Panel B: histogram of equilibrium slope margin N across scan
margins = np.array(margins)
c2.hist(np.clip(margins, -0.05, 0.4), bins=80, color="#54a24b", alpha=0.85)
c2.axvline(0, color="red", ls="--", lw=1.6)
c2.set_xlabel(r"equilibrium margin $N=\mu x(1-\lambda z)^2-2\lambda^2z^2$")
c2.set_ylabel("count")
c2.set_title(fr"(B) {100*prov/tot_c:.2f}% of equilibria have $N>0$ (provocation)")
c2.text(0.97,0.95, fr"$N>0$: {prov}/{tot_c}"+"\n"+fr"$N\leq0$: {disc}/{tot_c}",
        transform=c2.transAxes, ha="right", va="top", fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#999", alpha=0.9))
fig3.tight_layout()
f3 = os.path.join(FIGS, "laggard_BR.png"); fig3.savefig(f3, dpi=300, bbox_inches="tight")
log(f"figure written: {f3}")
log(f"  representative eq: z*={zeq:.4f}, x*={xeq:.4f}, margin N={mu_r*xeq*(1-lam_r*zeq)**2-2*lam_r**2*zeq**2:.4f} (>0 provocation)")

with open(os.path.join(NUM, "task4_results.txt"), "w") as f:
    f.write("\n".join(log_lines) + "\n")
print("\n[written]", os.path.join(NUM, "task4_results.txt"))
