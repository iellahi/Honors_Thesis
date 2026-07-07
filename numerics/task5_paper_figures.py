"""
Task 5 -- regenerate paper figures.
  1. beta11_sweep.png       : frontier-effort vs competitive-balance trade-off
                              (x11* falls, laggard x* rises, V11 rises as beta11 up).
  2. beta00_sweep.png       : (0,0) response to beta00 (x00*, P00=m x00*, V00).
  3. c_lambda_V11_sweep.png : V11 increasing in c, decreasing in lambda
                              (contrast: continuous-time model had V11 independent of c, lambda).
Already produced (referenced in the inventory): prop2_frontier.png (Task 2), laggard_BR.png (Task 4).

Solvers ported from verify_step5b.py.  One consistent baseline across sweeps.
Baseline: c=1, lam=mu=0.30, W=2.5, beta11=beta00=0.30   (satisfies A2 & A2').
Figures: PNG 300 dpi.  Outputs: figs/*.png, numerics/task5_results.txt
"""
import math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE   = os.path.dirname(os.path.abspath(__file__))
THESIS = os.path.dirname(HERE)
FIGS   = os.path.join(THESIS, "figs")
NUM    = os.path.join(THESIS, "numerics")
plt.rcParams.update({"axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 8.5})
log_lines = []
def log(s=""):
    print(s); log_lines.append(str(s))

# ---------------- solvers (verify_step5b.py) ----------------
def x11_and_V11(W, c, lam, b1):
    u = lam*lam*W*(1-b1*b1); Delta = 9*c*c - 2*c*u + u*u
    x11 = (3*c + u - math.sqrt(Delta))/(2*c*lam*(1+b1))
    V11 = c*x11/(2*lam*(1+b1)) + b1*W/(1+b1)
    return x11, V11
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
    return z, x
def values_10_01(z, x, W, c, lam, mu, V11):
    PL, PF = lam*z, mu*x; S = PL + PF - PL*PF
    V10 = (PL*W + (1-PL)*PF*V11 - c*z*z/2)/S
    V01 = ((1-PL)*PF*V11 - c*x*x/2)/S
    return V10, V01
def x00_V00(c, mu, b0, V10, V01, V11):
    m = mu*(1+b0); G = V10 + b0*V01 - (1+b0)*V11; D = V10 - V01
    A0, B0, C0 = c*m, -(3*c + 2*mu*m*G), 2*mu*(1-b0)*D
    x00 = (-B0 - math.sqrt(B0*B0 - 4*A0*C0))/(2*A0)
    P = m*x00
    V00 = (mu*(V10 + b0*V01) - m*P*(V10 + V01 - V11) - c*x00)/(m*(1-P))
    return x00, V00, P

# baseline
C0_, LAM0, MU0, W0, B11_0, B00_0 = 1.0, 0.30, 0.30, 2.5, 0.30, 0.30
BLUE, RED, GREEN, PURPLE = "#4c78c8", "#e45756", "#54a24b", "#b279a2"

# ==========================================================================
# FIGURE 1 -- beta11 sweep: the trade-off
# ==========================================================================
b11s = np.linspace(1e-3, 0.985, 300)
x11s, Vs, xstar, zstar = [], [], [], []
for b in b11s:
    x11, V11 = x11_and_V11(W0, C0_, LAM0, b)
    z, x = solve_fp(W0, C0_, LAM0, MU0, V11)
    x11s.append(x11); Vs.append(V11); xstar.append(x); zstar.append(z)
x11s, Vs, xstar, zstar = map(np.array, (x11s, Vs, xstar, zstar))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 5.0))
a1.plot(b11s, x11s, color=BLUE,  lw=2.4, label=r"development effort $x_{11}^*$ (frontier)")
a1.plot(b11s, xstar, color=GREEN, lw=2.4, label=r"laggard effort $x^*$ (catch-up)")
a1.plot(b11s, zstar, color=RED,   lw=2.0, ls="--", label=r"leader effort $z^*$")
a1.set_xlabel(r"development-phase spillover $\beta_{11}$"); a1.set_ylabel("effort")
a1.set_title(r"(A) frontier effort $\downarrow$, laggard catch-up $\uparrow$")
a1.legend(loc="best", framealpha=0.9); a1.set_xlim(0, 1)
a1.annotate("", xy=(0.85, x11s[int(0.85*len(b11s))]), xytext=(0.55, x11s[int(0.55*len(b11s))]),
            arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.2))
a1.annotate("", xy=(0.85, xstar[int(0.85*len(b11s))]), xytext=(0.55, xstar[int(0.55*len(b11s))]),
            arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.2))

a2.plot(b11s, Vs, color=PURPLE, lw=2.6, label=r"tied-state value $V_{11}$")
a2.axhline(W0/2, color="grey", ls=":", lw=1.3, label=r"$W/2$ (upper bound)")
a2.set_xlabel(r"development-phase spillover $\beta_{11}$"); a2.set_ylabel(r"$V_{11}$")
a2.set_title(r"(B) tied-state value $V_{11}$ rises (Prop 1$'$)")
a2.legend(loc="best", framealpha=0.9); a2.set_xlim(0, 1)
fig.suptitle(fr"$\beta_{{11}}$ sweep: frontier-effort vs competitive-balance trade-off "
             fr"($c={C0_:g},\ \lambda=\mu={LAM0:g},\ W={W0:g}$)", fontsize=11)
fig.tight_layout()
f1 = os.path.join(FIGS, "beta11_sweep.png"); fig.savefig(f1, dpi=300, bbox_inches="tight")
log(f"[beta11] x11*: {x11s[0]:.4f} -> {x11s[-1]:.4f} ({'down' if x11s[-1]<x11s[0] else 'up'});  "
    f"x*: {xstar[0]:.4f} -> {xstar[-1]:.4f} ({'up' if xstar[-1]>xstar[0] else 'down'});  "
    f"V11: {Vs[0]:.4f} -> {Vs[-1]:.4f} ({'up' if Vs[-1]>Vs[0] else 'down'})")
log(f"  monotone: x11* strictly down? {np.all(np.diff(x11s)<0)};  "
    f"x* strictly up? {np.all(np.diff(xstar)>0)};  V11 strictly up? {np.all(np.diff(Vs)>0)}")
log(f"figure: {f1}")

# ==========================================================================
# FIGURE 2 -- beta00 sweep: initial-state response
# ==========================================================================
# fix upstream values at baseline (they do not depend on beta00)
_x11b, V11b = x11_and_V11(W0, C0_, LAM0, B11_0)
zb, xb = solve_fp(W0, C0_, LAM0, MU0, V11b)
V10b, V01b = values_10_01(zb, xb, W0, C0_, LAM0, MU0, V11b)
b00s = np.linspace(0.0, 0.985, 300)
x00s, P00s, V00s = [], [], []
for b in b00s:
    x00, V00, P = x00_V00(C0_, MU0, b, V10b, V01b, V11b)
    x00s.append(x00); P00s.append(P); V00s.append(V00)
x00s, P00s, V00s = map(np.array, (x00s, P00s, V00s))

fig2, (b1, b2) = plt.subplots(1, 2, figsize=(12.4, 5.0))
b1.plot(b00s, x00s, color=BLUE, lw=2.4, label=r"initial effort $x_{00}^*$")
b1.plot(b00s, P00s, color=GREEN, lw=2.2, ls="--", label=r"breakthrough prob. $P_{00}=m x_{00}^*$")
b1.axhline(1.0, color="grey", ls=":", lw=1.0)
b1.set_xlabel(r"initial-phase spillover $\beta_{00}$"); b1.set_ylabel("effort / probability")
b1.set_title(r"(A) initial effort $x_{00}^*$ and $P_{00}$ vs $\beta_{00}$")
b1.legend(loc="best", framealpha=0.9); b1.set_xlim(0, 1)
b2.plot(b00s, V00s, color=PURPLE, lw=2.6, label=r"initial-state value $V_{00}$")
b2.set_xlabel(r"initial-phase spillover $\beta_{00}$"); b2.set_ylabel(r"$V_{00}$")
b2.set_title(r"(B) initial-state value $V_{00}$ vs $\beta_{00}$")
b2.legend(loc="best", framealpha=0.9); b2.set_xlim(0, 1)
fig2.suptitle(fr"$\beta_{{00}}$ sweep ($c={C0_:g},\ \lambda=\mu={LAM0:g},\ W={W0:g},\ \beta_{{11}}={B11_0:g}$)",
              fontsize=11)
fig2.tight_layout()
f2 = os.path.join(FIGS, "beta00_sweep.png"); fig2.savefig(f2, dpi=300, bbox_inches="tight")
log(f"[beta00] x00*: {x00s[0]:.4f} -> {x00s[-1]:.4f};  P00: {P00s[0]:.4f} -> {P00s[-1]:.4f};  "
    f"V00: {V00s[0]:.4f} -> {V00s[-1]:.4f}")
log(f"  x00* monotone decreasing? {np.all(np.diff(x00s)<0)}  (max P00={P00s.max():.4f}<1? {P00s.max()<1})")
log(f"figure: {f2}")

# ==========================================================================
# FIGURE 3 -- V11 vs c and vs lambda  (contrast with continuous-time model)
# ==========================================================================
#   Dependence is qualitatively new but modest; use beta11=0.1 (largest relative
#   effect), zoom each panel, overlay x11*, add a flat continuous-time counterfactual.
B11_F = 0.10
cs = np.linspace(0.76, 6.0, 300)                 # c >= lam W = 0.75
V11_c = np.array([x11_and_V11(W0, cc, LAM0, B11_F)[1] for cc in cs])
x11_c = np.array([x11_and_V11(W0, cc, LAM0, B11_F)[0] for cc in cs])
lams = np.linspace(0.02, 0.40, 300)              # lam <= c/W = 0.4
V11_l = np.array([x11_and_V11(W0, C0_, ll, B11_F)[1] for ll in lams])
x11_l = np.array([x11_and_V11(W0, C0_, ll, B11_F)[0] for ll in lams])

def _zoom(ax, y, frac=0.18):
    lo, hi = y.min(), y.max(); pad = (hi-lo)*frac; ax.set_ylim(lo-pad, hi+pad)

fig3, (c1, c2) = plt.subplots(1, 2, figsize=(12.4, 5.0))
# Panel A: V11 vs c (zoomed), with flat continuous-time counterfactual + x11* on twin axis
c1.plot(cs, V11_c, color=PURPLE, lw=2.8, label=r"$V_{11}$ (discrete-time)", zorder=3)
c1.axhline(V11_c[0], color="grey", ls="--", lw=1.5, label=r"continuous-time: $V_{11}\perp c$ (flat)")
_zoom(c1, V11_c); c1.set_xlabel(r"cost curvature $c$")
c1.set_ylabel(r"$V_{11}$", color=PURPLE); c1.tick_params(axis="y", labelcolor=PURPLE)
c1.set_title(r"(A) $V_{11}$ increasing in $c$")
c1.annotate(fr"+{100*(V11_c[-1]-V11_c[0])/V11_c[0]:.1f}% over range", xy=(0.52, 0.45),
            xycoords="axes fraction", ha="center", fontsize=9, color=PURPLE, style="italic")
c1x = c1.twinx(); c1x.plot(cs, x11_c, color=BLUE, lw=1.6, ls=":", label=r"$x_{11}^*$ (right axis)")
c1x.set_ylabel(r"$x_{11}^*$", color=BLUE); c1x.tick_params(axis="y", labelcolor=BLUE)
h, l = c1.get_legend_handles_labels(); hx, lx = c1x.get_legend_handles_labels()
c1.legend(h+hx, l+lx, loc="lower right", framealpha=0.9)
# Panel B: V11 vs lambda
c2.plot(lams, V11_l, color=PURPLE, lw=2.8, label=r"$V_{11}$ (discrete-time)", zorder=3)
c2.axhline(V11_l[0], color="grey", ls="--", lw=1.5, label=r"continuous-time: $V_{11}\perp\lambda$ (flat)")
_zoom(c2, V11_l); c2.set_xlabel(r"development rate $\lambda$")
c2.set_ylabel(r"$V_{11}$", color=PURPLE); c2.tick_params(axis="y", labelcolor=PURPLE)
c2.set_title(r"(B) $V_{11}$ decreasing in $\lambda$")
c2.annotate(fr"{100*(V11_l[-1]-V11_l[0])/V11_l[0]:.1f}% over range", xy=(0.52, 0.55),
            xycoords="axes fraction", ha="center", fontsize=9, color=PURPLE, style="italic")
c2x = c2.twinx(); c2x.plot(lams, x11_l, color=BLUE, lw=1.6, ls=":", label=r"$x_{11}^*$ (right axis)")
c2x.set_ylabel(r"$x_{11}^*$", color=BLUE); c2x.tick_params(axis="y", labelcolor=BLUE)
h, l = c2.get_legend_handles_labels(); hx, lx = c2x.get_legend_handles_labels()
c2.legend(h+hx, l+lx, loc="upper right", framealpha=0.9)
fig3.suptitle(fr"$V_{{11}}$ depends on $c$ and $\lambda$ (discrete-time) — qualitatively new vs "
              fr"continuous-time  ($W={W0:g},\ \beta_{{11}}={B11_F:g}$; modest magnitude)", fontsize=10.5)
fig3.tight_layout()
f3 = os.path.join(FIGS, "c_lambda_V11_sweep.png"); fig3.savefig(f3, dpi=300, bbox_inches="tight")
log(f"[c-sweep]   V11: {V11_c[0]:.4f} -> {V11_c[-1]:.4f} ({100*(V11_c[-1]-V11_c[0])/V11_c[0]:+.2f}%)  "
    f"strictly up? {np.all(np.diff(V11_c)>0)};  x11* down? {np.all(np.diff(x11_c)<0)}")
log(f"[lam-sweep] V11: {V11_l[0]:.4f} -> {V11_l[-1]:.4f} ({100*(V11_l[-1]-V11_l[0])/V11_l[0]:+.2f}%)  "
    f"strictly down? {np.all(np.diff(V11_l)<0)};  x11* up? {np.all(np.diff(x11_l)>0)}")
log(f"figure: {f3}  (beta11={B11_F} for visibility; monotone but modest)")

log("")
log("Referenced (already generated): figs/prop2_frontier.png (Task 2), figs/laggard_BR.png (Task 4)")
with open(os.path.join(NUM, "task5_results.txt"), "w") as f:
    f.write("\n".join(log_lines) + "\n")
print("\n[written]", os.path.join(NUM, "task5_results.txt"))
