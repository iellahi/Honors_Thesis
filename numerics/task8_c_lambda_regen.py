"""
Task 8 -- regenerate c_lambda_V11_sweep per Luke's 7/7 comment: strip ALL
continuous-time framing.  Two candidate variants (existing paper figure is NOT
overwritten; Luke picks):

  (1) figs/c_lambda_V11_sweep_norefline.png  -- no benchmark line at all.
  (2) figs/c_lambda_V11_sweep_limitline.png  -- dashed benchmark relabeled as
      the exact rare-breakthrough limit
          V11 -> (1+2 b11) W / (3 (1+b11))     (lambda -> 0),
      which is ALSO the exact c -> infinity saturation value, so the same
      constant is the honest benchmark in both panels.  No continuous-time
      machinery anywhere.

Config matches task5_paper_figures.py Figure 3 exactly (W=2.5, lam=0.30,
beta11=0.10, c-sweep [0.76, 6.0], lam-sweep [0.02, 0.40]).
Also numerically confirms the limit constant at both ends (printed below).
"""
import math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE   = os.path.dirname(os.path.abspath(__file__))
FIGS   = os.path.join(os.path.dirname(HERE), "figs")
plt.rcParams.update({"axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 8.5})

def x11_and_V11(W, c, lam, b1):
    u = lam*lam*W*(1-b1*b1); Delta = 9*c*c - 2*c*u + u*u
    x11 = (3*c + u - math.sqrt(Delta))/(2*c*lam*(1+b1))
    return x11, c*x11/(2*lam*(1+b1)) + b1*W/(1+b1)

W0, LAM0, B11_F, C0_ = 2.5, 0.30, 0.10, 1.0
BLUE, PURPLE = "#4c78c8", "#b279a2"
VLIM = (1 + 2*B11_F)*W0/(3*(1 + B11_F))     # rare-breakthrough limit

cs   = np.linspace(0.76, 6.0, 300)
V11_c = np.array([x11_and_V11(W0, cc, LAM0, B11_F)[1] for cc in cs])
x11_c = np.array([x11_and_V11(W0, cc, LAM0, B11_F)[0] for cc in cs])
lams  = np.linspace(0.02, 0.40, 300)
V11_l = np.array([x11_and_V11(W0, C0_, ll, B11_F)[1] for ll in lams])
x11_l = np.array([x11_and_V11(W0, C0_, ll, B11_F)[0] for ll in lams])

# confirm the limit constant numerically at both ends
_, v_smalllam = x11_and_V11(W0, C0_, 1e-7, B11_F)
_, v_bigc     = x11_and_V11(W0, 1e7, LAM0, B11_F)
print(f"limit constant (1+2b)W/(3(1+b)) = {VLIM:.6f}")
print(f"  V11 at lam=1e-7:  {v_smalllam:.6f}   V11 at c=1e7: {v_bigc:.6f}")

def _zoom(ax, y, extra=(), frac=0.18):
    vals = np.concatenate([y] + [np.atleast_1d(e) for e in extra])
    lo, hi = vals.min(), vals.max(); pad = (hi-lo)*frac
    ax.set_ylim(lo-pad, hi+pad)

def make(fname, with_limitline):
    fig, (c1, c2) = plt.subplots(1, 2, figsize=(12.4, 5.0))
    # Panel A: V11 vs c
    c1.plot(cs, V11_c, color=PURPLE, lw=2.8, label=r"$V_{11}$", zorder=3)
    if with_limitline:
        c1.axhline(VLIM, color="grey", ls="--", lw=1.5,
                   label=r"rare-breakthrough limit $\frac{(1+2\beta_{11})W}{3(1+\beta_{11})}$")
        _zoom(c1, V11_c, (VLIM,))
    else:
        _zoom(c1, V11_c)
    c1.set_xlabel(r"cost curvature $c$")
    c1.set_ylabel(r"$V_{11}$", color=PURPLE); c1.tick_params(axis="y", labelcolor=PURPLE)
    c1.set_title(r"(A) $V_{11}$ increasing in $c$")
    c1.annotate(fr"+{100*(V11_c[-1]-V11_c[0])/V11_c[0]:.1f}% over range", xy=(0.52, 0.45),
                xycoords="axes fraction", ha="center", fontsize=9, color=PURPLE, style="italic")
    c1x = c1.twinx(); c1x.plot(cs, x11_c, color=BLUE, lw=1.6, ls=":", label=r"$x_{11}^*$ (right axis)")
    c1x.set_ylabel(r"$x_{11}^*$", color=BLUE); c1x.tick_params(axis="y", labelcolor=BLUE)
    h, l = c1.get_legend_handles_labels(); hx, lx = c1x.get_legend_handles_labels()
    c1.legend(h+hx, l+lx, loc="lower right", framealpha=0.9)
    # Panel B: V11 vs lambda
    c2.plot(lams, V11_l, color=PURPLE, lw=2.8, label=r"$V_{11}$", zorder=3)
    if with_limitline:
        c2.axhline(VLIM, color="grey", ls="--", lw=1.5,
                   label=r"rare-breakthrough limit ($\lambda\to0$)")
        _zoom(c2, V11_l, (VLIM,))
    else:
        _zoom(c2, V11_l)
    c2.set_xlabel(r"development rate $\lambda$")
    c2.set_ylabel(r"$V_{11}$", color=PURPLE); c2.tick_params(axis="y", labelcolor=PURPLE)
    c2.set_title(r"(B) $V_{11}$ decreasing in $\lambda$")
    c2.annotate(fr"{100*(V11_l[-1]-V11_l[0])/V11_l[0]:.1f}% over range", xy=(0.52, 0.55),
                xycoords="axes fraction", ha="center", fontsize=9, color=PURPLE, style="italic")
    c2x = c2.twinx(); c2x.plot(lams, x11_l, color=BLUE, lw=1.6, ls=":", label=r"$x_{11}^*$ (right axis)")
    c2x.set_ylabel(r"$x_{11}^*$", color=BLUE); c2x.tick_params(axis="y", labelcolor=BLUE)
    h, l = c2.get_legend_handles_labels(); hx, lx = c2x.get_legend_handles_labels()
    c2.legend(h+hx, l+lx, loc="upper right", framealpha=0.9)
    fig.suptitle(fr"$V_{{11}}$ depends on $c$ and $\lambda$   ($W={W0:g},\ \beta_{{11}}={B11_F:g}$)",
                 fontsize=10.5)
    fig.tight_layout()
    out = os.path.join(FIGS, fname)
    fig.savefig(out, dpi=300, bbox_inches="tight"); plt.close(fig)
    print("figure written:", out)

make("c_lambda_V11_sweep_norefline.png", with_limitline=False)
make("c_lambda_V11_sweep_limitline.png", with_limitline=True)
