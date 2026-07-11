"""
Task 9 -- clean best-response figure for Section 4.2 (Luke: the asymmetric
equilibrium "is easy to demonstrate with best-response figures").
Decluttered variant of laggard_BR panel A: smooth curves, equilibrium point,
45-degree line showing leader dominance (z* > x*). Same representative config
as task4 (lam = mu = 0.30, c = 1, beta11 = 0.30, W = 0.85 c / lam).
Output: figs/best_response_clean.png (300 dpi).
"""
import math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figs")

def V11_of(W, c, lam, b1):
    u = lam*lam*W*(1-b1*b1); D = 9*c*c - 2*c*u + u*u
    x11 = (3*c + u - math.sqrt(D))/(2*c*lam*(1+b1))
    return c*x11/(2*lam*(1+b1)) + b1*W/(1+b1)
def br_x(z, W, c, lam, mu, V11):
    rad = c*c*lam*lam*z*z + 2*c*mu*mu*lam*z*(1-lam*z)**2*V11
    return (-c*lam*z + math.sqrt(rad))/(c*mu*(1-lam*z))
def br_z(x, W, c, lam, mu, V11):
    rad = c*c*mu*mu*x*x + 2*c*mu*lam*lam*x*(1-mu*x)*(W-V11)
    return (-c*mu*x + math.sqrt(rad))/(c*lam*(1-mu*x))

lam, mu, c, b1 = 0.30, 0.30, 1.0, 0.30
W = 0.85*c/max(lam, mu)
V11 = V11_of(W, c, lam, b1)
z, x = 0.5, 0.5
for _ in range(5000):
    zn = br_z(x, W, c, lam, mu, V11); xn = br_x(zn, W, c, lam, mu, V11)
    if abs(zn-z) < 1e-13 and abs(xn-x) < 1e-13: z, x = zn, xn; break
    z, x = zn, xn

zs  = np.linspace(1e-3, 1.0, 600)
Xz  = np.array([br_x(zz, W, c, lam, mu, V11) for zz in zs])
xs2 = np.linspace(1e-3, 1.0, 600)
Zx  = np.array([br_z(xx, W, c, lam, mu, V11) for xx in xs2])

fig, ax = plt.subplots(figsize=(6.4, 5.4))
ax.plot(zs, Xz, color="#54a24b", lw=2.4, label=r"laggard best response $x = X(z)$")
ax.plot(Zx, xs2, color="#4c78c8", lw=2.4, label=r"leader best response $z = Z(x)$")
ax.plot([0, 1], [0, 1], color="grey", ls="--", lw=1.2, label=r"$45^\circ$ line ($z = x$)")
ax.plot([z], [x], "k*", ms=16, zorder=5)
ax.annotate(fr"$(z^*, x^*) = ({z:.2f}, {x:.2f})$", xy=(z, x), xytext=(z+0.08, x-0.07),
            fontsize=10, arrowprops=dict(arrowstyle="-", color="grey", lw=0.8))
ax.set_xlim(0, 0.8); ax.set_ylim(0, 0.8)
ax.set_xlabel(r"leader effort $z$"); ax.set_ylabel(r"laggard effort $x$")
ax.set_title(fr"Best responses in the asymmetric state"
             + f"\n($\\lambda=\\mu={lam:g},\\ c={c:g},\\ \\beta_{{11}}={b1:g},\\ W={W:.2f}$)",
             fontsize=10.5)
ax.legend(fontsize=9, loc="upper left", framealpha=0.95)
fig.tight_layout()
out = os.path.join(FIGS, "best_response_clean.png")
fig.savefig(out, dpi=300, bbox_inches="tight")
print("written:", out, f"  eq: z*={z:.4f}, x*={x:.4f} (z*>x*: {z>x})")
