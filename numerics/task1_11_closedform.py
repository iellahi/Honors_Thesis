"""
Task 1 -- (1,1) symmetric development state.
  (a) Closed form (subtracted root) vs TWO independent numerical solves of the
      Bellman + FOC system (no quadratic used):
        * 1D bisection on the FOC residual with V = Bellman value  [primary]
        * best-response value iteration                            [secondary]
  (b) Corner-regime map outside (A2): where interior x* > 1 forces x* = 1.

No scipy (blocked); hand-rolled root-finders are used as the independent check.
Closed form ported from verify_scripts/verify_step5b.py.  Figures: PNG only, 300 dpi.
Outputs: numerics/task1_results.txt, numerics/task1_corner_reference.csv, figs/corner_regime.png
"""
import math, random, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE   = os.path.dirname(os.path.abspath(__file__))
THESIS = os.path.dirname(HERE)
FIGS   = os.path.join(THESIS, "figs")
NUM    = os.path.join(THESIS, "numerics")
os.makedirs(FIGS, exist_ok=True)

log_lines = []
def log(s=""):
    print(s); log_lines.append(str(s))

# ---------------- closed form (reference: verify_step5b.py) ----------------
def closed_form(W, c, lam, b1):
    u = lam*lam*W*(1.0-b1*b1)                 # u = lam^2 W (1-b^2) >= 0
    Delta = 9*c*c - 2*c*u + u*u               # = (u-c)^2 + 8c^2 > 0 always
    root = math.sqrt(Delta)
    denom = 2*c*lam*(1.0+b1)
    x_minus = (3*c + u - root)/denom          # subtracted root  -> interior
    x_plus  = (3*c + u + root)/denom          # additive root    -> P>1
    V = c*x_minus/(2*lam*(1.0+b1)) + b1*W/(1.0+b1)
    return x_minus, x_plus, V, Delta

def quad_resid(x, W, c, lam, b1):
    # A x^2 + B x + C  with A=c lam(1+b), B=-(3c+lam^2 W(1-b^2)), C=2 lam W(1-b)
    A = c*lam*(1+b1)
    B = -(3*c + lam*lam*W*(1-b1*b1))
    C = 2*lam*W*(1-b1)
    return A*x*x + B*x + C

# ---------------- primitives: Bellman + FOC (independent, NO quadratic) -----
def Vbell(x, W, c, lam, b1):
    M = lam*(1.0+b1)
    return (W*(2*M - M*M*x) - c*x)/(2*M*(2 - M*x))     # symmetric Bellman value

def foc_resid(x, W, c, lam, b1):
    # symmetric FOC (verify_step2 CHECK 4) with V replaced by the Bellman value:
    #   c x - [ lam W - M V - M^2 x (W/2 - V) ]
    M = lam*(1.0+b1)
    V = Vbell(x, W, c, lam, b1)
    return c*x - (lam*W - M*V - M*M*x*(W/2 - V))

def solve_bisect(W, c, lam, b1, tol=1e-15, itmax=300):
    """Independent interior solve: root of FOC residual on (0, 1/M)."""
    M = lam*(1.0+b1)
    a, b = 1e-13, (1.0-1e-13)/M
    fa, fb = foc_resid(a,W,c,lam,b1), foc_resid(b,W,c,lam,b1)
    if fa*fb > 0:
        return float('nan')
    for _ in range(itmax):
        mid = 0.5*(a+b); fm = foc_resid(mid,W,c,lam,b1)
        if fa*fm <= 0: b, fb = mid, fm
        else:          a, fa = mid, fm
        if b-a < tol: break
    return 0.5*(a+b)

def solve_br_iter(W, c, lam, b1, tol=1e-14, itmax=20000):
    """Best-response value iteration: fully independent of the closed form."""
    M = lam*(1.0+b1)
    V = W/4.0; x = min(0.3, 0.9/M)
    for _ in range(itmax):
        def foc_xi(xi):                       # FOC in xi given opponent x, continuation V
            Pi = lam*(xi + b1*x); Pj = lam*(x + b1*xi)
            return -c*xi + lam*(W - (Pj+b1*Pi)*W/2 - (1+b1)*V + (Pj+b1*Pi)*V)
        a, bb = 1e-13, (1.0-1e-13)/M
        fa, fb = foc_xi(a), foc_xi(bb)
        if fa*fb > 0:
            xi = a if abs(fa) < abs(fb) else bb
        else:
            for _ in range(300):
                mid = 0.5*(a+bb); fm = foc_xi(mid)
                if fa*fm <= 0: bb, fb = mid, fm
                else:          a, fa = mid, fm
                if bb-a < 1e-15: break
            xi = 0.5*(a+bb)
        Vn = Vbell(xi, W, c, lam, b1)
        if abs(xi-x) < tol and abs(Vn-V) < tol:
            x, V = xi, Vn; break
        x, V = xi, Vn
    return x, V

# ==========================================================================
# PART A -- validation sweep in the (A2) region
# ==========================================================================
log("="*74)
log("PART A: closed form vs independent Bellman+FOC solves  (A2 region: lamW<=c)")
log("="*74)
random.seed(11)
N = 60000
max_dx_bis = max_dV_bis = 0.0
max_dx_br  = max_dV_br  = 0.0
max_qres   = 0.0
n_br = 0
viol_interior = viol_prob = 0
for i in range(N):
    lam = random.uniform(1e-2, 0.4999)
    b1  = random.uniform(0.0, 0.9995)
    c   = random.uniform(5e-2, 10.0)
    W   = random.uniform(1e-3, c/lam)          # (A2): lam W <= c
    xm, xp, V, Delta = closed_form(W, c, lam, b1)
    # primary independent solve
    xb = solve_bisect(W, c, lam, b1)
    Vb = Vbell(xb, W, c, lam, b1)
    max_dx_bis = max(max_dx_bis, abs(xb-xm))
    max_dV_bis = max(max_dV_bis, abs(Vb-V))
    max_qres   = max(max_qres, abs(quad_resid(xm, W, c, lam, b1)))
    # interiority / probability under (A2)
    M = lam*(1+b1)
    if not (0 < xm <= 1.0+1e-12): viol_interior += 1
    if not (M*xm < 1.0):          viol_prob += 1
    # secondary independent solve on a subsample
    if i % 20 == 0:
        xi, Vi = solve_br_iter(W, c, lam, b1)
        max_dx_br = max(max_dx_br, abs(xi-xm))
        max_dV_br = max(max_dV_br, abs(Vi-V))
        n_br += 1

log(f"draws N = {N}   (seed 11, parity with verify_step5b)")
log(f"max |x_closed - x_bisect|          = {max_dx_bis:.3e}")
log(f"max |V_closed - V_bisect|          = {max_dV_bis:.3e}")
log(f"max |x_closed - x_BRiter| (n={n_br}) = {max_dx_br:.3e}")
log(f"max |V_closed - V_BRiter| (n={n_br}) = {max_dV_br:.3e}")
log(f"max |quadratic residual at x*|     = {max_qres:.3e}")
log(f"interiority x* in (0,1] violations = {viol_interior}   (expect 0 under A2)")
log(f"probability M x* < 1 violations    = {viol_prob}   (expect 0 under A2)")

# ==========================================================================
# PART B -- corner regime outside (A2)
# ==========================================================================
log("")
log("="*74)
log("PART B: corner regime outside (A2)  (interior x- > 1  =>  corner x* = 1)")
log("="*74)
# reference point from the task brief
lam0, c0, W0, b0 = 0.49, 0.1, 50.0, 0.0
xm0, xp0, V0, D0 = closed_form(W0, c0, lam0, b0)
log(f"reference lam={lam0}, c={c0}, W={W0}, b11={b0}:  lamW={lam0*W0:.3f} (>> c={c0}), "
    f"x- = {xm0:.4f}  ->  corner x* = {min(xm0,1.0):.1f}")
log(f"  (brief expected x- ~ 4.05; match = {abs(xm0-4.05) < 0.02})")

# does interior x- exceed 1 exactly at the (A2) line lamW=c?  -> NO: A2 is conservative.
log("")
log("Interior x- at the (A2) boundary lamW=c (b11=0) -- A2 is sufficient, not tight:")
for lam in (0.10, 0.25, 0.49):
    Wb = c0/lam
    xmb, *_ = closed_form(Wb, c0, lam, 0.0)
    log(f"  lam={lam:.2f}: at lamW=c (W={Wb:.3f})  x- = {xmb:.4f}  (< 1  => interior still holds)")

# per-lambda true corner onset: smallest W with x-(W) = 1  (bisection on W)
def x_of_W(W, lam, c=c0, b1=0.0):
    return closed_form(W, c, lam, b1)[0]
def corner_onset_W(lam, c=c0, b1=0.0):
    lo, hi = 1e-6, 1.0
    if x_of_W(lo,lam) > 1: return lo
    while x_of_W(hi,lam) < 1 and hi < 1e6: hi *= 2
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if x_of_W(mid,lam) < 1: lo = mid
        else: hi = mid
    return 0.5*(lo+hi)
log("")
log("True corner onset W (x-=1) vs (A2) boundary W=c/lam  (b11=0, c=0.1):")
onset_rows = []
for lam in (0.10, 0.25, 0.40, 0.49):
    Won = corner_onset_W(lam); Wa2 = c0/lam
    onset_rows.append((lam, Wa2, Won, lam*Won/c0))
    log(f"  lam={lam:.2f}:  A2 boundary W={Wa2:.3f}   corner onset W={Won:.3f}   "
        f"(lamW/c at onset = {lam*Won/c0:.3f} > 1)")

# save reference CSV
with open(os.path.join(NUM, "task1_corner_reference.csv"), "w") as f:
    f.write("lambda,c,W_A2_boundary,W_corner_onset,ratio_lamW_over_c_at_onset,x_at_A2_boundary\n")
    for lam, Wa2, Won, ratio in onset_rows:
        xa2 = closed_form(Wa2, c0, lam, 0.0)[0]
        f.write(f"{lam},{c0},{Wa2:.6f},{Won:.6f},{ratio:.6f},{xa2:.6f}\n")

# ==========================================================================
# FIGURE -- corner_regime.png  (2 panels)
# ==========================================================================
fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.4, 5.0))

# Panel A: x-(W) and clipped x* vs W for two lambdas, c=0.1, b11=0
c_fig = 0.1
for lam, col in [(0.25, "#1f77b4"), (0.49, "#d62728")]:
    Ws = np.linspace(1e-3, 2.0, 800)
    xs = np.array([closed_form(W, c_fig, lam, 0.0)[0] for W in Ws])
    xclip = np.minimum(xs, 1.0)
    axA.plot(Ws, xs, "--", color=col, lw=1.6, alpha=0.85,
             label=fr"$x^-$ (interior root), $\lambda={lam}$")
    axA.plot(Ws, xclip, "-", color=col, lw=2.4,
             label=fr"$x^*=\min(x^-,1)$, $\lambda={lam}$")
    Wa2 = c_fig/lam
    axA.axvline(Wa2, color=col, ls=":", lw=1.2, alpha=0.7)
    axA.text(Wa2-0.02, 0.45, fr"(A2) $\lambda W=c$  ($\lambda={lam}$)", color=col,
             rotation=90, ha="right", va="center", fontsize=7.5)
axA.axhline(1.0, color="grey", lw=1.0, alpha=0.6)
axA.set_xlim(0, 2.0); axA.set_ylim(0, 2.0)
axA.set_xlabel(r"prize $W$"); axA.set_ylabel(r"symmetric $(1,1)$ effort")
axA.set_title(r"(A) interior root $x^-$ vs feasible $x^*$  ($c=0.1,\ \beta_{11}=0$)")
axA.legend(fontsize=7.5, loc="upper right", framealpha=0.9)
axA.text(0.02, 0.03,
         r"ref: $\lambda{=}0.49,c{=}0.1,W{=}50\Rightarrow x^-\!\approx4.05$ (deep corner, off-scale)",
         transform=axA.transAxes, fontsize=7.3, style="italic", color="#333333")

# Panel B: regime map in (lambda, W), c=0.1, b11=0
lam_g = np.linspace(0.05, 0.49, 320)
W_g   = np.linspace(1e-3, 4.0, 320)
LAM, WW = np.meshgrid(lam_g, W_g)
XM = np.vectorize(lambda l, w: closed_form(w, c_fig, l, 0.0)[0])(LAM, WW)
# shaded interior root value (clipped for color), corner region hatched
pcm = axB.pcolormesh(LAM, WW, np.minimum(XM, 1.0), cmap="viridis",
                     shading="auto", vmin=0, vmax=1.0)
cb = fig.colorbar(pcm, ax=axB); cb.set_label(r"feasible effort $x^*=\min(x^-,1)$")
# corner region (x- >= 1)
axB.contourf(LAM, WW, (XM >= 1.0).astype(float), levels=[0.5, 1.5],
             colors="none", hatches=["////"])
axB.contour(LAM, WW, XM, levels=[1.0], colors="white", linewidths=2.2)
axB.plot(lam_g, c_fig/lam_g, color="red", lw=2.0)
axB.set_xlim(0.05, 0.49); axB.set_ylim(0, 4.0)
axB.set_xlabel(r"development rate $\lambda$"); axB.set_ylabel(r"prize $W$")
axB.set_title(r"(B) regime map: interior vs corner  ($c=0.1,\ \beta_{11}=0$)")
_proxies = [Line2D([0],[0], color="red",   lw=2.0, label=r"(A2) boundary $\lambda W=c$"),
            Line2D([0],[0], color="white", lw=2.2, label=r"corner onset $x^-=1$")]
axB.legend(handles=_proxies, fontsize=8, loc="upper right", framealpha=0.92)
axB.text(0.30, 3.3, "corner:  $x^*=1$\n(A2 fails)", ha="center", va="center",
         fontsize=10, color="black",
         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.82))
axB.text(0.40, 0.16, "interior:  $x^*\\in(0,1)$", ha="center", va="center",
         fontsize=9, color="black",
         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.82))

fig.tight_layout()
outpng = os.path.join(FIGS, "corner_regime.png")
fig.savefig(outpng, dpi=300, bbox_inches="tight")
log("")
log(f"figure written: {outpng}")

# note the gap between the two boundaries
log("")
log("NOTE: the white x-=1 contour lies ABOVE the red lamW=c line -> (A2) is a")
log("      conservative sufficient condition; interior equilibria persist in the")
log("      band between the two boundaries before the corner binds.")

with open(os.path.join(NUM, "task1_results.txt"), "w") as f:
    f.write("\n".join(log_lines) + "\n")
print("\n[written]", os.path.join(NUM, "task1_results.txt"))
