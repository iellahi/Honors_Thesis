"""
Task 6 -- consolidated confirmation scan.
One master scan (seed 11, N=60000, full assumption region, half sym / half asym)
that re-checks EVERY pointwise claim in one place, plus a finite-difference subsample
for the transmission claims. Emits a machine-readable confirmation table.

Pointwise claims:  (1,1) interiority & prob; value ordering; V11<W/2; Prop 2 (mu=lam);
z>x share (mu!=lam); Delta0>0; x00* in (0,1); m x00*<1; m x00*<1-b0; det J>0; provocation share.
Finite-difference (subsample): dx*/dbeta11>0; dz*/dbeta11 sign.

Solvers ported from verify_step5b.py / verify_step5d.py.
Outputs: numerics/task6_confirmation_table.csv, numerics/task6_results.txt
"""
import math, random, os
HERE = os.path.dirname(os.path.abspath(__file__)); NUM = HERE
log_lines = []
def log(s=""):
    print(s); log_lines.append(str(s))

def x11_V11(W, c, lam, b1):
    u = lam*lam*W*(1-b1*b1); D = 9*c*c - 2*c*u + u*u
    x = (3*c + u - math.sqrt(D))/(2*c*lam*(1+b1))
    return x, c*x/(2*lam*(1+b1)) + b1*W/(1+b1)
def br_x(z, W, c, lam, mu, V):
    r = c*c*lam*lam*z*z + 2*c*mu*mu*lam*z*(1-lam*z)**2*V
    return (-c*lam*z + math.sqrt(r))/(c*mu*(1-lam*z))
def br_z(x, W, c, lam, mu, V):
    r = c*c*mu*mu*x*x + 2*c*mu*lam*lam*x*(1-mu*x)*(W-V)
    return (-c*mu*x + math.sqrt(r))/(c*lam*(1-mu*x))
def solve_fp(W, c, lam, mu, V, tol=1e-13, itmax=5000):
    z = x = 0.5
    for _ in range(itmax):
        zn = br_z(x, W, c, lam, mu, V); xn = br_x(zn, W, c, lam, mu, V)
        if abs(zn-z) < tol and abs(xn-x) < tol: z, x = zn, xn; break
        z, x = zn, xn
    GL = (c*lam/2)*(1-mu*x)*z*z + c*mu*x*z - lam*mu*x*(W-V)
    GF = (c*mu/2)*(1-lam*z)*x*x + c*lam*z*x - mu*lam*z*(1-lam*z)*V
    return z, x, max(abs(GL), abs(GF))
def vals(z, x, W, c, lam, mu, V):
    PL, PF = lam*z, mu*x; S = PL+PF-PL*PF
    return (PL*W+(1-PL)*PF*V-c*z*z/2)/S, ((1-PL)*PF*V-c*x*x/2)/S
def x00_D0(c, mu, b0, V10, V01, V11):
    m = mu*(1+b0); G = V10+b0*V01-(1+b0)*V11; D = V10-V01
    A, B, C = c*m, -(3*c+2*mu*m*G), 2*mu*(1-b0)*D
    D0 = B*B-4*A*C
    x = (-B-math.sqrt(D0))/(2*A) if D0 > 0 else float('nan')
    return x, D0, m
def detJ(z, x, W, c, lam, mu, V, h=1e-7):
    GL = lambda z_, x_: (c*lam/2)*(1-mu*x_)*z_*z_ + c*mu*x_*z_ - lam*mu*x_*(W-V)
    GF = lambda z_, x_: (c*mu/2)*(1-lam*z_)*x_*x_ + c*lam*z_*x_ - mu*lam*z_*(1-lam*z_)*V
    a11 = (GL(z+h,x)-GL(z-h,x))/(2*h); a12 = (GL(z,x+h)-GL(z,x-h))/(2*h)
    a21 = (GF(z+h,x)-GF(z-h,x))/(2*h); a22 = (GF(z,x+h)-GF(z,x-h))/(2*h)
    return a11*a22 - a12*a21

random.seed(11); N = 60000
c_int=c_prob=0; ord_ok=0; vhalf=0
psym=fsym=0; pasym=fasym=0
D0_bad=x00_bad=mx00_bad=mxb0_bad=0; minD0=1e18
detbad=0; prov=disc=0
dxb_tot=dxb_pos=0; dzb_pos=dzb_neg=0
n_asym_conv=0; skip=0
worst = {}   # collect first violating tuple per claim
for i in range(N):
    lam = random.uniform(1e-2, 0.4999)
    mu  = random.uniform(1e-2, 0.4999) if i % 2 else lam
    b1  = random.uniform(1e-4, 0.9995); b0 = random.uniform(0.0, 0.9995)
    c   = random.uniform(5e-2, 10.0);   W  = random.uniform(1e-3, c/max(lam, mu))
    x11, V11 = x11_V11(W, c, lam, b1); M = lam*(1+b1)
    if 0 < x11 <= 1+1e-12: c_int += 1
    elif "int11" not in worst: worst["int11"] = (lam,mu,b1,b0,c,W)
    if M*x11 < 1: c_prob += 1
    if V11 < W/2 + 1e-12: vhalf += 1
    z, x, res = solve_fp(W, c, lam, mu, V11)
    if res > 1e-9: skip += 1; continue
    n_asym_conv += 1
    V10, V01 = vals(z, x, W, c, lam, mu, V11)
    if 0 <= V01 <= V11 <= V10 <= W+1e-12: ord_ok += 1
    elif "order" not in worst: worst["order"] = (lam,mu,b1,b0,c,W)
    if mu == lam:
        if z > x: psym += 1
        else: fsym += 1; worst.setdefault("prop2", (lam,mu,b1,b0,c,W))
    else:
        if z > x: pasym += 1
        else: fasym += 1
    Nmarg = mu*x*(1-lam*z)**2 - 2*lam*lam*z*z
    if Nmarg > 0: prov += 1
    else: disc += 1
    if detJ(z, x, W, c, lam, mu, V11) <= 0: detbad += 1; worst.setdefault("detJ",(lam,mu,b1,b0,c,W))
    x00, D0, m = x00_D0(c, mu, b0, V10, V01, V11)
    minD0 = min(minD0, D0)
    if D0 <= 0: D0_bad += 1; worst.setdefault("D0",(lam,mu,b1,b0,c,W)); continue
    if not (0 < x00 < 1): x00_bad += 1; worst.setdefault("x00",(lam,mu,b1,b0,c,W))
    if not (m*x00 < 1): mx00_bad += 1; worst.setdefault("mx00",(lam,mu,b1,b0,c,W))
    if not (m*x00 < 1-b0): mxb0_bad += 1; worst.setdefault("mxb0",(lam,mu,b1,b0,c,W))
    # transmission finite diff on a subsample
    if i % 15 == 0:
        h = 1e-6
        if b1+h < 1:
            _, V11h = x11_V11(W, c, lam, b1+h); zh, xh, rh = solve_fp(W, c, lam, mu, V11h)
            if rh < 1e-9:
                dxb_tot += 1
                if (xh-x)/h > 0: dxb_pos += 1
                if (zh-z)/h > 0: dzb_pos += 1
                else: dzb_neg += 1

def pct(a, b): return 100.0*a/b if b else float('nan')
log("="*74); log("TASK 6 -- consolidated confirmation scan (seed 11, N=60000)"); log("="*74)
log(f"converged asymmetric draws: {n_asym_conv}  (skipped: {skip})")
log("")
# build table rows: (claim, kind, checked, violations, key_stat)
rows = [
 ("(1,1) interiority 0<x11*<=1", "T-cross", N, N-c_int, f"{pct(c_int,N):.3f}% ok"),
 ("(1,1) probability M x11*<1", "T-cross", N, N-c_prob, f"{pct(c_prob,N):.3f}% ok"),
 ("V11 < W/2", "T-cross", N, N-vhalf, f"{pct(vhalf,N):.3f}% ok"),
 ("value ordering 0<=V01<=V11<=V10<=W", "T-cross", n_asym_conv, n_asym_conv-ord_ok, f"{pct(ord_ok,n_asym_conv):.3f}% ok"),
 ("Prop 2: z*>x* at mu=lam", "T-cross", psym+fsym, fsym, f"{fsym} fails"),
 ("z*>x* at mu!=lam (binding hyp.)", "N", pasym+fasym, fasym, f"{pct(fasym,pasym+fasym):.2f}% fail (~22%)"),
 ("Delta0 > 0", "N", n_asym_conv, D0_bad, f"min={minD0:.4e}"),
 ("x00* in (0,1)", "N", n_asym_conv, x00_bad, f"{x00_bad} viol"),
 ("m x00* < 1", "N", n_asym_conv, mx00_bad, f"{mx00_bad} viol"),
 ("m x00* < 1-b0  (=> dx00*/dV10>0)", "T-cross", n_asym_conv, mxb0_bad, f"{mxb0_bad} viol (unconditional)"),
 ("det J > 0", "N", n_asym_conv, detbad, f"{detbad} viol"),
 ("provocation region N>0", "N", prov+disc, disc, f"{pct(prov,prov+disc):.2f}% (~99.8%)"),
 ("dx*/dbeta11 > 0", "N", dxb_tot, dxb_tot-dxb_pos, f"{pct(dxb_pos,dxb_tot):.2f}% positive"),
 ("dz*/dbeta11 < 0 (ambiguous)", "N", dzb_pos+dzb_neg, dzb_pos, f"{pct(dzb_neg,dzb_pos+dzb_neg):.2f}% negative"),
]
log(f"{'claim':45} {'kind':8} {'checked':>8} {'viol':>6}  key")
for cl, kind, ch, vi, ks in rows:
    log(f"{cl:45} {kind:8} {ch:>8} {vi:>6}  {ks}")
log("")
tot_viol = sum(r[3] for r in rows if r[1] in ("T-cross",) ) + D0_bad+x00_bad+mx00_bad+detbad
log(f"HARD violations of proven/theorem cross-checks: "
    f"{(N-c_int)+(N-c_prob)+(N-vhalf)+(n_asym_conv-ord_ok)+fsym+mxb0_bad} (expect 0)")
log(f"HARD violations of (0,0) numerical claims (Delta0/x00/mx00) and det J: "
    f"{D0_bad+x00_bad+mx00_bad+detbad} (expect 0)")

with open(os.path.join(NUM, "task6_confirmation_table.csv"), "w") as f:
    f.write("claim,kind,checked,violations,key_stat\n")
    for cl, kind, ch, vi, ks in rows:
        f.write(f"\"{cl}\",{kind},{ch},{vi},\"{ks}\"\n")
if worst:
    log(""); log("first violating tuple per claim (empty => none):")
    for k, v in worst.items(): log(f"  {k}: {v}")
else:
    log(""); log("NO violations recorded for any claim across the scan.")

with open(os.path.join(NUM, "task6_results.txt"), "w") as f:
    f.write("\n".join(log_lines) + "\n")
print("\n[written]", os.path.join(NUM, "task6_confirmation_table.csv"))
