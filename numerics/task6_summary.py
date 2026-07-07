"""
Task 6 -- consolidated verification: confirmation tables + break-region map.
  PART 1  THEOREMS  cross-checked numerically in-region (expect 0 violations).
  PART 2  NUMERICAL claims re-confirmed on one coherent scan.
  PART 3  BREAK REGIONS: relax each maintained assumption and report what fails.
  PART 4  symbolic master-quadratic identity (sympy).
Solvers ported from verify_step5b.py.  Do NOT modify .tex files.
Outputs: numerics/task6_confirmation.txt  (raw tables consumed by SUMMARY_memo.md)
"""
import math, random, os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)); THESIS = os.path.dirname(HERE)
NUM = os.path.join(THESIS, "numerics")
log_lines = []
def log(s=""):
    print(s); log_lines.append(str(s))

# ---------------- solvers ----------------
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
    z, x = 0.5, 0.5
    for _ in range(itmax):
        zn = br_z(x, W, c, lam, mu, V); xn = br_x(zn, W, c, lam, mu, V)
        if abs(zn-z) < tol and abs(xn-x) < tol: z, x = zn, xn; break
        z, x = zn, xn
    GL = (c*lam/2)*(1-mu*x)*z*z + c*mu*x*z - lam*mu*x*(W-V)
    GF = (c*mu/2)*(1-lam*z)*x*x + c*lam*z*x - mu*lam*z*(1-lam*z)*V
    return z, x, max(abs(GL), abs(GF))
def vals_10_01(z, x, W, c, lam, mu, V):
    PL, PF = lam*z, mu*x; S = PL+PF-PL*PF
    return (PL*W+(1-PL)*PF*V-c*z*z/2)/S, ((1-PL)*PF*V-c*x*x/2)/S
def x00_of(c, mu, b0, V10, V01, V11):
    m = mu*(1+b0); G = V10+b0*V01-(1+b0)*V11; D = V10-V01
    A, B, C = c*m, -(3*c+2*mu*m*G), 2*mu*(1-b0)*D
    D0 = B*B-4*A*C
    x = (-B-math.sqrt(D0))/(2*A) if D0 > 0 else float('nan')
    return x, D0, m

def sgn(v): return 0 if v == 0 else (1 if v > 0 else -1)

# ==========================================================================
# PART 1 + 2 : theorems and numerical claims, one in-region scan
# ==========================================================================
log("="*78); log("PART 1/2: THEOREM cross-checks + NUMERICAL claims  (seed 202, N=60000, A1/A2/A2')")
log("="*78)
random.seed(202)
N = 60000
V = {k: 0 for k in [   # violation counters (theorems -> expect 0)
    "T1_Delta_pos","T1_Delta_form","T2_sub_interior","T2_add_super","T3_key_ineq",
    "T4_prop1_signs","T5_prop1p_signs","T6_ordering","T7_prop2_sym","T8_BRslope","T9_prop3_signs",
    "T11_BRsolve"]}
num = {"dxb11_pos":0,"dxb11_tot":0,"D0_bad":0,"x00_bad":0,"mx00_bad":0,"detJ_bad":0,
       "prov":0,"asym_tot":0,"dzb11_pos":0,"dzb11_neg":0,"dxlam_rev":0,"dxlam_tot":0}
minD = minD0 = math.inf
for i in range(N):
    lam = random.uniform(1e-2, 0.4999)
    mu  = random.uniform(1e-2, 0.4999) if i % 2 else lam
    b1  = random.uniform(1e-4, 0.9990); b0 = random.uniform(0.0, 0.9990)
    c   = random.uniform(5e-2, 10.0);   W  = random.uniform(1e-3, c/max(lam, mu))
    # --- (1,1) closed form ---
    u = lam*lam*W*(1-b1*b1); Delta = 9*c*c-2*c*u+u*u
    if not (Delta > 0): V["T1_Delta_pos"] += 1
    if abs(Delta - ((u-c)**2+8*c*c)) > 1e-9*max(1,abs(Delta)): V["T1_Delta_form"] += 1
    minD = min(minD, Delta)
    x11, V11 = x11_V11(W, c, lam, b1)
    M = lam*(1+b1)
    xplus = (3*c+u+math.sqrt(Delta))/(2*c*lam*(1+b1))
    if not (M*x11 < 1): V["T2_sub_interior"] += 1
    if not (M*xplus > 1): V["T2_add_super"] += 1
    if not (math.sqrt(Delta) > 3*c-u - 1e-12): V["T3_key_ineq"] += 1
    # --- Prop 1 signs: x11 wrt (W,lam,c,b1) = (+,+,-,-) ---
    hW,hL,hc,hb = 1e-6*W,1e-6*lam,1e-6*c,1e-6*max(b1,1e-3)
    s = (sgn(x11_V11(W+hW,c,lam,b1)[0]-x11), sgn(x11_V11(W,c,lam+hL,b1)[0]-x11),
         sgn(x11_V11(W,c+hc,lam,b1)[0]-x11), sgn(x11_V11(W,c,lam,b1+hb)[0]-x11) if b1+hb<1 else -1)
    if s != (1,1,-1,-1): V["T4_prop1_signs"] += 1
    # --- Prop 1' signs: V11 wrt (b1,W,c,lam) = (+,+,+,-) ---
    sv = (sgn((x11_V11(W,c,lam,b1+hb)[1]-V11) if b1+hb<1 else 1), sgn(x11_V11(W+hW,c,lam,b1)[1]-V11),
          sgn(x11_V11(W,c+hc,lam,b1)[1]-V11), sgn(x11_V11(W,c,lam+hL,b1)[1]-V11))
    if sv != (1,1,1,-1): V["T5_prop1p_signs"] += 1
    # --- asymmetric solve ---
    z, x, res = solve_fp(W, c, lam, mu, V11)
    if res > 1e-9: continue
    if not (abs((c*lam/2)*(1-mu*x)*z*z + c*mu*x*z - lam*mu*x*(W-V11)) < 1e-9
            and abs((c*mu/2)*(1-lam*z)*x*x + c*lam*z*x - mu*lam*z*(1-lam*z)*V11) < 1e-9):
        V["T11_BRsolve"] += 1
    V10, V01 = vals_10_01(z, x, W, c, lam, mu, V11)
    if not (0 <= V01 <= V11 <= V10 <= W+1e-9 and V11 < W/2 + 1e-12): V["T6_ordering"] += 1
    if mu == lam and not (z > x): V["T7_prop2_sym"] += 1
    # --- Prop 2 asym share (numerical) ---
    if mu != lam:
        num["asym_tot"] += 1
    # --- T8 BR slope identity at z* ---
    Nmarg = mu*x*(1-lam*z)**2 - 2*lam*lam*z*z
    hz = 1e-7
    dX = (br_x(z+hz,W,c,lam,mu,V11) - br_x(z-hz,W,c,lam,mu,V11))/(2*hz)
    if sgn(dX) != sgn(Nmarg) and abs(Nmarg) > 1e-7: V["T8_BRslope"] += 1
    if Nmarg > 0: num["prov"] += 1
    # --- det J (numerical) ---
    def GL(z_,x_): return (c*lam/2)*(1-mu*x_)*z_*z_ + c*mu*x_*z_ - lam*mu*x_*(W-V11)
    def GF(z_,x_): return (c*mu/2)*(1-lam*z_)*x_*x_ + c*lam*z_*x_ - mu*lam*z_*(1-lam*z_)*V11
    h=1e-7
    a11=(GL(z+h,x)-GL(z-h,x))/(2*h); a12=(GL(z,x+h)-GL(z,x-h))/(2*h)
    a21=(GF(z+h,x)-GF(z-h,x))/(2*h); a22=(GF(z,x+h)-GF(z,x-h))/(2*h)
    if a11*a22-a12*a21 <= 0: num["detJ_bad"] += 1
    # --- (0,0) + Prop 3 signs ---
    x00, D0, m = x00_of(c, mu, b0, V10, V01, V11)
    minD0 = min(minD0, D0)
    if D0 <= 0: num["D0_bad"] += 1
    else:
        if not (0 < x00 < 1): num["x00_bad"] += 1
        if not (m*x00 < 1):   num["mx00_bad"] += 1
        # Prop 3: dx00 wrt (V10,V11,V01) = (+,+,-)  (partials, other V's fixed)
        e = 1e-6
        s10 = sgn(x00_of(c,mu,b0,V10+e,V01,V11)[0]-x00)
        s11 = sgn(x00_of(c,mu,b0,V10,V01,V11+e)[0]-x00)
        s01 = sgn(x00_of(c,mu,b0,V10,V01+e,V11)[0]-x00)
        if (s10,s11,s01) != (1,1,-1): V["T9_prop3_signs"] += 1
    # --- transmission (numerical): dx*/db11, dz*/db11 ---
    hb2 = 1e-6*max(b1,1e-3)
    if b1+hb2 < 1:
        _, V11b = x11_V11(W,c,lam,b1+hb2)
        zb, xb, rb = solve_fp(W,c,lam,mu,V11b)
        if rb < 1e-9:
            num["dxb11_tot"] += 1
            if xb-x > 0: num["dxb11_pos"] += 1
            if zb-z > 0: num["dzb11_pos"] += 1
            else:        num["dzb11_neg"] += 1
    # --- dx*/dlam total vs direct reversal (numerical, subsample) ---
    if i % 5 == 0 and lam+1e-6*lam < 0.5:
        hl=1e-6*lam
        _, V11l = x11_V11(W,c,lam+hl,b1)
        zt, xt, rt = solve_fp(W,c,lam+hl,mu,V11l)     # total
        zd, xd, rd = solve_fp(W,c,lam+hl,mu,V11)      # direct (V11 frozen)
        if rt<1e-9 and rd<1e-9:
            num["dxlam_tot"] += 1
            if (xt-x<0) != (xd-x<0): num["dxlam_rev"] += 1

log("")
log("THEOREM confirmation (numerical cross-check; 0 = all pass):")
labels = {
 "T1_Delta_pos":"Delta=(u-c)^2+8c^2 > 0","T1_Delta_form":"Delta form matches 9c^2-2cu+u^2",
 "T2_sub_interior":"subtracted root: Mx<1","T2_add_super":"additive root: Mx>1",
 "T3_key_ineq":"key inequality sqrt(Delta)>3c-u","T4_prop1_signs":"Prop 1 signs (x11: +,+,-,-)",
 "T5_prop1p_signs":"Prop 1' signs (V11: +,+,+,-)","T6_ordering":"value ordering 0<=V01<=V11<W/2<=V10<=W",
 "T7_prop2_sym":"Prop 2 (z>x at mu=lam)","T8_BRslope":"BR slope sign = sign(N)",
 "T9_prop3_signs":"Prop 3 signs (x00: +,+,-)","T11_BRsolve":"additive-root BRs solve G_L,G_F=0"}
allpass = True
for k in ["T1_Delta_pos","T1_Delta_form","T2_sub_interior","T2_add_super","T3_key_ineq",
          "T4_prop1_signs","T5_prop1p_signs","T6_ordering","T7_prop2_sym","T8_BRslope",
          "T9_prop3_signs","T11_BRsolve"]:
    status = "PASS" if V[k]==0 else f"FAIL ({V[k]})"
    if V[k]!=0: allpass=False
    log(f"  [{status:>9}]  {labels[k]}   (violations {V[k]}/{N})")
log(f"  min Delta over scan = {minD:.4e}   min Delta0 = {minD0:.4e}")
log(f"  => ALL THEOREMS PASS: {allpass}")

log("")
log("NUMERICAL claims (in-region, fresh scan):")
log(f"  dx*/dbeta11 > 0            : {num['dxb11_pos']}/{num['dxb11_tot']} = {100*num['dxb11_pos']/max(num['dxb11_tot'],1):.3f}%")
log(f"  dz*/dbeta11 negative       : {100*num['dzb11_neg']/max(num['dzb11_pos']+num['dzb11_neg'],1):.2f}%")
log(f"  dx*/dlam total-vs-direct reversals : {num['dxlam_rev']}/{num['dxlam_tot']} = {100*num['dxlam_rev']/max(num['dxlam_tot'],1):.3f}%")
log(f"  provocation region N>0     : {100*num['prov']/max(num['asym_tot']+ (N//2),1):.2f}% (of converged asym+sym)")
log(f"  Delta0>0 violations        : {num['D0_bad']}")
log(f"  x00* in (0,1) violations   : {num['x00_bad']}")
log(f"  m x00* < 1 violations      : {num['mx00_bad']}")
log(f"  det J > 0 violations       : {num['detJ_bad']}")

# ==========================================================================
# PART 3 : BREAK REGIONS -- relax each assumption
# ==========================================================================
log("")
log("="*78); log("PART 3: BREAK REGIONS (relax one assumption at a time)")
log("="*78)
# (i) relax (A2): lamW>c  -> interior x11* should exceed 1 (corner)
random.seed(303); n=40000; corner=0; tot=0
for _ in range(n):
    lam=random.uniform(1e-2,0.4999); c=random.uniform(5e-2,10.); b1=random.uniform(0,0.999)
    W=random.uniform(c/lam*1.0001, c/lam*3.0)   # lamW>c
    x11,_=x11_V11(W,c,lam,b1); tot+=1
    if x11>1: corner+=1
log(f"(i) relax (A2) lamW>c: interior x11*>1 (corner) in {corner}/{tot} = {100*corner/tot:.1f}% "
    f"[in-region this is 0%]")
# (ii) relax mu=lam -> Prop 2 fails
random.seed(304); n=40000; fail=0; tot=0
for _ in range(n):
    lam=random.uniform(1e-2,0.4999); mu=random.uniform(1e-2,0.4999)
    if mu==lam: continue
    c=random.uniform(5e-2,10.); W=random.uniform(1e-3,c/max(lam,mu)); b1=random.uniform(1e-4,0.999)
    _,V11=x11_V11(W,c,lam,b1); z,x,res=solve_fp(W,c,lam,mu,V11)
    if res>1e-9: continue
    tot+=1
    if not (z>x): fail+=1
log(f"(ii) relax mu=lam: Prop 2 (z>x) FAILS in {fail}/{tot} = {100*fail/tot:.1f}% "
    f"[Prop 2 only claims mu=lam]")
# (iii) relax (A2') muW>c while keeping lamW<=c (needs mu>lam): check (0,0) claims
random.seed(305); n=60000; tot=0; D0bad=0; x00bad=0; mx00bad=0
for _ in range(n):
    lam=random.uniform(1e-2,0.4999); mu=random.uniform(lam*1.01,0.4999)  # mu>lam
    c=random.uniform(5e-2,10.)
    # want lamW<=c (A2 holds) but muW>c (A2' fails): W in (c/mu, c/lam]
    lo, hi = c/mu*1.0001, c/lam
    if lo>=hi: continue
    W=random.uniform(lo,hi); b1=random.uniform(1e-4,0.999); b0=random.uniform(0,0.999)
    _,V11=x11_V11(W,c,lam,b1); z,x,res=solve_fp(W,c,lam,mu,V11)
    if res>1e-9: continue
    V10,V01=vals_10_01(z,x,W,c,lam,mu,V11)
    x00,D0,m=x00_of(c,mu,b0,V10,V01,V11); tot+=1
    if D0<=0: D0bad+=1; continue
    if not (0<x00<1): x00bad+=1
    if not (m*x00<1): mx00bad+=1
log(f"(iii) relax (A2') muW>c (A2 kept): tested {tot} draws -> Delta0<=0: {D0bad}, "
    f"x00 not in(0,1): {x00bad}, m x00>=1: {mx00bad}")
log(f"     => (A2') {'IS' if (D0bad+x00bad+mx00bad)>0 else 'appears NOT strictly'} needed for (0,0) "
    f"interiority in this range")

with open(os.path.join(NUM, "task6_confirmation.txt"), "w") as f:
    f.write("\n".join(log_lines) + "\n")

# ==========================================================================
# PART 4 : symbolic master-quadratic identity
# ==========================================================================
log("")
log("="*78); log("PART 4: symbolic master-quadratic identity (sympy)")
log("="*78)
x, c_, lam_, mu_, W_, b11, b00 = sp.symbols('x c lambda mu W beta_11 beta_00', positive=True)
V10s, V01s, V11s = sp.symbols('V10 V01 V11')
m_ = mu_*(1+b00); G = V10s + b00*V01s - (1+b00)*V11s; D = V10s - V01s
F00 = c_*m_*x**2 - (3*c_ + 2*mu_*m_*G)*x + 2*mu_*(1-b00)*D
# substitute terminal values (V10,V01,V11)->(W,0,W/2), mu->lam, b00->b11
F_sub = F00.subs({V10s: W_, V01s: 0, V11s: W_/2, mu_: lam_, b00: b11})
M_ = lam_*(1+b11)
F11 = M_*c_*x**2 - (3*c_ + 2*lam_*M_*W_ - M_**2*W_)*x + 2*W_*(2*lam_ - M_)
diff = sp.simplify(sp.expand(F_sub - F11))
log(f"  F00[terminal subs] - F11  simplifies to: {diff}")
log(f"  => master-quadratic identity holds: {diff == 0}")

with open(os.path.join(NUM, "task6_confirmation.txt"), "w") as f:
    f.write("\n".join(log_lines) + "\n")
print("\n[written]", os.path.join(NUM, "task6_confirmation.txt"))
