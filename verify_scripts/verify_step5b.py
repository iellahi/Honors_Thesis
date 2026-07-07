"""
Step 5B: numerical closure of the full backward-induction system.
For random params under (A1), (A2) lam W <= c, (A2') mu W <= c:
  1. V11, x11* from Step 2 closed forms.
  2. Fixed point (z*, x*) of the BR quadratics (iteration + residual check).
  3. Checks: z* > x* (mu=lam AND mu!=lam); ordering 0<=V01<=V11<=V10<=W;
     z*,x* in (0,1); probabilities < 1.
  4. (0,0) closure: Delta0 > 0, x00* in (0,1), P00 < 1.  [Step 4 deferred items]
  5. Transmission: d x*/d beta11 > 0? d z*/d beta11 sign?  (finite differences, total)
  6. Luke's total-effect decomposition for c and lam on z*, x*:
     total = direct (V11 frozen) + indirect (via V11). Report sign patterns.
"""
import random, math

def V11_and_x11(W, c, lam, b1):
    u = lam*lam*W*(1-b1*b1)
    Delta = 9*c*c - 2*c*u + u*u
    x11 = (3*c + u - math.sqrt(Delta)) / (2*c*lam*(1+b1))
    V11 = c*x11/(2*lam*(1+b1)) + b1*W/(1+b1)
    return V11, x11

def br_x(z, W, c, lam, mu, V11):   # laggard BR (additive root)
    rad = c*c*lam*lam*z*z + 2*c*mu*mu*lam*z*(1-lam*z)**2*V11
    return (-c*lam*z + math.sqrt(rad)) / (c*mu*(1-lam*z))

def br_z(x, W, c, lam, mu, V11):   # leader BR (additive root)
    rad = c*c*mu*mu*x*x + 2*c*mu*lam*lam*x*(1-mu*x)*(W-V11)
    return (-c*mu*x + math.sqrt(rad)) / (c*lam*(1-mu*x))

def solve_fp(W, c, lam, mu, V11, tol=1e-13, itmax=5000):
    z, x = 0.5, 0.5
    for _ in range(itmax):
        zn = br_z(x, W, c, lam, mu, V11)
        xn = br_x(zn, W, c, lam, mu, V11)
        if abs(zn-z) < tol and abs(xn-x) < tol:
            z, x = zn, xn
            break
        z, x = zn, xn
    # residuals of the quadratics
    G_L = (c*lam/2)*(1-mu*x)*z*z + c*mu*x*z - lam*mu*x*(W-V11)
    G_F = (c*mu/2)*(1-lam*z)*x*x + c*lam*z*x - mu*lam*z*(1-lam*z)*V11
    return z, x, max(abs(G_L), abs(G_F))

def values_10_01(z, x, W, c, lam, mu, V11):
    PL, PF = lam*z, mu*x
    S = PL + PF - PL*PF
    V10 = (PL*W + (1-PL)*PF*V11 - c*z*z/2)/S
    V01 = ((1-PL)*PF*V11 - c*x*x/2)/S
    return V10, V01

def x00_and_D0(W, c, mu, b0, V10, V01, V11):
    m = mu*(1+b0)
    G = V10 + b0*V01 - (1+b0)*V11
    D = V10 - V01
    A0, B0, C0 = c*m, -(3*c + 2*mu*m*G), 2*mu*(1-b0)*D
    D0 = B0*B0 - 4*A0*C0
    x00 = (-B0 - math.sqrt(D0))/(2*A0) if D0 > 0 else float('nan')
    return x00, D0, m

def full_solve(W, c, lam, mu, b1):
    V11, x11 = V11_and_x11(W, c, lam, b1)
    z, x, res = solve_fp(W, c, lam, mu, V11)
    V10, V01 = values_10_01(z, x, W, c, lam, mu, V11)
    return V11, x11, z, x, res, V10, V01

random.seed(11)
N = 60000
bad = {k: 0 for k in ["res", "z_le_x_sym", "z_le_x_asym", "order", "interior_zx",
                      "prob", "D0", "x00_int", "trans_x", "V11_half"]}
mnD0 = 1e18; trans_z_pos = trans_z_neg = 0
tot = {k: [0, 0] for k in ["dz_dc", "dx_dc", "dz_dlam", "dx_dlam"]}  # [neg, pos] counts
dir_ind_flip = {k: 0 for k in ["dz_dc", "dx_dc", "dz_dlam", "dx_dlam"]}  # direct vs total sign flips

for i in range(N):
    lam_ = random.uniform(1e-2, 0.4999)
    mu_  = random.uniform(1e-2, 0.4999) if i % 2 else lam_   # half sym, half asym
    b1_  = random.uniform(1e-4, 0.9995)
    b0_  = random.uniform(0, 0.9995)
    c_   = random.uniform(5e-2, 10.0)
    W_   = random.uniform(1e-3, c_/max(lam_, mu_))            # A2 and A2'

    V11, x11, z, x, res, V10, V01 = full_solve(W_, c_, lam_, mu_, b1_)
    if res > 1e-9: bad["res"] += 1; continue
    if not (V11 < W_/2): bad["V11_half"] += 1
    if mu_ == lam_ and not (z > x): bad["z_le_x_sym"] += 1
    if mu_ != lam_ and not (z > x): bad["z_le_x_asym"] += 1
    if not (0 <= V01 <= V11 <= V10 <= W_ + 1e-12): bad["order"] += 1
    if not (0 < z < 1 and 0 <= x < 1): bad["interior_zx"] += 1
    if not (lam_*z < 1 and mu_*x < 1): bad["prob"] += 1

    x00, D0, m = x00_and_D0(W_, c_, mu_, b0_, V10, V01, V11)
    mnD0 = min(mnD0, D0)
    if D0 <= 0: bad["D0"] += 1
    elif not (0 < x00 < 1 and m*x00 < 1): bad["x00_int"] += 1

    # transmission wrt beta11 (total, finite difference)
    h = 1e-6 * max(b1_, 1e-3)
    if b1_ + h < 1:
        V11b, _, zb, xb, resb, _, _ = full_solve(W_, c_, lam_, mu_, b1_ + h)
        if resb < 1e-9:
            if (xb - x)/h <= 0: bad["trans_x"] += 1
            if (zb - z)/h > 0: trans_z_pos += 1
            else: trans_z_neg += 1

    # total vs direct effects for c and lam (subsample for speed)
    if i % 20 == 0:
        for par, key_z, key_x in [("c", "dz_dc", "dx_dc"), ("lam", "dz_dlam", "dx_dlam")]:
            if par == "c":
                h2 = 1e-6*c_
                V11p, _, zp, xp, rp, _, _ = full_solve(W_, c_+h2, lam_, mu_, b1_)
                zd, xd, rd = solve_fp(W_, c_+h2, lam_, mu_, V11)  # direct: V11 frozen
            else:
                h2 = 1e-6*lam_
                if lam_ + h2 >= 0.5: continue
                V11p, _, zp, xp, rp, _, _ = full_solve(W_, c_, lam_+h2, mu_, b1_)
                zd, xd, rd = solve_fp(W_, c_, lam_+h2, mu_, V11)
            if rp > 1e-9 or rd > 1e-9: continue
            tz, tx = (zp - z)/h2, (xp - x)/h2
            dz_dir, dx_dir = (zd - z)/h2, (xd - x)/h2
            tot[key_z][0 if tz < 0 else 1] += 1
            tot[key_x][0 if tx < 0 else 1] += 1
            if (tz < 0) != (dz_dir < 0): dir_ind_flip[key_z] += 1
            if (tx < 0) != (dx_dir < 0): dir_ind_flip[key_x] += 1

print(f"N = {N} parameter draws (half mu=lam, half mu!=lam)")
print("Violation counts:", bad)
print(f"min Delta0 over scan: {mnD0:.6e}")
print(f"d z*/d beta11 sign split:  positive {trans_z_pos}, negative {trans_z_neg}")
print("Total-effect sign counts [neg, pos]:")
for k, v in tot.items():
    print(f"  {k}: {v}   (direct-vs-total sign flips: {dir_ind_flip[k]})")
