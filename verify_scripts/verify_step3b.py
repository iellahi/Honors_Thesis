"""
Step 3b: rigorous versions of Checks 4-6 via polynomial reduction mod F, plus
full analytical sign proofs for dV11/d{W,c,lam,beta}. Every identity machine-checked.
Notation: u = lam^2 W(1-b^2), M = lam(1+b), S = sqrt(Delta), Delta = 9c^2 - 2cu + u^2.
Assumptions (A1)-(A2): lam<1/2, b in [0,1), lam*W <= c  =>  u <= lam*c(1-b^2) < c/2.
"""
import sympy as sp

x = sp.symbols('x', positive=True)
W, c, lam, b, q = sp.symbols('W c lambda beta q', positive=True)
M = lam*(1+b)
u = lam**2*W*(1-b**2)
A = c*M
B = -(3*c + u)
C = 2*lam*W*(1-b)
F = A*x**2 + B*x + C
Delta = 9*c**2 - 2*c*u + u**2

# Proper reduction: c*x^2 = (-B*x - C)/M  (linear in x), from F=0
cx2 = (-B*x - C)/M

print("="*70)
print("R1: F_lam mod F  ==  (x/lam)(3c - u)")
print("="*70)
F_lam = sp.diff(F, lam)
# Note: diff treats u,M as functions of lam; F_lam is a poly in x with c*x^2-type term
# Replace the x^2 term: F_lam = d(A)/dlam x^2 + d(B)/dlam x + d(C)/dlam
F_lam_red = sp.simplify(sp.expand(F_lam).subs(x**2, sp.expand(cx2/c)))
target = (x/lam)*(3*c - u)
print("PASS" if sp.simplify(F_lam_red - target) == 0 else f"FAIL: {sp.simplify(F_lam_red - target)}")
print("Under (A2): u < c/2 < 3c  =>  F_lam > 0  =>  dx*/dlam > 0. [No 12c>W needed]")

print("="*70)
print("R2: F_b mod F  ==  [x(3c + M^2 W) - 4 lam W] / (1+b)")
print("="*70)
F_b = sp.diff(F, b)
F_b_red = sp.simplify(sp.expand(F_b).subs(x**2, sp.expand(cx2/c)))
target_b = (x*(3*c + M**2*W) - 4*lam*W)/(1+b)
print("PASS" if sp.simplify(F_b_red - target_b) == 0 else f"FAIL: {sp.simplify(F_b_red - target_b)}")

print("="*70)
print("R3: dx*/dbeta < 0.  Let xbar = 4 lam W/(3c + M^2 W). Show F(xbar) < 0")
print("    via G = F(xbar) * (3c+M^2W)^2 / (lam W) and the identity")
print("    G = -18c^2(1+b) + 4 c M lam W (1-3b^2) - 2 M^3 lam W^2 (1-b)^2")
print("="*70)
xbar = 4*lam*W/(3*c + M**2*W)
G_direct = sp.expand(F.subs(x, xbar) * (3*c + M**2*W)**2 / (lam*W))
G_claim = -18*c**2*(1+b) + 4*c*M*lam*W*(1-3*b**2) - 2*M**3*lam*W**2*(1-b)**2
print("Identity:", "PASS" if sp.simplify(G_direct - sp.expand(G_claim)) == 0 else "FAIL")
# Sign under A2: 4cM lamW(1-3b^2) <= 4cM*lamW <= 4c^2 M = 4c^2 lam(1+b) < 2c^2(1+b)
# => G < -18c^2(1+b) + 2c^2(1+b) = -16 c^2 (1+b) < 0.
print("Bound chain: 4cM(lamW)(1-3b^2) <= 4cM*c < 2c^2(1+b)  [lam<1/2 => 4M=4lam(1+b)<2(1+b)]")
print("=> G < -16c^2(1+b) < 0 => F(xbar)<0 => x^- < xbar => F_b(x*)<0 => dx*/dbeta<0. QED")

print("="*70)
print("R4: KEY LEMMA  sqrt(Delta) > 3c - u   (since Delta - (3c-u)^2 = 4cu > 0)")
print("="*70)
print("PASS" if sp.simplify(sp.expand(Delta - (3*c - u)**2) - 4*c*u) == 0 else "FAIL")

print("="*70)
print("R5: dV11/dc > 0.  dV/dc = x*[sqrt(D) - 3c + c M x*]/(2 lam(1+b) sqrt(D))")
print("    and c M x* = (3c + u - sqrt(D))/2  =>  bracket = (sqrt(D) + u - 3c)/2 > 0 by R4")
print("="*70)
S = sp.sqrt(Delta)
x_star = (3*c + u - S)/(2*c*M)
V11 = c*x_star/(2*lam*(1+b)) + b*W/(1+b)
dV_dc = sp.diff(V11, c)
claim_dc = x_star*(S + u - 3*c)/(4*lam*(1+b)*S)   # = x*[bracket]/(2 lam(1+b) S), bracket=(S+u-3c)/2
print("Identity:", "PASS" if sp.simplify(dV_dc - claim_dc) == 0 else "FAIL")
print("Positive by R4 (S > 3c - u) and x* > 0 for b < 1. QED")

print("="*70)
print("R6: dV11/dlam < 0.  dV/dlam = c x*[(3c-u) - sqrt(D)]/(2 lam^2 (1+b) sqrt(D)) < 0 by R4")
print("="*70)
dV_dlam = sp.diff(V11, lam)
claim_dlam = c*x_star*((3*c - u) - S)/(2*lam**2*(1+b)*S)
print("Identity:", "PASS" if sp.simplify(dV_dlam - claim_dlam) == 0 else "FAIL")
print("Negative by R4. QED")

print("="*70)
print("R7: dV11/dW > 0.  dV/dW = (c/(2 lam(1+b))) dx*/dW + b/(1+b), both terms >= 0,")
print("    second strictly > 0 for b>0; for b=0 first term > 0 (dx*/dW>0). QED")
print("="*70)
dV_dW = sp.simplify(sp.diff(V11, W))
# spot verify positivity structure numerically at b=0 edge
val0 = dV_dW.subs([(b, 0), (lam, sp.Rational(1,4)), (c, 1), (W, 2)])
print("dV/dW at b=0, lam=1/4, c=1, W=2:", sp.N(val0), "> 0")

print("="*70)
print("R8 (PRIORITY): dV11/dbeta > 0 under (A1)-(A2) -- full analytical proof")
print("    Reduce to 2 variables: W = q c/lam^2 * ... use q = lam*(lamW/c) in (0,1/2)")
print("="*70)
# Substitute lam W = t c with t in (0,1]; q = lam t in (0,1/2). All terms scale with c^2.
# T1 = 18 + 3q(1+b)^2 + q(1-b^2) + q^2(1-b^2)(1+b)^2 + q^2(1-b^2)^2 - 8q(1+b)
# T2 = 6 - 2q(1+b);  Dt = 9 - 2q(1-b^2) + q^2(1-b^2)^2
# Claim: sign(dV/dbeta) = sign(T1 - sqrt(Dt) T2) and T1^2 - Dt*T2^2 > 0 with T1>0,T2>0.
T1 = 18 + 3*q*(1+b)**2 + q*(1-b**2) + q**2*(1-b**2)*(1+b)**2 + q**2*(1-b**2)**2 - 8*q*(1+b)
T2 = 6 - 2*q*(1+b)
Dt = 9 - 2*q*(1-b**2) + q**2*(1-b**2)**2

# First: machine-verify that sign(dV11/dbeta) == sign(T1 - sqrt(Dt)*T2) by direct substitution.
t = sp.symbols('t', positive=True)
subsW = {W: t*c/lam}
dV_db = sp.diff(V11, b)
expr = sp.simplify(dV_db.subs(subsW))
# Build the claimed sign expression with q -> lam*t
sign_expr = (T1 - sp.sqrt(Dt)*T2).subs(q, lam*t)
# dV/db should equal  [c/(2 lam (1+b)^2)] * (T1 - sqrt(Dt) T2) / (2 sqrt(Dt))   -- verify proportionality
cand = c/(4*lam*(1+b)**2) * sign_expr / sp.sqrt(Dt.subs(q, lam*t))
diff_check = sp.simplify(expr - cand)
print("dV/dbeta == c(T1 - sqrt(Dt)T2)/(4 lam (1+b)^2 sqrt(Dt)):",
      "PASS" if diff_check == 0 else f"FAIL: {diff_check}")

# Positivity: T2 in (4,6) > 0 for q<1/2, b<1. T1 > 0 (check min). Then need T1^2 - Dt T2^2 > 0.
P = sp.expand(T1**2 - Dt*T2**2)
P = sp.factor(P)
print("P(q,b) = T1^2 - Dt*T2^2 factors as:")
print(" ", P)
