"""Step 4c: small-mu limit of x00* and the equilibrium V00 expression."""
import sympy as sp
x = sp.symbols('x', positive=True)
c, mu, b0 = sp.symbols('c mu beta_00', positive=True)
V10, V01, V11 = sp.symbols('V_10 V_01 V_11', positive=True)
m = mu*(1+b0); G = V10 + b0*V01 - (1+b0)*V11; D = V10 - V01
A0 = c*m; B0 = -(3*c + 2*mu*m*G); C0 = 2*mu*(1-b0)*D
Delta0 = B0**2 - 4*A0*C0
x00 = (-B0 - sp.sqrt(Delta0))/(2*A0)   # subtracted root

print("L1: leading-order small-mu expansion of x00*")
ser = sp.series(x00, mu, 0, 2).removeO()
print("  x00* =", sp.simplify(ser), "+ O(mu^2)")
target = 2*mu*(1-b0)*D/(3*c)
print("  equals 2 mu (1-b0)(V10-V01)/(3c) at leading order:",
      "PASS" if sp.simplify(sp.expand(ser) - sp.expand(target)) == 0 else
      f"CHECK: {sp.simplify(ser - target)}")

print("L2: equilibrium V00 from the symmetric FOC (compact closed form)")
V00s = sp.symbols('V_00')
P = m*x
foc = -c*x + mu*(1-P)*(V10 + b0*V01) - mu*P*(b0*V10 + V01) + m*P*V11 - m*(1-P)*V00s
V00_foc = sp.solve(sp.Eq(foc, 0), V00s)[0]
print("  V00(x) =", sp.simplify(V00_foc))
# numeric consistency with Bellman-form at the root
import random
random.seed(4)
V00_bell = (P*(1-P)*(V10+V01) + P**2*V11 - sp.Rational(1,2)*c*x**2)/(P*(2-P))
f_root = sp.lambdify((c, mu, b0, V10, V01, V11), x00, 'math')
f_foc  = sp.lambdify((x, c, mu, b0, V10, V01, V11), V00_foc, 'math')
f_bell = sp.lambdify((x, c, mu, b0, V10, V01, V11), V00_bell, 'math')
worst = 0
for _ in range(5000):
    muv = random.uniform(1e-3, .4999); b0v = random.uniform(0, .9999)
    cv = random.uniform(1e-2, 10); Wv = random.uniform(1e-3, cv/muv)
    v01 = random.uniform(0, Wv/2); v11 = random.uniform(v01, Wv/2); v10 = random.uniform(v11, Wv)
    xr = f_root(cv, muv, b0v, v10, v01, v11)
    worst = max(worst, abs(f_foc(xr, cv, muv, b0v, v10, v01, v11)
                          - f_bell(xr, cv, muv, b0v, v10, v01, v11)))
print(f"  max |V00_FOC - V00_Bellman| at root over 5000 draws: {worst:.2e}  (expect ~0)")
