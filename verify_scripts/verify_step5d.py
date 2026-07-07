"""Step 5d: corrected laggard-BR slope analysis + what survives.
X'(z) sign = sign(N), N = mu x (1-lam z)^2 - 2 lam^2 z^2  (non-monotone BR).
Transmission dx*/dV11 > 0: direct channel proven; full equilibrium sign NUMERICAL.
det J sign checked numerically (uniqueness/stability support)."""
import sympy as sp

z, x = sp.symbols('z x', positive=True)
c, lam, mu, W, V11 = sp.symbols('c lambda mu W V_11', positive=True)
G_F = (c*mu/2)*(1-lam*z)*x**2 + c*lam*z*x - mu*lam*z*(1-lam*z)*V11
G_L = (c*lam/2)*(1-mu*x)*z**2 + c*mu*x*z - lam*mu*x*(W - V11)

print("D1: corrected identity for G_F_z on the BR")
V11_on = sp.solve(sp.Eq(G_F, 0), V11)[0]
G_F_z_on = sp.together(sp.simplify(sp.diff(G_F, z).subs(V11, V11_on)))
N = mu*x*(1-lam*z)**2 - 2*lam**2*z**2
cand = -c*x*N/(2*z*(1-lam*z))
print("  G_F_z|BR == -c x N / (2z(1-lam z)):",
      "PASS" if sp.simplify(G_F_z_on - cand) == 0 else "FAIL")
print("  => sign(X'(z)) = -sign(G_F_z) = sign(N): X increasing iff mu x(1-lam z)^2 > 2 lam^2 z^2")

print("D2: numerical -- det J > 0 and dx*/dV11 > 0 at equilibrium across the region")
import random, math
def V11f(Wv, cv, l, b1):
    u = l*l*Wv*(1-b1*b1); D = 9*cv*cv - 2*cv*u + u*u
    xs = (3*cv + u - math.sqrt(D))/(2*cv*l*(1+b1))
    return cv*xs/(2*l*(1+b1)) + b1*Wv/(1+b1)
def brx(zv, Wv, cv, l, m, V):
    r = cv*cv*l*l*zv*zv + 2*cv*m*m*l*zv*(1-l*zv)**2*V
    return (-cv*l*zv + math.sqrt(r))/(cv*m*(1-l*zv))
def brz(xv, Wv, cv, l, m, V):
    r = cv*cv*m*m*xv*xv + 2*cv*m*l*l*xv*(1-m*xv)*(Wv-V)
    return (-cv*m*xv + math.sqrt(r))/(cv*l*(1-m*xv))
def fp(Wv, cv, l, m, V):
    zv = xv = 0.5
    for _ in range(4000):
        zn = brz(xv, Wv, cv, l, m, V); xn = brx(zn, Wv, cv, l, m, V)
        if abs(zn-zv) < 1e-13 and abs(xn-xv) < 1e-13: return zn, xn
        zv, xv = zn, xn
    return zv, xv

fGLz = sp.lambdify((z,x,c,lam,mu,W,V11), sp.diff(G_L,z), 'math')
fGLx = sp.lambdify((z,x,c,lam,mu,W,V11), sp.diff(G_L,x), 'math')
fGFz = sp.lambdify((z,x,c,lam,mu,W,V11), sp.diff(G_F,z), 'math')
fGFx = sp.lambdify((z,x,c,lam,mu,W,V11), sp.diff(G_F,x), 'math')
fGLV = sp.lambdify((z,x,c,lam,mu,W,V11), sp.diff(G_L,V11), 'math')
fGFV = sp.lambdify((z,x,c,lam,mu,W,V11), sp.diff(G_F,V11), 'math')

random.seed(21)
badJ = badT = 0; n = 40000; Xup = Xdown = 0
for _ in range(n):
    l = random.uniform(1e-2,.4999); m = random.uniform(1e-2,.4999)
    cv = random.uniform(5e-2,10.); Wv = random.uniform(1e-3, cv/max(l,m))
    b1 = random.uniform(1e-4,.9995)
    V = V11f(Wv,cv,l,b1)
    zv,xv = fp(Wv,cv,l,m,V)
    a11,a12 = fGLz(zv,xv,cv,l,m,Wv,V), fGLx(zv,xv,cv,l,m,Wv,V)
    a21,a22 = fGFz(zv,xv,cv,l,m,Wv,V), fGFx(zv,xv,cv,l,m,Wv,V)
    det = a11*a22 - a12*a21
    if det <= 0: badJ += 1
    dxdV = -(a11*fGFV(zv,xv,cv,l,m,Wv,V) - a21*fGLV(zv,xv,cv,l,m,Wv,V))/det
    if dxdV <= 0: badT += 1
    if m*xv*(1-l*zv)**2 > 2*l*l*zv*zv: Xup += 1
    else: Xdown += 1
print(f"  n={n}: det J <= 0: {badJ};  dx*/dV11 <= 0: {badT}")
print(f"  equilibria in X-increasing region: {Xup}, X-decreasing region: {Xdown}")
