"""Step 6: (0,0) pseudo-statics for Prop 3 + frequency check of the V10-condition.
F00 = cm x^2 - (3c + 2 mu m G)x + 2 mu(1-b0)D,  G = V10 + b0 V01 - (1+b0)V11, D = V10-V01.
F_x(x00*) = -sqrt(Delta0) < 0  =>  sign(dx00*/dtheta) = sign(F_theta)."""
import sympy as sp
x = sp.symbols('x', positive=True)
c, mu, b0 = sp.symbols('c mu beta_00', positive=True)
V10, V01, V11 = sp.symbols('V_10 V_01 V_11', positive=True)
m = mu*(1+b0); G = V10 + b0*V01 - (1+b0)*V11; D = V10 - V01
F = c*m*x**2 - (3*c + 2*mu*m*G)*x + 2*mu*(1-b0)*D

print("P1: partials wrt continuation values")
for v, name in [(V10,"V10"), (V01,"V01"), (V11,"V11")]:
    print(f"  F_{name} =", sp.factor(sp.diff(F, v)))
# claims:
c1 = sp.simplify(sp.diff(F,V10) - 2*mu*((1-b0) - m*x))
c2 = sp.simplify(sp.diff(F,V01) - (-2*mu)*(b0*m*x + (1-b0)))
c3 = sp.simplify(sp.diff(F,V11) - 2*mu*m*(1+b0)*x)
print("  F_V10 == 2mu[(1-b0)-mx]:", c1==0, "| F_V01 == -2mu[b0 mx+(1-b0)]:", c2==0,
      "| F_V11 == 2mu m(1+b0)x:", c3==0)
print("  => dx00/dV11 > 0 (direct-jump channel), dx00/dV01 < 0,")
print("     dx00/dV10 > 0 iff m*x00* < 1-b0")

print("P2: frequency of m*x00* < 1-b0 at equilibrium (full closed system)")
import random, math
def V11f(W,cv,l,b1):
    u=l*l*W*(1-b1*b1); Dl=9*cv*cv-2*cv*u+u*u
    xs=(3*cv+u-math.sqrt(Dl))/(2*cv*l*(1+b1)); return cv*xs/(2*l*(1+b1))+b1*W/(1+b1)
def brx(z,W,cv,l,mm,V):
    r=cv*cv*l*l*z*z+2*cv*mm*mm*l*z*(1-l*z)**2*V; return (-cv*l*z+math.sqrt(r))/(cv*mm*(1-l*z))
def brz(xv,W,cv,l,mm,V):
    r=cv*cv*mm*mm*xv*xv+2*cv*mm*l*l*xv*(1-mm*xv)*(W-V); return (-cv*mm*xv+math.sqrt(r))/(cv*l*(1-mm*xv))
random.seed(31); n=30000; hold=0; fail=0
for _ in range(n):
    l=random.uniform(1e-2,.4999); mm=random.uniform(1e-2,.4999)
    cv=random.uniform(5e-2,10.); W=random.uniform(1e-3,cv/max(l,mm))
    b1=random.uniform(1e-4,.9995); b0v=random.uniform(0,.9995)
    V=V11f(W,cv,l,b1); z=xv=.5
    for _ in range(3000):
        zn=brz(xv,W,cv,l,mm,V); xn=brx(zn,W,cv,l,mm,V)
        if abs(zn-z)<1e-13 and abs(xn-xv)<1e-13: break
        z,xv=zn,xn
    PL,PF=l*z,mm*xv; S=PL+PF-PL*PF
    v10=(PL*W+(1-PL)*PF*V-cv*z*z/2)/S; v01=((1-PL)*PF*V-cv*xv*xv/2)/S
    mv=mm*(1+b0v); Gv=v10+b0v*v01-(1+b0v)*V; Dv=v10-v01
    D0=(3*cv+2*mm*mv*Gv)**2-8*cv*mv*mm*(1-b0v)*Dv
    x00=(3*cv+2*mm*mv*Gv-math.sqrt(D0))/(2*cv*mv)
    if mv*x00 < 1-b0v: hold+=1
    else: fail+=1
print(f"  n={n}: m*x00*<1-b0 holds {hold}, fails {fail}  ({100*fail/n:.1f}% ambiguous-sign region)")
