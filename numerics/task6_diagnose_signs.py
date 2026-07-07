"""
Diagnostic: are the Prop 1' / Prop 3 sign "failures" real, or finite-difference
artifacts near near-zero derivatives?  Re-check every float-flagged case with
50-digit mpmath central differences.  A theorem violation only counts if mpmath
(high precision) also reports the wrong sign.
"""
import math, random
import mpmath as mp
mp.mp.dps = 50

# float solvers (same as task6)
def x11_V11(W, c, lam, b1):
    u = lam*lam*W*(1-b1*b1); D = 9*c*c-2*c*u+u*u
    x = (3*c+u-math.sqrt(D))/(2*c*lam*(1+b1))
    return x, c*x/(2*lam*(1+b1)) + b1*W/(1+b1)
def br_x(z,W,c,lam,mu,V):
    r=c*c*lam*lam*z*z+2*c*mu*mu*lam*z*(1-lam*z)**2*V; return (-c*lam*z+math.sqrt(r))/(c*mu*(1-lam*z))
def br_z(x,W,c,lam,mu,V):
    r=c*c*mu*mu*x*x+2*c*mu*lam*lam*x*(1-mu*x)*(W-V); return (-c*mu*x+math.sqrt(r))/(c*lam*(1-mu*x))
def solve_fp(W,c,lam,mu,V,tol=1e-13,itmax=5000):
    z,x=0.5,0.5
    for _ in range(itmax):
        zn=br_z(x,W,c,lam,mu,V); xn=br_x(zn,W,c,lam,mu,V)
        if abs(zn-z)<tol and abs(xn-x)<tol: z,x=zn,xn; break
        z,x=zn,xn
    GL=(c*lam/2)*(1-mu*x)*z*z+c*mu*x*z-lam*mu*x*(W-V)
    GF=(c*mu/2)*(1-lam*z)*x*x+c*lam*z*x-mu*lam*z*(1-lam*z)*V
    return z,x,max(abs(GL),abs(GF))
def vals(z,x,W,c,lam,mu,V):
    PL,PF=lam*z,mu*x; S=PL+PF-PL*PF
    return (PL*W+(1-PL)*PF*V-c*z*z/2)/S, ((1-PL)*PF*V-c*x*x/2)/S
def x00_f(c,mu,b0,V10,V01,V11):
    m=mu*(1+b0); G=V10+b0*V01-(1+b0)*V11; D=V10-V01
    A,B,C=c*m,-(3*c+2*mu*m*G),2*mu*(1-b0)*D
    return (-B-math.sqrt(B*B-4*A*C))/(2*A)
def sgn(v): return 0 if v==0 else (1 if v>0 else -1)

# mpmath high-precision versions -- ALL arithmetic in mpf; perturbation added in mpf.
def V11_mp_raw(W,c,lam,b1):     # inputs are mpf
    u=lam*lam*W*(1-b1*b1); D=9*c*c-2*c*u+u*u
    x=(3*c+u-mp.sqrt(D))/(2*c*lam*(1+b1))
    return c*x/(2*lam*(1+b1))+b1*W/(1+b1)
def x00_mp_raw(c,mu,b0,V10,V01,V11):    # inputs are mpf
    m=mu*(1+b0); G=V10+b0*V01-(1+b0)*V11; D=V10-V01
    A,B,C=c*m,-(3*c+2*mu*m*G),2*mu*(1-b0)*D
    return (-B-mp.sqrt(B*B-4*A*C))/(2*A)

H = mp.mpf('1e-20')             # mpf step (base kept in mpf so it survives)
def hp_sign_V11(W,c,lam,b1,which):
    v=[mp.mpf(W),mp.mpf(c),mp.mpf(lam),mp.mpf(b1)]
    idx={'W':0,'c':1,'lam':2,'b1':3}[which]
    vp=v[:]; vp[idx]=vp[idx]+H; fp=V11_mp_raw(*vp)
    vm=v[:]; vm[idx]=vm[idx]-H; fm=V11_mp_raw(*vm)
    d=fp-fm
    return 0 if d==0 else (1 if d>0 else -1)
def hp_sign_x00(c,mu,b0,V10,V01,V11,which):
    cc,mm,bb=mp.mpf(c),mp.mpf(mu),mp.mpf(b0)
    v10,v01,v11=mp.mpf(V10),mp.mpf(V01),mp.mpf(V11)
    if which=='V10': fp,fm=x00_mp_raw(cc,mm,bb,v10+H,v01,v11), x00_mp_raw(cc,mm,bb,v10-H,v01,v11)
    if which=='V11': fp,fm=x00_mp_raw(cc,mm,bb,v10,v01,v11+H), x00_mp_raw(cc,mm,bb,v10,v01,v11-H)
    if which=='V01': fp,fm=x00_mp_raw(cc,mm,bb,v10,v01+H,v11), x00_mp_raw(cc,mm,bb,v10,v01-H,v11)
    d=fp-fm
    return 0 if d==0 else (1 if d>0 else -1)

random.seed(202); N=60000
p1p_flag=p1p_realbad=0; p3_flag=p3_realbad=0
p1p_examples=[]; p3_examples=[]
for i in range(N):
    lam=random.uniform(1e-2,0.4999); mu=random.uniform(1e-2,0.4999) if i%2 else lam
    b1=random.uniform(1e-4,0.9990); b0=random.uniform(0.0,0.9990)
    c=random.uniform(5e-2,10.0); W=random.uniform(1e-3,c/max(lam,mu))
    x11,V11=x11_V11(W,c,lam,b1)
    hW,hL,hc,hb=1e-6*W,1e-6*lam,1e-6*c,1e-6*max(b1,1e-3)
    # Prop 1' float signs (b1,W,c,lam) -> (+,+,+,-)
    sv=(sgn((x11_V11(W,c,lam,b1+hb)[1]-V11) if b1+hb<1 else 1), sgn(x11_V11(W+hW,c,lam,b1)[1]-V11),
        sgn(x11_V11(W,c+hc,lam,b1)[1]-V11), sgn(x11_V11(W,c,lam+hL,b1)[1]-V11))
    if sv!=(1,1,1,-1):
        p1p_flag+=1
        hp=(hp_sign_V11(W,c,lam,b1,'b1') if b1+hb<1 else 1, hp_sign_V11(W,c,lam,b1,'W'),
            hp_sign_V11(W,c,lam,b1,'c'), hp_sign_V11(W,c,lam,b1,'lam'))
        # real violation only if a hp sign is strictly wrong (nonzero and mismatched)
        exp=(1,1,1,-1); bad=any(hp[k]!=0 and hp[k]!=exp[k] for k in range(4))
        if bad: p1p_realbad+=1;
        if len(p1p_examples)<6: p1p_examples.append((round(lam,4),round(mu,4),round(c,4),round(W,4),round(b1,4),"float",sv,"hp",hp))
    # asymmetric for Prop 3
    z,x,res=solve_fp(W,c,lam,mu,V11)
    if res>1e-9: continue
    V10,V01=vals(z,x,W,c,lam,mu,V11)
    x00=x00_f(c,mu,b0,V10,V01,V11)
    e=1e-6
    s=(sgn(x00_f(c,mu,b0,V10+e,V01,V11)-x00), sgn(x00_f(c,mu,b0,V10,V01,V11+e)-x00),
       sgn(x00_f(c,mu,b0,V10,V01+e,V11)-x00))
    if s!=(1,1,-1):
        p3_flag+=1
        hp=(hp_sign_x00(c,mu,b0,V10,V01,V11,'V10'), hp_sign_x00(c,mu,b0,V10,V01,V11,'V11'),
            hp_sign_x00(c,mu,b0,V10,V01,V11,'V01'))
        exp=(1,1,-1); bad=any(hp[k]!=0 and hp[k]!=exp[k] for k in range(3))
        if bad: p3_realbad+=1
        if len(p3_examples)<6: p3_examples.append((round(lam,4),round(mu,4),round(c,4),round(W,4),"float",s,"hp",hp))

print(f"Prop 1' (V11 signs): float-flagged {p1p_flag}, HIGH-PRECISION real violations {p1p_realbad}")
for ex in p1p_examples: print("   ", ex)
print(f"Prop 3  (x00 signs): float-flagged {p3_flag}, HIGH-PRECISION real violations {p3_realbad}")
for ex in p3_examples: print("   ", ex)
print()
print("Interpretation: float flags with hp sign == 0 (or matching expected) are finite-difference")
print("artifacts where the true derivative is ~0 or below float resolution, NOT theorem violations.")
