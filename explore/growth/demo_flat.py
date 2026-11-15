"""Demonstration: unfold a flat point of real charges hierarchically and count. usage: demo_flat.py config.json rho_top sigma out_prefix"""
import numpy as np, json, sys, time
from fractions import Fraction as Fr
from flint import arb, fmpq, ctx
from harm import monos, Poly
from realize_jet import realize, jetvec, legcoef, pmul, r2pow
from flatsolve import top_poly
import finder
ctx.prec=400
cfg=json.load(open(sys.argv[1])); rho_top=Fr(sys.argv[2]); sig=Fr(sys.argv[3]); pref=sys.argv[4]
a0=np.array(cfg['sites']); q0=np.array(cfg['charges']); d=cfg['d']; n=len(q0)
def zon_exact(l, axis):
    lc=legcoef(l); out={}
    for k,ck in lc.items():
        zp={(0,0,l-2*k):Fr(1)} if axis=='z' else {(l-2*k,0,0):Fr(1)}
        for m,v in pmul(zp,r2pow(k)).items(): out[m]=out.get(m,0)+ck*v
    return out
# top polynomial of the base (degree d) for scale
ctop=top_poly(a0,q0,d); h=Fr(float(np.abs(ctop).max())).limit_denominator(10**6)
# levels: A_{d-1} about z, then alternate
coef={}; c=rho_top*h; rho=rho_top; axis='z'
levels=[]
for l in range(d-1,0,-1):
    levels.append((l,axis,c)); 
    rho=rho*sig; c=c*rho; axis='x' if axis=='z' else 'z'
target=[]; J={}
for l,ax,cl in levels:
    Z=zon_exact(l,ax); mx=max(abs(v) for v in Z.values())
    for m,v in Z.items(): J[m]=cl*v/mx
for l in range(1,d):
    L,_=monos(l); target+=[J.get(m,Fr(0)) for m in L if sum(m)==l]
print("levels (degree, axis, coefficient):",[(l,ax,float(cl)) for l,ax,cl in levels],flush=True)
t=time.time(); xa=realize(a0,q0,d,target); print(f"realized in {time.time()-t:.0f}s",flush=True)
sites=[[Fr(xa[3*i+k].str(100,radius=False)) for k in range(3)] for i in range(n)]
qs=[Fr(xa[3*n+i].exp().str(100,radius=False)) for i in range(n)]
# exact jet check of the rounded rational config
S=[[arb(fmpq(v.numerator,v.denominator)) for v in s] for s in sites]; Q=[arb(fmpq(v.numerator,v.denominator)) for v in qs]
jv=jetvec(S,Q,d,arb(1)); err=max(abs(float((jv[k]-arb(fmpq(target[k].numerator,target[k].denominator))).mid())) for k in range(len(target)))
print(f"rational config: max jet error {err:.2e}; smallest level coefficient {float(levels[-1][2]):.2e}; charges {[round(float(v),4) for v in qs]}",flush=True)
json.dump({"name":pref,"sites":[[str(v) for v in s] for s in sites],"charges":[str(v) for v in qs],"levels":[(l,ax,str(cl)) for l,ax,cl in levels]},open(f"configs/{pref}.json","w"))
if len(sys.argv)>5 and sys.argv[5]=='nofinder': sys.exit(0)
af=np.array([[float(v) for v in s] for s in sites]); qf=np.array([float(v) for v in qs])
t=time.time(); r,pts,idx=finder.report(af,qf,hp=True,a_exact=S,q_exact=Q)
print(f"finder: {r} ({time.time()-t:.0f}s)",flush=True)
rn=np.linalg.norm(pts,axis=1); print("radii of points near the flat point:",np.sort(rn[rn<0.2]))
np.save(f"out/{pref}_pts.npy",np.c_[pts,idx])
