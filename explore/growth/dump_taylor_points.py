"""Second, separate check of the points near the flat point: Taylor polynomial (degree 9) of the actual rational
configuration at the origin, computed in arb, then real critical points of that polynomial in |x| < 0.01 by
multi-scale Newton. Does not use finder.py."""
import numpy as np, json, sys
from fractions import Fraction as Fr
from flint import arb, fmpq, ctx
from realize_jet import jet
from harm import monos, Poly, newton_real
ctx.prec=400
cfg=json.load(open(sys.argv[1])); D=int(sys.argv[2])
S=[[arb(fmpq(Fr(v).numerator,Fr(v).denominator)) for v in s] for s in cfg['sites']]
Q=[arb(fmpq(Fr(v).numerator,Fr(v).denominator)) for v in cfg['charges']]
L,ix=monos(D); c=np.zeros(len(L))
for l in range(1,D+1):
    J=jet(S,Q,l,arb(1))
    for m,v in J.items(): c[ix[m]]=float(v.mid())
P=Poly(c,D); rng=np.random.default_rng(0)
lo=float(sys.argv[3]) if len(sys.argv)>3 else -7
X=np.vstack([rng.standard_normal((4000,3))*10**rng.uniform(lo,np.log10(float(sys.argv[4]) if len(sys.argv)>4 else 0.03),(4000,1)) for _ in range(int(sys.argv[5]) if len(sys.argv)>5 else 6)])
for _ in range(4): X=newton_real(P,X,iters=50)
g,H=P.gradhess(X); nx=np.linalg.norm(X,axis=1)
st=np.linalg.norm(np.linalg.solve(H,g[...,None])[...,0],axis=1); ok=np.isfinite(nx)&(st<1e-9*nx)&(nx<float(sys.argv[4]) if len(sys.argv)>4 else nx<1e-2)
X=X[ok]; keep=[]
for i in np.argsort(np.linalg.norm(X,axis=1)):
    if all(np.linalg.norm(X[i]-X[j])>1e-6*np.linalg.norm(X[i]) for j in keep): keep.append(i)
R=X[keep]; g,H=P.gradhess(R); ev=np.linalg.eigvalsh(H)
print(f"Taylor degree {D}: {len(R)} critical points with |x|<0.01; index1 {((ev<0).sum(1)==1).sum()}, index2 {((ev<0).sum(1)==2).sum()}")
np.save(sys.argv[6], R)
