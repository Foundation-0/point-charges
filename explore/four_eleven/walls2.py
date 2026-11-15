# Scale-aware wall tracer (numerical). From each index-2 point follow the two uphill branches of the
# unstable direction with RK4 on the normalised gradient; the step is a small fraction of the distance to
# the nearest balance point or charge, so clusters at 1e-5 are resolved. Usage: walls2.py indep_result.json config.json
import json, sys, collections, numpy as np
from fractions import Fraction as Fr
import finder
res=json.load(open(sys.argv[1])); cfg=json.load(open(sys.argv[2]))
P=np.array([[float(Fr(c)) for c in s] for s in cfg['sites']]); Q=np.array([float(Fr(c)) for c in cfg['charges']])
Z=[(np.array(z['centre_float']),int(z['index'])) for z in res['zeros']]
crit=np.array([z[0] for z in Z])
def g_unit(x):
    g,_=finder.gH(x,P,Q); return g/np.linalg.norm(g)
def climb(x0,v,sign,frac=0.02):
    others=np.array([c for c in crit if np.linalg.norm(c-x0)>0])
    dcrit=np.min(np.linalg.norm(others-x0,axis=1))
    x=x0+sign*v*dcrit*1e-3
    for it in range(200000):
        dch=np.linalg.norm(P-x,axis=1); i=np.argmin(dch)
        if dch[i]<0.05*np.min([np.linalg.norm(P[j]-P[k]) for j in range(4) for k in range(4) if j!=k]): return i,it
        dc=np.min(np.linalg.norm(crit-x,axis=1))
        h=frac*min(max(dc,1e-12),dch[i])
        k1=g_unit(x);k2=g_unit(x+h/2*k1);k3=g_unit(x+h/2*k2);k4=g_unit(x+h*k3)
        x=x+h/6*(k1+2*k2+2*k3+k4)
    return None,it
edges=[]; detail=[]
for x0,idx in Z:
    if idx!=2: continue
    _,H=finder.gH(x0,P,Q); w,V=np.linalg.eigh(H); v=V[:,2]
    a,ia=climb(x0,v,1); b,ib=climb(x0,v,-1)
    edges.append(tuple(sorted((a,b)))); detail.append((np.linalg.norm(x0),a,b,ia,ib))
for d in sorted(detail): print('index-2 point at r=%.3e: walls to charges %s and %s (steps %d, %d)'%d)
cnt=collections.Counter(edges); print('walls:',dict(cnt))
