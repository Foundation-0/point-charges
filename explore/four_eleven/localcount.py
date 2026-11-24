# Count balance points near the origin by multistart Newton at many scales (candidates, float64),
# then polish each in mpmath (40 digits) and report Hessian eigenvalues.
import numpy as np, mpmath as mp
import finder
def fib_dirs(n):
    i=np.arange(n)+0.5; phi=np.arccos(1-2*i/n); th=np.pi*(1+5**0.5)*i
    return np.stack([np.cos(th)*np.sin(phi),np.sin(th)*np.sin(phi),np.cos(phi)],1)
def newton_local(x,P,Q,Rmax,it=80):
    for k in range(it):
        g,H=finder.gH(x,P,Q)
        try: s=np.linalg.solve(H,g)
        except np.linalg.LinAlgError: return None
        ns=np.linalg.norm(s); cap=0.5*max(np.linalg.norm(x),1e-7)
        if ns>cap: s*=cap/ns
        x=x-s
        if np.linalg.norm(x)>Rmax: return None
        if ns<1e-15*max(np.linalg.norm(x),1e-12): break
    return x
def polish(x,P,Q,dps=40):
    mp.mp.dps=dps
    X=mp.matrix([mp.mpf(float(c)) for c in x])
    Pm=[[mp.mpf(float(c)) for c in p] for p in P]; Qm=[mp.mpf(float(q)) for q in Q]
    for k in range(60):
        g=mp.matrix(3,1); H=mp.matrix(3,3)
        for a,q in zip(Pm,Qm):
            d=[X[i]-a[i] for i in range(3)]; r2=sum(c*c for c in d); r=mp.sqrt(r2)
            for i in range(3):
                g[i]-=q*d[i]/r**3
                for j in range(3): H[i,j]+=q*(3*d[i]*d[j]/r**5-(1/r**3 if i==j else 0))
        s=mp.lu_solve(H,g); X=X-s
        if mp.norm(s)<mp.mpf(10)**(-dps+8)*(1+mp.norm(X)): break
    ev=sorted([float(e) for e in mp.eig(H)[0]],key=lambda e:e.real if hasattr(e,'real') else e)
    return np.array([float(c) for c in X]),float(mp.norm(g)),ev
def count_local(P,Q,Rmax=0.2,rmin=1e-7,nr=40,nd=200,extra=(),dedup=1e-3):
    P=np.asarray(P,float);Q=np.asarray(Q,float)
    D=fib_dirs(nd); sols=[]
    starts=[r*d for r in np.geomspace(rmin,Rmax,nr) for d in D]+list(extra)
    for s in starts:
        x=newton_local(np.array(s,float),P,Q,Rmax*1.5)
        if x is None: continue
        g,H=finder.gH(x,P,Q)
        if np.linalg.norm(g)>1e-12: continue
        if any(np.linalg.norm(x-u)<dedup*max(np.linalg.norm(u),1e-9) for u in sols): continue
        sols.append(x)
    return sols
