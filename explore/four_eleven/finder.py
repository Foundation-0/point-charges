# Robust numerical finder for zeros of grad V, V = sum q_i/|x-a_i| (candidates only, float64).
import numpy as np
def gH(x,P,Q):
    d=x[None,:]-P; r2=np.sum(d*d,1); r=np.sqrt(r2); r3=r2*r; r5=r3*r2
    g=-np.sum((Q/r3)[:,None]*d,0)
    H=np.sum(Q[:,None,None]*(3*d[:,:,None]*d[:,None,:]/r5[:,None,None]-np.eye(3)[None]/r3[:,None,None]),0)
    return g,H
def newton(x,P,Q,it=200,tol=1e-13,R=None):
    c0=P.mean(0); R=R or 3*np.max(np.linalg.norm(P-c0,axis=1))
    for k in range(it):
        if np.linalg.norm(x-c0)>R: return None
        g,H=gH(x,P,Q)
        try: s=np.linalg.solve(H,g)
        except np.linalg.LinAlgError: return None
        # damping: keep away from charges
        dmin=np.min(np.linalg.norm(x[None]-P,axis=1))
        ns=np.linalg.norm(s)
        if ns>0.5*dmin: s*=0.5*dmin/ns
        x=x-s
        if ns<tol*(1+np.linalg.norm(x)): return x
    return x if np.linalg.norm(gH(x,P,Q)[0])<1e-9 else None
def find_all(P,Q,nstart=4000,seed=0,extra=None,dedup=1e-7):
    P=np.asarray(P,float); Q=np.asarray(Q,float)
    rng=np.random.default_rng(seed)
    lo=P.min(0); hi=P.max(0); c=(lo+hi)/2; w=(hi-lo).max()
    starts=c+w*(rng.random((nstart,3))-0.5)*1.2
    # also starts near each pair segment and each triangle
    n=len(P); S=[]
    for i in range(n):
        for j in range(i+1,n):
            for tt in np.linspace(0.05,0.95,19): S.append(P[i]+tt*(P[j]-P[i]))
    starts=np.vstack([starts,np.array(S)+1e-3*rng.standard_normal((len(S),3))])
    if extra is not None: starts=np.vstack([starts,extra])
    sols=[]
    for s in starts:
        r=newton(s.copy(),P,Q)
        if r is None or not np.all(np.isfinite(r)): continue
        g,H=gH(r,P,Q)
        if np.linalg.norm(g)>1e-9*np.max(Q)*np.max(1/np.linalg.norm(r[None]-P,axis=1))**2: continue
        if any(np.linalg.norm(r-u)<dedup for u in sols): continue
        sols.append(r)
    out=[]
    for r in sols:
        g,H=gH(r,P,Q); ev=np.linalg.eigvalsh(H)
        out.append((r,ev))
    return out
