"""One level of a hierarchical unfolding: real nonzero critical points of A + B, A harmonic homogeneous of
degree p, B of degree p-1 (scale-free). Complex count is a^2+ab+b^2 with a=p-1, b=p-2 (rest sit at 0).
Real ones <-> critical points on RP^2 of (p-1)log|A| - p log|B|. Numerical."""
import numpy as np
from harm import monos, harmonic_basis, Poly, solve_all, evalmon

def emb(c_low, dlow, d):
    """embed coefficient vector of degree<=dlow into degree<=d monomial list"""
    L1,_=monos(dlow); L2,ix=monos(d); out=np.zeros(len(L2))
    for k,m in enumerate(L1): out[ix[m]]=c_low[k]
    return out

def basis_h(p):  # homogeneous harmonics degree p, in MON[p]
    return harmonic_basis(p,p)

def sectoral(p):
    L,ix=monos(p); c=np.zeros(len(L))
    from math import comb
    # Re (x+iy)^p = sum_k C(p,k) x^{p-k} (iy)^k real part
    for k in range(0,p+1,2):
        c[ix[(p-k,k,0)]]+=comb(p,k)*(-1)**(k//2)
    return c

def zonal(p):
    # r^p P_p(z/r) = sum_j coeff z^{p-2j} r^{2j}
    import sympy as sp
    x,y,z=sp.symbols('x y z'); r2=x*x+y*y+z*z
    t=sp.symbols('t'); P=sp.legendre(p,t)
    expr=sp.expand(sum(P.coeff(t,p-2*j)*z**(p-2*j)*r2**j for j in range(p//2+1)))
    L,ix=monos(p); c=np.zeros(len(L))
    for (i,j,k),v in sp.Poly(expr,x,y,z).terms(): c[ix[(i,j,k)]]=float(v)
    return c

def count_level(cA,cB,p,rng,nseed=4000,ret=False):
    c=cA+emb(cB,p-1,p)
    P=Poly(c,p)
    X=rng.standard_normal((nseed,3)); X*= (10**rng.uniform(-1.5,1.5,(nseed,1)))*np.abs(cB).max()/np.abs(cA).max()/np.linalg.norm(X,axis=1)[:,None]
    for _ in range(80):
        g,H=P.gradhess(X)
        det=np.linalg.det(H); ok=np.abs(det)>0
        st=np.zeros_like(X); st[ok]=np.linalg.solve(H[ok],g[ok][...,None])[...,0]
        nx=np.linalg.norm(X,axis=1); ns=np.linalg.norm(st,axis=1)
        f=np.minimum(1,0.3*nx/np.maximum(ns,1e-300)); X=X-st*f[:,None]
    g,H=P.gradhess(X); nx=np.linalg.norm(X,axis=1)
    det=np.linalg.det(H); ok=np.isfinite(nx)&(np.abs(det)>0)
    st=np.full(len(X),np.inf); st[ok]=np.linalg.norm(np.linalg.solve(H[ok],g[ok][...,None])[...,0],axis=1)
    ok&=(st<1e-10*nx)&(nx>1e-8)
    Y=X[ok]
    keep=[]
    for i in np.argsort(-np.linalg.norm(Y,axis=1)):
        if all(np.linalg.norm(Y[i]-Y[j])>1e-6*np.linalg.norm(Y[i]) for j in keep): keep.append(i)
    Y=Y[keep]
    g,H=P.gradhess(Y); ev=np.linalg.eigvalsh(H)
    if ret: return Y,(ev<0).sum(1),ev
    return len(Y),(ev<0).sum(1)

if __name__=="__main__":
    rng=np.random.default_rng(0)
    for p in range(3,9):
        a,b=p-1,p-2
        for nm,cA,cB in [("sect/zon",sectoral(p),zonal(p-1)),("zon/sect",zonal(p),sectoral(p-1))]:
            # small generic perturbation to break symmetry
            cA2=cA+1e-3*np.abs(cA).max()*(basis_h(p)@rng.standard_normal(2*p+1))
            cB2=cB+1e-3*np.abs(cB).max()*(basis_h(p-1)@rng.standard_normal(2*p-1))
            k,idx=count_level(cA2,cB2,p,rng,nseed=6000)
            forced = p*(p-2) if nm=="sect/zon" else (p-1)**2
            print(f"p={p} {nm}: real nonzero {k} (index1 {(idx==1).sum()}, index2 {(idx==2).sum()}); complex {a*a+a*b+b*b}; grid-forced {forced}",flush=True)
