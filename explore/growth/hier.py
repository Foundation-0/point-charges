"""Hierarchical polynomial P = sum_{k=1}^{d} c_k A_k, A_k zonal about z (k-d even) or about x (k-d odd), with
c_{k-1} = s * c_k * (scale). Prediction: #real critical points = sum_{k=2}^{d} (I_k - 1). Numerical check."""
import numpy as np, sys, time
from harm import Poly, solve_all, monos
from level import zonal, emb
from t_cross_util import rotate_zx, crossings
def hier_poly(d, s, axis_top='z'):
    L,ix=monos(d); c=np.zeros(len(L)); coef=1.0
    for k in range(d,0,-1):
        z=zonal(k); z=z/np.abs(z).max()
        ax = axis_top if (d-k)%2==0 else ('x' if axis_top=='z' else 'z')
        if ax=='x': z=rotate_zx(z,k)
        c+=coef*emb(z,k,d); coef*=s**(d-k+1)   # level k radius rho_k = c_{k-1}/c_k = s^(d-k+1)
    return c
if __name__=="__main__":
    s=float(sys.argv[1]); rng=np.random.default_rng(0)
    for d in range(3,int(sys.argv[2])+1):
        c=hier_poly(d,s); P=Poly(c,d)
        t=time.time()
        X=np.vstack([rng.standard_normal((3000,3))*(s**j)*10**rng.uniform(-1.3,1.3,(3000,1)) for j in range(1,d)])
        from harm import newton_real
        for _ in range(3):
            X=newton_real(P,X,iters=60)
        g,H=P.gradhess(X); nx=np.linalg.norm(X,axis=1)
        st=np.linalg.norm(np.linalg.solve(H,g[...,None])[...,0],axis=1); ok=np.isfinite(nx)&(st<1e-9*nx)&(nx>0)
        X=X[ok]; keep=[]
        for i in np.argsort(nx[ok]):
            if all(np.linalg.norm(X[i]-X[j])>1e-6*np.linalg.norm(X[i]) for j in keep): keep.append(i)
        R=X[keep]; g,H=P.gradhess(R); ev=np.linalg.eigvalsh(H); sols=R; lost=-1; r=np.ones(len(R),bool)
        pred=sum(crossings(k)-1 for k in range(2,d+1))
        print(f"d={d} s={s}: complex {len(sols)} (Bezout {(d-1)**3}, lost {lost}), real {r.sum()} (index1 {(ev<0).sum(1).tolist().count(1)}, index2 {(ev<0).sum(1).tolist().count(2)}), predicted {pred}, min|eig|/max {np.min(np.abs(ev).min(1)/np.abs(ev).max(1)):.1e}, |x| range {np.linalg.norm(R,axis=1).min():.1e}..{np.linalg.norm(R,axis=1).max():.1e}, {time.time()-t:.0f}s",flush=True)
