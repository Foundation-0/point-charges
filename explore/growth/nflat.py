"""Smallest n (found) with positive charges flat to order d-1 at a point (jets 1..d-1 vanish). Numerical search."""
import numpy as np, sys, time, json
from flatsolve import solve, eqs, jac
d=int(sys.argv[1]); ns=[int(v) for v in sys.argv[2].split(',')]; tries=int(sys.argv[3])
for n in ns:
    rng=np.random.default_rng(n*100+d); t=time.time(); hit=None
    for k in range(tries):
        r=solve(n,d,rng,iters=300,smin=-8)
        if r is not None:
            a,q=r; f,J=jac(a,np.log(q),d); rk=np.linalg.matrix_rank(J,tol=1e-8*np.abs(J).max())
            hit=(k,rk); json.dump({"d":d,"n":n,"sites":a.tolist(),"charges":q.tolist(),"jac_rank":int(rk),"conditions":d*d-1},open(f"configs/flat_d{d}_n{n}.json","w"))
            break
    print(f"d={d} n={n}: {'FOUND at try %d, jet Jacobian rank %d of %d conditions'%(hit[0],hit[1],d*d-1) if hit else 'none in %d tries'%tries} ({time.time()-t:.0f}s)",flush=True)
    if hit: break
