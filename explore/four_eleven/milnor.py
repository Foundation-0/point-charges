# Local multiplicity (Milnor number) of grad f at 0 via dim Q[x]/(J + m^k), exact over Q (flint fmpq_mat).
# The sequence stabilises at mu once m^k is inside the local primary component.
import sys, itertools
import sympy as sp
from flint import fmpq_mat, fmpq
x,y,z=sp.symbols('x y z')
X=(x,y,z)
def monos(deg):
    return [(a,b,deg-a-b) for a in range(deg+1) for b in range(deg+1-a)]
def local_dim(f,k):
    g=[sp.Poly(sp.diff(f,v),*X) for v in X]
    allm=[m for d in range(k) for m in monos(d)]
    idx={m:i for i,m in enumerate(allm)}
    rows=[]
    for gi in g:
        lowdeg=min(sum(m) for m in gi.as_dict())
        for d in range(0,k-lowdeg):
            for m in monos(d):
                row=[fmpq(0)]*len(allm)
                nz=False
                for mm,c in gi.as_dict().items():
                    e=(mm[0]+m[0],mm[1]+m[1],mm[2]+m[2])
                    if sum(e)<k:
                        c=sp.Rational(c); row[idx[e]]=fmpq(int(c.p),int(c.q)); nz=True
                if nz: rows.append(row)
    M=fmpq_mat(len(rows),len(allm),[v for r in rows for v in r])
    return len(allm)-M.rank()
K4=x**4+y**4+z**4-3*(x**2*y**2+y**2*z**2+z**2*x**2)
K6=2*x**6-15*x**4*y**2-15*x**4*z**2-15*x**2*y**4+180*x**2*y**2*z**2-15*x**2*z**4+2*y**6-15*y**4*z**2-15*y**2*z**4+2*z**6
# V/sqrt3 - const, exact: 20/27 xyz - 14/243 K4 + 4/729 K6
f4=sp.Rational(20,27)*x*y*z-sp.Rational(14,243)*K4
f6=f4+sp.Rational(4,729)*K6
for name,f in [('tet deg<=4',f4),('tet deg<=6',f6),('xyz+x4+y4+z4',x*y*z+x**4+y**4+z**4),('xyz only',x*y*z)]:
    print(name,[local_dim(f,k) for k in range(3,10)],flush=True)
