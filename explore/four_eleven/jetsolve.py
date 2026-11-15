# Solve for a four-charge configuration with a prescribed 3-jet of V at the origin.
# Unknowns: 12 site coordinates + 4 charges. Equations (16):
#   grad V(0) = lam (3); Hessian entries H_xx,H_yy,H_xy,H_xz,H_yz = A (5) (trace is 0 automatically);
#   cubic T1 parts tau = 0 (3, fixes rotation); cubic T2 parts nu (3); xyz coefficient = cxyz (1); q0 = 1 (1).
# Cubic C(x) = sum c_abc x^a y^b z^c = (1/6) T3 x x x.  Harmonic cubic basis:
#   xyz;  T1: x(y^2-z^2), y(z^2-x^2), z(x^2-y^2);  T2: x(2x^2-3y^2-3z^2) etc.
#   nu_x = c_xxx/2, tau_x = (c_xyy - c_xzz)/2 (cyclic).
import mpmath as mp
mp.mp.dps=40
def jets(P,Q):
    g=[mp.mpf(0)]*3; H=[[mp.mpf(0)]*3 for _ in range(3)]; T=[[[mp.mpf(0)]*3 for _ in range(3)] for _ in range(3)]
    for a,q in zip(P,Q):
        r2=sum(c*c for c in a); r=mp.sqrt(r2); r3=r*r2; r5=r3*r2; r7=r5*r2
        for k in range(3):
            g[k]+=q*a[k]/r3
            for l in range(3):
                H[k][l]+=q*(3*a[k]*a[l]-(r2 if k==l else 0))/r5
                for m in range(3):
                    d=lambda i,j: 1 if i==j else 0
                    T[k][l][m]+=q*(15*a[k]*a[l]*a[m]/r7-3*(d(k,l)*a[m]+d(k,m)*a[l]+d(l,m)*a[k])/r5)
    return g,H,T
def cubic_coeffs(T):
    # c_abc for monomial x^a y^b z^c: sum over index triples /6
    import itertools
    c={}
    for k,l,m in itertools.product(range(3),repeat=3):
        e=[0,0,0]; e[k]+=1;e[l]+=1;e[m]+=1; e=tuple(e)
        c[e]=c.get(e,0)+T[k][l][m]/6
    return c
def cubic_parts(T):
    c=cubic_coeffs(T)
    xyz=c[(1,1,1)]
    nu=[c[(3,0,0)]/2,c[(0,3,0)]/2,c[(0,0,3)]/2]
    tau=[(c[(1,2,0)]-c[(1,0,2)])/2,(c[(0,1,2)]-c[(2,1,0)])/2,(c[(2,0,1)]-c[(0,2,1)])/2]
    return xyz,nu,tau
def unpack(v): return [v[0:3],v[3:6],v[6:9],v[9:12]],list(v[12:16])
def F(v,tgt):
    P,Q=unpack(v); g,H,T=jets(P,Q); xyz,nu,tau=cubic_parts(T)
    lam,A,nut,cx=tgt
    out=[g[i]-lam[i] for i in range(3)]
    out+= [H[0][0]-A[0],H[1][1]-A[1],H[0][1]-A[2],H[0][2]-A[3],H[1][2]-A[4]]
    out+= tau
    out+= [nu[i]-nut[i] for i in range(3)]
    out+= [xyz-cx, Q[0]-1]
    return out
def solve(tgt,v0,it=40,tol=mp.mpf(10)**-32,verbose=False):
    v=[mp.mpf(x) for x in v0]
    for k in range(it):
        f=F(v,tgt); nf=mp.norm(mp.matrix(f))
        if verbose: print(k,mp.nstr(nf,5))
        if nf<tol: return v
        h=mp.mpf(10)**-20; J=mp.matrix(16,16)
        for j in range(16):
            w=list(v); w[j]+=h; fj=F(w,tgt)
            for i in range(16): J[i,j]=(fj[i]-f[i])/h
        dv=mp.lu_solve(J,mp.matrix(f))
        v=[v[i]-dv[i] for i in range(16)]
    raise RuntimeError('no convergence %s'%mp.nstr(nf,5))
TET=[1,1,1, 1,-1,-1, -1,1,-1, -1,-1,1, 1,1,1,1]
if __name__=='__main__':
    P,Q=unpack([mp.mpf(x) for x in TET]); g,H,T=jets(P,Q); print(cubic_parts(T))
    cx=cubic_parts(T)[0]
    v=solve(([0,0,0],[0]*5,[mp.mpf('0.005'),mp.mpf('0.006'),mp.mpf('0.004')],cx),TET,verbose=True)
    print([mp.nstr(x,12) for x in v])
