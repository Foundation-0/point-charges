import numpy as np, time
from level2d import *
from harm import monos
from numpy.polynomial import legendre
def rotate_zx(c,p):
    """zonal about z -> zonal about x: substitute (x,y,z)->(z,y,x)... (swap x and z)"""
    L,ix=monos(p); out=np.zeros(len(L))
    for k,(i,j,l) in enumerate(L): out[ix[(l,j,i)]]=c[k]
    return out
def crossings(p):
    # pairs (a,b) of zeros of P_p, P_{p-1} with a^2+b^2<1: each gives 2 points on S^2 -> 1 on RP^2? (x,y,z) and (x,-y,z): distinct on RP^2 unless...
    za=legendre.Legendre.basis(p).roots(); zb=legendre.Legendre.basis(p-1).roots()
    n=sum(1 for a in za for b in zb if a*a+b*b<1)
    return n  # points on S^2 = 2n; on RP^2 = n
