"""Four non-coplanar charges via the inverse map.

For sites A1..A4 (a tetrahedron) and x inside it, there is exactly one charge ray making x
an equilibrium: zeta_i = r_i^2 w_i with w the null vector of [u1 u2 u3 u4] (u_i unit vectors
from A_i to x). The Hessian H(x) of V is traceless, so det H < 0 <=> index 1, det H > 0 <=> index 2.
Total equilibria for a charge vector = 3 + 2 m1, m1 = #preimages in R1 = {det H < 0}.
Maxwell's bound 9 for four charges <=> m1 <= 3 for every charge vector.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np

def zeta_D(A, X):
    d = X[:, None, :] - A[None]
    r2 = (d * d).sum(2); r = np.sqrt(r2)
    u = d / r[..., None]                                   # (m,4,3)
    w = np.empty((X.shape[0], 4))
    for i in range(4):
        cols = [j for j in range(4) if j != i]
        w[:, i] = (-1) ** i * np.linalg.det(u[:, cols, :])
    w *= np.sign(w.sum(1, keepdims=True))
    z = w * r2
    z /= z.sum(1, keepdims=True)
    i3 = z / (r2 * r); i5 = z / (r2 * r2 * r)
    H = 3 * np.einsum('mn,mnk,mnl->mkl', i5, d, d) - i3.sum(1)[:, None, None] * np.eye(3)
    S = (z / r2).sum(1)
    return z, np.linalg.det(H) / S ** 3, (w > 0).all(1)

# Kuhn triangulation of the unit cube into 6 tetrahedra
_perm = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]

def coverage(A, N=64, P=64, box=None, return_cov=False):
    if box is None:
        box = (0, 1, 0, 1, 0, 1)
    sp = lambda lo, hi: lo + (hi - lo) * (0.5 - 0.5 * np.cos(np.pi * np.linspace(0, 1, N + 1)))
    s, t, v = sp(*box[0:2]), sp(*box[2:4]), sp(*box[4:6])
    S, T, Vv = np.meshgrid(s, t, v, indexing='ij')
    l = np.stack([1 - S, S * (1 - T), S * T * (1 - Vv), S * T * Vv], -1)
    X = (l @ A).reshape(-1, 3)
    inside = (l > 1e-14).all(-1).ravel()
    Z = np.full((X.shape[0], 3), np.nan); D = np.full(X.shape[0], np.nan)
    z, Dd, pos = zeta_D(A, X[inside])
    Z[inside] = z[:, :3]; D[inside] = Dd
    Z = Z.reshape(N + 1, N + 1, N + 1, 3); D = D.reshape(N + 1, N + 1, N + 1)
    neg = D < 0
    if not neg.any():
        return 0, None, None
    ii = np.argwhere(neg)
    lo_i = np.maximum(ii.min(0) - 2, 0); hi_i = np.minimum(ii.max(0) + 2, N)
    sbox = (s[lo_i[0]], s[hi_i[0]], t[lo_i[1]], t[hi_i[1]], v[lo_i[2]], v[hi_i[2]])
    zz = Z[neg]
    zlo = zz.min(0); w = np.maximum(zz.max(0) - zlo, 1e-300) * 1.02; zlo = zlo - 0.01 * w
    cov = np.zeros((P, P, P), np.int16)
    I = np.stack(np.meshgrid(np.arange(N), np.arange(N), np.arange(N), indexing='ij'), -1).reshape(-1, 3)
    cellneg = np.ones(len(I), bool)
    for o in np.ndindex(2, 2, 2):
        cellneg &= neg[I[:, 0] + o[0], I[:, 1] + o[1], I[:, 2] + o[2]]
    I = I[cellneg]
    for p in _perm:
        corners = [np.zeros(3, int)]
        c = np.zeros(3, int)
        for k in p:
            c = c.copy(); c[k] = 1; corners.append(c)
        Ps = [(Z[I[:, 0] + c[0], I[:, 1] + c[1], I[:, 2] + c[2]] - zlo) / w * P for c in corners]
        P0 = Ps[0]; M = np.stack([Ps[1] - P0, Ps[2] - P0, Ps[3] - P0], -1)   # columns
        det = np.linalg.det(M)
        good = np.isfinite(det) & (np.abs(det) > 1e-30)
        P0, M, Pall = P0[good], M[good], [q[good] for q in Ps]
        Minv = np.linalg.inv(M)
        lo = np.floor(np.minimum.reduce(Pall) - 0.5).astype(int)
        hi = np.ceil(np.maximum.reduce(Pall) - 0.5).astype(int)
        span = int((hi - lo).max()) if len(lo) else -1
        for o in np.ndindex(span + 1, span + 1, span + 1):
            px = lo + np.array(o)
            vv = px + 0.5 - P0
            b = np.einsum('mij,mj->mi', Minv, vv)
            ins = (b >= 0).all(1) & (b.sum(1) < 1) & (px >= 0).all(1) & (px < P).all(1)
            np.add.at(cov, (px[ins, 0], px[ins, 1], px[ins, 2]), 1)
    if return_cov:
        return int(cov.max()), sbox, (cov, zlo, w)
    return int(cov.max()), sbox, None

def max_m1(A, N1=40, N2=72, P=72):
    m, sbox, _ = coverage(A, N=N1, P=40)
    if sbox is None:
        return 0, 0
    m2, _, _ = coverage(A, N=N2, P=P, box=sbox)
    return m, m2

if __name__ == '__main__':
    import time
    s3 = np.sqrt(3) / 2
    for h in (0.05, 0.2, 0.5, 0.816):
        A = np.array([[1, 0, 0], [-0.5, s3, 0], [-0.5, -s3, 0], [0, 0, h]])
        t = time.time(); print('triangle + apex h=%g' % h, max_m1(A), '%.1fs' % (time.time() - t))
