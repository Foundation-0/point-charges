"""Three charges via the inverse map.

For a triangle A1 A2 A3 and a point x inside it there is exactly one charge vector (up to
scale) for which x is an equilibrium (Lami's theorem):
    zeta_i = r_i^2 * w_i,  w_1 : w_2 : w_3 = (u2 x u3) : (u3 x u1) : (u1 x u2),  u_i = (x - A_i)/r_i.
Let D(x) = det of the in-plane Hessian of V at x for that charge vector. x is a local
minimum (in the plane) iff D(x) > 0. Maxwell's bound 4 for three charges with this
triangle  <=>  Z: x -> [zeta(x)] is one-to-one on R_min = {D > 0}.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np

def cross(u, v):
    return u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]

def zeta_D(A, X):
    d = X[:, None, :] - A[None]
    r2 = (d * d).sum(2); r = np.sqrt(r2)
    u = d / r[..., None]
    w = np.stack([cross(u[:, 1], u[:, 2]), cross(u[:, 2], u[:, 0]), cross(u[:, 0], u[:, 1])], 1)
    w = w * np.sign(w.sum(1, keepdims=True))
    z = w * r2
    z = z / z.sum(1, keepdims=True)
    i3 = z / (r2 * r); i5 = z / (r2 * r2 * r)
    H = 3 * np.einsum('mn,mnk,mnl->mkl', i5, d, d) - i3.sum(1)[:, None, None] * np.eye(2)
    D = H[:, 0, 0] * H[:, 1, 1] - H[:, 0, 1] ** 2
    # normalise D by the field scale squared so it is comparable across x
    S = (z / r2).sum(1)
    return z, D / (S * S)

def grid_bary(N, edge_refine=True):
    """Points of the triangle in barycentric coordinates; refined near the edges."""
    t = np.linspace(0, 1, N + 1)
    if edge_refine:
        t = 0.5 - 0.5 * np.cos(np.pi * t)     # Chebyshev-like spacing
    L1, L2 = np.meshgrid(t, t, indexing='ij')
    return L1, L2

def analyse(A, N=400, P=600):
    """Max number of in-plane minima over all charge vectors, estimated by rasterising the
    image of R_min cells into a P x P grid on the charge simplex."""
    # square grid on (s, t) in [0,1]^2 mapped to triangle by x = A1 + s (A2-A1) + t(1-s)(A3-A1)... use Duffy map
    s = 0.5 - 0.5 * np.cos(np.pi * np.linspace(0, 1, N + 1))
    S, T = np.meshgrid(s, s, indexing='ij')
    l2 = S * (1 - T); l3 = S * T; l1 = 1 - S
    X = l1[..., None] * A[0] + l2[..., None] * A[1] + l3[..., None] * A[2]
    inside = (l1 > 0) & (l2 > 0) & (l3 > 0)
    Xf = X.reshape(-1, 2)
    z = np.full((Xf.shape[0], 3), np.nan); D = np.full(Xf.shape[0], np.nan)
    m = inside.ravel()
    z[m], D[m] = zeta_D(A, Xf[m])
    Z = z[:, :2].reshape(N + 1, N + 1, 2); D = D.reshape(N + 1, N + 1)
    cov = np.zeros((P, P), np.int32)
    # each grid square -> two triangles; keep triangles whose three vertices have D > 0
    tris = []
    for (i0, j0), (i1, j1), (i2, j2) in (((0, 0), (1, 0), (1, 1)), ((0, 0), (1, 1), (0, 1))):
        a_ = Z[i0:N + i0 - 1 + 1 - 0 if False else slice(i0, N + i0), slice(j0, N + j0)] if False else None
    idx = np.arange(N)
    I, J = np.meshgrid(idx, idx, indexing='ij')
    out = 0; maxcov = 0
    for (di0, dj0), (di1, dj1), (di2, dj2) in (((0, 0), (1, 0), (1, 1)), ((0, 0), (1, 1), (0, 1))):
        P0 = Z[I + di0, J + dj0]; P1 = Z[I + di1, J + dj1]; P2 = Z[I + di2, J + dj2]
        Dm = np.minimum(np.minimum(D[I + di0, J + dj0], D[I + di1, J + dj1]), D[I + di2, J + dj2])
        good = (Dm > 0) & np.isfinite(P0).all(-1) & np.isfinite(P1).all(-1) & np.isfinite(P2).all(-1)
        P0, P1, P2 = P0[good], P1[good], P2[good]
        # rasterise: pixel centres (k+0.5)/P in the unit square of (zeta1, zeta2)
        lo = np.floor(np.minimum(np.minimum(P0, P1), P2) * P - 0.5).astype(int)
        hi = np.ceil(np.maximum(np.maximum(P0, P1), P2) * P - 0.5).astype(int)
        span = (hi - lo).max() if len(lo) else 0
        for ox in range(span + 1):
            for oy in range(span + 1):
                px = lo[:, 0] + ox; py = lo[:, 1] + oy
                c = np.c_[(px + 0.5) / P, (py + 0.5) / P]
                # barycentric test with half-open edges
                v0 = P1 - P0; v1 = P2 - P0; v2 = c - P0
                den = v0[:, 0] * v1[:, 1] - v0[:, 1] * v1[:, 0]
                with np.errstate(divide='ignore', invalid='ignore'):
                    b1 = (v2[:, 0] * v1[:, 1] - v2[:, 1] * v1[:, 0]) / den
                    b2 = (v0[:, 0] * v2[:, 1] - v0[:, 1] * v2[:, 0]) / den
                ins = (b1 >= 0) & (b2 >= 0) & (b1 + b2 < 1) & (px >= 0) & (py >= 0) & (px < P) & (py < P) & (px <= hi[:, 0]) & (py <= hi[:, 1])
                np.add.at(cov, (px[ins], py[ins]), 1)
    return cov, D, Z

if __name__ == '__main__':
    s3 = np.sqrt(3) / 2
    for A in (np.array([[0, 0], [1, 0], [0.5, s3]]), np.array([[0, 0], [1, 0], [0.3, 0.4]]), np.array([[0, 0], [1, 0], [0.5, 0.15]])):
        cov, D, Z = analyse(A)
        print(A[2], 'max coverage', cov.max(), 'pixels covered', (cov > 0).sum(), 'pixels with >=2', (cov >= 2).sum(), 'D>0 frac', np.nanmean(D > 0))
