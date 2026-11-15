"""Robust count of real nonzero critical points of A+B (A degree p, B degree p-1, harmonic homogeneous):
they are x = t w with w on the sphere where grad A(w) || grad B(w), t = -(gA.gB)/|gA|^2.
Search: dense Fibonacci grid on the sphere, local minima of the parallelism residual, then 3D Newton. Numerical."""
import numpy as np
from scipy.spatial import cKDTree
from harm import Poly
from level import emb, basis_h, zonal, sectoral
import finder

def crit_level(cA, cB, p, ngrid=400000, k=8):
    P = Poly(cA + emb(cB, p - 1, p), p)
    PA = Poly(cA, p); PB = Poly(emb(cB, p - 1, p), p)
    W = finder.fib_sphere(ngrid); W = W[W[:, 2] >= 0]  # RP^2: upper hemisphere (+equator)
    gA = PA.grad(W); gB = PB.grad(W)
    res = np.linalg.norm(np.cross(gA, gB), axis=1) / (np.linalg.norm(gA, axis=1) * np.linalg.norm(gB, axis=1))
    tree = cKDTree(np.vstack([W, -W]))
    _, nb = tree.query(W, k=k + 1)
    resall = np.r_[res, res]
    ismin = np.all(res[:, None] <= resall[nb[:, 1:]], axis=1) & (res < 0.2)
    Ws = W[ismin]
    gA = PA.grad(Ws); gB = PB.grad(Ws)
    t = -(gA * gB).sum(1) / (gA * gA).sum(1)
    X = t[:, None] * Ws
    for _ in range(60):
        g, H = P.gradhess(X)
        st = np.linalg.solve(H, g[..., None])[..., 0]
        nx = np.linalg.norm(X, axis=1)
        f = np.minimum(1, 0.2 * nx / np.maximum(np.linalg.norm(st, axis=1), 1e-300))
        X = X - st * f[:, None]
    g, H = P.gradhess(X); nx = np.linalg.norm(X, axis=1)
    st = np.linalg.norm(np.linalg.solve(H, g[..., None])[..., 0], axis=1)
    ok = np.isfinite(nx) & (st < 1e-10 * nx) & (nx > 0)
    X = X[ok]
    return X

def dedup_pts(X, rel=1e-7):
    keep = []
    for i in range(len(X)):
        if all(np.linalg.norm(X[i] - X[j]) > rel * np.linalg.norm(X[i]) for j in keep):
            keep.append(i)
    return X[keep]
