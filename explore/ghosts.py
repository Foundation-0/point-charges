"""Near-miss detection: local minima of the normalised field |E|/sum(q/r^2).

Zeros are equilibria; positive local minima ("ghosts") mark where a small change of
the configuration can create a new pair of equilibria (a fold bifurcation).
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np
from charges import grad_hess, dedup


def resid_jac(x, a, q):
    d = x[:, None, :] - a[None, :, :]
    r2 = np.einsum('mnk,mnk->mn', d, d)
    r = np.sqrt(r2)
    g, H = grad_hess(x, a, q)
    S = (q[None, :] / r2).sum(1)
    dS = (-2 * q[None, :, None] * d / (r2 * r2)[:, :, None]).sum(1)
    R = g / S[:, None]
    J = H / S[:, None, None] - g[:, :, None] * dS[:, None, :] / (S * S)[:, None, None]
    return R, J


def lm(x, a, q, iters=80, scale=1.0):
    lam = np.full(len(x), 1e-3)
    R, J = resid_jac(x, a, q)
    f = (R * R).sum(1)
    for _ in range(iters):
        JT = np.transpose(J, (0, 2, 1))
        A = JT @ J + lam[:, None, None] * np.eye(3)[None]
        b = (JT @ R[..., None])[..., 0]
        step = np.linalg.solve(A, b[..., None])[..., 0]
        n = np.linalg.norm(step, axis=1)
        step *= np.minimum(1, 0.05 * scale / np.maximum(n, 1e-300))[:, None]
        xn = x - step
        Rn, Jn = resid_jac(xn, a, q)
        fn = (Rn * Rn).sum(1)
        better = fn < f
        x = np.where(better[:, None], xn, x)
        R = np.where(better[:, None], Rn, R)
        J = np.where(better[:, None, None], Jn, J)
        f = np.where(better, fn, f)
        lam = np.where(better, lam * 0.3, lam * 10)
        lam = np.clip(lam, 1e-12, 1e8)
    return x, np.sqrt(f)


def minima(a, q, nseed=3000, seed=0):
    """Return (points, |R| values) of local minima of |R| inside the bounding box."""
    rng = np.random.default_rng(seed)
    lo, hi = a.min(0), a.max(0)
    span = np.maximum(hi - lo, 1e-3 * np.max(hi - lo))
    sc = np.max(span)
    s1 = lo + span * rng.random((nseed, 3))
    s2 = rng.dirichlet(np.ones(len(a)) * 0.7, size=nseed) @ a
    x, v = lm(np.vstack([s1, s2]), a, q, scale=sc)
    d = x[:, None, :] - a[None, :, :]
    rmin = np.sqrt(np.einsum('mnk,mnk->mn', d, d)).min(1)
    ok = np.isfinite(v) & (rmin > 1e-4 * sc) & (v < 0.5)
    x, v = x[ok], v[ok]
    # dedup keeping the smaller value
    o = np.argsort(v); x, v = x[o], v[o]
    keep = []
    for i, p in enumerate(x):
        if all(np.linalg.norm(p - x[j]) > 1e-5 * sc for j in keep):
            keep.append(i)
    return x[keep], v[keep]
