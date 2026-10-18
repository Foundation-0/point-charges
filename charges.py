"""Numerical tools for the equilibria (critical points) of point-charge potentials.

V(x) = sum_i q_i / |x - a_i|,  x in R^3 minus the sites.  E = -grad V.
For positive charges every critical point lies in the convex hull of the sites.
This module is for exploration only; certificates live in certify.py.
"""
import numpy as np


def grad_hess(x, a, q):
    """Gradient and Hessian of V at points x (m,3). Returns g (m,3), H (m,3,3)."""
    d = x[:, None, :] - a[None, :, :]            # (m,n,3)
    r2 = np.einsum('mnk,mnk->mn', d, d)
    r = np.sqrt(r2)
    inv3 = q[None, :] / (r2 * r)                  # q/r^3
    inv5 = q[None, :] / (r2 * r2 * r)             # q/r^5
    g = -np.einsum('mn,mnk->mk', inv3, d)
    H = 3 * np.einsum('mn,mnk,mnl->mkl', inv5, d, d) - inv3.sum(1)[:, None, None] * np.eye(3)[None]
    return g, H


def potential(x, a, q):
    d = x[:, None, :] - a[None, :, :]
    return (q[None, :] / np.sqrt(np.einsum('mnk,mnk->mn', d, d))).sum(1)


def newton(x, a, q, iters=60, maxstep=None):
    """Damped Newton on grad V = 0 for many starting points at once."""
    scale = np.max(np.ptp(a, axis=0)) if len(a) > 1 else 1.0
    if maxstep is None:
        maxstep = 0.1 * scale
    x = x.copy()
    for _ in range(iters):
        g, H = grad_hess(x, a, q)
        try:
            step = np.linalg.solve(H, g[..., None])[..., 0]
        except np.linalg.LinAlgError:
            step = np.zeros_like(x)
            ok = np.abs(np.linalg.det(H)) > 1e-300
            step[ok] = np.linalg.solve(H[ok], g[ok][..., None])[..., 0]
        n = np.linalg.norm(step, axis=1)
        f = np.minimum(1.0, maxstep / np.maximum(n, 1e-300))
        x = x - step * f[:, None]
    return x


def find_critical_points(a, q, nseed=4000, rounds=3, seed=0, tol=1e-9, extra_seeds=None):
    """Multi-start Newton. Returns (points (k,3), index (k,), min |eig| (k,))."""
    a = np.asarray(a, float); q = np.asarray(q, float)
    rng = np.random.default_rng(seed)
    lo, hi = a.min(0), a.max(0)
    span = np.maximum(hi - lo, 1e-3 * np.max(hi - lo))
    lo = lo - 0.02 * span; hi = hi + 0.02 * span
    found = np.zeros((0, 3))
    for rnd in range(rounds):
        # seeds: random in bounding box, random convex combinations of sites, and jitter near found points
        s1 = lo + (hi - lo) * rng.random((nseed, 3))
        w = rng.dirichlet(np.ones(len(a)) * 0.7, size=nseed)
        s2 = w @ a
        seeds = [s1, s2]
        if len(found):
            sc = np.max(hi - lo)
            for eps in (1e-1, 1e-2, 1e-3, 1e-4):
                seeds.append(np.repeat(found, 20, 0) + eps * sc * rng.standard_normal((20 * len(found), 3)))
        if extra_seeds is not None and rnd == 0:
            seeds.append(np.asarray(extra_seeds, float))
        x = newton(np.vstack(seeds), a, q)
        g, H = grad_hess(x, a, q)
        # accept: small gradient relative to local field scale
        d = x[:, None, :] - a[None, :, :]
        rmin = np.sqrt(np.einsum('mnk,mnk->mn', d, d)).min(1)
        fieldscale = (np.abs(q)[None, :] / np.einsum('mnk,mnk->mn', d, d)).sum(1)
        ok = np.isfinite(x).all(1) & (np.linalg.norm(g, axis=1) < tol * fieldscale) & (rmin > 1e-6 * np.max(span))
        x = x[ok]
        allpts = np.vstack([found, x])
        found = dedup(allpts, 1e-7 * np.max(span))
    g, H = grad_hess(found, a, q)
    ev = np.linalg.eigvalsh(H)
    idx = (ev < 0).sum(1)
    return found, idx, np.abs(ev).min(1)


def dedup(pts, tol):
    if len(pts) == 0:
        return pts
    keep = []
    order = np.lexsort(pts.T[::-1])
    pts = pts[order]
    for p in pts:
        if not keep or np.min(np.linalg.norm(np.array(keep) - p, axis=1)) > tol:
            keep.append(p)
    return np.array(keep)


def summary(a, q, **kw):
    pts, idx, mev = find_critical_points(a, q, **kw)
    m1 = int((idx == 1).sum()); m2 = int((idx == 2).sum())
    return dict(n=len(pts), m1=m1, m2=m2, other=int(len(pts) - m1 - m2),
                morse_ok=(m2 - m1 == len(q) - 1), min_eig=float(mev.min()) if len(mev) else None), pts, idx


def _seeds(a, q, rng, nrand=400, nsite=60):
    n = len(a)
    S = [rng.dirichlet(np.full(n, 0.5), size=nrand) @ a]
    D = np.linalg.norm(a[:, None] - a[None], axis=2) + np.diag(np.full(n, np.inf))
    for i in range(n):
        dm = D[i].min()
        u = rng.standard_normal((nsite, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
        S.append(a[i] + u * (dm * 10 ** rng.uniform(-3, -0.05, nsite))[:, None])
        for j in range(i + 1, n):
            t0 = np.sqrt(q[i]) / (np.sqrt(q[i]) + np.sqrt(q[j]))
            t = np.clip(np.r_[t0, t0 * 10 ** rng.uniform(-1, 0, 8), 1 - (1 - t0) * 10 ** rng.uniform(-1, 0, 8)], 1e-6, 1 - 1e-6)
            S.append(a[i] + t[:, None] * (a[j] - a[i]))
    return np.vstack(S)


def newton_local(x, a, q, iters=50):
    """Damped Newton whose step never exceeds half the distance to the nearest site."""
    x = x.copy()
    for _ in range(iters):
        g, H = grad_hess(x, a, q)
        det = np.linalg.det(H)
        good = np.isfinite(det) & (np.abs(det) > 1e-300)
        step = np.zeros_like(x)
        step[good] = np.linalg.solve(H[good], g[good][..., None])[..., 0]
        rmin = np.sqrt(((x[:, None] - a[None]) ** 2).sum(2)).min(1)
        n = np.linalg.norm(step, axis=1)
        f = np.minimum(1.0, 0.5 * rmin / np.maximum(n, 1e-300))
        x = x - step * f[:, None]
    return x


def find_fast(a, q, seed=0, rounds=2, tol=1e-9):
    a = np.asarray(a, float); q = np.asarray(q, float)
    rng = np.random.default_rng(seed)
    sc = np.max(np.ptp(a, 0))
    found = np.zeros((0, 3))
    for rnd in range(rounds):
        X = _seeds(a, q, rng)
        if len(found):
            X = np.vstack([X] + [np.repeat(found, 10, 0) * 1 + e * rng.standard_normal((10 * len(found), 3)) *
                                 np.repeat(np.sqrt(((found[:, None] - a[None]) ** 2).sum(2)).min(1), 10)[:, None]
                                 for e in (0.3, 0.03)])
        X = newton_local(X, a, q)
        g, H = grad_hess(X, a, q)
        r2 = ((X[:, None] - a[None]) ** 2).sum(2)
        fs = (q[None] / r2).sum(1)
        ok = np.isfinite(X).all(1) & (np.linalg.norm(g, axis=1) < tol * fs) & (np.sqrt(r2).min(1) > 1e-9 * sc)
        found = dedup_r(np.vstack([found, X[ok]]), a)
    g, H = grad_hess(found, a, q)
    ev = np.linalg.eigvalsh(H) if len(found) else np.zeros((0, 3))
    return found, (ev < 0).sum(1), ev


def dedup_r(pts, a):
    """Deduplicate with a tolerance relative to the distance to the nearest site."""
    if len(pts) == 0:
        return pts
    rmin = np.sqrt(((pts[:, None] - a[None]) ** 2).sum(2)).min(1)
    keep = []
    for k in np.argsort(rmin):
        p = pts[k]
        if all(np.linalg.norm(p - pts[j]) > 1e-6 * min(rmin[k], rmin[j]) for j in keep):
            keep.append(k)
    return pts[keep]


def count_fast(a, q, seed=0):
    pts, idx, ev = find_fast(a, q, seed=seed)
    m1 = int((idx == 1).sum()); m2 = int((idx == 2).sum())
    ok = (m2 - m1 == len(q) - 1) and len(pts) == m1 + m2
    return len(pts), m1, m2, ok
