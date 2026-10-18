"""Three positive charges: equilibria lie in the open triangle, in the plane of the charges.
#equilibria = 2*m0 + 2 with m0 the number of local minima of the in-plane potential.
Maxwell's bound 4 <=> m0 <= 1. GNNS (2026) prove m0 <= 2."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np

def crit2d(A, z, nseed=600, seed=0, iters=60):
    rng = np.random.default_rng(seed)
    w = rng.dirichlet(np.full(3, 0.6), size=nseed)
    # extra seeds near each charge, scaled by charge ratios
    S = [w @ A]
    for i in range(3):
        for j in range(3):
            if i != j:
                t0 = np.sqrt(z[i]) / (np.sqrt(z[i]) + np.sqrt(z[j]))
                t = np.clip(t0 * 10 ** rng.uniform(-1.5, 0.3, 20), 1e-7, 1 - 1e-7)
                S.append(A[i] + t[:, None] * (A[j] - A[i]) + 1e-3 * rng.standard_normal((20, 2)) * np.linalg.norm(A[j] - A[i]) * t[:, None])
    X = np.vstack(S)
    for _ in range(iters):
        d = X[:, None, :] - A[None]
        r2 = (d * d).sum(2); r = np.sqrt(r2)
        i3 = z / (r2 * r); i5 = z / (r2 * r2 * r)
        g = -(i3[..., None] * d).sum(1)
        H = 3 * np.einsum('mn,mnk,mnl->mkl', i5, d, d) - i3.sum(1)[:, None, None] * np.eye(2)
        det = H[:, 0, 0] * H[:, 1, 1] - H[:, 0, 1] ** 2
        det = np.where(np.abs(det) < 1e-300, 1e-300, det)
        st = np.c_[H[:, 1, 1] * g[:, 0] - H[:, 0, 1] * g[:, 1], H[:, 0, 0] * g[:, 1] - H[:, 0, 1] * g[:, 0]] / det[:, None]
        n = np.linalg.norm(st, axis=1); f = np.minimum(1, 0.5 * r.min(1) / np.maximum(n, 1e-300))
        X = X - st * f[:, None]
    d = X[:, None, :] - A[None]; r2 = (d * d).sum(2); r = np.sqrt(r2)
    i3 = z / (r2 * r); i5 = z / (r2 * r2 * r)
    g = -(i3[..., None] * d).sum(1)
    ok = np.isfinite(X).all(1) & (np.linalg.norm(g, axis=1) < 1e-9 * (z / r2).sum(1)) & (r.min(1) > 1e-9)
    X = X[ok]; rm = r.min(1)[ok]
    pts = []; rr = []
    for k in np.argsort(rm):
        if all(np.linalg.norm(X[k] - pts[j]) > 1e-6 * min(rm[k], rr[j]) for j in range(len(pts))):
            pts.append(X[k]); rr.append(rm[k])
    pts = np.array(pts).reshape(-1, 2)
    d = pts[:, None, :] - A[None]; r2 = (d * d).sum(2); r = np.sqrt(r2)
    i3 = z / (r2 * r); i5 = z / (r2 * r2 * r)
    H = 3 * np.einsum('mn,mnk,mnl->mkl', i5, d, d) - i3.sum(1)[:, None, None] * np.eye(2)
    ev = np.linalg.eigvalsh(H) if len(pts) else np.zeros((0, 2))
    m0 = int((ev > 0).all(1).sum()); ms = int(((ev[:, 0] < 0) & (ev[:, 1] > 0)).sum())
    return pts, m0, ms

if __name__ == '__main__':
    import sys, time
    from collections import Counter
    rng = np.random.default_rng(int(sys.argv[2]))
    C = Counter(); t0 = time.time(); hits = []
    while time.time() - t0 < float(sys.argv[1]):
        a = rng.uniform(-1.5, 2.5); b = 10 ** rng.uniform(-3, 0.7)
        A = np.array([[0, 0], [1, 0], [a, b]])
        z = np.r_[10 ** rng.uniform(-4, 4, 2), 1.0]
        pts, m0, ms = crit2d(A, z, seed=int(rng.integers(1 << 30)))
        ok = ms - m0 == 2
        C[(m0, ok)] += 1
        if ok and m0 >= 2:
            print('HIT m0=%d' % m0, A.tolist(), z.tolist(), flush=True); hits.append((A.tolist(), z.tolist()))
    print(sorted(C.items()))
