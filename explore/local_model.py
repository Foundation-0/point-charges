"""Local model: one unit charge at the origin in a background b.X + X^T K X / 2.
Equilibria solve X/|X|^3 = K X + b. Count them over (K, b)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np
def count(K, b, nseed=4000, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((nseed, 3)) * 10 ** rng.uniform(-1.5, 1.5, (nseed, 1))
    for _ in range(100):
        r = np.linalg.norm(X, axis=1)
        F = X / r[:, None]**3 - X @ K.T - b
        J = (np.eye(3)[None] / r[:, None, None]**3 - 3 * X[:, :, None] * X[:, None, :] / r[:, None, None]**5) - K[None]
        st = np.linalg.solve(J, F[..., None])[..., 0]
        n = np.linalg.norm(st, axis=1); f = np.minimum(1, 0.3 * r / np.maximum(n, 1e-300))
        X = X - st * f[:, None]
    r = np.linalg.norm(X, axis=1)
    F = X / r[:, None]**3 - X @ K.T - b
    ok = np.isfinite(X).all(1) & (np.linalg.norm(F, axis=1) < 1e-10 * (1 / r**2 + 1))
    pts = []
    for p in X[ok]:
        if not pts or np.min(np.linalg.norm(np.array(pts) - p, axis=1)) > 1e-7 * (1 + np.linalg.norm(p)):
            pts.append(p)
    return len(pts)
if __name__ == '__main__':
    rng = np.random.default_rng(0)
    from collections import Counter
    C = Counter(); best = 0
    for t in range(3000):
        kap = rng.uniform(-1, 1)          # eigenvalues (1, kap, -1-kap): covers both signatures
        K = np.diag([1.0, kap, -1 - kap])
        b = rng.standard_normal(3) * 10 ** rng.uniform(-2, 0.5)
        c = count(K, b, nseed=600, seed=t)
        C[c] += 1
        if c > best:
            best = c; print('new best', c, kap, b, flush=True)
    print(C)
