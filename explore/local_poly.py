"""Local model behind ABK: Phi = L.X + Q(X) + C(X) + Z4(X), Z4 = 3 rho^4 - 24 rho^2 z^2 + 8 z^4.
grad Phi = 0 is a cubic system with no zeros at infinity, so 27 complex solutions (Bezout).
Count real ones; search (L, Q, C) to maximise."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np

def poly_grad_hess(X, L, Qm, Cc):
    x, y, z = X[..., 0], X[..., 1], X[..., 2]
    r2 = x * x + y * y
    # Z4 gradient and Hessian
    gx = 12 * r2 * x - 48 * x * z * z
    gy = 12 * r2 * y - 48 * y * z * z
    gz = -48 * r2 * z + 32 * z ** 3
    g = np.stack([gx, gy, gz], -1)
    H = np.zeros(X.shape + (3,), dtype=X.dtype)
    H[..., 0, 0] = 12 * r2 + 24 * x * x - 48 * z * z
    H[..., 1, 1] = 12 * r2 + 24 * y * y - 48 * z * z
    H[..., 2, 2] = -48 * r2 + 96 * z * z
    H[..., 0, 1] = H[..., 1, 0] = 24 * x * y
    H[..., 0, 2] = H[..., 2, 0] = -96 * x * z
    H[..., 1, 2] = H[..., 2, 1] = -96 * y * z
    # linear + quadratic (Qm symmetric traceless): grad = L + Qm X
    g = g + L + np.einsum('ij,...j->...i', Qm, X)
    H = H + Qm
    # cubic: C(X) = sum_{ijk} Cc[i,j,k] x_i x_j x_k with Cc symmetric
    g = g + 3 * np.einsum('ijk,...j,...k->...i', Cc, X, X)
    H = H + 6 * np.einsum('ijk,...k->...ij', Cc, X)
    return g, H

def all_roots(L, Qm, Cc, nseed=3000, seed=0):
    rng = np.random.default_rng(seed)
    X = (rng.standard_normal((nseed, 3)) + 1j * rng.standard_normal((nseed, 3))) * 10 ** rng.uniform(-1, 1, (nseed, 1))
    X = X.astype(complex)
    for _ in range(60):
        g, H = poly_grad_hess(X, L, Qm, Cc)
        try:
            st = np.linalg.solve(H, g[..., None])[..., 0]
        except np.linalg.LinAlgError:
            break
        X = X - st
    g, H = poly_grad_hess(X, L, Qm, Cc)
    ok = np.isfinite(X).all(1) & (np.abs(g).max(1) < 1e-9 * (1 + np.abs(X).max(1) ** 3))
    X = X[ok]
    roots = []
    for p in X:
        if not roots or np.min(np.abs(np.array(roots) - p).max(1)) > 1e-7 * (1 + np.abs(p).max()):
            roots.append(p)
    roots = np.array(roots)
    real = roots[np.abs(roots.imag).max(1) < 1e-9 * (1 + np.abs(roots).max(1))] if len(roots) else roots
    return len(roots), len(real), real.real

def harmonic_cubic(c):
    """Symmetric tensor for a cubic built from harmonic monomials with coefficients c (7)."""
    # basis of harmonic cubics in x,y,z (degree 3, dim 7)
    import itertools
    basis = [
        {(3, 0, 0): 1, (1, 2, 0): -3},            # x^3 - 3 x y^2
        {(2, 1, 0): 3, (0, 3, 0): -1},            # 3 x^2 y - y^3
        {(1, 0, 2): 4, (3, 0, 0): -1, (1, 2, 0): -1},   # x (4 z^2 - x^2 - y^2)
        {(0, 1, 2): 4, (2, 1, 0): -1, (0, 3, 0): -1},   # y (4 z^2 - x^2 - y^2)
        {(0, 0, 3): 2, (2, 0, 1): -3, (0, 2, 1): -3},   # z (2 z^2 - 3 x^2 - 3 y^2)
        {(2, 0, 1): 1, (0, 2, 1): -1},            # z (x^2 - y^2)
        {(1, 1, 1): 1},                           # x y z
    ]
    T = np.zeros((3, 3, 3))
    for ck, b in zip(c, basis):
        for (i, j, k), v in b.items():
            idx = [0] * i + [1] * j + [2] * k
            perms = set(itertools.permutations(idx))
            for p in perms:
                T[p] += ck * v / len(perms)
    return T

if __name__ == '__main__':
    # ABK local model check: Phi0 = 5/32 H2 + 15/8 H3 + 3/16 H4  (scale so Z4 coefficient is 1)
    Qm = (5 / 32) / (3 / 16) * 2 * np.diag([1, 1, -2]) * 1.0
    Cc = harmonic_cubic([15 / 8 / (3 / 16), 0, 0, 0, 0, 0, 0])
    print('ABK local model (complex, real):', all_roots(np.zeros(3), Qm, Cc)[:2])
