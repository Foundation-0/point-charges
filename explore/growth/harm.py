"""Harmonic polynomials in R^3: basis, gradient, and all complex critical points by homotopy.

Numerical exploration only (not proof). A polynomial of degree <= d is a coefficient vector over the
monomial list MON[d]. Critical points of P solve grad P = 0, a square system of degree d-1.
solve_all() tracks the (d-1)^3 paths of a total-degree homotopy and returns the finite endpoints.
"""
import itertools
import numpy as np

_cache = {}


def monos(d):
    if d not in _cache:
        L = [(i, j, k) for t in range(d + 1) for i in range(t + 1) for j in range(t + 1 - i) for k in [t - i - j]]
        _cache[d] = (L, {m: n for n, m in enumerate(L)})
    return _cache[d]


def deriv_mat(d, axis):
    L, ix = monos(d)
    D = np.zeros((len(L), len(L)))
    for n, m in enumerate(L):
        if m[axis] > 0:
            mm = list(m); mm[axis] -= 1
            D[ix[tuple(mm)], n] = m[axis]
    return D


def lap_mat(d):
    return sum(deriv_mat(d, a) @ deriv_mat(d, a) for a in range(3))


def harmonic_basis(d, l):
    """Orthonormal basis (columns, monomial coefficients in MON[d]) of homogeneous harmonics of degree l."""
    L, ix = monos(d)
    sel = [n for n, m in enumerate(L) if sum(m) == l]
    Lp = lap_mat(d)[:, sel]
    u, s, vt = np.linalg.svd(Lp)
    r = (s > 1e-9).sum()
    N = vt[r:].T
    B = np.zeros((len(L), N.shape[1]))
    B[sel] = N
    assert B.shape[1] == 2 * l + 1
    return B


def evalmon(d, X):
    """Monomial values, X (m,3) real or complex -> (m, M)."""
    L, _ = monos(d)
    E = np.array(L)
    pw = [np.stack([X[:, a] ** e for e in range(d + 1)], 1) for a in range(3)]
    return pw[0][:, E[:, 0]] * pw[1][:, E[:, 1]] * pw[2][:, E[:, 2]]


class Poly:
    def __init__(self, c, d):
        self.c = np.asarray(c); self.d = d
        self.D = [deriv_mat(d, a) for a in range(3)]
        self.g = [Dm @ self.c for Dm in self.D]
        self.h = [[self.D[a] @ self.g[b] for b in range(3)] for a in range(3)]

    def val(self, X):
        return evalmon(self.d, X) @ self.c

    def grad(self, X):
        M = evalmon(self.d, X)
        return np.stack([M @ gg for gg in self.g], 1)

    def gradhess(self, X):
        M = evalmon(self.d, X)
        g = np.stack([M @ gg for gg in self.g], 1)
        H = np.stack([np.stack([M @ self.h[a][b] for b in range(3)], 1) for a in range(3)], 1)
        return g, H


def solve_all(P, rng=None, steps0=200, tol_inf=1e7):
    """Total-degree homotopy for grad P = 0. Returns finite complex solutions (refined) and #lost paths."""
    rng = np.random.default_rng(0) if rng is None else rng
    e = P.d - 1
    roots = np.exp(2j * np.pi * np.arange(e) / e)
    S = np.array(list(itertools.product(roots, repeat=3)))  # starts of x_i^e = 1
    gam = np.exp(2j * np.pi * rng.random())

    def HJ(X, t):
        g, Hs = P.gradhess(X)
        G = X ** e - 1
        JG = np.zeros((len(X), 3, 3), complex)
        for a in range(3):
            JG[:, a, a] = e * X[:, a] ** (e - 1)
        Hval = (1 - t) * gam * G + t * g
        J = ((1 - t) * gam)[..., None] * JG + t[..., None] * Hs
        dHdt = -gam * G + g
        return Hval, J, dHdt

    X = S.astype(complex).copy()
    t = np.zeros(len(X)); dt = np.full(len(X), 1.0 / steps0)
    alive = np.ones(len(X), bool)
    for it in range(200000):
        act = alive & (t < 1)
        if not act.any():
            break
        ia = np.where(act)[0]
        x = X[ia]; tt = t[ia]; h = np.minimum(dt[ia], 1 - tt)
        # RK4-ish predictor via Euler with two half steps (Heun)
        Hv, J, Ht = HJ(x, tt[:, None])
        try:
            v1 = -np.linalg.solve(J, Ht[..., None])[..., 0]
        except np.linalg.LinAlgError:
            v1 = np.zeros_like(x)
        xp = x + h[:, None] * v1
        Hv, J, Ht = HJ(xp, (tt + h)[:, None])
        try:
            v2 = -np.linalg.solve(J, Ht[..., None])[..., 0]
        except np.linalg.LinAlgError:
            v2 = v1
        xp = x + h[:, None] * 0.5 * (v1 + v2)
        tn = tt + h
        ok = np.ones(len(ia), bool)
        for _ in range(3):
            Hv, J, Ht = HJ(xp, tn[:, None])
            try:
                dx = np.linalg.solve(J, Hv[..., None])[..., 0]
            except np.linalg.LinAlgError:
                dx = np.full_like(xp, np.inf)
            xp = xp - dx
        sz = np.maximum(1, np.abs(xp).max(1))
        ok = np.isfinite(xp).all(1) & (np.abs(dx).max(1) < 1e-6 * sz)
        X[ia[ok]] = xp[ok]; t[ia[ok]] = tn[ok]
        dt[ia[ok]] = np.minimum(dt[ia[ok]] * 1.5, 0.05)
        dt[ia[~ok]] *= 0.5
        dead = (dt < 1e-13) | (np.abs(X).max(1) > tol_inf)
        alive &= ~dead
    fin = alive & (t >= 1)
    Xf = X[fin]
    # final Newton polish on P itself
    for _ in range(8):
        g, H = P.gradhess(Xf)
        try:
            Xf = Xf - np.linalg.solve(H, g[..., None])[..., 0]
        except np.linalg.LinAlgError:
            break
    g, H = P.gradhess(Xf)
    good = np.isfinite(Xf).all(1) & (np.abs(g).max(1) < 1e-8 * np.maximum(1, np.abs(Xf).max(1)) ** (P.d - 1))
    Xf = Xf[good]
    # dedup
    keep = []
    for k in range(len(Xf)):
        if all(np.abs(Xf[k] - Xf[j]).max() > 1e-6 * max(1, np.abs(Xf[k]).max()) for j in keep):
            keep.append(k)
    return Xf[keep], int((~fin).sum())


def real_crit(P, sols, rtol=1e-7):
    r = np.abs(sols.imag).max(1) < rtol * np.maximum(1, np.abs(sols).max(1))
    R = sols[r].real
    if len(R) == 0:
        return R, np.zeros(0, int), np.zeros(0)
    g, H = P.gradhess(R)
    ev = np.linalg.eigvalsh(H)
    return R, (ev < 0).sum(1), np.abs(ev).min(1)


def newton_real(P, X, iters=60):
    for _ in range(iters):
        g, H = P.gradhess(X)
        det = np.linalg.det(H)
        okk = np.abs(det) > 1e-300
        st = np.zeros_like(X)
        st[okk] = np.linalg.solve(H[okk], g[okk][..., None])[..., 0]
        n = np.linalg.norm(st, axis=1)
        f = np.minimum(1, 0.3 * np.maximum(1, np.linalg.norm(X, axis=1)) / np.maximum(n, 1e-300))
        X = X - st * f[:, None]
    return X
