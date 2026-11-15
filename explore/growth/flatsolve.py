"""Positive charges whose potential is flat to order d-1 at the origin (all harmonic jets of degree 1..d-1 vanish).
Condition for degree l: sum_i q_i |a_i|^{-(2l+1)} p(a_i) = 0 for every harmonic polynomial p of degree l
(because the degree-l part of 1/|x-a| is |x|^l |a|^{-l-1} P_l(cos), and P_l(x.a) expands in degree-l harmonics).
Unknowns: sites a_i and log-charges. Underdetermined Gauss-Newton (least norm) from random starts. Numerical."""
import numpy as np
from harm import monos, evalmon, harmonic_basis

_B = {}


def basis(l):
    if l not in _B:
        B = harmonic_basis(l, l)
        L, _ = monos(l)
        _B[l] = B
    return _B[l]


def eqs(a, s, d):
    q = np.exp(s); r = np.linalg.norm(a, axis=1)
    out = []
    for l in range(1, d):
        M = evalmon(l, a) @ basis(l)        # (n, 2l+1): p(a_i)
        out.append((q * r ** (-(2 * l + 1))) @ M)
    return np.concatenate(out)


def jac(a, s, d, h=1e-7):
    x = np.r_[a.ravel(), s]
    f0 = eqs(a, s, d)
    J = np.zeros((len(f0), len(x)))
    for k in range(len(x)):
        xp = x.copy(); xp[k] += h
        xm = x.copy(); xm[k] -= h
        J[:, k] = (eqs(xp[:-len(s)].reshape(-1, 3), xp[-len(s):], d) - eqs(xm[:-len(s)].reshape(-1, 3), xm[-len(s):], d)) / (2 * h)
    return f0, J


def solve(n, d, rng, iters=200, a0=None, s0=None, smin=-6.0):
    a = rng.standard_normal((n, 3)) if a0 is None else a0.copy()
    a /= np.linalg.norm(a, axis=1)[:, None]
    a *= rng.uniform(0.7, 1.3, (n, 1)) if a0 is None else 1
    s = rng.uniform(-0.5, 0.5, n) if s0 is None else s0.copy()
    for it in range(iters):
        f, J = jac(a, s, d)
        nf = np.linalg.norm(f)
        if nf < 1e-13:
            break
        dx = np.linalg.lstsq(J, -f, rcond=None)[0]
        lam = 1.0
        for _ in range(20):
            an = a + lam * dx[:3 * n].reshape(-1, 3); sn = s + lam * dx[3 * n:]
            if np.linalg.norm(eqs(an, sn, d)) < nf:
                break
            lam *= 0.5
        a, s = an, sn
        # gauge: unit mean radius, max log-charge 0
        a /= np.linalg.norm(a, axis=1).mean()
        s -= s.max()
        if s.min() < smin or np.linalg.norm(a, axis=1).min() < 0.05:
            return None
    f = eqs(a, s, d)
    if np.linalg.norm(f) > 1e-11:
        return None
    return a, np.exp(s)


def top_poly(a, q, d):
    """Degree-d homogeneous part of sum q/|x - a| at 0, as monomial coefficients in MON[d] (exact up to rounding):
    projection of sum_i q_i |a_i|^{-(2d+1)} * (degree-d zonal kernel) computed by fitting on sample points."""
    L, ix = monos(d)
    sel = [k for k, m in enumerate(L) if sum(m) == d]
    rng = np.random.default_rng(5)
    X = rng.standard_normal((3 * len(sel), 3))
    from numpy.polynomial import legendre
    Pd = legendre.Legendre.basis(d)
    vals = np.zeros(len(X))
    rx = np.linalg.norm(X, axis=1)
    for ai, qi in zip(a, q):
        ra = np.linalg.norm(ai)
        cos = X @ ai / (rx * ra)
        vals += qi * rx ** d * Pd(cos) / ra ** (d + 1)
    M = evalmon(d, X)[:, sel]
    c = np.zeros(len(L))
    c[sel] = np.linalg.lstsq(M, vals, rcond=None)[0]
    return c
