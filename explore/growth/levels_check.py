"""Separate check of two counts claimed in the growth note (results/2026-11-15-growth-n-three-halves.md), with programs written apart from the search.

1. One level. A = zonal harmonic of degree p about the z axis, B = zonal harmonic of degree p - 1 about the x axis.
   Claim: A + B has exactly I_p - 1 real non-degenerate critical points away from the origin, with
   I_p = floor(p^2 / 2) + 2 [p = 3 mod 4].
   Method here: a critical point x = t w (w a unit vector) needs grad A(w) parallel to grad B(w), and then
   t = -(grad A . grad B) / |grad A|^2 is fixed. So the critical points away from 0 are the directions, up to sign,
   where the cross product of the two gradients vanishes. Those are found by Newton on the sphere from many starts.
2. Stacked levels. P = sum_{k=1..D} c_k Z_k with the axes alternating z, x and the levels at radius ratio s.
   Claim: P has exactly sum_{p=2..D} (I_p - 1) real critical points: 6, 13, 24, 41, 66 for D = 3..7.
   Method here: Newton in three dimensions in extended precision from starts spread over every scale.
   Claimed and checked here for D = 3, 4 and 5 only. At D = 6 this program's search reported 59 points against 41 (the
   extra points are not understood and are not claimed). The separate program hier.py finds 41 at D = 6.
Every mismatch within the claimed range stops the program with a nonzero exit status."""
import numpy as np, sys
def die(m): raise SystemExit('FAIL: ' + m)
from numpy.polynomial import legendre as Lg
rng = np.random.default_rng(0)
I = lambda p: p * p // 2 + (2 if p % 4 == 3 else 0)
def zonal_gh(X, k, axis):
    """gradient and Hessian of r^k P_k(x_axis / r) at the rows of X (k >= 1)."""
    r = np.linalg.norm(X, axis=1); u = X[:, axis] / r
    c = [0] * k + [1]; P = Lg.legval(u, c); dP = Lg.legval(u, Lg.legder(c)); d2P = Lg.legval(u, Lg.legder(c, 2)) if k >= 2 else 0 * u
    e = np.zeros(3); e[axis] = 1.0; n = X / r[:, None]
    du = (e[None] - u[:, None] * n) / r[:, None]                       # gradient of u
    g = (k * r ** (k - 1) * P)[:, None] * n + (r ** k * dP)[:, None] * du
    d2u = (-(du[:, :, None] * n[:, None, :]) - (n[:, :, None] * du[:, None, :]) - u[:, None, None] * (np.eye(3)[None] - n[:, :, None] * n[:, None, :]) / r[:, None, None]) / r[:, None, None]
    dn = (np.eye(3)[None] - n[:, :, None] * n[:, None, :]) / r[:, None, None]
    H = (k * (k - 1) * r ** (k - 2) * P)[:, None, None] * n[:, :, None] * n[:, None, :] + (k * r ** (k - 1) * P)[:, None, None] * dn \
        + (k * r ** (k - 1) * dP)[:, None, None] * (n[:, :, None] * du[:, None, :] + du[:, :, None] * n[:, None, :]) \
        + (r ** k * d2P)[:, None, None] * du[:, :, None] * du[:, None, :] + (r ** k * dP)[:, None, None] * d2u
    return g, H
print("check 1: one level")
for p in range(3, 13):
    W = rng.normal(size=(30000, 3)); W /= np.linalg.norm(W, axis=1)[:, None]
    def F(W):
        gA, HA = zonal_gh(W, p, 2); gB, HB = zonal_gh(W, p - 1, 0)
        return np.cross(gA, gB), gA, gB, HA, HB
    for it in range(60):
        c, gA, gB, HA, HB = F(W)
        # Jacobian of the cross product, then a least-squares step kept on the sphere
        J = np.cross(HA.transpose(0, 2, 1), gB[:, None, :]).transpose(0, 2, 1) + np.cross(gA[:, None, :], HB.transpose(0, 2, 1)).transpose(0, 2, 1)
        Jc = np.concatenate([J, W[:, None, :]], axis=1); rhs = np.concatenate([c, np.zeros((len(W), 1))], axis=1)
        step = np.einsum('mab,mb->ma', np.linalg.pinv(Jc), rhs)
        nrm = np.linalg.norm(step, axis=1); step *= np.minimum(1.0, 0.1 / np.maximum(nrm, 1e-300))[:, None]
        W = W - step; W /= np.linalg.norm(W, axis=1)[:, None]
    c, gA, gB, HA, HB = F(W)
    ok = np.linalg.norm(c, axis=1) < 1e-10 * np.linalg.norm(gA, axis=1) * np.linalg.norm(gB, axis=1) + 1e-13
    ok &= (np.linalg.norm(gA, axis=1) > 1e-7) & (np.linalg.norm(gB, axis=1) > 1e-7)
    dirs = []
    for w in W[ok]:
        if all(min(np.linalg.norm(w - d), np.linalg.norm(w + d)) > 1e-5 for d in dirs): dirs.append(w)
    dirs = np.array(dirs); gA, HA = zonal_gh(dirs, p, 2); gB, HB = zonal_gh(dirs, p - 1, 0)
    t = -(gA * gB).sum(1) / (gA * gA).sum(1); X = dirs * t[:, None]
    gA, HA = zonal_gh(X, p, 2); gB, HB = zonal_gh(X, p - 1, 0); H = HA + HB; ev = np.abs(np.linalg.eigvalsh(H))
    resid = np.linalg.norm(gA + gB, axis=1) / (np.linalg.norm(gA, axis=1) + 1e-300)
    idx = (np.linalg.eigvalsh(H) < 0).sum(1)
    if len(dirs) != I(p) - 1: die(f'p = {p}: found {len(dirs)} directions, claimed {I(p) - 1}')
    print(f"  p = {p:2d}: directions found {len(dirs)}; claimed I_p - 1 = {I(p) - 1}; largest relative residual {resid.max():.1e}; smallest eigenvalue ratio {(ev.min(1) / ev.max(1)).min():.3f}; index 2: {(idx == 2).sum()}, index 1: {(idx == 1).sum()}")
print("check 2: stacked levels (extended precision Newton)")
import mpmath as mp
mp.mp.dps = 60
def stacked(D, s):
    lv = [(k, 2 if (D - k) % 2 == 0 else 0, mp.mpf(s) ** ((D - k) * (D - k + 1) // 2)) for k in range(1, D + 1)]
    def P(x, y, z):
        r = mp.sqrt(x * x + y * y + z * z); tot = mp.mpf(0)
        for k, ax, c in lv: tot += c * r ** k * mp.legendre(k, (x, y, z)[ax] / r)
        return tot
    return P, lv
for D in range(3, 6):
    s = 0.01; P, lv = stacked(D, s)
    g = lambda x, y, z: [mp.diff(P, (x, y, z), d) for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    found = []
    for level in range(0, D - 1):
        rad = s ** level
        for trial in range(260):
            v = rng.normal(size=3); v *= rad * 10 ** rng.uniform(-1.2, 0.6) / np.linalg.norm(v)
            try: sol = mp.findroot(lambda x, y, z: g(x, y, z), tuple(mp.mpf(float(t)) for t in v), tol=1e-40, maxsteps=40)
            except Exception: continue
            x = np.array([float(t) for t in sol]); nx = np.linalg.norm(x)
            if not np.isfinite(nx) or nx < 1e-30 or nx > 50: continue
            if all(np.linalg.norm(x - f) > 1e-6 * max(nx, np.linalg.norm(f)) for f in found): found.append(x)
    radii = np.sort([np.linalg.norm(f) for f in found])
    groups = [int(((radii > s ** (j + 0.5)) & (radii <= s ** (j - 0.5))).sum()) for j in range(0, D)]
    claimed = sum(I(p) - 1 for p in range(2, D + 1))
    print(f"  D = {D}, ratio {s}: real critical points found {len(found)}; claimed {claimed}; by scale (outermost first) {groups}", flush=True)
    if len(found) != claimed: die(f'D = {D}: found {len(found)}, claimed {claimed}')
print('ALL CHECKS PASSED')
