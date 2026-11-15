"""Balance points of equal charges at the vertices of the Platonic solids, and the order of the centre.
By Kelvin inversion the Taylor terms of the potential at a point are the multipole moments of the inverted charges,
so the centre of a spherical k-design of equal charges has all derivatives up to order k equal to zero."""
import numpy as np, itertools
from numpy.polynomial import legendre as L
phi = (1 + 5 ** 0.5) / 2
def solids():
    tet = np.array([(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)], float)
    octa = np.array([(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)], float)
    cube = np.array(list(itertools.product((-1,1), repeat=3)), float)
    ico = np.array([p for s1 in (-1,1) for s2 in (-1,1) for p in ((0, s1, s2*phi), (s1, s2*phi, 0), (s2*phi, 0, s1))], float)
    dod = np.array(list(itertools.product((-1,1), repeat=3)) + [p for s1 in (-1,1) for s2 in (-1,1) for p in ((0, s1/phi, s2*phi), (s1/phi, s2*phi, 0), (s2*phi, 0, s1/phi))], float)
    return {"tetrahedron": tet, "octahedron": octa, "cube": cube, "icosahedron": ico, "dodecahedron": dod}
def field(X, A, q):
    D = X[:, None, :] - A[None]; r2 = (D * D).sum(-1); r = np.sqrt(r2)
    g = -(q[None, :, None] * D / (r2 * r)[..., None]).sum(1)
    H = (q[None, :, None, None] * (3 * D[..., :, None] * D[..., None, :] - r2[..., None, None] * np.eye(3)) / (r2 * r2 * r)[..., None, None]).sum(1)
    return g, H
rng = np.random.default_rng(0)
for name, A in solids().items():
    A = A / np.linalg.norm(A[0]); n = len(A); q = np.ones(n)
    # order of the centre: largest over random directions of |sum P_k(e.u_i)| for each degree k
    E = rng.normal(size=(200, 3)); E /= np.linalg.norm(E, axis=1)[:, None]
    orders = [float(np.abs(sum(L.legval(E @ A[i], [0] * k + [1]) for i in range(n))).max()) for k in range(1, 9)]
    first = next(k for k, v in zip(range(1, 9), orders) if v > 1e-9)
    X = rng.normal(size=(40000, 3)); X *= (0.999 * rng.uniform(0, 1, size=(40000, 1)) ** (1 / 3)) / np.linalg.norm(X, axis=1)[:, None]
    with np.errstate(all="ignore"):
        for it in range(150):
            g, H = field(X, A, q)
            d = np.einsum('mab,mb->ma', np.linalg.pinv(H, rcond=1e-13), g)
            st = np.linalg.norm(d, axis=1); d *= np.minimum(1.0, 0.03 / np.maximum(st, 1e-300))[:, None]; X = X - d
        g, H = field(X, A, q)
    w = np.linalg.eigvalsh(np.where(np.isfinite(H), H, 0))
    ok = np.isfinite(X).all(1) & (np.linalg.norm(g, axis=1) < 1e-11 * n) & (np.linalg.norm(X, axis=1) < 1) & (np.abs(w).min(1) > 1e-6 * np.abs(w).max(1)) & (np.linalg.norm(X, axis=1) > 0.02)
    pts = []
    for x, ww in zip(X[ok], w[ok]):
        if all(np.linalg.norm(x - p[0]) > 1e-6 for p in pts): pts.append((x, int((ww < 0).sum())))
    m2 = sum(1 for p in pts if p[1] == 2); m1 = sum(1 for p in pts if p[1] == 1)
    radii = sorted(set(np.round([np.linalg.norm(p[0]) for p in pts], 4)))
    print(f"{name}: n = {n}; first non-zero Taylor degree at the centre = {first}, so the centre is {first - 1}-flat and unfolds into at most {(first - 1) ** 3} points")
    print(f"   other balance points: {len(pts)} (index 2: {m2}, index 1: {m1}), at radii {radii}; index-2 minus index-1 = {m2 - m1}, so the local degree of the centre = {n - 1 - (m2 - m1)}")
    print(f"   Maxwell's bound (n-1)^2 = {(n - 1) ** 2}; best known 10n-26 or 10n-31 = {10 * n - (26 if n % 2 else 31)}; centre + others could reach {(first - 1) ** 3 + len(pts)}")
