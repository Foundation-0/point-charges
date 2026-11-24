"""Multi-scale numerical finder for the balance points of positive point charges (exploration, not proof).

V(x) = sum q_i/|x - a_i|.  Seeds are placed at every charge's own length scale (distance to its nearest
neighbour), along every pair of near neighbours, at triple and quadruple centroids of near neighbours,
and at random convex combinations.  Damped Newton never steps more than half the distance to the nearest
site.  Converged points are deduplicated with a tolerance relative to their distance to the nearest site
(so nested clusters at scales 1e-6 or 1e-12 are handled).  Reports the counts by index and two Morse checks:
  (M1) m2 - m1 = n - 1, and (M2) the wall network (uphill branches of each index-2 point) has
  walls - n + components independent loops, which must equal m1.
Failure modes: (a) a seed family may miss points that sit at a scale or place no seed reaches;
(b) near-degenerate points (min |eigenvalue| tiny) may be lost or merged; (c) cancellation in float64
when the configuration's own scales span more than about 1e-8 of its size; (d) both Morse checks can still
pass if points are missed in index-1/index-2 pairs (M1) and the missed index-2 points close no loop (M2 is
then violated only when the missing pair changes the loop count, so M2 catches most but not all).
"""
import numpy as np
from scipy.spatial import cKDTree
from flint import arb, ctx

HP = 212  # bits for refinement


def _gh_hp(x, A, Q):
    """Gradient and Hessian at one point in high precision (arb midpoints). A, Q: arb lists."""
    g = [arb(0)] * 3; H = [[arb(0)] * 3 for _ in range(3)]
    for (a0, a1, a2), qi in zip(A, Q):
        d = (x[0] - a0, x[1] - a1, x[2] - a2)
        r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
        r = r2.sqrt(); i3 = qi / (r2 * r); i5 = i3 / r2
        for k in range(3):
            g[k] -= i3 * d[k]
            for l in range(k, 3):
                H[k][l] += 3 * i5 * d[k] * d[l] - (i3 if k == l else 0)
    for k in range(3):
        for l in range(k):
            H[k][l] = H[l][k]
    return g, H


def refine_hp(X, a, q, iters=6, a_exact=None, q_exact=None):
    """Newton in high precision. Returns refined float points and the final relative Newton step."""
    old = ctx.prec; ctx.prec = HP
    try:
        A = a_exact if a_exact is not None else [[arb(float(v)) for v in row] for row in a]
        Q = q_exact if q_exact is not None else [arb(float(v)) for v in q]
        out = np.zeros_like(X); rel = np.zeros(len(X))
        for m, x0 in enumerate(X):
            x = [arb(float(v)) for v in x0]
            for it in range(iters):
                g, H = _gh_hp(x, A, Q)
                Hf = np.array([[float(H[i][j].mid()) for j in range(3)] for i in range(3)])
                try:
                    Hi = np.linalg.inv(Hf)
                except np.linalg.LinAlgError:
                    break
                # step = Hi g, carried in arb with float inverse (Newton-like, converges linearly then fast)
                st = [sum((arb(float(Hi[i, j])) * g[j] for j in range(3)), arb(0)).mid() for i in range(3)]
                x = [(x[i] - st[i]).mid() for i in range(3)]
            rmin = min(float(((x[0] - a0) ** 2 + (x[1] - a1) ** 2 + (x[2] - a2) ** 2).sqrt().mid()) for a0, a1, a2 in A)
            out[m] = [float(v.mid()) for v in x]
            rel[m] = float(sum(s * s for s in st).sqrt().mid()) / rmin
        return out, rel
    finally:
        ctx.prec = old


def grad_hess(x, a, q):
    d = x[:, None, :] - a[None, :, :]
    r2 = np.einsum('mnk,mnk->mn', d, d)
    r = np.sqrt(r2)
    inv3 = q[None, :] / (r2 * r)
    inv5 = inv3 / r2
    g = -np.einsum('mn,mnk->mk', inv3, d)
    H = 3 * np.einsum('mn,mnk,mnl->mkl', inv5, d, d) - inv3.sum(1)[:, None, None] * np.eye(3)[None]
    return g, H


def grad_only(x, a, q):
    d = x[:, None, :] - a[None, :, :]
    r2 = np.einsum('mnk,mnk->mn', d, d)
    inv3 = q[None, :] / (r2 * np.sqrt(r2))
    return -np.einsum('mn,mnk->mk', inv3, d), r2


def fib_sphere(k):
    i = np.arange(k) + 0.5
    phi = np.arccos(1 - 2 * i / k); th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)


def seeds(a, q, rng, ndir=40, nrad=8, knn=8, nrand=2000):
    n = len(a)
    S = []
    tree = cKDTree(a)
    k = min(knn + 1, n)
    dist, nb = tree.query(a, k=k)
    if n > 1:
        dn = dist[:, 1]
    else:
        dn = np.ones(1)
    U = fib_sphere(ndir)
    for i in range(n):
        R = rng.random((3, 3)); Qm, _ = np.linalg.qr(R)
        Ui = U @ Qm
        for s in np.geomspace(0.03, 1.5, nrad):
            S.append(a[i] + s * dn[i] * Ui)
    # pairs among near neighbours
    pairs = set()
    for i in range(n):
        for j in nb[i, 1:]:
            if j < n:
                pairs.add((min(i, j), max(i, j)))
    for i, j in pairs:
        t0 = np.sqrt(q[i]) / (np.sqrt(q[i]) + np.sqrt(q[j]))
        ts = np.r_[t0, np.linspace(0.05, 0.95, 7)]
        L = np.linalg.norm(a[j] - a[i])
        base = a[i] + ts[:, None] * (a[j] - a[i])
        S.append(base)
        off = rng.standard_normal((len(ts) * 4, 3)) * 0.15 * L
        S.append(np.repeat(base, 4, 0) + off)
    # triples / quadruples of near neighbours
    for i in range(n):
        nbs = [j for j in nb[i, 1:min(6, k)] if j < n]
        for x in range(len(nbs)):
            for y in range(x + 1, len(nbs)):
                c3 = (a[i] + a[nbs[x]] + a[nbs[y]]) / 3
                S.append(c3[None] + rng.standard_normal((3, 3)) * 0.1 * dn[i])
                for z in range(y + 1, len(nbs)):
                    c4 = (a[i] + a[nbs[x]] + a[nbs[y]] + a[nbs[z]]) / 4
                    S.append(c4[None] + rng.standard_normal((2, 3)) * 0.1 * dn[i])
    w = rng.dirichlet(np.full(n, 0.4), size=nrand)
    S.append(w @ a)
    return np.vstack(S)


def newton(x, a, q, iters=80):
    x = x.copy()
    for _ in range(iters):
        g, H = grad_hess(x, a, q)
        det = np.linalg.det(H)
        good = np.isfinite(det) & (np.abs(det) > 0)
        step = np.zeros_like(x)
        step[good] = np.linalg.solve(H[good], g[good][..., None])[..., 0]
        rmin = np.sqrt(((x[:, None] - a[None]) ** 2).sum(2)).min(1)
        nst = np.linalg.norm(step, axis=1)
        f = np.minimum(1.0, 0.5 * rmin / np.maximum(nst, 1e-300))
        x = x - step * f[:, None]
    return x


def accept(x, a, q, tol):
    """Scale-free test: the Newton step is below tol times the distance to the nearest site."""
    g, H = grad_hess(x, a, q)
    r2 = ((x[:, None] - a[None]) ** 2).sum(2)
    rmin = np.sqrt(r2.min(1))
    ok = np.isfinite(x).all(1) & (rmin > 0) & np.isfinite(g).all(1) & np.isfinite(H).all((1, 2))
    st = np.full(len(x), np.inf)
    if ok.any():
        Hk = H[ok]; det = np.linalg.det(Hk); good = np.abs(det) > 0
        sv = np.full(ok.sum(), np.inf)
        sv[good] = np.linalg.norm(np.linalg.solve(Hk[good], g[ok][good][..., None])[..., 0], axis=1)
        st[ok] = sv
    ok &= st < tol * rmin
    return ok, rmin


def dedup(pts, a, rel=1e-6):
    if len(pts) == 0:
        return pts
    rmin = np.sqrt(((pts[:, None] - a[None]) ** 2).sum(2)).min(1)
    tree = cKDTree(pts)
    par = np.arange(len(pts))

    def find(i):
        while par[i] != i:
            par[i] = par[par[i]]; i = par[i]
        return i
    for i, nbrs in enumerate(tree.query_ball_point(pts, rel * rmin)):
        for j in nbrs:
            ri, rj = find(i), find(j)
            if ri != rj:
                par[ri] = rj
    roots = np.array([find(i) for i in range(len(pts))])
    _, first = np.unique(roots, return_index=True)
    return pts[np.sort(first)]


def find(a, q, seed=0, rounds=3, tol=1e-9, extra=None, hp=None, a_exact=None, q_exact=None, **kw):
    """hp: None = decide automatically (use high-precision refinement when the charges span > 1e6 in size
    or sites span > 1e4 in spacing)."""
    a = np.asarray(a, float); q = np.asarray(q, float)
    rng = np.random.default_rng(seed)
    if hp is None:
        dn = cKDTree(a).query(a, k=2)[0][:, 1] if len(a) > 1 else np.ones(1)
        hp = (q.max() / q.min() > 1e6) or (dn.max() / dn.min() > 1e4)
    loose = 3e-2 if hp else tol
    found = np.zeros((0, 3))
    for rnd in range(rounds):
        X = seeds(a, q, rng, **kw)
        if extra is not None and rnd == 0:
            X = np.vstack([X, extra])
        if len(found):
            rm = np.sqrt(((found[:, None] - a[None]) ** 2).sum(2)).min(1)
            X = np.vstack([X] + [np.repeat(found, 6, 0) + e * rng.standard_normal((6 * len(found), 3)) * np.repeat(rm, 6)[:, None]
                                 for e in (0.5, 0.1, 0.02)])
        X = X[np.sqrt(((X[:, None] - a[None]) ** 2).sum(2)).min(1) > 0]
        with np.errstate(all='ignore'):
            X = newton(X, a, q)
            ok, _ = accept(X, a, q, loose)
        cand = dedup(X[ok], a, rel=1e-5)
        if hp and len(cand):
            new = dedup(np.vstack([found, cand]), a, rel=1e-5)
            # refine only candidates not already in found
            if len(found):
                tr = cKDTree(found); rmn = np.sqrt(((new[:, None] - a[None]) ** 2).sum(2)).min(1)
                dd, _ = tr.query(new)
                new = new[dd > 1e-5 * rmn]
            if len(new):
                Xr, rel = refine_hp(new, a, q, iters=10, a_exact=a_exact, q_exact=q_exact)
                cand = Xr[rel < 1e-12]
            else:
                cand = np.zeros((0, 3))
        found = dedup(np.vstack([found, cand]), a)
    if hp and len(found):
        old = ctx.prec; ctx.prec = HP
        A = a_exact if a_exact is not None else [[arb(float(v)) for v in row] for row in a]
        Q = q_exact if q_exact is not None else [arb(float(v)) for v in q]
        ev = []
        for x0 in found:
            _, H = _gh_hp([arb(float(v)) for v in x0], A, Q)
            ev.append(np.linalg.eigvalsh(np.array([[float(H[i][j].mid()) for j in range(3)] for i in range(3)])))
        ctx.prec = old
        ev = np.array(ev)
    else:
        g, H = grad_hess(found, a, q)
        ev = np.linalg.eigvalsh(H) if len(found) else np.zeros((0, 3))
    idx = (ev < 0).sum(1)
    return found, idx, ev


def dominance_radii(a, q):
    """Largest rho (by bisection) with q_i/rho^2 > sum_j q_j/(d_ij - rho)^2: inside, the flow goes straight to a_i."""
    n = len(a); out = np.zeros(n)
    D = np.sqrt(((a[:, None] - a[None]) ** 2).sum(2))
    for i in range(n):
        dj = np.delete(D[i], i); qj = np.delete(q, i)
        if len(dj) == 0:
            out[i] = 1.0; continue
        lo, hi = 0.0, 0.5 * dj.min()
        for _ in range(60):
            m = 0.5 * (lo + hi)
            if q[i] / m ** 2 > (qj / (dj - m) ** 2).sum():
                lo = m
            else:
                hi = m
        out[i] = 0.9 * lo
    return out


def trace_walls(pts, idx, ev_unused, a, q, maxsteps=4000):
    """Follow both uphill branches of each index-2 point. Returns list of (i, j) charge pairs (or None)."""
    sel = np.where(idx == 2)[0]
    if len(sel) == 0:
        return []
    P = pts[sel]
    g, H = grad_hess(P, a, q)
    w, V = np.linalg.eigh(H)
    v = V[:, :, 2]  # the positive eigenvalue's eigenvector (uphill)
    rmin0 = np.sqrt(((P[:, None] - a[None]) ** 2).sum(2)).min(1)
    X = np.vstack([P + 1e-3 * rmin0[:, None] * v, P - 1e-3 * rmin0[:, None] * v])
    m = len(X)
    dest = -np.ones(m, int)
    # capture radius: a ball where the charge's own field dominates, conservatively 0.2 x half nearest distance
    tree = cKDTree(a)
    dn = tree.query(a, k=2)[0][:, 1] if len(a) > 1 else np.ones(1)
    cap = dominance_radii(a, q)
    act = np.ones(m, bool)
    for it in range(maxsteps):
        ia = np.where(act)[0]
        if len(ia) == 0:
            break
        x = X[ia]
        dd, j = tree.query(x)
        hit = dd < cap[j]
        dest[ia[hit]] = j[hit]
        act[ia[hit]] = False
        ia = ia[~hit]; x = x[~hit]; dd = dd[~hit]
        if len(ia) == 0:
            break
        h = 0.05 * dd
        g, _ = grad_only(x, a, q)
        gn = np.linalg.norm(g, axis=1)[:, None]
        k1 = g / gn
        xm = x + 0.5 * h[:, None] * k1
        g2, _ = grad_only(xm, a, q)
        k2 = g2 / np.linalg.norm(g2, axis=1)[:, None]
        X[ia] = x + h[:, None] * k2
    half = len(P)
    out = []
    for k in range(half):
        i, j = dest[k], dest[k + half]
        out.append(None if i < 0 or j < 0 else (min(i, j), max(i, j)))
    return out


def wall_stats(edges, n):
    good = [e for e in edges if e is not None]
    par = list(range(n))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for i, j in good:
        par[f(i)] = f(j)
    comp = len({f(i) for i in range(n)})
    pairs = len({e for e in good if e[0] != e[1]})
    return dict(walls=len(good), failed=len(edges) - len(good), pairs=pairs, components=comp,
                loops=len(good) - n + comp)


def report(a, q, seed=0, rounds=3, walls=True, **kw):
    a = np.asarray(a, float); q = np.asarray(q, float)
    pts, idx, ev = find(a, q, seed=seed, rounds=rounds, **kw)
    n = len(q)
    m1 = int((idx == 1).sum()); m2 = int((idx == 2).sum())
    r = dict(n=n, total=len(pts), m1=m1, m2=m2, other=len(pts) - m1 - m2, M1=(m2 - m1 == n - 1),
             min_abs_eig_rel=float((np.abs(ev).min(1) / np.abs(ev).max(1)).min()) if len(pts) else None)
    if walls:
        e = trace_walls(pts, idx, ev, a, q)
        ws = wall_stats(e, n)
        r.update(ws)
        r['M2'] = (ws['failed'] == 0 and ws['loops'] == m1)
    return r, pts, idx
