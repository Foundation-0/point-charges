"""Rigorous count of the equilibria of a configuration of positive point charges.

Claim certified: inside R^3 minus the sites, grad V has exactly N zeros, each one
non-degenerate (non-singular Hessian), with the stated Morse indices.

Method (all arithmetic in ia.py, outward-rounded float64 intervals):
  * Every zero lies in the convex hull of the sites (positive charges), so in the
    bounding box B of the sites, padded.
  * B is bisected. A box is discarded when
      - (natural interval extension) some component of grad V excludes 0 on it, or
      - (dominance) for some site i, q_i / rmax_i^2 > sum_{j != i} q_j / rmin_j^2, so the
        field of site i outweighs all others everywhere in the box (this also handles
        boxes containing a site, where V is not defined).
  * A box that cannot be discarded is inflated and tested with the Krawczyk operator
    K = c - C g(c) + (I - C J(X))(X - c). K inside int(X) proves a unique zero in X;
    in addition ||I - C J(X)||_inf < 1 proves every Hessian in J(X) non-singular, so
    the zero is non-degenerate. Otherwise the box is bisected again.
  * Verified boxes that overlap are matched up through tight enclosures of their zeros.
  * The Morse index of each zero comes from Sylvester's law of inertia applied to
    Q^T J(X) Q with Gershgorin discs (Q = eigenvectors of mid J, any non-singular Q works).
"""
import numpy as np, json, sys, time
from fractions import Fraction
import ia

# ---------- interval kernels (boxes: lo, hi arrays of shape (m,3)) ----------

def _parts(lo, hi, a):
    d = [(ia.dn(lo[:, None, k] - a[None, :, k]), ia.up(hi[:, None, k] - a[None, :, k])) for k in range(3)]
    s = [ia.sqr(dk) for dk in d]
    r2 = ia.add(ia.add(s[0], s[1]), s[2])
    return d, s, r2


def _sum1(t):
    # Sum over sites with outward rounding: naive float sum error is bounded by
    # (n-1) * eps * sum|t|; we add that bound explicitly.
    n = t[0].shape[1]
    lo = t[0].sum(1); hi = t[1].sum(1)
    e = (n + 1) * 2.0 ** -52 * np.maximum(np.abs(t[0]), np.abs(t[1])).sum(1)
    return (ia.dn(lo - e), ia.up(hi + e))


def grad_iv(lo, hi, a, q):
    """Interval enclosure of grad V over boxes: list of 3 intervals of shape (m,)."""
    d, s, r2 = _parts(lo, hi, a)
    r = ia.sqrt(r2)
    inv3 = ia.scale(q[None, :], ia.recip(ia.mul(r2, r)))
    return [ia.neg(_sum1(ia.mul(inv3, d[k]))) for k in range(3)]


def hess_iv(lo, hi, a, q):
    """Interval Hessian of V: sum q (3 d d^T - r^2 I) / r^5. Returns 3x3 nested list of intervals."""
    d, s, r2 = _parts(lo, hi, a)
    r = ia.sqrt(r2)
    inv5 = ia.scale(q[None, :], ia.recip(ia.mul(ia.mul(r2, r2), r)))
    H = [[None] * 3 for _ in range(3)]
    for k in range(3):
        o = [m for m in range(3) if m != k]
        diag = ia.sub(ia.scale(2.0, s[k]), ia.add(s[o[0]], s[o[1]]))
        H[k][k] = _sum1(ia.mul(inv5, diag))
        for l in range(k + 1, 3):
            H[k][l] = _sum1(ia.mul(inv5, ia.scale(3.0, ia.mul(d[k], d[l]))))
            H[l][k] = H[k][l]
    return H


def dominance(lo, hi, a, q):
    """Boolean (m,): some site dominates the field everywhere in the box."""
    m, n = lo.shape[0], a.shape[0]
    # per-axis min and max distance from each site to the box
    dmin = np.zeros((m, n)); dmax = np.zeros((m, n))
    for k in range(3):
        L = ia.dn(lo[:, None, k] - a[None, :, k]); U = ia.up(hi[:, None, k] - a[None, :, k])
        inside = (L <= 0) & (U >= 0)
        mn = np.where(inside, 0.0, np.minimum(np.abs(L), np.abs(U)))
        mx = np.maximum(np.abs(L), np.abs(U))
        dmin = ia.dn(dmin + ia.dn(mn * mn))
        dmax = ia.up(dmax + ia.up(mx * mx))
    # dmin, dmax now hold squared distances (lower / upper bounds)
    with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
        strong = ia.dn(q[None, :] / dmax)                     # lower bound of q_i / r_i^2
        weak = np.where(dmin > 0, ia.up(q[None, :] / dmin), np.inf)   # upper bound of q_j / r_j^2
    # others[:, i] bounds sum_{j != i} q_j / r_j^2 from above. All terms are >= 0, so the
    # float sum has relative error at most (n - 2) * 2^-53; the factor below covers it.
    others = np.empty((m, n))
    for i in range(n):
        others[:, i] = ia.up(np.delete(weak, i, axis=1).sum(1) * (1 + (n + 2) * 2.0 ** -52))
    return (strong > others).any(1)


# ---------- small interval linear algebra ----------

def matvec_point(C, v):
    """C (m,3,3) exact floats times interval vector v (list of 3 intervals (m,))."""
    out = []
    for i in range(3):
        acc = None
        for j in range(3):
            t = ia.scale(C[:, i, j], v[j])
            acc = t if acc is None else ia.add(acc, t)
        out.append(acc)
    return out


def centred_grad(lo, hi, a, q):
    """Mean-value enclosure g(c) + J(X)(X - c) of grad V over boxes."""
    c = np.clip(0.5 * (lo + hi), lo, hi)
    gc = grad_iv(c, c, a, q)
    J = hess_iv(lo, hi, a, q)
    xc = [(ia.dn(lo[:, k] - c[:, k]), ia.up(hi[:, k] - c[:, k])) for k in range(3)]
    out = []
    for i in range(3):
        acc = gc[i]
        for j in range(3):
            acc = ia.add(acc, ia.mul(J[i][j], xc[j]))
        out.append(acc)
    return out


def krawczyk(lo, hi, a, q):
    """Krawczyk test on boxes. Returns (ok (m,), Klo, Khi)."""
    m = lo.shape[0]
    c = 0.5 * (lo + hi)
    c = np.clip(c, lo, hi)
    gc = grad_iv(c, c, a, q)
    J = hess_iv(lo, hi, a, q)
    Jm = np.zeros((m, 3, 3))
    for i in range(3):
        for j in range(3):
            Jm[:, i, j] = 0.5 * (J[i][j][0] + J[i][j][1])
    ok = np.isfinite(Jm).all((1, 2))
    Jm[~ok] = np.eye(3)
    det = np.linalg.det(Jm)
    ok &= np.abs(det) > 1e-300
    Jm[~ok] = np.eye(3)
    C = np.linalg.inv(Jm)
    Cg = matvec_point(C, gc)
    # M = I - C J  (interval matrix)
    M = [[None] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            acc = None
            for k in range(3):
                t = ia.scale(C[:, i, k], J[k][j])
                acc = t if acc is None else ia.add(acc, t)
            e = 1.0 if i == j else 0.0
            M[i][j] = ia.sub((np.full(m, e), np.full(m, e)), acc)
    xc = [(ia.dn(lo[:, k] - c[:, k]), ia.up(hi[:, k] - c[:, k])) for k in range(3)]
    Klo = np.zeros((m, 3)); Khi = np.zeros((m, 3))
    norm = np.zeros(m)
    for i in range(3):
        acc = ia.sub((c[:, i], c[:, i]), Cg[i])
        rs = np.zeros(m)
        for j in range(3):
            acc = ia.add(acc, ia.mul(M[i][j], xc[j]))
            rs = ia.up(rs + ia.mag(M[i][j]))
        norm = np.maximum(norm, rs)
        Klo[:, i], Khi[:, i] = acc
    ok &= np.isfinite(Klo).all(1) & np.isfinite(Khi).all(1)
    ok &= (Klo > lo).all(1) & (Khi < hi).all(1) & (norm < 1)
    return ok, Klo, Khi


def morse_index(lo, hi, a, q):
    """Index (number of negative eigenvalues) for every Hessian in J(box), or -1 if undecided."""
    J = hess_iv(lo[None], hi[None], a, q)
    Jm = np.array([[0.5 * (J[i][j][0][0] + J[i][j][1][0]) for j in range(3)] for i in range(3)])
    w, Q = np.linalg.eigh(Jm)
    # B = Q^T J Q in interval arithmetic
    B = [[None] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            acc = None
            for k in range(3):
                for l in range(3):
                    t = ia.scale(Q[k, i] * 1.0, ia.scale(Q[l, j] * 1.0, (J[k][l][0], J[k][l][1])))
                    acc = t if acc is None else ia.add(acc, t)
            B[i][j] = (acc[0][0], acc[1][0])
    neg = 0
    for i in range(3):
        rad = 0.0
        for j in range(3):
            if j != i:
                rad = ia.up(rad + max(abs(B[i][j][0]), abs(B[i][j][1])))
        if ia.up(B[i][i][1] + rad) < 0:
            neg += 1
        elif ia.dn(B[i][i][0] - rad) > 0:
            pass
        else:
            return -1
    # Q must be non-singular for Sylvester's law: |det Q - 1| <= ||Q^T Q - I|| bound, checked loosely
    G = Q.T @ Q
    if np.abs(G - np.eye(3)).max() > 1e-6:
        return -1
    return neg


# ---------- main subdivision ----------

def certify(a, q, pad=1e-3, max_boxes=5_000_000, verbose=True, batch_size=8192):
    a = np.asarray(a, float); q = np.asarray(q, float)
    assert (q > 0).all()
    span = np.max(np.ptp(a, 0))
    lo0 = a.min(0) - pad * span; hi0 = a.max(0) + pad * span
    if type(batch_size) is not int or batch_size < 1:
        raise ValueError("batch_size must be a positive integer")
    pending = [(lo0[None].copy(), hi0[None].copy())]
    # LIFO batches bound temporary interval arrays while retaining every child.
    # All exclusion and Krawczyk arithmetic below is unchanged.
    verified = []          # list of (lo, hi) verified boxes
    stats = dict(boxes=0, excl_interval=0, excl_dominance=0, krawczyk=0)
    t0 = time.time()
    while pending:
        lo, hi = pending.pop()
        if len(lo) > batch_size:
            pending.append((lo[batch_size:], hi[batch_size:]))
            lo, hi = lo[:batch_size], hi[:batch_size]
        m = len(lo)
        stats['boxes'] += m
        if stats['boxes'] > max_boxes:
            raise RuntimeError('box budget exceeded')
        dom = dominance(lo, hi, a, q)
        g = grad_iv(lo, hi, a, q)
        ex = np.zeros(m, bool)
        for k in range(3):
            ex |= ~ia.contains_zero(g[k])
        # mean-value form g(X) in g(c) + J(X)(X - c) for the boxes still undecided
        und = ~(ex | dom)
        if und.any():
            gm = centred_grad(lo[und], hi[und], a, q)
            exm = np.zeros(int(und.sum()), bool)
            for k in range(3):
                exm |= ~ia.contains_zero(gm[k])
            ex[np.where(und)[0][exm]] = True
        stats['excl_dominance'] += int((dom & ~ex).sum()); stats['excl_interval'] += int(ex.sum())
        keep = ~(ex | dom)
        lo, hi = lo[keep], hi[keep]
        if not len(lo):
            continue
        # Krawczyk on inflated boxes that are small compared with their distance to the sites
        w = hi - lo
        cen = 0.5 * (lo + hi)
        rmin = np.sqrt(((cen[:, None] - a[None]) ** 2).sum(2)).min(1)
        tryk = w.max(1) < 0.25 * rmin
        kv = np.zeros(len(lo), bool)
        if tryk.any():
            ilo = lo[tryk] - 0.1 * w[tryk]; ihi = hi[tryk] + 0.1 * w[tryk]
            ok, _, _ = krawczyk(ilo, ihi, a, q)
            idx = np.where(tryk)[0][ok]
            kv[idx] = True
            for l_, h_ in zip(ilo[ok], ihi[ok]):
                verified.append((l_.copy(), h_.copy()))
            stats['krawczyk'] += int(ok.sum())
        lo, hi = lo[~kv], hi[~kv]
        w = hi - lo
        if len(lo) and w.max(1).min() < 1e-13 * span:
            raise RuntimeError('box too small without resolution near %s' % (0.5 * (lo + hi))[np.argmin(w.max(1))])
        # bisect along widest dimension
        k = np.argmax(w, 1)
        mid = 0.5 * (lo[np.arange(len(lo)), k] + hi[np.arange(len(lo)), k])
        lo1, hi1 = lo.copy(), hi.copy(); lo2, hi2 = lo.copy(), hi.copy()
        hi1[np.arange(len(lo)), k] = mid; lo2[np.arange(len(lo)), k] = mid
        if len(lo):
            pending.append((np.vstack([lo1, lo2]), np.vstack([hi1, hi2])))
    if verbose:
        print(f"subdivision done: {stats['boxes']} boxes, {len(verified)} Krawczyk boxes, {time.time()-t0:.1f}s", flush=True)
    zeros = merge(verified, a, q)
    return zeros, stats


def tighten(lo, hi, a, q, iters=30):
    for _ in range(iters):
        ok, Klo, Khi = krawczyk(lo[None], hi[None], a, q)
        if not ok[0]:
            break
        nlo, nhi = np.maximum(lo, Klo[0]), np.minimum(hi, Khi[0])
        if np.all(nlo == lo) and np.all(nhi == hi):
            break
        lo, hi = nlo, nhi
    return lo, hi


def merge(verified, a, q):
    """Group verified boxes by the zero they contain."""
    tight = [tighten(l, h, a, q) for l, h in verified]
    n = len(verified)
    parent = list(range(n))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    L = np.array([v[0] for v in verified]); H = np.array([v[1] for v in verified])
    TL = np.array([t[0] for t in tight]); TH = np.array([t[1] for t in tight])
    for i in range(n):
        inter = ((L <= H[i]) & (H >= L[i])).all(1)
        inter[i] = False
        for j in np.where(inter)[0]:
            if j < i:
                continue
            a_in_b = ((TL[i] >= L[j]) & (TH[i] <= H[j])).all()
            b_in_a = ((TL[j] >= L[i]) & (TH[j] <= H[i])).all()
            a_out_b = ((TL[i] > H[j]) | (TH[i] < L[j])).any()
            b_out_a = ((TL[j] > H[i]) | (TH[j] < L[i])).any()
            if a_in_b or b_in_a:
                parent[find(i)] = find(j)
            elif a_out_b or b_out_a:
                pass
            else:
                raise RuntimeError('ambiguous overlap between verified boxes %d and %d' % (i, j))
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    zeros = []
    for g in groups.values():
        i = min(g, key=lambda k: np.max(TH[k] - TL[k]))
        idx = morse_index(TL[i], TH[i], a, q)
        zeros.append(dict(box=[L[i].tolist(), H[i].tolist()], enclosure=[TL[i].tolist(), TH[i].tolist()], index=int(idx), boxes=len(g)))
    zeros.sort(key=lambda z: z['enclosure'][0])
    return zeros


def exact(x):
    return str(Fraction(float(x)))


def run(a, q, name, out=None, **kw):
    a = np.asarray(a, float); q = np.asarray(q, float)
    zeros, stats = certify(a, q, **kw)
    idx = [z['index'] for z in zeros]
    res = dict(name=name,
               claim=f"grad V has exactly {len(zeros)} zeros in R^3 minus the sites; all non-degenerate",
               sites=[[exact(x) for x in p] for p in a], charges=[exact(x) for x in q],
               sites_float=a.tolist(), charges_float=q.tolist(),
               count=len(zeros), index_counts={str(k): idx.count(k) for k in sorted(set(idx))},
               morse_relation_m2_minus_m1=idx.count(2) - idx.count(1), n_minus_1=len(q) - 1,
               stats=stats,
               zeros=zeros)
    print(name, '->', res['count'], 'zeros; indices', res['index_counts'], flush=True)
    if out:
        json.dump(res, open(out, 'w'), indent=1)
    return res
