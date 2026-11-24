"""Independent verifier for certificates/*.json (no code shared with certify.py / ia.py).

Uses python-flint Arb balls for the field and Hessian and exact rationals (fmpq) for box
geometry and the dominance test. Checks, for a certificate claiming N zeros:
  1. each listed zero box passes a Krawczyk test with ||I - C J(X)||_inf < 1
     (a unique, non-degenerate zero) and its tight enclosure lies inside it;
  2. the tight enclosures are pairwise disjoint (so the N zeros are distinct);
  3. the Morse index of each zero (Gershgorin on Q^T J Q, Sylvester's law of inertia);
  4. completeness: an independent bisection of the padded bounding box of the sites in
     which every box is excluded (natural or mean-value enclosure of grad V misses 0, or one
     site dominates the field, decided in exact rational arithmetic) or lies inside one of
     the N verified zero boxes.
Zeros cannot lie outside the convex hull of positive charges, hence not outside the box.
"""
import json, sys, time
from fractions import Fraction
from flint import arb, fmpq, ctx

ctx.prec = 160


class VerificationError(Exception):
    """A proof step failed. Raised explicitly, so it is not removed by python -O (unlike assert)."""


def require(condition, message):
    if not condition:
        raise VerificationError(message)


def Q(x):
    f = Fraction(x) if not isinstance(x, str) else Fraction(x)
    return fmpq(f.numerator, f.denominator)


def ball(lo, hi):
    return arb((lo + hi) / 2, (hi - lo) / 2)


def sq(x):
    lo, hi = x.lower(), x.upper()
    s = (lo * lo).union(hi * hi)
    if lo <= 0 and hi >= 0 or x.contains(0):
        s = s.union(arb(0))
    return s


class Config:
    def __init__(self, sites, charges):
        require(type(sites) is list and type(charges) is list, 'sites and charges must be lists')
        require(len(sites) == len(charges) and len(charges) >= 1, 'sites and charges must be non-empty and same length')
        for p in sites:
            require(type(p) is list and len(p) == 3, 'each site must be a 3-element list')
        self.a = [[_safe_Q(c, 'site coordinate') for c in p] for p in sites]
        self.q = [_safe_Q(c, 'charge') for c in charges]
        require(all(c > 0 for c in self.q), 'all charges must be strictly positive')
        self.aa = [[arb(c) for c in p] for p in self.a]
        self.qa = [arb(c) for c in self.q]
        self.n = len(self.q)

    def grad(self, X):
        g = [arb(0), arb(0), arb(0)]
        for i in range(self.n):
            d = [X[k] - self.aa[i][k] for k in range(3)]
            r2 = sq(d[0]) + sq(d[1]) + sq(d[2])
            if not r2 > 0:
                return None
            inv3 = self.qa[i] / (r2 * r2.sqrt())
            for k in range(3):
                g[k] -= inv3 * d[k]
        return g

    def hess(self, X):
        H = [[arb(0), arb(0), arb(0)] for _ in range(3)]
        for i in range(self.n):
            d = [X[k] - self.aa[i][k] for k in range(3)]
            s = [sq(dk) for dk in d]
            r2 = s[0] + s[1] + s[2]
            if not r2 > 0:
                return None
            inv5 = self.qa[i] / (r2 * r2 * r2.sqrt())
            for k in range(3):
                o = [m for m in range(3) if m != k]
                H[k][k] += inv5 * (2 * s[k] - s[o[0]] - s[o[1]])
                for l in range(k + 1, 3):
                    H[k][l] += inv5 * 3 * d[k] * d[l]
        for k in range(3):
            for l in range(k):
                H[k][l] = H[l][k]
        return H

    def dominated(self, lo, hi):
        """Exact rational test: some site's field outweighs all others on the box."""
        dmin, dmax = [], []
        for i in range(self.n):
            mn = fmpq(0); mx = fmpq(0)
            for k in range(3):
                L = lo[k] - self.a[i][k]; U = hi[k] - self.a[i][k]
                if L <= 0 <= U:
                    m_ = fmpq(0)
                else:
                    m_ = min(abs(L), abs(U))
                M_ = max(abs(L), abs(U))
                mn += m_ * m_; mx += M_ * M_
            dmin.append(mn); dmax.append(mx)
        for i in range(self.n):
            if any(dmin[j] == 0 for j in range(self.n) if j != i):
                continue
            strong = self.q[i] / dmax[i]
            weak = sum((self.q[j] / dmin[j] for j in range(self.n) if j != i), fmpq(0))
            if strong > weak:
                return True
        return False


def mid_float(b):
    return float(b.mid())


def inv3x3(M):
    a, b, c = M[0]; d, e, f = M[1]; g, h, i = M[2]
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    if det == 0:
        return None
    return [[(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det],
            [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
            [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det]]


def krawczyk(cfg, lo, hi):
    X = [ball(lo[k], hi[k]) for k in range(3)]
    c = [(lo[k] + hi[k]) / 2 for k in range(3)]
    gc = cfg.grad([arb(ck) for ck in c])
    J = cfg.hess(X)
    if gc is None or J is None:
        return False, None
    Jm = [[fmpq(Fraction(mid_float(J[i][j])).numerator, Fraction(mid_float(J[i][j])).denominator) for j in range(3)] for i in range(3)]
    C = inv3x3(Jm)
    if C is None:
        return False, None
    Ca = [[arb(C[i][j]) for j in range(3)] for i in range(3)]
    K = []
    norm = arb(0)
    ok = True
    for i in range(3):
        acc = arb(c[i]) - sum((Ca[i][k] * gc[k] for k in range(3)), arb(0))
        row = arb(0)
        for j in range(3):
            Mij = (arb(1) if i == j else arb(0)) - sum((Ca[i][k] * J[k][j] for k in range(3)), arb(0))
            acc += Mij * (X[j] - arb(c[j]))
            row += abs(Mij)
        if not row < 1:
            ok = False
        K.append(acc)
        if not (acc.lower() > arb(lo[i]) and acc.upper() < arb(hi[i])):
            ok = False
    return ok, K


def morse_index(cfg, lo, hi):
    X = [ball(lo[k], hi[k]) for k in range(3)]
    J = cfg.hess(X)
    import numpy as np
    Jm = np.array([[mid_float(J[i][j]) for j in range(3)] for i in range(3)])
    w, V = np.linalg.eigh(Jm)
    Vq = [[fmpq(Fraction(float(V[i][j])).numerator, Fraction(float(V[i][j])).denominator) for j in range(3)] for i in range(3)]
    # Sylvester needs Vq non-singular: check det exactly
    det = (Vq[0][0] * (Vq[1][1] * Vq[2][2] - Vq[1][2] * Vq[2][1]) - Vq[0][1] * (Vq[1][0] * Vq[2][2] - Vq[1][2] * Vq[2][0])
           + Vq[0][2] * (Vq[1][0] * Vq[2][1] - Vq[1][1] * Vq[2][0]))
    if det == 0:
        return None
    Va = [[arb(x) for x in r] for r in Vq]
    B = [[sum((Va[k][i] * J[k][l] * Va[l][j] for k in range(3) for l in range(3)), arb(0)) for j in range(3)] for i in range(3)]
    neg = 0
    for i in range(3):
        rad = sum((abs(B[i][j]) for j in range(3) if j != i), arb(0))
        if (B[i][i] + rad) < 0:
            neg += 1
        elif (B[i][i] - rad) > 0:
            pass
        else:
            return None
    return neg


def verify(path):
    cert = _load_cert_json(path)
    cfg = Config(cert['sites'], cert['charges'])
    t0 = time.time()
    boxes = []
    for z in cert['zeros']:
        lo = [Q(x) for x in z['box'][0]]; hi = [Q(x) for x in z['box'][1]]
        ok, K = krawczyk(cfg, lo, hi)
        require(ok, 'Krawczyk failed for %s' % z)
        # Every zero in a box X lies in K(X); iterate X <- X cap K(X) (rounded outward)
        # to get our own tight enclosure E of the unique zero of the box.
        elo, ehi = lo, hi
        for _ in range(80):
            _, K = krawczyk(cfg, elo, ehi)
            if K is None:
                break
            nlo = [max(elo[k], Q(_dn(K[k]))) for k in range(3)]
            nhi = [min(ehi[k], Q(_up(K[k]))) for k in range(3)]
            require(all(nlo[k] <= nhi[k] for k in range(3)), 'empty intersection: no zero?')
            if nlo == elo and nhi == ehi:
                break
            elo, ehi = nlo, nhi
        idx = morse_index(cfg, elo, ehi)
        require(idx == z['index'], 'index mismatch %s vs %s' % (idx, z['index']))
        boxes.append((lo, hi, elo, ehi, idx))
    # distinct zeros: the enclosures are pairwise disjoint
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            e1, e2 = boxes[i], boxes[j]
            require(any(e1[3][k] < e2[2][k] or e2[3][k] < e1[2][k] for k in range(3)), 'enclosures %d, %d overlap' % (i, j))
    print('  %d zero boxes verified (%.1fs)' % (len(boxes), time.time() - t0), flush=True)
    # completeness
    import math
    pts = [[Q(c) for c in p] for p in cert['sites']]
    span = max(max(p[k] for p in pts) - min(p[k] for p in pts) for k in range(3))
    pad = span / 256
    lo0 = [min(p[k] for p in pts) - pad for k in range(3)]
    hi0 = [max(p[k] for p in pts) + pad for k in range(3)]
    stack = [(lo0, hi0)]
    leaves = dict(natural=0, centred=0, dominance=0, inside=0)
    nbox = 0
    while stack:
        lo, hi = stack.pop()
        nbox += 1
        if nbox > 3_000_000:
            raise RuntimeError('budget')
        if any(all(b[0][k] <= lo[k] and hi[k] <= b[1][k] for k in range(3)) for b in boxes):
            leaves['inside'] += 1; continue
        X = [ball(lo[k], hi[k]) for k in range(3)]
        g = cfg.grad(X)
        if g is not None and any(gk > 0 or gk < 0 for gk in g):
            leaves['natural'] += 1; continue
        if g is not None:
            c = [(lo[k] + hi[k]) / 2 for k in range(3)]
            gc = cfg.grad([arb(ck) for ck in c]); J = cfg.hess(X)
            if gc is not None and J is not None:
                gm = [gc[i] + sum((J[i][j] * (X[j] - arb(c[j])) for j in range(3)), arb(0)) for i in range(3)]
                if any(gk > 0 or gk < 0 for gk in gm):
                    leaves['centred'] += 1; continue
        if cfg.dominated(lo, hi):
            leaves['dominance'] += 1; continue
        w = [hi[k] - lo[k] for k in range(3)]
        k = max(range(3), key=lambda m: w[m])
        if w[k] < span / 2 ** 60:
            raise RuntimeError('cannot resolve box near %s' % [float(x) for x in lo])
        m = (lo[k] + hi[k]) / 2
        hi1 = list(hi); hi1[k] = m; lo2 = list(lo); lo2[k] = m
        stack.append((lo, hi1)); stack.append((lo2, hi))
    require(len(boxes) == cert['count'], 'certificate claims %s zeros but %d boxes were verified' % (cert['count'], len(boxes)))
    idx = [b[4] for b in boxes]
    m1, m2 = idx.count(1), idx.count(2)
    print('  completeness: %d boxes %s (%.1fs)' % (nbox, leaves, time.time() - t0))
    print('VERIFIED %s: exactly %d zeros, all non-degenerate; index 1: %d, index 2: %d; m2-m1 = %d = n-1: %s'
          % (cert['name'], len(boxes), m1, m2, m2 - m1, m2 - m1 == cfg.n - 1))
    return len(boxes), m1, m2


def _dn(b):
    import math
    return math.nextafter(float(b.lower()), -math.inf)


def _up(b):
    import math
    return math.nextafter(float(b.upper()), math.inf)


def _reject_dup_keys(pairs):
    out = {}
    for k, v in pairs:
        require(k not in out, 'duplicate JSON key: %r' % (k,))
        out[k] = v
    return out


def _reject_nonfinite(token):
    raise VerificationError('non-finite JSON constant forbidden: %r' % (token,))


def _safe_Q(x, label):
    require(type(x) is not bool, 'bool forbidden in %s' % label)
    try:
        return Q(x)
    except (ValueError, TypeError, OverflowError) as exc:
        raise VerificationError('invalid numeric %s: %s' % (label, exc)) from exc


def _load_cert_json(path):
    try:
        with open(path, 'r', encoding='utf-8') as f:
            cert = json.load(f, object_pairs_hook=_reject_dup_keys, parse_constant=_reject_nonfinite)
    except VerificationError:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, RecursionError, ValueError) as exc:
        raise VerificationError('failed to load certificate JSON: %s' % exc) from exc
    require(type(cert) is dict, 'certificate root must be a JSON object')
    for key in ('name', 'sites', 'charges', 'zeros', 'count'):
        require(key in cert, 'certificate missing required key: %s' % key)
    require(type(cert['name']) is str, 'certificate name must be str')
    require(type(cert['count']) is int and cert['count'] >= 0, 'certificate count must be a non-negative int')
    require(type(cert['zeros']) is list, 'certificate zeros must be a list')
    for z in cert['zeros']:
        require(type(z) is dict and 'box' in z and 'index' in z, 'each zero entry must be a dict with box and index')
        require(type(z['index']) is int, 'zero Morse index must be exact int')
        box = z['box']
        require(type(box) is list and len(box) == 2, 'zero box must be [lo, hi]')
        for end in box:
            require(type(end) is list and len(end) == 3, 'box bound must be a 3-element list')
        for k in range(3):
            require(_safe_Q(box[0][k], 'box bound') <= _safe_Q(box[1][k], 'box bound'), 'inverted zero box bound at axis %d' % k)
    return cert


if __name__ == '__main__':
    for p in sys.argv[1:]:
        print(p)
        verify(p)
