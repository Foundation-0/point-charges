#!/usr/bin/env python3
"""Independent count of the critical points of V(x) = sum_i q_i/|x - a_i| (all q_i > 0).

A third checker, written separately from certify.py / ia.py and from verify_arb.py: it shares no code
with them. Its method differs from verify_arb.py (see results/2026-11-15-third-checker.md).
Input: a JSON file with "sites" and "charges" as exact fraction strings. Nothing else in the
file is used for the proof, so a certificate's boxes are not read at all. (The "count" field, if present,
is compared at the end only.) Needs python-flint and NumPy.
Usage: python3 -I -O indep_count.py <certificate or configuration .json> [out.json]
Exit status: 0 the count agrees with the file's "count" (or the file has none), 2 it does not, 3 a proof step failed.

Proof structure (each step is proved in results/2026-11-15-count-lemmas.md; the derivative bounds are also
checked numerically by tests/check_lemmas.py):
  L1 (far field). If a half-space {u.x <= b} contains every site, then grad V != 0 wherever u.x > b.
      Half-spaces are built from triples of sites and checked in exact rationals.
      The root box is the bounding box of the sites, which contains the closed convex hull.
  L2 (near each charge). If q_i/rho^2 > sum_{j!=i} q_j/(d_ij - rho)^2 then grad V != 0 on
      0 < |x - a_i| <= rho. rho_i is a dyadic rational, checked in exact rationals using a
      rational lower bound for d_ij from an integer square root.
  T1 (first-order exclusion). |g(c)| > R * sum 2 q_i / delta_i^3, where R is the half-diagonal
      of the box, delta_i the distance from a_i to the box. Uses ||Hess(1/r)|| = 2/r^3.
  T2 (second-order exclusion). For some vector u: |u.g(c)| - R |H(c) u| > (R^2/2)|u| sum 6 q_i/delta_i^4.
      Uses |D^3(1/r)[v,v,.]| <= 6/r^4 and ||D^3(1/r)[v,.,.]|| <= 6/r^4 for unit v.
  Z (zero balls). Around a numerical zero x0 and radius rho, with A ~ H(x0)^{-1}:
      kappa = ||I - A H(x0)||_F + ||A||_F * rho * sum 6 q_i/(|x0-a_i| - rho)^4 < 1 and
      |A g(x0)| <= (1 - kappa) rho  ==>  T(x) = x - A g(x) is a kappa-contraction of the closed ball
      into itself, so g has exactly one zero there, and H is non-singular on the whole ball.
      Since trace H = 0 (V harmonic), the Morse index is 2 if det H > 0 and 1 if det H < 0.
  Balls are checked pairwise disjoint in exact rationals. An octree subdivision of the root box
  discards every box by L1, L2, T1, T2 or containment in a zero ball. Then the number of zeros
  of grad V in R^3 minus the sites equals the number of balls.

All proof inequalities are evaluated with python-flint arb (rigorous ball arithmetic; a comparison
returns True only when it holds for every point of the balls) or with exact fmpq/Fraction
arithmetic. Floats are used only to propose candidates (Newton) and approximate inverses; every
float that enters a proof is converted exactly to a rational first.
"""
import json
import math
import sys
import time
from fractions import Fraction

import numpy as np
from flint import arb, ctx, fmpq

ctx.prec = 160


class ProofFailure(Exception):
    pass


def frac_of(s):
    f = Fraction(s)
    return fmpq(f.numerator, f.denominator)


def fq_from_float(x):
    f = Fraction(float(x))  # exact: a double is a dyadic rational
    return fmpq(f.numerator, f.denominator)


def fq_to_frac(q):
    return Fraction(int(q.p), int(q.q))


# ----------------------------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------------------------
class Config:
    def __init__(self, sites, charges):
        self.n = len(sites)
        self.sites = [[frac_of(c) for c in s] for s in sites]
        self.q = [frac_of(c) for c in charges]
        for qi in self.q:
            if not qi > 0:
                raise ProofFailure("all charges must be positive")
        self.sites_arb = [[arb(c) for c in s] for s in self.sites]
        self.q_arb = [arb(c) for c in self.q]
        self.sites_f = np.array([[float(fq_to_frac(c)) for c in s] for s in self.sites])
        self.q_f = np.array([float(fq_to_frac(c)) for c in self.q])

    # float field and Hessian (candidate search only)
    def gH_float(self, x):
        d = x[None, :] - self.sites_f
        r2 = np.einsum("ij,ij->i", d, d)
        r = np.sqrt(r2)
        w3 = self.q_f / (r2 * r)
        g = -(w3[:, None] * d).sum(0)
        w5 = 3 * self.q_f / (r2 * r2 * r)
        H = np.einsum("i,ij,ik->jk", w5, d, d) - w3.sum() * np.eye(3)
        return g, H

    # rigorous field and Hessian at a point given as three arbs
    def gH_arb(self, c):
        cx, cy, cz = c
        gx = gy = gz = arb(0)
        hxx = hyy = hzz = hxy = hxz = hyz = arb(0)
        for (ax, ay, az), q in zip(self.sites_arb, self.q_arb):
            dx = cx - ax
            dy = cy - ay
            dz = cz - az
            r2 = dx * dx + dy * dy + dz * dz
            ir2 = 1 / r2
            w3 = q * r2.rsqrt() * ir2
            gx -= w3 * dx
            gy -= w3 * dy
            gz -= w3 * dz
            w5 = 3 * w3 * ir2
            hxx += w5 * dx * dx - w3
            hyy += w5 * dy * dy - w3
            hzz += w5 * dz * dz - w3
            hxy += w5 * dx * dy
            hxz += w5 * dx * dz
            hyz += w5 * dy * dz
        g = [gx, gy, gz]
        H = [[hxx, hxy, hxz], [hxy, hyy, hyz], [hxz, hyz, hzz]]
        return g, H


# ----------------------------------------------------------------------------------------------
# L1: half-spaces containing all sites (exact rationals)
# ----------------------------------------------------------------------------------------------
def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def sub(a, b):
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def halfspaces(cfg):
    S = [[fq_to_frac(c) for c in s] for s in cfg.sites]
    n = len(S)
    cand = []
    normals = []
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                N = cross(sub(S[j], S[i]), sub(S[k], S[i]))
                if N != [0, 0, 0]:
                    normals.append(N)
                    # in-plane edge normals (needed when all sites are coplanar)
                    for (p, r) in ((i, j), (j, k), (k, i)):
                        normals.append(cross(sub(S[r], S[p]), N))
    for N in normals:
        for u in (N, [-N[0], -N[1], -N[2]]):
            b = max(dot(u, s) for s in S)
            # checked: every site satisfies u.s <= b (b is the max, exact)
            if all(dot(u, s) <= b for s in S):
                cand.append((u, b))
    # deduplicate
    out = []
    seen = set()
    for u, b in cand:
        g = max(abs(x) for x in u)
        key = tuple(x / g for x in u) + (b / g,)
        if key not in seen:
            seen.add(key)
            out.append((u, b))
    return out


# ----------------------------------------------------------------------------------------------
# L2: exclusion ball around each charge (exact rationals)
# ----------------------------------------------------------------------------------------------
def sqrt_lower(x, bits=80):
    """Rational lower bound for sqrt(x), x a non-negative Fraction."""
    scale = 1 << (2 * bits)
    v = (x.numerator * scale) // x.denominator
    return Fraction(math.isqrt(v), 1 << bits)


def charge_ball_ok(cfg, i, rho):
    S = [[fq_to_frac(c) for c in s] for s in cfg.sites]
    Q = [fq_to_frac(c) for c in cfg.q]
    rhs = Fraction(0)
    for j in range(cfg.n):
        if j == i:
            continue
        d2 = sum((S[i][k] - S[j][k]) ** 2 for k in range(3))
        dl = sqrt_lower(d2)
        if not dl - rho > 0:
            return False
        rhs += Q[j] / (dl - rho) ** 2
    return Q[i] / rho ** 2 > rhs


def charge_balls(cfg):
    balls = []
    for i in range(cfg.n):
        a = cfg.sites_f[i]
        others = [j for j in range(cfg.n) if j != i]
        dmin = min(np.linalg.norm(a - cfg.sites_f[j]) for j in others)

        def f(r):
            return cfg.q_f[i] / r ** 2 - sum(
                cfg.q_f[j] / (np.linalg.norm(a - cfg.sites_f[j]) - r) ** 2 for j in others)

        lo, hi = 0.0, dmin
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(mid) > 0:
                lo = mid
            else:
                hi = mid
        # dyadic radius a little below the float root
        e = math.floor(math.log2(lo * 0.9))
        rho = Fraction(math.floor(lo * 0.9 / 2.0 ** (e - 20)), 1) * Fraction(2) ** (e - 20)
        if not charge_ball_ok(cfg, i, rho):
            raise ProofFailure(f"charge ball {i} failed exact check")
        balls.append(rho)
    return balls


# ----------------------------------------------------------------------------------------------
# Z: zero balls
# ----------------------------------------------------------------------------------------------
def newton(cfg, x, iters=80):
    """Candidate search only. x: float start. The iterate is kept as exact rationals and g, H are
    evaluated with arb at ctx.prec bits (midpoints used), so that clusters far below float64
    resolution of the coordinates can be resolved. Returns a list of 3 fmpq or None."""
    xq = [fq_from_float(v) for v in x]
    last = None
    for _ in range(iters):
        xa = [arb(v) for v in xq]
        g, H = cfg.gH_arb(xa)
        gm = np.array([float(v.mid()) for v in g])
        Hm = np.array([[float(H[i][j].mid()) for j in range(3)] for i in range(3)])
        if not (np.all(np.isfinite(gm)) and np.all(np.isfinite(Hm))):
            return None
        try:
            step = np.linalg.solve(Hm, gm)
        except np.linalg.LinAlgError:
            return None
        nrm = float(np.linalg.norm(step))
        if not np.isfinite(nrm):
            return None
        xf = np.array([float(fq_to_frac(v)) for v in xq])
        dmin = float(np.min(np.linalg.norm(cfg.sites_f - xf, axis=1)))
        if dmin == 0:
            return None
        if nrm > 0.5 * dmin:
            step *= 0.5 * dmin / nrm
        xq = [xq[k] - fq_from_float(step[k]) for k in range(3)]
        last = nrm
        if nrm < 1e-34 * max(1.0, float(np.linalg.norm(xf))):
            break
    if last is None or not last < 1e-25 * dmin:
        return None
    # keep the rationals short: round to 2^-120 (far below every scale used here)
    sc = 1 << 120
    xq = [fmpq(int((v * sc).floor()), sc) for v in xq]
    return xq


def frob_upper(M):
    s = arb(0)
    for row in M:
        for v in row:
            s += v * v
    return s.sqrt()


def verify_zero(cfg, x0, rho_cap=None):
    """x0: list of 3 fmpq. Return (x0, rho exact, index, det) or None."""
    xa = [arb(v) for v in x0]
    g, H = cfg.gH_arb(xa)
    Hm = np.array([[float(H[i][j].mid()) for j in range(3)] for i in range(3)])
    try:
        Af = np.linalg.inv(Hm)
    except np.linalg.LinAlgError:
        return None
    A = [[arb(fq_from_float(Af[i][j])) for j in range(3)] for i in range(3)]
    E = [[(arb(1) if i == j else arb(0)) - sum((A[i][k] * H[k][j] for k in range(3)), arb(0))
          for j in range(3)] for i in range(3)]
    nE = frob_upper(E)
    nA = frob_upper(A)
    y = [sum((A[i][k] * g[k] for k in range(3)), arb(0)) for i in range(3)]
    ny = (y[0] * y[0] + y[1] * y[1] + y[2] * y[2]).sqrt()
    d = []
    for s in cfg.sites_arb:
        dd = sub(xa, s)
        d.append(dot(dd, dd).sqrt())
    dmin = min(float(v.lower()) for v in d)
    rho_f = 0.5 * dmin if rho_cap is None else min(0.5 * dmin, float(rho_cap) * 0.999)
    floor_f = max(float(ny.upper()), 1e-300)
    best = None
    while rho_f > 2 * floor_f:
        e = math.floor(math.log2(rho_f))
        rho = fmpq(math.floor(rho_f / 2.0 ** (e - 30)), 1) * fmpq(2) ** (e - 30) if e - 30 >= 0 else \
            fmpq(math.floor(rho_f * 2.0 ** (30 - e)), 2 ** (30 - e))
        ra = arb(rho)
        ok = True
        L3 = arb(0)
        for di, q in zip(d, cfg.q_arb):
            gap = di - ra
            if not gap > 0:
                ok = False
                break
            L3 += 6 * q / (gap * gap * gap * gap)
        if ok:
            kappa = nE + nA * ra * L3
            if kappa < 1 and ny <= (1 - kappa) * ra:
                if rho_cap is None or rho < rho_cap:
                    best = rho
                    break
        rho_f *= 0.5
    if best is None:
        return None
    det = (H[0][0] * (H[1][1] * H[2][2] - H[1][2] * H[2][1])
           - H[0][1] * (H[1][0] * H[2][2] - H[1][2] * H[2][0])
           + H[0][2] * (H[1][0] * H[2][1] - H[1][1] * H[2][0]))
    tr = H[0][0] + H[1][1] + H[2][2]
    if not tr.contains(0):
        raise ProofFailure("trace of Hessian does not contain 0: evaluation bug")
    if det > 0:
        index = 2
    elif det < 0:
        index = 1
    else:
        return None
    return x0, best, index, det


class ZeroSet:
    def __init__(self, cfg):
        self.cfg = cfg
        self.balls = []  # (x0 fmpq[3], rho fmpq, index, x0 float, rho float)

    def contains_point_f(self, x):
        for b in self.balls:
            if np.linalg.norm(b[3] - x) <= b[4]:
                return True
        return False

    def try_add(self, xq):
        """Verify a new zero ball around xq (3 fmpq). Return True if a new ball was added."""
        for b in self.balls:
            d2 = sum(((xq[k] - b[0][k]) ** 2 for k in range(3)), fmpq(0))
            if d2 <= b[1] ** 2:
                return False
        cap = None
        for _ in range(4):
            res = verify_zero(self.cfg, xq, cap)
            if res is None:
                return False
            x0, rho, index, det = res
            clash = None
            for b in self.balls:
                d2 = sum(((x0[k] - b[0][k]) ** 2 for k in range(3)), fmpq(0))
                if d2 <= b[1] ** 2:
                    return False  # x0 lies in an existing ball: same zero
                if not d2 > (rho + b[1]) ** 2:
                    clash = b
                    break
            if clash is None:
                xf0 = np.array([float(fq_to_frac(v)) for v in x0])
                self.balls.append((x0, rho, index, xf0, float(fq_to_frac(rho))))
                return True
            dist = math.sqrt(float(fq_to_frac(sum(((x0[k] - clash[0][k]) ** 2 for k in range(3)), fmpq(0)))))
            cap_f = 0.999 * (dist - clash[4])
            if cap_f <= 0:
                return False
            cap = fq_from_float(cap_f)
        return False

    def check_pairwise_disjoint(self):
        B = self.balls
        for i in range(len(B)):
            for j in range(i + 1, len(B)):
                d2 = sum(((B[i][0][k] - B[j][0][k]) ** 2 for k in range(3)), fmpq(0))
                if not d2 > (B[i][1] + B[j][1]) ** 2:
                    raise ProofFailure(f"zero balls {i} and {j} overlap")


# ----------------------------------------------------------------------------------------------
# Octree subdivision
# ----------------------------------------------------------------------------------------------
def run(cfg, max_level=64, newton_level=4, log=None, box_budget=20_000_000):
    t0 = time.time()
    S = cfg.sites
    lo = [min(s[k] for s in S) for k in range(3)]
    hi = [max(s[k] for s in S) for k in range(3)]
    span = [hi[k] - lo[k] for k in range(3)]
    # degenerate span (all sites share a coordinate) cannot happen for these configs; guard anyway
    for k in range(3):
        if not span[k] > 0:
            raise ProofFailure("sites lie in a coordinate plane; root box degenerate")
    HS = [([fmpq(x.numerator, x.denominator) for x in u], fmpq(b.numerator, b.denominator))
          for u, b in halfspaces(cfg)]
    rho_ch = charge_balls(cfg)
    rho_ch_q = [fmpq(r.numerator, r.denominator) for r in rho_ch]
    zs = ZeroSet(cfg)
    stats = dict(boxes=0, L1=0, L2=0, T1=0, T2=0, Z=0, newton_tries=0, max_level=0)
    lo_f = np.array([float(fq_to_frac(v)) for v in lo])
    span_f = np.array([float(fq_to_frac(v)) for v in span])

    stack = [(0, 0, 0, 0)]
    while stack:
        L, i, j, k = stack.pop()
        stats["boxes"] += 1
        if stats["boxes"] > box_budget:
            raise ProofFailure("box budget exceeded")
        if L > stats["max_level"]:
            stats["max_level"] = L
        if log and stats["boxes"] % 200000 == 0:
            print(f"  boxes={stats['boxes']} stack={len(stack)} zeros={len(zs.balls)} "
                  f"t={time.time()-t0:.0f}s {stats}", file=log, flush=True)
        den = fmpq(2) ** L
        idx = (i, j, k)
        bl = [lo[a] + span[a] * idx[a] / den for a in range(3)]
        bh = [lo[a] + span[a] * (idx[a] + 1) / den for a in range(3)]
        h = [span[a] / (2 * den) for a in range(3)]
        c = [bl[a] + h[a] for a in range(3)]
        # L1
        out = False
        for u, b in HS:
            m = u[0] * c[0] + u[1] * c[1] + u[2] * c[2] - (abs(u[0]) * h[0] + abs(u[1]) * h[1] + abs(u[2]) * h[2])
            if m > b:
                out = True
                break
        if out:
            stats["L1"] += 1
            continue
        # distances box -> sites
        delta2 = []
        far2 = []
        for s in S:
            dsum = fmpq(0)
            fsum = fmpq(0)
            for a in range(3):
                t1 = bl[a] - s[a]
                t2 = s[a] - bh[a]
                if t1 > 0:
                    dsum += t1 * t1
                elif t2 > 0:
                    dsum += t2 * t2
                fsum += max(t1 * t1, t2 * t2)
            delta2.append(dsum)
            far2.append(fsum)
        # L2
        done = False
        for si in range(cfg.n):
            if far2[si] <= rho_ch_q[si] ** 2:
                done = True
                break
        if done:
            stats["L2"] += 1
            continue
        # Z containment
        cf = lo_f + span_f * (np.array(idx) + 0.5) / 2.0 ** L
        Rf = 0.5 * np.linalg.norm(span_f) / 2.0 ** L
        for b in zs.balls:
            if np.linalg.norm(b[3] - cf) + Rf < b[4] * 1.0000001 + 1e-300:
                corner_far = sum((max((bl[a] - b[0][a]) ** 2, (bh[a] - b[0][a]) ** 2) for a in range(3)), fmpq(0))
                if corner_far <= b[1] ** 2:
                    done = True
                    break
        if done:
            stats["Z"] += 1
            continue
        if all(d2 > 0 for d2 in delta2):
            ca = [arb(v) for v in c]
            g, H = cfg.gH_arb(ca)
            R = arb(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]).sqrt()
            M2 = arb(0)
            M3 = arb(0)
            for d2, q in zip(delta2, cfg.q_arb):
                da = arb(d2)
                dl = da.sqrt()
                M2 += 2 * q / (da * dl)
                M3 += 6 * q / (da * da)
            gn = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt()
            if gn > R * M2:
                stats["T1"] += 1
                continue
            rem = R * R * M3 / 2
            ok = False
            for a in range(3):
                Hr = (H[a][0] * H[a][0] + H[a][1] * H[a][1] + H[a][2] * H[a][2]).sqrt()
                if abs(g[a]) - R * Hr > rem:
                    ok = True
                    break
            if not ok:
                u = [arb(fq_from_float(float(v.mid()))) for v in g]
                un = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt()
                if un > 0:
                    ug = u[0] * g[0] + u[1] * g[1] + u[2] * g[2]
                    Hu = [H[a][0] * u[0] + H[a][1] * u[1] + H[a][2] * u[2] for a in range(3)]
                    Hun = (Hu[0] * Hu[0] + Hu[1] * Hu[1] + Hu[2] * Hu[2]).sqrt()
                    if abs(ug) - R * Hun > rem * un:
                        ok = True
            if ok:
                stats["T2"] += 1
                continue
        # undecided: maybe look for a zero, then subdivide
        if L >= newton_level and not zs.contains_point_f(cf):
            stats["newton_tries"] += 1
            x = newton(cfg, cf)
            if x is not None:
                zs.try_add(x)
        if L >= max_level:
            raise ProofFailure(f"undecided box at max level: L={L} idx={idx} centre={cf}")
        for di in (0, 1):
            for dj in (0, 1):
                for dk in (0, 1):
                    stack.append((L + 1, 2 * i + di, 2 * j + dj, 2 * k + dk))
    zs.check_pairwise_disjoint()
    stats["seconds"] = round(time.time() - t0, 1)
    return zs, stats, rho_ch, HS


def main(argv):
    path = argv[1]
    out = argv[2] if len(argv) > 2 else None
    D = json.load(open(path))
    cfg = Config(D["sites"], D["charges"])
    print(f"{path}: n={cfg.n}", flush=True)
    zs, stats, rho_ch, HS = run(cfg, log=sys.stdout)
    m1 = sum(1 for b in zs.balls if b[2] == 1)
    m2 = sum(1 for b in zs.balls if b[2] == 2)
    count = len(zs.balls)
    res = dict(file=path, n=cfg.n, count=count, index1=m1, index2=m2,
               morse_m2_minus_m1=m2 - m1, n_minus_1=cfg.n - 1, stats=stats,
               charge_ball_radii=[str(r) for r in rho_ch], n_halfspaces=len(HS),
               zeros=[dict(centre=[str(v) for v in b[0]], centre_float=b[3].tolist(),
                           rho=str(b[1]), rho_float=b[4], index=b[2]) for b in zs.balls])
    print(f"RESULT {path}: count={count} index1={m1} index2={m2} m2-m1={m2-m1} n-1={cfg.n-1} {stats}",
          flush=True)
    if out:
        json.dump(res, open(out, "w"), indent=1)
    claimed = D.get("count")
    if claimed is not None:
        try:
            claimed = int(claimed)          # a count stored as text is compared as a number
        except (TypeError, ValueError):
            print(f"MISMATCH: the file's count field is not a whole number: {claimed!r}", flush=True)
            return 2
        if claimed != count:
            print(f"MISMATCH: file claims {claimed}, independent count is {count}", flush=True)
            return 2
        print(f"AGREES with claimed count {claimed}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv))
    except ProofFailure as e:
        print(f"PROOF FAILURE: {e}", flush=True)
        sys.exit(3)
