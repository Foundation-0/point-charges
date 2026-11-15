#!/usr/bin/env python3
"""Lower-bound proof for the number of balance points of positive point charges. Written separately,
without reading indep_count.py or verify_arb.py, as a further check.

Input: JSON with "sites" (lists of three fraction strings) and "charges" (fraction strings). Nothing else is read.
Method: candidates are found here by Newton from random starts in floating point. Each candidate x0 is then tested by
the Krawczyk operator in ball arithmetic (python-flint arb): with B the cube of half-width r about x0 and
A an approximate inverse of the Hessian at x0,
        K = x0 - A g(x0) + (I - A H(B)) (B - x0).
If K lies strictly inside B, the gradient g has exactly one zero in B. If det H(B) has one sign on B, the zero is
non-degenerate, with Morse index 2 for a positive sign and 1 for a negative one (the Hessian is traceless).
Cubes are checked pairwise disjoint in exact rationals. The program proves 'at least N zeros'. It does not prove
that there are no others. Every failed proof step raises an error.
Usage: python3 -I -O lower_bound_check.py <config.json> <expected count>
Exit status 0 only if at least the expected number of zeros is proved. Needs python-flint and NumPy."""
import json, os, sys
from fractions import Fraction
import numpy as np
from flint import arb, ctx, fmpq
ctx.prec = 200
class CheckFailure(Exception): pass
def _hook(kind, value, tb):
    # print the failure without a traceback and stop with a nonzero status
    print(f"{kind.__name__}: {value}", file=sys.stderr, flush=True); sys.stdout.flush(); os._exit(1)
sys.excepthook = _hook
def Fq(s): f = Fraction(s); return fmpq(f.numerator, f.denominator)
cfg = json.load(open(sys.argv[1])); expected = int(sys.argv[2])
sites_f = [[Fraction(c) for c in s] for s in cfg["sites"]]; charges_f = [Fraction(c) for c in cfg["charges"]]
if len(sites_f) != len(charges_f): raise CheckFailure("sites and charges differ in number")
if any(q <= 0 for q in charges_f): raise CheckFailure("a charge is not positive")
n = len(sites_f)
A_np = np.array([[float(c) for c in s] for s in sites_f]); q_np = np.array([float(q) for q in charges_f])
def gH(X):
    D = X[:, None, :] - A_np[None]; r2 = (D * D).sum(-1); r = np.sqrt(r2)
    g = -(q_np[None, :, None] * D / (r2 * r)[..., None]).sum(1)
    H = (q_np[None, :, None, None] * (3 * D[..., :, None] * D[..., None, :] - r2[..., None, None] * np.eye(3)) / (r2 * r2 * r)[..., None, None]).sum(1)
    return g, H
rng = np.random.default_rng(20261115)
lo = A_np.min(0); hi = A_np.max(0); R = np.linalg.norm(hi - lo)
cands = []
for rnd in range(6):
    X = rng.uniform(lo, hi, size=(20000, 3))
    if rnd >= 2:                                    # also start densely around the centroid and the charge-weighted centre
        c = (A_np.mean(0), (q_np[:, None] * A_np).sum(0) / q_np.sum())[rnd % 2]
        X = c + rng.normal(size=(20000, 3)) * R * 10 ** rng.uniform(-4, -0.5, size=(20000, 1))
    with np.errstate(all="ignore"):
        for it in range(100):
            g, H = gH(X)
            d = np.linalg.solve(H, g[..., None])[..., 0]
            s = np.linalg.norm(d, axis=1); d *= np.minimum(1.0, 0.03 * R / np.maximum(s, 1e-300))[:, None]
            X = X - d
        g, H = gH(X)
    ok = np.isfinite(X).all(1) & (np.linalg.norm(g, axis=1) < 1e-10 * q_np.sum() / R**2) & (X >= lo - 1e-9).all(1) & (X <= hi + 1e-9).all(1)
    for x in X[ok]:
        if all(np.linalg.norm(x - c) > 1e-8 * R for c in cands): cands.append(x)
print("candidates found by float Newton:", len(cands))
sites = [[arb(Fq(str(c))) for c in s] for s in sites_f]; charges = [arb(Fq(str(q))) for q in charges_f]
def g_arb(x):
    g = [arb(0)] * 3
    for a, q in zip(sites, charges):
        d = [x[k] - a[k] for k in range(3)]; r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]; r3 = r2 * r2.sqrt()
        g = [g[k] - q * d[k] / r3 for k in range(3)]
    return g
def H_arb(x):
    H = [[arb(0)] * 3 for _ in range(3)]
    for a, q in zip(sites, charges):
        d = [x[k] - a[k] for k in range(3)]; r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
        if not (r2 > 0): raise CheckFailure("cube touches a charge")
        r5 = r2 * r2 * r2.sqrt()
        for i in range(3):
            for j in range(3):
                H[i][j] = H[i][j] + q * (3 * d[i] * d[j] - (r2 if i == j else 0)) / r5
    return H
def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
proved = []
for x0 in cands:
    # refine in floats, then fix an exact rational centre
    x = x0.copy()
    for it in range(5):
        g, H = gH(x[None]); x = x - np.linalg.solve(H[0], g[0])
    g, H = gH(x[None]); Ainv = np.linalg.inv(H[0])
    centre = [Fraction(float(v)) for v in x]
    done = False
    for rexp in range(8, 40):
        r = Fraction(1, 2 ** rexp)
        c_arb = [arb(Fq(str(v))) for v in centre]
        box = [arb(Fq(str(v)), float(r)) for v in centre]          # ball of radius r about each coordinate
        try:
            gc = g_arb(c_arb); HB = H_arb(box)
        except CheckFailure:
            continue
        Aq = [[arb(Fq(str(Fraction(float(Ainv[i][j]))))) for j in range(3)] for i in range(3)]
        delta = [arb(0, float(r)) for _ in range(3)]
        K = []
        for i in range(3):
            v = c_arb[i] - sum(Aq[i][j] * gc[j] for j in range(3))
            for j in range(3):
                m = (1 if i == j else 0) - sum(Aq[i][k] * HB[k][j] for k in range(3))
                v = v + m * delta[j]
            K.append(v)
        inside = all((K[i] - c_arb[i]).abs_upper() < arb(Fq(str(r))) * arb("0.999") for i in range(3))
        dH = det3(HB)
        if inside and (dH > 0 or dH < 0):
            proved.append((centre, r, 2 if dH > 0 else 1)); done = True; break
    if not done:
        print("  candidate not proved:", x)
for i in range(len(proved)):
    for j in range(i):
        ci, ri, _ = proved[i]; cj, rj, _ = proved[j]
        if not any(abs(ci[k] - cj[k]) > ri + rj for k in range(3)): raise CheckFailure("two cubes overlap")
m2 = sum(1 for p in proved if p[2] == 2); m1 = sum(1 for p in proved if p[2] == 1)
print(f"PROVED: at least {len(proved)} non-degenerate balance points; index 2: {m2}, index 1: {m1}; index-2 minus index-1 = {m2 - m1}, charges minus one = {n - 1}")
if len(proved) < expected: raise CheckFailure(f"proved {len(proved)}, expected at least {expected}")
if m2 - m1 != n - 1: print("note: the proved set does not satisfy the Morse relation, so it is not the full set")
print("OK")
