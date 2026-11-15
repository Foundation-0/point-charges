#!/usr/bin/env python3
"""Lower-bound proof for the number of balance points of positive point charges, for configurations whose balance
points sit at widely different scales. Written separately from the other checkers. It takes a list of candidate points from any source
and trusts none of them: each is refined by Newton in 260-digit arithmetic and then proved by the Krawczyk test in
ball arithmetic (python-flint arb, 1400 bits). Near a point where the potential is flat to high order the Hessian is many orders smaller than its separate terms, so the cubes there are made very small. A proved cube holds exactly one zero of the gradient, and the sign of
det H on the cube gives the index (the Hessian is traceless). Cubes are checked pairwise disjoint in exact rationals.
The result is 'at least N balance points'. It does not prove that there are no others.
Usage: python3 -I -O lower_bound_hp.py <config.json> <candidates .npy or .json> <expected count>
Exit status 0 only if at least the expected number of balance points is proved. Needs NumPy, mpmath and python-flint."""
import json, os, sys
from fractions import Fraction
import numpy as np, mpmath as mp
from flint import arb, ctx, fmpq
ctx.prec = 1400; mp.mp.dps = 260
class CheckFailure(Exception): pass
def _hook(kind, value, tb):
    print(f"{kind.__name__}: {value}", file=sys.stderr, flush=True); sys.stdout.flush(); os._exit(1)
sys.excepthook = _hook
cfg = json.load(open(sys.argv[1])); expected = int(sys.argv[3])
sites_f = [[Fraction(c) for c in s] for s in cfg["sites"]]; charges_f = [Fraction(c) for c in cfg["charges"]]
if any(q <= 0 for q in charges_f): raise CheckFailure("a charge is not positive")
n = len(sites_f)
cand = np.load(sys.argv[2]) if sys.argv[2].endswith(".npy") else np.array(json.load(open(sys.argv[2])), float)
cand = np.asarray(cand, float); cand = cand[:, :3] if cand.ndim == 2 else cand.reshape(-1, 3)   # any further columns (an index label) are ignored
mpF = lambda f: mp.mpf(f.numerator) / mp.mpf(f.denominator)
A_mp = [[mpF(c) for c in s] for s in sites_f]; q_mp = [mpF(q) for q in charges_f]
def gH_mp(x):
    g = [mp.mpf(0)] * 3; H = [[mp.mpf(0)] * 3 for _ in range(3)]
    for a, q in zip(A_mp, q_mp):
        d = [x[k] - a[k] for k in range(3)]; r2 = d[0] ** 2 + d[1] ** 2 + d[2] ** 2; r = mp.sqrt(r2)
        for i in range(3):
            g[i] -= q * d[i] / (r2 * r)
            for j in range(3): H[i][j] += q * (3 * d[i] * d[j] - (r2 if i == j else 0)) / (r2 * r2 * r)
    return mp.matrix(g), mp.matrix(H)
Fq = lambda f: fmpq(f.numerator, f.denominator)
sites = [[arb(Fq(c)) for c in s] for s in sites_f]; charges = [arb(Fq(q)) for q in charges_f]
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
            for j in range(3): H[i][j] = H[i][j] + q * (3 * d[i] * d[j] - (r2 if i == j else 0)) / r5
    return H
def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
def to_frac(v, digits=200): return Fraction(int(mp.nint(v * mp.mpf(10) ** digits)), 10 ** digits)
refined = []
for x0 in cand:
    x = mp.matrix([mp.mpf(float(v)) for v in x0]); good = True
    for it in range(60):
        g, H = gH_mp(x)
        try: dx = mp.lu_solve(H, g)
        except ZeroDivisionError: good = False; break
        x = x - dx
        if mp.norm(dx) < mp.mpf(10) ** (-180) * mp.norm(x): break
    if not good or mp.norm(x - mp.matrix([mp.mpf(float(v)) for v in x0])) > 0.5 * (mp.norm(x) + 1e-30) + 1e-3: continue
    if all(mp.norm(x - y) > mp.mpf(10) ** (-30) * mp.norm(x) for y in refined): refined.append(x)
print("candidates given:", len(cand), " distinct after refinement:", len(refined))
proved = []
for x in refined:
    others = [mp.norm(x - y) for y in refined if y is not x]
    sep = min(others + [min(mp.norm(x - mp.matrix(a)) for a in A_mp)])
    g, H = gH_mp(x); Ainv = H ** -1
    centre = [to_frac(x[k]) for k in range(3)]; c_arb = [arb(Fq(v)) for v in centre]
    Aq = [[arb(Fq(to_frac(Ainv[i, j], 40))) for j in range(3)] for i in range(3)]
    done = False; r = to_frac(sep / 4, 200)
    for attempt in range(70):
        r = r / (2 if attempt < 10 else 1000); rb = arb(Fq(r))     # near a flat point the cube must be tiny: the
        # interval Hessian cannot see the cancellation between charges, so its width must fall below the Hessian itself
        box = [c_arb[k] + arb(0, 1) * rb for k in range(3)]
        try: gc = g_arb(c_arb); HB = H_arb(box)
        except CheckFailure: continue
        K = []
        for i in range(3):
            v = -sum(Aq[i][j] * gc[j] for j in range(3))
            for j in range(3):
                m = (1 if i == j else 0) - sum(Aq[i][k] * HB[k][j] for k in range(3))
                v = v + m * (arb(0, 1) * rb)
            K.append(v)
        inside = all(K[i].abs_upper() < rb * arb("0.999") for i in range(3))
        dH = det3(HB)
        if inside and (dH > 0 or dH < 0):
            proved.append((centre, r, 2 if dH > 0 else 1)); done = True; break
    if not done: print("  not proved: candidate at distance %.3e from the origin" % float(mp.norm(x)))
for i in range(len(proved)):
    for j in range(i):
        ci, ri, _ = proved[i]; cj, rj, _ = proved[j]
        if not any(abs(ci[k] - cj[k]) > ri + rj for k in range(3)): raise CheckFailure("two cubes overlap")
m2 = sum(1 for p in proved if p[2] == 2); m1 = sum(1 for p in proved if p[2] == 1)
rad = sorted(float(mp.sqrt(sum(mpF(c) ** 2 for c in p[0]))) for p in proved)
print("distances of the proved points from the origin:", " ".join("%.2e" % v for v in rad))
print(f"PROVED: at least {len(proved)} non-degenerate balance points; index 2: {m2}, index 1: {m1}; index-2 minus index-1 = {m2 - m1}, charges minus one = {n - 1}")
if len(proved) < expected: raise CheckFailure(f"proved {len(proved)}, expected at least {expected}")
print("OK")
