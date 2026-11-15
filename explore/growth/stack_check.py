#!/usr/bin/env python3
"""Exact count of proved real critical points of a stack of zonal harmonics (step 3 of the growth note).
Written apart from the search programs hier.py and harm.py.

P = sum_{k=1..D} c_k N_k Z_k, where Z_k = r^k P_k(axis/r) is the zonal harmonic of degree k about z when D - k is even
and about x when D - k is odd, c_k = s^((D-k)(D-k+1)/2) with s = 1/100, and N_k is a normalisation:
  'legendre' : N_k = 1 (so Z_k = 1 at the pole);  'maxcoef' : N_k = 1 / (largest absolute monomial coefficient of Z_k).
All coefficients are exact rationals. Candidates come from Newton in 70-digit arithmetic with the exact gradient and
Hessian, from starts spread over every scale. Each candidate is then proved by a Krawczyk test in ball arithmetic,
which gives exactly one non-degenerate critical point in a small cube. Cubes are checked pairwise disjoint.
The output is the number proved, grouped by distance from the origin.
Usage: stack_check.py <D> <legendre|maxcoef> <starts per scale>"""
import sys, itertools
from fractions import Fraction
from math import comb, factorial
import numpy as np, mpmath as mp
from flint import arb, ctx, fmpq
ctx.prec = 400; mp.mp.dps = 70
D = int(sys.argv[1]); NORM = sys.argv[2]; NST = int(sys.argv[3]); s = Fraction(1, 100)
def legendre_coeffs(k):
    """P_k(u) = sum_j a_j u^(k-2j), exact."""
    return {k - 2 * j: Fraction((-1) ** j * comb(k, j) * comb(2 * k - 2 * j, k), 2 ** k) for j in range(k // 2 + 1)}
def zonal(k, axis):
    C = {}
    for e, a in legendre_coeffs(k).items():          # a * u^e * r^(k-e), with r^(k-e) = (x^2+y^2+z^2)^((k-e)/2)
        h = (k - e) // 2
        for i in range(h + 1):
            for j in range(h + 1 - i):
                l = h - i - j; m = [2 * i, 2 * j, 2 * l]; m[axis] += e
                C[tuple(m)] = C.get(tuple(m), Fraction(0)) + a * Fraction(factorial(h), factorial(i) * factorial(j) * factorial(l))
    return C
P = {}
for k in range(1, D + 1):
    Z = zonal(k, 2 if (D - k) % 2 == 0 else 0)
    nrm = Fraction(1) if NORM == "legendre" else 1 / max(abs(v) for v in Z.values())
    ck = s ** ((D - k) * (D - k + 1) // 2)
    for m, v in Z.items(): P[m] = P.get(m, Fraction(0)) + ck * nrm * v
lap = {}
for (i, j, k), c in P.items():
    for v, e in enumerate((i, j, k)):
        if e >= 2:
            t = [i, j, k]; t[v] -= 2; lap[tuple(t)] = lap.get(tuple(t), Fraction(0)) + c * e * (e - 1)
assert_ok = all(v == 0 for v in lap.values())
if not assert_ok: raise SystemExit("the polynomial is not harmonic")
terms = [(m, c) for m, c in P.items() if c != 0]
def grad_hess(X, num):
    g = [num(0)] * 3; H = [[num(0)] * 3 for _ in range(3)]
    pw = [[num(1)] for _ in range(3)]
    for v in range(3):
        for e in range(1, D + 1): pw[v].append(pw[v][-1] * X[v])
    for (e, c) in terms:
        cc = num(c)
        for u in range(3):
            if e[u] == 0: continue
            eu = list(e); eu[u] -= 1
            g[u] = g[u] + cc * e[u] * pw[0][eu[0]] * pw[1][eu[1]] * pw[2][eu[2]]
            for v in range(3):
                if eu[v] == 0: continue
                ev = list(eu); ev[v] -= 1
                H[u][v] = H[u][v] + cc * e[u] * eu[v] * pw[0][ev[0]] * pw[1][ev[1]] * pw[2][ev[2]]
    return g, H
mpn = lambda c: mp.mpf(c.numerator) / mp.mpf(c.denominator) if isinstance(c, Fraction) else mp.mpf(c)
arn = lambda c: arb(fmpq(c.numerator, c.denominator)) if isinstance(c, Fraction) else arb(c)
rng = np.random.default_rng(4); found = []
for level in range(0, D):
    rad = float(s) ** level
    for trial in range(NST):
        v = rng.normal(size=3); v *= rad * 10 ** rng.uniform(-1.3, 1.0) / np.linalg.norm(v)
        x = mp.matrix([mp.mpf(float(t)) for t in v]); ok = False
        for it in range(120):
            g, H = grad_hess(x, mpn)
            try: dx = mp.lu_solve(mp.matrix(H), mp.matrix(g))
            except ZeroDivisionError: break
            if mp.norm(dx) > 3 * mp.norm(x): dx = dx * (3 * mp.norm(x) / mp.norm(dx))     # stay within the scale
            x = x - dx
            if mp.norm(dx) < mp.mpf(10) ** -50 * mp.norm(x): ok = True; break
        if not ok or mp.norm(x) == 0: continue
        if all(mp.norm(x - y) > mp.mpf(10) ** -15 * (mp.norm(x) + mp.norm(y)) for y in found): found.append(x)
print(f"D = {D}, normalisation {NORM}: candidates from {NST * D} starts: {len(found)}", flush=True)
proved = []
for x in found:
    sep = min([mp.norm(x - y) for y in found if y is not x])
    centre = [Fraction(int(mp.nint(x[k] * mp.mpf(10) ** 60)), 10 ** 60) for k in range(3)]
    cA = [arn(c) for c in centre]
    g, H = grad_hess(x, mpn); Ainv = mp.matrix(H) ** -1
    Aq = [[arn(Fraction(int(mp.nint(Ainv[i, j] * mp.mpf(10) ** 20)), 10 ** 20)) for j in range(3)] for i in range(3)]
    r = Fraction(int(mp.nint(sep / 4 * mp.mpf(10) ** 60)), 10 ** 60); done = False
    for attempt in range(50):
        r = r / 2; rb = arn(r); box = [cA[k] + arb(0, 1) * rb for k in range(3)]
        gc, _ = grad_hess(cA, arn); _, HB = grad_hess(box, arn)
        K = []
        for i in range(3):
            v = -sum(Aq[i][j] * gc[j] for j in range(3))
            for j in range(3): v = v + ((1 if i == j else 0) - sum(Aq[i][k] * HB[k][j] for k in range(3))) * (arb(0, 1) * rb)
            K.append(v)
        det = (HB[0][0] * (HB[1][1] * HB[2][2] - HB[1][2] * HB[2][1]) - HB[0][1] * (HB[1][0] * HB[2][2] - HB[1][2] * HB[2][0]) + HB[0][2] * (HB[1][0] * HB[2][1] - HB[1][1] * HB[2][0]))
        if all(K[i].abs_upper() < rb * arb("0.999") for i in range(3)) and (det > 0 or det < 0):
            proved.append((centre, r, 2 if det > 0 else 1)); done = True; break
    if not done: print("  not proved: candidate at distance %.3e" % float(mp.norm(x)))
for i in range(len(proved)):
    for j in range(i):
        ci, ri, _ = proved[i]; cj, rj, _ = proved[j]
        if not any(abs(ci[k] - cj[k]) > ri + rj for k in range(3)): raise SystemExit("two cubes overlap")
rad = np.array(sorted(float(mp.sqrt(sum(mpn(c) ** 2 for c in p[0]))) for p in proved))
groups = [int(((rad > float(s) ** (j + 0.5)) & (rad <= float(s) ** (j - 0.5))).sum()) for j in range(0, D + 1)]
I = lambda p: p * p // 2 + (2 if p % 4 == 3 else 0)
print(f"PROVED: {len(proved)} real non-degenerate critical points (det Hessian positive at {sum(1 for p in proved if p[2] == 2)}, negative at {sum(1 for p in proved if p[2] == 1)}); the stacking rule predicts {sum(I(p) - 1 for p in range(2, D + 1))}; at most {(D - 1) ** 3} exist")
print("  by distance from the origin, outermost first:", groups, " rule:", [0] + [I(p) - 1 for p in range(D, 1, -1)])
if len(proved) < sum(I(p) - 1 for p in range(2, D + 1)): raise SystemExit("FAIL: fewer points proved than the rule predicts (try more starts)")
