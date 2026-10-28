import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
# Randomised soundness test of ia.py against exact rational arithmetic (fractions).
import numpy as np
from fractions import Fraction as F
import ia
rng = np.random.default_rng(1)
N = 20000
def rnd():
    lo = rng.standard_normal(N) * 10.0 ** rng.integers(-5, 5, N)
    w = np.abs(rng.standard_normal(N)) * 10.0 ** rng.integers(-12, 1, N)
    return (lo, lo + w)
a, b = rnd(), rnd()
ops = {'add': (ia.add, lambda x, y: x + y), 'sub': (ia.sub, lambda x, y: x - y),
       'mul': (ia.mul, lambda x, y: x * y), 'div': (ia.div, lambda x, y: x / y)}
bad = 0
for name, (fi, fe) in ops.items():
    r = fi(a, b)
    for k in range(0, N, 97):
        for ta in (a[0][k], a[1][k], (a[0][k] + a[1][k]) / 2):
            for tb in (b[0][k], b[1][k]):
                if name == 'div' and (b[0][k] <= 0 <= b[1][k]):
                    continue
                if name == 'div' and tb == 0:
                    continue
                e = fe(F(ta), F(tb))
                if not (F(r[0][k]) <= e <= F(r[1][k])):
                    bad += 1; print(name, ta, tb, r[0][k], r[1][k])
s = ia.sqr(a); sq = ia.sqrt((np.abs(a[0]), np.abs(a[0]) + (a[1] - a[0])))
for k in range(0, N, 97):
    for t in (a[0][k], a[1][k]):
        if not (F(s[0][k]) <= F(t) ** 2 <= F(s[1][k])): bad += 1
    t = abs(a[0][k])
    if not (F(sq[0][k]) ** 2 <= F(t) and F(t) <= F(sq[1][k]) ** 2): bad += 1
print('violations:', bad)
sys.exit(1 if bad else 0)
