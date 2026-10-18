"""Importance sampling for four charges through the inverse map.

Pick a random tetrahedron, pick random interior points x where det H < 0 (x is an index-1
equilibrium for its own charge vector zeta(x)), and count all equilibria for zeta(x).
Every sample therefore has at least one index-1 equilibrium. Any Morse-consistent count
above 9 is logged.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np, sys, json, time
from multiprocessing import Pool
from collections import Counter
from four_map import zeta_D
from charges import count_fast

def rand_tet(rng):
    k = rng.integers(4)
    if k == 0:
        A = rng.standard_normal((4, 3))
    elif k == 1:        # flattened
        A = rng.standard_normal((4, 3)); A[:, 2] *= 10 ** rng.uniform(-3, -0.5)
    elif k == 2:        # one vertex near the opposite face
        A = rng.standard_normal((4, 3)); A[:, 2] = 0
        w = rng.dirichlet(np.ones(3)); A[3] = w @ A[:3]; A[3, 2] = 10 ** rng.uniform(-3, 0)
    else:               # elongated
        A = rng.standard_normal((4, 3)) * np.r_[1, 10 ** rng.uniform(-2, 0), 10 ** rng.uniform(-2, 0)]
    return A

def work(wid):
    rng = np.random.default_rng(int(sys.argv[2]) * 1000 + wid)
    C = Counter(); hits = []; t0 = time.time(); nshape = 0
    while time.time() - t0 < float(sys.argv[1]):
        A = rand_tet(rng); nshape += 1
        L = rng.dirichlet(np.full(4, 0.4), size=4000)
        X = L @ A
        z, D, pos = zeta_D(A, X)
        ok = pos & np.isfinite(D) & (D < 0)
        if not ok.any():
            C['no R1'] += 1; continue
        for k in rng.choice(np.where(ok)[0], size=min(6, ok.sum()), replace=False):
            try:
                n, m1, m2, good = count_fast(A, z[k], seed=int(rng.integers(1 << 30)))
            except Exception:
                continue
            C[(n, good)] += 1
            if good and n >= 10:
                hits.append(dict(A=A.tolist(), q=z[k].tolist(), n=n)); print('HIT', n, m1, m2, flush=True)
    return C, hits, nshape

if __name__ == '__main__':
    with Pool(int(sys.argv[3]) if len(sys.argv) > 3 else 4) as p:
        res = p.map(work, range(int(sys.argv[3]) if len(sys.argv) > 3 else 4))
    C = Counter(); hits = []; ns = 0
    for c, h, n in res:
        C.update(c); hits += h; ns += n
    print('shapes', ns)
    print(sorted(((str(k), v) for k, v in C.items())))
    json.dump(dict(hist={str(k): v for k, v in C.items()}, hits=hits, shapes=ns), open('sample4_inverse_%s.json' % sys.argv[2], 'w'))
