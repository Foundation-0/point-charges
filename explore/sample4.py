"""Random multi-scale sampling of 4-charge configurations; logs histogram and any count >= 10."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np, sys, json, time
from multiprocessing import Pool
from collections import Counter
from charges import count_fast, find_fast

def rand_cfg(rng):
    kind = rng.integers(6)
    if kind == 0:   # comparable scales
        a = rng.standard_normal((4, 3)); q = 10 ** rng.uniform(-4, 0, 4)
    elif kind in (1, 2):  # 3 + 1 near a critical point of the 3-charge field
        a3 = rng.standard_normal((3, 3)); q3 = 10 ** rng.uniform(-2, 0, 3)
        if kind == 2:
            a3[:, 2] = 0; a3 = a3 - a3.mean(0)
        pts, idx, ev = find_fast(a3, q3, rounds=1)
        p = pts[rng.integers(len(pts))] if len(pts) else a3.mean(0)
        sc = np.min(np.linalg.norm(a3 - p, axis=1))
        off = rng.standard_normal(3) * sc * 10 ** rng.uniform(-3, 0)
        a = np.vstack([a3, p + off]); q = np.r_[q3, 10 ** rng.uniform(-8, -0.5)]
    elif kind == 3:  # 2 far + 2 near the far equilibrium
        q2 = np.r_[1.0, 10 ** rng.uniform(-1, 0)]
        a2 = np.array([[-1.0, 0, 0], [1.0, 0, 0]])
        t = np.sqrt(q2[0]) / (np.sqrt(q2[0]) + np.sqrt(q2[1])); p = a2[0] + t * (a2[1] - a2[0])
        s = 10 ** rng.uniform(-3, -0.3)
        u = rng.standard_normal(3); u /= np.linalg.norm(u)
        c = p + rng.standard_normal(3) * s * 10 ** rng.uniform(-2, 0.5)
        a = np.vstack([a2, c + s * u, c - s * u * 10 ** rng.uniform(-1, 1) ** 0])
        a[3] = c - s * u * rng.uniform(0.2, 1.5)
        q = np.r_[q2, 10 ** rng.uniform(-9, -1, 2)]
    elif kind == 4:  # 1 far + 3 near
        s = 10 ** rng.uniform(-3, -0.3)
        a3 = rng.standard_normal((3, 3)) * s; q3 = 10 ** rng.uniform(-2, 0, 3)
        a = np.vstack([a3, [3.0, 0, 0] + rng.standard_normal(3)]); q = np.r_[q3, 10 ** rng.uniform(-2, 3)]
    else:  # nested: 2 far, 1 mid, 1 near
        a = np.array([[-1.0, 0, 0], [1.0, 0, 0]] + [list(rng.standard_normal(3) * 10 ** rng.uniform(-1.5, 0))])
        pts, idx, ev = find_fast(a, np.r_[1, 10 ** rng.uniform(-1, 0), 10 ** rng.uniform(-4, 0)], rounds=1)
        q = np.r_[1, 10 ** rng.uniform(-1, 0), 10 ** rng.uniform(-4, 0)]
        pts, idx, ev = find_fast(a, q, rounds=1)
        p = pts[rng.integers(len(pts))] if len(pts) else a.mean(0)
        sc = np.min(np.linalg.norm(a - p, axis=1))
        a = np.vstack([a, p + rng.standard_normal(3) * sc * 10 ** rng.uniform(-3, 0)]); q = np.r_[q, 10 ** rng.uniform(-9, -1)]
    return kind, a, q

def work(wid):
    rng = np.random.default_rng(int(sys.argv[2]) * 100 + wid)
    C = Counter(); hits = []
    t0 = time.time()
    while time.time() - t0 < float(sys.argv[1]):
        try:
            kind, a, q = rand_cfg(rng)
            n, m1, m2, ok = count_fast(a, q, seed=int(rng.integers(1 << 30)))
        except Exception as e:
            continue
        C[(int(kind), int(n), bool(ok))] += 1
        if ok and n >= 10:
            hits.append(dict(kind=int(kind), n=n, a=a.tolist(), q=q.tolist()))
            print('HIT', kind, n, m1, m2, flush=True)
    return C, hits

if __name__ == '__main__':
    with Pool(4) as p:
        res = p.map(work, range(4))
    C = Counter(); hits = []
    for c, h in res:
        C.update(c); hits += h
    tot = Counter()
    for (k, n, ok), v in sorted(C.items()):
        if ok: tot[n] += v
    print('Morse-consistent histogram:', sorted(tot.items()))
    print('by kind max:', {k: max(n for (kk, n, ok) in C if kk == k and ok) for k in range(6) if any(kk == k and ok for (kk, n, ok) in C)})
    print('inconsistent:', sum(v for (k, n, ok), v in C.items() if not ok), 'total', sum(C.values()))
    json.dump(dict(hist=[[list(k), v] for k, v in C.items()], hits=hits), open(f'sample4_{sys.argv[2]}.json', 'w'))
