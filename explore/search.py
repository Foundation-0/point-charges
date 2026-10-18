"""Ghost-guided hill climbing over configurations of n positive charges.

State: positions a (n,3) and log-charges u (n). Score: (number of equilibria,
-smallest ghost value). A move is accepted when the Morse-consistent count does not
drop and the score does not get worse.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np, sys, json, time
from multiprocessing import Pool
from charges import summary
from ghosts import minima


def evaluate(a, q, seed):
    s, pts, idx = summary(a, q, nseed=1500, rounds=2, seed=seed)
    if not s['morse_ok'] or s['other']:
        return None
    x, v = minima(a, q, nseed=1500, seed=seed)
    gh = v[v > 1e-6]
    return s['n'], (float(gh.min()) if len(gh) else 1.0), s


def chain(args):
    cid, n, start, steps, sigma0 = args
    rng = np.random.default_rng(1000 + cid)
    a, u = start
    a = np.array(a, float); u = np.array(u, float)
    best = evaluate(a, np.exp(u), cid)
    while best is None:
        a = a + 0.01 * rng.standard_normal(a.shape); best = evaluate(a, np.exp(u), cid)
    sigma = sigma0
    log = []
    t0 = time.time()
    for it in range(steps):
        scale = np.max(np.ptp(a, 0))
        k = rng.integers(3)
        a2, u2 = a.copy(), u.copy()
        if k == 0:
            a2 += sigma * scale * rng.standard_normal(a.shape)
        elif k == 1:
            u2 += 5 * sigma * rng.standard_normal(u.shape)
        else:
            i = rng.integers(n); a2[i] += sigma * scale * rng.standard_normal(3); u2[i] += 5 * sigma * rng.standard_normal()
        r = evaluate(a2, np.exp(u2), it)
        if r is None:
            continue
        if (r[0] > best[0]) or (r[0] == best[0] and r[1] <= best[1]):
            if r[0] > best[0]:
                log.append(dict(it=it, n=r[0], a=a2.tolist(), u=u2.tolist(), s=r[2]))
                print(f'chain {cid} it {it}: count {r[0]} {r[2]} t={time.time()-t0:.0f}s', flush=True)
            a, u, best = a2, u2, r
            sigma = min(sigma * 1.2, 0.2)
        else:
            sigma = max(sigma * 0.9, 1e-4)
    return cid, best[0], best[1], a.tolist(), u.tolist(), log


if __name__ == '__main__':
    n = int(sys.argv[1]); steps = int(sys.argv[2]); nch = int(sys.argv[3])
    rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 0)
    starts = []
    s3 = np.sqrt(3)
    for c in range(nch):
        if c % 2 == 0 and n == 4:
            a = [[1, 0, 0], [-0.5, s3/2, 0], [-0.5, -s3/2, 0], [0, 0, 0]]
            u = [0, 0, 0, np.log(10 ** rng.uniform(-4, -2))]
            a = (np.array(a) + 0.02 * rng.standard_normal((4, 3))).tolist()
        else:
            a = rng.standard_normal((n, 3)).tolist(); u = rng.uniform(-3, 0, n).tolist()
        starts.append((c, n, (a, u), steps, 0.02))
    with Pool(4) as pool:
        res = pool.map(chain, starts)
    out = []
    for cid, cnt, g, a, u, log in res:
        print('chain', cid, 'final count', cnt, 'ghost', g)
        out.append(dict(cid=cid, count=cnt, ghost=g, a=a, u=u, log=log))
    json.dump(out, open(f'search_n{n}_{int(time.time())}.json', 'w'))
