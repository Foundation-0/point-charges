"""Look for four distinct index-1 equilibria with the same charges (=> >= 11 equilibria).

Unknowns: tetrahedron shape (a,b,c,d,e) and four interior points x_k (barycentric logits).
Residuals: log charge ratios of x_k minus those of x_1 (9 numbers), plus penalties keeping
det H(x_k) < 0 (index 1) and the points apart.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np, sys, time, json
from scipy.optimize import least_squares
from multiprocessing import Pool
from four_map import zeta_D

def unpack(p):
    a, b, c, d, e = p[:5]
    A = np.array([[0, 0, 0], [1, 0, 0], [a, b, 0], [c, d, e]], float)
    lg = p[5:].reshape(4, 4)
    l = np.exp(lg - lg.max(1, keepdims=True)); l /= l.sum(1, keepdims=True)
    X = l @ A
    return A, X, l

def resid(p, delta=0.05, sep=0.02):
    A, X, l = unpack(p)
    z, D, pos = zeta_D(A, X)
    L = np.log(np.abs(z[:, :3] / z[:, 3:4]))
    r = [(L[k] - L[0]) for k in range(1, 4)]
    pen = [np.logaddexp(0, (D[k] + delta) * 20) / 20 for k in range(4)]        # want D < -delta
    size = np.max(np.linalg.norm(A[:, None] - A[None], axis=2))
    pr = []
    for i in range(4):
        for j in range(i + 1, 4):
            dist = np.linalg.norm(X[i] - X[j]) / size
            pr.append(np.logaddexp(0, (sep - dist) * 200) / 200 * 10)
    # keep the tetrahedron from degenerating completely
    flat = [np.logaddexp(0, (0.02 - abs(p[4])) * 100) / 100, np.logaddexp(0, (0.02 - abs(p[1])) * 100) / 100]
    return np.r_[np.concatenate(r), pen, pr, flat]

def attempt(seed):
    rng = np.random.default_rng(seed)
    best = None
    for k in range(20):
        shape = np.r_[rng.uniform(-0.5, 1.5), rng.uniform(0.1, 1.2), rng.uniform(-0.5, 1.5), rng.uniform(-1, 1), rng.uniform(0.05, 1.2)]
        p0 = np.r_[shape, rng.normal(0, 1.5, 16)]
        try:
            sol = least_squares(resid, p0, method='trf', max_nfev=400, xtol=1e-14, ftol=1e-14, gtol=1e-14)
        except Exception:
            continue
        A, X, l = unpack(sol.x)
        z, D, pos = zeta_D(A, X)
        L = np.log(np.abs(z[:, :3] / z[:, 3:4]))
        err = max(np.abs(L[k] - L[0]).max() for k in range(1, 4))
        rec = dict(seed=seed, k=k, err=float(err), D=D.tolist(), cost=float(sol.cost), p=sol.x.tolist())
        if best is None or err < best['err'] or (err < 1e-8 and max(D) < max(best['D'])):
            if (D < 0).all():
                best = rec
        if err < 1e-9 and (D < -1e-3).all():
            print('CANDIDATE', seed, k, err, D, flush=True)
    return best

if __name__ == '__main__':
    t0 = time.time()
    with Pool(4) as p:
        res = p.map(attempt, range(int(sys.argv[1]), int(sys.argv[2])))
    res = [r for r in res if r]
    res.sort(key=lambda r: r['err'])
    for r in res[:10]:
        print('err %.3g  D %s' % (r['err'], np.round(r['D'], 4)))
    json.dump(res, open('four_lsq_%s.json' % sys.argv[1], 'w'))
    print('time %.0fs' % (time.time() - t0))
