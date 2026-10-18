import numpy as np, sys, json, time
from multiprocessing import Pool
from three_map import zeta_D

def coverage(A, N=500, P=500, box=None):
    """Rasterised count of R_min preimages over a P x P window of the charge simplex."""
    if box is None:
        s0, s1, t0, t1 = 0, 1, 0, 1
    else:
        s0, s1, t0, t1 = box
    s = s0 + (s1 - s0) * (0.5 - 0.5 * np.cos(np.pi * np.linspace(0, 1, N + 1)))
    t = t0 + (t1 - t0) * (0.5 - 0.5 * np.cos(np.pi * np.linspace(0, 1, N + 1)))
    S, T = np.meshgrid(s, t, indexing='ij')
    l2 = S * (1 - T); l3 = S * T; l1 = 1 - S
    X = (l1[..., None] * A[0] + l2[..., None] * A[1] + l3[..., None] * A[2]).reshape(-1, 2)
    inside = ((l1 > 1e-15) & (l2 > 1e-15) & (l3 > 1e-15)).ravel()
    z = np.full((X.shape[0], 3), np.nan); D = np.full(X.shape[0], np.nan)
    z[inside], D[inside] = zeta_D(A, X[inside])
    Z = z[:, :2].reshape(N + 1, N + 1, 2); D = D.reshape(N + 1, N + 1)
    pos = D > 0
    if not pos.any():
        return 0, 0, None, None
    ii, jj = np.where(pos)
    sbox = (s[max(ii.min() - 2, 0)], s[min(ii.max() + 2, N)], t[max(jj.min() - 2, 0)], t[min(jj.max() + 2, N)])
    zz = Z[pos]
    zlo = zz.min(0); zhi = zz.max(0); w = np.maximum(zhi - zlo, 1e-300) * 1.02
    zlo = zlo - 0.01 * w
    cov = np.zeros((P, P), np.int32)
    I, J = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
    for offs in (((0, 0), (1, 0), (1, 1)), ((0, 0), (1, 1), (0, 1))):
        Ps = [(Z[I + di, J + dj] - zlo) / w * P for di, dj in offs]
        Dm = np.minimum.reduce([D[I + di, J + dj] for di, dj in offs])
        good = (Dm > 0) & np.all([np.isfinite(p).all(-1) for p in Ps], 0)
        P0, P1, P2 = [p[good] for p in Ps]
        if not len(P0):
            continue
        lo = np.floor(np.minimum.reduce([P0, P1, P2]) - 0.5).astype(int)
        hi = np.ceil(np.maximum.reduce([P0, P1, P2]) - 0.5).astype(int)
        span = int((hi - lo).max())
        v0 = P1 - P0; v1 = P2 - P0
        den = v0[:, 0] * v1[:, 1] - v0[:, 1] * v1[:, 0]
        for ox in range(span + 1):
            for oy in range(span + 1):
                px = lo[:, 0] + ox; py = lo[:, 1] + oy
                v2 = np.c_[px + 0.5, py + 0.5] - P0
                with np.errstate(divide='ignore', invalid='ignore'):
                    b1 = (v2[:, 0] * v1[:, 1] - v2[:, 1] * v1[:, 0]) / den
                    b2 = (v0[:, 0] * v2[:, 1] - v0[:, 1] * v2[:, 0]) / den
                ins = (b1 >= 0) & (b2 >= 0) & (b1 + b2 < 1) & (px >= 0) & (py >= 0) & (px < P) & (py < P)
                np.add.at(cov, (px[ins], py[ins]), 1)
    return int(cov.max()), int((cov > 0).sum()), sbox, cov

def shape_job(ab):
    a, b = ab
    A = np.array([[0.0, 0.0], [1.0, 0.0], [a, b]])
    mc, npx, sbox, cov = coverage(A, N=200, P=200)
    if sbox is None:
        return a, b, 0, 0, 0
    mc2, npx2, _, cov2 = coverage(A, N=400, P=300, box=sbox)    # zoom on R_min
    return a, b, mc, mc2, npx2

if __name__ == '__main__':
    # fundamental domain of triangle shapes: |A1A2| = 1 is the longest side, A3 nearer A1
    bs = np.r_[10 ** np.linspace(-4, -1, 25), np.linspace(0.11, 0.87, 39)]
    jobs = []
    for b in bs:
        for a in np.linspace(0, 0.5, 26):
            if (a - 1) ** 2 + b ** 2 <= 1 + 1e-12:
                jobs.append((a, b))
    t0 = time.time()
    res = []
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as p, open('three_scan.log', 'w') as log:
        for r in p.imap_unordered(shape_job, jobs, chunksize=4):
            res.append(r)
            log.write(json.dumps([float(x) for x in r]) + '\n'); log.flush()
            if r[3] >= 2:
                print('DOUBLE', r, flush=True)
    print('shapes', len(res), 'time %.0fs' % (time.time() - t0))
    from collections import Counter
    print('max coverage (coarse):', Counter(r[2] for r in res))
    print('max coverage (zoom):', Counter(r[3] for r in res))
    print('no minimum region at all:', sum(1 for r in res if r[2] == 0))
    for r in res:
        if r[3] >= 2: print('DOUBLE', r)
    json.dump(res, open('three_scan.json', 'w'))
