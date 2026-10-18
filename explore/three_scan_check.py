"""Check the flagged shapes of three_scan.py: every pixel the rasteriser covered twice or more is
re-solved directly (crit2d counts the minima for that charge vector); a sample of singly covered
pixels is checked the same way. Writes ../results/three_scan_summary.txt."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import json, numpy as np
from collections import Counter
from three_scan import coverage
from three_map import zeta_D
from three import crit2d

def pixel_charges(A, sbox, N=400, P=300):
    s0, s1, t0, t1 = sbox
    sp = lambda a, b: a + (b - a) * (0.5 - 0.5 * np.cos(np.pi * np.linspace(0, 1, N + 1)))
    S, T = np.meshgrid(sp(s0, s1), sp(t0, t1), indexing='ij')
    l2 = S * (1 - T); l3 = S * T; l1 = 1 - S
    X = (l1[..., None] * A[0] + l2[..., None] * A[1] + l3[..., None] * A[2]).reshape(-1, 2)
    ins = ((l1 > 1e-15) & (l2 > 1e-15) & (l3 > 1e-15)).ravel()
    z = np.full((X.shape[0], 3), np.nan); D = np.full(X.shape[0], np.nan)
    z[ins], D[ins] = zeta_D(A, X[ins])
    zz = z[:, :2][D > 0]
    zlo = zz.min(0); w = np.maximum(zz.max(0) - zlo, 1e-300) * 1.02; zlo = zlo - 0.01 * w
    return lambda i, j: (lambda c: np.r_[c, 1 - c.sum()])(zlo + (np.array([i, j]) + 0.5) / P * w)

rows = [json.loads(l) for l in open('../results/three_scan.log')]
lines = []
zoom = Counter(int(r[3]) for r in rows)
bs = sorted(set(r[1] for r in rows))
res = [r for r in rows if r[3] >= 1]
lines.append('three_scan.py: %d triangle shapes (A1=(0,0), A2=(1,0), A3=(a,b); 0<=a<=1/2, |A2A3|<=1, b from 1e-4 to 0.87)' % len(rows))
lines.append('max pixel coverage of Z(R_min) at the fine resolution: %s' % dict(sorted(zoom.items())))
lines.append('shapes with no resolved minimum region: %d (smallest b with a resolved region: %.3g)' % (zoom[0], min(r[1] for r in res)))
rng = np.random.default_rng(0)
for r in rows:
    if r[3] < 2:
        continue
    a, b = r[0], r[1]
    A = np.array([[0, 0], [1, 0], [a, b]], float)
    _, _, sbox, _ = coverage(A, N=200, P=200)
    _, _, _, cov = coverage(A, N=400, P=300, box=sbox)
    zeta_at = pixel_charges(A, sbox)
    multi = np.argwhere(cov >= 2)
    out = Counter()
    for i, j in multi:
        out[crit2d(A, zeta_at(i, j), nseed=4000)[1]] += 1
    one = np.argwhere(cov == 1); sel = one[rng.choice(len(one), min(30, len(one)), replace=False)]
    chk = Counter(crit2d(A, zeta_at(i, j), nseed=4000)[1] for i, j in sel)
    lines.append('flagged shape a=%.3f b=%.4g: %d pixel(s) covered >= 2; direct solve finds minima %s there; 30 singly covered pixels: minima %s'
                 % (a, b, len(multi), dict(out), dict(chk)))
txt = '\n'.join(lines)
print(txt)
open('../results/three_scan_summary.txt', 'w').write(txt + '\n')
