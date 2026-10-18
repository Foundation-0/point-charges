"""C3v configurations: unit charges on the unit triangle in z=0, one charge q at (0,0,h).
Critical points in the mirror plane y=0 (x>0 toward a vertex, x<0 toward an edge)
count 1 each on the axis and 3 each off it."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np
from charges import grad_hess
s3 = np.sqrt(3)
tri = np.array([[1, 0, 0], [-0.5, s3/2, 0], [-0.5, -s3/2, 0]])

def slice_cps(h, q, nseed=3000, seed=0):
    a = np.vstack([tri, [0, 0, h]]); qq = np.array([1, 1, 1, q])
    rng = np.random.default_rng(seed)
    # multiscale seeds around the origin and the axial charge
    rad = 10 ** rng.uniform(-4, 0.1, nseed); ang = rng.uniform(0, 2*np.pi, nseed)
    cen = np.where(rng.random(nseed) < 0.5, 0.0, h)
    X = np.c_[rad*np.cos(ang), 0*rad, cen + rad*np.sin(ang)]
    for _ in range(80):
        g, H = grad_hess(X, a, qq)
        g2 = g[:, [0, 2]]; H2 = H[:, [0, 2]][:, :, [0, 2]]
        det = H2[:, 0, 0]*H2[:, 1, 1] - H2[:, 0, 1]*H2[:, 1, 0]
        det = np.where(np.abs(det) < 1e-300, 1e-300, det)
        sx = ( H2[:, 1, 1]*g2[:, 0] - H2[:, 0, 1]*g2[:, 1]) / det
        sz = (-H2[:, 1, 0]*g2[:, 0] + H2[:, 0, 0]*g2[:, 1]) / det
        st = np.c_[sx, sz]; n = np.linalg.norm(st, axis=1)
        f = np.minimum(1, 0.05 / np.maximum(n, 1e-300))
        X[:, 0] -= f*sx; X[:, 2] -= f*sz
    g, H = grad_hess(X, a, qq)
    d = X[:, None, :] - a[None]; r2 = (d*d).sum(2)
    fs = (qq/r2).sum(1)
    ok = np.isfinite(X).all(1) & (np.linalg.norm(g, axis=1) < 1e-10*fs) & (np.sqrt(r2).min(1) > 1e-7)
    X = X[ok]
    pts = []
    for p in X[np.lexsort(X.T[::-1])]:
        if not pts or min(np.linalg.norm(np.array(pts) - p, axis=1)) > 1e-8:
            pts.append(p)
    pts = np.array(pts)
    if len(pts) == 0:
        return 0, pts
    axis = np.abs(pts[:, 0]) < 1e-9
    return int(axis.sum() + 3*(~axis).sum()), pts
