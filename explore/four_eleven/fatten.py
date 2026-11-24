"""Make an 11-point configuration sturdier. Start from a certified configuration and its 11 balance points.
Propose small random changes of the sites and charges, follow the 11 points by Newton, and keep a change when all 11
survive with the same indices and the smallest Hessian eigenvalue among them (in units of total charge / size^3)
grows. Numerical only: the result has to go back through the exact provers.
Usage: fatten.py <config.json> <indep result .json> <iterations> <seed> <out.json>"""
import json, sys, numpy as np
from fractions import Fraction
cfg = json.load(open(sys.argv[1])); res = json.load(open(sys.argv[2])); NIT = int(sys.argv[3]); rng = np.random.default_rng(int(sys.argv[4]))
A = np.array([[float(Fraction(c)) for c in s] for s in cfg["sites"]]); q = np.array([float(Fraction(c)) for c in cfg["charges"]])
P = np.array([z["centre_float"] for z in res["zeros"]])
def gH(X, A, q):
    D = X[:, None, :] - A[None]; r2 = (D * D).sum(-1); r = np.sqrt(r2)
    g = -(q[None, :, None] * D / (r2 * r)[..., None]).sum(1)
    H = (q[None, :, None, None] * (3 * D[..., :, None] * D[..., None, :] - r2[..., None, None] * np.eye(3)) / (r2 * r2 * r)[..., None, None]).sum(1)
    return g, H
def track(P, A, q):
    X = P.copy()
    with np.errstate(all="ignore"):
        for it in range(25):
            g, H = gH(X, A, q)
            try: X = X - np.linalg.solve(H, g[..., None])[..., 0]
            except np.linalg.LinAlgError: return None
        g, H = gH(X, A, q)
    R = np.linalg.norm(A - A.mean(0), axis=1).max(); unit = q.sum() / R**3
    if not np.isfinite(X).all() or np.linalg.norm(g, axis=1).max() > 1e-11 * q.sum() / R**2: return None
    d = np.linalg.norm(X[:, None] - X[None], axis=-1) + np.eye(len(X)) * 1e9
    w = np.linalg.eigvalsh(H); idx = (w < 0).sum(1)
    return X, idx, np.abs(w).min() / unit, d.min() / R
X, idx0, J, sep = track(P, A, q); P = X
print("start: smallest eigenvalue %.3e, smallest separation %.3e, indices %s" % (J, sep, np.bincount(idx0).tolist()))
s = 1e-8; acc = 0
for it in range(NIT):
    R = np.linalg.norm(A - A.mean(0), axis=1).max()
    A2 = A + s * R * rng.normal(size=A.shape); q2 = q * (1 + s * rng.normal(size=4))
    t = track(P, A2, q2)
    if t is not None and (t[1] == idx0).all() and t[2] > J and t[3] > 0.3 * sep:
        A, q, P, J, sep = A2, q2, t[0], t[2], max(sep, t[3]) if False else t[3]; acc += 1; s = min(s * 1.3, 1e-2)
    else:
        s = max(s * 0.97, 1e-10)
    if it % 5000 == 0: print("  iteration %d: smallest eigenvalue %.3e, separation %.3e, step %.1e, accepted %d" % (it, J, sep, s, acc), flush=True)
print("end: smallest eigenvalue %.3e, smallest separation %.3e" % (J, sep))
json.dump({"sites": A.tolist(), "charges": (q / q.max()).tolist(), "points": P.tolist(), "indices": idx0.tolist(), "n": [0, 0, 0, 0], "eta": 0}, open(sys.argv[5], "w"), indent=1)
