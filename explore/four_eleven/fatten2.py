"""Sturdier 11-point configuration, second method: gradient ascent on a smooth minimum of the smallest Hessian
eigenvalues at the 11 tracked balance points (finite-difference gradient, backtracking line search).
Usage: fatten2.py <float config .json with sites, charges, points> <iterations> <out.json>"""
import json, sys, numpy as np
src = json.load(open(sys.argv[1])); NIT = int(sys.argv[2])
A = np.array(src["sites"], float); q = np.array(src["charges"], float); P = np.array(src["points"], float)
def gH(X, A, q):
    D = X[:, None, :] - A[None]; r2 = (D * D).sum(-1); r = np.sqrt(r2)
    g = -(q[None, :, None] * D / (r2 * r)[..., None]).sum(1)
    H = (q[None, :, None, None] * (3 * D[..., :, None] * D[..., None, :] - r2[..., None, None] * np.eye(3)) / (r2 * r2 * r)[..., None, None]).sum(1)
    return g, H
def unpack(v): return v[:12].reshape(4, 3), np.exp(v[12:])
def track(P, v):
    A, q = unpack(v); X = P.copy()
    with np.errstate(all="ignore"):
        for it in range(30):
            g, H = gH(X, A, q)
            try: X = X - np.linalg.solve(H, g[..., None])[..., 0]
            except np.linalg.LinAlgError: return None
        g, H = gH(X, A, q)
    R = np.linalg.norm(A - A.mean(0), axis=1).max(); unit = q.sum() / R**3
    if not np.isfinite(X).all() or np.linalg.norm(g, axis=1).max() > 1e-12 * q.sum() / R**2: return None
    d = np.linalg.norm(X[:, None] - X[None], axis=-1) + np.eye(len(X)) * 1e9
    w = np.linalg.eigvalsh(H); return X, (w < 0).sum(1), np.abs(w).min(1) / unit, d.min() / R
v = np.concatenate([A.ravel(), np.log(q)])
X, idx0, m, sep = track(P, v); P = X
def soft(m, beta): return -np.log(np.exp(-beta * (m - m.min())).sum()) / beta + m.min()
print("start: smallest eigenvalue %.3e, separation %.3e" % (m.min(), sep))
t = 1e-6
for it in range(NIT):
    beta = 8.0 / m.min(); h = 1e-7; g = np.zeros(16); ok = True
    for k in range(16):
        e = np.zeros(16); e[k] = h
        tp = track(P, v + e); tm = track(P, v - e)
        if tp is None or tm is None: ok = False; break
        g[k] = (soft(tp[2], beta) - soft(tm[2], beta)) / (2 * h)
    if not ok: print("gradient failed"); break
    g /= np.linalg.norm(g); f0 = soft(m, beta); moved = False
    for bt in range(40):
        tr = track(P, v + t * g)
        if tr is not None and (tr[1] == idx0).all() and soft(tr[2], beta) > f0 and tr[3] > 0.2 * sep:
            v = v + t * g; P, m, sep = tr[0], tr[2], tr[3]; t *= 1.5; moved = True; break
        t *= 0.5
    if not moved: print("no ascent step at iteration", it); break
    if it % 100 == 0: print("  iteration %d: smallest eigenvalue %.3e, separation %.3e, step %.1e" % (it, m.min(), sep, t), flush=True)
A, q = unpack(v)
print("end: smallest eigenvalue %.3e, separation %.3e; charges %s" % (m.min(), sep, np.round(q / q.max(), 4)))
json.dump({"sites": A.tolist(), "charges": (q / q.max()).tolist(), "points": P.tolist(), "indices": idx0.tolist(), "n": [0, 0, 0, 0], "eta": 0}, open(sys.argv[3], "w"), indent=1)
