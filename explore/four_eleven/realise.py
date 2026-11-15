"""Turn a hit of the cubic model (flat configuration n, traceless Q, vector l with 8 real local solutions)
into an actual four-charge configuration and count its real balance points numerically.
Usage: realise.py <hits.json> <how many hits to try>"""
import numpy as np, sys, json
from flat_param import flat_config, cubic_parts, field
import quadrics as qd
rng = np.random.default_rng(11)
hits = json.load(open(sys.argv[1])); NT = int(sys.argv[2])
def cond(u):
    A = u[:12].reshape(4, 3); q = u[12:]
    g, H = field(np.zeros((1, 3)), A, q); g = g[0]; H = H[0]
    return np.array([g[0], g[1], g[2], H[0, 0], H[1, 1], H[0, 1], H[0, 2], H[1, 2]])
def jac(u, h=1e-6):
    J = np.empty((8, 16))
    for k in range(16):
        e = np.zeros(16); e[k] = h; J[:, k] = (cond(u + e) - cond(u - e)) / (2 * h)
    return J
def project(u, target):
    for it in range(60):
        cnd = cond(u) - target
        if np.linalg.norm(cnd) < 1e-15: break
        u = u - np.linalg.lstsq(jac(u), cnd, rcond=None)[0]
    return u, np.linalg.norm(cond(u) - target)
def newton(X, A, q, iters=80, cap=None):
    for it in range(iters):
        g, H = field(X, A, q)
        d = np.linalg.solve(H, g[..., None])[..., 0]
        if cap is not None:
            s = np.linalg.norm(d, axis=1); d *= np.minimum(1.0, cap / np.maximum(s, 1e-300))[:, None]
        X = X - d
    return X
def all_points(A, q, seeds, nrand=6000):
    R = np.linalg.norm(A, axis=1).max()
    X0 = rng.normal(size=(nrand, 3)); X0 *= (R * rng.uniform(0, 1, size=(nrand, 1)) ** (1 / 3)) / np.linalg.norm(X0, axis=1)[:, None]
    X0 = np.vstack([seeds, X0])
    with np.errstate(all='ignore'):
        X = newton(X0, A, q, cap=0.05 * R)
        g, H = field(X, A, q)
    gs = np.abs(q).sum() / R**2
    ok = np.isfinite(X).all(1) & (np.linalg.norm(g, axis=1) < 1e-11 * gs) & (np.linalg.norm(X, axis=1) < 1.01 * R)
    pts = []
    for x, h in zip(X[ok], H[ok]):
        if any(np.linalg.norm(x - p[0]) < 1e-7 * R for p in pts): continue
        w = np.linalg.eigvalsh(h); pts.append((x, int((w < 0).sum()), float(np.abs(w).min() / np.abs(w).max())))
    return pts
# rank hits by how far the 8 real roots are from merging
def robustness(n, Qm, l):
    A, q, _ = flat_config(np.array(n)); U, c = cubic_parts(A, q); cmax = c.max(); c = c / cmax
    T = qd.tensor(U[None], c[None]); Qa = np.array(Qm)[None]; la = np.array(l)[None]
    X, good, gap = qd.solve_tql(T, Qa, la, rng); X, res = qd.polish(X, T, Qa, la)
    nreal, ok = qd.real_count(X, res)
    if not (ok[0] and nreal[0] == 8): return -1, None
    X = X[0].real; J = np.einsum('abc,kc->kab', T[0], X) + Qa[0]
    ev = np.abs(np.linalg.eigvalsh(J)); sep = min(np.linalg.norm(X[i] - X[j]) for i in range(8) for j in range(i))
    return float((ev.min(1) / ev.max(1)).min() * sep / np.abs(X).max()), X
ranked = sorted(((robustness(h[1], h[2], h[3])[0], i) for i, h in enumerate(hits)), reverse=True)
print("hits:", len(hits), " best robustness scores:", [round(r, 4) for r, _ in ranked[:5]])
best = None
for score, i in ranked[:NT]:
    _, n, Qm, l = hits[i]; n = np.array(n); Qm = np.array(Qm); l = np.array(l)
    _, Y = robustness(n, Qm, l)
    A0, q0, _ = flat_config(n); U, c = cubic_parts(A0, q0); cmax = c.max()
    u0 = np.concatenate([A0.ravel(), q0])
    for eta in (0.05, 0.02, 0.01, 0.005, 0.002):
        target = cmax * np.array([eta**2 * l[0], eta**2 * l[1], eta**2 * l[2], eta * Qm[0, 0], eta * Qm[1, 1], eta * Qm[0, 1], eta * Qm[0, 2], eta * Qm[1, 2]])
        u, err = project(u0.copy(), target)
        if err > 1e-13 or np.any(u[12:] <= 0): print("  projection failed", err); continue
        A = u[:12].reshape(4, 3); q = u[12:]
        pts = all_points(A, q, eta * Y)
        m1 = sum(1 for p in pts if p[1] == 1); m2 = sum(1 for p in pts if p[1] == 2); m0 = sum(1 for p in pts if p[1] not in (1, 2))
        inner = sum(1 for p in pts if np.linalg.norm(p[0]) < 3 * eta * np.abs(Y).max())
        print(f"hit {i} score {score:.4f} eta {eta}: total {len(pts)}  index-2 {m2}  index-1 {m1}  other {m0}  near flat point {inner}  smallest eigenvalue ratio {min(p[2] for p in pts):.2e}", flush=True)
        if m2 - m1 == 3 and m0 == 0 and (best is None or len(pts) > best[0] or (len(pts) == best[0] and min(p[2] for p in pts) > best[1])):
            best = (len(pts), min(p[2] for p in pts), dict(n=n.tolist(), eta=eta, sites=A.tolist(), charges=q.tolist(), points=[(p[0].tolist(), p[1]) for p in pts]))
if best:
    json.dump(best[2] | {"total": best[0]}, open("realise_best.json", "w"), indent=1); print("best total", best[0], "saved")
