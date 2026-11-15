import numpy as np, sys, json
from flat_param import flat_config, cubic_parts
import quadrics as qd
rng = np.random.default_rng(int(sys.argv[1])); NS = int(sys.argv[2]); SPREAD = float(sys.argv[3]); B = 20000
hist = np.zeros(9, int); incomplete = 0; hits = []
for b in range(0, NS, B):
    ns = np.abs(0.5 + SPREAD * rng.normal(size=(B, 4))) if SPREAD < 5 else np.abs(rng.normal(size=(B, 4)))
    ns /= np.linalg.norm(ns, axis=1)[:, None]
    t = ns / np.sqrt(1 - ns**2)
    G = -np.einsum('si,sj->sij', t**2, t**2); G[:, range(4), range(4)] = t**2
    w, V = np.linalg.eigh(G); A = V[:, :, 1:] * np.sqrt(np.maximum(w[:, None, 1:], 0))
    r = np.linalg.norm(A, axis=2); U = A / r[..., None]; q = t**3 / (1 + t**2); c = q / r**4; c /= c.max(1)[:, None]
    T = qd.tensor(U, c)
    Qm = rng.normal(size=(B, 3, 3)); Qm = (Qm + Qm.transpose(0, 2, 1)) / 2; Qm -= np.eye(3) * (np.trace(Qm, axis1=1, axis2=2) / 3)[:, None, None]
    Qm /= np.abs(np.linalg.eigvalsh(Qm)).max(-1)[:, None, None]
    l = rng.normal(size=(B, 3)); l *= (10 ** rng.uniform(-3, 1, size=(B, 1))) / np.linalg.norm(l, axis=1)[:, None]
    X, good, gap = qd.solve_tql(T, Qm, l, rng); X, res = qd.polish(X, T, Qm, l)
    nreal, ok = qd.real_count(X, res); ok &= good
    incomplete += (~ok).sum()
    for k in range(9): hist[k] += ((nreal == k) & ok).sum()
    for s in np.where(ok & (nreal >= 8))[0]: hits.append((int(nreal[s]), ns[s].tolist(), Qm[s].tolist(), l[s].tolist()))
print("spread", SPREAD, "samples", NS, "incomplete", incomplete)
print("real-solution counts 0..8:", hist.tolist())
json.dump(hits[:500], open("cubic_scan2_hits_%s_%s.json" % (sys.argv[1], sys.argv[3]), "w"))
