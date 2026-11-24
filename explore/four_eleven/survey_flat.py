"""Survey of the flat family: for random n in the positive part of the 3-sphere, count the real balance points of
the flat configuration other than the flat point itself. By the degree argument they satisfy
index-2 minus index-1 = 3, so the count is 3, 5, 7, ... A count of 5 or more, together with an all-real unfolding
of the flat point (8), would give 13 or more."""
import numpy as np, sys, json
from flat_param import flat_config, field
rng = np.random.default_rng(int(sys.argv[1])); NS = int(sys.argv[2]); NSTART = 3000
hist = {}; bad = 0; rich = []
for s in range(NS):
    n = np.abs(rng.normal(size=4)); n /= np.linalg.norm(n)
    if n.min() < 0.03: continue
    A, q, _ = flat_config(n); R = np.linalg.norm(A, axis=1).max()
    X = rng.normal(size=(NSTART, 3)); X *= (R * rng.uniform(0, 1, size=(NSTART, 1)) ** (1 / 3)) / np.linalg.norm(X, axis=1)[:, None]
    with np.errstate(all='ignore'):
        for it in range(80):
            g, H = field(X, A, q)
            d = np.einsum('mab,mb->ma', np.linalg.pinv(H), g)
            st = np.linalg.norm(d, axis=1); d *= np.minimum(1.0, 0.05 * R / np.maximum(st, 1e-300))[:, None]
            X = X - d
        g, H = field(X, A, q)
    gs = q.sum() / R**2
    w = np.linalg.eigvalsh(np.where(np.isfinite(H), H, 0.0))
    ok = np.isfinite(X).all(1) & (np.linalg.norm(g, axis=1) < 1e-11 * gs) & (np.linalg.norm(X, axis=1) < R) & (np.abs(w).min(1) > 1e-7 * np.abs(w).max(1))
    pts = []
    for x, ww in zip(X[ok], w[ok]):
        if np.linalg.norm(x) < 1e-5 * R: continue
        if all(np.linalg.norm(x - p[0]) > 1e-7 * R for p in pts): pts.append((x, int((ww < 0).sum())))
    m2 = sum(1 for p in pts if p[1] == 2); m1 = sum(1 for p in pts if p[1] == 1)
    if m2 - m1 != 3: bad += 1; continue
    hist[len(pts)] = hist.get(len(pts), 0) + 1
    if len(pts) >= 5: rich.append((len(pts), n.tolist()))
print("flat configurations surveyed:", sum(hist.values()), " inconsistent (Morse relation failed, skipped):", bad)
print("other balance points (count: configurations):", dict(sorted(hist.items())))
json.dump(rich, open("survey_flat_rich_%s.json" % sys.argv[1], "w"))
