"""Second ABK iteration (n=7): candidate search and Arb certificate builder.

Design (ABK Remark 3 / Proposition 1, arXiv:2607.27197): the centre of C5 has, in
circumradius-1 coordinates at the centre, alpha ~= 0.00116 and beta = 15/(8 sqrt(6)).
A second axial pair at normalised distance eta with charge alpha*eta^3 - gamma*eta^5,
gamma = beta^2/(30*alpha), bifurcates the centre equilibrium into the 21 critical
points of the ABK local model, for a total of 24 - 1 + 21 = 44.

This script uses per-coordinate offset s = 1/2048 (eta = s/sqrt(2)) and charge
q7 = 134713/2^60, the centre of the numerically found 44-count window
k in [134615, 134812] for q7 = k/2^60 at this offset. It finds the equilibria with
Newton iteration seeded from the C5 certificate points and the local-model cluster
(local_poly.all_roots), builds a Krawczyk box around each zero in Arb arithmetic,
tightens the enclosures exactly as verify_arb.py does, and writes
certificates/n7-abk-second-pair.json. The proof is then verify_arb.py run on that
file (per-box Krawczyk, disjoint enclosures, completeness subdivision).

The float certifier (certify.py) does not currently finish on this configuration
(box budget and memory; see the issue tracker), so unlike C3-C6 this certificate is
not reproduced by make_certs.py.

Run from the repository root:  python3 explore/second_pair.py
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # explore/ for local_poly
import numpy as np
from fractions import Fraction as F
from charges import grad_hess, newton_local, newton, dedup_r
import verify_arb
from verify_arb import Config, Q, _dn, _up
from local_poly import all_roots, harmonic_cubic

T = np.array([[3., 0, 0], [0, 3., 0], [0, 0, 3.]]); C = np.ones(3)
Q5V = 541 / 131072; S7 = 1 / 2048; Q7V = 134713 / 2 ** 60
A5 = np.vstack([T, [1.25] * 3, [0.75] * 3]); Q5 = np.array([1., 1, 1, Q5V, Q5V])
A7 = np.vstack([A5, C + S7, C - S7]); Q7 = np.r_[Q5, Q7V, Q7V]
LAM = 263.876  # 2*beta/(5*alpha) at the C5 centre (circumradius-1 frame)


def phi0_roots():
    Qm = (5 / 32) / (3 / 16) * 2 * np.diag([1, 1, -2.])
    Cc = harmonic_cubic([15 / 8 / (3 / 16), 0, 0, 0, 0, 0, 0])
    _, nr, real = all_roots(np.zeros(3), Qm, Cc, nseed=6000, seed=1)
    assert nr == 21, nr
    return real


def find44(seed=5):
    cert5 = json.load(open('certificates/n5-abk-bipyramid.json'))
    mids = np.array([(np.array(z['enclosure'][0]) + np.array(z['enclosure'][1])) / 2
                     for z in cert5['zeros']])
    E = np.array([[2, -1, -1], [0, 1, -1], [1, 1, 1]], float)
    E = (E / np.linalg.norm(E, axis=1)[:, None]).T
    eta = S7 / np.sqrt(2.0)
    seeds = [mids]
    for fac in (0.5, 0.8, 1.0, 1.3, 1.8):
        Y = (LAM * eta ** 2 * fac) * phi0_roots()
        seeds.append(C + (np.sqrt(6.0) * (E @ Y.T)).T)
    rng = np.random.default_rng(seed)
    for rad in (1e-4, 1e-3, 1e-2):
        u = rng.standard_normal((200, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
        seeds.append(C + u * rad * rng.uniform(0.2, 1.5, (200, 1)))
    X = newton_local(np.vstack(seeds), A7, Q7, iters=100)
    X = newton(X, A7, Q7, iters=30, maxstep=0.02)

    def accept(X):
        g, _ = grad_hess(X, A7, Q7)
        r2 = ((X[:, None] - A7[None]) ** 2).sum(2)
        fs = (Q7[None] / r2).sum(1)
        ok = np.isfinite(X).all(1) & (np.linalg.norm(g, axis=1) < 1e-10 * fs) & (np.sqrt(r2).min(1) > 1e-11)
        return X[ok]

    found = dedup_r(accept(X), A7)
    rmin = np.sqrt(((found[:, None] - A7[None]) ** 2).sum(2)).min(1)
    jit = np.repeat(found, 6, 0) + 0.1 * np.repeat(rmin, 6)[:, None] * rng.standard_normal((6 * len(found), 3))
    found = dedup_r(np.vstack([found, accept(newton_local(jit, A7, Q7, iters=80))]), A7)
    _, H = grad_hess(found, A7, Q7)
    idx = (np.linalg.eigvalsh(H) < 0).sum(1)
    return found, idx


def main():
    t0 = time.time()
    found, idx = find44()
    m1, m2 = int((idx == 1).sum()), int((idx == 2).sum())
    print(f'numerical: {len(found)} equilibria, index1={m1} index2={m2}', flush=True)
    assert len(found) == 44 and (m1, m2) == (19, 25), 'unexpected numerical count; not writing a certificate'

    def exact(x):
        return str(F(float(x)))

    sites = [[exact(x) for x in p] for p in A7]; charges = [exact(x) for x in Q7]
    cfg = Config(sites, charges)
    D = np.sqrt(((found[:, None] - found[None]) ** 2).sum(2)); np.fill_diagonal(D, np.inf)
    sep = D.min(1); dsite = np.sqrt(((found[:, None] - A7[None]) ** 2).sum(2)).min(1)
    zeros = []
    for i, p in enumerate(found):
        wmax = 0.3 * min(sep[i], dsite[i]); box = None
        for f in (0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.5, 1.0):
            w = wmax * f
            lo = [Q(v) for v in (p - w / 2)]; hi = [Q(v) for v in (p + w / 2)]
            try:
                ok, _ = verify_arb.krawczyk(cfg, lo, hi)
            except Exception:
                ok = False
            if ok:
                box = (lo, hi, [[float(v) for v in (p - w / 2)], [float(v) for v in (p + w / 2)]]); break
        assert box is not None, f'no Krawczyk box passed for zero {i}'
        lo, hi, boxf = box
        elo, ehi = lo, hi
        for _ in range(80):
            _, K = verify_arb.krawczyk(cfg, elo, ehi)
            if K is None:
                break
            nlo = [max(elo[k], Q(_dn(K[k]))) for k in range(3)]
            nhi = [min(ehi[k], Q(_up(K[k]))) for k in range(3)]
            assert all(nlo[k] <= nhi[k] for k in range(3))
            if nlo == elo and nhi == ehi:
                break
            elo, ehi = nlo, nhi
        assert verify_arb.morse_index(cfg, elo, ehi) == int(idx[i])
        zeros.append(dict(box=boxf, enclosure=[[float(v) for v in elo], [float(v) for v in ehi]],
                          index=int(idx[i]), boxes=1))
    zeros.sort(key=lambda z: z['enclosure'][0])
    cert = dict(name='n7-abk-second-pair',
                claim='grad V has exactly 44 zeros in R^3 minus the sites; all non-degenerate',
                sites=sites, charges=charges, sites_float=A7.tolist(), charges_float=Q7.tolist(),
                count=44, index_counts={'1': m1, '2': m2},
                morse_relation_m2_minus_m1=m2 - m1, n_minus_1=6,
                stats=dict(builder='explore/second_pair.py (targeted Newton + Arb Krawczyk boxes); '
                                   'proof by verify_arb.py completeness subdivision'),
                zeros=zeros)
    out = 'certificates/n7-abk-second-pair.json'
    json.dump(cert, open(out, 'w'), indent=1)
    print(f'wrote {out} in {time.time() - t0:.1f}s; now prove it with: python3 verify_arb.py {out}', flush=True)


if __name__ == '__main__':
    main()
