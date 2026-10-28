"""Wave-1 bounded search for a 4-charge configuration with >9 equilibria.

Issue: Foundation-0/point-charges#3. With non-degenerate equilibria the
count is 3+2*m1, so >9 means >=11 (m1>=4).

Method: screen with charges.count_fast (multi-start Newton + Morse check
m2-m1=3); any configuration screening >=9 is re-checked with two more
count_fast seeds and with the slower charges.summary. Only Morse-consistent
counts are recorded as counts; inconsistent screens are tallied separately
(they indicate solver misses, not equilibria).

Families (all deterministic, seed SEED):
  tri_interior   equilateral/isosceles/scalene triangle + interior 4th site
                 on a position x log-charge grid
  c3v_axial      equilateral triangle + 4th site on the symmetry axis
  symmetric      regular/stretched tetrahedron, square, rhombus, disphenoid
                 with structured charge patterns
  two_pair       two pairs at different scales/separations (issue suggestion)
  near_coplanar  flat tetrahedra: 3 sites in a plane + 4th at small height,
                 and 4 near-planar sites (least-squares drift direction)
  insertion      tiny 4th charge near an index-1 equilibrium of a 3-charge
                 configuration (Proposition-2 style, off-centre variants)
  random         multi-scale random positions and log-uniform charges
  hill_climb     mutation hill-climbing seeded from the certified C4 config
                 (9 equilibria) and from any 9 found in the sweep

Outputs: results/wave1_maxwell_n4.json and results/wave1_maxwell_n4.md
Usage: .venv/bin/python explore/wave1_maxwell_n4.py [--quick]
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import json, time
from collections import Counter
from multiprocessing import Pool
import numpy as np
from charges import count_fast, summary, find_fast

SEED = 20261018
QUICK = "--quick" in sys.argv


def screen(task):
    fam, a, q = task
    try:
        n, m1, m2, ok = count_fast(np.array(a), np.array(q), seed=0)
    except Exception as e:  # solver failure: tally, do not crash the sweep
        return fam, a, q, -1, 0, 0, False, str(e)[:80]
    return fam, a, q, n, m1, m2, ok, ""


# ---------------- configuration generators ----------------

def tri_equilateral(R=1.0):
    th = np.array([np.pi / 2, np.pi / 2 + 2 * np.pi / 3, np.pi / 2 + 4 * np.pi / 3])
    return np.c_[R * np.cos(th), R * np.sin(th), np.zeros(3)]


def gen_tri_interior(rng):
    tasks = []
    tris = [tri_equilateral(1.0)]
    # isosceles triangles: apex half-angle sweep
    for s in (0.3, 0.6, 1.5, 2.5):
        tris.append(np.array([[1.0, 0, 0], [-0.5, s, 0], [-0.5, -s, 0]]))
    # scalene triangles (seeded)
    for _ in range(4):
        tris.append(np.c_[rng.uniform(-1, 1, 3), rng.uniform(-1, 1, 3), np.zeros(3)])
    # interior positions: barycentric grid
    bary = []
    for i in range(5):
        for j in range(5 - i):
            w = np.array([i, j, 4 - i - j]) / 4.0
            bary.append(w)
    qs4 = 10 ** np.array([-6, -5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0.0])
    heights = [0.0, 0.05, -0.05, 0.2, -0.2, 0.5]
    for T in tris:
        for w in bary:
            p = w @ T
            for h in heights:
                a = np.vstack([T, p + [0, 0, h]])
                for q4 in qs4:
                    tasks.append(("tri_interior", a.tolist(), [1.0, 1.0, 1.0, float(q4)]))
    # unequal triangle charges on the equilateral triangle, centre-ish 4th site
    T = tri_equilateral(1.0)
    for qr in [(1, 1, 0.1), (1, 0.5, 0.1), (1, 0.1, 0.01), (1, 2, 4), (1, 0.3, 3)]:
        for w in bary[::3]:
            for q4 in qs4[::2]:
                a = np.vstack([T, w @ T])
                tasks.append(("tri_interior", a.tolist(), [*qr, float(q4)]))
    return tasks


def gen_c3v_axial(rng):
    tasks = []
    T = tri_equilateral(1.0)
    for h in np.r_[np.linspace(-1.5, 1.5, 25), 10 ** np.linspace(-3, -0.3, 8), -10 ** np.linspace(-3, -0.3, 8)]:
        for lq in np.linspace(-6, 1, 29):
            a = np.vstack([T, [0, 0, h]])
            tasks.append(("c3v_axial", a.tolist(), [1.0, 1.0, 1.0, float(10 ** lq)]))
    return tasks


def gen_symmetric(rng):
    tasks = []
    tet = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
    charge_patterns = [[1, 1, 1, 1]]
    for r in 10 ** np.linspace(-4, 1, 11):
        charge_patterns += [[1, 1, 1, r], [1, 1, r, r], [1, r, r, r], [1, r, r * r, r ** 3]]
    for s in (0.2, 0.5, 0.8, 1.0, 1.5, 2.5, 4.0):  # stretch along z
        a = tet * [1, 1, s]
        for q in charge_patterns:
            tasks.append(("symmetric_tet", a.tolist(), [float(x) for x in q]))
    for rho in np.linspace(0.2, 2.0, 10):  # rhombus / square in plane
        a = np.array([[1, 0, 0], [0, rho, 0], [-1, 0, 0], [0, -rho, 0]], float)
        for q in charge_patterns:
            tasks.append(("symmetric_rhombus", a.tolist(), [float(x) for x in q]))
    for s in np.linspace(0.1, 2.0, 10):  # disphenoid
        a = np.array([[1, 0, s], [-1, 0, s], [0, 1, -s], [0, -1, -s]], float)
        for q in charge_patterns[::2]:
            tasks.append(("symmetric_disphenoid", a.tolist(), [float(x) for x in q]))
    return tasks


def gen_two_pair(rng):
    """Two pairs at different scales: pair A wide, pair B tight, varied offset."""
    tasks = []
    for sepB in 10 ** np.linspace(-2.5, 0, 8):
        for off in (0.0, 0.2, 0.5, 1.0):
            for ang in (0.0, np.pi / 4, np.pi / 2):
                u = np.array([np.cos(ang), np.sin(ang), 0.3]); u /= np.linalg.norm(u)
                cB = np.array([off, 0.1 * off, 0.2 * off])
                a = np.vstack([[-1, 0, 0], [1, 0, 0], cB - sepB * u, cB + sepB * u])
                for lq in ((-4, -4), (-3, -1), (-2, -2), (-1, -3), (0, 0), (-5, -2)):
                    q = [1.0, 1.0, float(10 ** lq[0]), float(10 ** lq[1])]
                    tasks.append(("two_pair", a.tolist(), q))
                    q2 = [float(10 ** lq[0]), float(10 ** lq[1]), 1.0, 1.0]
                    tasks.append(("two_pair", a.tolist(), q2))
    # random two-pair: two random segments, scales differing by up to 100x
    n = 300 if QUICK else 1200
    for _ in range(n):
        c1, c2 = rng.standard_normal(3), rng.standard_normal(3)
        u1 = rng.standard_normal(3); u1 /= np.linalg.norm(u1)
        u2 = rng.standard_normal(3); u2 /= np.linalg.norm(u2)
        s1, s2 = 10 ** rng.uniform(-0.3, 0.3), 10 ** rng.uniform(-2.5, 0)
        a = np.vstack([c1 - s1 * u1, c1 + s1 * u1, c2 - s2 * u2, c2 + s2 * u2])
        q = 10 ** rng.uniform(-5, 0, 4)
        tasks.append(("two_pair_rand", a.tolist(), q.tolist()))
    return tasks


def gen_near_coplanar(rng):
    tasks = []
    n = 400 if QUICK else 1600
    for _ in range(n):
        a = rng.standard_normal((4, 3))
        a[:, 2] *= 10 ** rng.uniform(-4, -0.5)  # flatness
        q = 10 ** rng.uniform(-6, 0, 4)
        tasks.append(("near_coplanar", a.tolist(), q.tolist()))
    # exactly planar quadrilateral + one site lifted by tiny height
    for h in 10 ** np.linspace(-4, -0.5, 8):
        for lq in np.linspace(-5, 0, 11):
            a = np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, h]], float)
            tasks.append(("near_coplanar", a.tolist(), [1.0, 1.0, 1.0, float(10 ** lq)]))
    return tasks


def gen_insertion(rng):
    """Tiny 4th charge near an index-1 point of a 3-charge potential."""
    tasks = []
    bases = []
    T = tri_equilateral(1.0)
    bases.append((T, np.array([1.0, 1.0, 1.0])))
    for qr in ([1, 1, 0.3], [1, 0.5, 2.0], [1, 0.2, 0.2]):
        bases.append((T, np.array(qr, float)))
    for s in (0.5, 1.5):
        bases.append((np.array([[1.0, 0, 0], [-0.5, s, 0], [-0.5, -s, 0]]), np.ones(3)))
    for T3, q3 in bases:
        pts, idx, ev = find_fast(T3, q3, seed=0)
        saddles = pts[idx == 1]
        for p in saddles:
            for kappa in 10 ** np.linspace(-7, -1, 13):
                # on the saddle, and displaced along random directions
                tasks.append(("insertion", np.vstack([T3, p]).tolist(), [*q3.tolist(), float(kappa)]))
                for _ in range(3):
                    d = rng.standard_normal(3); d /= np.linalg.norm(d)
                    off = d * (kappa ** (1 / 3)) * rng.uniform(0.2, 3.0)
                    tasks.append(("insertion", np.vstack([T3, p + off]).tolist(), [*q3.tolist(), float(kappa)]))
    return tasks


def gen_random(rng):
    tasks = []
    n = 800 if QUICK else 4000
    for i in range(n):
        mode = i % 5
        if mode == 0:
            a = rng.standard_normal((4, 3))
        elif mode == 1:
            a = rng.uniform(-1, 1, (4, 3))
        elif mode == 2:  # one site near centroid of other three, small scale
            a = rng.standard_normal((4, 3)); a[3] = a[:3].mean(0) + 10 ** rng.uniform(-3, -0.3) * rng.standard_normal(3)
        elif mode == 3:  # hierarchical scales
            a = np.cumsum([rng.standard_normal(3) * 10 ** rng.uniform(-2, 0) for _ in range(4)], axis=0)
        else:  # near-collinear
            u = rng.standard_normal(3); u /= np.linalg.norm(u)
            a = np.sort(rng.uniform(-1, 1, 4))[:, None] * u + 10 ** rng.uniform(-3, -1) * rng.standard_normal((4, 3))
        q = 10 ** rng.uniform(-6, 0, 4)
        tasks.append(("random", a.tolist(), q.tolist()))
    return tasks


# ---------------- hill climbing ----------------

def hill_chain(args):
    cid, a0, q0, steps = args
    rng = np.random.default_rng(SEED + 777 + cid)
    a, u = np.array(a0, float), np.log(np.array(q0, float))
    cur = count_fast(a, np.exp(u), seed=cid)
    best_n, hist = cur[0] if cur[3] else -1, []
    sigma = 0.03
    evals = 1
    for it in range(steps):
        a2, u2 = a.copy(), u.copy()
        k = rng.integers(3)
        sc = max(np.max(np.ptp(a, 0)), 1e-9)
        if k == 0:
            a2 += sigma * sc * rng.standard_normal(a.shape)
        elif k == 1:
            u2 += 4 * sigma * rng.standard_normal(u.shape)
        else:
            i = rng.integers(4); a2[i] += sigma * sc * rng.standard_normal(3); u2[i] += 4 * sigma * rng.standard_normal()
        r = count_fast(a2, np.exp(u2), seed=int(rng.integers(1 << 30)))
        evals += 1
        if not r[3]:
            continue
        if r[0] >= best_n:
            if r[0] > best_n:
                hist.append(dict(it=it, n=int(r[0]), a=a2.tolist(), q=np.exp(u2).tolist()))
            a, u, best_n = a2, u2, r[0]
            sigma = min(sigma * 1.15, 0.25)
        else:
            sigma = max(sigma * 0.92, 1e-4)
    return dict(cid=cid, best_n=int(best_n), evals=evals, a=a.tolist(), q=np.exp(u).tolist(), improvements=hist)


def verify_candidate(a, q):
    """Re-check with fresh seeds and the slow solver. Returns dict of counts."""
    a, q = np.array(a, float), np.array(q, float)
    out = {}
    for s in (1, 2):
        out[f"count_fast_seed{s}"] = list(count_fast(a, q, seed=s))
    s, pts, idx = summary(a, q, nseed=2500, rounds=3, seed=7)
    out["summary"] = {k: s[k] for k in ("n", "m1", "m2", "other", "morse_ok", "min_eig")}
    return out


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    gens = [gen_tri_interior, gen_c3v_axial, gen_symmetric, gen_two_pair,
            gen_near_coplanar, gen_insertion, gen_random]
    tasks = []
    for g in gens:
        t = g(rng)
        print(f"{g.__name__}: {len(t)} configs", flush=True)
        tasks += t
    print(f"total sweep configs: {len(tasks)}", flush=True)

    with Pool(8) as pool:
        results = pool.map(screen, tasks, chunksize=32)
    print(f"sweep done in {time.time()-t0:.0f}s", flush=True)

    dist, inconsistent, best = Counter(), Counter(), []
    for fam, a, q, n, m1, m2, ok, err in results:
        if ok:
            dist[(fam.split('_')[0] if False else fam, n)] += 1
            if n >= 7:
                best.append(dict(family=fam, n=int(n), m1=int(m1), m2=int(m2), a=a, q=q))
        else:
            inconsistent[fam] += 1
    overall = Counter()
    for (fam, n), c in dist.items():
        overall[n] += c
    max_n = max(overall) if overall else -1
    print("overall Morse-consistent distribution:", dict(sorted(overall.items())), flush=True)
    print("Morse-inconsistent (solver misses, not counted):", dict(inconsistent), flush=True)

    # ---- hill climbing from certified C4 and the best sweep configs ----
    C4a = [[3, 0, 0], [0, 3, 0], [0, 0, 3], [1, 1, 1]]
    C4q = [1.0, 1.0, 1.0, 2 ** -10]
    seeds_hc = [(C4a, C4q)]
    nines = [b for b in best if b["n"] >= 9][:5]
    sevens = [b for b in best if b["n"] == 7][:3]
    for b in nines + sevens:
        seeds_hc.append((b["a"], b["q"]))
    steps = 150 if QUICK else 500
    args = [(i, a, q, steps) for i, (a, q) in enumerate(seeds_hc)]
    with Pool(min(8, len(args))) as pool:
        chains = pool.map(hill_chain, args)
    hc_evals = sum(c["evals"] for c in chains)
    hc_best = max(c["best_n"] for c in chains)
    print(f"hill-climb: {len(chains)} chains x {steps} steps, {hc_evals} evals, best {hc_best}", flush=True)
    for c in chains:
        if c["best_n"] >= 9:
            best.append(dict(family="hill_climb", n=c["best_n"], m1=None, m2=None, a=c["a"], q=c["q"]))
        for h in c["improvements"]:
            if h["n"] > 9:
                best.append(dict(family="hill_climb", n=h["n"], m1=None, m2=None, a=h["a"], q=h["q"]))

    # ---- verify every candidate screening >=9 (dedup by rounded coords) ----
    seen, verified = set(), []
    for b in sorted(best, key=lambda b: -b["n"]):
        if b["n"] < 9:
            continue
        key = (np.round(np.array(b["a"]) * 200).tobytes(), np.round(np.log10(np.array(b["q"])) * 20).tobytes())
        if key in seen:
            continue
        seen.add(key)
        v = verify_candidate(b["a"], b["q"])
        verified.append({**b, "verification": v})
        print(f"verify screen-n={b['n']} fam={b['family']}: {v}", flush=True)
        if len(verified) >= 40:
            break

    over9 = [v for v in verified
             if v["verification"]["summary"]["n"] > 9 and v["verification"]["summary"]["morse_ok"]]
    report = dict(
        seed=SEED, quick=QUICK, seconds=round(time.time() - t0, 1),
        sweep_configs=len(tasks), hill_climb_evals=int(hc_evals),
        overall_distribution={str(k): v for k, v in sorted(overall.items())},
        per_family_distribution={f"{f}|n={n}": c for (f, n), c in sorted(dist.items())},
        morse_inconsistent=dict(inconsistent),
        max_screened=int(max_n), hill_climb_best=int(hc_best),
        n_ge7=len(best), verified_candidates=verified,
        over9_verified=over9,
        hill_chains=[{k: c[k] for k in ("cid", "best_n", "evals")} for c in chains],
        best_ge7_sample=sorted(best, key=lambda b: -b["n"])[:25],
    )
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    jpath = os.path.join(root, "results", "wave1_maxwell_n4.json")
    json.dump(report, open(jpath, "w"), indent=1, default=str)

    fam_rows = {}
    for (f, n), c in dist.items():
        fam_rows.setdefault(f, Counter())[n] += c
    lines = [
        "# Wave-1 search for >9 equilibria with four charges (issue #3)",
        "",
        f"Script: `explore/wave1_maxwell_n4.py` (seed {SEED}). Screening: `charges.count_fast`;",
        "only Morse-consistent counts (m2-m1=3, no other indices) are tallied as counts.",
        f"{len(verified)} de-duplicated configurations screening >=9 were re-checked with fresh `count_fast` seeds and `charges.summary`; slow-solver counts: "
        + str(dict(sorted(Counter(v["verification"]["summary"]["n"] for v in verified).items()))) + ".",
        "Morse-inconsistent screens include degenerate tri_interior grid points where the 4th site coincides with a triangle vertex; they are excluded from all counts.",
        f"Runtime {report['seconds']}s. Sweep configurations: {len(tasks)}. Hill-climb evaluations: {hc_evals}.",
        "",
        f"**Maximum Morse-consistent count found: {max_n}.** "
        f"Configurations verified >9 by the slow solver: {len(over9)}.",
        "",
        "## Count distribution (Morse-consistent screens only)",
        "",
        "| family | " + " | ".join(f"n={k}" for k in sorted(overall)) + " | inconsistent |",
        "|---|" + "---|" * (len(overall) + 1),
    ]
    for f in sorted(fam_rows):
        lines.append(f"| {f} | " + " | ".join(str(fam_rows[f].get(k, 0)) for k in sorted(overall))
                     + f" | {inconsistent.get(f, 0)} |")
    lines.append(f"| **total** | " + " | ".join(str(overall[k]) for k in sorted(overall))
                 + f" | {sum(inconsistent.values())} |")
    lines += [
        "",
        f"Hill climbing ({len(chains)} chains x {steps} steps from the certified 9-equilibrium",
        f"configuration and the best sweep finds): best count {hc_best}.",
        "",
        "## Claim labels",
        "",
        "- VERIFIED (finite numerical check): the counts above, for the configurations actually run.",
        "- NOT a proof: no covering argument over the 8-parameter (shape+charge) space is made.",
        "  Absence of >9 here is evidence about the searched families, nothing more.",
        "",
        "Best configurations and all verified candidates: `results/wave1_maxwell_n4.json`.",
    ]
    open(os.path.join(root, "results", "wave1_maxwell_n4.md"), "w").write("\n".join(lines) + "\n")
    print(f"max={max_n} over9_verified={len(over9)} wrote results/wave1_maxwell_n4.{{json,md}}", flush=True)


if __name__ == "__main__":
    main()
