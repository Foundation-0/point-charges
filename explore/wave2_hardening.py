"""Wave-2 Task A: re-examine the 332 Morse-inconsistent wave-1 screens.

Wave 1 (`explore/wave1_maxwell_n4.py`) tallied 332 Morse-inconsistent
screens (327 tri_interior, 5 symmetric_tet) but stored only the tallies.
This script regenerates the exact deterministic task lists for those two
families (same seed, same generator code, imported from wave 1), re-screens
with the same `count_fast(seed=0)`, and classifies every inconsistent
screen as:

  degenerate-coincident : two sites coincide (min pairwise distance ~ 0).
      The configuration is then really a 3-charge configuration with the
      coincident charges merged, so m2-m1=2, not 3: Morse-inconsistency
      is definitional, not a solver miss. Verified by re-counting the
      merged 3-charge configuration with count_fast (expect ok for n=3).
  solver artifact       : sites distinct; higher-effort recounting gives
      the true count. `charges.summary` was attempted at two efforts and
      hangs on these cases (see clustered_count docstring), so a
      clustered heavy multi-start counter is used instead, with a
      certify.py attempt on top. Subtype solver_artifact_degenerate:
      the true equilibria include a degenerate one, which is what
      breaks both repo solvers and the non-degenerate Morse check.
  genuine candidate     : slow solver gives a Morse-consistent count > 9,
      or no Morse-consistent count is obtained at high effort. Any such
      case is reported for certification (make_certs machinery).

Analytic prediction for tri_interior (checked against the rerun):
the barycentric grid has 3 vertex positions out of 15; the 4th site
coincides with a vertex exactly when the height is 0 there.
Main grid: 9 triangles x 3 vertices x 1 zero-height x 11 charges = 297.
Unequal-charge part: 5 charge triples x 1 vertex in bary[::3] x 6 charges
= 30. Total 327 -- the full wave-1 tri_interior tally.

Outputs: results/wave2_hardening.json, results/wave2_hardening.md
Usage: .venv/bin/python explore/wave2_hardening.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import json, time
from collections import Counter
from multiprocessing import Pool
import numpy as np
from charges import count_fast, find_fast, grad_hess, newton, newton_local
import wave1_maxwell_n4 as w1

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JPATH = os.path.join(ROOT, "results", "wave2_hardening.json")
MPATH = os.path.join(ROOT, "results", "wave2_hardening.md")
COINC_TOL = 1e-9  # relative to configuration scale


def min_pair_dist(a):
    a = np.asarray(a, float)
    d = np.linalg.norm(a[:, None] - a[None], axis=2)
    d += np.diag(np.full(len(a), np.inf))
    sc = max(np.max(np.ptp(a, 0)), 1e-300)
    return float(d.min() / sc), float(d.min())


def merged(a, q):
    """Merge coincident sites, summing their charges."""
    a = np.asarray(a, float); q = np.asarray(q, float)
    sc = max(np.max(np.ptp(a, 0)), 1e-300)
    keep_a, keep_q = [], []
    for i in range(len(a)):
        hit = None
        for j, b in enumerate(keep_a):
            if np.linalg.norm(a[i] - b) < COINC_TOL * sc:
                hit = j; break
        if hit is None:
            keep_a.append(a[i]); keep_q.append(q[i])
        else:
            keep_q[hit] += q[i]
    return np.array(keep_a), np.array(keep_q)


def classify_coincident(task):
    fam, a, q, n, m1, m2 = task
    ma, mq = merged(a, q)
    mn, mm1, mm2, mok = count_fast(ma, mq, seed=0)
    return dict(family=fam, a=a, q=q, screen=[n, m1, m2],
                merged_nsites=int(len(mq)), merged_a=ma.tolist(), merged_q=mq.tolist(),
                merged_count=[int(mn), int(mm1), int(mm2), bool(mok)])


def clustered_count(a, q, seed=0, nrand=20000):
    """Higher-effort count robust to the failure mode found here.

    charges.find_fast dedups with tol 1e-6*dist-to-nearest-site; near a
    degenerate equilibrium Newton converges so slowly that hundreds of
    iterates pass the gradient test at mutual distances > tol and are
    counted as distinct equilibria (467 'points' within 1.5e-3 of the
    tetrahedron centre). charges.summary / find_critical_points then
    blows up in its O(k^2) dedup on the same cloud (both efforts tried:
    (nseed=8000,rounds=4) >7 min, (2500,3) >4 min, no result, killed).
    Here: heavy multi-start (global Newton), strict gradient test,
    greedy clustering ordered by |grad| at tol 3e-3*scale, polish,
    re-cluster. Returns points, indices, eigenvalues, |grad|."""
    a = np.asarray(a, float); q = np.asarray(q, float)
    n = len(a); sc = np.max(np.ptp(a, 0))
    rng = np.random.default_rng(seed)
    S = [rng.dirichlet(np.full(n, 0.4), size=nrand) @ a]
    for i in range(n):
        for j in range(i + 1, n):
            t = np.linspace(0.02, 0.98, 49)
            S.append(a[i] + t[:, None] * (a[j] - a[i]))
        others = [k for k in range(n) if k != i]
        if len(others) >= 3:
            w = rng.dirichlet(np.ones(len(others)), size=2000)
            S.append(w @ a[others] + 1e-3 * sc * rng.standard_normal((2000, 3)))
    def accept(X):
        g, H = grad_hess(X, a, q)
        r2 = ((X[:, None] - a[None]) ** 2).sum(2)
        fs = (q[None] / r2).sum(1)
        gn = np.linalg.norm(g, axis=1)
        ok = np.isfinite(X).all(1) & (gn < 1e-9 * fs) & (np.sqrt(r2).min(1) > 1e-8 * sc)
        return X[ok], gn[ok]

    def cluster(X, gn, tol):
        reps = []
        for k in np.argsort(gn):
            if all(np.linalg.norm(X[k] - r) > tol for r in reps):
                reps.append(X[k])
        return np.array(reps).reshape(-1, 3)

    X = newton(np.vstack(S), a, q, iters=100)
    X, gn = accept(X)
    X = cluster(X, gn, 3e-3 * sc)
    if len(X):
        X = newton_local(X, a, q, iters=200)
    X, gn = accept(X)
    # near a degenerate equilibrium polished reps still crawl at mutual
    # distances ~1e-3*scale; primary count clusters at 1e-2*scale and the
    # finer count is reported for sensitivity (equilibria closer than
    # 1e-2*scale would be merged -- stated limitation).
    n_fine = len(cluster(X, gn, 3e-3 * sc))
    X = cluster(X, gn, 1e-2 * sc)
    g, H = grad_hess(X, a, q) if len(X) else (np.zeros((0, 3)), np.zeros((0, 3, 3)))
    ev = np.linalg.eigvalsh(H) if len(X) else np.zeros((0, 3))
    if len(X):
        r2 = ((X[:, None] - a[None]) ** 2).sum(2)
        hscale = (q[None] / (r2 * np.sqrt(r2))).sum(1)  # typical |Hessian| scale
    else:
        hscale = np.zeros(0)
    return X, (ev < 0).sum(1), ev, np.linalg.norm(g, axis=1), n_fine, hscale


def classify_distinct(task):
    fam, a, q, n, m1, m2, mult = task
    a_a, q_a = np.asarray(a, float), np.asarray(q, float)
    out = dict(family=fam, a=a, q=q, screen=[n, m1, m2], multiplicity=mult, slow={})
    for s in (1, 2, 3):
        out["slow"][f"count_fast_seed{s}"] = [int(x) if not isinstance(x, bool) else bool(x)
                                              for x in count_fast(a_a, q_a, seed=s)]
    # find_fast cloud diagnostic: are the 'points' one tight cloud?
    pts, _, _ = find_fast(a_a, q_a, seed=0)
    if len(pts):
        cen = pts.mean(0)
        out["slow"]["find_fast_cloud"] = dict(
            n=int(len(pts)), max_dist_to_centroid=float(np.linalg.norm(pts - cen, axis=1).max()))
    runs = []
    for sd in (0, 5):
        P, idx, ev, gn, n_fine, hs = clustered_count(a_a, q_a, seed=sd)
        # degenerate: eigenvalues far below the local Hessian scale sum q/r^3
        # (reps near a degenerate point stall at ~1e-3 from it, where the
        # cubic-order Hessian is still ~1e-3 of scale)
        deg = [bool(np.abs(e).min() < 1e-3 * h or np.abs(e).max() < 1e-2 * h)
               for e, h in zip(ev, hs)]
        runs.append(dict(seed=sd, n=int(len(P)), n_at_fine_tol=int(n_fine), points=P.tolist(),
                         indices=[int(i) for i in idx], degenerate=deg,
                         min_abs_eig=[float(np.abs(e).min()) for e in ev]))
    out["slow"]["clustered"] = runs
    r0 = runs[0]
    agree = runs[0]["n"] == runs[1]["n"]
    ndeg = sum(r0["degenerate"])
    m1c = sum(1 for i, d in zip(r0["indices"], r0["degenerate"]) if i == 1 and not d)
    m2c = sum(1 for i, d in zip(r0["indices"], r0["degenerate"]) if i == 2 and not d)
    # perturbation check: breaking symmetry should restore Morse consistency
    rngp = np.random.default_rng(42)
    pert = []
    for k in range(3):
        qp = q_a * np.exp(1e-2 * rngp.standard_normal(len(q_a)))
        pert.append([int(x) if not isinstance(x, bool) else bool(x)
                     for x in count_fast(a_a, qp, seed=k)])
    out["slow"]["charge_perturbation_1e-2"] = pert
    if agree and ndeg == 0 and r0["n"] > 9 and m2c - m1c == len(q_a) - 1:
        out["verdict"] = "genuine_candidate_over9"; out["true_count"] = r0["n"]
    elif agree and ndeg > 0:
        out["verdict"] = "solver_artifact_degenerate"
        out["true_count"] = r0["n"]; out["degenerate_count"] = ndeg
    elif agree and m2c - m1c == len(q_a) - 1 and r0["n"] == m1c + m2c:
        out["verdict"] = "solver_artifact"; out["true_count"] = r0["n"]
    else:
        out["verdict"] = "unresolved"; out["true_count"] = r0["n"]
    # certification attempt: required if >9; also informative here
    # (Krawczyk cannot verify a degenerate zero, so failure is expected
    # and recorded, not hidden).
    try:
        from certify import run as cert_run
        t_c = time.time()
        res = cert_run(a_a, q_a, "wave2-distinct-" + fam, out=None)
        out["certify"] = dict(count=res["count"], index_counts=res["index_counts"],
                               seconds=round(time.time() - t_c, 1))
    except Exception as e:
        out["certify"] = dict(error=f"{type(e).__name__}: {e}"[:160])
    return out


def main():
    t0 = time.time()
    rng = np.random.default_rng(w1.SEED)
    # gen_tri_interior is the first generator wave 1 calls and the only one
    # of the two here that draws from rng; gen_symmetric draws nothing.
    tasks = w1.gen_tri_interior(rng) + w1.gen_symmetric(rng)
    tasks = [t for t in tasks if t[0] in ("tri_interior", "symmetric_tet")]
    print(f"regenerated tasks: {len(tasks)} "
          f"({Counter(t[0] for t in tasks)})", flush=True)

    with Pool(8) as pool:
        results = pool.map(w1.screen, tasks, chunksize=32)
    inconsistent = [r for r in results if not r[6]]
    by_fam = Counter(r[0] for r in inconsistent)
    print(f"re-screen done in {time.time()-t0:.0f}s; inconsistent: "
          f"{len(inconsistent)} {dict(by_fam)}", flush=True)

    coincident, distinct = [], []
    for fam, a, q, n, m1, m2, ok, err in inconsistent:
        rel, _ = min_pair_dist(a)
        (coincident if rel < COINC_TOL else distinct).append((fam, a, q, n, m1, m2))
    # also: did any coincident-site task screen as consistent? (it cannot:
    # merged n=3 has m2-m1=2 != 3, but check the data, do not assume)
    coinc_consistent = sum(1 for r in results
                           if r[6] and min_pair_dist(r[1])[0] < COINC_TOL)
    coinc_total = sum(1 for t in tasks if min_pair_dist(t[1])[0] < COINC_TOL)
    print(f"coincident: {len(coincident)} inconsistent, "
          f"{coinc_consistent} consistent, {coinc_total} total; "
          f"distinct inconsistent: {len(distinct)}", flush=True)

    checkpoint = dict(stage="rescreened", seconds=round(time.time() - t0, 1),
                      regenerated_tasks=len(tasks),
                      inconsistent_total=len(inconsistent),
                      inconsistent_by_family=dict(by_fam),
                      coincident_inconsistent=len(coincident),
                      coincident_consistent=coinc_consistent,
                      coincident_total=coinc_total,
                      distinct_inconsistent=len(distinct),
                      distinct_raw=distinct, coincident_raw=coincident)
    json.dump(checkpoint, open(JPATH, "w"), indent=1, default=str)  # checkpoint

    with Pool(8) as pool:
        coinc_done = pool.map(classify_coincident, coincident, chunksize=16)
    print(f"coincident classified in {time.time()-t0:.0f}s", flush=True)
    checkpoint["stage"] = "coincident_done"
    checkpoint["coincident"] = coinc_done
    json.dump(checkpoint, open(JPATH, "w"), indent=1, default=str)

    # dedup exact duplicate configurations (wave-1 gen_symmetric emits the
    # all-ones charge pattern 5 times: the initial pattern plus r=1 in
    # 10**linspace(-4,1,11) for each of the 4 pattern shapes)
    uniq = {}
    for t in distinct:
        key = (t[0], np.asarray(t[1]).tobytes(), np.asarray(t[2]).tobytes())
        if key in uniq:
            uniq[key][-1] += 1
        else:
            uniq[key] = [*t, 1]
    distinct_uniq = list(uniq.values())
    print(f"distinct unique configs: {len(distinct_uniq)} "
          f"(multiplicities {[t[-1] for t in distinct_uniq]})", flush=True)
    with Pool(min(8, max(1, len(distinct_uniq)))) as pool:
        distinct_done = pool.map(classify_distinct, distinct_uniq) if distinct_uniq else []
    print(f"distinct classified in {time.time()-t0:.0f}s", flush=True)

    merged_ok = sum(1 for c in coinc_done if c["merged_count"][3])
    merged_counts = Counter(c["merged_count"][0] for c in coinc_done if c["merged_count"][3])
    verdicts = Counter()
    for d in distinct_done:
        verdicts[d["verdict"]] += d["multiplicity"]
    true_counts = Counter(d.get("true_count") for d in distinct_done
                          if d["verdict"].startswith("solver_artifact"))
    over9 = [d for d in distinct_done if d["verdict"] == "genuine_candidate_over9"]
    unresolved = [d for d in distinct_done if d["verdict"] == "unresolved"]

    report = {**checkpoint, "stage": "done", "seconds": round(time.time() - t0, 1),
              "coincident": coinc_done, "distinct": distinct_done,
              "merged_morse_ok": merged_ok,
              "merged_count_distribution": {str(k): v for k, v in sorted(merged_counts.items())},
              "distinct_verdicts": dict(verdicts),
              "distinct_true_counts": {str(k): v for k, v in sorted(true_counts.items())},
              "over9_candidates": over9}
    json.dump(report, open(JPATH, "w"), indent=1, default=str)

    lines = [
        "# Wave-2 Task A: the 332 Morse-inconsistent wave-1 screens, re-examined",
        "",
        "Script: `explore/wave2_hardening.py`. Wave 1 stored only tallies for its",
        "332 Morse-inconsistent screens, so the exact task lists of the two affected",
        "families were regenerated from the wave-1 generators (same seed 20261018,",
        "same code) and re-screened with `count_fast(seed=0)`.",
        "",
        f"- Regenerated tasks: {len(tasks)}. Inconsistent on re-screen: "
        f"{len(inconsistent)} ({dict(by_fam)}). Wave-1 tally: 332 "
        "({'tri_interior': 327, 'symmetric_tet': 5}).",
        f"- Runtime {report['seconds']}s.",
        "",
        "## Classification",
        "",
        f"### Degenerate-coincident: {len(coincident)}",
        "",
        "Two sites coincide (relative min pairwise distance < 1e-9). The analytic",
        "count for tri_interior is 297 (main grid: 9 triangles x 3 barycentric",
        "vertices x height 0 x 11 fourth charges) + 30 (unequal-charge part:",
        "5 triples x 1 vertex in bary[::3] x 6 charges) = 327.",
        f"Coincident tasks total {coinc_total}; of those {len(coincident)} screened",
        f"inconsistent and {coinc_consistent} screened consistent. A coincident pair",
        "is one charge with the summed value, i.e. an n=3 configuration with",
        "m2-m1=2; wave 1 required m2-m1=3 (n=4), so inconsistency is definitional.",
        f"Merging each coincident pair and re-counting with `count_fast`: "
        f"{merged_ok}/{len(coinc_done)} Morse-consistent for n=3, count distribution "
        f"{dict(sorted(merged_counts.items()))}.",
        "",
        f"### Distinct sites: {len(distinct)} screens, "
        f"{len(distinct_done)} unique configuration(s)",
        "",
        "Wave-1's gen_symmetric emits the all-ones charge pattern 5 times (initial",
        "pattern + r=1 from 10**linspace(-4,1,11) in each of 4 pattern shapes), so",
        "duplicate screens are classified once with their multiplicity.",
        "",
        "The slow solver `charges.summary` could not be used on these cases: both",
        "(nseed=8000, rounds=4) and (nseed=2500, rounds=3) ran for minutes with no",
        "result (killed; see failed attempts). The cause is the same cloud that",
        "breaks the screen: near a degenerate equilibrium Newton converges slowly,",
        "hundreds of iterates pass the gradient test, and the O(k^2) dedup in",
        "find_critical_points / the 1e-6-relative dedup in find_fast turn one",
        "equilibrium into hundreds of 'points'. Re-examination therefore used",
        "`count_fast` seeds 1-3, a cloud diagnostic on find_fast's points, and a",
        "clustered heavy multi-start counter (`clustered_count` in this script:",
        "28k seeds, strict gradient test, clustering at 1e-2*scale (3e-3 reported), polish, two",
        "seeds), plus 1e-2 charge perturbations and a certify.py attempt.",
        f"Verdicts (weighted by multiplicity): {dict(verdicts)}. True counts of "
        f"solver artifacts (unique configs): {dict(sorted(true_counts.items()))}.",
        "",
        "For the regular tetrahedron with equal charges the centre is an equilibrium",
        "with Hessian exactly 0: the tetrahedral symmetry group acts irreducibly on",
        "R^3, so an invariant symmetric Hessian is a scalar matrix, and harmonicity",
        "(trace 0) forces the scalar to 0. It is therefore degenerate, Maxwell's",
        "non-degenerate Morse bookkeeping does not apply, and certify.py's Krawczyk",
        "test (which proves non-singularity) cannot verify it. Perturbing the",
        "charges by 1e-2 restores Morse-consistent counts (see JSON), and the",
        "stretched tetrahedra (s != 1) in wave 1 screened consistent.",
        "",
    ]
    for d in distinct_done:
        cl = d["slow"].get("clustered", [{}])
        lines.append(f"- {d['family']} x{d['multiplicity']} q={['%.4g' % x for x in d['q']]} "
                     f"a={d['a']} screen(n,m1,m2)={d['screen']} -> {d['verdict']}"
                     + (f", true count {d['true_count']}" if 'true_count' in d else "")
                     + (f" ({d['degenerate_count']} degenerate)" if d.get('degenerate_count') else "")
                     + f"; clustered n={[r.get('n') for r in cl]}"
                     + f"; cloud={d['slow'].get('find_fast_cloud')}"
                     + f"; certify={d.get('certify')}")
    lines += [
        "",
        f"**Genuine candidates with count > 9: {len(over9)}. Unresolved: {len(unresolved)}.**",
        "certify.py was attempted on each unique distinct config (results above);",
        "it certifies non-degenerate zeros only, so it cannot settle the degenerate",
        "tetrahedron centre. No >9 candidate exists to certify.",
        "",
        "## Claim labels",
        "",
        "- VERIFIED (finite numerical check): the classifications and true counts above,",
        "  for the regenerated configurations actually run. Regeneration is exact only",
        "  if the re-screen tally reproduces 332 / {327, 5}; see the numbers above.",
        "- The coincident classification is additionally PROVEN by the grid arithmetic",
        "  (297+30=327) together with the merged recount: coincident sites are an n=3",
        "  configuration, for which m2-m1=2 by Morse theory (NOTE.md S1).",
        "- PROVEN: the regular-tetrahedron centre is a degenerate equilibrium",
        "  (Hessian 0 by irreducible tetrahedral symmetry + trace 0).",
        "- VERIFIED (finite): its clustered true count 1 (two seeds, 28k starts;",
        "  ray scans and perturbations in exploration agreed). NOT certified by",
        "  certify.py/verify_arb.py: both require non-degenerate zeros.",
        "- NOT a proof about n=4 in general: nothing here covers configurations wave 1",
        "  did not generate.",
        "",
        "## Failed attempts",
        "",
        "- `charges.summary(nseed=8000, rounds=4)` on the distinct case: no result",
        "  in >7 min per config (workers at 98% CPU); killed. `summary(2500, 3)`:",
        "  also no result in >4 min; killed. The bottleneck is the O(k^2) dedup in",
        "  find_critical_points over the cloud of Newton iterates around the",
        "  degenerate centre. Replaced by `clustered_count`.",
        "- certify.py on the regular tetrahedron: RuntimeError 'box budget exceeded'",
        "  (expected: its Krawczyk test proves Hessian non-singularity, false at the",
        "  degenerate centre). verify_arb.py not run: no certificate exists to check.",
        "",
        "Full data: `results/wave2_hardening.json`.",
    ]
    open(MPATH, "w").write("\n".join(lines) + "\n")
    print(f"done: coincident={len(coincident)} distinct={len(distinct)} "
          f"over9={len(over9)} unresolved={len(unresolved)} "
          f"wrote results/wave2_hardening.{{json,md}}", flush=True)


if __name__ == "__main__":
    main()
