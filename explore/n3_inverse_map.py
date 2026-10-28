"""n=3 inverse-map scaffolding for a computer-assisted proof (issue #4).

Reformulation (NOTE.md S3): for a fixed triangle A1A2A3 and interior point
x, Lami's theorem gives the unique charge ratio Z(x) making x an
equilibrium: zeta_i = r_i^2 w_i, (w1:w2:w3) = (u2 x u3 : u3 x u1 : u1 x u2),
u_i = (x-A_i)/r_i. With D(x) the (normalised) determinant of the in-plane
Hessian at x for charges Z(x), R_min = {D>0} is the set of x that are
local minima for their own charge ratio, and

    Maxwell bound 4 for this triangle  <=>  Z is one-to-one on R_min.

This module implements, for a grid of triangle shapes:
  1. the inverse map Z and D (float64, via explore/three_map.zeta_D);
  2. rasterised coverage of Z(R_min) on the charge simplex (as in
     explore/three_scan.py) -- a finite numerical preimage count;
  3. forward spot-checks: sampled charge ratios in Z(R_min) are re-solved
     directly (explore/three.crit2d) and their minima counted;
  4. interval-arithmetic enclosures of D and Z over boxes in barycentric
     coordinates, using the repo's ia.py (outward-rounded float64
     intervals; ia.py is checked against exact rationals by
     tests/test_ia.py). Boxes are certified D>0 / D<0 where the enclosure
     decides the sign; straddling boxes are subdivided to a depth cap and
     the leftover is reported, not hidden. On certified D>0 boxes the
     Z-enclosures give an overlap multiplicity at sampled charge points.

What this is: VERIFIED finite checks (items 2-3) plus PROVEN sign facts on
the certified boxes of item 4 (modulo ia.py, which is independently
tested). What it is not: a proof of injectivity on any shape, let alone
all shapes. results/n3_inverse_map.md states the remaining gaps exactly.

Usage: .venv/bin/python explore/n3_inverse_map.py [--quick]
Outputs: results/n3_inverse_map.json, results/n3_inverse_map.md
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
import json, time
from collections import Counter
from multiprocessing import Pool
import numpy as np
import ia
from three_map import zeta_D
from three_scan import coverage
from three import crit2d

QUICK = "--quick" in sys.argv
JPATH = os.path.join(ROOT, "results", "n3_inverse_map.json")
MPATH = os.path.join(ROOT, "results", "n3_inverse_map.md")


def tri(a, b):
    return np.array([[0.0, 0.0], [1.0, 0.0], [a, b]])


# ---------------------------------------------------------------- shapes
def shape_grid():
    """Fundamental domain: side A1A2 = 1 longest, A3=(a,b), 0<=a<=1/2, b>0."""
    shapes = []
    bs = [0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.866]
    for b in bs:
        for a in (0.0, 0.15, 0.3, 0.5):
            if (a - 1) ** 2 + b ** 2 <= 1 + 1e-12 and a * a + b * b <= 1 + 1e-12:
                shapes.append((a, b))
    return shapes


def coverage_job(ab):
    a, b = ab
    A = tri(a, b)
    mc, npx, sbox, _ = coverage(A, N=150, P=150)
    if sbox is None:
        return dict(a=a, b=b, max_coverage=0, pixels=0, resolved=False)
    mc2, npx2, _, _ = coverage(A, N=300, P=250, box=sbox)
    # resolved = R_min image actually rasterised at zoom resolution;
    # flat shapes (b=0.02) have D>0 grid points but a sub-pixel image.
    return dict(a=a, b=b, max_coverage=int(mc2), pixels=int(npx2),
                resolved=bool(npx2 > 0))


# ------------------------------------------------- float inverse map grid
def inverse_grid(A, N=240):
    """Z and normalised D on a Duffy grid; returns interior points only."""
    s = 0.5 - 0.5 * np.cos(np.pi * np.linspace(0, 1, N + 1))
    S, T = np.meshgrid(s, s, indexing="ij")
    l2 = S * (1 - T); l3 = S * T; l1 = 1 - S
    X = (l1[..., None] * A[0] + l2[..., None] * A[1] + l3[..., None] * A[2]).reshape(-1, 2)
    inside = ((l1 > 1e-12) & (l2 > 1e-12) & (l3 > 1e-12)).ravel()
    z, D = zeta_D(A, X[inside])
    return X[inside], z, D


def jac_det_sign(A, X, z, D, h=1e-6):
    """det DZ (x -> (zeta1,zeta2)) by central differences, on D>0 points.
    NOTE.md predicts Z preserves orientation on R_min: det DZ > 0 there."""
    sel = np.where(D > 0)[0]
    if len(sel) > 400:
        sel = sel[np.linspace(0, len(sel) - 1, 400).astype(int)]
    dets = []
    for k in sel:
        x = X[k]
        cols = []
        bad = False
        for e in ((h, 0.0), (0.0, h)):
            zp, _ = zeta_D(A, (x + np.array(e))[None])
            zm, _ = zeta_D(A, (x - np.array(e))[None])
            if not (np.isfinite(zp).all() and np.isfinite(zm).all()):
                bad = True; break
            cols.append((zp[0, :2] - zm[0, :2]) / (2 * h))
        if not bad:
            dets.append(cols[0][0] * cols[1][1] - cols[0][1] * cols[1][0])
    dets = np.array(dets)
    return dict(n=int(len(dets)), n_positive=int((dets > 0).sum()),
                n_negative=int((dets < 0).sum()),
                min_abs=float(np.abs(dets).min()) if len(dets) else None)


def forward_checks(A, z, D, kmax=30, seed=0):
    """Re-solve sampled charge ratios from R_min directly; count minima."""
    rng = np.random.default_rng(seed)
    sel = np.where(D > 0)[0]
    if len(sel) == 0:
        return dict(n=0, m0_counts={}, morse_ok=0)
    pick = sel[rng.choice(len(sel), size=min(kmax, len(sel)), replace=False)]
    counts, ok = Counter(), 0
    for k in pick:
        m0s = []
        for s in (1, 2):
            pts, m0, ms = crit2d(A, z[k].copy(), seed=s)
            m0s.append(m0)
            ok += int(ms - m0 == 2)
        counts[int(max(m0s))] += 1
    return dict(n=int(len(pick)), m0_counts={str(k): v for k, v in sorted(counts.items())},
                morse_ok=int(ok), morse_total=int(2 * len(pick)))


# ------------------------------------------------------ interval machinery
def _cross(u, v):
    return ia.sub(ia.mul(u[0], v[1]), ia.mul(u[1], v[0]))


def interval_zeta_D(A, L2, L3):
    """Interval Z (3 x box) and normalised D over barycentric boxes.

    L2, L3 are ia intervals (lo,hi) arrays of shape (m,); x = A1 + l2*(A2-A1)
    + l3*(A3-A1). Returns (Z, Dn) as ia intervals; NaN = no information.
    Mirrors three_map.zeta_D operation by operation."""
    e2 = A[1] - A[0]; e3 = A[2] - A[0]
    X = []
    for c in range(2):
        xc = ia.add(ia.const(np.full(L2[0].shape, A[0][c])),
                    ia.add(ia.scale(e2[c], L2), ia.scale(e3[c], L3)))
        X.append(xc)
    d, r2, r, u = [], [], [], []
    for i in range(3):
        di = tuple(ia.sub(X[c], ia.const(np.full(L2[0].shape, A[i][c]))) for c in range(2))
        r2i = ia.add(ia.sqr(di[0]), ia.sqr(di[1]))
        ri = ia.sqrt(r2i)
        d.append(di); r2.append(r2i); r.append(ri)
        u.append(tuple(ia.div(di[c], ri) for c in range(2)))
    w = [_cross(u[1], u[2]), _cross(u[2], u[0]), _cross(u[0], u[1])]
    wsum = ia.add(ia.add(w[0], w[1]), w[2])
    # zeta_D flips sign so wsum>0; where the enclosure straddles 0 we keep
    # both signs by widening to the hull of w and -w (sound, may straddle).
    flip = wsum[1] < 0
    strad = (wsum[0] <= 0) & (wsum[1] >= 0)
    wx = []
    for wi in w:
        lo = np.where(flip, -wi[1], wi[0]); hi = np.where(flip, -wi[0], wi[1])
        lo = np.where(strad, np.minimum(wi[0], -wi[1]), lo)
        hi = np.where(strad, np.maximum(wi[1], -wi[0]), hi)
        wx.append((lo, hi))
    zun = [ia.mul(wx[i], r2[i]) for i in range(3)]
    zsum = ia.add(ia.add(zun[0], zun[1]), zun[2])
    Z = [ia.div(zun[i], zsum) for i in range(3)]
    H = [[ia.const(np.zeros(L2[0].shape)) for _ in range(2)] for _ in range(2)]
    i3sum = ia.const(np.zeros(L2[0].shape))
    Ssum = ia.const(np.zeros(L2[0].shape))
    for i in range(3):
        r3 = ia.mul(r2[i], r[i]); r5 = ia.mul(ia.mul(r2[i], r2[i]), r[i])
        i3 = ia.div(Z[i], r3); i5 = ia.div(Z[i], r5)
        i3sum = ia.add(i3sum, i3)
        Ssum = ia.add(Ssum, ia.div(Z[i], r2[i]))
        for a_ in range(2):
            for b_ in range(2):
                H[a_][b_] = ia.add(H[a_][b_], ia.scale(3.0, ia.mul(i5, ia.mul(d[i][a_], d[i][b_]))))
    H[0][0] = ia.sub(H[0][0], i3sum); H[1][1] = ia.sub(H[1][1], i3sum)
    D = ia.sub(ia.mul(H[0][0], H[1][1]), ia.sqr(H[0][1]))
    Dn = ia.div(D, ia.sqr(Ssum))
    return Z, Dn


def interval_analysis(A, delta=0.02, n0=20, max_depth=6, sample_pts=None,
                      max_boxes=400_000):
    """Certify the sign of D on boxes covering the delta-shrunk triangle
    {l1,l2,l3 >= delta}. Returns certified area fractions (barycentric
    measure) and, on D>0 boxes, the max number of Z-enclosures covering
    any sampled charge point."""
    # boxes: (l2lo,l2hi,l3lo,l3hi), kept only if all corners satisfy
    # l2>=d, l3>=d, l2+l3<=1-d (box inside shrunk triangle)
    edges = np.linspace(delta, 1 - 2 * delta, n0 + 1)
    boxes = []
    for i in range(n0):
        for j in range(n0):
            if edges[i + 1] + edges[j + 1] <= 1 - delta + 1e-15:
                boxes.append((edges[i], edges[i + 1], edges[j], edges[j + 1]))
    pos_boxes, neg_area, total_area = [], 0.0, 0.0
    unresolved_area, dropped_area = 0.0, 0.0
    depth = 0
    while boxes and depth <= max_depth:
        B = np.array(boxes)
        area = (B[:, 1] - B[:, 0]) * (B[:, 3] - B[:, 2])
        if depth == 0:
            total_area = float(area.sum())
        Z, Dn = interval_zeta_D(A, (B[:, 0], B[:, 1]), (B[:, 2], B[:, 3]))
        pos = Dn[0] > 0
        neg = Dn[1] < 0
        dec = pos | neg
        for k in np.where(pos)[0]:
            pos_boxes.append((B[k], area[k], (Z[0][0][k], Z[0][1][k], Z[1][0][k], Z[1][1][k])))
        neg_area += float(area[neg].sum())
        rest = B[~dec]
        if depth == max_depth or len(rest) == 0 or 4 * len(rest) > max_boxes:
            unresolved_area = float(area[~dec].sum()) if len(rest) else 0.0
            boxes = []
            break
        boxes = []
        for (a2, b2, a3, b3) in rest:
            m2, m3 = 0.5 * (a2 + b2), 0.5 * (a3 + b3)
            for sub in ((a2, m2, a3, m3), (m2, b2, a3, m3), (a2, m2, m3, b3), (m2, b2, m3, b3)):
                if sub[1] + sub[3] <= 1 - delta + 1e-15:
                    boxes.append(sub)
                else:  # sub-box pokes past l1=delta: dropped, counted unresolved
                    dropped_area += (sub[1] - sub[0]) * (sub[3] - sub[2])
        depth += 1
    pos_area = float(sum(a for _, a, _ in pos_boxes))
    max_mult, mult_hist = None, {}
    if pos_boxes and sample_pts is not None and len(sample_pts):
        enc = np.array([z for _, _, z in pos_boxes])  # z1lo,z1hi,z2lo,z2hi
        counts = []
        for p in sample_pts:
            inside = ((enc[:, 0] <= p[0]) & (p[0] <= enc[:, 1]) &
                      (enc[:, 2] <= p[1]) & (p[1] <= enc[:, 3]))
            counts.append(int(inside.sum()))
        counts = np.array(counts)
        # neighbouring boxes enclose overlapping Z-ranges by continuity, so
        # raw multiplicity counts boxes, not preimages; report it as the
        # enclosure-overlap diagnostic it is (see gap G3 in the writeup).
        max_mult = int(counts.max())
        mult_hist = dict(median=float(np.median(counts)),
                         frac_zero=float((counts == 0).mean()),
                         median_Z_width=[float(np.median(enc[:, 1] - enc[:, 0])),
                                         float(np.median(enc[:, 3] - enc[:, 2]))])
    return dict(delta=delta, total_area=total_area,
                pos_area=pos_area, neg_area=neg_area,
                unresolved_area=float(unresolved_area + dropped_area),
                n_pos_boxes=len(pos_boxes),
                enclosure_overlap_max=max_mult, enclosure_overlap_stats=mult_hist)


def detailed_shape(ab):
    a, b = ab
    A = tri(a, b)
    X, z, D = inverse_grid(A, N=160 if QUICK else 240)
    rmin_pts = z[D > 0]
    res = dict(a=a, b=b, grid_interior=int(len(X)),
               grid_Rmin=int((D > 0).sum()),
               jac=jac_det_sign(A, X, z, D),
               forward=forward_checks(A, z, D, kmax=12 if QUICK else 30))
    res["interval"] = interval_analysis(
        A, delta=0.02, n0=12 if QUICK else 20, max_depth=4 if QUICK else 6,
        sample_pts=rmin_pts[:, :2] if len(rmin_pts) else None)
    return res


GAPS = """## What a full computer-assisted proof for n=3 still needs

Claim to prove: for every non-degenerate triangle, Z is one-to-one on
R_min = {x in interior : D(x) > 0}. The items below are the gaps between
this module and that claim. None is filled here.

- **G1. Covering the shape space.** Shapes form a 2-parameter family
  (a,b) with b>0, non-compact at b=0. A proof needs (i) a finite box
  cover of a compact subset {b >= b0} x [0,1/2] with all estimates made
  uniform over each shape-box (interval arithmetic in (a,b) as well as
  in x -- this module fixes the shape as exact floats), and (ii) a
  separate asymptotic argument for flat triangles b < b0 (the equilibria
  then collapse towards the long side; the three-charge scan in NOTE.md
  already could not resolve R_min there).
- **G2. The boundary strips of the triangle.** This module certifies D
  only on {l1,l2,l3 >= delta}, delta=0.02. A proof must also handle
  x near the edges and vertices, where Z degenerates (one charge ratio
  tends to 0) and r_i -> 0 makes interval widths blow up. Needs
  rescaled local coordinates near each edge/vertex, or an analytic
  sign computation for D in those regimes.
- **G3. From local to global injectivity.** Certified D>0 plus Z
  enclosures do not give injectivity: neighbouring boxes have
  overlapping Z-enclosures by continuity, and one box could in
  principle contain two preimages of one charge ratio. A proof needs
  an extra topological step, e.g. (i) interval-enclose the Jacobian DZ
  (analytic derivatives, not finite differences) and certify det DZ > 0
  on every D>0 box, making Z a local diffeomorphism on R_min, then
  (ii) certify Z is one-to-one on the boundary of R_min (the fold
  curve D=0 and the triangle boundary), then conclude by a degree /
  monotonicity argument (a proper local homeomorphism injective on
  the boundary of a disc-like region is injective). This also needs
  control of the topology of R_min (connected components, and that
  Z maps different components to disjoint charge regions or counts
  them correctly).
- **G4. The fold curve itself.** Boxes straddling D=0 stay unresolved
  at every finite subdivision depth; their area is reported above and
  does not vanish by subdivision alone near cusp points of the fold.
  A proof needs the fold treated as a curve: e.g. certify D=0 is a
  smooth curve (interval gradient of D bounded away from 0 on
  straddling boxes, i.e. no cusp or a classified cusp) and track
  Z restricted to it.
- **G5. Charge-ratio ends.** Injectivity must hold for all charge
  ratios, including ratios tending to the simplex boundary (one charge
  -> 0 or dominating). These correspond under Z to x near edges/vertices
  (G2) and need the same asymptotic treatment; the raster and forward
  checks here only sample ratios attained at interior grid points.
- **G6. Independent checking.** The n=4,5,6 counts in this repo are
  certified twice (certify.py + verify_arb.py, sharing no code). A
  proof-grade n=3 argument should emit a certificate (box lists and
  enclosures) that a second, independent verifier (e.g. Arb-based,
  like verify_arb.py) can recheck. This module emits numbers, not a
  recheckable certificate.
- **G7. Termination bounds.** Subdivision depth here is capped by hand.
  A proof needs explicit width estimates for the interval extensions
  of D, Z and DZ (they are rational functions of x, A with square
  roots) guaranteeing the caps in G1-G4 are met, rather than observed
  to suffice on sampled shapes.

What would already be a publishable partial result (issue #4): G1(i)+G3+G4
on an explicit compact set of shapes, e.g. all triangles with b >= b0 for
some concrete b0, with delta-strips handled for those shapes. The sampled
shapes here suggest the margins (det DZ, D away from folds) are not
pathologically small away from flat triangles, but that is CONJECTURE
until the uniform box estimates exist.
"""


def main():
    t0 = time.time()
    shapes = shape_grid()
    if QUICK:
        shapes = shapes[::4]
    print(f"coverage on {len(shapes)} shapes", flush=True)
    with Pool(8) as pool:
        cov = pool.map(coverage_job, shapes)
    json.dump(dict(stage="coverage_done", coverage=cov), open(JPATH, "w"), indent=1)
    print(f"coverage done in {time.time()-t0:.0f}s", flush=True)

    detailed_shapes = [(0.5, 0.866), (0.3, 0.4), (0.0, 0.5), (0.15, 0.2),
                       (0.5, 0.1), (0.25, 0.65)]
    with Pool(min(6, len(detailed_shapes))) as pool:
        det = pool.map(detailed_shape, detailed_shapes)
    print(f"detailed done in {time.time()-t0:.0f}s", flush=True)

    resolved = [c for c in cov if c["resolved"]]
    maxcov = max((c["max_coverage"] for c in cov), default=0)
    report = dict(stage="done", quick=QUICK, seconds=round(time.time() - t0, 1),
                  coverage=cov, detailed=det,
                  max_coverage=maxcov, n_shapes=len(cov), n_resolved=len(resolved))
    json.dump(report, open(JPATH, "w"), indent=1, default=str)

    L = [
        "# n=3 inverse map: finite checks and the remaining proof gaps (issue #4)",
        "",
        "Script: `explore/n3_inverse_map.py`. Reformulation: NOTE.md S3. For each",
        "triangle, Maxwell's bound 4 holds iff the inverse map Z (Lami's theorem)",
        "is one-to-one on R_min = {D>0}.",
        "",
        f"Runtime {report['seconds']}s.",
        "",
        "## Finite checks run",
        "",
        f"- Rasterised coverage of Z(R_min) on the charge simplex, {len(cov)} shapes "
        f"in the fundamental domain (side 1 longest, A3=(a,b), a in [0,1/2]): "
        f"{len(resolved)} shapes have a resolved R_min; maximum coverage multiplicity "
        f"**{maxcov}**. Multiplicity 1 everywhere on resolved shapes is the finite",
        "  analogue of injectivity; >=2 anywhere would flag a candidate counterexample",
        "  (then re-solved directly, as in results/three_scan_summary.txt).",
        "- Detailed shapes: float inverse-map grid, det DZ sign by central differences",
        "  on R_min grid points, forward re-solves of sampled charge ratios",
        "  (`three.crit2d`, 2 seeds each), and interval certification with `ia.py`.",
        "",
        "| shape (a,b) | R_min grid pts | det DZ >0 / <0 (sampled) | forward m0 counts | interval D>0 area | interval D<0 area | unresolved area | Z-enclosure overlap max |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for d in det:
        iv = d["interval"]; j = d["jac"]; f = d["forward"]
        L.append(f"| ({d['a']},{d['b']}) | {d['grid_Rmin']} | {j['n_positive']} / {j['n_negative']} "
                 f"| {f['m0_counts']} (n={f['n']}) | {iv['pos_area']:.4f} | {iv['neg_area']:.4f} "
                 f"| {iv['unresolved_area']:.4f} | {iv['enclosure_overlap_max']} |")
    L += [
        "",
        "Areas are in barycentric (l2,l3) measure; the shrunk triangle "
        "{li >= 0.02} has area (1-3*0.02)^2/2 = 0.4418, of which the initial box "
        "cover captures `total_area` per shape (JSON). Unresolved area is the fold",
        "neighbourhood D~0 plus boxes dropped at the l1=0.02 cut (see G4, G2).",
        "Z-enclosure overlap counts *boxes* whose interval Z-range contains a",
        "sampled charge point; neighbouring boxes overlap by continuity, so this",
        "is a diagnostic of enclosure width, **not** a preimage count (see G3).",
        "",
        "## Claim labels",
        "",
        "- PROVEN (modulo ia.py, itself checked against exact rationals by",
        "  tests/test_ia.py): on every box reported D>0 above, D really is >0",
        "  throughout the box, for that exact float triangle; likewise D<0 boxes.",
        "  Same for the Z-enclosures containing the true Z on their boxes.",
        "- VERIFIED (finite numerical check): raster coverage multiplicities,",
        "  det DZ signs at sampled points, forward minima counts m0 at sampled",
        "  charge ratios. These concern only the shapes/points actually run.",
        "- CONJECTURE: that the margins seen on sampled shapes persist uniformly",
        "  over compact shape sets away from flat triangles.",
        "- NOT proven: injectivity of Z on R_min for any single shape, and",
        "  Maxwell's bound for n=3 in any regime. See the gap list.",
        "",
        GAPS,
        "Full data: `results/n3_inverse_map.json`.",
    ]
    open(MPATH, "w").write("\n".join(L) + "\n")
    print(f"max_coverage={maxcov} resolved={len(resolved)}/{len(cov)} "
          f"wrote results/n3_inverse_map.{{json,md}}", flush=True)


if __name__ == "__main__":
    main()

