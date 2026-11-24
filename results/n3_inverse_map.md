# n=3 inverse map: finite checks and the remaining proof gaps (issue #4)

Script: `explore/n3_inverse_map.py`. Reformulation: NOTE.md S3. For each
triangle, Maxwell's bound 4 holds iff the inverse map Z (Lami's theorem)
is one-to-one on R_min = {D>0}.

Runtime 5.6s.

## Finite checks run

- Rasterised coverage of Z(R_min) on the charge simplex, 22 shapes in the fundamental domain (side 1 longest, A3=(a,b), a in [0,1/2]): 19 shapes have a resolved R_min; maximum coverage multiplicity **1**. Multiplicity 1 everywhere on resolved shapes is the finite
  analogue of injectivity; >=2 anywhere would flag a candidate counterexample
  (then re-solved directly, as in results/three_scan_summary.txt).
- Detailed shapes: float inverse-map grid, det DZ sign by central differences
  on R_min grid points, forward re-solves of sampled charge ratios
  (`three.crit2d`, 2 seeds each), and interval certification with `ia.py`.

| shape (a,b) | R_min grid pts | det DZ >0 / <0 (sampled) | forward m0 counts | interval D>0 area | interval D<0 area | unresolved area | Z-enclosure overlap max |
|---|---|---|---|---|---|---|---|
| (0.5,0.866) | 1333 | 400 / 0 | {'1': 30} (n=30) | 0.0278 | 0.3763 | 0.0156 | 1431 |
| (0.3,0.4) | 1125 | 400 / 0 | {'1': 30} (n=30) | 0.0216 | 0.3852 | 0.0129 | 2850 |
| (0.0,0.5) | 1192 | 400 / 0 | {'1': 30} (n=30) | 0.0207 | 0.3862 | 0.0127 | 2476 |
| (0.15,0.2) | 840 | 400 / 0 | {'1': 30} (n=30) | 0.0112 | 0.3971 | 0.0113 | 2332 |
| (0.5,0.1) | 351 | 351 / 0 | {'1': 30} (n=30) | 0.0062 | 0.4056 | 0.0079 | 3046 |
| (0.25,0.65) | 1308 | 400 / 0 | {'1': 30} (n=30) | 0.0261 | 0.3794 | 0.0142 | 2189 |

Areas are in barycentric (l2,l3) measure; the shrunk triangle {li >= 0.02} has area (1-3*0.02)^2/2 = 0.4418, of which the initial box cover captures `total_area` per shape (JSON). Unresolved area is the fold
neighbourhood D~0 plus boxes dropped at the l1=0.02 cut (see G4, G2).
Z-enclosure overlap counts *boxes* whose interval Z-range contains a
sampled charge point; neighbouring boxes overlap by continuity, so this
is a diagnostic of enclosure width, **not** a preimage count (see G3).

## Claim labels

- PROVEN (modulo ia.py, itself checked against exact rationals by
  tests/test_ia.py): on every box reported D>0 above, D really is >0
  throughout the box, for that exact float triangle; likewise D<0 boxes.
  Same for the Z-enclosures containing the true Z on their boxes.
- VERIFIED (finite numerical check): raster coverage multiplicities,
  det DZ signs at sampled points, forward minima counts m0 at sampled
  charge ratios. These concern only the shapes/points actually run.
- CONJECTURE: that the margins seen on sampled shapes persist uniformly
  over compact shape sets away from flat triangles.
- NOT proven: injectivity of Z on R_min for any single shape, and
  Maxwell's bound for n=3 in any regime. See the gap list.

## What a full computer-assisted proof for n=3 still needs

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

Full data: `results/n3_inverse_map.json`.
