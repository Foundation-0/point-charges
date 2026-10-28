# Wave-2 Task A: the 332 Morse-inconsistent wave-1 screens, re-examined

Script: `explore/wave2_hardening.py`. Wave 1 stored only tallies for its
332 Morse-inconsistent screens, so the exact task lists of the two affected
families were regenerated from the wave-1 generators (same seed 20261018,
same code) and re-screened with `count_fast(seed=0)`.

- Regenerated tasks: 9375. Inconsistent on re-screen: 332 ({'tri_interior': 327, 'symmetric_tet': 5}). Wave-1 tally: 332 ({'tri_interior': 327, 'symmetric_tet': 5}).
- Runtime 260.5s.

## Classification

### Degenerate-coincident: 327

Two sites coincide (relative min pairwise distance < 1e-9). The analytic
count for tri_interior is 297 (main grid: 9 triangles x 3 barycentric
vertices x height 0 x 11 fourth charges) + 30 (unequal-charge part:
5 triples x 1 vertex in bary[::3] x 6 charges) = 327.
Coincident tasks total 327; of those 327 screened
inconsistent and 0 screened consistent. A coincident pair
is one charge with the summed value, i.e. an n=3 configuration with
m2-m1=2; wave 1 required m2-m1=3 (n=4), so inconsistency is definitional.
Merging each coincident pair and re-counting with `count_fast`: 327/327 Morse-consistent for n=3, count distribution {2: 298, 4: 29}.

### Distinct sites: 5 screens, 1 unique configuration(s)

Wave-1's gen_symmetric emits the all-ones charge pattern 5 times (initial
pattern + r=1 from 10**linspace(-4,1,11) in each of 4 pattern shapes), so
duplicate screens are classified once with their multiplicity.

The slow solver `charges.summary` could not be used on these cases: both
(nseed=8000, rounds=4) and (nseed=2500, rounds=3) ran for minutes with no
result (killed; see failed attempts). The cause is the same cloud that
breaks the screen: near a degenerate equilibrium Newton converges slowly,
hundreds of iterates pass the gradient test, and the O(k^2) dedup in
find_critical_points / the 1e-6-relative dedup in find_fast turn one
equilibrium into hundreds of 'points'. Re-examination therefore used
`count_fast` seeds 1-3, a cloud diagnostic on find_fast's points, and a
clustered heavy multi-start counter (`clustered_count` in this script:
28k seeds, strict gradient test, clustering at 1e-2*scale (3e-3 reported), polish, two
seeds), plus 1e-2 charge perturbations and a certify.py attempt.
Verdicts (weighted by multiplicity): {'solver_artifact_degenerate': 5}. True counts of solver artifacts (unique configs): {1: 1}.

For the regular tetrahedron with equal charges the centre is an equilibrium
with Hessian exactly 0: the tetrahedral symmetry group acts irreducibly on
R^3, so an invariant symmetric Hessian is a scalar matrix, and harmonicity
(trace 0) forces the scalar to 0. It is therefore degenerate, Maxwell's
non-degenerate Morse bookkeeping does not apply, and certify.py's Krawczyk
test (which proves non-singularity) cannot verify it. Perturbing the
charges by 1e-2 restores Morse-consistent counts (see JSON), and the
stretched tetrahedra (s != 1) in wave 1 screened consistent.

- symmetric_tet x5 q=['1', '1', '1', '1'] a=[[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]] screen(n,m1,m2)=[467, 3, 464] -> solver_artifact_degenerate, true count 1 (1 degenerate); clustered n=[1, 1]; cloud={'n': 467, 'max_dist_to_centroid': 0.0014627524972960428}; certify={'error': 'RuntimeError: box budget exceeded'}

**Genuine candidates with count > 9: 0. Unresolved: 0.**
certify.py was attempted on each unique distinct config (results above);
it certifies non-degenerate zeros only, so it cannot settle the degenerate
tetrahedron centre. No >9 candidate exists to certify.

## Claim labels

- VERIFIED (finite numerical check): the classifications and true counts above,
  for the regenerated configurations actually run. Regeneration is exact only
  if the re-screen tally reproduces 332 / {327, 5}; see the numbers above.
- The coincident classification is additionally PROVEN by the grid arithmetic
  (297+30=327) together with the merged recount: coincident sites are an n=3
  configuration, for which m2-m1=2 by Morse theory (NOTE.md S1).
- PROVEN: the regular-tetrahedron centre is a degenerate equilibrium
  (Hessian 0 by irreducible tetrahedral symmetry + trace 0).
- VERIFIED (finite): its clustered true count 1 (two seeds, 28k starts;
  ray scans and perturbations in exploration agreed). NOT certified by
  certify.py/verify_arb.py: both require non-degenerate zeros.
- NOT a proof about n=4 in general: nothing here covers configurations wave 1
  did not generate.

## Failed attempts

- `charges.summary(nseed=8000, rounds=4)` on the distinct case: no result
  in >7 min per config (workers at 98% CPU); killed. `summary(2500, 3)`:
  also no result in >4 min; killed. The bottleneck is the O(k^2) dedup in
  find_critical_points over the cloud of Newton iterates around the
  degenerate centre. Replaced by `clustered_count`.
- certify.py on the regular tetrahedron: RuntimeError 'box budget exceeded'
  (expected: its Krawczyk test proves Hessian non-singularity, false at the
  degenerate centre). verify_arb.py not run: no certificate exists to check.

Full data: `results/wave2_hardening.json`.
