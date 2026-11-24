# 15 November 2026: balance points can grow faster than the number of charges

This note adds to issue 5 (eight or more charges). It follows `results/2026-11-15-growth-per-charge.md`, which describes the
constructions known before this one. Marks: **[F]** fact from a run, **[P]** a published result, **[I]** inference,
**[C]** conjecture. Everything numerical is uncertified except where a proof is named.

## Result

[I] With n positive charges the number of balance points can grow like n^(3/2). Every construction known before this grows like 10n.

Status: an argument in four steps with a proof sketch, numerical checks of every step, and two exact lower bounds in actual charges.
It is not a finished proof. No count for any particular number of charges is beaten, so the table of certified counts in the README is unchanged.
On the numbers below this route passes the known constructions only near 160 charges.

## An inversion rule

[F] Take charges q_i at sites a_i, at distances r_i from a point O. Reflect them in the unit sphere about O: charge m_i = q_i / r_i at b_i = a_i / r_i^2.
Then the term of degree k in the Taylor series of the potential at O is the multipole of order k of the reflected charges. The dipole gives the gradient, the
quadrupole the Hessian, and so on.

The proof is the Legendre expansion. The degree-k term is the sum of q_i P_k(cos) |x|^k / r_i^(k+1), and q_i / r_i^(k+1) is m_i |b_i|^k.
`explore/growth/kelvin_check.py` confirms it for seven random charges to 12 digits, for k from 0 to 6 (`results/growth/kelvin_check.out`).

So O is flat to order k, meaning every derivative of order 1 to k vanishes, exactly when the reflected charges have no multipole of order 1 to k.

[P] Equal weights at the points of a spherical k-design have no multipole of order 1 to k about the centre. The regular solids are designs: the tetrahedron of order 2,
the octahedron and cube of order 3, the icosahedron and dodecahedron of order 5. `explore/growth/platonic.py` checks the order for equal charges at the vertices
(`results/growth/platonic.out`).

## The mechanism

**Step 1. A point that is flat to high order.** By the inversion rule, charges can be placed so that the potential at a point has no terms of degree 1 to d - 1.

- [P] Equal charges at the points of a spherical design do this, and designs exist with a number of points that grows like the square of their order (Bondarenko, Radchenko and Viazovska, 2013).
- [F] Solving for sites and charges together finds such points with fewer charges: 6, 9, 12, 16 and 21 charges for d = 4 to 8. That is about (d^2 - 1)/3. None was found with 5 charges for d = 4 or with 8 for d = 5. `explore/growth/nflat.py` runs the search. The saved runs in `results/growth/` cover d = 4, 5 and 6 only.

**Step 2. One level.** Let A be the zonal harmonic of degree p about the z axis and B the zonal harmonic of degree p - 1 about the x axis.

[F] A + B has exactly I(p) - 1 real, non-degenerate critical points away from the origin, where I(p) is p^2/2 rounded down, plus 2 when p leaves remainder 3 on division by 4.

| p | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|
| I(p) - 1 | 5 | 7 | 11 | 17 | 25 | 31 | 39 | 49 | 61 | 71 |

- [F] The search found these counts for p = 3 to 14 by scanning a grid over the sphere, and for p = 3 to 9 by a homotopy method. Those programs are not part of this branch.
- [F] A program written apart from both finds the same counts for p = 3 to 12 (`explore/growth/levels_check.py`, output `results/growth/levels_check.out`). It stops with an error on any mismatch. It uses a different fact: a critical point away from the origin lies in a direction where the gradients of A and B are parallel, and its distance is then fixed.
- [I] Proof sketch. These directions are the critical points, on the sphere with opposite points identified, of a function that goes to minus infinity on the zero set of A and plus infinity on the zero set of B. The two zero sets cut the surface into faces. A face with 2m sides forces m - 1 critical points, and Euler's formula sums these to I(p) - 1.

**Step 3. Stacking levels.** Take the sum of zonal harmonics of degrees 1 to D, with axes alternating between z and x and with coefficients chosen so that each degree rules at its own distance from the origin. Each level then behaves as in step 2.

[F] The sum has exactly S(D) real critical points, where S(D) adds I(p) - 1 for p from 2 to D.

| D | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|
| S(D) | 6 | 13 | 24 | 41 | 66 | 97 | 136 | 185 |

- [F] The search finds these exactly for D = 3 to 10 with the levels a factor of 100 apart (`explore/growth/hier.py`, `results/growth/hier_s0.01.out`). With a factor of 33 a few points are lost from D = 7 up.
- [F] A separate program with exact coefficients, `explore/growth/stack_check.py`, proves 6, 13, 24 and 41 real non-degenerate critical points for D = 3 to 6, each split by scale as predicted. It encloses every point in a cube by a ball-arithmetic test. Outputs are `results/growth/stack_check_D3.out` to `stack_check_D6.out`.
- [F] An earlier version of that check used numerical derivatives and reported 59 at D = 6. It had counted poorly converged points twice. It is withdrawn, and the exact program replaces it.
- [I] S(D) grows like D^3/6.

**Step 4. Charges.** Start from a point that is flat to order d - 1. Small changes of the charges set the terms of degree 1 to d - 1 to any small values, here the stacked levels of step 3. The balance points near the flat point are then the critical points counted in step 3, once the scale is small enough.

[I] With about d^2/3 charges and about d^3/6 balance points, the count grows like 0.87 n^(3/2).

## A proof of step 2 as a lower bound

The sketch in step 2 needs every face to be a disc. The argument below does not, and it is complete as a lower bound. It was written for this review and has not been checked by a second person.

**Claim.** Let A and B be harmonic forms of degrees p and p - 1 in three variables. Suppose their zero sets on the unit sphere are smooth curves that cross transversally, at 2I points. Suppose also that the function phi below has only non-degenerate critical points. Then A + B has at least I - 1 non-degenerate real critical points away from the origin.

**1. Reduction to the sphere.** Write a point as x = t w with w a unit vector. Then (A + B)(t w) = t^(p-1) (t A(w) + B(w)). The derivative in t vanishes only at t = -(p - 1) B(w) / (p A(w)). Where A(w) is not zero this is a non-degenerate zero in t.

Putting that t back in gives a function of w alone, a constant times B^p / A^(p-1). Its critical points are those of phi = (p - 1) log|A| - p log|B|, and they are non-degenerate together. The directions w and -w give the same point x. So the critical points of A + B away from the origin, with A and B not zero there, match the critical points of phi on the sphere with opposite points identified.

**2. A one-form with isolated zeros.** Put omega = (p - 1) B dA - p A dB on the sphere. It is smooth everywhere. Where A and B are both non-zero, omega is A B times d phi, so its zeros there are the critical points of phi with the same index: -1 at a saddle and +1 at a maximum or minimum.

Where exactly one of A and B vanishes, omega is a non-zero multiple of dA or dB, which does not vanish on a smooth zero curve. Where both vanish, at a crossing, take A and B as local coordinates. Then omega is the linear field with matrix rows (0, p - 1) and (-p, 0), whose determinant p(p - 1) is positive. So each crossing is a zero of index +1.

**3. Counting.** The indices of the zeros of a one-form on the sphere add up to 2 (Poincaré and Hopf). The 2I crossings give 2I. The critical points of phi come in opposite pairs, so they give twice their sum over the sphere with opposite points identified. That sum is therefore 1 - I.

So the number of saddles of phi, minus the number of its maxima and minima, is I - 1. There are at least I - 1 saddles, each giving one non-degenerate critical point of A + B.

**4. The hypotheses hold for a dense set of pairs.** Transversal crossing is stable under small changes of A and B. Harmonic forms of degree p can take any value and any tangential gradient at a given point of the sphere, so by the standard transversality argument phi is non-degenerate for a dense set of pairs near any given one.

**For the crossed zonal harmonics.** The zero sets are circles z = a and x = b, with a and b zeros of the Legendre polynomials of degrees p and p - 1. Two such circles cross, transversally, exactly when a² + b² is below 1. So I counts the pairs of zeros with a² + b² below 1.

- [F] That count is the I(p) of step 2 for every p up to 400.
- [I] It grows like p²/2, because the zeros are spread like the cosines of equally spaced angles. Classical bounds on Legendre zeros turn this into a proof that I is at least a fixed multiple of p².
- [F] In the runs the count of critical points is exactly I(p) - 1, so phi has no maximum or minimum for p up to 14. [I] A minimum is ruled out in general: at a critical point the Laplacian of phi on the sphere is negative. A maximum is not ruled out by this argument, and the lower bound does not need it.

**Steps 3 and 4, what a full proof still needs.** Each level's critical points are non-degenerate, so they survive the small terms from the other levels once the levels are far enough apart. That is the implicit function theorem applied level by level, with the scale ratio chosen after the levels.

For step 4, take equal charges on a design of order 2(d - 1). The harmonics of degree up to d - 1 are then orthogonal on its points, so the charges alone reach every small set of terms of degree 1 to d - 1, and small changes keep the charges positive.

Terms of degree d and above change too, by amounts of the same size, and are smaller by a factor of the distance at the scales used. Neither step is written out in full here.

## The same thing in actual charges

Two configurations were built with exact rational sites and charges: `explore/configs/growth-flat-d5-n9.json` and `explore/configs/growth-flat-d8-n21.json`.
They are **not certificates** in this repository's format and are not under `certificates/`, because no exact count was completed. For each, a lower bound is proved by
`lower_bound_hp.py`. It trusts no candidate: each is refined in 260-digit arithmetic and then enclosed in a cube where a ball-arithmetic test proves exactly one balance point.
The candidates are in `explore/growth/data/`. `explore/growth/demo_flat.py` builds the configurations from the flat starting points there. Rerun here, it reproduces the 21-charge configuration exactly and the 9-charge one to 1e-34. The committed 9-charge configuration has the smaller jet error, 5e-43. The output is `results/growth/demo_flat_reproduces_configs.out`. The proved bounds are for the committed files.
`explore/growth/dump_taylor_points.py` lists the candidates near the flat point. For the 21-charge configuration it gave 87 points, 40 of index 1 and 47 of index 2 (`results/growth/dump_taylor_points_d8_n21.out`).

| Charges | Flat to order | Proved balance points | Of index 2, of index 1 | Groups by distance from the flat point |
|---|---|---|---|---|
| 9 | 4 | at least 26 | 17, 9 | 1, 5 and 7 at three scales, then 13 further out |
| 21 | 7 | at least 87 near the flat point | 47, 40 | 1, 5, 7, 11, 17 and 25 at six scales, then 21 at the top level |

- [F] The proved lower bounds are checked by `tests/test_lower_bound_hp.py`. The 9-charge case takes a few seconds. The 21-charge case takes several minutes and runs only with `--slow`. Its output is `results/growth/lower_bound_hp_n21.out`.
- [F] The groups are the predicted I(p) - 1 for p = 2 to 4 and p = 2 to 7.
- [F] For 9 charges the 26 satisfy the index relation, so they are plausibly all of them. That is not proved: the exact count prover did not finish in an hour.
- [F] For 21 charges the 87 are only those within 0.06 of the flat point. The index relation shows at least 13 more further out.
- [F] These are needles. In the 21-charge example the level coefficients run from 1.8 down to 1.8e-54, and the innermost balance point is 1e-14 from the flat point.

## What it means

- [I] The largest number of real critical points of a harmonic polynomial of degree d in three variables is at least of order d^3/6. In two variables the answer is d - 1 [P], so three variables do not behave like two.
- [I] Walls between charges, the index-2 points, also grow faster than the number of charges. So no bound of the form "a fixed number of walls per charge" can hold.
- [I] Maxwell's bound of (n - 1)^2 is not threatened by this for large n. The route here gives n^(3/2), and at most 8 n^(3/2) even if every complex solution were real.
- [F] It does not beat any known count yet. The best known counts are 10n - 26 and 10n - 31.

| d | Charges needed, about | Balance points from the levels | Best known for that many charges |
|---|---|---|---|
| 8 | 21 | 66 | 184 |
| 12 | 48 | 246 | 449 |
| 20 | 133 | 1,222 | 1,304 |
| 22 | 161 | 1,640 | 1,584 |
| 28 | 261 | 3,446 | 2,584 |

## Other findings of the search (numerical, programs not included)

- [F] Random clusters of 5 to 30 charges, chains and lattice blocks have n - 1 balance points, the least possible. Rings have n + 1.
- [F] Equal charges on the regular solids have 3, 5, 17, 11 and 43 balance points once the degenerate centre is split. None has a balance point on a face.
- [F] A tiny charge placed at a balance point gives new points at the positive critical points of the leading harmonic on the sphere. This is confirmed for the cube, where it adds 20.
- [F] Small clusters placed at a flat point gave 4 to 7 new balance points per charge, below the 10 of the known constructions.
- [F] Climbing at random reached 19 of 27 real critical points for a harmonic quartic and 26 of 64 for a quintic. The stacked levels give 13 and 24, so they are not the best for small degree.

## Not verified

- The proof of step 2 below is complete only as a lower bound and for a dense set of pairs of forms. For the exact crossed zonal pair, non-degeneracy is numerical, for p up to 14. Nobody else has checked the proof.
- The scaling argument that carries step 3 into step 4 is standard and is not written out.
- That changing the charges alone reaches every small set of low-degree terms is argued from the orthogonality of harmonics on a design. It was checked by rank for 6, 9 and 21 charges only.
- No exact count was completed for either configuration. The lower bounds are proved.
- No earlier statement of growth faster than linear was found in arXiv:2607.27197, arXiv:2501.05315 or arXiv:2609.09637, or by search. That is not proof that it is new.

## Next

1. Write out step 2 and the scaling step as full proofs.
2. Find the best real count for a harmonic polynomial of each small degree. The stacked levels are a lower bound and random climbs already beat them at degree 4 and 5.
3. Look for a level pair that fills more of the available solutions without spreading over many scales. That would bring the crossing point down from 160 charges.
4. Whether a fixed share of the d^3 critical points of the leading form of a design can be made real is the same question for the regular solids. For the icosahedron, 12 charges, the leading form can have up to 125 critical points, and the best count known for 12 charges is 89. This is open.

## Rerun

With NumPy, SciPy, mpmath and python-flint, from the repository root:

```
python3 explore/growth/levels_check.py
python3 -I -O lower_bound_hp.py explore/configs/growth-flat-d5-n9.json explore/growth/data/n9_candidates.npy 26
python3 -I -O lower_bound_hp.py explore/configs/growth-flat-d8-n21.json explore/growth/data/n21_candidates.npy 87
python3 tests/test_lower_bound_hp.py            # add --slow for the 21-charge case
sh results/growth/run_programs.sh               # reruns the programs of explore/growth/, about 15 minutes
```
