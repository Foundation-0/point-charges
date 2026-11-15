# 15 November 2026: four charges with eleven balance points

This note adds to issue 3 (can four charges have more than 9 equilibria?) and issue 5. It follows the head `a2fd159`.
Marks: **[F]** fact from a run, **[P]** a published result, **[I]** inference, **[C]** conjecture.

## Result

[F] Four positive point charges can have exactly 11 balance points, all non-degenerate: 7 of index 2 and 4 of index 1.
Maxwell's bound for four charges is 9.

[F] Two configurations were built, by two different methods and two separately written sets of programs. For each one,
two exact programs prove the count is exactly 11, and a third program proves it is at least 11.

Status: claimed in this repository's terms. In the terms of the review that found it, VERIFIED, meaning separately written
programs accept it. Nobody outside the project has checked it.

## A short certificate

[F] The flat-route configuration was made sturdier and then rounded to five decimal places. It still has exactly 11 balance points, by `verify_arb.py`, `indep_count.py` and `lower_bound_check.py`. Certificate: `certificates/n4-eleven-short.json`. Plain sites and charges: `explore/configs/n4-eleven-short.json`.

| Charge | x | y | z |
|---|---|---|---|
| 0.26712 | -0.35180 | -0.27373 | -0.20165 |
| 0.28010 | -0.05642 | 0.44608 | -0.21436 |
| 1 | -0.00274 | 0.00705 | 0.84585 |
| 0.26274 | 0.40862 | -0.17404 | -0.19749 |

The 11 balance points, to five decimals, with their index:

```
(-0.00021, +0.00022, +0.00155)  2      (-0.00540, +0.00212, -0.00514)  2
(+0.00014, -0.00055, -0.00251)  1      (+0.00583, -0.00243, +0.00455)  1
(+0.00041, -0.00328, -0.00342)  2      (+0.01945, -0.00805, +0.01529)  2
(+0.00395, +0.00287, -0.00454)  2      (-0.01747, -0.01312, +0.01561)  2
(-0.00487, -0.00360, +0.00427)  1      (-0.00307, +0.02304, +0.01606)  2
(-0.00088, +0.00615, +0.00431)  1
```

Anyone can test a listed point with a few lines of code: the field there is zero to the digits shown, and the signs of the three Hessian eigenvalues give the index. Rounding to four decimals leaves 9 and to three leaves 5.

[F] It was made by two numerical climbs, `explore/four_eleven/fatten.py` and `fatten2.py`, which move the sites and charges so that the smallest Hessian eigenvalue among the 11 points grows, and then `short_cert.py`. The climbs are a search only. The certificate stands on its own. It survives relative changes of about 1e-6, ten times more than the flat-route configuration.

## The two configurations

| | Jet route | Flat route |
|---|---|---|
| Certificate | `certificates/n4-eleven-jet.json` | `certificates/n4-eleven-flat.json` |
| Plain sites and charges | `explore/configs/n4-eleven-jet.json` | `explore/configs/n4-eleven-flat.json` |
| Charges | 1, 0.851, 0.864, 0.838 | 1, 0.280, 0.267, 0.263 |
| Shape | a tetrahedron about 0.05 from the regular one | one large charge at 0.846 from the flat point, three small ones at 0.49 |
| Balance points | 8 within 0.001 of the design point, 3 at 0.011 to 0.018 | 8 within 0.0026 of the flat point, 3 at 0.027 to 0.030 |
| Coordinates | 12-digit decimals | 14-digit decimals |

Both certificates are in the repository's format. They were written by `explore/to_repo_format.py` from the zeros found by
`indep_count.py`, as described in `results/2026-11-15-third-checker.md`.

## Checks

| Program | What it proves | Jet route | Flat route |
|---|---|---|---|
| `indep_count.py`, run with `-O` | exactly N, by covering the whole region with boxes that hold no zero | 11, from 36,897 boxes | 11, from 33,873 boxes |
| `verify_arb.py` (with the fix for issue 21), plain and with `-O` | exactly N | VERIFIED 11 | VERIFIED 11 |
| `lower_bound_check.py`, new, run with `-O` | at least N | 11 | 11 |

- `indep_count.py` reads only the sites and the charges. It is the third checker of `results/2026-11-15-third-checker.md`.
- `verify_arb.py` is the repository's own verifier, run both ways because the fix replaced its `assert` statements.
- `lower_bound_check.py` was written for this result without reading the other two. It finds its own candidates, then proves one zero in each small cube with the Krawczyk test in ball arithmetic, and fixes the index from the sign of the determinant. It proves "at least N" only.
- The exact commands and outputs are in `results/four-eleven/checks.log`.

[F] Controls, all of which fail as they should (`tests/test_four_eleven.py`):

- `indep_count.py` with 13 claimed on the jet-route configuration reports a mismatch.
- One charge of the jet-route configuration multiplied by 1.000001 gives a count of 5 in `indep_count.py`, and `verify_arb.py` rejects the certificate.
- `lower_bound_check.py` asked to prove 12 stops with an error.
- The jet-route recipe with a linear term of -1e-10 in place of -3e-10 gave 9 in the review's own count. That control is not in the tests.

## How narrow it is

[F] Both configurations are needles. `explore/four_eleven/width.py` multiplies every coordinate and charge by a random
factor and counts again, numerically (saved in `results/four-eleven/width_flat.out` and `width_jet.out`).

| Relative change | Flat route, 12 trials | Jet route, 12 trials |
|---|---|---|
| 1e-8 | 11 in all 12 | 11 in 1 |
| 1e-7 | 11 in all 12 | 11 in 2 |
| 1e-6 | 11 in 7 | 11 in none |
| 1e-5 | 11 in 2 | 11 in none |
| 1e-4 and above | 11 in none | 11 in none |

[I] This is why random search never found it. In the repository's sample of 25,108 random four-charge configurations (`results/SUMMARY.md`), 25,054 have the least possible count of 3.
The table is from one run of the review. The saved outputs of a rerun are in `results/four-eleven/` and can differ in the trials that count as 11.

## How it was found

### The starting point

[F] At the centre of a regular tetrahedron of equal charges the potential is 4/sqrt(3) + (20 sqrt(3)/27) xyz - (14 sqrt(3)/243) K4 + higher terms, with K4 = x^4 + y^4 + z^4 - 3(x^2 y^2 + y^2 z^2 + z^2 x^2). The linear and quadratic terms vanish.

[F] That centre is a balance point of multiplicity exactly 11. The dimension count that gives the multiplicity runs 7, 10, 11, 11 and stays there, for the true series and not only for a truncation. `explore/four_eleven/milnor.py` computes it (`results/four-eleven/milnor.out`).

[F] Numerically the centre is the only balance point of the regular tetrahedron, from 20,000 starts. So 11 is the most this one point could give, and the question was whether all 11 can be made real at once.

### Jet route

`explore/four_eleven/jetsolve.py` solves 16 equations for the sites and charges that give a chosen gradient, Hessian and cubic part at the origin. The values that work are a gradient of -3e-10 in each coordinate, a Hessian with diagonal 2e-6, 2e-6 and -4e-6, and cubic coefficients of 1.0e-3, 1.2e-3 and 0.8e-3. `explore/four_eleven/make_cert.py` rounds the result to 12-digit fractions and recounts. The command is in `results/four-eleven/run_programs.sh`, and its output says the configuration is reproduced exactly.

[I] Along each axis the reduced equation has three roots at well-separated scales. Two more points sit near the centre. Three axes times three, plus two, is 11. A first try with a different Hessian gave 7, as the reality condition for the outer pair predicts.

### Flat route

Call a point *flat* when both the gradient and the whole Hessian of the potential vanish there. That is 8 conditions.

[F] **Closed form.** Every four-charge configuration with a flat point at the origin comes from a unit vector n in four dimensions with all entries positive. Put t = n / sqrt(1 - n^2) entry by entry. Then site i lies at distance t_i from the origin, sites i and j have dot product -t_i^2 t_j^2, and charge i is t_i^3 / (1 + t_i^2). Rotation and overall scale are free. The regular tetrahedron is n = (1, 1, 1, 1)/2.

The derivation is in the header of `explore/four_eleven/flat_param.py`. The same program checks it on random n: gradient and Hessian come out below 1e-13.

[I] In other words, reflect the four sites in a sphere about the flat point. The four images form a tetrahedron whose four altitudes meet at the flat point.

**What a flat point is worth.** Near a flat point the potential starts with a harmonic cubic C. Moving the charges slightly adds a linear term l and a quadratic term with a traceless matrix Q. The balance points near the flat point are then the real solutions of grad C(x) + Qx + l = 0. These are three quadrics, so there are at most 8 solutions.

[I] The equation keeps its form when x, Q and l are scaled by s, s and s^2. So if the cubic model has 8 real solutions, the true potential has them too once the scale is small enough.

[F] All 8 can be real. `explore/four_eleven/cubic_scan2.py` draws flat configurations and pairs (Q, l) at random and solves for all 8 roots.

| Where n is drawn | Samples | With all 8 real |
|---|---|---|
| near the regular tetrahedron, spread 0.05 | 200,000 | 2,197 |
| spread 0.15 | 200,000 | 641 |
| the whole flat family | 200,000 | 139 |

[I] **The rest of the configuration.** The gradient of a cubic is a quadratic map. Its local algebra has graded pieces of size 1, 3, 3, 1 and pairs degree k with degree 3 - k, so its signature is 0. The local degree of a flat point is therefore 0. The other balance points must then satisfy index-2 minus index-1 = 3 on their own, so there are at least 3 of them.

[F] So an all-real unfolding of one flat point gives at least 8 + 3 = 11. `explore/four_eleven/realise.py` turns a hit of the cubic model into an actual configuration by Newton steps on the 8 conditions, and `rationalise.py` writes it as 14-digit fractions. The first hit tried gave 11 at three scales.

[F] In the certified one, 4 points of each index lie near the flat point and 3 of index 2 lie further out. The chain from the scan to the configuration is reproduced exactly by `results/four-eleven/run_programs.sh`.

## Walls

[F] In both configurations all six pairs of charges share a wall, and one pair shares two. That is seven walls, and 2 x 7 - 3 = 11. The trace is numerical (`explore/four_eleven/walls2.py`, outputs `walls2_jet.out` and `walls2_flat.out`). It held when the step was cut from 2% to 0.5%.

## What is published

- [P] arXiv:2607.27197, 29 July 2026: five charges with at least 24 balance points, and 3 + 2m charges with at least 4 + 20m. It does not treat four charges.
- [P] arXiv:2607.28785: three positive charges have at most 6.
- [P] A 2022 paper gives nine balance points for four charges in a plane, which meets the bound.
- [P] arXiv:2501.05315 treats the centre of the regular tetrahedron as degenerate and says only that it splits into three saddles.

[F] Searches on 15 November 2026 found no report of more than 9 for four charges, and none of the multiplicity 11. That is "no earlier report found". It is not proof that the result is new.

[I] If the result holds, four is the smallest number of charges at which Maxwell's bound is known to fail. For three charges the bound is 4 and the proved limit is 6.

## Not verified

- Nothing here shows that 11 is the largest count for four charges.
- The two derivative bounds that `indep_count.py` relies on are proved by hand, not by machine.
- That the regular tetrahedron has only one balance point is numerical.
- The degree argument is a sketch. The certified counts do not depend on it.
- The wall traces and the width table are numerical.
- The published papers were read from their abstracts and cached texts. No author was contacted.

## Open

1. A flat configuration with five or more other balance points would give 13 by the same unfolding. `explore/four_eleven/survey_flat.py` counts them across the flat family. In the saved run (200 random flat configurations, 155 consistent with the Morse relation) every one had exactly 3 other balance points, so nothing above 11 came from this route.
2. Flat points whose cubic is itself degenerate have multiplicity above 8.
3. A sturdier example, one that survives changes of 1e-3, would be easier for others to check.
4. For more charges, one flat point is worth 8 balance points where a plain one is worth 1. Whether that raises the growth per charge is untested.

## Rerun

With python-flint, NumPy, SciPy, SymPy and mpmath installed, from the repository root:

```
python3 -O verify_arb.py certificates/n4-eleven-flat.json
python3 -I -O indep_count.py certificates/n4-eleven-flat.json
python3 -I -O lower_bound_check.py certificates/n4-eleven-flat.json 11
python3 -I -O lower_bound_check.py certificates/n4-eleven-jet.json 11
python3 tests/test_four_eleven.py
sh results/four-eleven/run_programs.sh      # rebuilds both configurations, about 15 minutes
```

`results/four-eleven/run_checks.sh` runs all the checks and the whole test set, and its output is `results/four-eleven/checks.log`.
