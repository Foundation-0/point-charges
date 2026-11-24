# 15 November 2026: how the best known counts grow with the number of charges, and a search for three and four charges

This note adds to issue 5 (eight or more charges), issue 4 (three charges) and issue 3 (four charges). A later note, `2026-11-15-four-charges-eleven.md`, gives four charges with 11 equilibria, above the figures of 9 quoted here for four charges. It follows the
head `a2fd159`. Marks: **[C]** certified by the programs in this repository, **[N]** numerical, not certified,
**[I]** inference.

## What the new certificates show

- [C] Eight charges can have exactly 49 equilibria and nine charges exactly 64. Both equal Maxwell's bound (n-1)^2 and
  neither beats it. The certificates are `certificates/n8-c7-plus-centre.json` and `certificates/n9-abk-third-pair.json`.
  Each is accepted by `verify_arb.py` and by `indep_count.py`.
- [C] The eight-charge configuration is the seven-charge one with a charge of 2^-86 added at the centre of the triangle.
  This is the construction of Proposition 2(b) in NOTE.md, which is only a sketch. This certificate checks one case of it.
  The nine-charge configuration adds a third axial pair at distance 2^-21, with the charge from the Arathoon, Ball and
  Kvalheim formula (arXiv:2607.27197, Proposition 1).
- [N] The nine-charge window is narrow. Changing the last charge by half a percent gives 52 or 56 instead of 64.
- [C] A second family with four charges and exactly 9 equilibria: `certificates/n4-second-family.json`. It is an
  irregular tetrahedron with no symmetry and charges near 1, 1.90, 1.20 and 8.30. All nine equilibria lie within about 0.05
  of one point. The first four-charge example in this repository is a triangle with a charge near its centre, so this is a
  different family. It does not beat 9. It was found by a numerical search (below) and rounded to nine decimals.

## Growth per charge

Let M(n) be the most non-degenerate equilibria known for n positive charges. `explore/growth.py` takes the best of the
known rules, and `results/growth_rules.txt` is its output for n up to 40:

1. ABK, pairs on the axis: M(3 + 2m) is at least 4 + 20m (ABK Proposition 1, proved there).
2. A tiny charge at the centre: M(4 + 2m) is at least 9 + 20m. This is NOTE.md Proposition 2(b), a sketch. It is certified here for m = 0, 1, 2 (n = 4, 6, 8).
3. Adding one charge on an index-1 point raises the count by at least 3 (NOTE.md Proposition 2(a), a sketch).
4. Clustering: M(a + b - 1) is at least M(a) + M(b). Shrink a configuration of a charges to a point and use it in place of one charge of a configuration of b charges. [I] This is a written argument, not a formal proof.

[I] Every one of these grows by about 10 equilibria per charge. The best of them is 10n - 26 for odd n and 10n - 31 for even n.
Maxwell's bound grows with the square.

| charges | best known | Maxwell (n-1)^2 | |
|---|---|---|---|
| 5, 6, 7 | 24, 29, 44 | 16, 25, 36 | beaten |
| 8, 9 | 49, 64 | 49, 64 | met, not beaten |
| 10 | 69 | 81 | below |
| 20 | 169 | 361 | below |
| 40 | 369 | 1521 | below |

[I] Note added later: the statement that every known construction grows by about 10 equilibria per charge is about the constructions known before `results/2026-11-15-growth-n-three-halves.md`. That note argues, with a proof sketch, that the growth can be like n^(3/2). It beats no count yet.

[I] So with what is known today, the bound is beaten only for 5, 6 and 7 charges and met for 8 and 9.
From 10 charges up, no known construction reaches it. Clustering cannot change the rate, because it never gives more per charge than
the best single configuration it uses (superadditivity and Fekete's lemma).

[I] Whether the true maximum grows faster than a straight line is open. One route to test: how many real critical points can a
harmonic polynomial of degree d on R^3 have? If that grows faster than d squared, nested constructions might beat a straight line.

No earlier report of the explicit counts 49 and 64 was found. The best published count for eight charges found is 25 (numerical).
Seven and nine charges follow from the ABK formula, but no explicit configuration with an exact count was found. This is
"no earlier report found" after reading the arXiv listing and three closely related papers.

## Search for three and four charges (numerical)

A push method was tried on the open cases. It tracks all complex solutions of the equilibrium equations and moves the
charges and sites so that complex pairs approach the real axis and become real equilibria, while keeping real ones apart.
The search code is not part of this branch. These counts are reported, not checkable here. Only the one configuration
that was then certified (`certificates/n4-second-family.json`) is a result of this repository.

| charges | pushes | steps | best count reached |
|---|---|---|---|
| 3 | 311 | 46,063 | 4 in 197 pushes, 2 in 114 |
| 4 | 130 | 16,495 | 9 in 17, 7 in 23, 5 in 56, 3 in 34 |

- [N] The push found nothing above 4 for three charges and nothing above 9 for four charges. For four charges this is now known to be a limit of the push, not of the problem: a different method found 11 (`results/2026-11-15-four-charges-eleven.md`). For three charges the earlier searches in `results/SUMMARY.md` agree, and none of this is a proof.
- [N] For three charges, in 333 of 344 configurations that reach 4, no complex pair lies in the all-positive region. Where one does, the nearest stays far from real. So a fifth point has no nearby source in the configurations tried. [I] A route to a proof: show that at count 4 no complex solution lies in the all-positive region.
- [N] For four charges, at 9 there are still one or two candidate complex pairs, and one comes close. Moving toward it first merges two existing equilibria, so the count falls from 9 to 7. When that merge is forbidden, the push stops. [I] A tenth point is close and blocked by a merge. That is the pattern expected near a configuration where the potential is flat to second order at an equilibrium. The five-charge counterexample was built around such a point.
- [N] The complex solution counts of the system were 68 for three charges (34 points) and 420 for four charges (210 points), by two solvers each. They are not proved.

## Not verified

- The push results above, and the complex counts 68 and 420.
- The clustering argument and the claim about the growth rate.
- Proposition 2 of NOTE.md for general m. Only m = 0, 1, 2 are certified here.
