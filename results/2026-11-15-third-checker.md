# 15 November 2026: a third checker for the equilibrium counts, `indep_count.py`

This note describes `indep_count.py`, a program that proves the number of equilibria of positive point charges without
reading any certificate box. It follows the head `a2fd159` and adds to issues 2 and 5. Marks: **[C]** checked by a program in
this repository, **[I]** inference.

## What it proves

For sites a_i and positive charges q_i, it proves the exact number of zeros of grad V in R^3 minus the sites, that every
zero is non-degenerate, and the Morse index of each. It reads only `sites` and `charges` from the file. The `count`
field is compared with the result at the end and changes nothing in the proof. Exit status 0 means the counts agree,
2 means they differ, 3 means a proof step failed.

[C] It accepts all eight committed certificates. The saved run is `results/indep_count_all.log`.

| certificate | zeros (index 1 / index 2) | boxes in `indep_count.py` | boxes in `verify_arb.py` |
|---|---|---|---|
| n3-triangle | 4 (1 / 3) | 1,801 | 3,563 |
| n4-triangle-plus-centre | 9 (3 / 6) | 6,009 | 15,379 |
| n4-second-family | 9 (3 / 6) | 23,065 | 28,983 |
| n5-abk-bipyramid | 24 (10 / 14) | 39,617 | 104,271 |
| n6-abk-plus-centre | 29 (12 / 17) | 58,561 | 152,925 |
| n7-abk-second-pair | 44 (19 / 25) | 128,817 | 394,265 |
| n8-c7-plus-centre | 49 (21 / 28) | 162,225 | 480,647 |
| n9-abk-third-pair | 64 (28 / 36) | 279,385 | 915,205 |

The `verify_arb.py` counts are from `results/verify_arb_fixed_all.log`, a run of the fixed verifier under `python3 -O`.

## How the method differs from `verify_arb.py`

| | `verify_arb.py` | `indep_count.py` |
|---|---|---|
| Where zeros come from | the boxes listed in the certificate | its own Newton search, never the certificate boxes |
| Existence and uniqueness of a zero | Krawczyk test on a box | contraction on a ball, using a bound on the third derivative |
| Morse index | Gershgorin discs and Sylvester's law of inertia | sign of det H, because trace H = 0 (V is harmonic) |
| Far field | not used | a half-space lemma in exact rationals |
| Near each charge | a dominance test on boxes | one exclusion ball per charge, checked in exact rationals |
| Elsewhere | natural and mean-value enclosures | two Taylor tests at box centres |
| Subdivision | bisection | octree |

The two derivative bounds used by the Taylor tests are ||D^2(1/r)|| = 2/r^3 and ||D^3(1/r)[v,.,.]|| at most 6/r^4. They were
derived by hand. `tests/check_lemmas.py` checks the closed forms against finite differences and samples the bounds
(saved in `results/check_lemmas.txt`). That is a numerical check, not a proof.

Both programs use python-flint for ball arithmetic. Every inequality in `indep_count.py` is evaluated with ball arithmetic
or exact rationals. Floats only propose candidates and approximate inverses, and each float that enters a proof is first
converted exactly to a rational. Proof steps raise errors, never `assert`, so the results are the same under `python3 -O`.

[I] The count of index-2 points minus index-1 points equals the number of charges minus one in every row. The program does not use this relation,
so it is a further check. It follows from Morse theory (NOTE.md section 1).

## Tests

`tests/test_indep_count.py` runs the program under `-O` on corrupted inputs. A claimed count that is too high, a changed charge and a negative charge
are rejected, with exit status 2 or 3. Garbage zero boxes change nothing, because the boxes are never read.

## How the new certificates were made

The configurations are in `explore/configs/`. Each file holds exact sites and charges. For each, `indep_count.py` was run with
an output file, which lists a ball around each zero. `explore/to_repo_format.py` then wrote the certificate in this repository's format:
for each zero it takes a cube around the ball centre, halves it until `verify_arb.krawczyk` accepts it, and computes the tight enclosure and
the index with the code of `verify_arb.py`. The boxes are only a choice. The proof is the separate run of `verify_arb.py` on the written file, and the run of `indep_count.py`, which ignores the boxes.

```
python3 -I -O indep_count.py explore/configs/n8-c7-plus-centre.json /tmp/n8.json
python3 explore/to_repo_format.py explore/configs/n8-c7-plus-centre.json /tmp/n8.json n8-c7-plus-centre certificates/n8-c7-plus-centre.json
python3 -O verify_arb.py certificates/n8-c7-plus-centre.json
```

## Rerun

```
python3 -I -O indep_count.py certificates/n7-abk-second-pair.json     # about a minute; n9 takes about three minutes
python3 tests/test_indep_count.py                                      # about 15 s
python3 tests/check_lemmas.py
```

Needs python-flint and NumPy.

## Not verified

- The two derivative bounds are proved on paper, not by machine. The proofs, and the proof of each step of the count, are written out in `results/2026-11-15-count-lemmas.md`.
- `indep_count.py` was written before its author opened `certify.py`, `ia.py` and `verify_arb.py`, but after reading NOTE.md's description of their method. It is a separately written program, not a blind one.
- Nobody outside the project has run it.
