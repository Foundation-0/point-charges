# point-charges

*A [physics0](https://github.com/Foundation-0/physics0) project of Foundation0.*

## Maxwell's point charges

How many places of perfect balance can a handful of electric charges have? Put positive charges at fixed points
in space. A test charge feels no force at the equilibria of the field, the critical points of the potential. In
1873 James Clerk Maxwell asserted that $n$ charges have at most $(n-1)^2$ such points (when none is degenerate).
In July 2026 a five-charge counterexample with at least 24 equilibria was published (Arathoon, Ball and
Kvalheim, with the idea suggested by an OpenAI model). The three- and four-charge cases are open.

This study adds four computer-certified results and evidence on the open cases. The full write-up is in
[NOTE.md](NOTE.md).

| charges | Maxwell's bound | this study (exact count, all non-degenerate) | certificate |
|---|---|---|---|
| 3 | 4 | 4 (sanity check: the equilateral triangle) | [n3-triangle.json](certificates/n3-triangle.json) |
| 4 | 9 | **9**: the bound is reached (EFO 2026 knew of no example for $n>3$) | [n4-triangle-plus-centre.json](certificates/n4-triangle-plus-centre.json) |
| 5 | 16 | **24**: an explicit ABK counterexample, every equilibrium enclosed | [n5-abk-bipyramid.json](certificates/n5-abk-bipyramid.json) |
| 6 | 25 | **29**: Maxwell's bound fails for six charges too | [n6-abk-plus-centre.json](certificates/n6-abk-plus-centre.json) |
| 7 | 36 | **44** †: the second ABK iteration, completely counted | [n7-abk-second-pair.json](certificates/n7-abk-second-pair.json) |

† Seven charges: proved by the independent verifier only for now (float certifier in issue #12).

Evidence, not proof: no three-charge configuration with more than 4 equilibria, and no four-charge
configuration with more than 9, turned up in the searches described in NOTE.md §3.

## Check it yourself

```
pip install -r requirements.txt
python3 make_certs.py                      # interval-arithmetic certificates for 3-6 charges, about 5 s
python3 explore/second_pair.py             # the 7-charge certificate (Arb Krawczyk boxes), seconds
python3 verify_arb.py certificates/*.json  # independent re-check with Arb balls and exact rationals, about 2 min
python3 tests/test_ia.py                   # the float interval library against exact rationals
python3 tests/test_verifier_rejects.py     # corrupted certificates are rejected
```

Each certificate lists the sites and charges as exact fractions, and gives for every equilibrium a box that
contains exactly one, its Morse index, and a tight enclosure. For 3 to 6 charges two verifiers that share
no code prove the counts. For 7 charges the independent verifier proves the count from the boxes alone;
the float certifier does not yet finish on that configuration (issue #12). In every case the independent
verifier reads only the configuration and the boxes, and re-proves both that each box holds a
non-degenerate equilibrium and that nothing exists outside the boxes.

## Files

| path | what |
|---|---|
| `NOTE.md` | statements, constructions, proof sketch of the insertion lemma, method, evidence, references |
| `certify.py`, `ia.py` | primary certifier: outward-rounded float64 intervals, bisection, Krawczyk test, Gershgorin indices |
| `verify_arb.py` | independent verifier: python-flint Arb balls and exact rationals |
| `make_certs.py` | the 3- to 6-charge configurations, as exact rationals |
| `explore/second_pair.py` | the 7-charge configuration, its numerical search and Arb certificate builder |
| `charges.py` | fast numerical equilibrium finder used for exploration |
| `export_page_data.py` | writes the data file for the Foundation0 page from the certificates |
| `explore/` | the searches and local models behind NOTE.md §3 |
| `results/` | logs and summaries of those runs |
| `tests/` | interval arithmetic soundness, numerical reproductions, negative controls |

## Status

Done by Claude (an AI model by Anthropic) on 18 October 2026, at the request of Foundation0 founding member
Ahraz Arifuddin. Reviewed by an independent AI agent before merge (PR #1); not reviewed by a person. Novelty
was checked by web search only. Corrections and pointers to earlier work are welcome as issues.

This study started in [Foundation-0/physics0](https://github.com/Foundation-0/physics0), the hub for Foundation0's
physics work. It moved here with its history once it outgrew the hub. The open problems it leaves are tracked
as issues in this repository.

## Licence

Copyright 2026 Foundation0 contributors.

Code (`*.py`, including `explore/` and `tests/`): [Apache License 2.0](LICENSE).
Text and data (`NOTE.md`, this README, `certificates/` and `results/`):
[Creative Commons Attribution 4.0 International](LICENSE-CC-BY-4.0) (CC BY 4.0).
