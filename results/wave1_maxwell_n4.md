# Wave-1 search for >9 equilibria with four charges (issue #3)

Script: `explore/wave1_maxwell_n4.py` (seed 20261018). Screening: `charges.count_fast`;
only Morse-consistent counts (m2-m1=3, no other indices) are tallied as counts.
Up to 40 de-duplicated configurations screening >=9 were re-checked with fresh `count_fast` seeds and `charges.summary`; all re-checked as 9.
Morse-inconsistent screens include degenerate tri_interior grid points where the 4th site coincides with a triangle vertex; they are excluded from all counts.
Runtime 470.6s. Sweep configurations: 19336. Hill-climb evaluations: 4509.

**Maximum Morse-consistent count found: 9.** Configurations verified >9 by the slow solver: 0.

## Count distribution (Morse-consistent screens only)

| family | n=3 | n=5 | n=7 | n=9 | inconsistent |
|---|---|---|---|---|---|
| c3v_axial | 416 | 598 | 0 | 175 | 0 |
| insertion | 7 | 32 | 3 | 10 | 0 |
| near_coplanar | 1680 | 8 | 0 | 0 | 0 |
| random | 4000 | 0 | 0 | 0 | 0 |
| symmetric_disphenoid | 215 | 11 | 4 | 0 | 0 |
| symmetric_rhombus | 403 | 32 | 15 | 0 | 0 |
| symmetric_tet | 285 | 25 | 0 | 0 | 5 |
| tri_interior | 7999 | 731 | 3 | 0 | 327 |
| two_pair | 1147 | 4 | 1 | 0 | 0 |
| two_pair_rand | 1200 | 0 | 0 | 0 | 0 |
| **total** | 17352 | 1441 | 26 | 185 | 332 |

Hill climbing (9 chains x 500 steps from the certified 9-equilibrium
configuration and the best sweep finds): best count 9.

## Claim labels

- VERIFIED (finite numerical check): the counts above, for the configurations actually run.
- NOT a proof: no covering argument over the 8-parameter (shape+charge) space is made.
  Absence of >9 here is evidence about the searched families, nothing more.

Best configurations and all verified candidates: `results/wave1_maxwell_n4.json`.
