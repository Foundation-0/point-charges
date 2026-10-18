# Runs behind NOTE.md

All runs on 18 October 2026, 4 CPU cores. Counts are numerical (multi-start Newton with a Morse-relation check),
except the certificates, which are rigorous.

| run | script | configurations | outcome |
|---|---|---|---|
| certificates | `make_certs.py` | 4 | exactly 4, 9, 24, 29 equilibria, all non-degenerate (`make_certs.log`) |
| seven-charge certificate | `explore/second_pair.py` | 1 (second ABK pair; charge scan over $k/2^{60}$ at offset $1/2048$) | exactly 44 equilibria (19 index 1, 25 index 2), all non-degenerate; 44-count window $k \in [134615, 134812]$, centre $k=134713$ used; larger offsets $1/512$, $1/1024$ reach only 38 numerically. Float certifier does not finish on this configuration (issue #12) |
| independent check | `verify_arb.py` | 5 | all five verified; seven-charge completeness 394,265 boxes (`verify_arb.log`) |
| negative controls | `tests/test_verifier_rejects.py` | 4 corrupted | all rejected |
| C3v coarse scan | `explore/scan_c3v.py` | 899 (h × log q grid) | counts 3: 258, 4: 239, 5: 369, 9: 33; maximum 9 |
| C3v fine scan | `explore/scan_c3v_fine.py` | 2,928 (near the H2-cancelling charge) | counts 3: 1,272, 5: 942, 9: 714; maximum 9 |
| random four-charge sampling | `explore/sample4.py` | 25,108 | Morse-consistent counts 3: 25,054, 5: 53, 7: 1; maximum 7 (`sample4_1.log`) |
| ghost-guided hill climbing | `explore/search.py` | 8 chains × 300 steps | maximum 9 (`search_n4_ghost.json`) |
| four index-1 points, least squares | `explore/four_lsq.py` | 160 restarts | best residual 3.1e-3 in log-charge space, never 0; best runs flatten the tetrahedron (`four_lsq_0.json`) |
| one-charge local model | `explore/local_model.py` | 3,000 random backgrounds | counts 2: 2,119, 4: 881; maximum 4 |
| ABK local model | `explore/local_poly.py` | 150 random (L, Q, C) | finite complex solutions 21 in 143 cases (fewer in 7: solver misses); ABK's coefficients give 21 real |
| polygon + axial pair | `explore/kgon_pair.py` | k = 3..6 | best counts 24 (k=3, n=5), 13 (n=6), 16 (n=7), 19 (n=8) |
| random three-charge sampling | `explore/three.py` | 16,049 | 0 minima: 16,047; 1 minimum: 2; 2 minima: none |
| three-charge shape scan | `explore/three_scan.py`, `explore/three_scan_check.py` (log: `three_scan.log`) | 1,289 triangle shapes | no charge ratio covered twice by the minimum region in the 799 resolved shapes; 3 flagged edge pixels re-solved directly: 0 or 1 minimum (artifacts); 487 mostly flat shapes show no minimum region at grid resolution (`three_scan_summary.txt`) |
| four-charge importance sampling | `explore/sample4_inverse.py` | 2,703 random tetrahedra, 16,179 charge vectors, each with at least one 1-saddle | Morse-consistent counts 5: 16,078, 7: 98, 9: 1; one run failed the Morse check and was discarded; maximum 9 (`sample4_inverse.log`) |
