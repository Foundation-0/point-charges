# Runs behind NOTE.md

All runs on 18 October 2026, 4 CPU cores. Counts are numerical (multi-start Newton with a Morse-relation check),
except the certificates, which are rigorous.

| run | script | configurations | outcome |
|---|---|---|---|
| certificates | `make_certs.py` | 4 | exactly 4, 9, 24, 29 equilibria, all non-degenerate (`make_certs.log`) |
| independent check | `verify_arb.py` | 4 | all four re-verified (`verify_arb.log`) |
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
