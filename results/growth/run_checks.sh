#!/bin/sh
# Checks for the n^(3/2) note and the existing test set. Run from the repository root:
#   sh results/growth/run_checks.sh > results/growth/checks.log 2>&1
# Needs NumPy, SciPy, mpmath and python-flint. Total time: about 10 minutes.
PY=${PY:-python3}
run() { echo; echo "\$ $*" | sed "s|^\$ $PY|\$ python3|"; "$@"; echo "[exit status $?]"; }
run $PY -I -O lower_bound_hp.py explore/configs/growth-flat-d5-n9.json explore/growth/data/n9_candidates.npy 26
run $PY -I -O lower_bound_hp.py explore/configs/growth-flat-d8-n21.json explore/growth/data/n21_candidates.npy 87
run $PY explore/growth/levels_check.py
run $PY -O explore/growth/stack_check.py 6 maxcoef 200
run $PY tests/test_lower_bound_hp.py --slow
run $PY tests/test_four_eleven.py
run $PY tests/test_indep_count.py
run $PY tests/test_verifier_rejects.py
run $PY tests/test_verifier_rejects_O.py
run $PY tests/test_ia.py
run $PY tests/test_known.py
run $PY -O verify_arb.py certificates/n4-eleven-short.json certificates/n4-eleven-jet.json certificates/n4-second-family.json
