#!/bin/sh
# Checks of the three four-charge certificates with 11 balance points, and the whole test set.
# Run it from the repository root: sh results/four-eleven/run_checks.sh > results/four-eleven/checks.log 2>&1
# Needs python-flint, NumPy and SciPy. Total time: about 20 minutes.
PY=${PY:-python3}
run() { echo; echo "\$ $*" | sed "s|^\$ $PY|\$ python3|"; "$@"; echo "[exit status $?]"; }
for n in n4-eleven-jet n4-eleven-flat n4-eleven-short; do
  run $PY verify_arb.py certificates/$n.json
  run $PY -O verify_arb.py certificates/$n.json
  run $PY -I -O indep_count.py certificates/$n.json
  run $PY -I -O lower_bound_check.py certificates/$n.json 11
done
run $PY -I -O lower_bound_check.py certificates/n4-eleven-flat.json 12
run $PY tests/test_four_eleven.py
run $PY tests/test_indep_count.py
run $PY tests/test_verifier_rejects.py
run $PY tests/test_verifier_rejects_O.py
run $PY tests/test_ia.py
run $PY tests/test_known.py
run $PY -O verify_arb.py certificates/*.json
