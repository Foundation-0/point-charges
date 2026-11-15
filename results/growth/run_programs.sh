#!/bin/sh
# Runs the programs under explore/growth/ that support results/2026-11-15-growth-n-three-halves.md.
# Run it from the repository root. Scratch files go to .work/ (not committed).
# The output of each program is saved in results/growth/<name>.out. Total time: about 15 minutes.
# Needs NumPy, SciPy, mpmath and python-flint.
set -e
PY=${PY:-python3}
E=../explore/growth
O=../results/growth
C=../explore/configs
mkdir -p .work/configs .work/out
cd .work
$PY $E/kelvin_check.py > $O/kelvin_check.out
$PY $E/platonic.py > $O/platonic.out
$PY $E/levels_check.py > $O/levels_check.out
$PY $E/hier.py 0.01 10 > $O/hier_s0.01.out
$PY $E/level.py > $O/level.out
$PY $E/nflat.py 4 6 40 > $O/nflat_d4.out
$PY $E/nflat.py 5 9 40 > $O/nflat_d5.out
$PY $E/nflat.py 6 12 40 > $O/nflat_d6.out
cp $E/data/flat_d5_n9.json $E/data/flat_d8_n21.json configs/
$PY $E/demo_flat.py configs/flat_d5_n9.json 0.02 1/20 demo_d5_n9 > $O/demo_flat_d5_n9.out
$PY $E/demo_flat.py configs/flat_d8_n21.json 0.01 1/100 demo_d8_n21 nofinder > $O/demo_flat_d8_n21.out
$PY - <<'PYEND' > $O/demo_flat_reproduces_configs.out
import json
from fractions import Fraction as F
for mine, kept in (("demo_d5_n9", "growth-flat-d5-n9"), ("demo_d8_n21", "growth-flat-d8-n21")):
    a = json.load(open("configs/%s.json" % mine)); b = json.load(open("../explore/configs/%s.json" % kept))
    same = a["sites"] == b["sites"] and a["charges"] == b["charges"]
    flat = lambda d: [x for s in d["sites"] for x in s] + d["charges"]
    diff = max(abs(F(x) - F(y)) for x, y in zip(flat(a), flat(b)))
    print(kept, "regenerated sites and charges identical:", same, "; largest difference: %.1e" % float(diff))
PYEND
$PY $E/check_taylor.py $C/growth-flat-d5-n9.json 9 > $O/check_taylor_d5_n9.out
$PY $E/dump_taylor_points.py $C/growth-flat-d8-n21.json 18 -16 0.06 6 taylor_n21.npy > $O/dump_taylor_points_d8_n21.out
