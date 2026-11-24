#!/bin/sh
# Runs every program under explore/four_eleven/ that gives the two 11-point certificates.
# Run it from the repository root. Scratch files go to .work/ (not committed).
# The output of each program is saved in results/four-eleven/<name>.out.
# Needs python-flint, NumPy, SciPy, SymPy and mpmath. Total time: about 15 minutes.
set -e
PY=${PY:-python3}
E=../explore/four_eleven
O=../results/four-eleven
mkdir -p .work
cd .work
$PY $E/flat_param.py > $O/flat_param.out
$PY $E/cubic_scan2.py 1 200000 0.05 > $O/cubic_scan2_0.05.out
$PY $E/cubic_scan2.py 1 200000 0.15 > $O/cubic_scan2_0.15.out
$PY $E/cubic_scan2.py 1 200000 9.0 > $O/cubic_scan2_9.0.out
$PY $E/realise.py cubic_scan2_hits_1_0.05.json 4 > $O/realise.out
$PY $E/rationalise.py realise_best.json flat-config.json n4-eleven-flat 11
cmp flat-config.json ../explore/configs/n4-eleven-flat.json && echo "flat configuration reproduced exactly" > $O/rationalise.out
$PY $E/jetsolve.py > $O/jetsolve.out
$PY $E/try11.py > $O/try11.out
$PY $E/make_cert.py 1e-3,1.2e-3,0.8e-3 2e-6,2e-6,0,0,0 -3e-10,-3e-10,-3e-10 n4-eleven-jet 12 > $O/make_cert.out
cmp n4-eleven-jet.json ../explore/configs/n4-eleven-jet.json && echo "jet configuration reproduced exactly" >> $O/make_cert.out
$PY $E/milnor.py > $O/milnor.out
$PY -I -O ../indep_count.py ../explore/configs/n4-eleven-jet.json jet.indep.json > $O/indep_count_jet_config.out
$PY -I -O ../indep_count.py ../explore/configs/n4-eleven-flat.json flat.indep.json > $O/indep_count_flat_config.out
$PY $E/walls2.py jet.indep.json ../explore/configs/n4-eleven-jet.json > $O/walls2_jet.out
$PY $E/walls2.py flat.indep.json ../explore/configs/n4-eleven-flat.json > $O/walls2_flat.out
$PY $E/width.py ../explore/configs/n4-eleven-flat.json flat.indep.json > $O/width_flat.out
$PY $E/width.py ../explore/configs/n4-eleven-jet.json jet.indep.json > $O/width_jet.out
$PY $E/survey_flat.py 1 200 > $O/survey_flat.out
