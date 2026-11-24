"""Tests for the four-charge certificates with 11 balance points (n4-eleven-jet.json, n4-eleven-flat.json and
n4-eleven-short.json).

Every program runs under python -O in a fresh process. Checked here:
  - lower_bound_check.py proves at least 11 on both certificates and stops with an error when asked for 12;
  - indep_count.py rejects a claimed count of 13, and a charge multiplied by 1.000001 (the count falls to 5);
  - verify_arb.py rejects the jet certificate when one charge is changed.
Usage: python3 tests/test_four_eleven.py   (about 2 minutes; needs python-flint and NumPy)
"""
import copy
import json
import os
import subprocess
import sys
import tempfile
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NAMES = ['n4-eleven-jet', 'n4-eleven-flat', 'n4-eleven-short']
ok = True


def load(name):
    return json.load(open(os.path.join(ROOT, 'certificates', name + '.json')))


def run(script, cert, *extra):
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
        json.dump(cert, f)
    try:
        return subprocess.run([sys.executable, '-I', '-O', os.path.join(ROOT, script), f.name] + list(extra),
                              capture_output=True, text=True)
    finally:
        os.unlink(f.name)


def report(label, good, r):
    global ok
    ok &= good
    tail = ((r.stdout + r.stderr).strip().splitlines() or [''])[-1]
    print('%-62s exit %d %s | %s' % (label, r.returncode, 'OK' if good else 'WRONG', tail[:60]))


for name in NAMES:
    cert = load(name)
    r = run('lower_bound_check.py', cert, '11')
    report('%s: lower_bound_check proves 11' % name, r.returncode == 0 and 'PROVED: at least 11' in r.stdout, r)
r = run('lower_bound_check.py', load(NAMES[1]), '12')
report('n4-eleven-flat: lower_bound_check asked for 12 stops', r.returncode != 0 and 'expected at least 12' in r.stderr, r)

jet = load(NAMES[0])
c = copy.deepcopy(jet); c['count'] = 13
r = run('indep_count.py', c)
report('n4-eleven-jet: indep_count, claimed 13', r.returncode == 2, r)
c = copy.deepcopy(jet)
q = Fraction(c['charges'][1]) * (1 + Fraction(1, 10 ** 6)); c['charges'][1] = str(q)
r = run('indep_count.py', c)
report('n4-eleven-jet: indep_count, a charge times 1.000001', r.returncode == 2 and 'independent count is 5' in r.stdout, r)
r = run('verify_arb.py', c)
report('n4-eleven-jet: verify_arb, a charge times 1.000001', r.returncode != 0 and 'VERIFIED' not in r.stdout, r)
print('ALL CASES BEHAVED' if ok else 'FAILURE')
sys.exit(0 if ok else 1)
