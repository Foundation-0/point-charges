"""Tests for lower_bound_hp.py, the lower-bound prover for configurations with many scales.

Default case (a few seconds): on the 9-charge configuration (flat to order 4) it proves at least 26 balance points
under python -O, it stops with an error when asked for 27, and it fails when one charge is changed by one part in a
million (the candidates then no longer hold). Optional second case, off by default because it takes several minutes:
  python3 tests/test_lower_bound_hp.py --slow
proves at least 87 on the 21-charge configuration (flat to order 7).
Needs NumPy, mpmath and python-flint. Exit status 0 only if every case behaves."""
import json
import os
import subprocess
import sys
import tempfile
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCRIPT = os.path.join(ROOT, 'lower_bound_hp.py')
CFG9 = os.path.join(ROOT, 'explore', 'configs', 'growth-flat-d5-n9.json')
CAND9 = os.path.join(ROOT, 'explore', 'growth', 'data', 'n9_candidates.npy')
CFG21 = os.path.join(ROOT, 'explore', 'configs', 'growth-flat-d8-n21.json')
CAND21 = os.path.join(ROOT, 'explore', 'growth', 'data', 'n21_candidates.npy')
ok = True


def run(cfg, cand, n):
    return subprocess.run([sys.executable, '-I', '-O', SCRIPT, cfg, cand, str(n)], capture_output=True, text=True)


def report(label, good, r):
    global ok
    ok &= good
    tail = ((r.stdout + r.stderr).strip().splitlines() or [''])[-1]
    print('%-58s exit %d %s | %s' % (label, r.returncode, 'OK' if good else 'WRONG', tail[:60]))


r = run(CFG9, CAND9, 26)
report('9 charges: proves at least 26', r.returncode == 0 and 'PROVED: at least 26' in r.stdout and 'index 2: 17, index 1: 9' in r.stdout, r)
r = run(CFG9, CAND9, 27)
report('9 charges: asked for 27 stops with an error', r.returncode != 0 and 'expected at least 27' in r.stderr, r)
cfg = json.load(open(CFG9))
cfg['charges'][0] = str(Fraction(cfg['charges'][0]) * (1 + Fraction(1, 10 ** 6)))
with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
    json.dump(cfg, f)
try:
    r = run(f.name, CAND9, 26)
finally:
    os.unlink(f.name)
report('9 charges: one charge changed by 1e-6, 26 not proved', r.returncode != 0 and 'PROVED: at least 26' not in r.stdout, r)
if '--slow' in sys.argv:
    r = run(CFG21, CAND21, 87)
    report('21 charges: proves at least 87 (slow)', r.returncode == 0 and 'PROVED: at least 87' in r.stdout, r)
print('ALL CASES BEHAVED' if ok else 'FAILURE')
sys.exit(0 if ok else 1)
