"""Tests for indep_count.py, the third checker: it accepts committed certificates and rejects corrupted ones.

Every run is `python -O indep_count.py` in a fresh process. Exit status 0 means the independent count agrees with the
file's "count", 2 means it differs, 3 means a proof step failed. Each corrupted configuration must give a nonzero status.
Usage: python3 tests/test_indep_count.py   (about 15 s; needs python-flint and NumPy)
"""
import copy
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCRIPT = os.path.join(ROOT, 'indep_count.py')


def load(name):
    return json.load(open(os.path.join(ROOT, 'certificates', name + '.json')))


def run(cert):
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
        json.dump(cert, f)
    try:
        r = subprocess.run([sys.executable, '-O', SCRIPT, f.name], capture_output=True, text=True)
    finally:
        os.unlink(f.name)
    return r


ok = True


def expect(label, cert, want):
    global ok
    r = run(cert)
    good = r.returncode == want if want == 0 else (r.returncode != 0)
    ok &= good
    last = [ln for ln in r.stdout.splitlines() if ln.startswith(('AGREES', 'MISMATCH', 'PROOF FAILURE'))]
    print('%-46s exit %d %s | %s' % (label, r.returncode, 'OK' if good else 'WRONG', (last or [''])[0][:70]))


n4 = load('n4-triangle-plus-centre')
expect('n4 unchanged', n4, 0)
c = copy.deepcopy(n4); c['count'] = 10
expect('n4 claimed count 10 (true count is 9)', c, 2)
c = copy.deepcopy(n4)
for z in c['zeros']:
    z['box'] = [[0, 0, 0], [1e-9, 1e-9, 1e-9]]; z['index'] = 1
expect('n4 boxes replaced by garbage (never read)', c, 0)
c = copy.deepcopy(n4); c['charges'][3] = '1/32'
expect('n4 centre charge 1/32 (count changes)', c, 2)
c = copy.deepcopy(n4); c['charges'][3] = '-1/1024'
expect('n4 negative charge', c, 3)
n5 = load('n5-abk-bipyramid')
expect('n5 unchanged (24 zeros)', n5, 0)
c = copy.deepcopy(n5); c['charges'][3] = c['charges'][4] = '441/131072'
expect('n5 axial charges changed (count changes)', c, 2)
print('ALL CASES BEHAVED' if ok else 'FAILURE')
sys.exit(0 if ok else 1)
