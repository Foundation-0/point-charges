"""Negative controls under python -O: the independent verifier must still reject corrupted certificates.

Python removes `assert` statements when run with -O, so a verifier whose proof steps are asserts prints VERIFIED
for a corrupted certificate. Each case below runs `python -O verify_arb.py` in a fresh process and requires a
nonzero exit status and no VERIFIED line. The unchanged certificate must still be accepted under -O.
Usage: python3 tests/test_verifier_rejects_O.py   (exit status 0 only if every case behaves)
"""
import copy
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
VERIFIER = os.path.join(ROOT, 'verify_arb.py')
base = json.load(open(os.path.join(ROOT, 'certificates', 'n4-triangle-plus-centre.json')))


def run(cert, optimize):
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
        json.dump(cert, f)
    try:
        cmd = [sys.executable] + (['-O'] if optimize else []) + [VERIFIER, f.name]
        return subprocess.run(cmd, capture_output=True, text=True)
    finally:
        os.unlink(f.name)


def accepted(r):
    return r.returncode == 0 and 'VERIFIED' in r.stdout


cases = []
c = copy.deepcopy(base); c['zeros'].pop(3); c['count'] -= 1
cases.append(('one zero removed', c))
c = copy.deepcopy(base); c['zeros'].append(copy.deepcopy(c['zeros'][0])); c['count'] += 1
cases.append(('one zero duplicated', c))
c = copy.deepcopy(base); c['zeros'][0]['index'] = 3 - c['zeros'][0]['index']
cases.append(('wrong Morse index', c))
c = copy.deepcopy(base); c['charges'][3] = '1/32'
cases.append(('centre charge changed to 1/32', c))
c = copy.deepcopy(base); c['count'] += 1
cases.append(('claimed count one too high', c))
c = copy.deepcopy(base); c['zeros'][0]['box'][0][0] = c['zeros'][0]['box'][1][0]
cases.append(('zero box squeezed to a plane', c))

ok = True
for optimize in (False, True):
    tag = 'python -O' if optimize else 'python   '
    r = run(base, optimize)
    good = accepted(r)
    ok &= good
    print('%s control (unchanged certificate): %s %s' % (tag, 'accepted' if good else 'NOT ACCEPTED', 'OK' if good else 'WRONG'))
    for label, cert in cases:
        r = run(cert, optimize)
        good = (r.returncode != 0) and ('VERIFIED' not in r.stdout)
        ok &= good
        last = (r.stderr.strip().splitlines() or [''])[-1]
        print('%s %-32s %s %s | %s' % (tag, label, 'rejected' if r.returncode != 0 else 'ACCEPTED', 'OK' if good else 'WRONG', last[:70]))
print('ALL NEGATIVE CONTROLS PASSED UNDER -O' if ok else 'FAILURE')
sys.exit(0 if ok else 1)
