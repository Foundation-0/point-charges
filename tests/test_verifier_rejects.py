"""Negative controls: the independent verifier must reject corrupted certificates."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import json, copy, tempfile
import verify_arb

base = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'certificates', 'n4-triangle-plus-centre.json')))

def rejects(cert, label):
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
        json.dump(cert, f)
    try:
        verify_arb.verify(f.name)
    except (AssertionError, RuntimeError) as e:
        print('rejected as expected (%s): %s' % (label, str(e)[:80]))
        return True
    finally:
        os.unlink(f.name)
    print('NOT REJECTED:', label)
    return False

ok = True
c = copy.deepcopy(base); c['zeros'].pop(3); c['count'] -= 1
ok &= rejects(c, 'one zero removed')
c = copy.deepcopy(base); c['zeros'].append(copy.deepcopy(c['zeros'][0])); c['count'] += 1
ok &= rejects(c, 'one zero duplicated')
c = copy.deepcopy(base); c['zeros'][0]['index'] = 3 - c['zeros'][0]['index']
ok &= rejects(c, 'wrong Morse index')
c = copy.deepcopy(base); c['charges'][3] = '1/32'
ok &= rejects(c, 'centre charge changed to 1/32')
print('ALL NEGATIVE CONTROLS PASSED' if ok else 'FAILURE')
sys.exit(0 if ok else 1)
