"""Certify the four explicit configurations (exact dyadic rational data)."""
import numpy as np, time, sys
from certify import run
T = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
c = [1, 1, 1]
t = 0.25
q5 = 541 / 131072            # 0.00412750244140625
q6 = 2.0 ** -44
cfgs = {
    'n3-triangle': (T, [1, 1, 1]),
    'n4-triangle-plus-centre': (T + [c], [1, 1, 1, 2.0 ** -10]),
    'n5-abk-bipyramid': (T + [[1 + t] * 3, [1 - t] * 3], [1, 1, 1, q5, q5]),
    'n6-abk-plus-centre': (T + [[1 + t] * 3, [1 - t] * 3, c], [1, 1, 1, q5, q5, q6]),
}
for name in (sys.argv[1:] or cfgs):
    a, q = cfgs[name]
    t0 = time.time()
    run(np.array(a, float), np.array(q, float), name, out=f'certificates/{name}.json')
    print('  %.1fs' % (time.time() - t0), flush=True)
