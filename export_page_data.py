"""Write the data file for the Foundation0 page research/maxwell-charges/ from the certificates.

Frame: the triangle's centre c = (1,1,1) at the origin, its axis along z, a vertex on +x, unit circumradius:
  X = ((p - c).e1, (p - c).e2, (p - c).e3) / sqrt(6),  e1 = (2,-1,-1)/sqrt6, e2 = (0,1,-1)/sqrt2, e3 = (1,1,1)/sqrt3.
Equilibrium positions are the midpoints of the certified enclosures (width below 1e-9).
"""
import json, sys
import numpy as np

COMMIT = sys.argv[1] if len(sys.argv) > 1 else 'main'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'equilibria.json'
BASE = 'https://github.com/Foundation-0/point-charges/blob/%s/certificates/' % COMMIT
c = np.ones(3)
E = np.array([[2, -1, -1], [0, 1, -1], [1, 1, 1]], float)
E /= np.linalg.norm(E, axis=1)[:, None]
R = np.sqrt(6)
names = [('n3', 'n3-triangle', 'Three equal charges on a triangle'),
         ('n4', 'n4-triangle-plus-centre', 'The triangle plus a tiny charge at its centre'),
         ('n5', 'n5-abk-bipyramid', 'The triangle plus two small charges on its axis (Arathoon, Ball and Kvalheim)'),
         ('n6', 'n6-abk-plus-centre', 'The five-charge example plus a tiny charge at its centre')]
out = dict(source='Foundation-0/point-charges, certificates at ' + COMMIT,
           frame='centre of the triangle at the origin, its axis along z, a vertex on +x, unit circumradius',
           configs=[])
for key, fn, title in names:
    C = json.load(open('certificates/%s.json' % fn))
    to = lambda p: [float(v) for v in ((np.asarray(p, float) - c) @ E.T) / R]
    eq = []
    for z in C['zeros']:
        mid = [(a + b) / 2 for a, b in zip(*z['enclosure'])]
        w = max(b - a for a, b in zip(*z['enclosure']))
        assert w < 1e-9, w
        eq.append(to(mid) + [z['index']])
    eq.sort(key=lambda e: (round(float(np.linalg.norm(e[:3])), 9), e[3], round(e[0], 9), round(e[1], 9), round(e[2], 9)))
    n = len(C['charges'])
    out['configs'].append(dict(id=key, title=title, n=n, maxwell=(n - 1) ** 2, count=C['count'],
                               m1=C['index_counts'].get('1', 0), m2=C['index_counts'].get('2', 0),
                               sites=[to(p) for p in C['sites_float']], charges=C['charges'],
                               charges_float=C['charges_float'], equilibria=eq, cert=BASE + fn + '.json'))
json.dump(out, open(OUT, 'w'), indent=1)
print('wrote', OUT)
