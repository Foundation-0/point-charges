import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np
from charges import summary
s3 = np.sqrt(3)
tri = np.array([[1, 0, 0], [-0.5, s3/2, 0], [-0.5, -s3/2, 0]])
print('triangle', summary(tri, np.ones(3))[0])
for q4 in (0.01, 0.1, 0.5, 1.0):
    a = np.vstack([tri, [0, 0, 0]])
    print('triangle+centre q=%g' % q4, summary(a, np.array([1, 1, 1, q4]))[0])
for eps in (1/6, 0.1):
    qe = 0.75*eps**3 - 5/32*eps**5
    a = np.vstack([tri, [0, 0, eps], [0, 0, -eps]])
    print('ABK eps=%g' % eps, summary(a, np.array([1, 1, 1, qe, qe]), nseed=20000, rounds=4)[0])
