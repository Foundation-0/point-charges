import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np, sys
from multiprocessing import Pool
from charges import summary
s3 = np.sqrt(3)
tri = np.array([[1, 0, 0], [-0.5, s3/2, 0], [-0.5, -s3/2, 0]])
def job(p):
    h, lq = p
    a = np.vstack([tri, [0, 0, h]]); q = np.array([1, 1, 1, 10**lq])
    s = summary(a, q, nseed=1500, rounds=2)[0]
    return h, lq, s['n'], s['m1'], s['morse_ok']
if __name__ == '__main__':
    P = [(h, lq) for h in np.linspace(0.0, 1.5, 31) for lq in np.linspace(-6, 1, 29)]
    with Pool(4) as pool:
        res = pool.map(job, P)
    from collections import Counter
    print(Counter(r[2] for r in res))
    for r in res:
        if r[2] >= 9: print(r)
