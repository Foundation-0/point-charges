import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np
from charges import count_fast
def cfg(k, eps, qq, h=0.0):
    th = 2 * np.pi * np.arange(k) / k
    a = np.c_[np.cos(th), np.sin(th), np.zeros(k)]
    a = np.vstack([a, [0, 0, eps + h], [0, 0, -eps + h]])
    return a, np.r_[np.ones(k), qq, qq]
for k in (3, 4, 5, 6):
    # H2 coefficient of the k-gon at its centre: V ~ k + (k/4)(r^2 - 2 z^2) + ... ; pair gives -(2q/eps^3)(r^2-2z^2)/2
    alpha = k / 4.0
    best = (0, None)
    for eps in (0.25, 0.15, 0.1, 0.06):
        for f in np.linspace(0.90, 1.02, 61):
            qq = f * alpha * eps ** 3
            n, m1, m2, ok = count_fast(*cfg(k, eps, qq))
            if ok and n > best[0]:
                best = (n, (eps, f, m1, m2))
    print('k=%d (n=%d charges): best %s' % (k, k + 2, best), flush=True)
