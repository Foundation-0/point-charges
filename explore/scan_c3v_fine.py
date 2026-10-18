import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np
from c3v_slice import slice_cps
from collections import Counter
res = []
for h in np.r_[np.linspace(0.005, 0.1, 20), np.linspace(0.12, 1.2, 28)]:
    for lr in np.linspace(-1.5, 1.5, 61):       # q = 1.5 h^3 * 10^lr  (around the H2-cancelling value)
        q = 1.5 * h**3 * 10**lr
        c, pts = slice_cps(h, q, nseed=1500)
        res.append((h, lr, c))
print(Counter(r[2] for r in res))
best = max(r[2] for r in res)
for r in res:
    if r[2] >= 9: print(r)
