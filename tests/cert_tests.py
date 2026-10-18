import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E702
import numpy as np, time
from certify import run
s3 = 0.8660254037844386
tri = np.array([[1, 0, 0], [-0.5, s3, 0], [-0.5, -s3, 0]])
t = time.time(); run(tri, np.ones(3), 'triangle'); print(time.time() - t)
t = time.time(); run(np.vstack([tri, [0, 0, 0]]), np.array([1, 1, 1, 2.0**-10]), 'triangle+centre'); print(time.time() - t)
eps = 1/8; qe = 0.75*eps**3 - 5/32*eps**5
t = time.time(); run(np.vstack([tri, [0, 0, eps], [0, 0, -eps]]), np.array([1, 1, 1, qe, qe]), 'ABK eps=1/8'); print(time.time() - t)
