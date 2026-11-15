"""Check of the inversion rule: the degree-k Taylor term of V = sum q_i/|x - a_i| at the origin equals the
degree-k multipole polynomial of the inverted charges m_i = q_i/r_i at b_i = a_i/r_i^2.
Both sides are evaluated at random points for random charges; the Taylor term is taken from a high-order
finite-difference-free formula (Legendre expansion) and compared with direct numerical differentiation of V."""
import numpy as np
from numpy.polynomial import legendre as L
import mpmath as mp
rng = np.random.default_rng(1); n = 7
A = rng.normal(size=(n, 3)) + 2.0; q = rng.uniform(0.2, 2, size=n)
r = np.linalg.norm(A, axis=1); B = A / r[:, None] ** 2; m = q / r
x = rng.normal(size=3); x /= np.linalg.norm(x)
def V_line(s): return sum(mp.mpf(q[i]) / mp.sqrt(sum((s * mp.mpf(x[j]) - mp.mpf(A[i, j])) ** 2 for j in range(3))) for i in range(n))
mp.mp.dps = 40
print(" k   Taylor coefficient of V along x (numerical derivative)   multipole of the inverted charges")
for k in range(0, 7):
    taylor = mp.diff(V_line, 0, k) / mp.factorial(k)
    multipole = sum(m[i] * np.linalg.norm(B[i]) ** k * L.legval(x @ B[i] / np.linalg.norm(B[i]), [0] * k + [1]) for i in range(n))
    print(" %d   %+.12f   %+.12f" % (k, float(taylor), multipole))
