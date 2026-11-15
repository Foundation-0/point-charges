"""Closed form for every four-charge configuration that is flat at the origin
(gradient and Hessian of the potential both vanish there), and the cubic term there.

Derivation. Write a_i = r_i u_i. Hessian zero  <=>  sum v_i u_i u_i^T = lambda I with v_i = q_i / r_i^3,
so the 3x4 matrix with columns sqrt(v_i) u_i has orthogonal rows; let n (unit, in R^4) span the missing direction.
Then v_i = lambda (1 - n_i^2) and u_i.u_j = -n_i n_j / sqrt((1-n_i^2)(1-n_j^2)).  Gradient zero <=> sum v_i r_i u_i = 0
<=> r_i proportional to t_i = n_i / sqrt(1 - n_i^2), which needs every n_i > 0.
So: flat configurations <-> points n of the positive part of the unit 3-sphere, with
     a_i.a_j = -t_i^2 t_j^2 (i != j),  |a_i| = t_i,  q_i = t_i^3 / (1 + t_i^2).
"""
import numpy as np
def flat_config(n):
    n = np.asarray(n, float); n = n / np.linalg.norm(n)
    t = n / np.sqrt(1 - n**2)
    G = -np.outer(t**2, t**2); np.fill_diagonal(G, t**2)
    w, V = np.linalg.eigh(G)
    A = V[:, 1:] * np.sqrt(np.maximum(w[1:], 0))
    q = t**3 / (1 + t**2)
    return A, q, w[0]
def field(X, A, Q):
    Dm = X[:, None, :] - A[None, :, :]; r2 = (Dm * Dm).sum(-1); r = np.sqrt(r2)
    g = -(Q[None, :, None] * Dm / (r2 * r)[:, :, None]).sum(1)
    H = (Q[None, :, None, None] * (3 * Dm[:, :, :, None] * Dm[:, :, None, :] - r2[:, :, None, None] * np.eye(3)) / (r2 * r2 * r)[:, :, None, None]).sum(1)
    return g, H
def cubic_parts(A, q):
    r = np.linalg.norm(A, axis=1); U = A / r[:, None]; c = q / r**4
    return U, c
if __name__ == "__main__":
    rng = np.random.default_rng(0)
    for trial in range(5):
        n = np.abs(rng.normal(size=4)) if trial else np.ones(4)
        A, q, w0 = flat_config(n)
        g, H = field(np.zeros((1, 3)), A, q)
        print("n =", np.round(n / np.linalg.norm(n), 4), " smallest Gram eigenvalue %.1e" % w0, " |grad| %.1e  |Hess| %.1e" % (np.abs(g).max(), np.abs(H).max()), " charges", np.round(q / q.max(), 4))
