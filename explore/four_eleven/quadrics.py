"""All 8 solutions of three quadrics in three unknowns, batched, by the null space of the degree-4 Macaulay matrix
(eigenvectors of a multiplication matrix). Used for  grad C(x) + Q x + l = 0."""
import numpy as np, itertools
MON4 = [e for d in range(5) for e in itertools.product(range(d + 1), repeat=3) if sum(e) == d]
IDX4 = {e: i for i, e in enumerate(MON4)}
MON2 = [e for e in MON4 if sum(e) <= 2]
ADD = np.array([[IDX4[tuple(np.add(m, e))] for e in MON2] for m in MON2])       # [multiplier, monomial] -> column
BASIS = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)]
ROWS_B = np.array([IDX4[b] for b in BASIS])
ROWS_XB = [np.array([IDX4[tuple(np.add(b, np.eye(3, dtype=int)[a]))] for b in BASIS]) for a in range(3)]
def coefficients(T, Qm, l):
    """f_a = (1/2) T_abc x_b x_c + Q_ab x_b + l_a as coefficients on MON2. T [S,3,3,3], Qm [S,3,3], l [S,3]."""
    S = T.shape[0]; co = np.zeros((S, 3, len(MON2)))
    for j, e in enumerate(MON2):
        d = sum(e)
        if d == 0: co[:, :, j] = l
        elif d == 1: co[:, :, j] = Qm[:, :, e.index(1)]
        else:
            idx = [a for a in range(3) for _ in range(e[a])]
            co[:, :, j] = T[:, :, idx[0], idx[1]] * (0.5 if idx[0] == idx[1] else 1.0)
    return co
LMIX = np.array([[1.0, 0.37, -0.52], [0.29, 1.0, 0.61], [-0.43, 0.23, 1.0]])
def solve_tql(T, Qm, l, rng):
    """Solve in mixed coordinates x = LMIX y. For a harmonic cubic the x^2, y^2, z^2 coefficients of the three
    quadrics are linearly dependent (each row of them sums to zero), so the standard monomial basis fails in the
    original coordinates. A fixed non-orthogonal change of variables removes that."""
    Tm = np.einsum('sabc,bd,ce->sade', T, LMIX, LMIX); Qmm = np.einsum('sab,bd->sad', Qm, LMIX)
    Y, good, gap = solve(coefficients(Tm, Qmm, l), rng)
    return np.einsum('bd,skd->skb', LMIX, Y), good, gap
def solve(co, rng):
    S = co.shape[0]; M = np.zeros((S, 30, 35))
    for a in range(3):
        for m in range(10):
            M[:, a * 10 + m, ADD[m]] = co[:, a, :]
    u, sv, vt = np.linalg.svd(M)
    N = vt[:, 27:, :].transpose(0, 2, 1)                       # [S,35,8] null space
    gap = sv[:, 26] / sv[:, 0]
    NB = N[:, ROWS_B, :]
    r = rng.normal(size=3)
    NX = sum(r[a] * N[:, ROWS_XB[a], :] for a in range(3))
    cond = np.linalg.cond(NB)
    good = np.isfinite(cond) & (cond < 1e10)
    X = np.full((S, 8, 3), np.nan, complex)
    if good.any():
        Mx = NX[good] @ np.linalg.inv(NB[good])
        w, V = np.linalg.eig(Mx)                               # columns of V: values of BASIS at each root
        V = V / V[:, 0:1, :]
        X[good] = V[:, 1:4, :].transpose(0, 2, 1)
    return X, good, gap
def polish(X, T, Qm, l, iters=8):
    for it in range(iters):
        F = 0.5 * np.einsum('sabc,skb,skc->ska', T, X, X) + np.einsum('sab,skb->ska', Qm, X) + l[:, None, :]
        J = np.einsum('sabc,skc->skab', T, X) + Qm[:, None, :, :]
        with np.errstate(all='ignore'):
            try: X = X - np.linalg.solve(J, F[..., None])[..., 0]
            except np.linalg.LinAlgError: pass
    F = 0.5 * np.einsum('sabc,skb,skc->ska', T, X, X) + np.einsum('sab,skb->ska', Qm, X) + l[:, None, :]
    return X, np.abs(F).max(-1)
def tensor(U, c):
    """Third-derivative tensor of the cubic term. U [S,4,3], c [S,4]."""
    I = np.eye(3); cu = np.einsum('si,sia->sa', c, U)
    T = 15 * np.einsum('si,sia,sib,sic->sabc', c, U, U, U)
    T -= 3 * (I[None, :, :, None] * cu[:, None, None, :] + I[None, :, None, :] * cu[:, None, :, None] + I[None, None, :, :] * cu[:, :, None, None])
    return T
def real_count(X, res, tol=1e-9):
    """Number of real solutions per sample and whether the set is complete (8 distinct, small residual)."""
    S = X.shape[0]; scale = np.abs(X).max(axis=(1, 2)) + 1e-300
    ok = (res < 1e-8 * (1 + scale[:, None]**2)).all(1) & np.isfinite(X).all(axis=(1, 2))
    d = np.abs(X[:, :, None, :] - X[:, None, :, :]).max(-1) + np.eye(8)[None] * 1e300
    ok &= d.min(axis=(1, 2)) > 1e-7 * scale
    nreal = (np.abs(X.imag).max(-1) < tol * scale[:, None]).sum(1)
    return nreal, ok
if __name__ == "__main__":
    rng = np.random.default_rng(0)
    # self-test: x^2 = 1, y^2 = 4, z^2 = 9 written as T, Q, l
    T = np.zeros((1, 3, 3, 3)); 
    for a in range(3): T[0, a, a, a] = 2.0
    Qm = np.zeros((1, 3, 3)); l = -np.array([[1.0, 4.0, 9.0]])
    X, good, gap = solve_tql(T, Qm, l, rng); X, res = polish(X, T, Qm, l)
    print("self-test roots:", np.round(np.sort_complex(X[0][:, 0] * 100 + X[0][:, 1] * 10 + X[0][:, 2]), 6), real_count(X, res))
