#!/usr/bin/env python3
"""Numerical sanity checks (not proofs) of the derivative-norm lemmas used by indep_count.py.

For f = 1/|x| at |x| = 1 (the bounds scale as r^-(k+1)):
  ||D^2 f||_2                 = 2
  max_{|v|=1} |D^3 f[v,v,.]|  = 6
  max_{|v|=1} ||D^3 f[v,.,.]||_2 = 6
The proofs are written in the comments of indep_count.py. This script samples random unit vectors and also checks
the closed-form tensors against finite differences.
"""
import numpy as np

rng = np.random.default_rng(1)


def T3(n):
    I = np.eye(3)
    return (3 * (np.einsum("ab,c->abc", I, n) + np.einsum("ac,b->abc", I, n) + np.einsum("bc,a->abc", I, n))
            - 15 * np.einsum("a,b,c->abc", n, n, n))


def H2(n):
    return 3 * np.outer(n, n) - np.eye(3)


def f(x):
    return 1 / np.linalg.norm(x)


# finite-difference check of the tensors at a random point
x = rng.normal(size=3)
r = np.linalg.norm(x)
n = x / r
h = 1e-4
H_fd = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        ei = np.eye(3)[i] * h
        ej = np.eye(3)[j] * h
        H_fd[i, j] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * h * h)
print("Hessian formula vs finite differences, max abs err:", np.abs(H_fd - H2(n) / r**3).max())


def grad(x):
    return -x / np.linalg.norm(x) ** 3


def hess(x):
    r = np.linalg.norm(x)
    return H2(x / r) / r**3


T_fd = np.zeros((3, 3, 3))
for k in range(3):
    ek = np.eye(3)[k] * h
    T_fd[:, :, k] = (hess(x + ek) - hess(x - ek)) / (2 * h)
print("third-derivative formula vs finite differences, max abs err:", np.abs(T_fd - T3(n) / r**4).max())

n = np.array([0, 0, 1.0])
T = T3(n)
print("||D^2|| at r=1:", np.abs(np.linalg.eigvalsh(H2(n))).max())
m1 = m2 = 0
for _ in range(200000):
    v = rng.normal(size=3)
    v /= np.linalg.norm(v)
    y = np.einsum("abc,a,b->c", T, v, v)
    m1 = max(m1, np.linalg.norm(y))
    M = np.einsum("abc,a->bc", T, v)
    m2 = max(m2, np.abs(np.linalg.eigvalsh(M)).max())
print("sampled max |D^3[v,v,.]| at r=1:", m1, "(bound 6)")
print("sampled max ||D^3[v,.,.]|| at r=1:", m2, "(bound 6)")
print("at v=n:", np.linalg.norm(np.einsum("abc,a,b->c", T, n, n)), np.abs(np.linalg.eigvalsh(np.einsum("abc,a->bc", T, n))).max())
