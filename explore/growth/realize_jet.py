"""Move charges (positions and values) so that the potential's jet at the origin, degrees 1..d-1, equals a target
harmonic jet J exactly (to ~1e-40 relative). Jets are computed in arb (python-flint) from
  V_l(x) = sum_i q_i sum_k c_{lk} (x.a_i)^{l-2k} (x.x)^k / |a_i|^{2l-2k+1},   P_l(t) = sum_k c_{lk} t^{l-2k}.
Newton with a float Jacobian (finite differences of float jets) and arb residuals. Output: rational configuration."""
import numpy as np, json, math, itertools
from fractions import Fraction
from flint import arb, ctx, fmpq
from numpy.polynomial import legendre
from harm import monos

ctx.prec = 400


def legcoef(l):
    c = legendre.leg2poly(np.eye(l + 1)[l])  # power-basis coefficients (float); exact rationals below
    out = {}
    for k in range(l // 2 + 1):
        # P_l(t) coefficient of t^{l-2k}: (-1)^k (2l-2k)! / (2^l k! (l-k)! (l-2k)!)
        out[k] = Fraction((-1) ** k * math.factorial(2 * l - 2 * k), 2 ** l * math.factorial(k) * math.factorial(l - k) * math.factorial(l - 2 * k))
    return out


def mono_pow_dot(a, m):
    """(x.a)^m as dict monomial -> coefficient (works for arb or float a)."""
    out = {}
    for i in range(m + 1):
        for j in range(m + 1 - i):
            k = m - i - j
            c = math.factorial(m) // (math.factorial(i) * math.factorial(j) * math.factorial(k))
            out[(i, j, k)] = c * a[0] ** i * a[1] ** j * a[2] ** k
    return out


def r2pow(k):
    out = {}
    for i in range(k + 1):
        for j in range(k + 1 - i):
            kk = k - i - j
            c = math.factorial(k) // (math.factorial(i) * math.factorial(j) * math.factorial(kk))
            out[(2 * i, 2 * j, 2 * kk)] = c
    return out


def pmul(A, B):
    out = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            m = (ma[0] + mb[0], ma[1] + mb[1], ma[2] + mb[2])
            out[m] = out.get(m, 0) + ca * cb
    return out


def jet(sites, q, l, one):
    """degree-l part, dict monomial -> value; `one` is arb(1) or 1.0 to fix the arithmetic."""
    lc = legcoef(l)
    tot = {}
    for a, qi in zip(sites, q):
        rr = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
        r = rr.sqrt() if hasattr(rr, 'sqrt') else math.sqrt(rr)
        for k, ck in lc.items():
            coef = qi * (one * ck.numerator) / ck.denominator / r ** (2 * l - 2 * k + 1)
            term = pmul(mono_pow_dot(a, l - 2 * k), r2pow(k))
            for m, v in term.items():
                tot[m] = tot.get(m, 0) + coef * v
    return tot


def jetvec(sites, q, d, one):
    out = []
    for l in range(1, d):
        L, _ = monos(l)
        J = jet(sites, q, l, one)
        out += [J.get(m, 0 * one) for m in L if sum(m) == l]
    return out


def realize(sites0, q0, d, target, iters=30, verbose=True):
    """target: list (same layout as jetvec) of Fractions/floats. Unknowns: sites (3n) and log charges (n)."""
    n = len(q0)
    x = np.r_[np.asarray(sites0, float).ravel(), np.log(np.asarray(q0, float))]
    tgt_arb = [arb(fmpq(Fraction(t).numerator, Fraction(t).denominator)) if not isinstance(t, arb) else t for t in target]
    tgt_f = np.array([float(Fraction(t)) for t in target])

    def F_float(x):
        s = x[:3 * n].reshape(-1, 3); q = np.exp(x[3 * n:])
        return np.array([float(v) for v in jetvec(s, q, d, 1.0)]) - tgt_f

    def F_arb(xa):
        s = [xa[3 * i:3 * i + 3] for i in range(n)]; q = [xa[3 * n + i].exp() for i in range(n)]
        v = jetvec(s, q, d, arb(1))
        return [vi - ti for vi, ti in zip(v, tgt_arb)]

    xa = [arb(float(v)) for v in x]
    scale = np.maximum(np.abs(tgt_f), 1e-300)
    for it in range(iters * 3):
        r = F_arb(xa)
        rf = np.array([float(v.mid()) for v in r])
        xf = np.array([float(v.mid()) for v in xa])
        h = 1e-7
        J = np.zeros((len(rf), len(xf)))
        for k in range(len(xf)):
            e = np.zeros(len(xf)); e[k] = h
            J[:, k] = (F_float(xf + e) - F_float(xf - e)) / (2 * h)
        dx = np.linalg.lstsq(J, -rf, rcond=1e-12)[0]
        lam = 1.0
        if np.max(np.abs(rf)) > 1e-12 * np.max(np.abs(tgt_f)):   # damped phase (float is accurate enough here)
            for _ in range(30):
                if np.max(np.abs(F_float(xf + lam * dx))) < np.max(np.abs(rf)):
                    break
                lam *= 0.5
        xa = [(xa[k] + arb(float(lam * dx[k]))).mid() for k in range(len(xa))]
        rel = np.max(np.abs(rf)) / np.max(np.abs(tgt_f))
        if verbose:
            print(f"  realize it {it}: max residual {np.max(np.abs(rf)):.3e} (relative to max target {rel:.2e}); smallest target {np.min(np.abs(tgt_f[tgt_f!=0])):.2e}", flush=True)
        if np.max(np.abs(rf)) < 1e-10 * np.min(np.abs(tgt_f[tgt_f != 0])):
            break
    return xa
