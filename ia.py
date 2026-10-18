"""Vectorised interval arithmetic on float64 with outward rounding.

Every elementary operation (+ - * / sqrt) is IEEE-754 correctly rounded to nearest, so
moving each computed endpoint one ulp outward (numpy.nextafter) gives a rigorous
enclosure. Intervals are pairs (lo, hi) of numpy arrays of the same shape.
NaN endpoints mean "no information" and must be treated as failure by callers.
"""
import numpy as np

_inf = np.inf


def dn(x):
    return np.nextafter(x, -_inf)


def up(x):
    return np.nextafter(x, _inf)


def const(c, shape=()):
    c = np.broadcast_to(np.asarray(c, float), shape) if shape else np.asarray(c, float)
    return (c, c)


def add(a, b):
    return (dn(a[0] + b[0]), up(a[1] + b[1]))


def sub(a, b):
    return (dn(a[0] - b[1]), up(a[1] - b[0]))


def neg(a):
    return (-a[1], -a[0])


def mul(a, b):
    p = np.stack([a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]])
    return (dn(np.min(p, 0)), up(np.max(p, 0)))


def scale(c, a):
    """Exact float c times interval a."""
    c = np.asarray(c, float)
    lo, hi = c * a[0], c * a[1]
    return (dn(np.minimum(lo, hi)), up(np.maximum(lo, hi)))


def sqr(a):
    lo2, hi2 = a[0] * a[0], a[1] * a[1]
    lo = np.where(a[0] >= 0, lo2, np.where(a[1] <= 0, hi2, 0.0))
    hi = np.maximum(lo2, hi2)
    return (np.where(lo > 0, dn(lo), 0.0), up(hi))


def sqrt(a):
    lo = np.maximum(a[0], 0.0)
    return (np.maximum(dn(np.sqrt(lo)), 0.0), up(np.sqrt(a[1])))


def recip(a):
    """1/a; returns NaN interval where a contains 0."""
    bad = (a[0] <= 0) & (a[1] >= 0)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        lo, hi = dn(1.0 / a[1]), up(1.0 / a[0])
    return (np.where(bad, np.nan, lo), np.where(bad, np.nan, hi))


def div(a, b):
    return mul(a, recip(b))


def contains_zero(a):
    """True where the interval may contain 0 (NaN counts as 'may')."""
    return ~((a[0] > 0) | (a[1] < 0))


def mag(a):
    return np.maximum(np.abs(a[0]), np.abs(a[1]))


def hull(a, b):
    return (np.minimum(a[0], b[0]), np.maximum(a[1], b[1]))
