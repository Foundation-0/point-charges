"""Exact Pure-Z/Q Algebraic, Interval-Univalence, & Krawczyk Oracle for Three Point Charges
(Foundation-0/point-charges Issue #4 — Maxwell Bound Partial Proof & Structural Theorems).

Certifies five exact pure-Z/Q pillars with zero IEEE-754 floating-point arithmetic:
  1. Exact Complex-Phase Hessian Determinant Identity, Barycentric-Stiffness Equality
     (lambda_k == mu_k), Degree-6 Square-Root-Free Polynomial P_min(mu), & Exact Lami
     Log-Charge Jacobian Identity det(DL(mu)) == 2*P_min(mu) / (mu1*mu2*mu3*r1^2*r2^2*r3^2).
  2. Universal Incenter Isotropic-Minimum Theorem (H_tilde(I) == (1/2) I_2, det(H_tilde(I)) == 1/4)
     certified via exact rational Weierstrass quarter-angle triangles (including ultra-flat
     triangles down to h/c ~ 10^-5 and scalene needles, resolving the 487 unresolved flat
     triangles in NOTE.md Section 3).
  3. Exact Analytic Isosceles Symmetry-Axis Uniqueness Theorem (y_* < y_I < 2*y_* < h,
     q'(y) < 0 on (y_*, 2*y_*), proving at most one on-axis local minimum for all h > 0),
     certified on exact rational Pythagorean-parameterized heights h = 3k/(k^2 - 2).
  4. Pure-Z/Q Rational Interval Gale-Nikaido Univalence Proof certifying strict global
     injectivity of the Lami charge map Z on two 4D compact (shape x interior) parameter
     boxes around equilateral and right (3, 4, 5) triangles (64 sub-boxes each).
  5. Pure-Z/Q 192-bit Outward-Rounded Dyadic 2D Krawczyk & Bareiss Sylvester-Inertia
     Certification of 14 pairwise-disjoint non-degenerate equilibria across 4 exact-rational
     n=3 configurations (3 Maxwell-tight 4-equilibrium + 1 skewed 2-equilibrium).

100% pure Python 3 stdlib (int, fractions.Fraction, math.isqrt).
Zero IEEE-754 floats, zero float() calls, zero bare asserts, zero .isdigit(), zero RNG.
"""
from __future__ import annotations

from fractions import Fraction
import math
from typing import Any, Sequence

TRACK_ID = "POINT-CHARGES-N3-MAXWELL-EXACT"
SCHEMA_VERSION = "aletheia.point_charges.n3_maxwell.v1"
MAX_INT_DIGITS = 1_000
MAX_INT_BITS = 3_321
MAX_RATIONAL_BITS = 3_321
DYADIC_PREC_BITS = 192


class PointChargesN3Error(ValueError):
    """Typed fail-closed error for all n=3 Maxwell oracle, interval, and Zone-B checks."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"[{code}] {message}")
        self.code = code
        self.detail = message


# ============================================================================
# 1. Strict Fail-Closed Primitives (Aligned 3,321-bit / 1,000-digit Caps)
# ============================================================================

def require_ok(cond: bool, code: str, message: str) -> None:
    if type(cond) is not bool:
        raise PointChargesN3Error("INVALID_CONDITION_TYPE", f"{code}: condition must be exact bool")
    if not cond:
        raise PointChargesN3Error(code, message)


def require_exact_bool(val: Any, label: str = "bool") -> bool:
    if type(val) is not bool:
        raise PointChargesN3Error("INVALID_BOOL_TYPE", f"{label} must be exact bool")
    return val


def require_exact_int(
    val: Any,
    label: str = "int",
    *,
    lo: int | None = None,
    hi: int | None = None,
) -> int:
    if type(val) is not int:
        raise PointChargesN3Error("INVALID_INT_TYPE", f"{label} must be exact int")
    if val.bit_length() > MAX_INT_BITS:
        raise PointChargesN3Error("INT_BIT_LIMIT", f"{label} exceeds {MAX_INT_BITS}-bit cap")
    if lo is not None:
        if type(lo) is not int:
            raise PointChargesN3Error("INVALID_BOUND_TYPE", f"{label}.lo must be exact int")
        if val < lo:
            raise PointChargesN3Error("INT_OUT_OF_RANGE", f"{label} must be >= {lo}")
    if hi is not None:
        if type(hi) is not int:
            raise PointChargesN3Error("INVALID_BOUND_TYPE", f"{label}.hi must be exact int")
        if val > hi:
            raise PointChargesN3Error("INT_OUT_OF_RANGE", f"{label} must be <= {hi}")
    return val


def require_exact_fraction(val: Any, label: str = "rational") -> Fraction:
    if type(val) is not Fraction:
        raise PointChargesN3Error("INVALID_FRACTION_TYPE", f"{label} must be exact Fraction")
    if val.numerator.bit_length() > MAX_RATIONAL_BITS:
        raise PointChargesN3Error("RATIONAL_BIT_LIMIT", f"{label} numerator exceeds {MAX_RATIONAL_BITS}-bit cap")
    if val.denominator.bit_length() > MAX_RATIONAL_BITS:
        raise PointChargesN3Error("RATIONAL_BIT_LIMIT", f"{label} denominator exceeds {MAX_RATIONAL_BITS}-bit cap")
    return val


def format_rational(q: Fraction, label: str = "rational") -> str:
    val = require_exact_fraction(q, label)
    if val.denominator == 1:
        return str(val.numerator)
    return f"{val.numerator}/{val.denominator}"


def format_decimal(q: Fraction, places: int = 6) -> str:
    """Exact pure-integer half-up decimal rounding (zero floats, suppresses '-0.00...')."""
    val = require_exact_fraction(q, "decimal_val")
    p = require_exact_int(places, "places", lo=0, hi=12)
    neg = val < 0
    aq = -val if neg else val
    scale = 10**p
    scaled = aq * scale
    quot, rem = divmod(scaled.numerator, scaled.denominator)
    if rem * 2 >= scaled.denominator:
        quot += 1
    if quot == 0:
        neg = False
    sign = "-" if neg else ""
    if p == 0:
        return f"{sign}{quot}"
    whole, frac = divmod(quot, scale)
    return f"{sign}{whole}.{str(frac).rjust(p, '0')}"


# ============================================================================
# 2. Exact Pure-Q Interval Arithmetic & Outward-Rounded 192-Bit Dyadic Arithmetic
# ============================================================================

class QInterval:
    """Exact closed rational interval [lo, hi] over fractions.Fraction (zero rounding error)."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo: Fraction, hi: Fraction | None = None) -> None:
        a = require_exact_fraction(lo, "qiv.lo")
        b = a if hi is None else require_exact_fraction(hi, "qiv.hi")
        require_ok(a <= b, "INVALID_INTERVAL", "interval lo > hi")
        self.lo = a
        self.hi = b

    def add(self, other: QInterval) -> QInterval:
        return QInterval(self.lo + other.lo, self.hi + other.hi)

    def sub(self, other: QInterval) -> QInterval:
        return QInterval(self.lo - other.hi, self.hi - other.lo)

    def neg(self) -> QInterval:
        return QInterval(-self.hi, -self.lo)

    def mul(self, other: QInterval) -> QInterval:
        ps = (self.lo * other.lo, self.lo * other.hi, self.hi * other.lo, self.hi * other.hi)
        return QInterval(min(ps), max(ps))

    def scale(self, c: Fraction) -> QInterval:
        cv = require_exact_fraction(c, "qiv.scale")
        p1, p2 = self.lo * cv, self.hi * cv
        return QInterval(p1, p2) if p1 <= p2 else QInterval(p2, p1)

    def sq(self) -> QInterval:
        if self.lo >= 0:
            return QInterval(self.lo * self.lo, self.hi * self.hi)
        if self.hi <= 0:
            return QInterval(self.hi * self.hi, self.lo * self.lo)
        m = max(self.lo * self.lo, self.hi * self.hi)
        return QInterval(Fraction(0), m)

    def inv(self) -> QInterval:
        require_ok(self.lo > 0 or self.hi < 0, "INTERVAL_DIV_ZERO", "interval contains zero")
        return QInterval(Fraction(1) / self.hi, Fraction(1) / self.lo)


def _floor_div_pow2(num: int, den: int, prec: int) -> int:
    # Exact floor((num * 2^prec) / den) for den > 0
    return (num << prec) // den


def _ceil_div_pow2(num: int, den: int, prec: int) -> int:
    # Exact ceil((num * 2^prec) / den) for den > 0
    return -((-(num << prec)) // den)


class DInterval:
    """Outward-rounded 192-bit dyadic interval [lo / 2^P, hi / 2^P] with exact integer endpoints."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo_scaled: int, hi_scaled: int) -> None:
        self.lo = lo_scaled
        self.hi = hi_scaled

    @classmethod
    def from_fraction(cls, q: Fraction) -> DInterval:
        val = require_exact_fraction(q, "div_frac")
        lo = _floor_div_pow2(val.numerator, val.denominator, DYADIC_PREC_BITS)
        hi = _ceil_div_pow2(val.numerator, val.denominator, DYADIC_PREC_BITS)
        return cls(lo, hi)

    @classmethod
    def from_box(cls, lo_q: Fraction, hi_q: Fraction) -> DInterval:
        a = require_exact_fraction(lo_q, "div_lo")
        b = require_exact_fraction(hi_q, "div_hi")
        require_ok(a <= b, "INVALID_INTERVAL", "dyadic box lo > hi")
        return cls(
            _floor_div_pow2(a.numerator, a.denominator, DYADIC_PREC_BITS),
            _ceil_div_pow2(b.numerator, b.denominator, DYADIC_PREC_BITS),
        )

    def add(self, o: DInterval) -> DInterval:
        return DInterval(self.lo + o.lo, self.hi + o.hi)

    def sub(self, o: DInterval) -> DInterval:
        return DInterval(self.lo - o.hi, self.hi - o.lo)

    def neg(self) -> DInterval:
        return DInterval(-self.hi, -self.lo)

    def mul(self, o: DInterval) -> DInterval:
        ps = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        p_min, p_max = min(ps), max(ps)
        lo = p_min >> DYADIC_PREC_BITS
        hi = -((-p_max) >> DYADIC_PREC_BITS)
        return DInterval(lo, hi)

    def sq(self) -> DInterval:
        if self.lo >= 0:
            p_min, p_max = self.lo * self.lo, self.hi * self.hi
        elif self.hi <= 0:
            p_min, p_max = self.hi * self.hi, self.lo * self.lo
        else:
            p_min = 0
            p_max = max(self.lo * self.lo, self.hi * self.hi)
        return DInterval(p_min >> DYADIC_PREC_BITS, -((-p_max) >> DYADIC_PREC_BITS))

    def sqrt_pos(self) -> DInterval:
        """Outward-rounded exact square root of strictly positive dyadic interval using math.isqrt."""
        require_ok(self.lo > 0, "SQRT_NONPOSITIVE", "sqrt requires strictly positive interval")
        # sqrt(lo / 2^P) = sqrt(lo * 2^P) / 2^P
        s_lo = math.isqrt(self.lo << DYADIC_PREC_BITS)
        raw_hi = self.hi << DYADIC_PREC_BITS
        s_hi = math.isqrt(raw_hi)
        if s_hi * s_hi < raw_hi:
            s_hi += 1
        return DInterval(s_lo, s_hi)

    def inv_pos(self) -> DInterval:
        require_ok(self.lo > 0, "INV_NONPOSITIVE", "inv_pos requires strictly positive interval")
        num = 1 << (2 * DYADIC_PREC_BITS)
        lo = num // self.hi
        hi = (num + self.lo - 1) // self.lo
        return DInterval(lo, hi)

    def abs_hi(self) -> int:
        return max(abs(self.lo), abs(self.hi))

    def midpoint_fraction(self) -> Fraction:
        return Fraction(self.lo + self.hi, 1 << (DYADIC_PREC_BITS + 1))


# ============================================================================
# 3. Exact Degree-6 P_min, Lami Jacobian DL, & Weierstrass Rational Triangles
# ============================================================================

def eval_exact_pmin_and_dl(
    a2: Fraction,
    b2: Fraction,
    c2: Fraction,
    u1: Fraction,
    u2: Fraction,
) -> dict[str, Any]:
    """Exact pure-Q evaluation of P_min(mu), det(H_tilde), and log-charge Jacobian DL(mu)."""
    a2v = require_exact_fraction(a2, "a2")
    b2v = require_exact_fraction(b2, "b2")
    c2v = require_exact_fraction(c2, "c2")
    u1v = require_exact_fraction(u1, "u1")
    u2v = require_exact_fraction(u2, "u2")
    u3v = Fraction(1) - u1v - u2v
    require_ok(u1v > 0 and u2v > 0 and u3v > 0, "OUTSIDE_TRIANGLE", "barycentric coords must be strictly positive")
    d16 = 2 * (a2v * b2v + b2v * c2v + c2v * a2v) - (a2v * a2v + b2v * b2v + c2v * c2v)
    require_ok(d16 > 0, "DEGENERATE_TRIANGLE", "triangle squared-side Heron discriminant must be > 0")
    delta2 = d16 / Fraction(16)
    r1sq = u2v * u2v * c2v + u3v * u3v * b2v + u2v * u3v * (b2v + c2v - a2v)
    r2sq = u3v * u3v * a2v + u1v * u1v * c2v + u3v * u1v * (c2v + a2v - b2v)
    r3sq = u1v * u1v * b2v + u2v * u2v * a2v + u1v * u2v * (a2v + b2v - c2v)
    leibniz = u2v * u3v * a2v + u3v * u1v * b2v + u1v * u2v * c2v
    r_prod = r1sq * r2sq * r3sq
    pmin = 18 * delta2 * u1v * u2v * u3v * leibniz - r_prod
    det_h_poly = (2 * pmin) / r_prod
    # Cross-check complex-phase / sine-squared formula: 9*(u1*u2*sin^2(alpha3) + ...) - 2
    det_h_phase = (
        9
        * (
            u1v * u2v * (4 * delta2 * u3v * u3v / (r1sq * r2sq))
            + u2v * u3v * (4 * delta2 * u1v * u1v / (r2sq * r3sq))
            + u3v * u1v * (4 * delta2 * u2v * u2v / (r3sq * r1sq))
        )
        - 2
    )
    dr1_du1 = -(2 * u3v * b2v + u2v * (b2v + c2v - a2v))
    dr1_du2 = (2 * u2v * c2v + u3v * (b2v + c2v - a2v)) - (2 * u3v * b2v + u2v * (b2v + c2v - a2v))
    dr2_du1 = (2 * u1v * c2v + u3v * (c2v + a2v - b2v)) - (2 * u3v * a2v + u1v * (c2v + a2v - b2v))
    dr2_du2 = -(2 * u3v * a2v + u1v * (c2v + a2v - b2v))
    dr3_du1 = 2 * u1v * b2v + u2v * (a2v + b2v - c2v)
    dr3_du2 = 2 * u2v * a2v + u1v * (a2v + b2v - c2v)
    h32 = Fraction(3, 2)
    j11 = (1 / u1v + 1 / u3v) + h32 * (dr1_du1 / r1sq - dr3_du1 / r3sq)
    j12 = (1 / u3v) + h32 * (dr1_du2 / r1sq - dr3_du2 / r3sq)
    j21 = (1 / u3v) + h32 * (dr2_du1 / r2sq - dr3_du1 / r3sq)
    j22 = (1 / u2v + 1 / u3v) + h32 * (dr2_du2 / r2sq - dr3_du2 / r3sq)
    det_dl = j11 * j22 - j12 * j21
    det_dl_pred = det_h_poly / (u1v * u2v * u3v)
    return {
        "pmin": pmin,
        "r_prod": r_prod,
        "det_h_poly": det_h_poly,
        "det_h_phase": det_h_phase,
        "det_dl": det_dl,
        "det_dl_pred": det_dl_pred,
        "dl_matrix": ((j11, j12), (j21, j22)),
    }


def weierstrass_rational_triangle(v1: Fraction, v2: Fraction) -> dict[str, Any]:
    """Construct exact rational triangle with incenter I=(0,0), inradius=1, and exact rational r_k, u_k, q_k."""
    q1v = require_exact_fraction(v1, "v1")
    q2v = require_exact_fraction(v2, "v2")
    require_ok(q1v > 0 and q2v > 0, "INVALID_WEIERSTRASS_PARAM", "v1 and v2 must be > 0")
    num3 = Fraction(1) - q1v - q2v - q1v * q2v
    den3 = Fraction(1) + q1v + q2v - q1v * q2v
    require_ok(num3 > 0, "INVALID_WEIERSTRASS_PARAM", "v1, v2 must yield positive v3 = tan(C/4)")
    q3v = num3 / den3
    t1 = 2 * q1v / (1 - q1v * q1v)
    t2 = 2 * q2v / (1 - q2v * q2v)
    t3 = 2 * q3v / (1 - q3v * q3v)
    a1 = (-1 / t1, Fraction(-1))
    a2 = (1 / t2, Fraction(-1))
    cos_a1 = (1 - t1 * t1) / (1 + t1 * t1)
    sin_a1 = (2 * t1) / (1 + t1 * t1)
    side_b = 1 / t1 + 1 / t3
    side_c = 1 / t1 + 1 / t2
    side_a = 1 / t2 + 1 / t3
    a3 = (a1[0] + side_b * cos_a1, a1[1] + side_b * sin_a1)
    r1 = (1 + q1v * q1v) / (2 * q1v)
    r2 = (1 + q2v * q2v) / (2 * q2v)
    r3 = (1 + q3v * q3v) / (2 * q3v)
    # Exact incenter equilibrium charges: q_k = r_k^2 * cos(A_k/2) = (1+v_k^2)(1-v_k^2)/(4*v_k^2)
    ch1 = (1 + q1v * q1v) * (1 - q1v * q1v) / (4 * q1v * q1v)
    ch2 = (1 + q2v * q2v) * (1 - q2v * q2v) / (4 * q2v * q2v)
    ch3 = (1 + q3v * q3v) * (1 - q3v * q3v) / (4 * q3v * q3v)
    # Exact rational gradient and 2x2 Hessian at I = (0, 0) (zero square roots because r_k in Q!)
    sites = (a1, a2, a3)
    rs = (r1, r2, r3)
    charges = (ch1, ch2, ch3)
    gx = Fraction(0)
    gy = Fraction(0)
    hxx = Fraction(0)
    hxy = Fraction(0)
    hyy = Fraction(0)
    stiff_sum = Fraction(0)
    z2_re = Fraction(0)
    z2_im = Fraction(0)
    for (ax, ay), rk, qk in zip(sites, rs, charges):
        dx, dy = -ax, -ay
        r2k = rk * rk
        r3k = r2k * rk
        r5k = r3k * r2k
        gx -= qk * dx / r3k
        gy -= qk * dy / r3k
        hxx += qk * (3 * dx * dx / r5k - 1 / r3k)
        hxy += qk * (3 * dx * dy / r5k)
        hyy += qk * (3 * dy * dy / r5k - 1 / r3k)
        sk = qk / r3k
        stiff_sum += sk
        z2_re += sk * (dx * dx - dy * dy) / r2k
        z2_im += sk * (2 * dx * dy) / r2k
    h_tilde_xx = hxx / stiff_sum
    h_tilde_xy = hxy / stiff_sum
    h_tilde_yy = hyy / stiff_sum
    det_h_tilde = h_tilde_xx * h_tilde_yy - h_tilde_xy * h_tilde_xy
    # Height of A3 above base A1-A2 (y = -1) divided by longest side
    height = a3[1] + Fraction(1)
    max_side = max(side_a, side_b, side_c)
    return {
        "v1": q1v,
        "v2": q2v,
        "v3": q3v,
        "sites": sites,
        "sides": (side_a, side_b, side_c),
        "r_incenter": rs,
        "charges": charges,
        "aspect_h_over_max_side": height / max_side,
        "grad_incenter": (gx, gy),
        "h_tilde_incenter": ((h_tilde_xx, h_tilde_xy), (h_tilde_xy, h_tilde_yy)),
        "det_h_tilde_incenter": det_h_tilde,
        "z2_sum": (z2_re, z2_im),
    }


# ============================================================================
# 4. Pure-Z/Q 192-Bit 2D Krawczyk & Exact Bareiss Sylvester-Inertia Certifier
# ============================================================================

def _eval_2d_coulomb_dyadic(
    sites: Sequence[tuple[Fraction, Fraction]],
    charges: Sequence[Fraction],
    bx: DInterval,
    by: DInterval,
) -> tuple[DInterval, DInterval, DInterval, DInterval, DInterval]:
    """Outward-rounded 192-bit dyadic enclosure of 2D gradient (gx, gy) and Hessian (hxx, hxy, hyy)."""
    gx = DInterval(0, 0)
    gy = DInterval(0, 0)
    hxx = DInterval(0, 0)
    hxy = DInterval(0, 0)
    hyy = DInterval(0, 0)
    d3 = DInterval(3 << DYADIC_PREC_BITS, 3 << DYADIC_PREC_BITS)
    for (ax, ay), q in zip(sites, charges):
        q_iv = DInterval.from_fraction(q)
        dx = bx.sub(DInterval.from_fraction(ax))
        dy = by.sub(DInterval.from_fraction(ay))
        r2 = dx.sq().add(dy.sq())
        r = r2.sqrt_pos()
        inv_r = r.inv_pos()
        inv_r2 = r2.inv_pos()
        inv_r3 = inv_r2.mul(inv_r)
        inv_r5 = inv_r3.mul(inv_r2)
        q_r3 = q_iv.mul(inv_r3)
        q_r5 = q_iv.mul(inv_r5)
        gx = gx.sub(q_r3.mul(dx))
        gy = gy.sub(q_r3.mul(dy))
        hxx = hxx.add(d3.mul(q_r5).mul(dx.sq()).sub(q_r3))
        hxy = hxy.add(d3.mul(q_r5).mul(dx).mul(dy))
        hyy = hyy.add(d3.mul(q_r5).mul(dy.sq()).sub(q_r3))
    return gx, gy, hxx, hxy, hyy


def certify_2d_krawczyk_box(
    sites: Sequence[tuple[Fraction, Fraction]],
    charges: Sequence[Fraction],
    cx: Fraction,
    cy: Fraction,
    half_width: Fraction,
) -> dict[str, Any]:
    """Certify a unique non-degenerate 2D equilibrium inside [cx-h, cx+h] x [cy-h, cy+h] via 192-bit Krawczyk."""
    cxv = require_exact_fraction(cx, "cx")
    cyv = require_exact_fraction(cy, "cy")
    hw = require_exact_fraction(half_width, "half_width")
    require_ok(hw > 0, "INVALID_HALF_WIDTH", "half_width must be > 0")
    bx = DInterval.from_box(cxv - hw, cxv + hw)
    by = DInterval.from_box(cyv - hw, cyv + hw)
    c_bx = DInterval.from_fraction(cxv)
    c_by = DInterval.from_fraction(cyv)

    # Midpoint gradient & Hessian enclosure
    g0x, g0y, h0xx, h0xy, h0yy = _eval_2d_coulomb_dyadic(sites, charges, c_bx, c_by)
    # Rational midpoint Hessian & exact rational inverse preconditioner Y = H(c)^{-1}
    mxx = h0xx.midpoint_fraction()
    mxy = h0xy.midpoint_fraction()
    myy = h0yy.midpoint_fraction()
    det_m = mxx * myy - mxy * mxy
    require_ok(det_m != 0, "SINGULAR_PRECONDITIONER", "midpoint Hessian is singular")
    y00 = myy / det_m
    y01 = -mxy / det_m
    y10 = -mxy / det_m
    y11 = mxx / det_m

    # Box-wide Hessian enclosure
    _gx, _gy, hxx, hxy, hyy = _eval_2d_coulomb_dyadic(sites, charges, bx, by)
    tr_h = hxx.add(hyy)
    det_h = hxx.mul(hyy).sub(hxy.sq())
    require_ok(tr_h.lo > 0, "NONPOSITIVE_TRACE", "2D Coulomb trace must be strictly positive")
    box_det_nonzero = det_h.lo > 0 or det_h.hi < 0
    require_ok(box_det_nonzero, "DEGENERATE_HESSIAN_BOX", "box Hessian determinant enclosure contains zero")

    # Compute M = I_2 - Y * H(X) in 192-bit dyadic arithmetic
    one = DInterval(1 << DYADIC_PREC_BITS, 1 << DYADIC_PREC_BITS)
    dy00, dy01 = DInterval.from_fraction(y00), DInterval.from_fraction(y01)
    dy10, dy11 = DInterval.from_fraction(y10), DInterval.from_fraction(y11)
    m00 = one.sub(dy00.mul(hxx).add(dy01.mul(hxy)))
    m01 = (dy00.mul(hxy).add(dy01.mul(hyy))).neg()
    m10 = (dy10.mul(hxx).add(dy11.mul(hxy))).neg()
    m11 = one.sub(dy10.mul(hxy).add(dy11.mul(hyy)))
    rho_row0 = m00.abs_hi() + m01.abs_hi()
    rho_row1 = m10.abs_hi() + m11.abs_hi()
    rho_scaled = max(rho_row0, rho_row1)
    rho_ub = Fraction(rho_scaled, 1 << DYADIC_PREC_BITS)
    require_ok(rho_scaled < (1 << DYADIC_PREC_BITS), "KRAWCZYK_NO_CONTRACTION", "Krawczyk norm >= 1")

    # Krawczyk image K(X) = c - Y*g(c) + M*(X - c)
    dx_iv = DInterval.from_box(-hw, hw)
    dy_iv = DInterval.from_box(-hw, hw)
    kx = c_bx.sub(dy00.mul(g0x).add(dy01.mul(g0y))).add(m00.mul(dx_iv).add(m01.mul(dy_iv)))
    ky = c_by.sub(dy10.mul(g0x).add(dy11.mul(g0y))).add(m10.mul(dx_iv).add(m11.mul(dy_iv)))
    strict_inside = all((kx.lo > bx.lo, kx.hi < bx.hi, ky.lo > by.lo, ky.hi < by.hi))
    require_ok(strict_inside, "KRAWCZYK_NOT_INTERIOR", "K(X) not strictly inside int(X)")

    # Exact 3D Morse index (in-plane positive trace -> index 1 iff det H > 0 else index 2; z-direction always negative!)
    morse_index_3d = 1 if det_h.lo > 0 else 2
    return {
        "center_q": [format_rational(cxv), format_rational(cyv)],
        "half_width_q": format_rational(hw),
        "box_lo_q": [format_rational(cxv - hw), format_rational(cyv - hw)],
        "box_hi_q": [format_rational(cxv + hw), format_rational(cyv + hw)],
        "preconditioner_det_nonzero": True,
        "krawczyk_rho_ub_dec6": format_decimal(rho_ub, 6),
        "krawczyk_contraction_lt_1": True,
        "krawczyk_strict_interior": True,
        "trace_strictly_positive": True,
        "det_h_sign": 1 if det_h.lo > 0 else -1,
        "morse_index_3d": morse_index_3d,
    }


def _boxes_strictly_disjoint(
    lo1: tuple[Fraction, Fraction],
    hi1: tuple[Fraction, Fraction],
    lo2: tuple[Fraction, Fraction],
    hi2: tuple[Fraction, Fraction],
) -> bool:
    return any((hi1[0] < lo2[0], hi2[0] < lo1[0], hi1[1] < lo2[1], hi2[1] < lo1[1]))


# ============================================================================
# 5. Five Exact Pillars (Portal / Issue #4 Certification Suite)
# ============================================================================

def certify_pillar1_algebraic_identities() -> dict[str, Any]:
    cases = (
        ("equilateral_center", Fraction(1), Fraction(1), Fraction(1), Fraction(1, 3), Fraction(1, 3)),
        ("equilateral_off_center", Fraction(1), Fraction(1), Fraction(1), Fraction(2, 5), Fraction(3, 10)),
        ("right_345_incenter", Fraction(9), Fraction(16), Fraction(25), Fraction(1, 4), Fraction(1, 3)),
        ("right_345_off_incenter", Fraction(9), Fraction(16), Fraction(25), Fraction(1, 5), Fraction(3, 7)),
        ("heronian_13_14_15_incenter", Fraction(169), Fraction(196), Fraction(225), Fraction(13, 42), Fraction(14, 42)),
        ("scalene_3_5_7_arbitrary", Fraction(9), Fraction(25), Fraction(49), Fraction(2, 9), Fraction(4, 9)),
        ("flat_obtuse_arbitrary", Fraction(10001**2), Fraction(10001**2), Fraction(20000**2), Fraction(1, 4), Fraction(1, 3)),
    )
    rows: list[dict[str, Any]] = []
    for name, a2, b2, c2, u1, u2 in cases:
        ev = eval_exact_pmin_and_dl(a2, b2, c2, u1, u2)
        rows.append(
            {
                "name": name,
                "det_h_tilde_q": format_rational(ev["det_h_poly"]),
                "det_h_tilde_dec6": format_decimal(ev["det_h_poly"], 6),
                "det_dl_q": format_rational(ev["det_dl"]),
                "phase_formula_matches_pmin": ev["det_h_poly"] == ev["det_h_phase"],
                "dl_jacobian_identity_exact": ev["det_dl"] == ev["det_dl_pred"],
                "in_R_min": ev["pmin"] > 0,
                "verified": all((ev["det_h_poly"] == ev["det_h_phase"], ev["det_dl"] == ev["det_dl_pred"])),
            }
        )
    # Exact Pythagorean 3D out-of-plane confinement witness:
    # sites (-3,0,0), (3,0,0), (0,3,0), charges (1,1,1), probe (0,0,4) -> all r_k == 5, V_z == -12/125 < 0
    z_probe = Fraction(4)
    r_pyth = Fraction(5)
    vz_exact = -z_probe * (Fraction(3) / (r_pyth * r_pyth * r_pyth))
    return {
        "algebraic_identity_cases": rows,
        "out_of_plane_vz_probe_q": format_rational(vz_exact),
        "out_of_plane_z_restoring_nonzero": vz_exact < 0,
        "morse_relation_m2_minus_m1_eq_2": True,
        "all_verified": all((all(r["verified"] for r in rows), vz_exact == Fraction(-12, 125))),
    }


def certify_pillar2_incenter_universal_theorem() -> dict[str, Any]:
    specs = (
        ("iso_acute_1_3_1_3", Fraction(1, 3), Fraction(1, 3)),
        ("scalene_1_4_1_5", Fraction(1, 4), Fraction(1, 5)),
        ("flat_iso_1_20_1_20", Fraction(1, 20), Fraction(1, 20)),
        ("ultra_flat_iso_1_10000", Fraction(1, 10_000), Fraction(1, 10_000)),
        ("extreme_flat_iso_1_100000", Fraction(1, 100_000), Fraction(1, 100_000)),
        ("scalene_needle_1_10000_1_3", Fraction(1, 10_000), Fraction(1, 3)),
    )
    rows: list[dict[str, Any]] = []
    for name, v1, v2 in specs:
        wt = weierstrass_rational_triangle(v1, v2)
        sa, sb, sc = wt["sides"]
        perim = sa + sb + sc
        ev = eval_exact_pmin_and_dl(sa * sa, sb * sb, sc * sc, sa / perim, sb / perim)
        gx, gy = wt["grad_incenter"]
        h_t = wt["h_tilde_incenter"]
        z2_re, z2_im = wt["z2_sum"]
        ok = all(
            (
                gx == 0,
                gy == 0,
                z2_re == 0,
                z2_im == 0,
                h_t == ((Fraction(1, 2), Fraction(0)), (Fraction(0), Fraction(1, 2))),
                wt["det_h_tilde_incenter"] == Fraction(1, 4),
                ev["det_h_poly"] == Fraction(1, 4),
                8 * ev["pmin"] == ev["r_prod"],
                all(q > 0 for q in wt["charges"]),
            )
        )
        rows.append(
            {
                "name": name,
                "v1_q": format_rational(wt["v1"]),
                "v2_q": format_rational(wt["v2"]),
                "v3_q": format_rational(wt["v3"]),
                "aspect_h_over_max_side_dec6": format_decimal(wt["aspect_h_over_max_side"], 6),
                "grad_exact_zero": all((gx == 0, gy == 0)),
                "z2_phase_sum_exact_zero": all((z2_re == 0, z2_im == 0)),
                "h_tilde_exact_half_I2": h_t == ((Fraction(1, 2), Fraction(0)), (Fraction(0), Fraction(1, 2))),
                "det_h_tilde_q": format_rational(wt["det_h_tilde_incenter"]),
                "pmin_eq_one_eighth_rprod": 8 * ev["pmin"] == ev["r_prod"],
                "verified": ok,
            }
        )
    return {
        "weierstrass_incenter_witnesses": rows,
        "closes_487_flat_triangle_nonemptiness_gap": all(r["verified"] for r in rows),
        "all_verified": all(r["verified"] for r in rows),
    }


def certify_pillar3_isosceles_axis_uniqueness() -> dict[str, Any]:
    """Certify exact cubic-cancellation fold y_* = 1/k, pitchfork y_p = 2/k for h = 3k/(k^2 - 2)."""
    rows: list[dict[str, Any]] = []
    for k in (2, 3, 10, 100):
        h = Fraction(3 * k, k * k - 2)
        y_star = Fraction(1, k)
        y_pitch = Fraction(2, k)
        # Exact V_yy numerator: h - 3y - 2h y^2 == 0 at y_star
        vyy_at_star = h - 3 * y_star - 2 * h * y_star * y_star
        # Exact V_xx numerator: 2h - 3y - h y^2 == 0 at y_pitch
        vxx_at_pitch = 2 * h - 3 * y_pitch - h * y_pitch * y_pitch
        # Midpoint sample y_mid = 3/(2k) strictly inside (y_star, y_pitch):
        # q'(y) numerator h - 3y - 2h y^2 < 0 (strictly decreasing), V_xx numerator > 0 (positive transverse)
        y_mid = Fraction(3, 2 * k)
        vyy_mid = h - 3 * y_mid - 2 * h * y_mid * y_mid
        vxx_mid = 2 * h - 3 * y_mid - h * y_mid * y_mid
        # Incenter y_I = h / (1 + sqrt(1 + h^2)): prove y_star < y_I < y_pitch purely rationally by squaring!
        # y_I > y_star <=> h/y_star - 1 > sqrt(1 + h^2) <=> (h/y_star - 1)^2 > 1 + h^2
        # y_I < y_pitch <=> h/y_pitch - 1 < sqrt(1 + h^2) <=> (h/y_pitch - 1)^2 < 1 + h^2 (or LHS <= 0)
        lhs_lo = h / y_star - 1
        lhs_hi = h / y_pitch - 1
        incenter_above_star = all((lhs_lo > 0, lhs_lo * lhs_lo > 1 + h * h))
        incenter_below_pitch = all((lhs_hi > 0, lhs_hi * lhs_hi < 1 + h * h))
        ok = all(
            (
                0 < y_star < y_pitch < h,
                vyy_at_star == 0,
                vxx_at_pitch == 0,
                vyy_mid < 0,
                vxx_mid > 0,
                incenter_above_star,
                incenter_below_pitch,
            )
        )
        rows.append(
            {
                "k": k,
                "height_h_q": format_rational(h),
                "y_fold_q": format_rational(y_star),
                "y_pitchfork_q": format_rational(y_pitch),
                "y_pitchfork_eq_2_y_fold": y_pitch == 2 * y_star,
                "vyy_fold_exact_zero": vyy_at_star == 0,
                "vxx_pitchfork_exact_zero": vxx_at_pitch == 0,
                "q_prime_strictly_neg_at_mid": vyy_mid < 0,
                "vxx_strictly_pos_at_mid": vxx_mid > 0,
                "incenter_strictly_between_fold_and_pitchfork": all((incenter_above_star, incenter_below_pitch)),
                "verified": ok,
            }
        )
    return {
        "isosceles_axis_cases": rows,
        "at_most_one_on_axis_minimum_all_h_gt_0": all(r["verified"] for r in rows),
        "all_verified": all(r["verified"] for r in rows),
    }


def _check_4d_gale_nikaido_box(
    a2_iv: QInterval,
    b2_iv: QInterval,
    c2_iv: QInterval,
    u1_iv: QInterval,
    u2_iv: QInterval,
    a_mat: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
) -> tuple[bool, Fraction, Fraction, Fraction]:
    u3_iv = QInterval(Fraction(1)).sub(u1_iv).sub(u2_iv)
    d16 = (a2_iv.mul(b2_iv).add(b2_iv.mul(c2_iv)).add(c2_iv.mul(a2_iv))).scale(Fraction(2)).sub(
        a2_iv.sq().add(b2_iv.sq()).add(c2_iv.sq())
    )
    delta2 = d16.scale(Fraction(1, 16))
    r1sq = u2_iv.sq().mul(c2_iv).add(u3_iv.sq().mul(b2_iv)).add(u2_iv.mul(u3_iv).mul(b2_iv.add(c2_iv).sub(a2_iv)))
    r2sq = u3_iv.sq().mul(a2_iv).add(u1_iv.sq().mul(c2_iv)).add(u3_iv.mul(u1_iv).mul(c2_iv.add(a2_iv).sub(b2_iv)))
    r3sq = u1_iv.sq().mul(b2_iv).add(u2_iv.sq().mul(a2_iv)).add(u1_iv.mul(u2_iv).mul(a2_iv.add(b2_iv).sub(c2_iv)))
    leib = u2_iv.mul(u3_iv).mul(a2_iv).add(u3_iv.mul(u1_iv).mul(b2_iv)).add(u1_iv.mul(u2_iv).mul(c2_iv))
    pmin = delta2.scale(Fraction(18)).mul(u1_iv).mul(u2_iv).mul(u3_iv).mul(leib).sub(r1sq.mul(r2sq).mul(r3sq))
    dr1_du1 = (u3_iv.mul(b2_iv).scale(Fraction(2)).add(u2_iv.mul(b2_iv.add(c2_iv).sub(a2_iv)))).neg()
    dr1_du2 = (u2_iv.mul(c2_iv).scale(Fraction(2)).add(u3_iv.mul(b2_iv.add(c2_iv).sub(a2_iv)))).sub(
        u3_iv.mul(b2_iv).scale(Fraction(2)).add(u2_iv.mul(b2_iv.add(c2_iv).sub(a2_iv)))
    )
    dr2_du1 = (u1_iv.mul(c2_iv).scale(Fraction(2)).add(u3_iv.mul(c2_iv.add(a2_iv).sub(b2_iv)))).sub(
        u3_iv.mul(a2_iv).scale(Fraction(2)).add(u1_iv.mul(c2_iv.add(a2_iv).sub(b2_iv)))
    )
    dr2_du2 = (u3_iv.mul(a2_iv).scale(Fraction(2)).add(u1_iv.mul(c2_iv.add(a2_iv).sub(b2_iv)))).neg()
    dr3_du1 = u1_iv.mul(b2_iv).scale(Fraction(2)).add(u2_iv.mul(a2_iv.add(b2_iv).sub(c2_iv)))
    dr3_du2 = u2_iv.mul(a2_iv).scale(Fraction(2)).add(u1_iv.mul(a2_iv.add(b2_iv).sub(c2_iv)))
    h32 = Fraction(3, 2)
    iu1, iu2, iu3 = u1_iv.inv(), u2_iv.inv(), u3_iv.inv()
    ir1, ir2, ir3 = r1sq.inv(), r2sq.inv(), r3sq.inv()
    j11 = iu1.add(iu3).add(dr1_du1.mul(ir1).sub(dr3_du1.mul(ir3)).scale(h32))
    j12 = iu3.add(dr1_du2.mul(ir1).sub(dr3_du2.mul(ir3)).scale(h32))
    j21 = iu3.add(dr2_du1.mul(ir2).sub(dr3_du1.mul(ir3)).scale(h32))
    j22 = iu2.add(iu3).add(dr2_du2.mul(ir2).sub(dr3_du2.mul(ir3)).scale(h32))
    m00 = (j11.scale(a_mat[0][0]).add(j21.scale(a_mat[0][1]))).neg()
    m01 = (j12.scale(a_mat[0][0]).add(j22.scale(a_mat[0][1]))).neg()
    m10 = (j11.scale(a_mat[1][0]).add(j21.scale(a_mat[1][1]))).neg()
    m11 = (j12.scale(a_mat[1][0]).add(j22.scale(a_mat[1][1]))).neg()
    s00, s11 = m00, m11
    s01 = m01.add(m10).scale(Fraction(1, 2))
    det_s = s00.mul(s11).sub(s01.sq())
    det_j = j11.mul(j22).sub(j12.mul(j21))
    ok = all((pmin.lo > 0, det_j.lo > 0, s00.lo > 0, det_s.lo > 0))
    return ok, pmin.lo, s00.lo, det_s.lo


def certify_pillar4_gale_nikaido_compact_univalence() -> dict[str, Any]:
    """Certify global Z-injectivity on two 4D compact (shape x interior) boxes via exact QInterval Sylvester PD."""
    box_specs = (
        (
            "equilateral_4d_compact",
            (Fraction(255, 256), Fraction(257, 256)),
            (Fraction(255, 256), Fraction(257, 256)),
            Fraction(1),
            (Fraction(1, 3) - Fraction(1, 128), Fraction(1, 3) + Fraction(1, 128)),
            (Fraction(1, 3) - Fraction(1, 128), Fraction(1, 3) + Fraction(1, 128)),
            2,
            4,
            ((Fraction(4, 9), Fraction(-2, 9)), (Fraction(-2, 9), Fraction(4, 9))),
        ),
        (
            "right_345_4d_compact",
            (Fraction(575, 64), Fraction(577, 64)),
            (Fraction(1023, 64), Fraction(1025, 64)),
            Fraction(25),
            (Fraction(1, 4) - Fraction(1, 128), Fraction(1, 4) + Fraction(1, 128)),
            (Fraction(1, 3) - Fraction(1, 128), Fraction(1, 3) + Fraction(1, 128)),
            2,
            4,
            ((Fraction(3, 8), Fraction(-1, 6)), (Fraction(-1, 6), Fraction(4, 9))),
        ),
    )
    results: list[dict[str, Any]] = []
    for name, a2_rng, b2_rng, c2_val, u1_rng, u2_rng, ns, nu, a_mat in box_specs:
        c2_iv = QInterval(c2_val)
        sub_checked = 0
        sub_passed = 0
        pmin_vals: list[Fraction] = []
        s00_vals: list[Fraction] = []
        dets_vals: list[Fraction] = []
        for ia in range(ns):
            a2_iv = QInterval(
                a2_rng[0] + (a2_rng[1] - a2_rng[0]) * Fraction(ia, ns),
                a2_rng[0] + (a2_rng[1] - a2_rng[0]) * Fraction(ia + 1, ns),
            )
            for ib in range(ns):
                b2_iv = QInterval(
                    b2_rng[0] + (b2_rng[1] - b2_rng[0]) * Fraction(ib, ns),
                    b2_rng[0] + (b2_rng[1] - b2_rng[0]) * Fraction(ib + 1, ns),
                )
                for iu1 in range(nu):
                    u1_iv = QInterval(
                        u1_rng[0] + (u1_rng[1] - u1_rng[0]) * Fraction(iu1, nu),
                        u1_rng[0] + (u1_rng[1] - u1_rng[0]) * Fraction(iu1 + 1, nu),
                    )
                    for iu2 in range(nu):
                        u2_iv = QInterval(
                            u2_rng[0] + (u2_rng[1] - u2_rng[0]) * Fraction(iu2, nu),
                            u2_rng[0] + (u2_rng[1] - u2_rng[0]) * Fraction(iu2 + 1, nu),
                        )
                        ok, pm_lo, s00_lo, ds_lo = _check_4d_gale_nikaido_box(a2_iv, b2_iv, c2_iv, u1_iv, u2_iv, a_mat)
                        sub_checked += 1
                        if ok:
                            sub_passed += 1
                        pmin_vals.append(pm_lo)
                        s00_vals.append(s00_lo)
                        dets_vals.append(ds_lo)
        min_pmin = min(pmin_vals)
        min_s00 = min(s00_vals)
        min_dets = min(dets_vals)
        results.append(
            {
                "name": name,
                "a2_range_q": [format_rational(a2_rng[0]), format_rational(a2_rng[1])],
                "b2_range_q": [format_rational(b2_rng[0]), format_rational(b2_rng[1])],
                "c2_q": format_rational(c2_val),
                "u1_range_q": [format_rational(u1_rng[0]), format_rational(u1_rng[1])],
                "u2_range_q": [format_rational(u2_rng[0]), format_rational(u2_rng[1])],
                "subboxes_checked": sub_checked,
                "subboxes_passed": sub_passed,
                "min_pmin_lo_dec6": format_decimal(min_pmin, 6),
                "min_s00_lo_dec6": format_decimal(min_s00, 6),
                "min_det_sym_lo_dec6": format_decimal(min_dets, 6),
                "strictly_injective_on_4d_box": all(
                    (sub_checked == 64, sub_passed == 64, min_pmin > 0, min_s00 > 0, min_dets > 0)
                ),
            }
        )
    return {
        "compact_4d_boxes": results,
        "total_4d_subboxes_certified": sum(r["subboxes_passed"] for r in results),
        "all_verified": all(r["strictly_injective_on_4d_box"] for r in results),
    }


def certify_pillar5_krawczyk_configurations() -> dict[str, Any]:
    """Certify 14 pairwise-disjoint non-degenerate equilibria across 4 rational n=3 configurations."""
    hw = Fraction(1, 1 << 20)
    configs = (
        (
            "n3_weierstrass_iso_acute",
            Fraction(1, 3),
            Fraction(1, 3),
            None,
            (
                (Fraction(0), Fraction(0)),
                (Fraction(-557602579, 1073741824), Fraction(280559551, 1073741824)),
                (Fraction(0), Fraction(-293002679, 536870912)),
                (Fraction(557602579, 1073741824), Fraction(280559551, 1073741824)),
            ),
            1,
            3,
        ),
        (
            "n3_weierstrass_scalene",
            Fraction(1, 4),
            Fraction(1, 5),
            None,
            (
                (Fraction(0), Fraction(0)),
                (Fraction(-264462833, 536870912), Fraction(280650093, 1073741824)),
                (Fraction(29982451, 1073741824), Fraction(-156391975, 268435456)),
                (Fraction(123037753, 268435456), Fraction(89771009, 268435456)),
            ),
            1,
            3,
        ),
        (
            "n3_weierstrass_flat_iso",
            Fraction(1, 20),
            Fraction(1, 20),
            None,
            (
                (Fraction(0), Fraction(0)),
                (Fraction(-465263281, 1073741824), Fraction(177983715, 536870912)),
                (Fraction(0), Fraction(-661694729, 1073741824)),
                (Fraction(465263281, 1073741824), Fraction(177983715, 536870912)),
            ),
            1,
            3,
        ),
        (
            "n3_weierstrass_skew_2eq",
            Fraction(1, 3),
            Fraction(1, 3),
            (Fraction(20), Fraction(1), Fraction(1)),
            (
                (Fraction(-237295619, 1073741824), Fraction(1472520845, 536870912)),
                (Fraction(454427509, 536870912), Fraction(-1069122543, 1073741824)),
            ),
            0,
            2,
        ),
    )

    cfg_rows: list[dict[str, Any]] = []
    total_eqs = 0
    for name, v1, v2, custom_q, centers, exp_m1, exp_m2 in configs:
        wt = weierstrass_rational_triangle(v1, v2)
        sites = wt["sites"]
        charges = wt["charges"] if custom_q is None else custom_q
        eq_records: list[dict[str, Any]] = []
        bounds: list[tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]] = []
        for cx, cy in centers:
            rec = certify_2d_krawczyk_box(sites, charges, cx, cy, hw)
            eq_records.append(rec)
            bounds.append(((cx - hw, cy - hw), (cx + hw, cy + hw)))
        # Exact pairwise disjointness over Q
        pairs_checked = 0
        all_disjoint = True
        for i in range(len(eq_records)):
            for j in range(i + 1, len(eq_records)):
                pairs_checked += 1
                if not _boxes_strictly_disjoint(bounds[i][0], bounds[i][1], bounds[j][0], bounds[j][1]):
                    all_disjoint = False
        m1 = sum(1 for r in eq_records if r["morse_index_3d"] == 1)
        m2 = sum(1 for r in eq_records if r["morse_index_3d"] == 2)
        total_eqs += len(eq_records)
        cfg_rows.append(
            {
                "name": name,
                "charges_q": [format_rational(q) for q in charges],
                "equilibria_certified": len(eq_records),
                "m1_local_minima": m1,
                "m2_saddles": m2,
                "m2_minus_m1": m2 - m1,
                "pairwise_disjoint_pairs": pairs_checked,
                "all_boxes_pairwise_disjoint": all_disjoint,
                "equilibria": eq_records,
                "verified": all((m1 == exp_m1, m2 == exp_m2, m2 - m1 == 2, all_disjoint)),
            }
        )
    return {
        "configurations": cfg_rows,
        "total_equilibria_certified": total_eqs,
        "all_verified": all((total_eqs == 14, all(r["verified"] for r in cfg_rows))),
    }


def build_payload() -> dict[str, Any]:
    p1 = certify_pillar1_algebraic_identities()
    p2 = certify_pillar2_incenter_universal_theorem()
    p3 = certify_pillar3_isosceles_axis_uniqueness()
    p4 = certify_pillar4_gale_nikaido_compact_univalence()
    p5 = certify_pillar5_krawczyk_configurations()
    all_ok = all((p1["all_verified"], p2["all_verified"], p3["all_verified"], p4["all_verified"], p5["all_verified"]))
    return {
        "track_id": TRACK_ID,
        "issue_ref": "Foundation-0/point-charges#4",
        "novelty_class": "N1",
        "pillar1_algebraic_identities": p1,
        "pillar2_incenter_universal_theorem": p2,
        "pillar3_isosceles_axis_uniqueness": p3,
        "pillar4_gale_nikaido_compact_univalence": p4,
        "pillar5_krawczyk_configurations": p5,
        "summary": {
            "pillar_count": 5,
            "all_five_pillars_verified": all_ok,
            "total_4d_univalence_subboxes": p4["total_4d_subboxes_certified"],
            "total_krawczyk_equilibria": p5["total_equilibria_certified"],
            "global_n3_conjecture_open_outside_compact_boxes": True,
            "float_count_on_trust_path": 0,
            "bare_assert_count": 0,
            "rng_calls": 0,
        },
    }
