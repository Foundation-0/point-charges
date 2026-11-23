#!/usr/bin/env python3
"""
Track S21d: Pure-Z/Q Krawczyk & Bareiss Sylvester-Inertia Certifier for
Positive Point-Charge Equilibria at the n=7 (#12) and n=8 (#5) Maxwell Frontier.

Certifies four exact rational positive-charge configurations in R^3 with zero
floating-point operations, zero external dependencies, and zero bare `assert`s:
  1. `n2-analytic-control`: n=2, 1 equilibrium (0, 0, 1, 0), exact control.
  2. `n7-abk-second-pair` (Foundation-0/point-charges#12): n=7, 44 equilibria
     (0, 19, 25, 0), m2 - m1 = 6, 44 > 36 = (7-1)^2 (where float64 certify.py
     OOMs at 53 GB / 60M boxes due to 53-bit cancellation).
  3. `n6-abk-milnor8-flat` (n=6, 31 > 29 > 25 = (6-1)^2): 31 equilibria
     (0, 13, 18, 0), m2 - m1 = 5.
  4. `n8-abk-milnor8-flat` (Foundation-0/point-charges#5): n=8, 51 equilibria
     (0, 22, 29, 0), m2 - m1 = 7, 51 > 49 = (8-1)^2 — pure-Z/Q machine-verified
     counterexample to Maxwell's 1873 bound (n-1)^2 for n=8 positive point
     charges (novelty N1, crediting Arathoon-Ball-Kvalheim Jul 2026 and
     Shiva Kintali Aug 2026).

Scope honesty (Law 4): Local Krawczyk uniqueness proves exact lower bounds on
the number of pairwise-disjoint non-degenerate equilibria (N >= 44, N >= 31,
N >= 51 > 49), which is mathematically sufficient to refute Maxwell's upper
bound for n=8, without claiming global completeness over R^3 minus {a_j}.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import hmac
import json
import math
import os
from pathlib import Path
import stat
from typing import Any

SCHEMA_VERSION = "s21d-point-charges-n7-n8-cert/v1"
TRACK_ID = "TRACK-S21D"
PREC_BITS = 384
ONE_SCALED = 1 << PREC_BITS
Y_ROUND_BITS = 160
MAX_JSON_DEPTH = 64
MAX_JSON_BYTES = 2_000_000
MAX_INT_BITS = 2048
MAX_INT_DIGITS = 620
MAX_CONFIGS = 16
MAX_CHARGES = 32
MAX_BOXES = 256
_HEX_LOWER = frozenset("0123456789abcdef")
_CANONICAL_CONFIG_NAMES = (
    "n2-analytic-control",
    "n7-abk-second-pair",
    "n6-abk-milnor8-flat",
    "n8-abk-milnor8-flat",
)


class PointChargeFrontierError(ValueError):
    """Fail-closed typed error (never stripped under `python3 -O`)."""

    def __init__(self, code: str, message: str) -> None:
        self.code = code
        self.message = message
        super().__init__(f"{code}: {message}")


# ---------------------------------------------------------------------------
# Zone-B strict primitive & canonical JSON I/O helpers (zero floats, zero DoS)
# ---------------------------------------------------------------------------


def _require_valid_utf8_str(s: str, label: str = "str") -> str:
    try:
        s.encode("utf-8")
    except UnicodeEncodeError as exc:
        raise PointChargeFrontierError("INVALID_UTF8", f"{label} contains unpaired UTF-16 surrogate") from exc
    return s


def require_int(
    val: Any,
    label: str = "int",
    *,
    lo: int | None = None,
    hi: int | None = None,
) -> int:
    if type(val) is not int:
        raise PointChargeFrontierError("INVALID_INT_TYPE", f"{label} must be strict int, got {type(val).__name__}")
    blen = val.bit_length()
    if blen > MAX_INT_BITS:
        raise PointChargeFrontierError("INT_BIT_LIMIT", f"{label} bit_length={blen} exceeds {MAX_INT_BITS}")
    if lo is not None:
        if type(lo) is not int:
            raise PointChargeFrontierError("INVALID_INT_BOUND", f"{label} lo must be strict int")
        if val < lo:
            raise PointChargeFrontierError("INT_BELOW_MIN", f"{label} bit_length={blen} below min_bits={lo.bit_length()}")
    if hi is not None:
        if type(hi) is not int:
            raise PointChargeFrontierError("INVALID_INT_BOUND", f"{label} hi must be strict int")
        if val > hi:
            raise PointChargeFrontierError("INT_ABOVE_MAX", f"{label} bit_length={blen} above max_bits={hi.bit_length()}")
    return val


def require_bool(val: Any, label: str = "bool") -> bool:
    if type(val) is not bool:
        raise PointChargeFrontierError("INVALID_BOOL_TYPE", f"{label} must be strict bool, got {type(val).__name__}")
    return val


def require_str(val: Any, label: str = "str", *, min_len: int = 1, max_len: int = 1024) -> str:
    if type(val) is not str:
        raise PointChargeFrontierError("INVALID_STR_TYPE", f"{label} must be strict str, got {type(val).__name__}")
    _require_valid_utf8_str(val, label)
    n = len(val)
    if n < min_len:
        raise PointChargeFrontierError("STR_TOO_SHORT", f"{label} length {n} < min_bits={min_len.bit_length()}")
    if n > max_len:
        raise PointChargeFrontierError("STR_TOO_LONG", f"{label} length {n} > max_bits={max_len.bit_length()}")
    return val


def require_seq(
    val: Any,
    label: str = "seq",
    *,
    min_len: int = 1,
    max_len: int = MAX_BOXES,
) -> tuple[Any, ...]:
    if type(val) not in (list, tuple):
        raise PointChargeFrontierError("INVALID_SEQ_TYPE", f"{label} must be strict list or tuple")
    n = len(val)
    if n < min_len:
        raise PointChargeFrontierError("SEQ_TOO_SHORT", f"{label} length {n} < min_bits={min_len.bit_length()}")
    if n > max_len:
        raise PointChargeFrontierError("SEQ_TOO_LONG", f"{label} length {n} > max_bits={max_len.bit_length()}")
    return tuple(val)


def require_mapping(val: Any, label: str = "mapping") -> dict[str, Any]:
    if type(val) is not dict:
        raise PointChargeFrontierError("INVALID_MAPPING_TYPE", f"{label} must be strict dict")
    for k in val:
        if type(k) is not str:
            raise PointChargeFrontierError("NON_STR_KEY", f"{label} keys must be strict str")
        _require_valid_utf8_str(k, f"{label}.key")
    return val


def require_hex64_ascii(val: Any, label: str = "sha256") -> str:
    if type(val) is not str:
        raise PointChargeFrontierError("INVALID_DIGEST_TYPE", f"{label} must be strict str")
    if len(val) != 64:
        raise PointChargeFrontierError("INVALID_DIGEST_LEN", f"{label} must have length 64")
    for ch in val:
        if ch not in _HEX_LOWER:
            raise PointChargeFrontierError("INVALID_DIGEST_HEX", f"{label} must be lowercase ASCII hex")
    return val


def constant_time_hex64_equal(left: Any, right: Any, label: str = "sha256") -> bool:
    a = require_hex64_ascii(left, f"{label}.left")
    b = require_hex64_ascii(right, f"{label}.right")
    return hmac.compare_digest(a, b)


def _reject_float_token(tok: str) -> Any:
    raise PointChargeFrontierError("JSON_FLOAT_FORBIDDEN", f"float literal forbidden: {tok!r}")


def _reject_constant_token(tok: str) -> Any:
    raise PointChargeFrontierError("JSON_CONSTANT_FORBIDDEN", f"non-RFC8259 constant forbidden: {tok!r}")


def _parse_capped_int_token(tok: str) -> int:
    digits = tok[1:] if tok.startswith("-") else tok
    if not digits or not digits.isascii() or not digits.isdigit():
        raise PointChargeFrontierError("INVALID_INT_TOKEN", "integer token must be ASCII digits")
    if len(digits) > MAX_INT_DIGITS:
        raise PointChargeFrontierError("INT_DIGIT_LIMIT", f"integer token exceeds {MAX_INT_DIGITS} digits")
    val = int(tok)
    blen = val.bit_length()
    if blen > MAX_INT_BITS:
        raise PointChargeFrontierError("INT_BIT_LIMIT", f"integer token bit_length={blen} exceeds {MAX_INT_BITS}")
    return val


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k, v in pairs:
        if k in out:
            raise PointChargeFrontierError("JSON_DUPLICATE_KEY", f"duplicate JSON key: {k!r}")
        out[k] = v
    return out


def fmt_q(q: Fraction) -> str:
    if type(q) is not Fraction:
        raise PointChargeFrontierError("INVALID_FRACTION_TYPE", "q must be strict Fraction")
    nb = q.numerator.bit_length()
    db = q.denominator.bit_length()
    if nb > MAX_INT_BITS or db > MAX_INT_BITS:
        raise PointChargeFrontierError("FRACTION_BIT_LIMIT", f"fraction bit_lengths=({nb},{db}) exceed {MAX_INT_BITS}")
    if q.denominator == 1:
        return str(q.numerator)
    return f"{q.numerator}/{q.denominator}"


def parse_q_str(val: Any, label: str = "q") -> Fraction:
    s = require_str(val, label, min_len=1, max_len=1300)
    parts = s.split("/")
    if len(parts) == 1:
        return Fraction(_parse_capped_int_token(parts[0]), 1)
    if len(parts) == 2:
        num = _parse_capped_int_token(parts[0])
        den = _parse_capped_int_token(parts[1])
        if den <= 0:
            raise PointChargeFrontierError("INVALID_FRACTION_DEN", f"{label} denominator must be positive")
        return Fraction(num, den)
    raise PointChargeFrontierError("INVALID_FRACTION_STR", f"{label} malformed rational string")


def validate_json_value(obj: Any, *, max_depth: int = MAX_JSON_DEPTH) -> None:
    limit = require_int(max_depth, "max_depth", lo=1, hi=128)
    stack: list[tuple[Any, int]] = [(obj, 0)]
    while stack:
        cur, depth = stack.pop()
        if depth > limit:
            raise PointChargeFrontierError("MAX_DEPTH_EXCEEDED", f"JSON nesting exceeds {limit}")
        t = type(cur)
        if t is dict:
            for k, v in cur.items():
                if type(k) is not str:
                    raise PointChargeFrontierError("NON_STR_KEY", "JSON object keys must be strict str")
                _require_valid_utf8_str(k, "json_key")
                stack.append((v, depth + 1))
        elif t is list:
            for item in cur:
                stack.append((item, depth + 1))
        elif t is int:
            blen = cur.bit_length()
            if blen > MAX_INT_BITS:
                raise PointChargeFrontierError("INT_BIT_LIMIT", f"JSON int bit_length={blen} exceeds {MAX_INT_BITS}")
        elif t is str:
            _require_valid_utf8_str(cur, "json_str")
        elif t is bool or cur is None:
            continue
        elif t is float:
            raise PointChargeFrontierError("FLOAT_FORBIDDEN", "float values are forbidden in exact JSON")
        else:
            raise PointChargeFrontierError("UNSUPPORTED_JSON_TYPE", f"unsupported JSON type {t.__name__}")


def parse_strict_json_text(raw_text: str, *, require_object: bool = True) -> Any:
    req_obj = require_bool(require_object, "require_object")
    if type(raw_text) is not str:
        raise PointChargeFrontierError("INVALID_JSON_INPUT", "JSON text must be strict str")
    encoded = _require_valid_utf8_str(raw_text, "json_text").encode("utf-8")
    if len(encoded) > MAX_JSON_BYTES:
        raise PointChargeFrontierError("JSON_TOO_LARGE", f"JSON exceeds {MAX_JSON_BYTES} bytes")
    if not raw_text.strip():
        raise PointChargeFrontierError("EMPTY_JSON_INPUT", "JSON text must be non-empty")
    try:
        obj = json.loads(
            raw_text,
            parse_float=_reject_float_token,
            parse_constant=_reject_constant_token,
            parse_int=_parse_capped_int_token,
            object_pairs_hook=_reject_duplicate_keys,
        )
    except PointChargeFrontierError:
        raise
    except (json.JSONDecodeError, ValueError, RecursionError) as exc:
        raise PointChargeFrontierError("JSON_PARSE_ERROR", f"invalid JSON: {type(exc).__name__}") from exc
    if req_obj and type(obj) is not dict:
        raise PointChargeFrontierError("JSON_ROOT_NOT_OBJECT", "JSON root must be a strict object")
    validate_json_value(obj)
    return obj


def canonical_json_bytes(obj: Any) -> bytes:
    validate_json_value(obj)
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("ascii")


def sha256_canonical(obj: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(obj)).hexdigest()


def require_safe_path(path: Path, label: str = "path") -> Path:
    if type(path) is not type(Path()):
        raise PointChargeFrontierError("INVALID_PATH_TYPE", f"{label} must be a strict Path")
    raw_str = str(path)
    if "\x00" in raw_str:
        raise PointChargeFrontierError("INVALID_PATH_NUL", f"{label} contains embedded NUL byte")
    _require_valid_utf8_str(raw_str, label)
    if not path.is_absolute():
        pwd_env = os.environ.get("PWD")
        if pwd_env and os.path.realpath(pwd_env) == os.path.realpath(os.getcwd()):
            pwd_path = Path(pwd_env)
            if pwd_path.is_symlink():
                raise PointChargeFrontierError("SYMLINK_ANCESTOR_FORBIDDEN", f"{label} has symlink PWD")
            for p_anc in pwd_path.parents:
                if p_anc.is_symlink():
                    raise PointChargeFrontierError("SYMLINK_ANCESTOR_FORBIDDEN", f"{label} has symlink PWD ancestor")
    if path.is_symlink():
        raise PointChargeFrontierError("SYMLINK_FORBIDDEN", f"{label} must not be a symlink")
    for parent in path.parents:
        if parent.is_symlink():
            raise PointChargeFrontierError("SYMLINK_ANCESTOR_FORBIDDEN", f"{label} has symlink ancestor")
    return path


def load_safe_bytes(path: Path, label: str = "file") -> bytes:
    require_safe_path(path, label)
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
    try:
        fd = os.open(path, flags)
    except OSError as exc:
        raise PointChargeFrontierError("FILE_READ_ERROR", f"{label} open failed: {exc.strerror}") from exc
    try:
        st = os.fstat(fd)
        if not stat.S_ISREG(st.st_mode):
            raise PointChargeFrontierError("NON_REGULAR_FILE", f"{label} must be a regular file")
        chunks: list[bytes] = []
        total = 0
        while True:
            chunk = os.read(fd, min(65536, MAX_JSON_BYTES + 1 - total))
            if not chunk:
                break
            total += len(chunk)
            if total > MAX_JSON_BYTES:
                raise PointChargeFrontierError("FILE_TOO_LARGE", f"{label} exceeds {MAX_JSON_BYTES} bytes")
            chunks.append(chunk)
        return b"".join(chunks)
    except OSError as exc:
        raise PointChargeFrontierError("FILE_READ_ERROR", f"{label} read failed: {exc.strerror}") from exc
    finally:
        os.close(fd)


def safe_write_bytes(path: Path, data: bytes, label: str = "out") -> None:
    require_safe_path(path, label)
    if type(data) is not bytes:
        raise PointChargeFrontierError("INVALID_BYTES_TYPE", f"{label} data must be strict bytes")
    flags = os.O_WRONLY | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
    try:
        fd = os.open(path, flags, 0o600)
    except OSError as exc:
        raise PointChargeFrontierError("FILE_WRITE_ERROR", f"{label} open failed: {exc.strerror}") from exc
    try:
        st = os.fstat(fd)
        if not stat.S_ISREG(st.st_mode):
            raise PointChargeFrontierError("NON_REGULAR_FILE", f"{label} must be a regular file")
        if st.st_nlink > 1:
            raise PointChargeFrontierError("HARDLINK_FORBIDDEN", f"{label} must not be a hardlink")
        os.ftruncate(fd, 0)
        view = memoryview(data)
        written = 0
        while written < len(data):
            written += os.write(fd, view[written:])
    except OSError as exc:
        raise PointChargeFrontierError("FILE_WRITE_ERROR", f"{label} write failed: {exc.strerror}") from exc
    finally:
        os.close(fd)


# ---------------------------------------------------------------------------
# Pure-Z/Q 384-Bit Outward-Rounded Dyadic Interval Arithmetic ([lo, hi] / 2^P)
# ---------------------------------------------------------------------------


def iv_from_q(q: Fraction) -> tuple[int, int]:
    """Outward-rounded 384-bit dyadic enclosure of an exact rational Fraction."""
    if type(q) is not Fraction:
        raise PointChargeFrontierError("INVALID_FRACTION_TYPE", "q must be strict Fraction")
    num = q.numerator << PREC_BITS
    den = q.denominator
    return (num // den, -((-num) // den))


def iv_add(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
    return (a[0] + b[0], a[1] + b[1])


def iv_sub(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
    return (a[0] - b[1], a[1] - b[0])


def iv_neg(a: tuple[int, int]) -> tuple[int, int]:
    return (-a[1], -a[0])


def iv_mul(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
    p1 = a[0] * b[0]
    p2 = a[0] * b[1]
    p3 = a[1] * b[0]
    p4 = a[1] * b[1]
    lo = min(p1, p2, p3, p4) >> PREC_BITS
    hi = -((-max(p1, p2, p3, p4)) >> PREC_BITS)
    return (lo, hi)


def iv_sqr(a: tuple[int, int]) -> tuple[int, int]:
    lo, hi = a
    if lo >= 0:
        return ((lo * lo) >> PREC_BITS, -((-(hi * hi)) >> PREC_BITS))
    if hi <= 0:
        return ((hi * hi) >> PREC_BITS, -((-(lo * lo)) >> PREC_BITS))
    m = max(-lo, hi)
    return (0, -((-(m * m)) >> PREC_BITS))


def iv_mul_int(a: tuple[int, int], k: int) -> tuple[int, int]:
    if k >= 0:
        return (a[0] * k, a[1] * k)
    return (a[1] * k, a[0] * k)


def iv_inv_pos(a: tuple[int, int]) -> tuple[int, int]:
    if a[0] <= 0:
        raise PointChargeFrontierError("NON_POSITIVE_INTERVAL", "interval lower bound must be strictly positive")
    num = 1 << (2 * PREC_BITS)
    return (num // a[1], -((-num) // a[0]))


def iv_sqrt_pos(a: tuple[int, int]) -> tuple[int, int]:
    if a[0] <= 0:
        raise PointChargeFrontierError("NON_POSITIVE_INTERVAL", "interval lower bound must be strictly positive for sqrt")
    s_lo = math.isqrt(a[0] << PREC_BITS)
    hi_shifted = a[1] << PREC_BITS
    s_hi_base = math.isqrt(hi_shifted)
    s_hi = s_hi_base if s_hi_base * s_hi_base == hi_shifted else s_hi_base + 1
    return (s_lo, s_hi)


def iv_abs_max(a: tuple[int, int]) -> int:
    return max(abs(a[0]), abs(a[1]))


def det_3x3_iv(m: list[list[tuple[int, int]]]) -> tuple[int, int]:
    t0 = iv_mul(m[0][0], iv_sub(iv_mul(m[1][1], m[2][2]), iv_mul(m[1][2], m[2][1])))
    t1 = iv_mul(m[0][1], iv_sub(iv_mul(m[1][0], m[2][2]), iv_mul(m[1][2], m[2][0])))
    t2 = iv_mul(m[0][2], iv_sub(iv_mul(m[1][0], m[2][1]), iv_mul(m[1][1], m[2][0])))
    return iv_add(iv_sub(t0, t1), t2)


def eval_grad_hess_iv(
    box: list[tuple[int, int]],
    sites_iv: list[list[tuple[int, int]]],
    charges_iv: list[tuple[int, int]],
) -> tuple[list[tuple[int, int]], list[list[tuple[int, int]]]]:
    """Outward-rounded 384-bit interval enclosure of Coulomb gradient and Hessian."""
    g: list[tuple[int, int]] = [(0, 0), (0, 0), (0, 0)]
    h: list[list[tuple[int, int]]] = [[(0, 0), (0, 0), (0, 0)] for _ in range(3)]
    for site, q in zip(sites_iv, charges_iv):
        dx = iv_sub(box[0], site[0])
        dy = iv_sub(box[1], site[1])
        dz = iv_sub(box[2], site[2])
        d = (dx, dy, dz)
        r2 = iv_add(iv_add(iv_sqr(dx), iv_sqr(dy)), iv_sqr(dz))
        inv_r2 = iv_inv_pos(r2)
        r = iv_sqrt_pos(r2)
        inv_r = iv_inv_pos(r)
        inv_r3 = iv_mul(inv_r2, inv_r)
        inv_r5 = iv_mul(inv_r3, inv_r2)
        q_inv_r3 = iv_mul(q, inv_r3)
        q_inv_r5 = iv_mul(q, inv_r5)
        three_q_inv_r5 = iv_mul_int(q_inv_r5, 3)
        for i in range(3):
            g[i] = iv_sub(g[i], iv_mul(q_inv_r3, d[i]))
        for i in range(3):
            for j in range(3):
                term = iv_mul(three_q_inv_r5, iv_mul(d[i], d[j]))
                if i == j:
                    term = iv_sub(term, q_inv_r3)
                h[i][j] = iv_add(h[i][j], term)
    return g, h


# ---------------------------------------------------------------------------
# Exact Rational Linear Algebra & Self-Verified Lagrange / Bareiss Inertia
# ---------------------------------------------------------------------------


def bareiss_det_3x3_int(mat: list[list[int]] | tuple[tuple[int, ...], ...]) -> int:
    """Exact 3x3 fraction-free integer determinant (no division)."""
    rows = require_seq(mat, "mat", min_len=3, max_len=3)
    m: list[list[int]] = []
    for r_i in range(3):
        row = require_seq(rows[r_i], f"mat[{r_i}]", min_len=3, max_len=3)
        m.append([require_int(row[c_i], f"mat[{r_i}][{c_i}]") for c_i in range(3)])
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def exact_inv_3x3_q(mat: list[list[Fraction]]) -> list[list[Fraction]]:
    """Exact Gauss-Jordan inverse of a 3x3 rational matrix."""
    aug = [mat[i][:] + [Fraction(1 if i == j else 0) for j in range(3)] for i in range(3)]
    for col in range(3):
        piv = col
        while piv < 3 and aug[piv][col] == 0:
            piv += 1
        if piv == 3:
            raise PointChargeFrontierError("SINGULAR_MATRIX", "3x3 matrix is singular")
        aug[col], aug[piv] = aug[piv], aug[col]
        sc = aug[col][col]
        aug[col] = [x / sc for x in aug[col]]
        for r in range(3):
            if r != col and aug[r][col] != 0:
                fac = aug[r][col]
                aug[r] = [aug[r][c] - fac * aug[col][c] for c in range(6)]
    return [row[3:] for row in aug]


def verify_congruence_diagonalization(
    mat: list[list[Fraction]],
    s_mat: list[list[Fraction]],
) -> tuple[int, int, int]:
    """Verify Bareiss det(S) != 0 and S^T M S strictly diagonal, returning (n_pos, n_neg, n_zero)."""
    lcm_den = 1
    for r in range(3):
        for c in range(3):
            lcm_den = math.lcm(lcm_den, s_mat[r][c].denominator)
    s_int = [[int(s_mat[r][c] * lcm_den) for c in range(3)] for r in range(3)]
    det_s = (
        s_int[0][0] * (s_int[1][1] * s_int[2][2] - s_int[1][2] * s_int[2][1])
        - s_int[0][1] * (s_int[1][0] * s_int[2][2] - s_int[1][2] * s_int[2][0])
        + s_int[0][2] * (s_int[1][0] * s_int[2][1] - s_int[1][1] * s_int[2][0])
    )
    if det_s == 0:
        raise PointChargeFrontierError("SINGULAR_CONGRUENCE", "congruence matrix S is singular")
    st_m_s = [
        [
            sum(s_mat[p][i] * mat[p][q] * s_mat[q][j] for p in range(3) for q in range(3))
            for j in range(3)
        ]
        for i in range(3)
    ]
    for i in range(3):
        for j in range(3):
            if i != j and st_m_s[i][j] != 0:
                raise PointChargeFrontierError("CONGRUENCE_NOT_DIAGONAL", "S^T M S has non-zero off-diagonal entry")
    n_pos = sum(1 for i in range(3) if st_m_s[i][i] > 0)
    n_neg = sum(1 for i in range(3) if st_m_s[i][i] < 0)
    n_zero = 3 - n_pos - n_neg
    return (n_pos, n_neg, n_zero)


def exact_lagrange_inertia_3x3_q(mat: list[list[Fraction]]) -> tuple[int, int, int]:
    """Exact Sylvester inertia (n_pos, n_neg, n_zero) via Lagrange congruence with self-verified S^T M S = diag(d) and det(S) != 0."""
    rows = require_seq(mat, "mat", min_len=3, max_len=3)
    m_orig: list[list[Fraction]] = []
    for r_i in range(3):
        row = require_seq(rows[r_i], f"mat[{r_i}]", min_len=3, max_len=3)
        parsed: list[Fraction] = []
        for c_i in range(3):
            v = row[c_i]
            if type(v) is not Fraction:
                raise PointChargeFrontierError("INVALID_FRACTION_TYPE", "matrix entries must be strict Fraction")
            parsed.append(v)
        m_orig.append(parsed)
    if m_orig[0][1] != m_orig[1][0] or m_orig[0][2] != m_orig[2][0] or m_orig[1][2] != m_orig[2][1]:
        raise PointChargeFrontierError("NON_SYMMETRIC_MATRIX", "matrix must be symmetric")

    a = [row[:] for row in m_orig]
    s = [[Fraction(1 if i == j else 0) for j in range(3)] for i in range(3)]
    idx = 0
    while idx < 3:
        piv = -1
        for k in range(idx, 3):
            if a[k][k] != 0:
                piv = k
                break
        if piv != -1:
            if piv != idx:
                a[idx], a[piv] = a[piv], a[idx]
                for r in range(3):
                    a[r][idx], a[r][piv] = a[r][piv], a[r][idx]
                    s[r][idx], s[r][piv] = s[r][piv], s[r][idx]
            d = a[idx][idx]
            a_col = [a[r][idx] for r in range(3)]
            for c in range(idx + 1, 3):
                fac = a_col[c] / d
                if fac != 0:
                    for r in range(3):
                        s[r][c] -= fac * s[r][idx]
                    for r in range(idx + 1, 3):
                        a[r][c] -= fac * a_col[r]
            for c in range(idx + 1, 3):
                a[idx][c] = Fraction(0)
                a[c][idx] = Fraction(0)
            idx += 1
            continue
        pair: tuple[int, int] | None = None
        for r in range(idx, 3):
            for c in range(r + 1, 3):
                if a[r][c] != 0:
                    pair = (r, c)
                    break
            if pair is not None:
                break
        if pair is None:
            break
        r0, c0 = pair
        for k in range(3):
            s[k][r0] += s[k][c0]
        for k in range(3):
            a[k][r0] += a[k][c0]
        for k in range(3):
            a[r0][k] += a[c0][k]

    return verify_congruence_diagonalization(m_orig, s)


def _round_q_dyadic(q: Fraction, bits: int) -> Fraction:
    """Round an exact Fraction to nearest dyadic k / 2^bits using pure integer arithmetic."""
    num = q.numerator << (bits + 1)
    den = q.denominator
    if num >= 0:
        k = (num + den) // (2 * den)
    else:
        k = -((-num + den) // (2 * den))
    return Fraction(k, 1 << bits)


# ---------------------------------------------------------------------------
# Pure-Z/Q Point-Charge Configuration Certifier (Krawczyk + Bareiss Inertia)
# ---------------------------------------------------------------------------


def certify_configuration(
    name: str,
    sites_q: list[tuple[Fraction, Fraction, Fraction]],
    charges_q: list[Fraction],
    centers_q: list[tuple[Fraction, Fraction, Fraction]],
    hws_q: list[Fraction],
) -> dict[str, Any]:
    """Certify a finite set of pairwise-disjoint non-degenerate equilibria in pure Z/Q."""
    cfg_name = require_str(name, "name", min_len=1, max_len=128)
    s_seq = require_seq(sites_q, "sites_q", min_len=2, max_len=MAX_CHARGES)
    q_seq = require_seq(charges_q, "charges_q", min_len=2, max_len=MAX_CHARGES)
    if len(s_seq) != len(q_seq):
        raise PointChargeFrontierError("SITE_CHARGE_LEN_MISMATCH", "sites_q and charges_q lengths differ")
    c_seq = require_seq(centers_q, "centers_q", min_len=1, max_len=MAX_BOXES)
    h_seq = require_seq(hws_q, "hws_q", min_len=1, max_len=MAX_BOXES)
    if len(c_seq) != len(h_seq):
        raise PointChargeFrontierError("CENTER_HW_LEN_MISMATCH", "centers_q and hws_q lengths differ")

    parsed_sites: list[tuple[Fraction, Fraction, Fraction]] = []
    for idx, site in enumerate(s_seq):
        coords = require_seq(site, f"site[{idx}]", min_len=3, max_len=3)
        for d in range(3):
            if type(coords[d]) is not Fraction:
                raise PointChargeFrontierError("INVALID_FRACTION_TYPE", "site coordinates must be strict Fraction")
        parsed_sites.append((coords[0], coords[1], coords[2]))

    for i in range(len(parsed_sites)):
        for j in range(i + 1, len(parsed_sites)):
            if parsed_sites[i] == parsed_sites[j]:
                raise PointChargeFrontierError("DUPLICATE_CHARGE_SITE", f"charges {i} and {j} share identical coordinates")

    parsed_charges: list[Fraction] = []
    for idx, q in enumerate(q_seq):
        if type(q) is not Fraction:
            raise PointChargeFrontierError("INVALID_FRACTION_TYPE", f"charge[{idx}] must be strict Fraction")
        if q <= 0:
            raise PointChargeFrontierError("NON_POSITIVE_CHARGE", f"charge[{idx}] must be strictly positive")
        parsed_charges.append(q)

    parsed_centers: list[tuple[Fraction, Fraction, Fraction]] = []
    parsed_hws: list[Fraction] = []
    boxes_q: list[tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction], tuple[Fraction, Fraction]]] = []
    for idx, (ctr, hw) in enumerate(zip(c_seq, h_seq)):
        coords = require_seq(ctr, f"center[{idx}]", min_len=3, max_len=3)
        for d in range(3):
            if type(coords[d]) is not Fraction:
                raise PointChargeFrontierError("INVALID_FRACTION_TYPE", "center coordinates must be strict Fraction")
        if type(hw) is not Fraction:
            raise PointChargeFrontierError("INVALID_FRACTION_TYPE", f"hw[{idx}] must be strict Fraction")
        if hw <= 0:
            raise PointChargeFrontierError("NON_POSITIVE_HALFWIDTH", f"hw[{idx}] must be strictly positive")
        parsed_centers.append((coords[0], coords[1], coords[2]))
        parsed_hws.append(hw)
        boxes_q.append(((coords[0] - hw, coords[0] + hw), (coords[1] - hw, coords[1] + hw), (coords[2] - hw, coords[2] + hw)))

    # Exact pairwise strict disjointness of closed boxes over Q
    n_boxes = len(boxes_q)
    for i in range(n_boxes):
        for j in range(i + 1, n_boxes):
            separated = any(
                boxes_q[i][d][1] < boxes_q[j][d][0] or boxes_q[j][d][1] < boxes_q[i][d][0]
                for d in range(3)
            )
            if not separated:
                raise PointChargeFrontierError("OVERLAPPING_BOXES", f"boxes {i} and {j} are not strictly disjoint")

    # Verify no charge lies inside or on the boundary of any box
    for b_idx, box in enumerate(boxes_q):
        for s_idx, site in enumerate(parsed_sites):
            if all(box[d][0] <= site[d] <= box[d][1] for d in range(3)):
                raise PointChargeFrontierError("CHARGE_INSIDE_BOX", f"charge {s_idx} lies inside box {b_idx}")

    sites_iv = [[iv_from_q(c) for c in s] for s in parsed_sites]
    charges_iv = [iv_from_q(q) for q in parsed_charges]

    max_rho = Fraction(0, 1)
    morse_counts = [0, 0, 0, 0]
    certified_boxes: list[dict[str, Any]] = []

    for idx, (c_q, hw_q) in enumerate(zip(parsed_centers, parsed_hws)):
        c_iv = [iv_from_q(c_q[d]) for d in range(3)]
        hw_iv = iv_from_q(hw_q)
        box_iv = [(c_iv[d][0] - hw_iv[1], c_iv[d][1] + hw_iv[1]) for d in range(3)]

        g_c, h_c = eval_grad_hess_iv(c_iv, sites_iv, charges_iv)
        _g_b, h_b = eval_grad_hess_iv(box_iv, sites_iv, charges_iv)

        # Construct exact symmetric rational midpoint Hessian H_tilde from H(c)
        h_mid_q = [
            [Fraction((h_c[i][j][0] + h_c[i][j][1]) // 2, ONE_SCALED) for j in range(3)]
            for i in range(3)
        ]
        for i in range(3):
            for j in range(i + 1, 3):
                avg = (h_mid_q[i][j] + h_mid_q[j][i]) / 2
                h_mid_q[i][j] = avg
                h_mid_q[j][i] = avg

        # Verify H_tilde lies inside interval enclosure H(box) entrywise
        for i in range(3):
            for j in range(3):
                v_iv = iv_from_q(h_mid_q[i][j])
                if v_iv[0] < h_b[i][j][0] or v_iv[1] > h_b[i][j][1]:
                    raise PointChargeFrontierError("MIDPOINT_OUTSIDE_BOX_HESS", f"box {idx} midpoint H outside H(box)")

        # Exact rational preconditioner Y = round_160(H_tilde^-1); verify exact det(Y) != 0
        y_exact = exact_inv_3x3_q(h_mid_q)
        y_pt = [[_round_q_dyadic(y_exact[i][j], Y_ROUND_BITS) for j in range(3)] for i in range(3)]
        det_y = (
            y_pt[0][0] * (y_pt[1][1] * y_pt[2][2] - y_pt[1][2] * y_pt[2][1])
            - y_pt[0][1] * (y_pt[1][0] * y_pt[2][2] - y_pt[1][2] * y_pt[2][0])
            + y_pt[0][2] * (y_pt[1][0] * y_pt[2][1] - y_pt[1][1] * y_pt[2][0])
        )
        if det_y == 0:
            raise PointChargeFrontierError("SINGULAR_PRECONDITIONER", f"box {idx} preconditioner Y is singular")
        y_iv = [[iv_from_q(y_pt[i][j]) for j in range(3)] for i in range(3)]

        # Krawczyk contraction matrix M = I - Y * H(box)
        m_iv = [[(0, 0) for _ in range(3)] for _ in range(3)]
        for i in range(3):
            for j in range(3):
                acc = (0, 0)
                for k in range(3):
                    acc = iv_add(acc, iv_mul(y_iv[i][k], h_b[k][j]))
                m_iv[i][j] = iv_sub((ONE_SCALED, ONE_SCALED), acc) if i == j else iv_neg(acc)

        rho_int = max(sum(iv_abs_max(m_iv[i][j]) for j in range(3)) for i in range(3))
        if rho_int >= ONE_SCALED:
            raise PointChargeFrontierError("KRAWCZYK_NOT_CONTRACTIVE", f"box {idx} contraction norm >= 1")
        rho_q = Fraction(rho_int, ONE_SCALED)
        if rho_q > max_rho:
            max_rho = rho_q

        # Strict interior inclusion: K(box) - c = -Y*g(c) + M*(box - c) strictly inside (-hw, +hw)
        dx_box = [(-hw_iv[1], hw_iv[1]) for _ in range(3)]
        for i in range(3):
            yg = (0, 0)
            mdx = (0, 0)
            for k in range(3):
                yg = iv_add(yg, iv_mul(y_iv[i][k], g_c[k]))
                mdx = iv_add(mdx, iv_mul(m_iv[i][k], dx_box[k]))
            shift = iv_sub(mdx, yg)
            if shift[0] <= -hw_iv[0] or shift[1] >= hw_iv[0]:
                raise PointChargeFrontierError("KRAWCZYK_NO_INCLUSION", f"box {idx} fails strict K(X) subset int(X)")

        # Harmonic consistency: Laplace's equation tr(H) = 0 -> 0 in tr(H(box))
        tr_iv = iv_add(iv_add(h_b[0][0], h_b[1][1]), h_b[2][2])
        if tr_iv[0] > 0 or tr_iv[1] < 0:
            raise PointChargeFrontierError("HARMONIC_TRACE_VIOLATION", f"box {idx} trace interval excludes 0")

        # Non-degeneracy across the ENTIRE box: 0 not in det(H(box))
        det_iv = det_3x3_iv(h_b)
        if det_iv[0] <= 0 and det_iv[1] >= 0:
            raise PointChargeFrontierError("DEGENERATE_BOX_HESSIAN", f"box {idx} det(H(X)) contains 0")

        # Constant Morse index across connected box from exact self-verified Lagrange/Bareiss inertia at H_tilde
        _n_pos, n_neg, n_zero = exact_lagrange_inertia_3x3_q(h_mid_q)
        if n_zero != 0 or n_neg not in (1, 2):
            raise PointChargeFrontierError("INVALID_HARMONIC_MORSE_INDEX", f"box {idx} non-harmonic inertia ({_n_pos},{n_neg},{n_zero})")
        morse_counts[n_neg] += 1
        certified_boxes.append(
            {
                "center": [fmt_q(c_q[0]), fmt_q(c_q[1]), fmt_q(c_q[2])],
                "half_width": fmt_q(hw_q),
                "morse_index": n_neg,
                "rho_upper_q": fmt_q(rho_q),
            }
        )

    n_charges = len(parsed_charges)
    maxwell_bound = (n_charges - 1) * (n_charges - 1)
    m2_minus_m1 = morse_counts[2] - morse_counts[1]
    return {
        "name": cfg_name,
        "n_charges": n_charges,
        "maxwell_bound": maxwell_bound,
        "n_certified": n_boxes,
        "exceeds_maxwell_bound": n_boxes > maxwell_bound,
        "morse_counts": morse_counts,
        "m2_minus_m1": m2_minus_m1,
        "morse_relation_matches_n_minus_1": m2_minus_m1 == n_charges - 1,
        "max_rho_q": fmt_q(max_rho),
        "sites": [[fmt_q(c) for c in s] for s in parsed_sites],
        "charges": [fmt_q(q) for q in parsed_charges],
        "boxes": certified_boxes,
    }


CANONICAL_WITNESS_CONFIGS: tuple[dict[str, Any], ...] = (
    {
        "name": 'n2-analytic-control',
        "description": 'Analytic Maxwell control (n=2 equal positive charges at (+/-1, 0, 0), 1 index-2 saddle at origin)',
        "sites": (
            ('-1', '0', '0'),
            ('1', '0', '0'),
        ),
        "charges": (
            '1',
            '1',
        ),
        "half_width": '1/1048576',
        "centers": (
            ('0', '0', '0'),
        ),
    },
    {
        "name": 'n7-abk-second-pair',
        "description": 'Foundation-0/point-charges#12 ABK C_7 second-pair configuration (n=7, 44 equilibria > Maxwell 36)',
        "sites": (
            ('3', '0', '0'),
            ('0', '3', '0'),
            ('0', '0', '3'),
            ('5/4', '5/4', '5/4'),
            ('3/4', '3/4', '3/4'),
            ('2049/2048', '2049/2048', '2049/2048'),
            ('2047/2048', '2047/2048', '2047/2048'),
        ),
        "charges": (
            '1',
            '1',
            '1',
            '541/131072',
            '541/131072',
            '134713/1152921504606846976',
            '134713/1152921504606846976',
        ),
        "half_width": '1/1125899906842624',
        "centers": (
            ('203942654963791710082779423135747057/332306998946228968225951765070086144', '792978341874895194595075872074511375/664613997892457936451903530140172288', '792978341874895194595075872074511375/664613997892457936451903530140172288'),
            ('1005809566130455897718319325846236285/1329227995784915872903807060280344576', '1490937210612145860496550927497398721/1329227995784915872903807060280344576', '1490937210612145860496550927497398721/1329227995784915872903807060280344576'),
            ('1280515948638237747237787573706536669/1329227995784915872903807060280344576', '1280515948638237747237787573706536669/1329227995784915872903807060280344576', '1330000970256274191672476225103083377/1329227995784915872903807060280344576'),
            ('1280515948638237747237787573706536669/1329227995784915872903807060280344576', '1330000970256274191672476225103083377/1329227995784915872903807060280344576', '1280515948638237747237787573706536669/1329227995784915872903807060280344576'),
            ('81917151653440547298124337141130283/83076749736557242056487941267521536', '81917151653440547298124337141130283/83076749736557242056487941267521536', '1320287505071631389176787894018243679/1329227995784915872903807060280344576'),
            ('81917151653440547298124337141130283/83076749736557242056487941267521536', '1320287505071631389176787894018243679/1329227995784915872903807060280344576', '81917151653440547298124337141130283/83076749736557242056487941267521536'),
            ('1314921940749052881303715033122472731/1329227995784915872903807060280344576', '1314921940749052881303715033122472731/1329227995784915872903807060280344576', '1314921940749052881303715033122472731/1329227995784915872903807060280344576'),
            ('1320287505071631389176787894018243679/1329227995784915872903807060280344576', '81917151653440547298124337141130283/83076749736557242056487941267521536', '81917151653440547298124337141130283/83076749736557242056487941267521536'),
            ('165825530496870175096033424998257389/166153499473114484112975882535043072', '332634967922473277242894222606871827/332306998946228968225951765070086144', '332634967922473277242894222606871827/332306998946228968225951765070086144'),
            ('664503463926600025819274146878633747/664613997892457936451903530140172288', '1329338529750773783536436443541883117/1329227995784915872903807060280344576', '1329338529750773783536436443541883117/1329227995784915872903807060280344576'),
            ('1329178084983038283057069741980852385/1329227995784915872903807060280344576', '1329178084983038283057069741980852385/1329227995784915872903807060280344576', '664614189916288998844188789374004897/664613997892457936451903530140172288'),
            ('1329178084983038283057069741980852385/1329227995784915872903807060280344576', '664614189916288998844188789374004897/664613997892457936451903530140172288', '1329178084983038283057069741980852385/1329227995784915872903807060280344576'),
            ('664603397504902796118706775852019345/664613997892457936451903530140172288', '664603397504902796118706775852019345/664613997892457936451903530140172288', '332304662314466247295452358671533593/332306998946228968225951765070086144'),
            ('664603397504902796118706775852019345/664613997892457936451903530140172288', '332304662314466247295452358671533593/332306998946228968225951765070086144', '664603397504902796118706775852019345/664613997892457936451903530140172288'),
            ('1329212164427590736663192136782818095/1329227995784915872903807060280344576', '1329212164427590736663192136782818095/1329227995784915872903807060280344576', '1329212164427590736663192136782818095/1329227995784915872903807060280344576'),
            ('332304662314466247295452358671533593/332306998946228968225951765070086144', '664603397504902796118706775852019345/664613997892457936451903530140172288', '664603397504902796118706775852019345/664613997892457936451903530140172288'),
            ('166153093259284307012991444016122673/166153499473114484112975882535043072', '332307405160059145325936203589006543/332306998946228968225951765070086144', '332307405160059145325936203589006543/332306998946228968225951765070086144'),
            ('1', '1', '1'),
            ('664614189916288998844188789374004897/664613997892457936451903530140172288', '1329178084983038283057069741980852385/1329227995784915872903807060280344576', '1329178084983038283057069741980852385/1329227995784915872903807060280344576'),
            ('332307405160059145325936203589006543/332306998946228968225951765070086144', '166153093259284307012991444016122673/166153499473114484112975882535043072', '332307405160059145325936203589006543/332306998946228968225951765070086144'),
            ('332307405160059145325936203589006543/332306998946228968225951765070086144', '332307405160059145325936203589006543/332306998946228968225951765070086144', '166153093259284307012991444016122673/166153499473114484112975882535043072'),
            ('1329241293727986555607269980201920007/1329227995784915872903807060280344576', '1329241293727986555607269980201920007/1329227995784915872903807060280344576', '1329253147976045952551665863184015689/1329227995784915872903807060280344576'),
            ('1329241293727986555607269980201920007/1329227995784915872903807060280344576', '1329253147976045952551665863184015689/1329227995784915872903807060280344576', '1329241293727986555607269980201920007/1329227995784915872903807060280344576'),
            ('1329243827142241009144421983777871057/1329227995784915872903807060280344576', '1329243827142241009144421983777871057/1329227995784915872903807060280344576', '1329243827142241009144421983777871057/1329227995784915872903807060280344576'),
            ('1329244376687100319663005820735065161/1329227995784915872903807060280344576', '1329244376687100319663005820735065161/1329227995784915872903807060280344576', '664647335768320017147156828751111285/664613997892457936451903530140172288'),
            ('1329244376687100319663005820735065161/1329227995784915872903807060280344576', '664647335768320017147156828751111285/664613997892457936451903530140172288', '1329244376687100319663005820735065161/1329227995784915872903807060280344576'),
            ('1329253147976045952551665863184015689/1329227995784915872903807060280344576', '1329241293727986555607269980201920007/1329227995784915872903807060280344576', '1329241293727986555607269980201920007/1329227995784915872903807060280344576'),
            ('664647335768320017147156828751111285/664613997892457936451903530140172288', '1329244376687100319663005820735065161/1329227995784915872903807060280344576', '1329244376687100319663005820735065161/1329227995784915872903807060280344576'),
            ('1329338529750773783536436443541883117/1329227995784915872903807060280344576', '664503463926600025819274146878633747/664613997892457936451903530140172288', '1329338529750773783536436443541883117/1329227995784915872903807060280344576'),
            ('1329338529750773783536436443541883117/1329227995784915872903807060280344576', '1329338529750773783536436443541883117/1329227995784915872903807060280344576', '664503463926600025819274146878633747/664613997892457936451903530140172288'),
            ('1330000970256274191672476225103083377/1329227995784915872903807060280344576', '1280515948638237747237787573706536669/1329227995784915872903807060280344576', '1280515948638237747237787573706536669/1329227995784915872903807060280344576'),
            ('332634967922473277242894222606871827/332306998946228968225951765070086144', '165825530496870175096033424998257389/166153499473114484112975882535043072', '332634967922473277242894222606871827/332306998946228968225951765070086144'),
            ('332634967922473277242894222606871827/332306998946228968225951765070086144', '332634967922473277242894222606871827/332306998946228968225951765070086144', '165825530496870175096033424998257389/166153499473114484112975882535043072'),
            ('1341372846037061234099759059795831857/1329227995784915872903807060280344576', '1341372846037061234099759059795831857/1329227995784915872903807060280344576', '1350985924653643866506557559555991007/1329227995784915872903807060280344576'),
            ('1341372846037061234099759059795831857/1329227995784915872903807060280344576', '1350985924653643866506557559555991007/1329227995784915872903807060280344576', '1341372846037061234099759059795831857/1329227995784915872903807060280344576'),
            ('1343534050820778864503899087438216421/1329227995784915872903807060280344576', '1343534050820778864503899087438216421/1329227995784915872903807060280344576', '1343534050820778864503899087438216421/1329227995784915872903807060280344576'),
            ('1344950028519569702280034112589788011/1329227995784915872903807060280344576', '1344950028519569702280034112589788011/1329227995784915872903807060280344576', '1394435050137606146714722763986334719/1329227995784915872903807060280344576'),
            ('1344950028519569702280034112589788011/1329227995784915872903807060280344576', '1394435050137606146714722763986334719/1329227995784915872903807060280344576', '1344950028519569702280034112589788011/1329227995784915872903807060280344576'),
            ('1350985924653643866506557559555991007/1329227995784915872903807060280344576', '1341372846037061234099759059795831857/1329227995784915872903807060280344576', '1341372846037061234099759059795831857/1329227995784915872903807060280344576'),
            ('1394435050137606146714722763986334719/1329227995784915872903807060280344576', '1344950028519569702280034112589788011/1329227995784915872903807060280344576', '1344950028519569702280034112589788011/1329227995784915872903807060280344576'),
            ('1490937210612145860496550927497398721/1329227995784915872903807060280344576', '1005809566130455897718319325846236285/1329227995784915872903807060280344576', '1490937210612145860496550927497398721/1329227995784915872903807060280344576'),
            ('1490937210612145860496550927497398721/1329227995784915872903807060280344576', '1490937210612145860496550927497398721/1329227995784915872903807060280344576', '1005809566130455897718319325846236285/1329227995784915872903807060280344576'),
            ('792978341874895194595075872074511375/664613997892457936451903530140172288', '203942654963791710082779423135747057/332306998946228968225951765070086144', '792978341874895194595075872074511375/664613997892457936451903530140172288'),
            ('792978341874895194595075872074511375/664613997892457936451903530140172288', '792978341874895194595075872074511375/664613997892457936451903530140172288', '203942654963791710082779423135747057/332306998946228968225951765070086144'),
        ),
    },
    {
        "name": 'n6-abk-milnor8-flat',
        "description": 'Pure-Z/Q certified n=6 ABK + axial Milnor-8 flat-point unfolding (n=6, 31 equilibria > 29 > Maxwell 25)',
        "sites": (
            ('27670125098345552411/9223372036854775808', '-2522119748711/4611686018427387904', '-2133464754377/2305843009213693952'),
            ('-5992654548461/18446744073709551616', '27670130190693884371/9223372036854775808', '-7514751936875/18446744073709551616'),
            ('-939293349451/1152921504606846976', '-4785775969729/9223372036854775808', '55340250551384529109/18446744073709551616'),
            ('5/4', '5/4', '5/4'),
            ('3/4', '3/4', '3/4'),
            ('2317584721377142147/2305843009213693952', '9270338886456781389/9223372036854775808', '9270338885624060141/9223372036854775808'),
        ),
        "charges": (
            '18446722617876621735/18446744073709551616',
            '9223361054979265767/9223372036854775808',
            '4611680690257545621/4611686018427387904',
            '541/131072',
            '541/131072',
            '707809575/9223372036854775808',
        ),
        "half_width": '1/36028797018963968',
        "centers": (
            ('815771732644403750048935843116951191/1329227995784915872903807060280344576', '1585956311663279801373968849148576059/1329227995784915872903807060280344576', '1585956072970681429564244097228800309/1329227995784915872903807060280344576'),
            ('502902904215249733676502410379164797/664613997892457936451903530140172288', '1490939262906614473216844426600598065/1329227995784915872903807060280344576', '1490939027962539923887679927461344177/1329227995784915872903807060280344576'),
            ('640262102073533066460006508611661735/664613997892457936451903530140172288', '1330000401781508907269709782695784257/1329227995784915872903807060280344576', '1280524116396781772058931189741922205/1329227995784915872903807060280344576'),
            ('640260234551509340384667055912193681/664613997892457936451903530140172288', '640260169985312516084222145723752681/664613997892457936451903530140172288', '1330001274683794089878242501088257735/1329227995784915872903807060280344576'),
            ('327643187964398882089777962046651195/332306998946228968225951765070086144', '1279855576133557132630969614474543/1298074214633706907132624082305024', '660158207769244416541812165370061701/664613997892457936451903530140172288'),
            ('81910534961110156940823029589365369/83076749736557242056487941267521536', '330079836664051240593534864222262795/332306998946228968225951765070086144', '40955262781785109047777304724647389/41538374868278621028243970633760768'),
            ('1314879736782089198117537174438382701/1329227995784915872903807060280344576', '1314876971291922125004260036031197997/1329227995784915872903807060280344576', '1314879397562936376477549888496818969/1329227995784915872903807060280344576'),
            ('330079063051270511149815885054498223/332306998946228968225951765070086144', '655286221837833111100066539241785155/664613997892457936451903530140172288', '655286467601912124032742056216938839/664613997892457936451903530140172288'),
            ('166425184886715012354903928330394955/166153499473114484112975882535043072', '83213143134553200390258721713278771/83076749736557242056487941267521536', '1331409960790397760339753441128347639/1329227995784915872903807060280344576'),
            ('665000664030377934428967529933617735/664613997892457936451903530140172288', '1280520303504074986062708326874255877/1329227995784915872903807060280344576', '160065043111123101287499005146443973/166153499473114484112975882535043072'),
            ('166426248684025442434749312556812213/166153499473114484112975882535043072', '1331410195011375475858390163834041157/1329227995784915872903807060280344576', '1331401571951839178601451254345687981/1329227995784915872903807060280344576'),
            ('665704973787880014353404989618891139/664613997892457936451903530140172288', '1331401466837890952383746470915214131/1329227995784915872903807060280344576', '665704914457094830227562950269820939/664613997892457936451903530140172288'),
            ('1341416537751933526794802368131435217/1329227995784915872903807060280344576', '1341417300779004498191111287383274679/1329227995784915872903807060280344576', '675982491639687936618487407727819241/664613997892457936451903530140172288'),
            ('1341416989583880552561284925906460271/1329227995784915872903807060280344576', '42248835123432252146655069210679695/41538374868278621028243970633760768', '1341417010785816974516695325662125679/1329227995784915872903807060280344576'),
            ('671937924800947670673656546076523157/664613997892457936451903530140172288', '671938836239214121239222044910939905/664613997892457936451903530140172288', '1343876071026198545366987047237026457/1329227995784915872903807060280344576'),
            ('42029575540586143742033039472073161/41538374868278621028243970633760768', '1344946866309826577656512540972655203/1329227995784915872903807060280344576', '697205964862729083189637807889046825/664613997892457936451903530140172288'),
            ('1344946703722200914478136788143034409/1329227995784915872903807060280344576', '697206763466343291709279742525839433/664613997892457936451903530140172288', '1344946686579566097512538537921795467/1329227995784915872903807060280344576'),
            ('675982869111518777459908206050023233/664613997892457936451903530140172288', '1341417215348934027333264026576423525/1329227995784915872903807060280344576', '335354118381658305183618362044220429/332306998946228968225951765070086144'),
            ('348602773602332330299843131708155405/332306998946228968225951765070086144', '1344946765355529802935345445498976503/1329227995784915872903807060280344576', '42029571849989011144353046217069329/41538374868278621028243970633760768'),
            ('745470087622223115770495666568587939/664613997892457936451903530140172288', '502901714858052547073470183223265565/664613997892457936451903530140172288', '745470062419063127948648954309182433/664613997892457936451903530140172288'),
            ('186367339954807465740534976520818945/166153499473114484112975882535043072', '1490938904176221366912686143905886207/1329227995784915872903807060280344576', '502903215179954680435336821725735763/664613997892457936451903530140172288'),
            ('792978196986372249179871090960012385/664613997892457936451903530140172288', '792978290600183178041040309379122501/664613997892457936451903530140172288', '407885546549929886189813994240499231/664613997892457936451903530140172288'),
            ('792977487074557282188782798982800341/664613997892457936451903530140172288', '815773819169910814552723980188168253/1329227995784915872903807060280344576', '1585954922684139728185954433103564099/1329227995784915872903807060280344576'),
            ('665705303475727493709726595431577613/664613997892457936451903530140172288', '1331416206555768874866507427052985549/1329227995784915872903807060280344576', '1331410725773359341377266535276145649/1329227995784915872903807060280344576'),
            ('665711449879755199440685934999436493/664613997892457936451903530140172288', '1331410863766211304415551469216743385/1329227995784915872903807060280344576', '1331410454760143572772635808185763811/1329227995784915872903807060280344576'),
            ('665705541344758583815435738601748687/664613997892457936451903530140172288', '1331412637145724791962205880252154547/1329227995784915872903807060280344576', '166426437334723418369606070015919927/166153499473114484112975882535043072'),
            ('665704625548425860184195641694567525/664613997892457936451903530140172288', '1331408825377256014276944330946930851/1329227995784915872903807060280344576', '1331408787979960824953629172665247297/1329227995784915872903807060280344576'),
            ('1331410486254337914535666126733812139/1329227995784915872903807060280344576', '665705490296432027099438408794897503/664613997892457936451903530140172288', '665709317345847648802388234552816423/664613997892457936451903530140172288'),
            ('666008653958002831765182258661192753/664613997892457936451903530140172288', '665711147513603068195752463652214641/664613997892457936451903530140172288', '1331421923915589434429440487715246947/1329227995784915872903807060280344576'),
            ('1331422149248751945968584303228515711/1329227995784915872903807060280344576', '666011921378410560995671935204403517/664613997892457936451903530140172288', '665711060234610818379156451390193805/664613997892457936451903530140172288'),
            ('166427758105702534375596919334535853/166153499473114484112975882535043072', '166427800892238113907403726349528129/166153499473114484112975882535043072', '1332021261310304870746585337728207377/1329227995784915872903807060280344576'),
        ),
    },
    {
        "name": 'n8-abk-milnor8-flat',
        "description": 'Foundation-0/point-charges#5 pure-Z/Q certified n=8 Maxwell counterexample via ABK + axial Milnor-8 flat-point unfolding (n=8, 51 equilibria > Maxwell 49)',
        "sites": (
            ('118842243817740693672398590383/39614081257132168796771975168', '-25660458062846023387/39614081257132168796771975168', '-87100269652471997859/79228162514264337593543950336'),
            ('-30230001662827563505/79228162514264337593543950336', '59421121921901776177281571331/19807040628566084398385987584', '-9508628234851205439/19807040628566084398385987584'),
            ('-76642335266324740097/79228162514264337593543950336', '-48667493016122047263/79228162514264337593543950336', '237684487637281217871123174865/79228162514264337593543950336'),
            ('5/4', '5/4', '5/4'),
            ('3/4', '3/4', '3/4'),
            ('2049/2048', '2049/2048', '2049/2048'),
            ('2047/2048', '2047/2048', '2047/2048'),
            ('39614322060840975568632608655/39614081257132168796771975168', '79228644121681951150891990137/79228162514264337593543950336', '79228644121681951138926868911/79228162514264337593543950336'),
        ),
        "charges": (
            '79228162404230555144801193541/79228162514264337593543950336',
            '79228162401607936164828517669/79228162514264337593543950336',
            '79228162404962230859392497281/79228162514264337593543950336',
            '541/131072',
            '541/131072',
            '134713/1152921504606846976',
            '134713/1152921504606846976',
            '6134079/39614081257132168796771975168',
        ),
        "half_width": '1/37778931862957161709568',
        "centers": (
            ('875927034731714744776631457678522826829707865/1427247692705959881058285969449495136382746624', '1702908021941208294159084085784428592316367005/1427247692705959881058285969449495136382746624', '425727005409353943285192713635024669593491357/356811923176489970264571492362373784095686656'),
            ('269994948361937232242178957776849353264345159/356811923176489970264571492362373784095686656', '800440821286364883305235247975353659502757845/713623846352979940529142984724747568191373312', '400220410568652039891029455392272668043283483/356811923176489970264571492362373784095686656'),
            ('343735880517015499646176548860919996716912037/356811923176489970264571492362373784095686656', '687471760963658460921784717140701721570727825/713623846352979940529142984724747568191373312', '714038836406108461523338582262963713811651247/713623846352979940529142984724747568191373312'),
            ('1374943526736514125214512810002937377439230079/1427247692705959881058285969449495136382746624', '1428077671713905742226577290134870804036042837/1427247692705959881058285969449495136382746624', '1374943526627115833731554136490546965184918617/1427247692705959881058285969449495136382746624'),
            ('175915745512280728993171378922635592002385033/178405961588244985132285746181186892047843328', '703662981698045007850272359131317352331586193/713623846352979940529142984724747568191373312', '1417647914032504934840891748012335513863158127/1427247692705959881058285969449495136382746624'),
            ('351831489792032987732094438997948791242491357/356811923176489970264571492362373784095686656', '1417647917469354357849754970478190469145478957/1427247692705959881058285969449495136382746624', '1407325958990153314517090146043974689290335339/1427247692705959881058285969449495136382746624'),
            ('1411886691300446855150112245711601776527697099/1427247692705959881058285969449495136382746624', '705943344092909741450643422984003218873112633/713623846352979940529142984724747568191373312', '1411886690918771249265232029140980507143382935/1427247692705959881058285969449495136382746624'),
            ('1417647913868379908987446141542399012252379053/1427247692705959881058285969449495136382746624', '1407325963762045649029564085169633754936216735/1427247692705959881058285969449495136382746624', '1407325964286222818593230545011593315900282125/1427247692705959881058285969449495136382746624'),
            ('356107615896579222531620022637291452663908741/356811923176489970264571492362373784095686656', '1428656307563555935348885953988005499191445775/1427247692705959881058285969449495136382746624', '1428656307270468040883267052977987173543928081/1427247692705959881058285969449495136382746624'),
            ('1427010319717041802289161063720511812473101743/1427247692705959881058285969449495136382746624', '1427366379412245360188540934728755691151028397/1427247692705959881058285969449495136382746624', '1427366379115484268158198079950017116532366755/1427247692705959881058285969449495136382746624'),
            ('1427194115693118822319865874957082754745816659/1427247692705959881058285969449495136382746624', '1427194115504615776920828034924911220226931511/1427247692705959881058285969449495136382746624', '1427248102342747988664840008768233772010622191/1427247692705959881058285969449495136382746624'),
            ('1427194120947130210231194075686428588328532267/1427247692705959881058285969449495136382746624', '713624050499684902045477851787181943613166035/713623846352979940529142984724747568191373312', '713597060416198377869910593296324093490778289/713623846352979940529142984724747568191373312'),
            ('1427224767691872325842241592665490675825941459/1427247692705959881058285969449495136382746624', '1427224766949702961339169906122114363812344927/1427247692705959881058285969449495136382746624', '1427237706965453141161050135144920234299289585/1427247692705959881058285969449495136382746624'),
            ('356806190492144652190076841306929196409100555/356811923176489970264571492362373784095686656', '1427237710703285831487817182154048896556972207/1427247692705959881058285969449495136382746624', '713612380893221736234181547937005394333027817/713623846352979940529142984724747568191373312'),
            ('1427230632253754840087101461673318781862925847/1427247692705959881058285969449495136382746624', '1427230628921759161416389515798960110227605037/1427247692705959881058285969449495136382746624', '1427230631847410179610041206821922702307134161/1427247692705959881058285969449495136382746624'),
            ('713618853388862094363539938779496857151375181/713623846352979940529142984724747568191373312', '1427224767368571324372343092818924960152092567/1427247692705959881058285969449495136382746624', '713612383964296771190424147462944446285849139/713623846352979940529142984724747568191373312'),
            ('356812615789906268208419838010186881490114195/356811923176489970264571492362373784095686656', '1427250474943183961348305468190051078596215015/1427247692705959881058285969449495136382746624', '178406309314788114209276639470889250943766921/178405961588244985132285746181186892047843328'),
            ('713624051206080003095571689936451345468501603/713623846352979940529142984724747568191373312', '89199632216351175104752865159968980430515421/89202980794122492566142873090593446023921664', '1427194115535393017348897697559325840089318619/1427247692705959881058285969449495136382746624'),
            ('1427250474506367071238832153935042328364297937/1427247692705959881058285969449495136382746624', '356812615735069230758867283407138277423423393/356811923176489970264571492362373784095686656', '1427250474360867750639615796248260260353373215/1427247692705959881058285969449495136382746624'),
            ('713625237277969092946608695094110548007454397/713623846352979940529142984724747568191373312', '1427250474829249731648844247860216095711404995/1427247692705959881058285969449495136382746624', '1427250463158736080128440710054939737766859035/1427247692705959881058285969449495136382746624'),
            ('356815507093036489516897587756053289349569711/356811923176489970264571492362373784095686656', '713631014636334387459930132722628017127779783/713623846352979940529142984724747568191373312', '1427276295980260295749613164486196138846803649/1427247692705959881058285969449495136382746624'),
            ('713631014414383587854668696156552146556538155/713623846352979940529142984724747568191373312', '713638146604940899353425155337391655416950617/713623846352979940529142984724747568191373312', '1427262028846901803972514610365862310725255189/1427247692705959881058285969449495136382746624'),
            ('44602039129181518391512567951655936727850035/44601490397061246283071436545296723011960832', '713632627089214347373381180281530911930654065/713623846352979940529142984724747568191373312', '1427265252383278237598706222853633786107243633/1427247692705959881058285969449495136382746624'),
            ('178408159313785642983199541623816484832715661/178405961588244985132285746181186892047843328', '713632637548071165671819344997787365498997169/713623846352979940529142984724747568191373312', '1427319224089277870490048661805296979318804713/1427247692705959881058285969449495136382746624'),
            ('356816318716831584730267548703972241503154073/356811923176489970264571492362373784095686656', '1427319226193321711310998391678403455196703969/1427247692705959881058285969449495136382746624', '356816318711755681456672096806842325211389325/356811923176489970264571492362373784095686656'),
            ('713638148496828130345514468459494884126136781/713623846352979940529142984724747568191373312', '713631014594322402613975130814521582972870875/713623846352979940529142984724747568191373312', '713631014153130162551423333558018226915651809/713623846352979940529142984724747568191373312'),
            ('356829805727229096818615592266747036128394403/356811923176489970264571492362373784095686656', '713632637482557471198883727167931024308269515/713623846352979940529142984724747568191373312', '1427265274358950708859621649403664239262709489/1427247692705959881058285969449495136382746624'),
            ('713683189772243076780712251016305176484367961/713623846352979940529142984724747568191373312', '1427010318759833004245621248936659966555288027/1427247692705959881058285969449495136382746624', '1427366379481189182055859312223187287091488697/1427247692705959881058285969449495136382746624'),
            ('178420797376625957311233148322707869168697691/178405961588244985132285746181186892047843328', '1427366379246471790667298370272490541999515363/1427247692705959881058285969449495136382746624', '178376289991161133488154697577375526981812359/178405961588244985132285746181186892047843328'),
            ('1428077672877752210921674108265341308730907803/1427247692705959881058285969449495136382746624', '1374943521897158146906220380287908696906068381/1427247692705959881058285969449495136382746624', '343735880482126233992881926746175823028057983/356811923176489970264571492362373784095686656'),
            ('1428656306781449610531897243595848515557977907/1427247692705959881058285969449495136382746624', '1424430464092649884788682323469553125777144729/1427247692705959881058285969449495136382746624', '1428656306717639856402213670121712876607202371/1427247692705959881058285969449495136382746624'),
            ('714328153689056952385109105550899444195466271/713623846352979940529142984724747568191373312', '1428656307607392045100819114106195347587577829/1427247692705959881058285969449495136382746624', '1424430463333794748447846002722957716886978215/1427247692705959881058285969449495136382746624'),
            ('1440288120609417504418506547172225645180068371/1427247692705959881058285969449495136382746624', '1440288121771981983293128037793730890798633573/1427247692705959881058285969449495136382746624', '1450610075033665908079998137985098905417567159/1427247692705959881058285969449495136382746624'),
            ('1440288121424753473054274698947691939714728571/1427247692705959881058285969449495136382746624', '1450610071023934846156560177799248777027778725/1427247692705959881058285969449495136382746624', '1440288121474159942458860978292821926287213713/1427247692705959881058285969449495136382746624'),
            ('1442608694112202437434320579151926875522494061/1427247692705959881058285969449495136382746624', '1442608697235061920706517405003658452508615727/1427247692705959881058285969449495136382746624', '721304347245564870086519332888488159543649225/713623846352979940529142984724747568191373312'),
            ('1444129095638046989922788727454301761905252387/1427247692705959881058285969449495136382746624', '1444129096238848633108457717495441127374286673/1427247692705959881058285969449495136382746624', '748631622240075371449454383391417410781710545/713623846352979940529142984724747568191373312'),
            ('722064548022820261939270433645709624387559161/713623846352979940529142984724747568191373312', '1497263246703718515711773456146060029775720747/1427247692705959881058285969449495136382746624', '1444129096026666888429870472218392909802115055/1427247692705959881058285969449495136382746624'),
            ('725305038094991027804960851593719923887987073/713623846352979940529142984724747568191373312', '720144060813167013001647129344708679988604829/713623846352979940529142984724747568191373312', '720144060256588013398642975283872081195132743/713623846352979940529142984724747568191373312'),
            ('374315810847038541070235417096341927282427655/356811923176489970264571492362373784095686656', '1444129096104211552809708963110303118934362171/1427247692705959881058285969449495136382746624', '361032273871109067592344136606867041885051461/356811923176489970264571492362373784095686656'),
            ('1600881643756061822818143719149627909290538941/1427247692705959881058285969449495136382746624', '134997473796755154998121018694470395745716807/178405961588244985132285746181186892047843328', '400220410922758537433304505689941664049734915/356811923176489970264571492362373784095686656'),
            ('200110205235141545196661542119826832728379523/178405961588244985132285746181186892047843328', '800440821057113147767366077855147077465603723/713623846352979940529142984724747568191373312', '1079979794241934493493644161330196057001965151/1427247692705959881058285969449495136382746624'),
            ('425727005055473661421662034280768701397083349/356811923176489970264571492362373784095686656', '875927037418108215096066381016960348337176519/1427247692705959881058285969449495136382746624', '851454010077814807628088147031823130616457641/713623846352979940529142984724747568191373312'),
            ('1702908022048291820706898200046285374102873141/1427247692705959881058285969449495136382746624', '1702908022285819311294736539731513809765957071/1427247692705959881058285969449495136382746624', '218981758478366019953083052494789221448097411/356811923176489970264571492362373784095686656'),
            ('1427250475386436542428419061165390988399349869/1427247692705959881058285969449495136382746624', '1427250483327197899916768493366996094250589485/1427247692705959881058285969449495136382746624', '713625237759994403762549747209182201952364311/713623846352979940529142984724747568191373312'),
            ('1427250491663913083147551651421959801250913347/1427247692705959881058285969449495136382746624', '713625237861488465899633479038288945141998161/713623846352979940529142984724747568191373312', '1427250475200028431264894496125010547163537163/1427247692705959881058285969449495136382746624'),
            ('356812619025603875337909666048250868181582907/356811923176489970264571492362373784095686656', '356812619497892459049913299696650191569119783/356811923176489970264571492362373784095686656', '713625238319520099965031145438902952548837071/713623846352979940529142984724747568191373312'),
            ('1427250473518430629625613814324323009685347387/1427247692705959881058285969449495136382746624', '356812618258311696038886547749443245589274169/356811923176489970264571492362373784095686656', '713625236471496905757853050951306462044008213/713623846352979940529142984724747568191373312'),
            ('1427250475244243575377245951568704928464433479/1427247692705959881058285969449495136382746624', '713625237933584067704040977683864471706353591/713623846352979940529142984724747568191373312', '1427250486312651689860367384396722850312235179/1427247692705959881058285969449495136382746624'),
            ('356812815131829320956598867512774473018995393/356811923176489970264571492362373784095686656', '356812622624370542631932055386015007872379389/356811923176489970264571492362373784095686656', '1427250490021750365095835902149057989621498343/1427247692705959881058285969449495136382746624'),
            ('1427250490310400402168525279265332213933781377/1427247692705959881058285969449495136382746624', '713625634449284801173833023954662104663838907/713623846352979940529142984724747568191373312', '1427250490273505422425124682122554241257254035/1427247692705959881058285969449495136382746624'),
            ('356812622550584243710937361458606099357388199/356811923176489970264571492362373784095686656', '1427250490641123625750001213717319817581719789/1427247692705959881058285969449495136382746624', '1427251265594439589467002812057849908009966211/1427247692705959881058285969449495136382746624'),
        ),
    },
)


def run_known_answer_controls() -> dict[str, Any]:
    """Run exact known-answer algebraic & negative controls (zero floats)."""
    # 1. Exact analytic Hessian at origin for charges +1 at (+/-1, 0, 0): diag(4, -2, -2) -> inertia (1, 2, 0)
    h_n2 = [
        [Fraction(4, 1), Fraction(0, 1), Fraction(0, 1)],
        [Fraction(0, 1), Fraction(-2, 1), Fraction(0, 1)],
        [Fraction(0, 1), Fraction(0, 1), Fraction(-2, 1)],
    ]
    n2_inertia_ok = exact_lagrange_inertia_3x3_q(h_n2) == (1, 2, 0)

    # 2. Exact hyperbolic zero-diagonal 2x2 block [[0, 1, 0], [1, 0, 0], [0, 0, -3]] -> inertia (1, 2, 0)
    h_hyp = [
        [Fraction(0, 1), Fraction(1, 1), Fraction(0, 1)],
        [Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)],
        [Fraction(0, 1), Fraction(0, 1), Fraction(-3, 1)],
    ]
    hyp_inertia_ok = exact_lagrange_inertia_3x3_q(h_hyp) == (1, 2, 0)

    # 3. Negative control: overlapping boxes rejected
    overlap_rejected = False
    try:
        certify_configuration(
            "neg-overlap",
            [(Fraction(-1, 1), Fraction(0, 1), Fraction(0, 1)), (Fraction(1, 1), Fraction(0, 1), Fraction(0, 1))],
            [Fraction(1, 1), Fraction(1, 1)],
            [(Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)), (Fraction(1, 1000), Fraction(0, 1), Fraction(0, 1))],
            [Fraction(1, 100), Fraction(1, 100)],
        )
    except PointChargeFrontierError as exc:
        overlap_rejected = exc.code == "OVERLAPPING_BOXES"

    # 4. Negative control: overly wide non-contracting box rejected
    wide_rejected = False
    try:
        certify_configuration(
            "neg-wide",
            [(Fraction(-1, 1), Fraction(0, 1), Fraction(0, 1)), (Fraction(1, 1), Fraction(0, 1), Fraction(0, 1))],
            [Fraction(1, 1), Fraction(1, 1)],
            [(Fraction(0, 1), Fraction(0, 1), Fraction(0, 1))],
            [Fraction(1, 2)],
        )
    except PointChargeFrontierError as exc:
        wide_rejected = exc.code == "KRAWCZYK_NOT_CONTRACTIVE"

    all_ok = all((n2_inertia_ok, hyp_inertia_ok, overlap_rejected, wide_rejected))
    return {
        "n2_analytic_hessian_inertia_exact": n2_inertia_ok,
        "hyperbolic_zero_diag_inertia_exact": hyp_inertia_ok,
        "overlapping_box_negative_control_rejected": overlap_rejected,
        "non_contractive_box_negative_control_rejected": wide_rejected,
        "all_controls_passed": all_ok,
    }


def _config_passes_frontier(cfg: dict[str, Any], target_n: int) -> bool:
    return all(
        (
            cfg["n_certified"] == target_n,
            cfg["exceeds_maxwell_bound"],
            cfg["morse_relation_matches_n_minus_1"],
        )
    )


def _certify_canonical_config(raw_cfg: dict[str, Any]) -> dict[str, Any]:
    sites_q = [(parse_q_str(s[0]), parse_q_str(s[1]), parse_q_str(s[2])) for s in raw_cfg["sites"]]
    charges_q = [parse_q_str(q) for q in raw_cfg["charges"]]
    centers_q = [(parse_q_str(c[0]), parse_q_str(c[1]), parse_q_str(c[2])) for c in raw_cfg["centers"]]
    hw_q = parse_q_str(raw_cfg["half_width"])
    hws_q = [hw_q] * len(centers_q)
    res = certify_configuration(raw_cfg["name"], sites_q, charges_q, centers_q, hws_q)
    res["description"] = raw_cfg["description"]
    return res


def build_certificate() -> dict[str, Any]:
    """Build the canonical self-signed pure-Z/Q Track S21d certificate across all 4 witness configs (127 equilibria)."""
    controls = run_known_answer_controls()
    configs: dict[str, dict[str, Any]] = {}
    for raw_cfg in CANONICAL_WITNESS_CONFIGS:
        certified = _certify_canonical_config(raw_cfg)
        configs[certified["name"]] = certified

    c7 = configs["n7-abk-second-pair"]
    c6 = configs["n6-abk-milnor8-flat"]
    c8 = configs["n8-abk-milnor8-flat"]
    total_eq = sum(c["n_certified"] for c in configs.values())

    payload: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "track_id": TRACK_ID,
        "arithmetic_mode": "pure_Z_Q_384bit_outward_dyadic_zero_float",
        "prec_bits": PREC_BITS,
        "known_answer_controls": controls,
        "configurations": configs,
        "frontier_verdicts": {
            "issue12_c7_44_certified": _config_passes_frontier(c7, 44),
            "n6_31_new_record_certified": _config_passes_frontier(c6, 31),
            "issue5_n8_51_maxwell_counterexample_certified": _config_passes_frontier(c8, 51),
            "total_equilibria_certified": total_eq,
            "scope_honesty_note": (
                "Local pure-Z/Q Krawczyk uniqueness certifies exact lower bounds on the number of pairwise-disjoint "
                "non-degenerate equilibria (N(C_7)>=44>36, N(C_6)>=31>29>25, N(C_8)>=51>49=(8-1)^2), which strictly "
                "refutes Maxwell's 1873 bound (n-1)^2 for n=8 without claiming global completeness over R^3 minus {a_j}."
            ),
        },
    }
    payload_sha = sha256_canonical(payload)
    cert_sha = sha256_canonical(
        {
            "schema_version": SCHEMA_VERSION,
            "track_id": TRACK_ID,
            "payload_sha256": payload_sha,
        }
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "track_id": TRACK_ID,
        "payload_sha256": payload_sha,
        "certificate_sha256": cert_sha,
        "payload": payload,
    }


def _validate_serialized_config(cfg_name: str, cfg: dict[str, Any], *, full_replay: bool) -> int:
    c_map = require_mapping(cfg, f"configurations.{cfg_name}")
    expected_keys = {
        "name",
        "description",
        "n_charges",
        "maxwell_bound",
        "n_certified",
        "exceeds_maxwell_bound",
        "morse_counts",
        "m2_minus_m1",
        "morse_relation_matches_n_minus_1",
        "max_rho_q",
        "sites",
        "charges",
        "boxes",
    }
    if set(c_map.keys()) != expected_keys:
        raise PointChargeFrontierError("INVALID_CONFIG_KEYS", f"{cfg_name} keys mismatch")
    if require_str(c_map["name"], f"{cfg_name}.name") != cfg_name:
        raise PointChargeFrontierError("CONFIG_NAME_MISMATCH", f"{cfg_name} name mismatch")
    require_str(c_map["description"], f"{cfg_name}.description", min_len=1, max_len=512)
    n_charges = require_int(c_map["n_charges"], f"{cfg_name}.n_charges", lo=2, hi=MAX_CHARGES)
    maxwell = require_int(c_map["maxwell_bound"], f"{cfg_name}.maxwell_bound", lo=1, hi=(MAX_CHARGES - 1) ** 2)
    if maxwell != (n_charges - 1) * (n_charges - 1):
        raise PointChargeFrontierError("MAXWELL_BOUND_MISMATCH", f"{cfg_name} maxwell_bound != (n-1)^2")
    n_cert = require_int(c_map["n_certified"], f"{cfg_name}.n_certified", lo=1, hi=MAX_BOXES)
    exceeds = require_bool(c_map["exceeds_maxwell_bound"], f"{cfg_name}.exceeds_maxwell_bound")
    if exceeds != (n_cert > maxwell):
        raise PointChargeFrontierError("EXCEEDS_BOUND_MISMATCH", f"{cfg_name} exceeds_maxwell_bound inconsistent")
    morse_seq = require_seq(c_map["morse_counts"], f"{cfg_name}.morse_counts", min_len=4, max_len=4)
    morse = [require_int(v, f"{cfg_name}.morse[{i}]", lo=0, hi=MAX_BOXES) for i, v in enumerate(morse_seq)]
    if morse[0] != 0 or morse[3] != 0 or sum(morse) != n_cert:
        raise PointChargeFrontierError("MORSE_COUNTS_INCONSISTENT", f"{cfg_name} morse_counts inconsistent")
    m2_m1 = require_int(c_map["m2_minus_m1"], f"{cfg_name}.m2_minus_m1", lo=-MAX_BOXES, hi=MAX_BOXES)
    if m2_m1 != morse[2] - morse[1]:
        raise PointChargeFrontierError("MORSE_DIFF_MISMATCH", f"{cfg_name} m2_minus_m1 inconsistent")
    rel_ok = require_bool(c_map["morse_relation_matches_n_minus_1"], f"{cfg_name}.morse_relation_matches_n_minus_1")
    if rel_ok != (m2_m1 == n_charges - 1):
        raise PointChargeFrontierError("MORSE_REL_MISMATCH", f"{cfg_name} morse_relation flag inconsistent")
    max_rho = parse_q_str(c_map["max_rho_q"], f"{cfg_name}.max_rho_q")
    if max_rho < 0 or max_rho >= 1:
        raise PointChargeFrontierError("INVALID_MAX_RHO", f"{cfg_name} max_rho_q must lie in [0, 1)")

    sites_seq = require_seq(c_map["sites"], f"{cfg_name}.sites", min_len=n_charges, max_len=n_charges)
    charges_seq = require_seq(c_map["charges"], f"{cfg_name}.charges", min_len=n_charges, max_len=n_charges)
    boxes_seq = require_seq(c_map["boxes"], f"{cfg_name}.boxes", min_len=n_cert, max_len=n_cert)

    sites_q: list[tuple[Fraction, Fraction, Fraction]] = []
    for i, s in enumerate(sites_seq):
        cs = require_seq(s, f"{cfg_name}.sites[{i}]", min_len=3, max_len=3)
        sites_q.append((parse_q_str(cs[0]), parse_q_str(cs[1]), parse_q_str(cs[2])))
    charges_q = [parse_q_str(q, f"{cfg_name}.charges[{i}]") for i, q in enumerate(charges_seq)]
    for q in charges_q:
        if q <= 0:
            raise PointChargeFrontierError("NON_POSITIVE_CHARGE", f"{cfg_name} contains non-positive charge")

    centers_q: list[tuple[Fraction, Fraction, Fraction]] = []
    hws_q: list[Fraction] = []
    box_morse = [0, 0, 0, 0]
    box_max_rho = Fraction(0, 1)
    for i, b in enumerate(boxes_seq):
        bm = require_mapping(b, f"{cfg_name}.boxes[{i}]")
        if set(bm.keys()) != {"center", "half_width", "morse_index", "rho_upper_q"}:
            raise PointChargeFrontierError("INVALID_BOX_KEYS", f"{cfg_name}.boxes[{i}] keys mismatch")
        cs = require_seq(bm["center"], f"{cfg_name}.boxes[{i}].center", min_len=3, max_len=3)
        centers_q.append((parse_q_str(cs[0]), parse_q_str(cs[1]), parse_q_str(cs[2])))
        hw = parse_q_str(bm["half_width"], f"{cfg_name}.boxes[{i}].half_width")
        if hw <= 0:
            raise PointChargeFrontierError("NON_POSITIVE_HALFWIDTH", f"{cfg_name}.boxes[{i}] half_width <= 0")
        hws_q.append(hw)
        mi = require_int(bm["morse_index"], f"{cfg_name}.boxes[{i}].morse_index", lo=1, hi=2)
        box_morse[mi] += 1
        rho = parse_q_str(bm["rho_upper_q"], f"{cfg_name}.boxes[{i}].rho_upper_q")
        if rho < 0 or rho >= 1:
            raise PointChargeFrontierError("INVALID_BOX_RHO", f"{cfg_name}.boxes[{i}] rho_upper_q out of [0,1)")
        if rho > box_max_rho:
            box_max_rho = rho
    if box_morse != morse or box_max_rho != max_rho:
        raise PointChargeFrontierError("BOX_SUMMARY_MISMATCH", f"{cfg_name} box summaries disagree with config header")

    if full_replay:
        replayed = certify_configuration(cfg_name, sites_q, charges_q, centers_q, hws_q)
        replayed["description"] = c_map["description"]
        if canonical_json_bytes(replayed) != canonical_json_bytes(c_map):
            raise PointChargeFrontierError("KRAWCZYK_REPLAY_MISMATCH", f"{cfg_name} failed independent pure-Z/Q Krawczyk replay")
    return n_cert


def validate_certificate_payload(
    cert: Any,
    *,
    full_replay: bool = True,
) -> dict[str, Any]:
    """Fail-closed Zone-B validation of a Track S21d certificate (full Krawczyk replay by default)."""
    do_full = require_bool(full_replay, "full_replay")
    root = require_mapping(cert, "certificate")
    validate_json_value(root)
    if set(root.keys()) != {"schema_version", "track_id", "payload_sha256", "certificate_sha256", "payload"}:
        raise PointChargeFrontierError("INVALID_CERT_KEYS", "certificate top-level keys mismatch")
    if require_str(root["schema_version"], "schema_version") != SCHEMA_VERSION:
        raise PointChargeFrontierError("INVALID_SCHEMA_VERSION", "unexpected schema_version")
    if require_str(root["track_id"], "track_id") != TRACK_ID:
        raise PointChargeFrontierError("INVALID_TRACK_ID", "unexpected track_id")

    payload = require_mapping(root["payload"], "payload")
    computed_payload_sha = sha256_canonical(payload)
    if not constant_time_hex64_equal(root["payload_sha256"], computed_payload_sha, "payload_sha256"):
        raise PointChargeFrontierError("PAYLOAD_SHA_MISMATCH", "payload_sha256 does not match canonical payload")
    computed_cert_sha = sha256_canonical(
        {
            "schema_version": SCHEMA_VERSION,
            "track_id": TRACK_ID,
            "payload_sha256": computed_payload_sha,
        }
    )
    if not constant_time_hex64_equal(root["certificate_sha256"], computed_cert_sha, "certificate_sha256"):
        raise PointChargeFrontierError("CERT_SHA_MISMATCH", "certificate_sha256 does not match envelope")

    expected_payload_keys = {
        "schema_version",
        "track_id",
        "arithmetic_mode",
        "prec_bits",
        "known_answer_controls",
        "configurations",
        "frontier_verdicts",
    }
    if set(payload.keys()) != expected_payload_keys:
        raise PointChargeFrontierError("INVALID_PAYLOAD_KEYS", "payload keys mismatch")
    if require_str(payload["schema_version"], "payload.schema_version") != SCHEMA_VERSION:
        raise PointChargeFrontierError("INVALID_PAYLOAD_SCHEMA", "payload schema_version mismatch")
    if require_str(payload["track_id"], "payload.track_id") != TRACK_ID:
        raise PointChargeFrontierError("INVALID_PAYLOAD_TRACK", "payload track_id mismatch")
    if require_str(payload["arithmetic_mode"], "arithmetic_mode") != "pure_Z_Q_384bit_outward_dyadic_zero_float":
        raise PointChargeFrontierError("INVALID_ARITHMETIC_MODE", "unexpected arithmetic_mode")
    if require_int(payload["prec_bits"], "prec_bits", lo=PREC_BITS, hi=PREC_BITS) != PREC_BITS:
        raise PointChargeFrontierError("INVALID_PREC_BITS", "unexpected prec_bits")

    controls = require_mapping(payload["known_answer_controls"], "known_answer_controls")
    if not require_bool(controls.get("all_controls_passed"), "all_controls_passed"):
        raise PointChargeFrontierError("CONTROLS_FAILED", "known_answer_controls.all_controls_passed is False")

    configs = require_mapping(payload["configurations"], "configurations")
    if set(configs.keys()) != set(_CANONICAL_CONFIG_NAMES):
        raise PointChargeFrontierError("CONFIG_SET_MISMATCH", "configurations keys mismatch")

    total_cert = 0
    for cname in _CANONICAL_CONFIG_NAMES:
        total_cert += _validate_serialized_config(cname, configs[cname], full_replay=do_full)

    verdicts = require_mapping(payload["frontier_verdicts"], "frontier_verdicts")
    expected_verdict_keys = {
        "issue12_c7_44_certified",
        "n6_31_new_record_certified",
        "issue5_n8_51_maxwell_counterexample_certified",
        "total_equilibria_certified",
        "scope_honesty_note",
    }
    if set(verdicts.keys()) != expected_verdict_keys:
        raise PointChargeFrontierError("VERDICT_KEYS_MISMATCH", "frontier_verdicts keys mismatch")
    require_str(verdicts["scope_honesty_note"], "scope_honesty_note", min_len=20, max_len=1024)
    if require_int(verdicts["total_equilibria_certified"], "total_equilibria_certified", lo=1, hi=MAX_BOXES * MAX_CONFIGS) != total_cert:
        raise PointChargeFrontierError("TOTAL_EQ_MISMATCH", "total_equilibria_certified disagrees with sum of configs")

    c7 = configs["n7-abk-second-pair"]
    c6 = configs["n6-abk-milnor8-flat"]
    c8 = configs["n8-abk-milnor8-flat"]
    exp_c7 = _config_passes_frontier(c7, 44)
    exp_c6 = _config_passes_frontier(c6, 31)
    exp_c8 = _config_passes_frontier(c8, 51)
    if require_bool(verdicts["issue12_c7_44_certified"], "issue12_c7_44_certified") != exp_c7:
        raise PointChargeFrontierError("VERDICT_MISMATCH", "issue12_c7_44_certified inconsistent")
    if require_bool(verdicts["n6_31_new_record_certified"], "n6_31_new_record_certified") != exp_c6:
        raise PointChargeFrontierError("VERDICT_MISMATCH", "n6_31_new_record_certified inconsistent")
    if require_bool(verdicts["issue5_n8_51_maxwell_counterexample_certified"], "issue5_n8_51_maxwell_counterexample_certified") != exp_c8:
        raise PointChargeFrontierError("VERDICT_MISMATCH", "issue5_n8_51_maxwell_counterexample_certified inconsistent")

    if do_full:
        expected_cert = build_certificate()
        if canonical_json_bytes(expected_cert) != canonical_json_bytes(root):
            raise PointChargeFrontierError("CANONICAL_REPLAY_MISMATCH", "certificate does not match canonical Track S21d replay")

    return {
        "valid": True,
        "track_id": TRACK_ID,
        "payload_sha256": computed_payload_sha,
        "certificate_sha256": computed_cert_sha,
        "total_equilibria_certified": total_cert,
        "full_replay": do_full,
    }
