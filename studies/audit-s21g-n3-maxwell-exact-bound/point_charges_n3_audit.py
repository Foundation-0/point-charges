#!/usr/bin/env python3
"""Zone-B R2-grade hardened standalone CLI, full pure-Z/Q replay & deep secondary auditor
for Foundation-0/point-charges Issue #4 (Three-Charge Maxwell Bound Exact Study).

100% pure Python 3 stdlib, zero IEEE-754 floats, zero bare asserts, zero .isdigit(),
zero mathematical RNG (secrets.token_hex used strictly for atomic temp filename uniqueness).
"""
from __future__ import annotations

import errno
from fractions import Fraction
import hashlib
import hmac
import json
import math
import os
from pathlib import Path
import secrets
import stat
import sys
import threading
from typing import Any, Sequence

_THIS_DIR = str(Path(__file__).resolve().parent)
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import point_charges_n3_oracle as oracle  # noqa: E402

MAX_JSON_BYTES = 1_048_576
MAX_ALLOWED_CAP_BYTES = 64 * 1024 * 1024
_READ_CHUNK_BYTES = 65_536
MAX_JSON_DEPTH = 64
MAX_JSON_INT_DIGITS = oracle.MAX_INT_DIGITS
_EXACT_PATH_TYPES = (Path, type(Path()))
_HEX_LOWER = frozenset("0123456789abcdef")
_ASCII_DIGITS = frozenset("0123456789")
_CERT_KEYS = frozenset({"schema_version", "track_id", "payload_sha256", "certificate_sha256", "payload"})
_PAYLOAD_KEYS = frozenset(
    {
        "track_id",
        "issue_ref",
        "novelty_class",
        "pillar1_algebraic_identities",
        "pillar2_incenter_universal_theorem",
        "pillar3_isosceles_axis_uniqueness",
        "pillar4_gale_nikaido_compact_univalence",
        "pillar5_krawczyk_configurations",
        "summary",
    }
)
_P1_KEYS = frozenset(
    {
        "algebraic_identity_cases",
        "out_of_plane_vz_probe_q",
        "out_of_plane_z_restoring_nonzero",
        "morse_relation_m2_minus_m1_eq_2",
        "all_verified",
    }
)
_P2_KEYS = frozenset({"weierstrass_incenter_witnesses", "closes_487_flat_triangle_nonemptiness_gap", "all_verified"})
_P3_KEYS = frozenset({"isosceles_axis_cases", "at_most_one_on_axis_minimum_all_h_gt_0", "all_verified"})
_P4_KEYS = frozenset({"compact_4d_boxes", "total_4d_subboxes_certified", "all_verified"})
_P5_KEYS = frozenset({"configurations", "total_equilibria_certified", "all_verified"})
_P1_ROW_KEYS = frozenset(
    {
        "name",
        "det_h_tilde_q",
        "det_h_tilde_dec6",
        "det_dl_q",
        "phase_formula_matches_pmin",
        "dl_jacobian_identity_exact",
        "in_R_min",
        "verified",
    }
)
_P2_ROW_KEYS = frozenset(
    {
        "name",
        "v1_q",
        "v2_q",
        "v3_q",
        "aspect_h_over_max_side_dec6",
        "grad_exact_zero",
        "z2_phase_sum_exact_zero",
        "h_tilde_exact_half_I2",
        "det_h_tilde_q",
        "pmin_eq_one_eighth_rprod",
        "verified",
    }
)
_P3_ROW_KEYS = frozenset(
    {
        "k",
        "height_h_q",
        "y_fold_q",
        "y_pitchfork_q",
        "y_pitchfork_eq_2_y_fold",
        "vyy_fold_exact_zero",
        "vxx_pitchfork_exact_zero",
        "q_prime_strictly_neg_at_mid",
        "vxx_strictly_pos_at_mid",
        "incenter_strictly_between_fold_and_pitchfork",
        "verified",
    }
)
_P4_ROW_KEYS = frozenset(
    {
        "name",
        "a2_range_q",
        "b2_range_q",
        "c2_q",
        "u1_range_q",
        "u2_range_q",
        "subboxes_checked",
        "subboxes_passed",
        "min_pmin_lo_dec6",
        "min_s00_lo_dec6",
        "min_det_sym_lo_dec6",
        "strictly_injective_on_4d_box",
    }
)
_P5_CFG_KEYS = frozenset(
    {
        "name",
        "charges_q",
        "equilibria_certified",
        "m1_local_minima",
        "m2_saddles",
        "m2_minus_m1",
        "pairwise_disjoint_pairs",
        "all_boxes_pairwise_disjoint",
        "equilibria",
        "verified",
    }
)
_P5_EQ_KEYS = frozenset(
    {
        "center_q",
        "half_width_q",
        "box_lo_q",
        "box_hi_q",
        "preconditioner_det_nonzero",
        "krawczyk_rho_ub_dec6",
        "krawczyk_contraction_lt_1",
        "krawczyk_strict_interior",
        "trace_strictly_positive",
        "det_h_sign",
        "morse_index_3d",
    }
)
_SUMMARY_KEYS = frozenset(
    {
        "pillar_count",
        "all_five_pillars_verified",
        "total_4d_univalence_subboxes",
        "total_krawczyk_equilibria",
        "global_n3_conjecture_open_outside_compact_boxes",
        "float_count_on_trust_path",
        "bare_assert_count",
        "rng_calls",
    }
)
_DIR_OPEN_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
_TEMP_NAME_LOCK = threading.Lock()
_TEMP_COUNTER = 0


# ============================================================================
# 1. Zone-B Path, CWD, Atomic dir_fd Traversal & Durable 0o600 Atomic Write
# ============================================================================

def _require_utf8_no_surrogate(text: Any, label: str) -> str:
    if type(text) is not str:
        raise oracle.PointChargesN3Error("INVALID_STR", f"{label}: exact str required")
    oracle.require_ok("\x00" not in text, "NUL_BYTE_FORBIDDEN", f"{label}: NUL byte forbidden")
    try:
        text.encode("utf-8")
    except UnicodeEncodeError as exc:
        raise oracle.PointChargesN3Error("INVALID_UTF8_SURROGATE", f"{label}: lone UTF-16 surrogate forbidden") from exc
    return text


def _require_concrete_path(path: Any, label: str) -> Path:
    oracle.require_ok(type(path) in _EXACT_PATH_TYPES, "INVALID_PATH_TYPE", f"{label}: concrete pathlib.Path required")
    _require_utf8_no_surrogate(str(path), label)
    return path


def _safe_cwd(label: str) -> Path:
    try:
        return Path.cwd()
    except OSError as exc:
        raise oracle.PointChargesN3Error("CWD_UNAVAILABLE", f"{label}: CWD unavailable ({type(exc).__name__})") from exc


def _check_pwd_no_symlinks(label: str) -> None:
    pwd = os.environ.get("PWD", "")
    if pwd.startswith("/"):
        _require_utf8_no_surrogate(pwd, f"{label} PWD")
        pwd_path = Path(pwd)
        try:
            dfd = os.open("/", _DIR_OPEN_FLAGS)
        except OSError as exc:
            raise oracle.PointChargesN3Error("CWD_UNAVAILABLE", f"{label}: cannot open root for PWD") from exc
        try:
            for part in pwd_path.parts[1:]:
                if part == ".":
                    continue
                try:
                    next_fd = os.open(part, _DIR_OPEN_FLAGS, dir_fd=dfd)
                except OSError as exc:
                    if exc.errno in (errno.ELOOP, errno.ENOTDIR):
                        raise oracle.PointChargesN3Error(
                            "SYMLINK_ANCESTOR_FORBIDDEN", f"{label}: symlink PWD ancestor forbidden"
                        ) from exc
                    raise oracle.PointChargesN3Error(
                        "CWD_UNAVAILABLE", f"{label}: PWD ancestor unavailable ({type(exc).__name__})"
                    ) from exc
                os.close(dfd)
                dfd = next_fd
        finally:
            os.close(dfd)


def _open_dir_chain(target: Path, label: str, *, create_parents: bool) -> tuple[int, str]:
    """Open ancestor directory chain atomically via dir_fd + O_DIRECTORY | O_NOFOLLOW."""
    if target.is_absolute():
        parts = target.parts[1:]
    else:
        cwd = _safe_cwd(label)
        _check_pwd_no_symlinks(label)
        parts = cwd.parts[1:] + target.parts
    oracle.require_ok(
        len(parts) >= 1 and parts[-1] not in ("", ".", ".."),
        "NOT_A_REGULAR_FILE",
        f"{label}: path does not name a regular file leaf",
    )
    parent_parts = parts[:-1]
    leaf = parts[-1]
    err_code = "FILE_WRITE_ERROR" if create_parents else "FILE_READ_ERROR"
    try:
        dfd = os.open("/", _DIR_OPEN_FLAGS)
    except OSError as exc:
        raise oracle.PointChargesN3Error(err_code, f"{label}: cannot open root") from exc
    try:
        for part in parent_parts:
            if part == ".":
                continue
            try:
                next_fd = os.open(part, _DIR_OPEN_FLAGS, dir_fd=dfd)
            except OSError as exc:
                if create_parents and exc.errno == errno.ENOENT:
                    try:
                        os.mkdir(part, 0o755, dir_fd=dfd)
                    except FileExistsError:
                        pass
                    except OSError as mk_exc:
                        raise oracle.PointChargesN3Error(
                            "FILE_WRITE_ERROR", f"{label}: mkdir failed ({type(mk_exc).__name__})"
                        ) from mk_exc
                    try:
                        next_fd = os.open(part, _DIR_OPEN_FLAGS, dir_fd=dfd)
                    except OSError as re_exc:
                        raise oracle.PointChargesN3Error(
                            "FILE_WRITE_ERROR", f"{label}: cannot open created dir ({type(re_exc).__name__})"
                        ) from re_exc
                elif exc.errno in (errno.ELOOP, errno.ENOTDIR):
                    raise oracle.PointChargesN3Error(
                        "SYMLINK_FORBIDDEN", f"{label}: symlink or non-directory ancestor forbidden"
                    ) from exc
                else:
                    raise oracle.PointChargesN3Error(
                        err_code, f"{label}: cannot traverse ancestor ({type(exc).__name__})"
                    ) from exc
            os.close(dfd)
            dfd = next_fd
        return dfd, leaf
    except BaseException:
        os.close(dfd)
        raise


def read_regular_bytes_capped(path: Path, label: str = "file", max_bytes: int = MAX_JSON_BYTES) -> bytes:
    """Read regular file bytes under strict cap using atomic dir_fd openat + fstat S_ISREG + 64 KiB chunks."""
    target = _require_concrete_path(path, label)
    oracle.require_exact_int(max_bytes, "max_bytes", lo=1, hi=MAX_ALLOWED_CAP_BYTES)
    dfd, leaf = _open_dir_chain(target, label, create_parents=False)
    flags = os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK
    try:
        try:
            fd = os.open(leaf, flags, dir_fd=dfd)
        except OSError as exc:
            if exc.errno == errno.ELOOP:
                raise oracle.PointChargesN3Error("SYMLINK_FORBIDDEN", f"{label}: symlink leaf forbidden") from exc
            raise oracle.PointChargesN3Error(
                "FILE_READ_ERROR", f"{label}: cannot open ({type(exc).__name__})"
            ) from exc
    finally:
        os.close(dfd)

    chunks: list[bytes] = []
    total = 0
    try:
        st_f = os.fstat(fd)
        oracle.require_ok(stat.S_ISREG(st_f.st_mode), "NOT_A_REGULAR_FILE", f"{label}: not a regular file")
        while total <= max_bytes:
            chunk = os.read(fd, min(_READ_CHUNK_BYTES, max_bytes + 1 - total))
            if not chunk:
                break
            chunks.append(chunk)
            total += len(chunk)
    except oracle.PointChargesN3Error:
        raise
    except OSError as exc:
        raise oracle.PointChargesN3Error("FILE_READ_ERROR", f"{label}: read failed ({type(exc).__name__})") from exc
    finally:
        os.close(fd)
    oracle.require_ok(total <= max_bytes, "FILE_TOO_LARGE", f"{label}: exceeds {max_bytes} bytes")
    return b"".join(chunks)


def write_regular_bytes_safe(path: Path, payload: bytes, label: str = "output") -> None:
    """Durable atomic dir_fd write: 0o600 temp, lock+counter+CSPRNG+EEXIST retry, fsync, st_nlink==1, os.replace, dir fsync."""
    global _TEMP_COUNTER
    oracle.require_ok(type(payload) is bytes, "INVALID_PAYLOAD_BYTES", f"{label}: payload must be exact bytes")
    target = _require_concrete_path(path, label)
    dfd, leaf = _open_dir_chain(target, label, create_parents=True)
    tmp_name: str | None = None
    try:
        try:
            st_existing = os.lstat(leaf, dir_fd=dfd)
        except FileNotFoundError:
            st_existing = None
        except OSError as exc:
            raise oracle.PointChargesN3Error(
                "FILE_WRITE_ERROR", f"{label}: lstat failed ({type(exc).__name__})"
            ) from exc
        if st_existing is not None:
            oracle.require_ok(not stat.S_ISLNK(st_existing.st_mode), "SYMLINK_FORBIDDEN", f"{label}: symlink leaf forbidden")
            oracle.require_ok(stat.S_ISREG(st_existing.st_mode), "NOT_A_REGULAR_FILE", f"{label}: non-regular target")
            oracle.require_ok(st_existing.st_nlink <= 1, "HARDLINK_FORBIDDEN", f"{label}: hardlinked target forbidden")

        pid = os.getpid()
        for _ in range(64):
            with _TEMP_NAME_LOCK:
                _TEMP_COUNTER += 1
                seq_num = _TEMP_COUNTER
                cand = f".n3.{pid}.{seq_num}.{secrets.token_hex(8)}.tmp"
            try:
                fd = os.open(
                    cand,
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | os.O_NONBLOCK,
                    0o600,
                    dir_fd=dfd,
                )
                tmp_name = cand
                break
            except FileExistsError:
                continue
            except OSError as exc:
                raise oracle.PointChargesN3Error(
                    "FILE_WRITE_ERROR", f"{label}: create temp failed ({type(exc).__name__})"
                ) from exc
        else:
            raise oracle.PointChargesN3Error("FILE_WRITE_ERROR", f"{label}: exhausted temp file slots")

        try:
            view = memoryview(payload)
            written = 0
            while written < len(view):
                written += os.write(fd, view[written:])
            os.fsync(fd)
            st_tmp = os.fstat(fd)
            oracle.require_ok(st_tmp.st_nlink == 1, "HARDLINK_FORBIDDEN", f"{label}: temp link count altered")
        except oracle.PointChargesN3Error:
            raise
        except OSError as exc:
            raise oracle.PointChargesN3Error(
                "FILE_WRITE_ERROR", f"{label}: write failed ({type(exc).__name__})"
            ) from exc
        finally:
            os.close(fd)

        try:
            os.replace(tmp_name, leaf, src_dir_fd=dfd, dst_dir_fd=dfd)
            tmp_name = None
            os.fsync(dfd)
        except OSError as exc:
            raise oracle.PointChargesN3Error(
                "FILE_WRITE_ERROR", f"{label}: atomic replace failed ({type(exc).__name__})"
            ) from exc
    finally:
        if tmp_name is not None:
            try:
                os.unlink(tmp_name, dir_fd=dfd)
            except OSError:
                pass
        os.close(dfd)


# ============================================================================
# 2. Strict ASCII Parsers, Aligned 3,321-Bit Validation & Canonical JSON
# ============================================================================

def require_hex64_lower_ascii(val: Any, label: str = "sha256") -> str:
    s = _require_utf8_no_surrogate(val, label)
    oracle.require_ok(len(s) == 64, "INVALID_HEX64", f"{label}: expected 64-char lowercase ASCII hex")
    for ch in s:
        oracle.require_ok(ch in _HEX_LOWER, "INVALID_HEX64", f"{label}: non-hex-lowercase char")
    return s


def _parse_bounded_int(token: str) -> int:
    s = _require_utf8_no_surrogate(token, "json_int")
    oracle.require_ok(s != "-0", "INVALID_JSON_INT", "negative zero forbidden")
    digits = s[1:] if s.startswith("-") else s
    oracle.require_ok(len(digits) >= 1, "INVALID_JSON_INT", "JSON int has no digits")
    oracle.require_ok(len(digits) <= MAX_JSON_INT_DIGITS, "INT_DIGIT_LIMIT", "JSON int exceeds digit cap")
    oracle.require_ok(digits == "0" or not digits.startswith("0"), "INVALID_JSON_INT", "leading zero forbidden")
    for ch in digits:
        oracle.require_ok(ch in _ASCII_DIGITS, "INVALID_JSON_INT", "non-ASCII digit forbidden")
    val = int(s)
    oracle.require_ok(val.bit_length() <= oracle.MAX_INT_BITS, "INT_BIT_LIMIT", "JSON int exceeds 3,321-bit cap")
    return val


def parse_exact_rational(token: Any, label: str = "rational") -> Fraction:
    s = _require_utf8_no_surrogate(token, label)
    parts = s.split("/")
    oracle.require_ok(len(parts) in (1, 2), "INVALID_RATIONAL", f"{label}: expected 'n' or 'n/d'")
    num = _parse_bounded_int(parts[0])
    if len(parts) == 1:
        return oracle.require_exact_fraction(Fraction(num, 1), label)
    den = _parse_bounded_int(parts[1])
    oracle.require_ok(den >= 2, "INVALID_RATIONAL", f"{label}: denominator must be >= 2")
    oracle.require_ok(math.gcd(abs(num), den) == 1, "NON_CANONICAL_RATIONAL", f"{label}: not in lowest terms")
    return oracle.require_exact_fraction(Fraction(num, den), label)


def _validate_json_tree(root: Any, label: str = "json") -> None:
    stack: list[tuple[Any, int]] = [(root, 0)]
    while stack:
        node, depth = stack.pop()
        oracle.require_ok(depth <= MAX_JSON_DEPTH, "JSON_DEPTH_LIMIT", f"{label}: nesting exceeds {MAX_JSON_DEPTH}")
        t = type(node)
        if node is None or t is bool:
            continue
        if t is int:
            oracle.require_ok(node.bit_length() <= oracle.MAX_INT_BITS, "INT_BIT_LIMIT", f"{label}: int exceeds 3,321 bits")
            continue
        if t is str:
            _require_utf8_no_surrogate(node, label)
            continue
        if t is list:
            nxt = depth + 1
            for item in node:
                stack.append((item, nxt))
            continue
        if t is dict:
            nxt = depth + 1
            for k, v in node.items():
                _require_utf8_no_surrogate(k, f"{label}.key")
                stack.append((v, nxt))
            continue
        raise oracle.PointChargesN3Error("FORBIDDEN_JSON_TYPE", f"{label}: unsupported type {t.__name__}")


def _reject_float(token: str) -> Any:
    raise oracle.PointChargesN3Error("FLOAT_LITERAL_FORBIDDEN", f"float literal forbidden: {token}")


def _reject_constant(token: str) -> Any:
    raise oracle.PointChargesN3Error("NONFINITE_LITERAL_FORBIDDEN", f"non-finite constant forbidden: {token}")


def _strict_pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k, v in pairs:
        _require_utf8_no_surrogate(k, "json_key")
        oracle.require_ok(k not in out, "DUPLICATE_JSON_KEY", f"duplicate key {k!r}")
        out[k] = v
    return out


def parse_strict_json_text(text: Any, *, require_object: bool = True, max_bytes: int = MAX_JSON_BYTES) -> Any:
    s = _require_utf8_no_surrogate(text, "json_text")
    oracle.require_exact_bool(require_object, "require_object")
    oracle.require_exact_int(max_bytes, "max_bytes", lo=1, hi=MAX_ALLOWED_CAP_BYTES)
    oracle.require_ok(len(s.encode("utf-8")) <= max_bytes, "JSON_TOO_LARGE", "JSON text exceeds byte cap")
    try:
        parsed = json.loads(
            s,
            parse_float=_reject_float,
            parse_constant=_reject_constant,
            parse_int=_parse_bounded_int,
            object_pairs_hook=_strict_pairs_hook,
        )
    except oracle.PointChargesN3Error:
        raise
    except (json.JSONDecodeError, RecursionError, ValueError) as exc:
        raise oracle.PointChargesN3Error("JSON_PARSE_ERROR", f"invalid JSON ({type(exc).__name__})") from exc
    _validate_json_tree(parsed, "parsed_json")
    if require_object and type(parsed) is not dict:
        raise oracle.PointChargesN3Error("JSON_ROOT_NOT_OBJECT", "JSON root must be an object")
    return parsed


def canonical_json_bytes(obj: Any) -> bytes:
    _validate_json_tree(obj, "canonical_json")
    return json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def _sha256_canonical(obj: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(obj)).hexdigest()


# ============================================================================
# 3. Deep Schema-Safe Secondary Auditor (Full Enforcement Under full_replay=False)
# ============================================================================

def _req_dict(val: Any, label: str, expected_keys: frozenset[str] | None = None) -> dict[str, Any]:
    oracle.require_ok(type(val) is dict, "INVALID_SCHEMA_DICT", f"{label}: exact dict required")
    if expected_keys is not None:
        oracle.require_ok(frozenset(val.keys()) == expected_keys, "INVALID_SCHEMA_KEYS", f"{label}: key set mismatch")
    return val


def _req_list(val: Any, label: str, exact_len: int) -> list[Any]:
    oracle.require_ok(type(val) is list and len(val) == exact_len, "INVALID_SCHEMA_LIST", f"{label}: list len {exact_len} required")
    return val


def _req_int_f(obj: dict[str, Any], key: str, label: str) -> int:
    oracle.require_ok(key in obj, "MISSING_FIELD", f"{label}: missing {key!r}")
    return oracle.require_exact_int(obj[key], f"{label}.{key}")


def _req_bool_f(obj: dict[str, Any], key: str, label: str) -> bool:
    oracle.require_ok(key in obj, "MISSING_FIELD", f"{label}: missing {key!r}")
    return oracle.require_exact_bool(obj[key], f"{label}.{key}")


_P1_EXPECTED_CASES = (
    ("equilateral_center", Fraction(1), Fraction(1), Fraction(1), Fraction(1, 3), Fraction(1, 3)),
    ("equilateral_off_center", Fraction(1), Fraction(1), Fraction(1), Fraction(2, 5), Fraction(3, 10)),
    ("right_345_incenter", Fraction(9), Fraction(16), Fraction(25), Fraction(1, 4), Fraction(1, 3)),
    ("right_345_off_incenter", Fraction(9), Fraction(16), Fraction(25), Fraction(1, 5), Fraction(3, 7)),
    ("heronian_13_14_15_incenter", Fraction(169), Fraction(196), Fraction(225), Fraction(13, 42), Fraction(14, 42)),
    ("scalene_3_5_7_arbitrary", Fraction(9), Fraction(25), Fraction(49), Fraction(2, 9), Fraction(4, 9)),
    ("flat_obtuse_arbitrary", Fraction(10001**2), Fraction(10001**2), Fraction(20000**2), Fraction(1, 4), Fraction(1, 3)),
)
_P2_EXPECTED_SPECS = (
    ("iso_acute_1_3_1_3", Fraction(1, 3), Fraction(1, 3)),
    ("scalene_1_4_1_5", Fraction(1, 4), Fraction(1, 5)),
    ("flat_iso_1_20_1_20", Fraction(1, 20), Fraction(1, 20)),
    ("ultra_flat_iso_1_10000", Fraction(1, 10_000), Fraction(1, 10_000)),
    ("extreme_flat_iso_1_100000", Fraction(1, 100_000), Fraction(1, 100_000)),
    ("scalene_needle_1_10000_1_3", Fraction(1, 10_000), Fraction(1, 3)),
)
_P4_EXPECTED_BOXES = (
    (
        "equilateral_4d_compact",
        ("255/256", "257/256"),
        ("255/256", "257/256"),
        "1",
        ("125/384", "131/384"),
        ("125/384", "131/384"),
        "0.000576",
        "0.753628",
        "0.594345",
    ),
    (
        "right_345_4d_compact",
        ("575/64", "577/64"),
        ("1023/64", "1025/64"),
        "25",
        ("31/128", "33/128"),
        ("125/384", "131/384"),
        "2.782278",
        "0.772429",
        "0.670078",
    ),
)
_P5_EXPECTED_CONFIGS = (
    ("n3_weierstrass_iso_acute", Fraction(1, 3), Fraction(1, 3), None),
    ("n3_weierstrass_scalene", Fraction(1, 4), Fraction(1, 5), None),
    ("n3_weierstrass_flat_iso", Fraction(1, 20), Fraction(1, 20), None),
    ("n3_weierstrass_skew_2eq", Fraction(1, 3), Fraction(1, 3), (Fraction(20), Fraction(1), Fraction(1))),
)


def secondary_audit_payload(payload: Any) -> dict[str, Any]:
    """Schema-safe fail-closed exact algebraic & witness cross-check across ALL 5 pillars."""
    pl = _req_dict(payload, "payload", _PAYLOAD_KEYS)
    p1 = _req_dict(pl["pillar1_algebraic_identities"], "p1", _P1_KEYS)
    p2 = _req_dict(pl["pillar2_incenter_universal_theorem"], "p2", _P2_KEYS)
    p3 = _req_dict(pl["pillar3_isosceles_axis_uniqueness"], "p3", _P3_KEYS)
    p4 = _req_dict(pl["pillar4_gale_nikaido_compact_univalence"], "p4", _P4_KEYS)
    p5 = _req_dict(pl["pillar5_krawczyk_configurations"], "p5", _P5_KEYS)
    summ = _req_dict(pl["summary"], "summary", _SUMMARY_KEYS)

    # P1: 7 algebraic identity cases verified live via eval_exact_pmin_and_dl
    p1_rows = _req_list(p1["algebraic_identity_cases"], "p1.cases", 7)
    p1_in_rmin: list[bool] = []
    p1_rows_ok = True
    for i, (row, (exp_name, a2, b2, c2, u1, u2)) in enumerate(zip(p1_rows, _P1_EXPECTED_CASES)):
        rd = _req_dict(row, f"p1.cases[{i}]", _P1_ROW_KEYS)
        ev = oracle.eval_exact_pmin_and_dl(a2, b2, c2, u1, u2)
        det_h_q = parse_exact_rational(rd["det_h_tilde_q"], f"p1[{i}].detH")
        det_dl_q = parse_exact_rational(rd["det_dl_q"], f"p1[{i}].detDL")
        in_rmin = _req_bool_f(rd, "in_R_min", f"p1[{i}]")
        p1_in_rmin.append(in_rmin)
        if not all(
            (
                rd["name"] == exp_name,
                det_h_q == ev["det_h_poly"],
                rd["det_h_tilde_dec6"] == oracle.format_decimal(det_h_q, 6),
                det_dl_q == ev["det_dl"],
                in_rmin is (ev["pmin"] > 0),
                _req_bool_f(rd, "phase_formula_matches_pmin", f"p1[{i}]"),
                _req_bool_f(rd, "dl_jacobian_identity_exact", f"p1[{i}]"),
                _req_bool_f(rd, "verified", f"p1[{i}]"),
            )
        ):
            p1_rows_ok = False
    if p1_in_rmin != [True, True, True, True, True, True, False]:
        p1_rows_ok = False

    # P2: 6 Weierstrass incenter witnesses verified live via weierstrass_rational_triangle
    p2_rows = _req_list(p2["weierstrass_incenter_witnesses"], "p2.witnesses", 6)
    p2_rows_ok = True
    for i, (row, (exp_name, exp_v1, exp_v2)) in enumerate(zip(p2_rows, _P2_EXPECTED_SPECS)):
        rd = _req_dict(row, f"p2[{i}]", _P2_ROW_KEYS)
        wt = oracle.weierstrass_rational_triangle(exp_v1, exp_v2)
        v1_q = parse_exact_rational(rd["v1_q"], f"p2[{i}].v1")
        v2_q = parse_exact_rational(rd["v2_q"], f"p2[{i}].v2")
        v3_q = parse_exact_rational(rd["v3_q"], f"p2[{i}].v3")
        det_h_q = parse_exact_rational(rd["det_h_tilde_q"], f"p2[{i}].detH")
        if not all(
            (
                rd["name"] == exp_name,
                v1_q == exp_v1,
                v2_q == exp_v2,
                v3_q == wt["v3"],
                rd["aspect_h_over_max_side_dec6"] == oracle.format_decimal(wt["aspect_h_over_max_side"], 6),
                det_h_q == Fraction(1, 4),
                _req_bool_f(rd, "grad_exact_zero", f"p2[{i}]"),
                _req_bool_f(rd, "z2_phase_sum_exact_zero", f"p2[{i}]"),
                _req_bool_f(rd, "h_tilde_exact_half_I2", f"p2[{i}]"),
                _req_bool_f(rd, "pmin_eq_one_eighth_rprod", f"p2[{i}]"),
                _req_bool_f(rd, "verified", f"p2[{i}]"),
            )
        ):
            p2_rows_ok = False

    # P3: 4 isosceles symmetry-axis cases verified against exact Pythagorean formulas
    p3_rows = _req_list(p3["isosceles_axis_cases"], "p3.cases", 4)
    p3_ks: list[int] = []
    p3_rows_ok = True
    for i, row in enumerate(p3_rows):
        rd = _req_dict(row, f"p3[{i}]", _P3_ROW_KEYS)
        kv = _req_int_f(rd, "k", f"p3[{i}]")
        p3_ks.append(kv)
        h_q = parse_exact_rational(rd["height_h_q"], f"p3[{i}].h")
        y_f = parse_exact_rational(rd["y_fold_q"], f"p3[{i}].yf")
        y_p = parse_exact_rational(rd["y_pitchfork_q"], f"p3[{i}].yp")
        exp_h = Fraction(3 * kv, kv * kv - 2) if kv >= 2 else Fraction(0)
        exp_yf = Fraction(1, kv) if kv >= 2 else Fraction(0)
        if not all(
            (
                kv >= 2,
                h_q == exp_h,
                y_f == exp_yf,
                y_p == 2 * y_f,
                _req_bool_f(rd, "y_pitchfork_eq_2_y_fold", f"p3[{i}]"),
                _req_bool_f(rd, "vyy_fold_exact_zero", f"p3[{i}]"),
                _req_bool_f(rd, "vxx_pitchfork_exact_zero", f"p3[{i}]"),
                _req_bool_f(rd, "q_prime_strictly_neg_at_mid", f"p3[{i}]"),
                _req_bool_f(rd, "vxx_strictly_pos_at_mid", f"p3[{i}]"),
                _req_bool_f(rd, "incenter_strictly_between_fold_and_pitchfork", f"p3[{i}]"),
                _req_bool_f(rd, "verified", f"p3[{i}]"),
            )
        ):
            p3_rows_ok = False

    # P4: 2 compact 4D boxes x 64 sub-boxes = 128 sub-boxes with exact range & certified bound match
    p4_rows = _req_list(p4["compact_4d_boxes"], "p4.boxes", 2)
    p4_sub_total = 0
    p4_rows_ok = True
    for i, (row, spec) in enumerate(zip(p4_rows, _P4_EXPECTED_BOXES)):
        rd = _req_dict(row, f"p4[{i}]", _P4_ROW_KEYS)
        exp_name, exp_a2, exp_b2, exp_c2, exp_u1, exp_u2, exp_pm, exp_s00, exp_ds = spec
        parse_exact_rational(rd["c2_q"], f"p4[{i}].c2")
        for rng_k in ("a2_range_q", "b2_range_q", "u1_range_q", "u2_range_q"):
            rng_l = _req_list(rd[rng_k], f"p4[{i}].{rng_k}", 2)
            parse_exact_rational(rng_l[0], f"p4[{i}].{rng_k}[0]")
            parse_exact_rational(rng_l[1], f"p4[{i}].{rng_k}[1]")
        chk = _req_int_f(rd, "subboxes_checked", f"p4[{i}]")
        pas = _req_int_f(rd, "subboxes_passed", f"p4[{i}]")
        p4_sub_total += pas
        if not all(
            (
                rd["name"] == exp_name,
                tuple(rd["a2_range_q"]) == exp_a2,
                tuple(rd["b2_range_q"]) == exp_b2,
                rd["c2_q"] == exp_c2,
                tuple(rd["u1_range_q"]) == exp_u1,
                tuple(rd["u2_range_q"]) == exp_u2,
                rd["min_pmin_lo_dec6"] == exp_pm,
                rd["min_s00_lo_dec6"] == exp_s00,
                rd["min_det_sym_lo_dec6"] == exp_ds,
                chk == 64,
                pas == 64,
                _req_bool_f(rd, "strictly_injective_on_4d_box", f"p4[{i}]"),
            )
        ):
            p4_rows_ok = False

    # P5: 4 configurations (4, 4, 4, 2 -> 14 equilibria) verified LIVE via 192-bit Krawczyk & exact box geometry
    p5_rows = _req_list(p5["configurations"], "p5.configs", 4)
    p5_counts: list[int] = []
    p5_rows_ok = True
    exp_hw = Fraction(1, 1 << 20)
    for i, (row, (exp_name, v1, v2, custom_q)) in enumerate(zip(p5_rows, _P5_EXPECTED_CONFIGS)):
        rd = _req_dict(row, f"p5[{i}]", _P5_CFG_KEYS)
        wt = oracle.weierstrass_rational_triangle(v1, v2)
        exp_sites = wt["sites"]
        exp_charges = wt["charges"] if custom_q is None else custom_q
        ch_list = _req_list(rd["charges_q"], f"p5[{i}].charges", 3)
        parsed_charges = tuple(parse_exact_rational(tok, f"p5[{i}].charges[{ci}]") for ci, tok in enumerate(ch_list))
        eq_n = _req_int_f(rd, "equilibria_certified", f"p5[{i}]")
        m1 = _req_int_f(rd, "m1_local_minima", f"p5[{i}]")
        m2 = _req_int_f(rd, "m2_saddles", f"p5[{i}]")
        diff_m = _req_int_f(rd, "m2_minus_m1", f"p5[{i}]")
        disj_n = _req_int_f(rd, "pairwise_disjoint_pairs", f"p5[{i}]")
        p5_counts.append(eq_n)
        eqs = _req_list(rd["equilibria"], f"p5[{i}].eqs", eq_n)
        eqs_ok = True
        obs_m1 = 0
        obs_m2 = 0
        box_bounds: list[tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]] = []
        for j, eq in enumerate(eqs):
            ed = _req_dict(eq, f"p5[{i}].eq[{j}]", _P5_EQ_KEYS)
            ctr_l = _req_list(ed["center_q"], f"p5[{i}].eq[{j}].center", 2)
            lo_l = _req_list(ed["box_lo_q"], f"p5[{i}].eq[{j}].lo", 2)
            hi_l = _req_list(ed["box_hi_q"], f"p5[{i}].eq[{j}].hi", 2)
            cx = parse_exact_rational(ctr_l[0], f"p5[{i}].eq[{j}].cx")
            cy = parse_exact_rational(ctr_l[1], f"p5[{i}].eq[{j}].cy")
            hw = parse_exact_rational(ed["half_width_q"], f"p5[{i}].eq[{j}].hw")
            lx = parse_exact_rational(lo_l[0], f"p5[{i}].eq[{j}].lx")
            ly = parse_exact_rational(lo_l[1], f"p5[{i}].eq[{j}].ly")
            ux = parse_exact_rational(hi_l[0], f"p5[{i}].eq[{j}].ux")
            uy = parse_exact_rational(hi_l[1], f"p5[{i}].eq[{j}].uy")
            box_bounds.append(((lx, ly), (ux, uy)))
            idx = _req_int_f(ed, "morse_index_3d", f"p5[{i}].eq[{j}]")
            sign_h = _req_int_f(ed, "det_h_sign", f"p5[{i}].eq[{j}]")
            if idx == 1:
                obs_m1 += 1
            elif idx == 2:
                obs_m2 += 1
            else:
                eqs_ok = False
            live_eq = oracle.certify_2d_krawczyk_box(exp_sites, exp_charges, cx, cy, hw)
            if not all(
                (
                    hw == exp_hw,
                    lx == cx - hw,
                    ly == cy - hw,
                    ux == cx + hw,
                    uy == cy + hw,
                    ed["krawczyk_rho_ub_dec6"] == live_eq["krawczyk_rho_ub_dec6"],
                    idx == live_eq["morse_index_3d"],
                    sign_h == live_eq["det_h_sign"],
                    _req_bool_f(ed, "preconditioner_det_nonzero", f"p5[{i}].eq[{j}]"),
                    _req_bool_f(ed, "krawczyk_contraction_lt_1", f"p5[{i}].eq[{j}]"),
                    _req_bool_f(ed, "krawczyk_strict_interior", f"p5[{i}].eq[{j}]"),
                    _req_bool_f(ed, "trace_strictly_positive", f"p5[{i}].eq[{j}]"),
                )
            ):
                eqs_ok = False
        disjoint_ok = True
        for a_i in range(len(box_bounds)):
            for b_i in range(a_i + 1, len(box_bounds)):
                if not oracle._boxes_strictly_disjoint(box_bounds[a_i][0], box_bounds[a_i][1], box_bounds[b_i][0], box_bounds[b_i][1]):
                    disjoint_ok = False
        if not all(
            (
                rd["name"] == exp_name,
                parsed_charges == exp_charges,
                m1 == obs_m1,
                m2 == obs_m2,
                m1 + m2 == eq_n,
                m2 - m1 == 2,
                diff_m == 2,
                disj_n == eq_n * (eq_n - 1) // 2,
                disjoint_ok,
                _req_bool_f(rd, "all_boxes_pairwise_disjoint", f"p5[{i}]"),
                _req_bool_f(rd, "verified", f"p5[{i}]"),
                eqs_ok,
            )
        ):
            p5_rows_ok = False

    vz_probe_q = parse_exact_rational(p1.get("out_of_plane_vz_probe_q"), "p1.vz_probe")
    ok_all = all(
        (
            pl["track_id"] == oracle.TRACK_ID,
            pl["issue_ref"] == "Foundation-0/point-charges#4",
            pl["novelty_class"] == "N1",
            _req_bool_f(p1, "all_verified", "p1") is True,
            vz_probe_q == Fraction(-12, 125),
            _req_bool_f(p1, "out_of_plane_z_restoring_nonzero", "p1") is True,
            _req_bool_f(p1, "morse_relation_m2_minus_m1_eq_2", "p1") is True,
            p1_rows_ok,
            _req_bool_f(p2, "all_verified", "p2") is True,
            _req_bool_f(p2, "closes_487_flat_triangle_nonemptiness_gap", "p2") is True,
            p2_rows_ok,
            _req_bool_f(p3, "all_verified", "p3") is True,
            _req_bool_f(p3, "at_most_one_on_axis_minimum_all_h_gt_0", "p3") is True,
            p3_ks == [2, 3, 10, 100],
            p3_rows_ok,
            _req_bool_f(p4, "all_verified", "p4") is True,
            _req_int_f(p4, "total_4d_subboxes_certified", "p4") == 128,
            p4_sub_total == 128,
            p4_rows_ok,
            _req_bool_f(p5, "all_verified", "p5") is True,
            _req_int_f(p5, "total_equilibria_certified", "p5") == 14,
            p5_counts == [4, 4, 4, 2],
            p5_rows_ok,
            _req_int_f(summ, "pillar_count", "summary") == 5,
            _req_bool_f(summ, "all_five_pillars_verified", "summary") is True,
            _req_int_f(summ, "total_4d_univalence_subboxes", "summary") == 128,
            _req_int_f(summ, "total_krawczyk_equilibria", "summary") == 14,
            _req_bool_f(summ, "global_n3_conjecture_open_outside_compact_boxes", "summary") is True,
            _req_int_f(summ, "float_count_on_trust_path", "summary") == 0,
            _req_int_f(summ, "bare_assert_count", "summary") == 0,
            _req_int_f(summ, "rng_calls", "summary") == 0,
        )
    )
    oracle.require_ok(ok_all, "SECONDARY_AUDIT_FAILED", "secondary algebraic cross-check failed")
    return {"secondary_audit_passed": True, "pillars_audited": 5}


def build_standalone_certificate() -> dict[str, Any]:
    payload = oracle.build_payload()
    payload_sha = _sha256_canonical(payload)
    outer = {
        "schema_version": oracle.SCHEMA_VERSION,
        "track_id": oracle.TRACK_ID,
        "payload_sha256": payload_sha,
        "payload": payload,
    }
    cert_sha = _sha256_canonical(outer)
    return {
        "schema_version": oracle.SCHEMA_VERSION,
        "track_id": oracle.TRACK_ID,
        "payload_sha256": payload_sha,
        "certificate_sha256": cert_sha,
        "payload": payload,
    }


def verify_certificate_dict(cert: Any, *, full_replay: bool = True) -> dict[str, Any]:
    oracle.require_exact_bool(full_replay, "full_replay")
    oracle.require_ok(type(cert) is dict, "INVALID_CERT_TYPE", "certificate must be exact dict")
    _validate_json_tree(cert, "cert")
    oracle.require_ok(frozenset(cert.keys()) == _CERT_KEYS, "INVALID_CERT_KEYS", "certificate keys mismatch")
    oracle.require_ok(cert["schema_version"] == oracle.SCHEMA_VERSION, "INVALID_SCHEMA", "schema_version mismatch")
    oracle.require_ok(cert["track_id"] == oracle.TRACK_ID, "INVALID_TRACK_ID", "track_id mismatch")
    oracle.require_ok(type(cert["payload"]) is dict, "INVALID_PAYLOAD_TYPE", "payload must be exact dict")

    claimed_p = require_hex64_lower_ascii(cert["payload_sha256"], "payload_sha256")
    claimed_c = require_hex64_lower_ascii(cert["certificate_sha256"], "certificate_sha256")

    actual_p = _sha256_canonical(cert["payload"])
    oracle.require_ok(hmac.compare_digest(claimed_p, actual_p), "PAYLOAD_SHA256_MISMATCH", "payload_sha256 mismatch")

    outer = {
        "schema_version": cert["schema_version"],
        "track_id": cert["track_id"],
        "payload_sha256": claimed_p,
        "payload": cert["payload"],
    }
    actual_c = _sha256_canonical(outer)
    oracle.require_ok(hmac.compare_digest(claimed_c, actual_c), "CERT_SHA256_MISMATCH", "certificate_sha256 mismatch")

    sec = secondary_audit_payload(cert["payload"])
    if full_replay:
        ref = build_standalone_certificate()
        oracle.require_ok(hmac.compare_digest(actual_p, ref["payload_sha256"]), "REPLAY_MISMATCH", "replay mismatch")

    return {
        "valid": True,
        "schema_version": oracle.SCHEMA_VERSION,
        "track_id": oracle.TRACK_ID,
        "payload_sha256": actual_p,
        "certificate_sha256": actual_c,
        "full_replay": full_replay,
        "secondary_audit": sec,
    }


def check_certificate_file(path: Path, *, full_replay: bool = True) -> dict[str, Any]:
    raw = read_regular_bytes_capped(path, "certificate", MAX_JSON_BYTES)
    oracle.require_ok(b"\r" not in raw, "CRLF_FORBIDDEN", "CRLF line endings forbidden")
    oracle.require_ok(raw.endswith(b"\n"), "MISSING_TRAILING_LF", "certificate must end with LF")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise oracle.PointChargesN3Error("INVALID_UTF8", "certificate is not valid UTF-8") from exc
    cert = parse_strict_json_text(text, require_object=True)
    canon = canonical_json_bytes(cert)
    oracle.require_ok(raw == canon + b"\n", "NON_CANONICAL_BYTES", "certificate is not canonical ASCII JSON + LF")
    res = verify_certificate_dict(cert, full_replay=full_replay)
    res["canonical_sha256"] = hashlib.sha256(canon).hexdigest()
    res["file_sha256"] = hashlib.sha256(raw).hexdigest()
    return res


def write_certificate_file(path: Path) -> dict[str, Any]:
    cert = build_standalone_certificate()
    canon = canonical_json_bytes(cert)
    raw = canon + b"\n"
    write_regular_bytes_safe(path, raw, "certificate")
    return {
        "written": True,
        "schema_version": oracle.SCHEMA_VERSION,
        "track_id": oracle.TRACK_ID,
        "payload_sha256": cert["payload_sha256"],
        "certificate_sha256": cert["certificate_sha256"],
        "file_sha256": hashlib.sha256(raw).hexdigest(),
    }


# ============================================================================
# 4. Strict Fail-Closed CLI (Space & Equals Flag-Like Rejection)
# ============================================================================

def _reject_bad_path_arg(tok: str, val: str) -> Path:
    oracle.require_ok(bool(val) and not val.startswith("-"), "INVALID_PATH_ARG", f"{tok} path empty or flag-like")
    return Path(val)


def _parse_cli_argv(argv: Sequence[str]) -> tuple[str, Path | None, bool]:
    oracle.require_ok(type(argv) in (list, tuple), "INVALID_ARGV", "argv must be list or tuple")
    for item in argv:
        _require_utf8_no_surrogate(item, "argv_item")

    mode: str | None = None
    target: Path | None = None
    as_json = False
    i = 0
    n = len(argv)
    while i < n:
        tok = argv[i]
        if tok == "--json":
            oracle.require_ok(not as_json, "DUPLICATE_FLAG", "duplicate --json flag")
            as_json = True
            i += 1
        elif tok == "--standalone":
            oracle.require_ok(mode is None, "CONFLICTING_MODE", "multiple mode flags forbidden")
            mode = "standalone"
            i += 1
        elif tok in ("--check", "--write"):
            oracle.require_ok(mode is None, "CONFLICTING_MODE", "multiple mode flags forbidden")
            oracle.require_ok(i + 1 < n, "MISSING_PATH_ARG", f"{tok} requires a path argument")
            target = _reject_bad_path_arg(tok, argv[i + 1])
            mode = "check" if tok == "--check" else "write"
            i += 2
        elif tok.startswith("--check=") or tok.startswith("--write="):
            oracle.require_ok(mode is None, "CONFLICTING_MODE", "multiple mode flags forbidden")
            flag_name, val = tok.split("=", 1)
            target = _reject_bad_path_arg(flag_name, val)
            mode = "check" if flag_name == "--check" else "write"
            i += 1
        else:
            raise oracle.PointChargesN3Error("UNKNOWN_FLAG", f"unrecognized CLI argument: {tok}")

    if mode is None:
        mode = "standalone"
    return mode, target, as_json


def main(argv: Sequence[str] | None = None) -> int:
    raw_argv: Any = sys.argv[1:] if argv is None else argv
    try:
        mode, target, as_json = _parse_cli_argv(raw_argv)
        if mode == "standalone":
            out: dict[str, Any] = build_standalone_certificate()
        elif mode == "write":
            out = write_certificate_file(target)  # type: ignore[arg-type]
        else:
            out = check_certificate_file(target)  # type: ignore[arg-type]
    except oracle.PointChargesN3Error as exc:
        err = {"valid": False, "error_code": exc.code, "error_message": exc.detail}
        sys.stderr.write(canonical_json_bytes(err).decode("ascii") + "\n")
        return 2

    if as_json:
        sys.stdout.write(canonical_json_bytes(out).decode("ascii") + "\n")
    else:
        sys.stdout.write(
            f"OK {oracle.SCHEMA_VERSION} payload_sha256={out['payload_sha256']} "
            f"certificate_sha256={out['certificate_sha256']}\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
