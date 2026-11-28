#!/usr/bin/env python3
"""
Pure-Init Lean 4.36.0 Exact-Rational Certificate Compiler & Verifier
for Foundation-0/point-charges (Refs #6).

Pure Python 3 standard library only (`int`, `fractions.Fraction`, `math.isqrt`, `math.lcm`).
Zero IEEE-754 floats, zero `ast.Div`, zero bare `assert`s, zero external Lean packages.

CLI modes:
  --compile [--offline] [--json]
      Validate `witness_pack.json`, emit `PointChargesExactCert.lean`, `AxiomAudit.lean`,
      `lakefile.toml`, `lake-manifest.json`, and `lean-toolchain`, run cold Lean 4.36.0 +
      `leanchecker --fresh` (unless `--offline`), and write self-bound
      `point_charges_lean_cert.json`.
  --check [<cert>] [--offline] [--json]
      Validate witness pack, byte-exact compiled Lean/Lake/manifest parity, static Lean
      trust-escape scan, live Lean 4.36.0 + `leanchecker --fresh` (unless `--offline`),
      and canonical self-bound certificate integrity.
  --standalone
      Validate witness pack and on-disk Lean/Lake/manifest parity, then print the
      canonical certificate JSON to stdout.
"""
from __future__ import annotations

import ast
from fractions import Fraction
import hashlib
import hmac
import json
import math
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
from typing import Any, Sequence

SCHEMA_VERSION = "point_charges.formal_pure_init.lean_cert.v1"
WITNESS_SCHEMA = "point_charges.formal_pure_init.witness_pack.v1"
UPSTREAM_ISSUE = "Foundation-0/point-charges#6"
EXPECTED_LEAN_TOOLCHAIN = "leanprover/lean4:v4.36.0"
LEAN_MODULE_NAME = "PointChargesExactCert"
LEAN_FILENAME = "PointChargesExactCert.lean"
AXIOM_AUDIT_FILENAME = "AxiomAudit.lean"
LAKEFILE_FILENAME = "lakefile.toml"
MANIFEST_FILENAME = "lake-manifest.json"
TOOLCHAIN_FILENAME = "lean-toolchain"
WITNESS_PACK_FILENAME = "witness_pack.json"
DEFAULT_CERT_FILENAME = "point_charges_lean_cert.json"

MAX_FILE_BYTES = 1024 * 1024
MAX_INT_DIGITS = 120
_HEX64_RE = re.compile(r"\A[0-9a-f]{64}\Z")
_CANON_Q_RE = re.compile(r"\A(?:0|-?[1-9][0-9]*)(?:/[1-9][0-9]*)?\Z")
_ALLOWED_LEAN_AXIOMS: frozenset[str] = frozenset({"propext", "Quot.sound"})
_AXIOM_NONE_RE = re.compile(r"\A'([A-Za-z0-9_.]+)' does not depend on any axioms\Z")
_AXIOM_DEP_RE = re.compile(r"\A'([A-Za-z0-9_.]+)' depends on axioms: \[([^\]]*)\]\Z")
_HEX_DIGITS: frozenset[str] = frozenset("0123456789abcdefABCDEF")

_ALLOWED_MATH_ATTRS: frozenset[str] = frozenset({"isqrt", "lcm", "gcd"})
_FORBIDDEN_IMPORTS: frozenset[str] = frozenset(
    {"builtins", "cmath", "decimal", "numpy", "scipy", "statistics", "torch", "operator", "importlib"}
)
_FORBIDDEN_BUILTINS: frozenset[str] = frozenset(
    {
        "float",
        "complex",
        "pow",
        "getattr",
        "setattr",
        "delattr",
        "eval",
        "exec",
        "compile",
        "__import__",
        "vars",
        "globals",
        "locals",
        "dir",
        "builtins",
        "__builtins__",
    }
)
_FORBIDDEN_ATTRS: frozenset[str] = frozenset(
    {
        "from_float",
        "from_decimal",
        "fromhex",
        "__truediv__",
        "__rtruediv__",
        "__itruediv__",
        "__pow__",
        "__rpow__",
        "__ipow__",
        "__float__",
        "__complex__",
        "__dict__",
        "__builtins__",
        "__globals__",
        "__subclasses__",
        "__getattribute__",
        "modules",
        "pow",
        "float",
        "complex",
    }
)

_LEAN4_BANNED_CORE_RE = re.compile(
    r"(?:"
    r"[«»]"
    r"|\b(?:sorry|admit|axiom|unsafe|partial|opaque|noncomputable|prelude"
    r"|native_decide|bv_decide|bv_check|decide_cbv|Classical"
    r"|sorryAx|ofReduceBool|ofReduceNat|trustCompiler|mkSorry"
    r"|extern|implemented_by|csimp|initialize|builtin_initialize"
    r"|attribute|instance"
    r"|macro|macro_rules|elab|elab_rules|syntax|declare_syntax_cat"
    r"|notation|prefix|infix|infixl|infixr|postfix|run_cmd|run_elab|run_meta"
    r"|throwError|interpolatedStr|set_option"
    r"|dbg_trace|dbg_sleep|dbg_stack_trace"
    r"|dbgTrace|dbgTraceVal|dbgTraceValId|dbgTraceIfShared|dbgSleep|dbgStackTrace"
    r"|panicCore|panicWithPos|panicWithPosWithDecl)\b"
    r"|\b(?:panic|unreachable|assert|debug_assert)!"
    r"|\b(?:debug|warn)\s*\.\s*(?:skipKernelTC|byAsSorry|proofAsSorry|sorryAsError|terminalTacticsAsSorry|sorry)\b"
    r"|@\s*\[\s*(?:(?:local|scoped)\s+)?(?:[^\]\n]*,\s*)?(?:builtin_)?(?:init|export|term_elab|command_elab|tactic|macro(?:_rules)?)\b"
    r"|#\s*(?:eval|reduce|exit|check|synth|guard(?:_msgs)?)\b"
    r"|#\s*print\s+(?!axioms\s+PointChargesExactCert\.[A-Za-z0-9_]+\b)"
    r"|(?<!\w)\+\s*native\b"
    r"|\bnative\s*:=(?!\s*(?:Bool\.)?false\b(?:\s*[},)]|\Z))"
    r"|\bDecide(?:Config|\.Config)\b"
    r"|\bconfig\s*:=(?!\s*\{\s*native\s*:=\s*(?:Bool\.)?false\s*\})"
    r")"
)
_LEAN4_IMPORT_RE = re.compile(r"\bimport\s+(?:all\s+)?([A-Za-z0-9_.]+)")

EXPECTED_THEOREM_NAMES: tuple[str, ...] = (
    "PointChargesExactCert.t01_coulomb_single_hessian_trace_zero",
    "PointChargesExactCert.t02_earnshaw_total_hessian_trace_zero",
    "PointChargesExactCert.t03_earnshaw_no_strict_definite_quadform",
    "PointChargesExactCert.t04_planar_offplane_force_and_z_curvature",
    "PointChargesExactCert.t05_lagrange_congruence_quadform_identity",
    "PointChargesExactCert.t06_lagrange_diagonal_morse_index2_witness",
    "PointChargesExactCert.t07_lagrange_diagonal_morse_index1_witness",
    "PointChargesExactCert.t08_disjoint_boxes_no_common_point",
    "PointChargesExactCert.t09_inv_sqrt_bracket_order_preservation",
    "PointChargesExactCert.t10_stereographic_s2q_unit_sphere",
    "PointChargesExactCert.t11_stereographic_s2q_chord_distance_sq",
    "PointChargesExactCert.t12_cayley_s3q_unit_norm_and_equatorial_bridge",
    "PointChargesExactCert.t13_n3_isosceles_cubic_cancellation_and_pythagorean_roots",
    "PointChargesExactCert.t14_regular_tetrahedron_s2q_sqrt3_obstruction",
    "PointChargesExactCert.t15_compiled_n2_analytic_equilibrium_and_hessian",
    "PointChargesExactCert.t16_compiled_n3_weierstrass_incenter_iso_acute",
    "PointChargesExactCert.t17_compiled_n3_weierstrass_incenter_scalene",
    "PointChargesExactCert.t18_compiled_n3_pythagorean_off_plane_confinement",
    "PointChargesExactCert.t19_compiled_n4_eleven_disjointness_and_morse",
    "PointChargesExactCert.t20_compiled_n8_fiftyone_disjointness_and_morse",
    "PointChargesExactCert.t21_compiled_thomson_s2q_energy_dominance",
    "PointChargesExactCert.t22_compiled_stereographic_s2q_and_cayley_s3q_bridges",
)

EXPECTED_LAKEFILE_TEXT = (
    'name = "PointChargesExactCert"\n'
    'version = "0.1.0"\n'
    'defaultTargets = ["PointChargesExactCert"]\n\n'
    "[[lean_lib]]\n"
    'name = "PointChargesExactCert"\n'
)
EXPECTED_MANIFEST_TEXT = (
    '{"lakeDir":".lake","name":"PointChargesExactCert",'
    '"packages":[],"packagesDir":".lake/packages","version":"1.1.0"}\n'
)
EXPECTED_TOOLCHAIN_TEXT = f"{EXPECTED_LEAN_TOOLCHAIN}\n"


class CertError(ValueError):
    """Typed fail-closed error for the pure-Init Lean 4.36.0 certificate compiler."""

    def __init__(self, code: str, detail: str) -> None:
        self.code = code
        self.detail = detail
        super().__init__(f"{code}: {detail}")


def require_ok(cond: bool, code: str, detail: str) -> None:
    if not cond:
        raise CertError(code, detail)


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fmt_q(val: Fraction) -> str:
    return str(val.numerator) if val.denominator == 1 else f"{val.numerator}/{val.denominator}"


def _q_div(a: Fraction, b: Fraction) -> Fraction:
    require_ok(b.numerator != 0, "ZERO_DIVISOR", "exact rational division by zero")
    return Fraction(a.numerator * b.denominator, a.denominator * b.numerator)


def _q_pow_pos(base: Fraction, exp: int) -> Fraction:
    require_ok(exp >= 1, "INVALID_EXP", "positive integer exponent required")
    num_acc = 1
    den_acc = 1
    for _ in range(exp):
        num_acc *= base.numerator
        den_acc *= base.denominator
    return Fraction(num_acc, den_acc)


def canonical_json_bytes(obj: Any) -> bytes:
    return (
        json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False, ensure_ascii=True) + "\n"
    ).encode("ascii")


def _collect_pure_annotation_ids(node: ast.AST, out: set[int]) -> None:
    """Collect node ids strictly inside pure static type syntax (never descending into calls or walrus exprs)."""
    if isinstance(node, (ast.Name, ast.Attribute, ast.Subscript, ast.Tuple, ast.List, ast.BinOp, ast.Constant)):
        out.add(id(node))
        for child in ast.iter_child_nodes(node):
            _collect_pure_annotation_ids(child, out)


def verify_pure_q_ast(source: str, label: str = "source") -> list[str]:
    """
    Fail-closed two-pass AST purity audit for exact rational Python sources.

    Rejects: syntax errors, bare `assert`, IEEE-754 `float` / `complex` literals or constructors,
    `ast.Div` / `ast.FloorDiv` (`/` and `//`), negative or non-literal exponents in `ast.Pow` or
    `pow()`, float f-string specifiers (`:.3f`, `:.2e`, `:.1g`, `:.1%`), dunder numeric/reflection
    escapes (`__pow__`, `__float__`, `__truediv__`, `__dict__`, `fromhex`, `vars`, `globals`,
    `locals`, `getattr`, `builtins`), non-call aliasing of `Fraction`/`math`/`fractions` (including
    tuple unpacks, attribute RHS, and walrus expressions), `Fraction(str)` / starred / non-2-arg
    `Fraction(...)` / `Fraction.from_float` / `Fraction.from_decimal`, imports of floating-point
    modules, and any `math.*` attribute outside `{isqrt, lcm, gcd}`.
    """
    violations: list[str] = []
    try:
        tree = ast.parse(source, filename=label)
    except SyntaxError as exc:
        err_ln = exc.lineno if exc.lineno is not None else 0
        return [f"L{err_ln}: syntax error: {exc.msg}"]

    all_nodes = list(ast.walk(tree))
    math_aliases: set[str] = {"math"}
    frac_aliases: set[str] = {"Fraction"}
    frac_mod_aliases: set[str] = {"fractions"}
    call_func_ids: set[int] = set()
    attr_value_ids: set[int] = set()
    annotation_ids: set[int] = set()

    # Pass 1: collect complete module/constructor alias sets & parent roles (order-independent).
    for node in all_nodes:
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "math":
                    math_aliases.add(alias.asname if alias.asname else "math")
                if alias.name == "fractions":
                    frac_mod_aliases.add(alias.asname if alias.asname else "fractions")
        elif isinstance(node, ast.ImportFrom):
            mod_name = node.module if node.module else ""
            if mod_name.split(".", 1)[0] == "fractions":
                for alias in node.names:
                    if alias.name == "Fraction":
                        frac_aliases.add(alias.asname if alias.asname else "Fraction")
        elif isinstance(node, ast.Call):
            call_func_ids.add(id(node.func))
        elif isinstance(node, ast.Attribute):
            attr_value_ids.add(id(node.value))
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.returns is not None:
                _collect_pure_annotation_ids(node.returns, annotation_ids)
        elif isinstance(node, ast.arg):
            if node.annotation is not None:
                _collect_pure_annotation_ids(node.annotation, annotation_ids)
        elif isinstance(node, ast.AnnAssign):
            _collect_pure_annotation_ids(node.annotation, annotation_ids)

    # Pass 2: validate every AST node against complete alias sets.
    for node in all_nodes:
        ln = node.lineno if hasattr(node, "lineno") else 0
        if isinstance(node, ast.Assert):
            violations.append(f"L{ln}: bare assert forbidden")
        elif isinstance(node, ast.Constant):
            if type(node.value).__name__ in ("float", "complex"):
                violations.append(f"L{ln}: float/complex literal forbidden")
        elif isinstance(node, (ast.Div, ast.FloorDiv)):
            violations.append(f"L{ln}: division operator forbidden")
        elif isinstance(node, ast.FormattedValue):
            if node.format_spec is not None:
                for sub in ast.walk(node.format_spec):
                    if isinstance(sub, ast.Constant):
                        if type(sub.value) is str:
                            if any(ch in sub.value for ch in "eEfFgG%"):
                                violations.append(f"L{ln}: float format specifier forbidden")
        elif isinstance(node, ast.Name):
            if node.id in _FORBIDDEN_BUILTINS:
                violations.append(f"L{ln}: forbidden builtin/name reference {node.id}")
            elif node.id in frac_aliases:
                if id(node) not in call_func_ids and id(node) not in annotation_ids:
                    violations.append(f"L{ln}: non-call Fraction reference/alias forbidden")
            elif node.id in (math_aliases | frac_mod_aliases):
                if id(node) not in attr_value_ids:
                    violations.append(f"L{ln}: bare module reference/alias {node.id} forbidden")
        elif isinstance(node, ast.Attribute):
            if node.attr in _FORBIDDEN_ATTRS:
                violations.append(f"L{ln}: forbidden attribute .{node.attr}")
            elif isinstance(node.value, ast.Name):
                if node.value.id in math_aliases and node.attr not in _ALLOWED_MATH_ATTRS:
                    violations.append(f"L{ln}: forbidden math.{node.attr}")
                elif node.value.id in frac_mod_aliases and node.attr == "Fraction":
                    if id(node) not in call_func_ids and id(node) not in annotation_ids:
                        violations.append(f"L{ln}: non-call fractions.Fraction reference/alias forbidden")
        elif isinstance(node, ast.Import):
            for alias in node.names:
                root_mod = alias.name.split(".", 1)[0]
                if root_mod in _FORBIDDEN_IMPORTS:
                    violations.append(f"L{ln}: forbidden import {alias.name}")
        elif isinstance(node, ast.ImportFrom):
            if node.level > 0:
                violations.append(f"L{ln}: relative import forbidden")
            mod_name = node.module if node.module else ""
            root_mod = mod_name.split(".", 1)[0]
            if root_mod in _FORBIDDEN_IMPORTS:
                violations.append(f"L{ln}: forbidden from-import {mod_name}")
            if root_mod == "math":
                for alias in node.names:
                    if alias.name not in _ALLOWED_MATH_ATTRS:
                        violations.append(f"L{ln}: forbidden math import {alias.name}")
            if root_mod == "fractions":
                for alias in node.names:
                    if alias.name == "*":
                        violations.append(f"L{ln}: wildcard fractions import forbidden")
        elif isinstance(node, ast.BinOp):
            if isinstance(node.op, ast.Pow):
                rhs = node.right
                is_safe_const = isinstance(rhs, ast.Constant) and type(rhs.value) is int and rhs.value >= 0
                if not is_safe_const:
                    violations.append(f"L{ln}: unsafe Pow exponent")
        elif isinstance(node, ast.AugAssign):
            if isinstance(node.op, ast.Pow):
                violations.append(f"L{ln}: augmented Pow forbidden")
        elif isinstance(node, ast.Call):
            fn = node.func
            if isinstance(fn, ast.Name):
                if fn.id in frac_aliases:
                    if len(node.args) != 2 or bool(node.keywords) or any(isinstance(a, ast.Starred) for a in node.args):
                        violations.append(f"L{ln}: Fraction() requires exact 2 non-starred positional int args")
            elif isinstance(fn, ast.Attribute):
                if isinstance(fn.value, ast.Name):
                    if fn.value.id in frac_mod_aliases and fn.attr == "Fraction":
                        if len(node.args) != 2 or bool(node.keywords) or any(isinstance(a, ast.Starred) for a in node.args):
                            violations.append(f"L{ln}: fractions.Fraction() requires exact 2 non-starred positional int args")
    return violations


def parse_q(raw: Any, label: str) -> Fraction:
    require_ok(type(raw) is str, "INVALID_Q_TYPE", f"{label} must be exact str")
    require_ok(_CANON_Q_RE.match(raw) is not None, "INVALID_Q_FORMAT", f"{label} invalid rational: {raw!r}")
    if "/" not in raw:
        require_ok(len(raw.lstrip("-")) <= MAX_INT_DIGITS, "INT_TOO_LARGE", label)
        return Fraction(int(raw), 1)
    n_s, d_s = raw.split("/", 1)
    require_ok(len(n_s.lstrip("-")) <= MAX_INT_DIGITS, "INT_TOO_LARGE", f"{label}.num")
    require_ok(len(d_s) <= MAX_INT_DIGITS, "INT_TOO_LARGE", f"{label}.den")
    num, den = int(n_s), int(d_s)
    require_ok(den > 1, "NON_CANONICAL_Q", f"{label}: denominator must be > 1")
    fr = Fraction(num, den)
    require_ok(fr.denominator == den, "NON_REDUCED_Q", f"{label}: fraction not reduced")
    return fr


def require_exact_int(val: Any, label: str, *, min_val: int | None = None, max_val: int | None = None) -> int:
    require_ok(type(val) is int, "INVALID_INT", f"{label} must be exact int")
    if min_val is not None:
        require_ok(val >= min_val, "INT_OUT_OF_RANGE", f"{label} < {min_val}")
    if max_val is not None:
        require_ok(val <= max_val, "INT_OUT_OF_RANGE", f"{label} > {max_val}")
    return val


def require_hex64(val: Any, label: str) -> str:
    require_ok(type(val) is str, "INVALID_HEX64", f"{label} must be str")
    require_ok(_HEX64_RE.match(val) is not None, "INVALID_HEX64", f"{label} must be 64 lower-hex")
    return val


def _reject_float(tok: str) -> Any:
    raise CertError("FLOAT_FORBIDDEN", f"IEEE-754 float forbidden: {tok[:32]}")


def _reject_constant(tok: str) -> Any:
    raise CertError("JSON_CONSTANT_FORBIDDEN", f"JSON constant forbidden: {tok[:32]}")


def _no_dup_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k, v in pairs:
        require_ok(k not in out, "DUPLICATE_JSON_KEY", f"duplicate key {k!r}")
        out[k] = v
    return out


def parse_strict_json_bytes(raw: bytes, label: str) -> dict[str, Any]:
    require_ok(len(raw) > 0, "BAD_FILE_SIZE", f"{label} empty")
    require_ok(len(raw) <= MAX_FILE_BYTES, "BAD_FILE_SIZE", f"{label} exceeds size cap")
    require_ok(b"\x00" not in raw, "BAD_BYTES", f"{label} contains NUL")
    require_ok(b"\r" not in raw, "BAD_BYTES", f"{label} contains CR")
    try:
        text = raw.decode("ascii", "strict")
    except UnicodeDecodeError as exc:
        raise CertError("NON_ASCII_JSON", f"{label} must be strict ASCII") from exc
    try:
        obj = json.loads(
            text,
            parse_float=_reject_float,
            parse_constant=_reject_constant,
            object_pairs_hook=_no_dup_keys,
        )
    except CertError:
        raise
    except Exception as exc:
        raise CertError("INVALID_JSON", f"{label}: {exc}") from exc
    require_ok(type(obj) is dict, "JSON_NOT_OBJECT", f"{label} root must be JSON object")
    return obj


def read_regular_bytes(path: Path, label: str) -> bytes:
    try:
        st = os.lstat(path)
    except OSError as exc:
        raise CertError("IO_ERROR", f"{label} unreadable: {exc}") from exc
    require_ok(stat.S_ISREG(st.st_mode), "NON_REGULAR_FILE", f"{label} must be a regular non-symlink file")
    require_ok(st.st_size <= MAX_FILE_BYTES, "BAD_FILE_SIZE", f"{label} exceeds size cap")
    o_nofollow = os.O_NOFOLLOW if hasattr(os, "O_NOFOLLOW") else 0
    o_nonblock = os.O_NONBLOCK if hasattr(os, "O_NONBLOCK") else 0
    fd = os.open(path, os.O_RDONLY | o_nofollow | o_nonblock)
    try:
        fst = os.fstat(fd)
        require_ok(stat.S_ISREG(fst.st_mode), "NON_REGULAR_FILE", f"{label} fd not regular")
        chunks: list[bytes] = []
        total = 0
        while True:
            chunk = os.read(fd, 65536)
            if not chunk:
                break
            total += len(chunk)
            require_ok(total <= MAX_FILE_BYTES, "BAD_FILE_SIZE", f"{label} stream exceeds size cap")
            chunks.append(chunk)
        return b"".join(chunks)
    finally:
        os.close(fd)


def write_regular_atomic_0644(path: Path, raw: bytes) -> None:
    parent = path.parent
    require_ok(parent.is_dir(), "PARENT_NOT_DIR", f"parent missing: {parent}")
    fd, tmp_name = tempfile.mkstemp(prefix=".pc_lean_", dir=str(parent))
    tmp_p = Path(tmp_name)
    try:
        os.fchmod(fd, 0o644)
        view = memoryview(raw)
        written = 0
        while written < len(raw):
            n = os.write(fd, view[written:])
            written += n
        os.fsync(fd)
        os.close(fd)
        fd = -1
        os.replace(tmp_p, path)
    except Exception:
        if fd >= 0:
            os.close(fd)
        if tmp_p.exists():
            tmp_p.unlink()
        raise


def _parse_vec3(raw: Any, label: str) -> tuple[Fraction, Fraction, Fraction]:
    require_ok(type(raw) is list, "INVALID_VEC3", f"{label} must be list")
    require_ok(len(raw) == 3, "INVALID_VEC3", f"{label} must have length 3")
    return (parse_q(raw[0], f"{label}[0]"), parse_q(raw[1], f"{label}[1]"), parse_q(raw[2], f"{label}[2]"))


def _det3_q(m: Sequence[Sequence[Fraction]]) -> Fraction:
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def _floor_mul_q(q: Fraction, scale: int) -> int:
    return divmod(q.numerator * scale, q.denominator)[0]


def _ceil_mul_q(q: Fraction, scale: int) -> int:
    return -divmod(-q.numerator * scale, q.denominator)[0]


def _lcm_many(values: Sequence[int]) -> int:
    acc = 1
    for v in values:
        acc = math.lcm(acc, abs(v))
    return acc


def validate_equilibrium_witness(w: Any) -> dict[str, Any]:
    """Validate an exact-rational equilibrium & Earnshaw Hessian witness in pure Fraction arithmetic."""
    require_ok(type(w) is dict, "INVALID_WITNESS", "equilibrium witness must be dict")
    name = w.get("name")
    require_ok(type(name) is str, "INVALID_WITNESS", "equilibrium witness name must be str")
    require_ok(bool(name.strip()), "INVALID_WITNESS", "missing equilibrium witness name")
    sites_raw = w.get("sites")
    charges_raw = w.get("charges")
    radii_raw = w.get("radii")
    require_ok(type(sites_raw) is list, "INVALID_WITNESS", f"{name}: sites must be list")
    require_ok(type(charges_raw) is list, "INVALID_WITNESS", f"{name}: charges must be list")
    require_ok(type(radii_raw) is list, "INVALID_WITNESS", f"{name}: radii must be list")
    require_ok(len(sites_raw) >= 2, "INVALID_WITNESS", f"{name}: need >= 2 sites")
    require_ok(len(sites_raw) == len(charges_raw), "INVALID_WITNESS", f"{name}: charges length mismatch")
    require_ok(len(sites_raw) == len(radii_raw), "INVALID_WITNESS", f"{name}: radii length mismatch")

    sites = [_parse_vec3(s, f"{name}.sites[{i}]") for i, s in enumerate(sites_raw)]
    require_ok(len(set(sites)) == len(sites), "DUPLICATE_SITE", f"{name}: duplicate charge site")
    charges = [parse_q(c, f"{name}.charges[{i}]") for i, c in enumerate(charges_raw)]
    require_ok(all(q > 0 for q in charges), "NONPOSITIVE_CHARGE", f"{name}: charges must be strictly positive")
    point = _parse_vec3(w.get("point"), f"{name}.point")
    radii = [parse_q(r, f"{name}.radii[{i}]") for i, r in enumerate(radii_raw)]
    require_ok(all(r > 0 for r in radii), "NONPOSITIVE_RADIUS", f"{name}: radii must be strictly positive")
    expected_morse = require_exact_int(w.get("expected_morse_index"), f"{name}.expected_morse_index", min_val=1, max_val=2)

    disps: list[tuple[Fraction, Fraction, Fraction]] = []
    for i, (site, r) in enumerate(zip(sites, radii)):
        d = (point[0] - site[0], point[1] - site[1], point[2] - site[2])
        d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
        require_ok(d2 == r * r, "RADIUS_SQ_MISMATCH", f"{name}: ||x - a_{i}||^2 != r_{i}^2")
        disps.append(d)

    zero_q = Fraction(0, 1)
    grad = [zero_q, zero_q, zero_q]
    hess = [[zero_q, zero_q, zero_q] for _ in range(3)]
    for q, d, r in zip(charges, disps, radii):
        r2 = r * r
        r3 = r2 * r
        r5 = r3 * r2
        for a in range(3):
            grad[a] -= _q_div(q * d[a], r3)
            for b in range(3):
                hess[a][b] += _q_div(q * (3 * d[a] * d[b] - (r2 if a == b else zero_q)), r5)

    require_ok(grad == [zero_q, zero_q, zero_q], "NONZERO_GRADIENT", f"{name}: grad V != 0")
    tr = hess[0][0] + hess[1][1] + hess[2][2]
    det = _det3_q(hess)
    require_ok(det != 0, "DEGENERATE_HESSIAN", f"{name}: det(H) == 0")
    require_ok(hess[0][1] == 0, "NONDIAGONAL_WITNESS", f"{name}: H_xy != 0")
    require_ok(hess[0][2] == 0, "NONDIAGONAL_WITNESS", f"{name}: H_xz != 0")
    require_ok(hess[1][2] == 0, "NONDIAGONAL_WITNESS", f"{name}: H_yz != 0")
    diag = [hess[0][0], hess[1][1], hess[2][2]]
    n_pos = sum(1 for v in diag if v > 0)
    n_neg = sum(1 for v in diag if v < 0)
    require_ok(n_neg == expected_morse, "MORSE_MISMATCH", f"{name}: Morse index mismatch")

    d_pos = _lcm_many([c.denominator for d in disps for c in d] + [r.denominator for r in radii])
    int_disps = [(int(d[0] * d_pos), int(d[1] * d_pos), int(d[2] * d_pos)) for d in disps]
    int_radii = [int(r * d_pos) for r in radii]
    raw_weights_h = [_q_div(q, _q_pow_pos(r, 5)) for q, r in zip(charges, radii)]
    scaled_w_h = [_q_div(w_h, Fraction(d_pos * d_pos, 1)) for w_h in raw_weights_h]
    l_h = _lcm_many([sw.denominator for sw in scaled_w_h])
    int_weights_h = [int(sw * l_h) for sw in scaled_w_h]
    scaled_w_g = [_q_div(q, _q_pow_pos(r, 3) * Fraction(d_pos, 1)) for q, r in zip(charges, radii)]
    l_g = _lcm_many([sw.denominator for sw in scaled_w_g])
    int_weights_g = [int(sw * l_g) for sw in scaled_w_g]

    return {
        "name": name,
        "n_charges": len(sites),
        "hessian_diag_q": [fmt_q(v) for v in diag],
        "hessian_trace_q": fmt_q(tr),
        "hessian_det_q": fmt_q(det),
        "morse_index": n_neg,
        "inertia": [n_pos, n_neg, 0],
        "cleared_grad_scale": l_g,
        "cleared_hess_scale": l_h,
        "int_disps": [list(d) for d in int_disps],
        "int_radii": int_radii,
        "int_weights_g": int_weights_g,
        "int_weights_h": int_weights_h,
        "cleared_hess_diag": [int(v * l_h) for v in diag],
    }


def validate_off_plane_witness(w: Any) -> dict[str, Any]:
    require_ok(type(w) is dict, "INVALID_WITNESS", "off_plane witness must be dict")
    sites_raw = w.get("sites")
    charges_raw = w.get("charges")
    radii_raw = w.get("radii")
    require_ok(type(sites_raw) is list, "INVALID_WITNESS", "off_plane sites must be list")
    require_ok(type(charges_raw) is list, "INVALID_WITNESS", "off_plane charges must be list")
    require_ok(type(radii_raw) is list, "INVALID_WITNESS", "off_plane radii must be list")
    require_ok(len(sites_raw) >= 1, "INVALID_WITNESS", "off_plane requires >= 1 site")
    require_ok(len(sites_raw) == len(charges_raw), "INVALID_WITNESS", "off_plane charges length mismatch")
    require_ok(len(sites_raw) == len(radii_raw), "INVALID_WITNESS", "off_plane radii length mismatch")

    sites = [_parse_vec3(s, f"off_plane.sites[{i}]") for i, s in enumerate(sites_raw)]
    require_ok(all(s[2] == 0 for s in sites), "SITE_NOT_PLANAR", "off_plane sites must lie on z = 0")
    charges = [parse_q(c, f"off_plane.charges[{i}]") for i, c in enumerate(charges_raw)]
    require_ok(all(q > 0 for q in charges), "NONPOSITIVE_CHARGE", "off_plane charges must be > 0")
    point = _parse_vec3(w.get("point"), "off_plane.point")
    require_ok(point[2] != 0, "POINT_ON_PLANE", "off_plane point must have z != 0")
    radii = [parse_q(r, f"off_plane.radii[{i}]") for i, r in enumerate(radii_raw)]
    require_ok(all(r > 0 for r in radii), "NONPOSITIVE_RADIUS", "off_plane radii must be > 0")
    for s, r in zip(sites, radii):
        d2 = (point[0] - s[0]) ** 2 + (point[1] - s[1]) ** 2 + (point[2] - s[2]) ** 2
        require_ok(d2 == r * r, "RADIUS_SQ_MISMATCH", "off_plane ||x - a_i||^2 != r_i^2")
    vz = sum((_q_div(-q * point[2], _q_pow_pos(r, 3)) for q, r in zip(charges, radii)), Fraction(0, 1))
    z_den = point[2].denominator
    w_q = [_q_div(q, _q_pow_pos(r, 3)) for q, r in zip(charges, radii)]
    l_w = _lcm_many([wq.denominator for wq in w_q])
    return {
        "name": w["name"],
        "z_q": fmt_q(point[2]),
        "vz_q": fmt_q(vz),
        "z_num": int(point[2] * z_den),
        "z_den": z_den,
        "weight_scale": l_w,
        "int_weights": [int(wq * l_w) for wq in w_q],
    }


def _boxes_strictly_separated(a: list[int], b: list[int]) -> bool:
    for axis in (0, 2, 4):
        if a[axis + 1] < b[axis]:
            return True
        if b[axis + 1] < a[axis]:
            return True
    return False


def validate_box_collection_witness(w: Any) -> dict[str, Any]:
    require_ok(type(w) is dict, "INVALID_WITNESS", "box collection witness must be dict")
    name = w.get("name")
    require_ok(type(name) is str, "INVALID_WITNESS", "box collection name must be str")
    require_ok(bool(name.strip()), "INVALID_WITNESS", "missing box collection name")
    n_charges = require_exact_int(w.get("n_charges"), f"{name}.n_charges", min_val=2, max_val=64)
    grid_bits = require_exact_int(w.get("grid_shift_bits"), f"{name}.grid_shift_bits", min_val=8, max_val=128)
    sites_raw = w.get("sites")
    charges_raw = w.get("charges")
    boxes_raw = w.get("boxes")
    require_ok(type(sites_raw) is list, "INVALID_WITNESS", f"{name}: sites must be list")
    require_ok(len(sites_raw) == n_charges, "INVALID_WITNESS", f"{name}: sites count mismatch")
    require_ok(type(charges_raw) is list, "INVALID_WITNESS", f"{name}: charges must be list")
    require_ok(len(charges_raw) == n_charges, "INVALID_WITNESS", f"{name}: charges count mismatch")
    require_ok(type(boxes_raw) is list, "INVALID_WITNESS", f"{name}: boxes must be list")
    require_ok(len(boxes_raw) >= 1, "INVALID_WITNESS", f"{name}: boxes empty")

    sites = [_parse_vec3(s, f"{name}.sites[{i}]") for i, s in enumerate(sites_raw)]
    require_ok(len(set(sites)) == n_charges, "DUPLICATE_SITE", f"{name}: duplicate charge site")
    charges = [parse_q(c, f"{name}.charges[{i}]") for i, c in enumerate(charges_raw)]
    require_ok(all(q > 0 for q in charges), "NONPOSITIVE_CHARGE", f"{name}: charges must be > 0")

    scale = 1 << grid_bits
    # Outward-rounded 2^grid_bits site enclosures: exact dyadic sites (denominator | 2^grid_bits)
    # yield zero-width point intervals [k, k], while non-dyadic sites yield width-1 intervals [k, k+1].
    grid_sites = [
        [_floor_mul_q(s[0], scale), _ceil_mul_q(s[0], scale),
         _floor_mul_q(s[1], scale), _ceil_mul_q(s[1], scale),
         _floor_mul_q(s[2], scale), _ceil_mul_q(s[2], scale)]
        for s in sites
    ]
    grid_boxes: list[list[int]] = []
    morse_indices: list[int] = []
    for i, b in enumerate(boxes_raw):
        require_ok(type(b) is dict, "INVALID_BOX", f"{name}.boxes[{i}] must be dict")
        ctr = _parse_vec3(b.get("center"), f"{name}.boxes[{i}].center")
        hw = parse_q(b.get("half_width"), f"{name}.boxes[{i}].half_width")
        require_ok(hw > 0, "NONPOSITIVE_HALF_WIDTH", f"{name}.boxes[{i}]: half_width must be > 0")
        m_idx = require_exact_int(b.get("morse_index"), f"{name}.boxes[{i}].morse_index", min_val=1, max_val=2)
        morse_indices.append(m_idx)
        grid_boxes.append([
            _floor_mul_q(ctr[0] - hw, scale), _ceil_mul_q(ctr[0] + hw, scale),
            _floor_mul_q(ctr[1] - hw, scale), _ceil_mul_q(ctr[1] + hw, scale),
            _floor_mul_q(ctr[2] - hw, scale), _ceil_mul_q(ctr[2] + hw, scale),
        ])

    n_boxes = len(grid_boxes)
    for i in range(n_boxes):
        bi = grid_boxes[i]
        for j in range(i + 1, n_boxes):
            require_ok(_boxes_strictly_separated(bi, grid_boxes[j]), "OVERLAPPING_BOXES", f"{name}: boxes {i} and {j} overlap")
        for s_idx, gs in enumerate(grid_sites):
            require_ok(_boxes_strictly_separated(gs, bi), "SITE_IN_BOX", f"{name}: box {i} touches site {s_idx}")

    m1 = sum(1 for m in morse_indices if m == 1)
    m2 = sum(1 for m in morse_indices if m == 2)
    require_ok(m2 - m1 == n_charges - 1, "POINCARE_HOPF_MISMATCH", f"{name}: m2 - m1 != n - 1")
    maxwell_bound = (n_charges - 1) ** 2
    return {
        "name": name,
        "n_charges": n_charges,
        "grid_shift_bits": grid_bits,
        "n_boxes": n_boxes,
        "pairwise_disjoint_pairs": divmod(n_boxes * (n_boxes - 1), 2)[0],
        "site_avoidance_pairs": n_boxes * n_charges,
        "morse_counts": [0, m1, m2, 0],
        "m2_minus_m1": m2 - m1,
        "maxwell_bound": maxwell_bound,
        "exceeds_maxwell_bound": n_boxes > maxwell_bound,
        "morse_indices": morse_indices,
        "grid_sites": grid_sites,
        "grid_boxes": grid_boxes,
    }


def _sqrt_bracket_q(r2: Fraction, bits: int) -> tuple[Fraction, Fraction]:
    scale = 1 << bits
    fl = math.isqrt(divmod(r2.numerator * scale * scale, r2.denominator)[0])
    require_ok(fl > 0, "BRACKET_UNDERFLOW", "sqrt lower bracket reached 0")
    return Fraction(fl, scale), Fraction(fl + 1, scale)


def _eval_thomson_s2_config(sites: list[tuple[Fraction, Fraction, Fraction]], bits: int, label: str) -> dict[str, Any]:
    require_ok(len(sites) >= 2, "INVALID_THOMSON_SITES", f"{label}: need >= 2 sites")
    require_ok(len(set(sites)) == len(sites), "INVALID_THOMSON_SITES", f"{label}: duplicate sites")
    for i, p in enumerate(sites):
        norm2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
        require_ok(norm2 == 1, "NOT_ON_UNIT_SPHERE", f"{label}: site[{i}] norm^2 != 1")
    e2 = Fraction(0, 1)
    e1_lo = Fraction(0, 1)
    e1_hi = Fraction(0, 1)
    one_q = Fraction(1, 1)
    pair_brackets: list[dict[str, Any]] = []
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            d2 = (sites[i][0] - sites[j][0]) ** 2 + (sites[i][1] - sites[j][1]) ** 2 + (sites[i][2] - sites[j][2]) ** 2
            e2 += _q_div(one_q, d2)
            lo, hi = _sqrt_bracket_q(d2, bits)
            e1_lo += _q_div(one_q, hi)
            e1_hi += _q_div(one_q, lo)
            pair_brackets.append({
                "i": i, "j": j,
                "r2_n": d2.numerator, "r2_d": d2.denominator,
                "l_n": lo.numerator, "l_d": lo.denominator,
                "u_n": hi.numerator, "u_d": hi.denominator,
            })
    return {
        "n_charges": len(sites),
        "e2": e2,
        "e1_lo": e1_lo,
        "e1_hi": e1_hi,
        "e2_q": fmt_q(e2),
        "e1_lo_q": fmt_q(e1_lo),
        "e1_hi_q": fmt_q(e1_hi),
        "sites_cleared": [[p[0].numerator, p[0].denominator, p[1].numerator, p[1].denominator, p[2].numerator, p[2].denominator] for p in sites],
        "pair_brackets": pair_brackets,
    }


def validate_thomson_pair_witness(w: Any) -> dict[str, Any]:
    require_ok(type(w) is dict, "INVALID_WITNESS", "thomson witness must be dict")
    name = w.get("name")
    require_ok(type(name) is str, "INVALID_WITNESS", "thomson witness name must be str")
    require_ok(bool(name.strip()), "INVALID_WITNESS", "missing thomson witness name")
    n_charges = require_exact_int(w.get("n_charges"), f"{name}.n_charges", min_val=2, max_val=32)
    bits = require_exact_int(w.get("sqrt_bracket_bits"), f"{name}.sqrt_bracket_bits", min_val=4, max_val=32)
    lo_raw = w.get("lower_energy_sites")
    hi_raw = w.get("higher_energy_sites")
    require_ok(type(lo_raw) is list, "INVALID_WITNESS", f"{name}: lower_energy_sites must be list")
    require_ok(type(hi_raw) is list, "INVALID_WITNESS", f"{name}: higher_energy_sites must be list")
    require_ok(len(lo_raw) == n_charges, "INVALID_WITNESS", f"{name}: lower site count mismatch")
    require_ok(len(hi_raw) == n_charges, "INVALID_WITNESS", f"{name}: higher site count mismatch")

    lo_sites = [_parse_vec3(s, f"{name}.lo[{i}]") for i, s in enumerate(lo_raw)]
    hi_sites = [_parse_vec3(s, f"{name}.hi[{i}]") for i, s in enumerate(hi_raw)]
    lo_ev = _eval_thomson_s2_config(lo_sites, bits, f"{name}.lo")
    hi_ev = _eval_thomson_s2_config(hi_sites, bits, f"{name}.hi")
    require_ok(lo_ev["e2"] < hi_ev["e2"], "THOMSON_E2_NOT_STRICT", f"{name}: E2(lo) not < E2(hi)")
    require_ok(lo_ev["e1_hi"] < hi_ev["e1_lo"], "THOMSON_E1_NOT_STRICT", f"{name}: E1_hi(lo) not < E1_lo(hi)")
    return {
        "name": name,
        "n_charges": n_charges,
        "lower_label": w["lower_energy_label"],
        "higher_label": w["higher_energy_label"],
        "lower": lo_ev,
        "higher": hi_ev,
    }


def validate_sphere_bridges(w: Any) -> dict[str, Any]:
    require_ok(type(w) is dict, "INVALID_WITNESS", "sphere_bridge_witnesses must be dict")
    pairs_raw = w.get("stereographic_s2q_pairs")
    s3_raw = w.get("cayley_s3q_points")
    require_ok(type(pairs_raw) is list, "INVALID_WITNESS", "stereographic_s2q_pairs must be list")
    require_ok(len(pairs_raw) >= 2, "INVALID_WITNESS", "need >= 2 stereographic pairs")
    require_ok(type(s3_raw) is list, "INVALID_WITNESS", "cayley_s3q_points must be list")
    require_ok(len(s3_raw) >= 2, "INVALID_WITNESS", "need >= 2 cayley S3 points")

    one_q = Fraction(1, 1)
    stereo_out: list[dict[str, Any]] = []
    for idx, item in enumerate(pairs_raw):
        require_ok(type(item) is dict, "INVALID_WITNESS", f"stereo[{idx}] must be dict")
        u1 = parse_q(item.get("u1"), f"stereo[{idx}].u1")
        v1 = parse_q(item.get("v1"), f"stereo[{idx}].v1")
        u2 = parse_q(item.get("u2"), f"stereo[{idx}].u2")
        v2 = parse_q(item.get("v2"), f"stereo[{idx}].v2")
        d1 = one_q + u1 * u1 + v1 * v1
        d2 = one_q + u2 * u2 + v2 * v2
        p1 = (_q_div(2 * u1, d1), _q_div(2 * v1, d1), _q_div(u1 * u1 + v1 * v1 - one_q, d1))
        p2 = (_q_div(2 * u2, d2), _q_div(2 * v2, d2), _q_div(u2 * u2 + v2 * v2 - one_q, d2))
        chord_direct = (p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2
        stereo_out.append({
            "u1": fmt_q(u1), "v1": fmt_q(v1), "u2": fmt_q(u2), "v2": fmt_q(v2),
            "u1_nd": [u1.numerator, u1.denominator], "v1_nd": [v1.numerator, v1.denominator],
            "u2_nd": [u2.numerator, u2.denominator], "v2_nd": [v2.numerator, v2.denominator],
            "p1": [fmt_q(c) for c in p1], "p2": [fmt_q(c) for c in p2],
            "chord_sq_q": fmt_q(chord_direct),
            "p1_nd": [[c.numerator, c.denominator] for c in p1],
            "p2_nd": [[c.numerator, c.denominator] for c in p2],
            "chord_nd": [chord_direct.numerator, chord_direct.denominator],
        })

    s3_out: list[dict[str, Any]] = []
    for idx, item in enumerate(s3_raw):
        require_ok(type(item) is dict, "INVALID_WITNESS", f"cayley[{idx}] must be dict")
        w1 = parse_q(item.get("w1"), f"cayley[{idx}].w1")
        w2 = parse_q(item.get("w2"), f"cayley[{idx}].w2")
        w3 = parse_q(item.get("w3"), f"cayley[{idx}].w3")
        s = w1 * w1 + w2 * w2 + w3 * w3
        denom = one_q + s
        q4 = (_q_div(one_q - s, denom), _q_div(2 * w1, denom), _q_div(2 * w2, denom), _q_div(2 * w3, denom))
        s3_out.append({
            "w": [fmt_q(w1), fmt_q(w2), fmt_q(w3)],
            "w_nd": [[w1.numerator, w1.denominator], [w2.numerator, w2.denominator], [w3.numerator, w3.denominator]],
            "q4": [fmt_q(c) for c in q4],
            "q4_nd": [[c.numerator, c.denominator] for c in q4],
            "equatorial_matches_stereographic": w3 == 0,
        })
    return {"stereographic_pairs": stereo_out, "cayley_s3_points": s3_out}


def validate_witness_pack(pack: Any) -> dict[str, Any]:
    require_ok(type(pack) is dict, "INVALID_PACK", "witness pack must be dict")
    require_ok(pack.get("schema") == WITNESS_SCHEMA, "SCHEMA_MISMATCH", "witness pack schema mismatch")
    require_ok(pack.get("upstream_issue") == UPSTREAM_ISSUE, "ISSUE_MISMATCH", "upstream_issue mismatch")

    eq_raw = pack.get("equilibrium_witnesses")
    box_raw = pack.get("box_collection_witnesses")
    th_raw = pack.get("thomson_pair_witnesses")
    require_ok(type(eq_raw) is list, "INVALID_PACK", "equilibrium_witnesses must be list")
    require_ok(len(eq_raw) == 3, "INVALID_PACK", "expected 3 equilibrium witnesses")
    require_ok(type(box_raw) is list, "INVALID_PACK", "box_collection_witnesses must be list")
    require_ok(len(box_raw) == 2, "INVALID_PACK", "expected 2 box collection witnesses")
    require_ok(type(th_raw) is list, "INVALID_PACK", "thomson_pair_witnesses must be list")
    require_ok(len(th_raw) == 2, "INVALID_PACK", "expected 2 thomson pair witnesses")

    eq_v = [validate_equilibrium_witness(w) for w in eq_raw]
    off_v = validate_off_plane_witness(pack.get("off_plane_confinement_witness"))
    box_v = [validate_box_collection_witness(w) for w in box_raw]
    th_v = [validate_thomson_pair_witness(w) for w in th_raw]
    br_v = validate_sphere_bridges(pack.get("sphere_bridge_witnesses"))

    require_ok(
        [e["name"] for e in eq_v] == [
            "n2_analytic_control",
            "n3_weierstrass_incenter_iso_acute",
            "n3_weierstrass_incenter_scalene",
        ],
        "WITNESS_SET_MISMATCH",
        "unexpected equilibrium witness names",
    )
    require_ok(
        [b["name"] for b in box_v] == ["n4_eleven_short", "n8_abk_milnor8_flat"],
        "WITNESS_SET_MISMATCH",
        "unexpected box collection witness names",
    )
    require_ok(
        [t["name"] for t in th_v] == [
            "thomson_n4_s2q_near_tet_vs_square",
            "thomson_n6_s2q_octahedron_vs_equatorial6",
        ],
        "WITNESS_SET_MISMATCH",
        "unexpected thomson witness names",
    )
    return {
        "equilibria": eq_v,
        "off_plane": off_v,
        "box_collections": box_v,
        "thomson_pairs": th_v,
        "sphere_bridges": br_v,
    }


# ---------------------------------------------------------------------------
# Deterministic Pure-Init Lean 4.36.0 Code Generator
# ---------------------------------------------------------------------------

def _lean_int(v: int) -> str:
    return f"({v} : Int)" if v < 0 else str(v)


def _emit_eq_terms_lean(eq: dict[str, Any]) -> tuple[str, str]:
    g_items = [
        f"⟨{_lean_int(w)}, ⟨{_lean_int(d[0])}, {_lean_int(d[1])}, {_lean_int(d[2])}⟩⟩"
        for w, d in zip(eq["int_weights_g"], eq["int_disps"])
    ]
    h_items = [
        f"⟨{_lean_int(w)}, {_lean_int(d[0])}, {_lean_int(d[1])}, {_lean_int(d[2])}⟩"
        for w, d in zip(eq["int_weights_h"], eq["int_disps"])
    ]
    return "[" + ", ".join(g_items) + "]", "[" + ", ".join(h_items) + "]"


def _emit_box_list_lean(boxes: list[list[int]]) -> str:
    rows = [
        f"  ⟨{_lean_int(b[0])}, {_lean_int(b[1])}, {_lean_int(b[2])}, {_lean_int(b[3])}, {_lean_int(b[4])}, {_lean_int(b[5])}⟩"
        for b in boxes
    ]
    return "[\n" + ",\n".join(rows) + "\n]"


def _emit_thomson_brackets_lean(cfg: dict[str, Any]) -> str:
    rows = [
        f"  ⟨{b['r2_n']}, {b['r2_d']}, {b['l_n']}, {b['l_d']}, {b['u_n']}, {b['u_d']}⟩"
        for b in cfg["pair_brackets"]
    ]
    return "[\n" + ",\n".join(rows) + "\n]"


def _emit_thomson_sites_lean(cfg: dict[str, Any]) -> str:
    rows = [
        f"  ⟨{_lean_int(s[0])}, {s[1]}, {_lean_int(s[2])}, {s[3]}, {_lean_int(s[4])}, {s[5]}⟩"
        for s in cfg["sites_cleared"]
    ]
    return "[\n" + ",\n".join(rows) + "\n]"


def compile_axiom_audit_source() -> str:
    """Deterministically compile `AxiomAudit.lean` importing `PointChargesExactCert`."""
    lines = [
        "-- Pure-Init Axiom Audit for PointChargesExactCert (Lean 4.36.0)",
        "import PointChargesExactCert",
        "",
    ]
    for thm in EXPECTED_THEOREM_NAMES:
        lines.append(f"#print axioms {thm}")
    lines.append("")
    return "\n".join(lines)


def compile_lean4_source(validated: dict[str, Any]) -> str:
    """Deterministically compile validated exact-rational witnesses into pure-Init Lean 4.36.0 source."""
    eq_n2, eq_iso, eq_sca = validated["equilibria"]
    off_p = validated["off_plane"]
    box_n4, box_n8 = validated["box_collections"]
    th_n4, th_n6 = validated["thomson_pairs"]
    br = validated["sphere_bridges"]

    n2_g, n2_h = _emit_eq_terms_lean(eq_n2)
    iso_g, iso_h = _emit_eq_terms_lean(eq_iso)
    sca_g, sca_h = _emit_eq_terms_lean(eq_sca)

    off_weights = "[" + ", ".join(str(x) for x in off_p["int_weights"]) + "]"
    n4_morse = "[" + ", ".join(str(x) for x in box_n4["morse_indices"]) + "]"
    n8_morse = "[" + ", ".join(str(x) for x in box_n8["morse_indices"]) + "]"

    sp0 = br["stereographic_pairs"][0]
    sp1 = br["stereographic_pairs"][1]
    cs0 = br["cayley_s3_points"][0]
    cs1 = br["cayley_s3_points"][1]

    lines: list[str] = [
        "-- Pure-Init Lean 4.36.0 Point-Charge Exact Certificate Module",
        "-- Deterministically compiled from exact rational witnesses (Foundation-0/point-charges#6).",
        "-- Zero external imports, zero Mathlib, zero floats, zero sorry, zero Classical.choice,",
        "-- zero bv_decide, zero native_decide, zero warnings.",
        "",
        "namespace PointChargesExactCert",
        "",
        "-- ============================================================================",
        "-- 0. Core Exact Integer / Rational Structures (Pure Init)",
        "-- ============================================================================",
        "",
        "structure Vec3Z where",
        "  x : Int",
        "  y : Int",
        "  z : Int",
        "  deriving DecidableEq",
        "",
        "def vecZero : Vec3Z := ⟨0, 0, 0⟩",
        "def e1 : Vec3Z := ⟨1, 0, 0⟩",
        "def e2 : Vec3Z := ⟨0, 1, 0⟩",
        "def e3 : Vec3Z := ⟨0, 0, 1⟩",
        "",
        "def vecAdd (a b : Vec3Z) : Vec3Z := ⟨a.x + b.x, a.y + b.y, a.z + b.z⟩",
        "def vecScale (s : Int) (v : Vec3Z) : Vec3Z := ⟨s * v.x, s * v.y, s * v.z⟩",
        "def vecNormSq (v : Vec3Z) : Int := v.x * v.x + v.y * v.y + v.z * v.z",
        "",
        "structure Sym3Z where",
        "  h11 : Int",
        "  h22 : Int",
        "  h33 : Int",
        "  h12 : Int",
        "  h13 : Int",
        "  h23 : Int",
        "  deriving DecidableEq",
        "",
        "def symZero : Sym3Z := ⟨0, 0, 0, 0, 0, 0⟩",
        "def diagSym (d1 d2 d3 : Int) : Sym3Z := ⟨d1, d2, d3, 0, 0, 0⟩",
        "def symAdd (a b : Sym3Z) : Sym3Z :=",
        "  ⟨a.h11 + b.h11, a.h22 + b.h22, a.h33 + b.h33, a.h12 + b.h12, a.h13 + b.h13, a.h23 + b.h23⟩",
        "def symTrace (h : Sym3Z) : Int := h.h11 + h.h22 + h.h33",
        "def symDet (h : Sym3Z) : Int :=",
        "  h.h11 * (h.h22 * h.h33 - h.h23 * h.h23)",
        "  - h.h12 * (h.h12 * h.h33 - h.h23 * h.h13)",
        "  + h.h13 * (h.h12 * h.h23 - h.h22 * h.h13)",
        "",
        "def quadFormSym (h : Sym3Z) (v : Vec3Z) : Int :=",
        "  h.h11 * v.x * v.x + h.h22 * v.y * v.y + h.h33 * v.z * v.z",
        "  + 2 * h.h12 * v.x * v.y + 2 * h.h13 * v.x * v.z + 2 * h.h23 * v.y * v.z",
        "",
        "structure Mat3Z where",
        "  a11 : Int",
        "  a12 : Int",
        "  a13 : Int",
        "  a21 : Int",
        "  a22 : Int",
        "  a23 : Int",
        "  a31 : Int",
        "  a32 : Int",
        "  a33 : Int",
        "  deriving DecidableEq",
        "",
        "def matVec (s : Mat3Z) (u : Vec3Z) : Vec3Z :=",
        "  ⟨s.a11 * u.x + s.a12 * u.y + s.a13 * u.z,",
        "   s.a21 * u.x + s.a22 * u.y + s.a23 * u.z,",
        "   s.a31 * u.x + s.a32 * u.y + s.a33 * u.z⟩",
        "",
        "def bilSym (h : Sym3Z) (u v : Vec3Z) : Int :=",
        "  h.h11 * u.x * v.x + h.h22 * u.y * v.y + h.h33 * u.z * v.z",
        "  + h.h12 * (u.x * v.y + u.y * v.x)",
        "  + h.h13 * (u.x * v.z + u.z * v.x)",
        "  + h.h23 * (u.y * v.z + u.z * v.y)",
        "",
        "def congruenceSym (s : Mat3Z) (h : Sym3Z) : Sym3Z :=",
        "  let c1 : Vec3Z := ⟨s.a11, s.a21, s.a31⟩",
        "  let c2 : Vec3Z := ⟨s.a12, s.a22, s.a32⟩",
        "  let c3 : Vec3Z := ⟨s.a13, s.a23, s.a33⟩",
        "  ⟨quadFormSym h c1, quadFormSym h c2, quadFormSym h c3,",
        "   bilSym h c1 c2, bilSym h c1 c3, bilSym h c2 c3⟩",
        "",
        "structure CoulombTermZ where",
        "  w : Int",
        "  dx : Int",
        "  dy : Int",
        "  dz : Int",
        "  deriving DecidableEq",
        "",
        "structure GradTermZ where",
        "  w : Int",
        "  d : Vec3Z",
        "  deriving DecidableEq",
        "",
        "def coulombSingleHessian (w dx dy dz : Int) : Sym3Z :=",
        "  let r2 := dx * dx + dy * dy + dz * dz",
        "  ⟨w * (3 * dx * dx - r2), w * (3 * dy * dy - r2), w * (3 * dz * dz - r2),",
        "   w * (3 * dx * dy), w * (3 * dx * dz), w * (3 * dy * dz)⟩",
        "",
        "def totalCoulombHessian : List CoulombTermZ → Sym3Z",
        "  | [] => symZero",
        "  | t :: ts => symAdd (coulombSingleHessian t.w t.dx t.dy t.dz) (totalCoulombHessian ts)",
        "",
        "def totalCoulombGrad : List GradTermZ → Vec3Z",
        "  | [] => vecZero",
        "  | t :: ts => vecAdd (vecScale (-t.w) t.d) (totalCoulombGrad ts)",
        "",
        "def sumPosList : List Int → Int",
        "  | [] => 0",
        "  | w :: ws => w + sumPosList ws",
        "",
        "def allPosNonEmpty : List Int → Prop",
        "  | [] => False",
        "  | [w] => 0 < w",
        "  | w :: w2 :: ws => 0 < w ∧ allPosNonEmpty (w2 :: ws)",
        "",
        "def planarForceZ (z : Int) (ws : List Int) : Int := -z * sumPosList ws",
        "def planarHessZZ (ws : List Int) : Int := -sumPosList ws",
        "def planarInPlaneTrace (ws : List Int) : Int := sumPosList ws",
        "",
        "structure Box3Z where",
        "  lx : Int",
        "  ux : Int",
        "  ly : Int",
        "  uy : Int",
        "  lz : Int",
        "  uz : Int",
        "  deriving DecidableEq",
        "",
        "def inBox (x : Vec3Z) (b : Box3Z) : Prop :=",
        "  b.lx ≤ x.x ∧ x.x ≤ b.ux ∧ b.ly ≤ x.y ∧ x.y ≤ b.uy ∧ b.lz ≤ x.z ∧ x.z ≤ b.uz",
        "",
        "def boxDisjoint (a b : Box3Z) : Bool :=",
        "  decide (a.ux < b.lx) || decide (b.ux < a.lx) ||",
        "  decide (a.uy < b.ly) || decide (b.uy < a.ly) ||",
        "  decide (a.uz < b.lz) || decide (b.uz < a.lz)",
        "",
        "def allPairsDisjoint : List Box3Z → Bool",
        "  | [] => true",
        "  | b :: bs => bs.all (fun c => boxDisjoint b c) && allPairsDisjoint bs",
        "",
        "def allBoxesAvoidSites (boxes sites : List Box3Z) : Bool :=",
        "  boxes.all (fun b => sites.all (fun s => boxDisjoint b s))",
        "",
        "def countVal (k : Nat) (xs : List Nat) : Nat :=",
        "  xs.foldl (fun acc x => if x == k then acc + 1 else acc) 0",
        "",
        "structure RationalSite3 where",
        "  xn : Int",
        "  xd : Nat",
        "  yn : Int",
        "  yd : Nat",
        "  zn : Int",
        "  zd : Nat",
        "  deriving DecidableEq",
        "",
        "def siteOnUnitSphere (p : RationalSite3) : Bool :=",
        "  let xd2 : Int := (p.xd * p.xd : Nat)",
        "  let yd2 : Int := (p.yd * p.yd : Nat)",
        "  let zd2 : Int := (p.zd * p.zd : Nat)",
        "  decide (0 < p.xd ∧ 0 < p.yd ∧ 0 < p.zd) &&",
        "  decide (p.xn * p.xn * yd2 * zd2 + p.yn * p.yn * xd2 * zd2 + p.zn * p.zn * xd2 * yd2 = xd2 * yd2 * zd2)",
        "",
        "structure SqrtBracket where",
        "  r2n : Nat",
        "  r2d : Nat",
        "  ln : Nat",
        "  ld : Nat",
        "  un : Nat",
        "  ud : Nat",
        "  deriving DecidableEq",
        "",
        "def validSqrtBracket (b : SqrtBracket) : Bool :=",
        "  decide (0 < b.r2n ∧ 0 < b.r2d ∧ 0 < b.ln ∧ 0 < b.ld ∧ 0 < b.un ∧ 0 < b.ud) &&",
        "  decide (b.ln * b.ln * b.r2d ≤ b.r2n * b.ld * b.ld) &&",
        "  decide (b.r2n * b.ud * b.ud ≤ b.un * b.un * b.r2d)",
        "",
        "structure QNat where",
        "  n : Nat",
        "  d : Nat",
        "  deriving DecidableEq",
        "",
        "def qAdd (a b : QNat) : QNat := ⟨a.n * b.d + b.n * a.d, a.d * b.d⟩",
        "def qLt (a b : QNat) : Bool := decide (a.n * b.d < b.n * a.d)",
        "def sumE2 (bs : List SqrtBracket) : QNat :=",
        "  bs.foldl (fun acc b => qAdd acc ⟨b.r2d, b.r2n⟩) ⟨0, 1⟩",
        "def sumE1Lo (bs : List SqrtBracket) : QNat :=",
        "  bs.foldl (fun acc b => qAdd acc ⟨b.ud, b.un⟩) ⟨0, 1⟩",
        "def sumE1Hi (bs : List SqrtBracket) : QNat :=",
        "  bs.foldl (fun acc b => qAdd acc ⟨b.ld, b.ln⟩) ⟨0, 1⟩",
        "",
        "def checkPairDistSq (p q : RationalSite3) (r2n r2d : Nat) : Bool :=",
        "  let dxn := p.xn * (q.xd : Int) - q.xn * (p.xd : Int)",
        "  let dxd : Int := (p.xd * q.xd : Nat)",
        "  let dyn := p.yn * (q.yd : Int) - q.yn * (p.yd : Int)",
        "  let dyd : Int := (p.yd * q.yd : Nat)",
        "  let dzn := p.zn * (q.zd : Int) - q.zn * (p.zd : Int)",
        "  let dzd : Int := (p.zd * q.zd : Nat)",
        "  let num := dxn * dxn * dyd * dyd * dzd * dzd + dyn * dyn * dxd * dxd * dzd * dzd + dzn * dzn * dxd * dxd * dyd * dyd",
        "  let den := dxd * dxd * dyd * dyd * dzd * dzd",
        "  decide (num * (r2d : Int) = (r2n : Int) * den)",
        "",
        "def pairDistsMatch : List RationalSite3 → List SqrtBracket → Bool",
        "  | [], [] => true",
        "  | p :: ps, bs =>",
        "    let k := ps.length",
        "    let head := bs.take k",
        "    let tail := bs.drop k",
        "    (head.length == k) &&",
        "    ((ps.zip head).all (fun (q, b) => checkPairDistSq p q b.r2n b.r2d)) &&",
        "    pairDistsMatch ps tail",
        "  | _, _ => false",
        "",
        "/-- Rational point on the 3-sphere S^3(Q) stored as four (numerator, positive denominator) pairs. -/",
        "structure RationalSite4 where",
        "  wn : Int",
        "  wd : Nat",
        "  xn : Int",
        "  xd : Nat",
        "  yn : Int",
        "  yd : Nat",
        "  zn : Int",
        "  zd : Nat",
        "  deriving DecidableEq",
        "",
        "/-- Exact rational equality with positive-denominator guards: aN / aD = bN / bD. -/",
        "def ratEq (aN : Int) (aD : Nat) (bN : Int) (bD : Nat) : Bool :=",
        "  decide (0 < aD ∧ 0 < bD ∧ aN * (bD : Int) = bN * (aD : Int))",
        "",
        "/-- Cleared-denominator S^3(Q) unit-norm check: w^2 + x^2 + y^2 + z^2 = 1. -/",
        "def siteOnUnitSphere4 (q : RationalSite4) : Bool :=",
        "  let wd2 : Int := (q.wd * q.wd : Nat)",
        "  let xd2 : Int := (q.xd * q.xd : Nat)",
        "  let yd2 : Int := (q.yd * q.yd : Nat)",
        "  let zd2 : Int := (q.zd * q.zd : Nat)",
        "  decide (0 < q.wd ∧ 0 < q.xd ∧ 0 < q.yd ∧ 0 < q.zd) &&",
        "  decide (",
        "    q.wn * q.wn * xd2 * yd2 * zd2 +",
        "    q.xn * q.xn * wd2 * yd2 * zd2 +",
        "    q.yn * q.yn * wd2 * xd2 * zd2 +",
        "    q.zn * q.zn * wd2 * xd2 * yd2 = wd2 * xd2 * yd2 * zd2",
        "  )",
        "",
        "/-- Verify that `p` equals the exact inverse stereographic image of `(un/ud, vn/vd)` on S^2(Q). -/",
        "def stereographicMatches (un : Int) (ud : Nat) (vn : Int) (vd : Nat) (p : RationalSite3) : Bool :=",
        "  let ud2 : Int := (ud * ud : Nat)",
        "  let vd2 : Int := (vd * vd : Nat)",
        "  let sNum : Int := un * un * vd2 + vn * vn * ud2",
        "  let dDen : Int := ud2 * vd2",
        "  let qDen : Int := sNum + dDen",
        "  decide (0 < ud ∧ 0 < vd ∧ 0 < p.xd ∧ 0 < p.yd ∧ 0 < p.zd ∧ 0 < qDen) &&",
        "  decide (p.xn * qDen = 2 * un * (ud : Int) * vd2 * (p.xd : Int)) &&",
        "  decide (p.yn * qDen = 2 * vn * (vd : Int) * ud2 * (p.yd : Int)) &&",
        "  decide (p.zn * qDen = (sNum - dDen) * (p.zd : Int))",
        "",
        "/-- Verify that `q` equals the exact Cayley S^3(Q) map of `(w1n/w1d, w2n/w2d, w3n/w3d)`. -/",
        "def cayleyS3Matches (w1n : Int) (w1d : Nat) (w2n : Int) (w2d : Nat) (w3n : Int) (w3d : Nat) (q : RationalSite4) : Bool :=",
        "  let d1 : Int := (w1d * w1d : Nat)",
        "  let d2 : Int := (w2d * w2d : Nat)",
        "  let d3 : Int := (w3d * w3d : Nat)",
        "  let sNum : Int := w1n * w1n * d2 * d3 + w2n * w2n * d1 * d3 + w3n * w3n * d1 * d2",
        "  let dDen : Int := d1 * d2 * d3",
        "  let qDen : Int := dDen + sNum",
        "  decide (0 < w1d ∧ 0 < w2d ∧ 0 < w3d ∧ 0 < q.wd ∧ 0 < q.xd ∧ 0 < q.yd ∧ 0 < q.zd ∧ 0 < qDen) &&",
        "  decide (q.wn * qDen = (dDen - sNum) * (q.wd : Int)) &&",
        "  decide (q.xn * qDen = 2 * w1n * (w1d : Int) * d2 * d3 * (q.xd : Int)) &&",
        "  decide (q.yn * qDen = 2 * w2n * (w2d : Int) * d1 * d3 * (q.yd : Int)) &&",
        "  decide (q.zn * qDen = 2 * w3n * (w3d : Int) * d1 * d2 * (q.zd : Int))",
        "",
        "/-- Verify the 4-coordinate equatorial S^3(Q) <-> S^2(Q) bridge: (q_w, q_x, q_y, q_z) = (-p_z, p_x, p_y, 0). -/",
        "def equatorialCayleyMatchesStereo (q : RationalSite4) (p : RationalSite3) : Bool :=",
        "  ratEq q.wn q.wd (-p.zn) p.zd &&",
        "  ratEq q.xn q.xd p.xn p.xd &&",
        "  ratEq q.yn q.yd p.yn p.yd &&",
        "  ratEq q.zn q.zd 0 1",
        "",
        "-- ============================================================================",
        "-- 1. Fourteen Universal (forall) Pure-Init Electrostatic, Box, & Sphere Theorems",
        "-- ============================================================================",
        "",
        "private theorem two_mul_int (x : Int) : 2 * x = x + x := by omega",
        "private theorem three_mul_int (x : Int) : 3 * x = x + x + x := by omega",
        "private theorem four_mul_int (x : Int) : 4 * x = x + x + x + x := by omega",
        "private theorem mul_four_int (x : Int) : x * 4 = x + x + x + x := by omega",
        "",
        "private theorem symTrace_add (a b : Sym3Z) :",
        "    symTrace (symAdd a b) = symTrace a + symTrace b := by",
        "  simp only [symTrace, symAdd]",
        "  omega",
        "",
        "theorem t01_coulomb_single_hessian_trace_zero (w dx dy dz : Int) :",
        "    symTrace (coulombSingleHessian w dx dy dz) = 0 := by",
        "  simp only [symTrace, coulombSingleHessian]",
        "  rw [← Int.mul_add, ← Int.mul_add]",
        "  have h0 : 3 * dx * dx - (dx * dx + dy * dy + dz * dz)",
        "      + (3 * dy * dy - (dx * dx + dy * dy + dz * dz))",
        "      + (3 * dz * dz - (dx * dx + dy * dy + dz * dz)) = 0 := by",
        "    rw [Int.mul_assoc 3 dx dx, Int.mul_assoc 3 dy dy, Int.mul_assoc 3 dz dz]",
        "    omega",
        "  rw [h0, Int.mul_zero]",
        "",
        "theorem t02_earnshaw_total_hessian_trace_zero (terms : List CoulombTermZ) :",
        "    symTrace (totalCoulombHessian terms) = 0 := by",
        "  induction terms with",
        "  | nil => rfl",
        "  | cons t ts ih =>",
        "    simp only [totalCoulombHessian]",
        "    rw [symTrace_add, t01_coulomb_single_hessian_trace_zero, ih]",
        "    omega",
        "",
        "private theorem quadFormSym_e1 (h : Sym3Z) : quadFormSym h e1 = h.h11 := by",
        "  simp only [quadFormSym, e1, Int.mul_one, Int.mul_zero, Int.add_zero]",
        "",
        "private theorem quadFormSym_e2 (h : Sym3Z) : quadFormSym h e2 = h.h22 := by",
        "  simp only [quadFormSym, e2, Int.mul_one, Int.mul_zero, Int.zero_add, Int.add_zero]",
        "",
        "private theorem quadFormSym_e3 (h : Sym3Z) : quadFormSym h e3 = h.h33 := by",
        "  simp only [quadFormSym, e3, Int.mul_one, Int.mul_zero, Int.zero_add, Int.add_zero]",
        "",
        "theorem t03_earnshaw_no_strict_definite_quadform (h : Sym3Z) (htr : symTrace h = 0) :",
        "    (¬ ∀ v : Vec3Z, v ≠ vecZero → 0 < quadFormSym h v) ∧",
        "    (¬ ∀ v : Vec3Z, v ≠ vecZero → quadFormSym h v < 0) := by",
        "  have he1 : e1 ≠ vecZero := by decide",
        "  have he2 : e2 ≠ vecZero := by decide",
        "  have he3 : e3 ≠ vecZero := by decide",
        "  constructor",
        "  · intro hpd",
        "    have h1 := hpd e1 he1; rw [quadFormSym_e1] at h1",
        "    have h2 := hpd e2 he2; rw [quadFormSym_e2] at h2",
        "    have h3 := hpd e3 he3; rw [quadFormSym_e3] at h3",
        "    simp only [symTrace] at htr",
        "    omega",
        "  · intro hnd",
        "    have h1 := hnd e1 he1; rw [quadFormSym_e1] at h1",
        "    have h2 := hnd e2 he2; rw [quadFormSym_e2] at h2",
        "    have h3 := hnd e3 he3; rw [quadFormSym_e3] at h3",
        "    simp only [symTrace] at htr",
        "    omega",
        "",
        "private theorem sumPosList_pos (ws : List Int) (h : allPosNonEmpty ws) : 0 < sumPosList ws := by",
        "  match ws with",
        "  | [] => exact False.elim h",
        "  | [w] => simp only [allPosNonEmpty, sumPosList] at h ⊢; omega",
        "  | w :: w2 :: rest =>",
        "    simp only [allPosNonEmpty] at h",
        "    have ih := sumPosList_pos (w2 :: rest) h.2",
        "    change 0 < w + sumPosList (w2 :: rest)",
        "    omega",
        "",
        "theorem t04_planar_offplane_force_and_z_curvature (z : Int) (ws : List Int) (hpos : allPosNonEmpty ws) :",
        "    (0 < z → planarForceZ z ws < 0) ∧",
        "    (z < 0 → 0 < planarForceZ z ws) ∧",
        "    planarHessZZ ws < 0 ∧",
        "    0 < planarInPlaneTrace ws := by",
        "  have hs : 0 < sumPosList ws := sumPosList_pos ws hpos",
        "  refine ⟨?_, ?_, ?_, ?_⟩",
        "  · intro hz",
        "    simp only [planarForceZ]",
        "    have hp : 0 < z * sumPosList ws := Int.mul_pos hz hs",
        "    rw [Int.neg_mul]; omega",
        "  · intro hz",
        "    simp only [planarForceZ]",
        "    have hnz : 0 < -z := by omega",
        "    exact Int.mul_pos hnz hs",
        "  · simp only [planarHessZZ]; omega",
        "  · simp only [planarInPlaneTrace]; exact hs",
        "",
        "theorem t05_lagrange_congruence_quadform_identity (s : Mat3Z) (h : Sym3Z) (u : Vec3Z) :",
        "    quadFormSym h (matVec s u) = quadFormSym (congruenceSym s h) u := by",
        "  simp only [quadFormSym, matVec, congruenceSym, bilSym]",
        "  simp only [two_mul_int, Int.add_mul, Int.mul_add]",
        "  simp only [Int.mul_comm, Int.mul_left_comm]",
        "  omega",
        "",
        "private theorem quadFormSym_diag (d1 d2 d3 : Int) (u : Vec3Z) :",
        "    quadFormSym (diagSym d1 d2 d3) u = d1 * u.x * u.x + d2 * u.y * u.y + d3 * u.z * u.z := by",
        "  simp only [quadFormSym, diagSym, Int.mul_zero, Int.zero_mul, Int.add_zero]",
        "",
        "theorem t06_lagrange_diagonal_morse_index2_witness (s : Mat3Z) (h : Sym3Z) (d1 d2 d3 : Int)",
        "    (hcong : congruenceSym s h = diagSym d1 d2 d3) (hd1 : 0 < d1) (hd2 : d2 < 0) (hd3 : d3 < 0) :",
        "    0 < quadFormSym h (matVec s e1) ∧ quadFormSym h (matVec s e2) < 0 ∧ quadFormSym h (matVec s e3) < 0 := by",
        "  rw [t05_lagrange_congruence_quadform_identity, t05_lagrange_congruence_quadform_identity,",
        "      t05_lagrange_congruence_quadform_identity, hcong]",
        "  rw [quadFormSym_diag, quadFormSym_diag, quadFormSym_diag]",
        "  simp only [e1, e2, e3, Int.mul_one, Int.mul_zero, Int.add_zero, Int.zero_add]",
        "  exact ⟨hd1, hd2, hd3⟩",
        "",
        "theorem t07_lagrange_diagonal_morse_index1_witness (s : Mat3Z) (h : Sym3Z) (d1 d2 d3 : Int)",
        "    (hcong : congruenceSym s h = diagSym d1 d2 d3) (hd1 : 0 < d1) (hd2 : 0 < d2) (hd3 : d3 < 0) :",
        "    0 < quadFormSym h (matVec s e1) ∧ 0 < quadFormSym h (matVec s e2) ∧ quadFormSym h (matVec s e3) < 0 := by",
        "  rw [t05_lagrange_congruence_quadform_identity, t05_lagrange_congruence_quadform_identity,",
        "      t05_lagrange_congruence_quadform_identity, hcong]",
        "  rw [quadFormSym_diag, quadFormSym_diag, quadFormSym_diag]",
        "  simp only [e1, e2, e3, Int.mul_one, Int.mul_zero, Int.add_zero, Int.zero_add]",
        "  exact ⟨hd1, hd2, hd3⟩",
        "",
        "theorem t08_disjoint_boxes_no_common_point (b1 b2 : Box3Z) (hd : boxDisjoint b1 b2 = true)",
        "    (x : Vec3Z) : ¬(inBox x b1 ∧ inBox x b2) := by",
        "  intro ⟨h1, h2⟩",
        "  simp only [inBox] at h1 h2",
        "  simp only [boxDisjoint, Bool.or_eq_true, decide_eq_true_eq] at hd",
        "  omega",
        "",
        "theorem t09_inv_sqrt_bracket_order_preservation (ln ld un ud r2n r2d : Nat)",
        "    (hr2d : 0 < r2d)",
        "    (hlo : ln * ln * r2d ≤ r2n * ld * ld) (hhi : r2n * ud * ud ≤ un * un * r2d) :",
        "    ln * ud ≤ un * ld := by",
        "  have hmul1 : (ln * ln * r2d) * (ud * ud) ≤ (r2n * ld * ld) * (ud * ud) := Nat.mul_le_mul_right _ hlo",
        "  have hmul2 : (r2n * ud * ud) * (ld * ld) ≤ (un * un * r2d) * (ld * ld) := Nat.mul_le_mul_right _ hhi",
        "  have hmid : (r2n * ld * ld) * (ud * ud) = (r2n * ud * ud) * (ld * ld) := by",
        "    simp only [Nat.mul_left_comm, Nat.mul_comm]",
        "  have hchain : (ln * ln * r2d) * (ud * ud) ≤ (un * un * r2d) * (ld * ld) := by omega",
        "  have hlhs : (ln * ln * r2d) * (ud * ud) = r2d * ((ln * ud) * (ln * ud)) := by",
        "    simp only [Nat.mul_left_comm, Nat.mul_comm]",
        "  have hrhs : (un * un * r2d) * (ld * ld) = r2d * ((un * ld) * (un * ld)) := by",
        "    simp only [Nat.mul_left_comm, Nat.mul_comm]",
        "  rw [hlhs, hrhs] at hchain",
        "  have hsq : (ln * ud) * (ln * ud) ≤ (un * ld) * (un * ld) := Nat.le_of_mul_le_mul_left hchain hr2d",
        "  if hle : ln * ud ≤ un * ld then",
        "    exact hle",
        "  else",
        "    have hgt : un * ld < ln * ud := by omega",
        "    have h1 : (un * ld) * (un * ld) < (ln * ud) * (ln * ud) := Nat.mul_lt_mul_of_lt_of_lt hgt hgt",
        "    omega",
        "",
        "theorem t10_stereographic_s2q_unit_sphere (un ud vn vd : Int) :",
        "    let a := 2 * un * ud * vd * vd",
        "    let b := 2 * vn * vd * ud * ud",
        "    let c := (un * un * vd * vd + vn * vn * ud * ud) - ud * ud * vd * vd",
        "    let d := (un * un * vd * vd + vn * vn * ud * ud) + ud * ud * vd * vd",
        "    a * a + b * b + c * c = d * d := by",
        "  dsimp only",
        "  simp only [two_mul_int, Int.sub_mul, Int.mul_sub, Int.add_mul, Int.mul_add]",
        "  simp only [Int.mul_comm, Int.mul_left_comm]",
        "  omega",
        "",
        "theorem t11_stereographic_s2q_chord_distance_sq (u1 v1 u2 v2 : Int) :",
        "    let d1 := 1 + u1 * u1 + v1 * v1",
        "    let d2 := 1 + u2 * u2 + v2 * v2",
        "    let dx := 2 * u1 * d2 - 2 * u2 * d1",
        "    let dy := 2 * v1 * d2 - 2 * v2 * d1",
        "    let dz := (u1 * u1 + v1 * v1 - 1) * d2 - (u2 * u2 + v2 * v2 - 1) * d1",
        "    dx * dx + dy * dy + dz * dz = 4 * ((u1 - u2) * (u1 - u2) + (v1 - v2) * (v1 - v2)) * d1 * d2 := by",
        "  dsimp only",
        "  simp only [two_mul_int, four_mul_int, Int.sub_mul, Int.mul_sub, Int.add_mul, Int.mul_add, Int.mul_one, Int.one_mul]",
        "  simp only [Int.mul_comm, Int.mul_left_comm]",
        "  omega",
        "",
        "theorem t12_cayley_s3q_unit_norm_and_equatorial_bridge (w1 w2 w3 : Int) :",
        "    let s := w1 * w1 + w2 * w2 + w3 * w3",
        "    ((1 - s) * (1 - s) + (2 * w1) * (2 * w1) + (2 * w2) * (2 * w2) + (2 * w3) * (2 * w3)",
        "      = (1 + s) * (1 + s)) ∧",
        "    (let s2 := w1 * w1 + w2 * w2; (1 - s2) = -(s2 - 1)) := by",
        "  dsimp only",
        "  constructor",
        "  · simp only [two_mul_int, Int.sub_mul, Int.mul_sub, Int.add_mul, Int.mul_add, Int.mul_one, Int.one_mul]",
        "    simp only [Int.mul_comm, Int.mul_left_comm]",
        "    omega",
        "  · omega",
        "",
        "theorem t13_n3_isosceles_cubic_cancellation_and_pythagorean_roots (h y k : Int) :",
        "    ((h - y) * (h - 3 * y) * (1 + y * y) - 3 * y * (h - y) * (h - y) * y",
        "      = (h - y) * (h - 3 * y - 2 * h * y * y)) ∧",
        "    (let hNum := 3 * k; let hDen := k * k - 2;",
        "      hNum * k * k - 3 * hDen * k - 2 * hNum = 0 ∧",
        "      2 * hNum * k * k - 3 * hDen * (2 * k) - hNum * 4 = 0) := by",
        "  dsimp only",
        "  refine ⟨?_, ?_, ?_⟩",
        "  · simp only [two_mul_int, three_mul_int, Int.sub_mul, Int.mul_sub, Int.add_mul, Int.mul_add, Int.mul_one]",
        "    simp only [Int.mul_comm, Int.mul_left_comm]",
        "    omega",
        "  · simp only [two_mul_int, three_mul_int, Int.sub_mul, Int.mul_sub, Int.add_mul, Int.mul_add]",
        "    simp only [Int.mul_comm]",
        "    omega",
        "  · simp only [two_mul_int, three_mul_int, mul_four_int, Int.sub_mul, Int.mul_sub, Int.add_mul, Int.mul_add]",
        "    simp only [Int.mul_comm]",
        "    omega",
        "",
        "theorem t14_regular_tetrahedron_s2q_sqrt3_obstruction : ∀ (a b : Nat), 0 < b → 3 * (b * b) ≠ a * a := by",
        "  intro a",
        "  induction a using Nat.strongRecOn with",
        "  | ind a ih =>",
        "    intro b hb hEq",
        "    have ha3 : a % 3 = 0 := by",
        "      have hmod : (a * a) % 3 = 0 := by omega",
        "      have hm := Nat.mul_mod a a 3",
        "      have hlt : a % 3 < 3 := Nat.mod_lt a (by decide)",
        "      rcases (show a % 3 = 0 ∨ a % 3 = 1 ∨ a % 3 = 2 by omega) with h0 | h1 | h2",
        "      · exact h0",
        "      · rw [h1] at hm; omega",
        "      · rw [h2] at hm; omega",
        "    have ha_eq : a = 3 * (a / 3) := by omega",
        "    have hb_eq : b * b = 3 * ((a / 3) * (a / 3)) := by",
        "      have h1 : 3 * (b * b) = 9 * ((a / 3) * (a / 3)) := by",
        "        rw [hEq]",
        "        conv => lhs; rw [ha_eq]",
        "        rw [Nat.mul_assoc, Nat.mul_left_comm (a / 3) 3 (a / 3), ← Nat.mul_assoc]",
        "      omega",
        "    have hb3 : b % 3 = 0 := by",
        "      have hmod : (b * b) % 3 = 0 := by omega",
        "      have hm := Nat.mul_mod b b 3",
        "      have hlt : b % 3 < 3 := Nat.mod_lt b (by decide)",
        "      rcases (show b % 3 = 0 ∨ b % 3 = 1 ∨ b % 3 = 2 by omega) with h0 | h1 | h2",
        "      · exact h0",
        "      · rw [h1] at hm; omega",
        "      · rw [h2] at hm; omega",
        "    have hba : a / 3 < a := by",
        "      have ha_pos : 0 < a := by",
        "        rcases Nat.eq_zero_or_pos a with rfl | hpos",
        "        · have : 0 < b * b := Nat.mul_pos hb hb",
        "          omega",
        "        · exact hpos",
        "      omega",
        "    have hb_div_pos : 0 < b / 3 := by omega",
        "    have hb_sub : 3 * ((b / 3) * (b / 3)) = (a / 3) * (a / 3) := by",
        "      have hb_eq2 : b = 3 * (b / 3) := by omega",
        "      have h2 : b * b = 9 * ((b / 3) * (b / 3)) := by",
        "        conv => lhs; rw [hb_eq2]",
        "        rw [Nat.mul_assoc, Nat.mul_left_comm (b / 3) 3 (b / 3), ← Nat.mul_assoc]",
        "      omega",
        "    exact ih (a / 3) hba (b / 3) hb_div_pos hb_sub",
        "",
        "-- ============================================================================",
        "-- 2. Eight Compiled Exact-Rational Witness Theorems (Emitted From Witness Pack)",
        "-- ============================================================================",
        "",
        f"def w_n2_grad : List GradTermZ := {n2_g}",
        f"def w_n2_hess : List CoulombTermZ := {n2_h}",
        "",
        "theorem t15_compiled_n2_analytic_equilibrium_and_hessian :",
        "    totalCoulombGrad w_n2_grad = vecZero ∧",
        f"    totalCoulombHessian w_n2_hess = diagSym {_lean_int(eq_n2['cleared_hess_diag'][0])} {_lean_int(eq_n2['cleared_hess_diag'][1])} {_lean_int(eq_n2['cleared_hess_diag'][2])} ∧",
        "    symTrace (totalCoulombHessian w_n2_hess) = 0 ∧",
        "    symDet (totalCoulombHessian w_n2_hess) ≠ 0 ∧",
        "    0 < quadFormSym (totalCoulombHessian w_n2_hess) e1 ∧",
        "    quadFormSym (totalCoulombHessian w_n2_hess) e2 < 0 ∧",
        "    quadFormSym (totalCoulombHessian w_n2_hess) e3 < 0 := by",
        "  decide",
        "",
        f"def w_n3_iso_grad : List GradTermZ := {iso_g}",
        f"def w_n3_iso_hess : List CoulombTermZ := {iso_h}",
        "",
        "theorem t16_compiled_n3_weierstrass_incenter_iso_acute :",
        "    totalCoulombGrad w_n3_iso_grad = vecZero ∧",
        f"    totalCoulombHessian w_n3_iso_hess = diagSym {_lean_int(eq_iso['cleared_hess_diag'][0])} {_lean_int(eq_iso['cleared_hess_diag'][1])} {_lean_int(eq_iso['cleared_hess_diag'][2])} ∧",
        "    symTrace (totalCoulombHessian w_n3_iso_hess) = 0 ∧",
        "    symDet (totalCoulombHessian w_n3_iso_hess) ≠ 0 ∧",
        "    0 < quadFormSym (totalCoulombHessian w_n3_iso_hess) e1 ∧",
        "    0 < quadFormSym (totalCoulombHessian w_n3_iso_hess) e2 ∧",
        "    quadFormSym (totalCoulombHessian w_n3_iso_hess) e3 < 0 := by",
        "  decide",
        "",
        f"def w_n3_sca_grad : List GradTermZ := {sca_g}",
        f"def w_n3_sca_hess : List CoulombTermZ := {sca_h}",
        "",
        "theorem t17_compiled_n3_weierstrass_incenter_scalene :",
        "    totalCoulombGrad w_n3_sca_grad = vecZero ∧",
        f"    totalCoulombHessian w_n3_sca_hess = diagSym {_lean_int(eq_sca['cleared_hess_diag'][0])} {_lean_int(eq_sca['cleared_hess_diag'][1])} {_lean_int(eq_sca['cleared_hess_diag'][2])} ∧",
        "    symTrace (totalCoulombHessian w_n3_sca_hess) = 0 ∧",
        "    symDet (totalCoulombHessian w_n3_sca_hess) ≠ 0 ∧",
        "    0 < quadFormSym (totalCoulombHessian w_n3_sca_hess) e1 ∧",
        "    0 < quadFormSym (totalCoulombHessian w_n3_sca_hess) e2 ∧",
        "    quadFormSym (totalCoulombHessian w_n3_sca_hess) e3 < 0 := by",
        "  decide",
        "",
        f"def w_offplane_weights : List Int := {off_weights}",
        "",
        "theorem t18_compiled_n3_pythagorean_off_plane_confinement :",
        f"    planarForceZ {_lean_int(off_p['z_num'])} w_offplane_weights < 0 ∧",
        "    planarHessZZ w_offplane_weights < 0 ∧",
        "    0 < planarInPlaneTrace w_offplane_weights := by",
        "  decide",
        "",
        f"def w_n4_sites : List Box3Z := {_emit_box_list_lean(box_n4['grid_sites'])}",
        f"def w_n4_boxes : List Box3Z := {_emit_box_list_lean(box_n4['grid_boxes'])}",
        f"def w_n4_morse : List Nat := {n4_morse}",
        "",
        "theorem t19_compiled_n4_eleven_disjointness_and_morse :",
        "    w_n4_boxes.length = 11 ∧",
        "    allPairsDisjoint w_n4_boxes = true ∧",
        "    allBoxesAvoidSites w_n4_boxes w_n4_sites = true ∧",
        "    countVal 1 w_n4_morse = 4 ∧ countVal 2 w_n4_morse = 7 ∧",
        "    countVal 2 w_n4_morse - countVal 1 w_n4_morse = 4 - 1 ∧",
        "    (4 - 1) * (4 - 1) < w_n4_boxes.length := by",
        "  decide",
        "",
        f"def w_n8_sites : List Box3Z := {_emit_box_list_lean(box_n8['grid_sites'])}",
        f"def w_n8_boxes : List Box3Z := {_emit_box_list_lean(box_n8['grid_boxes'])}",
        f"def w_n8_morse : List Nat := {n8_morse}",
        "",
        "theorem t20_compiled_n8_fiftyone_disjointness_and_morse :",
        "    w_n8_boxes.length = 51 ∧",
        "    allPairsDisjoint w_n8_boxes = true ∧",
        "    allBoxesAvoidSites w_n8_boxes w_n8_sites = true ∧",
        "    countVal 1 w_n8_morse = 22 ∧ countVal 2 w_n8_morse = 29 ∧",
        "    countVal 2 w_n8_morse - countVal 1 w_n8_morse = 8 - 1 ∧",
        "    (8 - 1) * (8 - 1) < w_n8_boxes.length := by",
        "  decide",
        "",
        f"def w_th_n4_lo_sites : List RationalSite3 := {_emit_thomson_sites_lean(th_n4['lower'])}",
        f"def w_th_n4_hi_sites : List RationalSite3 := {_emit_thomson_sites_lean(th_n4['higher'])}",
        f"def w_th_n4_lo_br : List SqrtBracket := {_emit_thomson_brackets_lean(th_n4['lower'])}",
        f"def w_th_n4_hi_br : List SqrtBracket := {_emit_thomson_brackets_lean(th_n4['higher'])}",
        f"def w_th_n6_lo_sites : List RationalSite3 := {_emit_thomson_sites_lean(th_n6['lower'])}",
        f"def w_th_n6_hi_sites : List RationalSite3 := {_emit_thomson_sites_lean(th_n6['higher'])}",
        f"def w_th_n6_lo_br : List SqrtBracket := {_emit_thomson_brackets_lean(th_n6['lower'])}",
        f"def w_th_n6_hi_br : List SqrtBracket := {_emit_thomson_brackets_lean(th_n6['higher'])}",
        "",
        "theorem t21_compiled_thomson_s2q_energy_dominance :",
        "    w_th_n4_lo_sites.all siteOnUnitSphere = true ∧",
        "    w_th_n4_hi_sites.all siteOnUnitSphere = true ∧",
        "    pairDistsMatch w_th_n4_lo_sites w_th_n4_lo_br = true ∧",
        "    pairDistsMatch w_th_n4_hi_sites w_th_n4_hi_br = true ∧",
        "    w_th_n4_lo_br.all validSqrtBracket = true ∧",
        "    w_th_n4_hi_br.all validSqrtBracket = true ∧",
        "    qLt (sumE2 w_th_n4_lo_br) (sumE2 w_th_n4_hi_br) = true ∧",
        "    qLt (sumE1Hi w_th_n4_lo_br) (sumE1Lo w_th_n4_hi_br) = true ∧",
        "    qLt ⟨9, 4⟩ (sumE2 w_th_n4_lo_br) = true ∧",
        "    w_th_n6_lo_sites.all siteOnUnitSphere = true ∧",
        "    w_th_n6_hi_sites.all siteOnUnitSphere = true ∧",
        "    pairDistsMatch w_th_n6_lo_sites w_th_n6_lo_br = true ∧",
        "    pairDistsMatch w_th_n6_hi_sites w_th_n6_hi_br = true ∧",
        "    w_th_n6_lo_br.all validSqrtBracket = true ∧",
        "    w_th_n6_hi_br.all validSqrtBracket = true ∧",
        "    qLt (sumE2 w_th_n6_lo_br) (sumE2 w_th_n6_hi_br) = true ∧",
        "    qLt (sumE1Hi w_th_n6_lo_br) (sumE1Lo w_th_n6_hi_br) = true := by",
        "  decide",
        "",
        f"def w_stereo_p0a : RationalSite3 := ⟨{_lean_int(sp0['p1_nd'][0][0])}, {sp0['p1_nd'][0][1]}, {_lean_int(sp0['p1_nd'][1][0])}, {sp0['p1_nd'][1][1]}, {_lean_int(sp0['p1_nd'][2][0])}, {sp0['p1_nd'][2][1]}⟩",
        f"def w_stereo_p0b : RationalSite3 := ⟨{_lean_int(sp0['p2_nd'][0][0])}, {sp0['p2_nd'][0][1]}, {_lean_int(sp0['p2_nd'][1][0])}, {sp0['p2_nd'][1][1]}, {_lean_int(sp0['p2_nd'][2][0])}, {sp0['p2_nd'][2][1]}⟩",
        f"def w_stereo_p1a : RationalSite3 := ⟨{_lean_int(sp1['p1_nd'][0][0])}, {sp1['p1_nd'][0][1]}, {_lean_int(sp1['p1_nd'][1][0])}, {sp1['p1_nd'][1][1]}, {_lean_int(sp1['p1_nd'][2][0])}, {sp1['p1_nd'][2][1]}⟩",
        f"def w_stereo_p1b : RationalSite3 := ⟨{_lean_int(sp1['p2_nd'][0][0])}, {sp1['p2_nd'][0][1]}, {_lean_int(sp1['p2_nd'][1][0])}, {sp1['p2_nd'][1][1]}, {_lean_int(sp1['p2_nd'][2][0])}, {sp1['p2_nd'][2][1]}⟩",
        f"def w_cayley_q0 : RationalSite4 := ⟨{_lean_int(cs0['q4_nd'][0][0])}, {cs0['q4_nd'][0][1]}, {_lean_int(cs0['q4_nd'][1][0])}, {cs0['q4_nd'][1][1]}, {_lean_int(cs0['q4_nd'][2][0])}, {cs0['q4_nd'][2][1]}, {_lean_int(cs0['q4_nd'][3][0])}, {cs0['q4_nd'][3][1]}⟩",
        f"def w_cayley_q1 : RationalSite4 := ⟨{_lean_int(cs1['q4_nd'][0][0])}, {cs1['q4_nd'][0][1]}, {_lean_int(cs1['q4_nd'][1][0])}, {cs1['q4_nd'][1][1]}, {_lean_int(cs1['q4_nd'][2][0])}, {cs1['q4_nd'][2][1]}, {_lean_int(cs1['q4_nd'][3][0])}, {cs1['q4_nd'][3][1]}⟩",
        "",
        "/-- Compiled rational stereographic S^2(Q), Cayley S^3(Q), and 4-coordinate equatorial bridge theorem. -/",
        "theorem t22_compiled_stereographic_s2q_and_cayley_s3q_bridges :",
        f"    stereographicMatches {_lean_int(sp0['u1_nd'][0])} {sp0['u1_nd'][1]} {_lean_int(sp0['v1_nd'][0])} {sp0['v1_nd'][1]} w_stereo_p0a = true ∧",
        f"    stereographicMatches {_lean_int(sp0['u2_nd'][0])} {sp0['u2_nd'][1]} {_lean_int(sp0['v2_nd'][0])} {sp0['v2_nd'][1]} w_stereo_p0b = true ∧",
        "    siteOnUnitSphere w_stereo_p0a = true ∧",
        "    siteOnUnitSphere w_stereo_p0b = true ∧",
        f"    checkPairDistSq w_stereo_p0a w_stereo_p0b {sp0['chord_nd'][0]} {sp0['chord_nd'][1]} = true ∧",
        f"    stereographicMatches {_lean_int(sp1['u1_nd'][0])} {sp1['u1_nd'][1]} {_lean_int(sp1['v1_nd'][0])} {sp1['v1_nd'][1]} w_stereo_p1a = true ∧",
        f"    stereographicMatches {_lean_int(sp1['u2_nd'][0])} {sp1['u2_nd'][1]} {_lean_int(sp1['v2_nd'][0])} {sp1['v2_nd'][1]} w_stereo_p1b = true ∧",
        "    siteOnUnitSphere w_stereo_p1a = true ∧",
        "    siteOnUnitSphere w_stereo_p1b = true ∧",
        f"    checkPairDistSq w_stereo_p1a w_stereo_p1b {sp1['chord_nd'][0]} {sp1['chord_nd'][1]} = true ∧",
        f"    cayleyS3Matches {_lean_int(cs0['w_nd'][0][0])} {cs0['w_nd'][0][1]} {_lean_int(cs0['w_nd'][1][0])} {cs0['w_nd'][1][1]} {_lean_int(cs0['w_nd'][2][0])} {cs0['w_nd'][2][1]} w_cayley_q0 = true ∧",
        "    siteOnUnitSphere4 w_cayley_q0 = true ∧",
        f"    cayleyS3Matches {_lean_int(cs1['w_nd'][0][0])} {cs1['w_nd'][0][1]} {_lean_int(cs1['w_nd'][1][0])} {cs1['w_nd'][1][1]} {_lean_int(cs1['w_nd'][2][0])} {cs1['w_nd'][2][1]} w_cayley_q1 = true ∧",
        "    siteOnUnitSphere4 w_cayley_q1 = true ∧",
        "    equatorialCayleyMatchesStereo w_cayley_q1 w_stereo_p0a = true := by",
        "  decide",
        "",
        "end PointChargesExactCert",
        "",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Char-Literal-Safe Static Lean 4 Gate + Scrubbed Cold Live Lean & leanchecker
# ---------------------------------------------------------------------------

def _match_lean_char_literal_end(source: str, i: int, n: int) -> int | None:
    """Return end index of a Lean 4 Char literal starting at `source[i] == \"'\"`, or None."""
    if i + 2 < n:
        if source[i + 1] not in ("\\", "'"):
            if source[i + 2] == "'":
                return i + 3
    if i + 3 < n:
        if source[i + 1] == "\\":
            if source[i + 3] == "'":
                if source[i + 2] != "\n":
                    return i + 4
            esc = source[i + 2]
            if esc == "x":
                if i + 5 < n:
                    if source[i + 5] == "'":
                        if all(c in _HEX_DIGITS for c in source[i + 3 : i + 5]):
                            return i + 6
            if esc == "u":
                if i + 7 < n:
                    if source[i + 7] == "'":
                        if all(c in _HEX_DIGITS for c in source[i + 3 : i + 7]):
                            return i + 8
    return None


def strip_lean_comments_and_strings(source: str) -> str:
    """
    Strip Lean 4 `--` line comments, nested `/- ... -/` block comments, single-quoted
    char literals (including `'\"'`, `'\\\\''`, `'\\x22'`, `'\\u0022'`), and double-quoted
    string literals while preserving line numbers. Fails closed on unclosed comments or strings.
    """
    out: list[str] = []
    i, n = 0, len(source)
    depth = 0
    in_ident = False
    while i < n:
        ch = source[i]
        nxt = source[i + 1] if i + 1 < n else ""
        if depth > 0:
            if ch == "/" and nxt == "-":
                depth += 1
                out.append("  ")
                i += 2
            elif ch == "-" and nxt == "/":
                depth -= 1
                out.append("  ")
                i += 2
            else:
                out.append("\n" if ch == "\n" else " ")
                i += 1
            continue
        if ch == "/" and nxt == "-":
            depth = 1
            in_ident = False
            if i + 2 < n and source[i + 2] in ("-", "!"):
                out.append("   ")
                i += 3
            else:
                out.append("  ")
                i += 2
            continue
        if ch == "-" and nxt == "-":
            in_ident = False
            while i < n and source[i] != "\n":
                out.append(" ")
                i += 1
            continue
        if ch == "'" and not in_ident:
            char_end = _match_lean_char_literal_end(source, i, n)
            if char_end is not None:
                out.append(" " * (char_end - i))
                i = char_end
                continue
        if ch == '"':
            in_ident = False
            out.append(" ")
            i += 1
            closed = False
            while i < n:
                sc = source[i]
                if sc == "\\" and i + 1 < n:
                    out.append(" \n" if source[i + 1] == "\n" else "  ")
                    i += 2
                elif sc == '"':
                    out.append(" ")
                    i += 1
                    closed = True
                    break
                else:
                    out.append("\n" if sc == "\n" else " ")
                    i += 1
            require_ok(closed, "LEAN4_UNTERMINATED_STRING", "unterminated string literal in Lean source")
            continue
        in_ident = ch.isalnum() or ch in ("_", "'")
        out.append(ch)
        i += 1
    require_ok(depth == 0, "LEAN4_UNTERMINATED_COMMENT", "unterminated block comment in Lean source")
    return "".join(out)


def static_gate_lean4_source(source: str, *, mode: str = "library") -> None:
    require_ok("\r" not in source, "LEAN4_BAD_CHARS", "Lean source contains CR")
    require_ok("\x00" not in source, "LEAN4_BAD_CHARS", "Lean source contains NUL")
    stripped = strip_lean_comments_and_strings(source)
    m = _LEAN4_BANNED_CORE_RE.search(stripped)
    require_ok(m is None, "LEAN4_TRUST_ESCAPE", f"forbidden Lean token: {m.group(0)!r}" if m else "")
    imports = _LEAN4_IMPORT_RE.findall(stripped)
    if mode == "library":
        require_ok(len(imports) == 0, "LEAN4_FORBIDDEN_IMPORT", f"pure-Init library forbids imports: {imports}")
        for full_name in EXPECTED_THEOREM_NAMES:
            short = full_name.split(".", 1)[1]
            require_ok(re.search(rf"\btheorem\s+{re.escape(short)}\b", stripped) is not None, "LEAN4_MISSING_THEOREM", f"missing {short}")
    else:
        require_ok(imports == [LEAN_MODULE_NAME], "LEAN4_FORBIDDEN_IMPORT", f"AxiomAudit must import only {LEAN_MODULE_NAME}: {imports}")
        for full_name in EXPECTED_THEOREM_NAMES:
            require_ok(re.search(rf"#\s*print\s+axioms\s+{re.escape(full_name)}\b", stripped) is not None, "LEAN4_MISSING_AXIOMS", f"missing #print axioms {full_name}")


def resolve_lean_binaries() -> tuple[str, str, str]:
    lean_bin = os.environ.get("POINT_CHARGES_LEAN_BIN")
    if lean_bin:
        lp = Path(lean_bin)
        require_ok(lp.is_file(), "LEAN_UNAVAILABLE", f"invalid POINT_CHARGES_LEAN_BIN: {lean_bin}")
        require_ok(os.access(lp, os.X_OK), "LEAN_UNAVAILABLE", f"non-executable POINT_CHARGES_LEAN_BIN: {lean_bin}")
    else:
        found = shutil.which("lean")
        if found is None:
            home = os.environ.get("HOME", "")
            cand = Path(home).joinpath(".elan", "bin", "lean") if home else None
            if cand is not None and cand.is_file():
                if os.access(cand, os.X_OK):
                    found = str(cand)
        require_ok(found is not None, "LEAN_UNAVAILABLE", "lean binary not found")
        lp = Path(str(found))
    try:
        pr = subprocess.run(
            [str(lp), "--print-prefix"],
            cwd=str(Path(__file__).resolve().parent),
            capture_output=True, text=True, timeout=30, check=False,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise CertError("LEAN_EXEC_ERROR", f"lean --print-prefix failed: {exc}") from exc
    require_ok(pr.returncode == 0, "LEAN_EXEC_ERROR", "lean --print-prefix nonzero")
    require_ok(bool(pr.stdout.strip()), "LEAN_EXEC_ERROR", "lean --print-prefix empty")
    prefix = Path(pr.stdout.strip())
    real_lean = prefix.joinpath("bin", "lean")
    checker = prefix.joinpath("bin", "leanchecker")
    require_ok(real_lean.is_file(), "LEAN_UNAVAILABLE", "real lean binary missing")
    require_ok(os.access(real_lean, os.X_OK), "LEAN_UNAVAILABLE", "real lean binary not executable")
    require_ok(checker.is_file(), "LEANCHECKER_UNAVAILABLE", "leanchecker binary missing")
    require_ok(os.access(checker, os.X_OK), "LEANCHECKER_UNAVAILABLE", "leanchecker binary not executable")
    return str(real_lean), str(checker), str(prefix)


def verify_lean4_live_and_fresh(project_files: dict[str, bytes]) -> dict[str, Any]:
    """
    Cold-verify the complete pure-Init Lean 4.36.0 package inside an isolated temp workspace
    (copying `.lean`, `lakefile.toml`, `lake-manifest.json`, and `lean-toolchain`), check
    `AxiomAudit.lean` `#print axioms`, and replay with `leanchecker --fresh`.
    """
    lean_raw = project_files[LEAN_FILENAME]
    audit_raw = project_files[AXIOM_AUDIT_FILENAME]
    try:
        lean_src = lean_raw.decode("utf-8", "strict")
        audit_src = audit_raw.decode("utf-8", "strict")
    except UnicodeDecodeError as exc:
        raise CertError("INVALID_UTF8", "Lean file not valid UTF-8") from exc
    static_gate_lean4_source(lean_src, mode="library")
    static_gate_lean4_source(audit_src, mode="audit")

    lean_bin, checker_bin, prefix = resolve_lean_binaries()
    clean_env = {"PATH": f"{Path(lean_bin).parent}:/usr/bin:/bin"}

    try:
        vr = subprocess.run([lean_bin, "--version"], capture_output=True, text=True, timeout=30, env=clean_env, check=False)
    except (OSError, subprocess.SubprocessError) as exc:
        raise CertError("LEAN_EXEC_ERROR", f"lean --version failed: {exc}") from exc
    require_ok(vr.returncode == 0, "LEAN_VERSION_MISMATCH", "lean --version nonzero")
    require_ok("4.36.0" in vr.stdout, "LEAN_VERSION_MISMATCH", f"got {vr.stdout.strip()!r}")

    with tempfile.TemporaryDirectory(prefix="pc_lean_cold_") as td:
        tmp = Path(td)
        for fname, fbytes in project_files.items():
            tmp.joinpath(fname).write_bytes(fbytes)
        olean = tmp.joinpath(f"{LEAN_MODULE_NAME}.olean")
        try:
            comp = subprocess.run(
                [lean_bin, "-R", str(tmp), "-o", str(olean), str(tmp.joinpath(LEAN_FILENAME))],
                capture_output=True, text=True, timeout=180, env=clean_env, check=False,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise CertError("LEAN_EXEC_ERROR", f"lean compile failed: {exc}") from exc
        require_ok(comp.returncode == 0, "LEAN_COMPILE_FAIL", (comp.stdout + comp.stderr)[:200])
        require_ok(comp.stderr == "", "LEAN_STDERR_NONEMPTY", comp.stderr[:200])
        require_ok(comp.stdout == "", "LEAN_WARNING", comp.stdout[:200])

        audit_env = {
            "PATH": f"{Path(lean_bin).parent}:/usr/bin:/bin",
            "LEAN_PATH": f"{tmp}:{prefix}/lib/lean",
        }
        try:
            aud = subprocess.run(
                [lean_bin, "-R", str(tmp), str(tmp.joinpath(AXIOM_AUDIT_FILENAME))],
                capture_output=True, text=True, timeout=180, env=audit_env, check=False,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise CertError("LEAN_EXEC_ERROR", f"lean AxiomAudit failed: {exc}") from exc
        require_ok(aud.returncode == 0, "LEAN_AUDIT_FAIL", (aud.stdout + aud.stderr)[:200])
        require_ok(aud.stderr == "", "LEAN_STDERR_NONEMPTY", aud.stderr[:200])
        require_ok("warning:" not in aud.stdout.lower(), "LEAN_WARNING", aud.stdout[:200])

        folded_lines: list[str] = []
        for raw_line in aud.stdout.splitlines():
            s = raw_line.strip()
            if not s:
                continue
            if raw_line[:1].isspace() and folded_lines:
                folded_lines[-1] = folded_lines[-1] + " " + s
            else:
                folded_lines.append(s)

        seen: dict[str, list[str]] = {}
        for s in folded_lines:
            m_none = _AXIOM_NONE_RE.match(s)
            if m_none is not None:
                nm = m_none.group(1)
                require_ok(nm not in seen, "LEAN_DUP_THM", nm)
                seen[nm] = []
                continue
            m_dep = _AXIOM_DEP_RE.match(s)
            require_ok(m_dep is not None, "LEAN_UNEXPECTED_OUT", s[:100])
            nm = m_dep.group(1)
            require_ok(nm not in seen, "LEAN_DUP_THM", nm)
            axs = [a.strip() for a in m_dep.group(2).split(",") if a.strip()]
            for ax in axs:
                require_ok(ax in _ALLOWED_LEAN_AXIOMS, "LEAN_FORBIDDEN_AXIOM", f"{nm}: {ax}")
            seen[nm] = sorted(axs)
        require_ok(set(seen.keys()) == set(EXPECTED_THEOREM_NAMES), "LEAN_THM_SET_MISMATCH", f"got {len(seen)} theorems")

        try:
            chk = subprocess.run(
                [checker_bin, "--fresh", LEAN_MODULE_NAME],
                capture_output=True, text=True, timeout=180, env=audit_env, check=False,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise CertError("LEANCHECKER_EXEC_ERROR", f"leanchecker failed: {exc}") from exc
        require_ok(chk.returncode == 0, "LEANCHECKER_FAIL", (chk.stdout + chk.stderr)[:200])
        require_ok(chk.stderr == "", "LEANCHECKER_FAIL", chk.stderr[:200])

    return {
        "lean_version": vr.stdout.strip(),
        "lean_source_sha256": sha256_hex(lean_raw),
        "axiom_audit_sha256": sha256_hex(audit_raw),
        "theorems_verified": len(seen),
        "allowed_axioms": sorted(_ALLOWED_LEAN_AXIOMS),
        "leanchecker_fresh_ec": 0,
    }


# ---------------------------------------------------------------------------
# Certificate Construction, Self-Bound Audit, and CLI
# ---------------------------------------------------------------------------

def build_certificate(
    *,
    witness_pack_sha256: str,
    lean_source_sha256: str,
    axiom_audit_sha256: str,
    lakefile_sha256: str,
    manifest_sha256: str,
    toolchain_sha256: str,
    validated: dict[str, Any],
) -> dict[str, Any]:
    payload = {
        "schema_version": SCHEMA_VERSION,
        "upstream_issue": UPSTREAM_ISSUE,
        "oracle_class": "exact_pure_q_plus_lean4_pure_init_compiler",
        "pinned_artifacts": {
            "witness_pack_file": WITNESS_PACK_FILENAME,
            "witness_pack_sha256": witness_pack_sha256,
            "lean_source_file": LEAN_FILENAME,
            "lean_source_sha256": lean_source_sha256,
            "axiom_audit_file": AXIOM_AUDIT_FILENAME,
            "axiom_audit_sha256": axiom_audit_sha256,
            "lakefile_sha256": lakefile_sha256,
            "lake_manifest_sha256": manifest_sha256,
            "lean_toolchain": EXPECTED_LEAN_TOOLCHAIN,
            "lean_toolchain_sha256": toolchain_sha256,
            "theorem_count": len(EXPECTED_THEOREM_NAMES),
            "universal_theorem_count": 14,
            "compiled_witness_theorem_count": 8,
            "theorem_names": list(EXPECTED_THEOREM_NAMES),
            "allowed_axioms": sorted(_ALLOWED_LEAN_AXIOMS),
        },
        "equilibrium_summaries": [
            {
                "name": e["name"],
                "n_charges": e["n_charges"],
                "hessian_diag_q": e["hessian_diag_q"],
                "hessian_trace_q": e["hessian_trace_q"],
                "hessian_det_q": e["hessian_det_q"],
                "morse_index": e["morse_index"],
                "inertia": e["inertia"],
            }
            for e in validated["equilibria"]
        ],
        "off_plane_summary": {
            "name": validated["off_plane"]["name"],
            "z_q": validated["off_plane"]["z_q"],
            "vz_q": validated["off_plane"]["vz_q"],
        },
        "box_collection_summaries": [
            {
                "name": b["name"],
                "n_charges": b["n_charges"],
                "grid_shift_bits": b["grid_shift_bits"],
                "n_boxes": b["n_boxes"],
                "pairwise_disjoint_pairs": b["pairwise_disjoint_pairs"],
                "site_avoidance_pairs": b["site_avoidance_pairs"],
                "morse_counts": b["morse_counts"],
                "m2_minus_m1": b["m2_minus_m1"],
                "maxwell_bound": b["maxwell_bound"],
                "exceeds_maxwell_bound": b["exceeds_maxwell_bound"],
            }
            for b in validated["box_collections"]
        ],
        "thomson_summaries": [
            {
                "name": t["name"],
                "n_charges": t["n_charges"],
                "lower_label": t["lower_label"],
                "lower_e2_q": t["lower"]["e2_q"],
                "lower_e1_hi_q": t["lower"]["e1_hi_q"],
                "higher_label": t["higher_label"],
                "higher_e2_q": t["higher"]["e2_q"],
                "higher_e1_lo_q": t["higher"]["e1_lo_q"],
            }
            for t in validated["thomson_pairs"]
        ],
        "sphere_bridge_summaries": {
            "stereographic_s2q_pairs": [
                {"u1": b["u1"], "v1": b["v1"], "u2": b["u2"], "v2": b["v2"], "chord_sq_q": b["chord_sq_q"]}
                for b in validated["sphere_bridges"]["stereographic_pairs"]
            ],
            "cayley_s3q_points": [
                {"w": h["w"], "q4": h["q4"], "equatorial_matches_stereographic": h["equatorial_matches_stereographic"]}
                for h in validated["sphere_bridges"]["cayley_s3_points"]
            ],
        },
    }
    payload_sha = sha256_hex(canonical_json_bytes(payload))
    env_body = {
        "schema_version": SCHEMA_VERSION,
        "payload_sha256": payload_sha,
        "payload": payload,
    }
    cert_sha = sha256_hex(canonical_json_bytes(env_body))
    return {**env_body, "certificate_sha256": cert_sha}


def _load_and_verify_project_files(root: Path, validated: dict[str, Any]) -> dict[str, bytes]:
    """Read each Lean/Lake artifact once from disk and assert byte-exact parity against deterministic compilation."""
    expected_lean = compile_lean4_source(validated).encode("utf-8")
    expected_audit = compile_axiom_audit_source().encode("utf-8")
    expected_lake = EXPECTED_LAKEFILE_TEXT.encode("ascii")
    expected_manifest = EXPECTED_MANIFEST_TEXT.encode("ascii")
    expected_tc = EXPECTED_TOOLCHAIN_TEXT.encode("ascii")

    disk_lean = read_regular_bytes(root.joinpath(LEAN_FILENAME), LEAN_FILENAME)
    require_ok(disk_lean == expected_lean, "LEAN_DRIFT", "on-disk Lean file differs from compiled witness pack")
    static_gate_lean4_source(disk_lean.decode("utf-8", "strict"), mode="library")

    disk_audit = read_regular_bytes(root.joinpath(AXIOM_AUDIT_FILENAME), AXIOM_AUDIT_FILENAME)
    require_ok(disk_audit == expected_audit, "AUDIT_DRIFT", "on-disk AxiomAudit.lean differs from expected")
    static_gate_lean4_source(disk_audit.decode("utf-8", "strict"), mode="audit")

    disk_lake = read_regular_bytes(root.joinpath(LAKEFILE_FILENAME), LAKEFILE_FILENAME)
    require_ok(disk_lake == expected_lake, "LAKEFILE_DRIFT", "lakefile.toml mismatch")

    disk_manifest = read_regular_bytes(root.joinpath(MANIFEST_FILENAME), MANIFEST_FILENAME)
    require_ok(disk_manifest == expected_manifest, "MANIFEST_DRIFT", "lake-manifest.json mismatch")

    disk_tc = read_regular_bytes(root.joinpath(TOOLCHAIN_FILENAME), TOOLCHAIN_FILENAME)
    require_ok(disk_tc == expected_tc, "TOOLCHAIN_DRIFT", "lean-toolchain mismatch")

    return {
        LEAN_FILENAME: disk_lean,
        AXIOM_AUDIT_FILENAME: disk_audit,
        LAKEFILE_FILENAME: disk_lake,
        MANIFEST_FILENAME: disk_manifest,
        TOOLCHAIN_FILENAME: disk_tc,
    }


def compile_and_verify_bundle(bundle_dir: Path, *, live_lean: bool = True) -> dict[str, Any]:
    wp_raw = read_regular_bytes(bundle_dir.joinpath(WITNESS_PACK_FILENAME), WITNESS_PACK_FILENAME)
    pack = parse_strict_json_bytes(wp_raw, WITNESS_PACK_FILENAME)
    validated = validate_witness_pack(pack)

    lean_src = compile_lean4_source(validated)
    audit_src = compile_axiom_audit_source()
    static_gate_lean4_source(lean_src, mode="library")
    static_gate_lean4_source(audit_src, mode="audit")

    project_files = {
        LEAN_FILENAME: lean_src.encode("utf-8"),
        AXIOM_AUDIT_FILENAME: audit_src.encode("utf-8"),
        LAKEFILE_FILENAME: EXPECTED_LAKEFILE_TEXT.encode("ascii"),
        MANIFEST_FILENAME: EXPECTED_MANIFEST_TEXT.encode("ascii"),
        TOOLCHAIN_FILENAME: EXPECTED_TOOLCHAIN_TEXT.encode("ascii"),
    }
    for fname, fbytes in project_files.items():
        write_regular_atomic_0644(bundle_dir.joinpath(fname), fbytes)

    lean_live = verify_lean4_live_and_fresh(project_files) if live_lean else None
    cert = build_certificate(
        witness_pack_sha256=sha256_hex(wp_raw),
        lean_source_sha256=sha256_hex(project_files[LEAN_FILENAME]),
        axiom_audit_sha256=sha256_hex(project_files[AXIOM_AUDIT_FILENAME]),
        lakefile_sha256=sha256_hex(project_files[LAKEFILE_FILENAME]),
        manifest_sha256=sha256_hex(project_files[MANIFEST_FILENAME]),
        toolchain_sha256=sha256_hex(project_files[TOOLCHAIN_FILENAME]),
        validated=validated,
    )
    cert_bytes = canonical_json_bytes(cert)
    write_regular_atomic_0644(bundle_dir.joinpath(DEFAULT_CERT_FILENAME), cert_bytes)
    return {
        "compiled": True,
        "lean_source_sha256": sha256_hex(project_files[LEAN_FILENAME]),
        "axiom_audit_sha256": sha256_hex(project_files[AXIOM_AUDIT_FILENAME]),
        "certificate_file_sha256": sha256_hex(cert_bytes),
        "certificate_sha256": cert["certificate_sha256"],
        "payload_sha256": cert["payload_sha256"],
        "theorems_compiled": len(EXPECTED_THEOREM_NAMES),
        "live_lean_checked": live_lean,
        "lean_live": lean_live,
    }


def check_certificate_bundle(
    cert_path: Path,
    *,
    bundle_dir: Path | None = None,
    live_lean: bool = True,
) -> dict[str, Any]:
    root = cert_path.parent if bundle_dir is None else bundle_dir
    raw = read_regular_bytes(cert_path, cert_path.name)
    cert = parse_strict_json_bytes(raw, cert_path.name)
    require_ok(
        frozenset(cert.keys()) == frozenset({"schema_version", "payload_sha256", "payload", "certificate_sha256"}),
        "ENVELOPE_KEYS",
        "unexpected certificate envelope keys",
    )
    require_ok(cert["schema_version"] == SCHEMA_VERSION, "SCHEMA_MISMATCH", "envelope schema_version")

    p_sha = require_hex64(cert["payload_sha256"], "payload_sha256")
    c_sha = require_hex64(cert["certificate_sha256"], "certificate_sha256")
    exp_p_sha = sha256_hex(canonical_json_bytes(cert["payload"]))
    require_ok(hmac.compare_digest(p_sha, exp_p_sha), "PAYLOAD_SHA_MISMATCH", "payload_sha256 mismatch")
    env_body = {
        "schema_version": cert["schema_version"],
        "payload_sha256": p_sha,
        "payload": cert["payload"],
    }
    exp_c_sha = sha256_hex(canonical_json_bytes(env_body))
    require_ok(hmac.compare_digest(c_sha, exp_c_sha), "CERT_SHA_MISMATCH", "certificate_sha256 mismatch")

    wp_raw = read_regular_bytes(root.joinpath(WITNESS_PACK_FILENAME), WITNESS_PACK_FILENAME)
    pack = parse_strict_json_bytes(wp_raw, WITNESS_PACK_FILENAME)
    validated = validate_witness_pack(pack)
    project_files = _load_and_verify_project_files(root, validated)

    rebuilt = build_certificate(
        witness_pack_sha256=sha256_hex(wp_raw),
        lean_source_sha256=sha256_hex(project_files[LEAN_FILENAME]),
        axiom_audit_sha256=sha256_hex(project_files[AXIOM_AUDIT_FILENAME]),
        lakefile_sha256=sha256_hex(project_files[LAKEFILE_FILENAME]),
        manifest_sha256=sha256_hex(project_files[MANIFEST_FILENAME]),
        toolchain_sha256=sha256_hex(project_files[TOOLCHAIN_FILENAME]),
        validated=validated,
    )
    require_ok(hmac.compare_digest(canonical_json_bytes(rebuilt), raw), "FULL_REPLAY_MISMATCH", "certificate disagrees with deterministic replay")

    lean_live = verify_lean4_live_and_fresh(project_files) if live_lean else None
    return {
        "valid": True,
        "file_sha256": sha256_hex(raw),
        "certificate_sha256": exp_c_sha,
        "payload_sha256": exp_p_sha,
        "theorems_verified": len(EXPECTED_THEOREM_NAMES),
        "live_lean_checked": live_lean,
        "lean_live": lean_live,
    }


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:]) if argv is None else list(argv)
    root = Path(__file__).resolve().parent
    try:
        mode = "check"
        target = root.joinpath(DEFAULT_CERT_FILENAME)
        want_json = False
        offline = False
        i = 0
        while i < len(args):
            a = args[i]
            if a == "--compile":
                mode = "compile"
                i += 1
            elif a == "--standalone":
                mode = "standalone"
                i += 1
            elif a == "--check":
                mode = "check"
                if i + 1 < len(args) and not args[i + 1].startswith("-"):
                    target = Path(args[i + 1])
                    i += 2
                else:
                    i += 1
            elif a == "--offline":
                offline = True
                i += 1
            elif a == "--json":
                want_json = True
                i += 1
            else:
                raise CertError("UNKNOWN_FLAG", f"unsupported flag: {a!r}")

        if mode == "compile":
            res = compile_and_verify_bundle(root, live_lean=not offline)
            sys.stdout.write(canonical_json_bytes(res).decode("ascii") if want_json else f"COMPILED {res['certificate_sha256']}\n")
            return 0
        if mode == "standalone":
            wp_raw = read_regular_bytes(root.joinpath(WITNESS_PACK_FILENAME), WITNESS_PACK_FILENAME)
            validated = validate_witness_pack(parse_strict_json_bytes(wp_raw, WITNESS_PACK_FILENAME))
            project_files = _load_and_verify_project_files(root, validated)
            cert = build_certificate(
                witness_pack_sha256=sha256_hex(wp_raw),
                lean_source_sha256=sha256_hex(project_files[LEAN_FILENAME]),
                axiom_audit_sha256=sha256_hex(project_files[AXIOM_AUDIT_FILENAME]),
                lakefile_sha256=sha256_hex(project_files[LAKEFILE_FILENAME]),
                manifest_sha256=sha256_hex(project_files[MANIFEST_FILENAME]),
                toolchain_sha256=sha256_hex(project_files[TOOLCHAIN_FILENAME]),
                validated=validated,
            )
            sys.stdout.write(canonical_json_bytes(cert).decode("ascii"))
            return 0
        res = check_certificate_bundle(target, bundle_dir=root, live_lean=not offline)
        sys.stdout.write(canonical_json_bytes(res).decode("ascii") if want_json else f"OK {res['certificate_sha256']}\n")
        return 0
    except CertError as exc:
        err = {"error": True, "code": exc.code, "detail": exc.detail}
        sys.stderr.write(json.dumps(err, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
