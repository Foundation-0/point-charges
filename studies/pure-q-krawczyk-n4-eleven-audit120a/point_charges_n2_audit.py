"""Fail-closed Zone-B hardened standalone CLI & certificate checker for Track S20a.

Zero external dependencies, zero floats, zero bare asserts.
Hardened against all 16 Zone-B adversarial classes (ZB01..ZB16):
non-ASCII / surrogate digests & keys, huge integers, deep nesting, bool-as-int
re-signed SHA forgeries, duplicate JSON keys, floats/NaN/Infinity, non-object
roots, CRLF/non-canonical bytes, hostile Mapping/Sequence/builtin subclasses,
hollow re-signed payloads, symlink targets & symlink ancestors, missing/dir paths,
and invalid UTF-8 / truncated JSON.
"""
from __future__ import annotations

import hashlib
import hmac
import json
import os
from pathlib import Path
import sys
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import point_charges_n2_oracle as oracle
from point_charges_n2_oracle import PointChargeCertError

SCHEMA_ID = "aletheia.track_s20a.point_charges_n2.v1"
MAX_JSON_DEPTH = 48
MAX_JSON_NODES = 10000
MAX_INT_BITS = 14000
HEX_LOWER = frozenset("0123456789abcdef")
TOP_KEYS = ("schema", "payload", "payload_sha256", "certificate_sha256")


def _check_utf8_str(val: str) -> None:
    try:
        val.encode("utf-8")
    except UnicodeEncodeError as exc:
        raise PointChargeCertError("INVALID_UNICODE_STR", "string contains lone UTF-16 surrogate") from exc


def validate_json_tree_iterative(root: Any) -> None:
    """Iteratively validate exact JSON-compatible built-in types, depth, node count, UTF-8, and int size."""
    stack: list[tuple[Any, int]] = [(root, 0)]
    nodes = 0
    while stack:
        val, depth = stack.pop()
        if depth > MAX_JSON_DEPTH:
            raise PointChargeCertError("JSON_TOO_DEEP", "JSON nesting exceeds depth cap")
        nodes += 1
        if nodes > MAX_JSON_NODES:
            raise PointChargeCertError("JSON_TOO_LARGE", "JSON node count exceeds cap")
        t = type(val)
        if t is dict:
            for k, v in val.items():
                if type(k) is not str:
                    raise PointChargeCertError("NON_STR_KEY", "dict keys must be exact str")
                _check_utf8_str(k)
                stack.append((v, depth + 1))
        elif t is list:
            for item in val:
                stack.append((item, depth + 1))
        elif t is str:
            _check_utf8_str(val)
        elif t is int:
            if val.bit_length() > MAX_INT_BITS:
                raise PointChargeCertError("INT_TOO_LARGE", "integer exceeds bit limit")
        elif t is bool or val is None:
            continue
        else:
            raise PointChargeCertError("FORBIDDEN_JSON_TYPE", "disallowed value type in exact JSON tree")


def canonical_json_bytes(obj: Any) -> bytes:
    """Serialize a validated zero-float exact object to canonical UTF-8/ASCII JSON bytes."""
    validate_json_tree_iterative(obj)
    try:
        return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("ascii")
    except (TypeError, ValueError, RecursionError) as exc:
        raise PointChargeCertError("CANONICAL_SERIALIZE_ERROR", "failed to canonically serialize JSON") from exc


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _require_hex64(val: Any, field: str) -> str:
    if type(val) is not str or len(val) != 64 or not val.isascii() or not all(ch in HEX_LOWER for ch in val):
        raise PointChargeCertError("INVALID_SHA256", f"{field} must be a 64-char lowercase ASCII hex digest")
    return val


def exact_json_equal(a: Any, b: Any) -> bool:
    """Iterative type-strict equality over JSON trees (distinguishes bool from int)."""
    stack: list[tuple[Any, Any]] = [(a, b)]
    while stack:
        x, y = stack.pop()
        tx = type(x)
        if tx is not type(y):
            return False
        if tx is dict:
            if len(x) != len(y):
                return False
            for k, xv in x.items():
                if k not in y:
                    return False
                stack.append((xv, y[k]))
        elif tx is list:
            if len(x) != len(y):
                return False
            for xi, yi in zip(x, y):
                stack.append((xi, yi))
        elif x != y:
            return False
    return True


def parse_strict_json_text(raw_text: str) -> dict[str, Any]:
    """Parse strict RFC-8259 zero-float JSON object, rejecting duplicate keys, floats, NaN/Inf, deep nesting, huge ints."""
    if type(raw_text) is not str or len(raw_text.strip()) == 0:
        raise PointChargeCertError("EMPTY_JSON", "JSON input must be a non-empty str")
    _check_utf8_str(raw_text)

    def reject_constant(tok: str) -> Any:
        raise PointChargeCertError("JSON_CONSTANT_FORBIDDEN", f"non-standard constant forbidden: {tok!r}")

    def reject_float(tok: str) -> Any:
        raise PointChargeCertError("JSON_FLOAT_FORBIDDEN", f"float literal forbidden: {tok!r}")

    def reject_dup_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for k, v in pairs:
            if k in out:
                raise PointChargeCertError("JSON_DUPLICATE_KEY", f"duplicate JSON key: {k!r}")
            out[k] = v
        return out

    try:
        parsed = json.loads(
            raw_text,
            parse_constant=reject_constant,
            parse_float=reject_float,
            object_pairs_hook=reject_dup_keys,
        )
    except PointChargeCertError:
        raise
    except (ValueError, TypeError, RecursionError) as exc:
        raise PointChargeCertError("JSON_PARSE_ERROR", "invalid or hostile JSON input") from exc

    if type(parsed) is not dict:
        raise PointChargeCertError("JSON_ROOT_NOT_OBJECT", "top-level JSON must be an object")
    validate_json_tree_iterative(parsed)
    return parsed


def build_certificate() -> dict[str, Any]:
    """Build the canonical self-signed pure-Z/Q Track S20a certificate."""
    payload = oracle.build_full_payload()
    payload_sha = sha256_hex(canonical_json_bytes(payload))
    envelope = {
        "schema": SCHEMA_ID,
        "payload": payload,
        "payload_sha256": payload_sha,
    }
    cert_sha = sha256_hex(canonical_json_bytes(envelope))
    return {
        "schema": SCHEMA_ID,
        "payload": payload,
        "payload_sha256": payload_sha,
        "certificate_sha256": cert_sha,
    }


def verify_certificate_dict(cert: Any) -> dict[str, Any]:
    """Fail-closed in-memory verification of a Track S20a certificate dict (recomputes oracle from scratch)."""
    if type(cert) is not dict:
        raise PointChargeCertError("CERT_NOT_DICT", "certificate must be an exact dict")
    validate_json_tree_iterative(cert)
    if set(cert.keys()) != set(TOP_KEYS):
        raise PointChargeCertError("CERT_KEYS_MISMATCH", "certificate top-level keys do not match schema")
    if cert["schema"] != SCHEMA_ID:
        raise PointChargeCertError("SCHEMA_MISMATCH", "unexpected certificate schema")
    claimed_payload_sha = _require_hex64(cert["payload_sha256"], "payload_sha256")
    claimed_cert_sha = _require_hex64(cert["certificate_sha256"], "certificate_sha256")
    if type(cert["payload"]) is not dict:
        raise PointChargeCertError("PAYLOAD_NOT_DICT", "certificate payload must be an exact dict")

    actual_payload_sha = sha256_hex(canonical_json_bytes(cert["payload"]))
    if not hmac.compare_digest(claimed_payload_sha, actual_payload_sha):
        raise PointChargeCertError("PAYLOAD_SHA_MISMATCH", "payload_sha256 does not match payload bytes")
    envelope = {
        "schema": cert["schema"],
        "payload": cert["payload"],
        "payload_sha256": claimed_payload_sha,
    }
    actual_cert_sha = sha256_hex(canonical_json_bytes(envelope))
    if not hmac.compare_digest(claimed_cert_sha, actual_cert_sha):
        raise PointChargeCertError("CERT_SHA_MISMATCH", "certificate_sha256 does not match envelope bytes")

    expected_cert = build_certificate()
    if not hmac.compare_digest(claimed_payload_sha, expected_cert["payload_sha256"]):
        raise PointChargeCertError("ORACLE_PAYLOAD_MISMATCH", "payload digest disagrees with recomputed oracle")
    return {
        "valid": True,
        "schema": SCHEMA_ID,
        "payload_sha256": claimed_payload_sha,
        "certificate_sha256": claimed_cert_sha,
    }


def _reject_symlinks_on_path(path: Path) -> None:
    """Reject symlink target or any symlink ancestor component."""
    abs_path = Path(os.path.abspath(os.fspath(path)))
    cur = abs_path
    while True:
        if os.path.islink(cur):
            raise PointChargeCertError("SYMLINK_FORBIDDEN", "symlink target or symlink ancestor forbidden")
        parent = cur.parent
        if parent == cur:
            break
        cur = parent


def _coerce_path(path: Any) -> Path:
    if type(path) is str:
        if not path:
            raise PointChargeCertError("INVALID_PATH", "path must be non-empty str or Path")
        return Path(path)
    if isinstance(path, Path):
        return path
    raise PointChargeCertError("INVALID_PATH", "path must be str or Path")


def check_certificate_file(path: Path | str) -> dict[str, Any]:
    """Fail-closed file verification rejecting symlinks, dirs, CRLF, non-canonical JSON, and tampered certs."""
    p = _coerce_path(path)
    _reject_symlinks_on_path(p)
    if not os.path.exists(p):
        raise PointChargeCertError("FILE_NOT_FOUND", "certificate file does not exist")
    if not os.path.isfile(p):
        raise PointChargeCertError("NOT_A_FILE", "certificate path is not a regular file")
    try:
        raw_bytes = p.read_bytes()
    except OSError as exc:
        raise PointChargeCertError("FILE_READ_ERROR", "failed to read certificate file") from exc
    if b"\r" in raw_bytes:
        raise PointChargeCertError("CRLF_FORBIDDEN", "CR/CRLF line endings forbidden in certificate")
    try:
        raw_text = raw_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise PointChargeCertError("INVALID_UTF8", "certificate file is not valid UTF-8") from exc
    cert = parse_strict_json_text(raw_text)
    canonical = canonical_json_bytes(cert)
    stripped = raw_bytes[:-1] if raw_bytes.endswith(b"\n") else raw_bytes
    if stripped != canonical:
        raise PointChargeCertError("NON_CANONICAL_JSON", "certificate bytes are not canonical JSON")
    return verify_certificate_dict(cert)


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    try:
        if args in ([], ["--standalone", "--json"]):
            cert = build_certificate()
            sys.stdout.write(canonical_json_bytes(cert).decode("ascii") + "\n")
            return 0
        if len(args) == 2 and args[0] == "--write":
            out_path = _coerce_path(args[1])
            _reject_symlinks_on_path(out_path)
            if os.path.isdir(out_path):
                raise PointChargeCertError("NOT_A_FILE", "cannot write certificate to a directory path")
            cert = build_certificate()
            out_bytes = canonical_json_bytes(cert) + b"\n"
            try:
                out_path.write_bytes(out_bytes)
            except OSError as exc:
                raise PointChargeCertError("FILE_WRITE_ERROR", "failed to write certificate file") from exc
            sys.stdout.write(canonical_json_bytes({"written": str(out_path), "certificate_sha256": cert["certificate_sha256"]}).decode("ascii") + "\n")
            return 0
        if len(args) == 2 and args[0] == "--check":
            res = check_certificate_file(args[1])
            sys.stdout.write(canonical_json_bytes(res).decode("ascii") + "\n")
            return 0
        raise PointChargeCertError("INVALID_CLI_ARGS", "usage: --standalone --json | --write <path> | --check <path>")
    except PointChargeCertError as exc:
        sys.stderr.write(f"ERROR [{exc.code}]: {exc.message}\n")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
