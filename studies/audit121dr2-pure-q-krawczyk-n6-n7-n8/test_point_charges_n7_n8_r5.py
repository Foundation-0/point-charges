#!/usr/bin/env python3
"""
R5 Pure-Z/Q Adversarial & Mutation-Complete Test Suite for Track S21d
(Point Charges n=7 #12 & n=8 #5 Maxwell Frontier Certifier).

Zero bare `assert`s; independent zero-shared-code exact rational cross-checks;
full Zone-B adversarial hardening coverage across oracle & CLI.
"""
from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import copy
from fractions import Fraction
import hashlib
import inspect
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import point_charges_n7_n8_audit as audit_cli
import point_charges_n7_n8_oracle as oracle


class _IntSub(int):
    pass


class _StrSub(str):
    pass


class _ListSub(list):
    pass


class _DictSub(dict):
    pass


class _FracSub(Fraction):
    pass


def _rehash_cert(cert: dict[str, object]) -> dict[str, object]:
    out = copy.deepcopy(cert)
    payload_sha = oracle.sha256_canonical(out["payload"])
    out["payload_sha256"] = payload_sha
    out["certificate_sha256"] = oracle.sha256_canonical(
        {
            "schema_version": out["schema_version"],
            "track_id": out["track_id"],
            "payload_sha256": payload_sha,
        }
    )
    return out


def _brute_force_coulomb_exact_pythagorean(
    x: tuple[Fraction, Fraction, Fraction],
    sites: list[tuple[Fraction, Fraction, Fraction]],
    charges: list[Fraction],
    dists_exact: list[Fraction],
) -> tuple[list[Fraction], list[list[Fraction]]]:
    """Independent zero-shared-code exact rational Coulomb grad & Hess when ||x - a_j|| in Q."""
    g = [Fraction(0), Fraction(0), Fraction(0)]
    h = [[Fraction(0)] * 3 for _ in range(3)]
    for site, q, r in zip(sites, charges, dists_exact):
        d = [x[i] - site[i] for i in range(3)]
        r2 = r * r
        inv_r3 = Fraction(1) / (r2 * r)
        inv_r5 = inv_r3 / r2
        for i in range(3):
            g[i] -= q * inv_r3 * d[i]
            for j in range(3):
                h[i][j] += 3 * q * inv_r5 * d[i] * d[j]
            h[i][i] -= q * inv_r3
    return g, h


class TestA_ConstantsAndZoneBPrimitives(unittest.TestCase):
    """Fast constant pins and full Zone-B primitive hardening tests."""

    def test_constants_and_defaults_pinned(self) -> None:
        self.assertEqual(oracle.SCHEMA_VERSION, "s21d-point-charges-n7-n8-cert/v1")
        self.assertEqual(oracle.TRACK_ID, "TRACK-S21D")
        self.assertEqual(oracle.PREC_BITS, 384)
        self.assertEqual(oracle.ONE_SCALED, 1 << 384)
        self.assertEqual(oracle.Y_ROUND_BITS, 160)
        self.assertEqual(oracle.MAX_JSON_DEPTH, 64)
        self.assertEqual(oracle.MAX_JSON_BYTES, 2_000_000)
        self.assertEqual(oracle.MAX_INT_BITS, 2048)
        self.assertEqual(oracle.MAX_INT_DIGITS, 620)
        self.assertEqual(oracle.MAX_CONFIGS, 16)
        self.assertEqual(oracle.MAX_CHARGES, 32)
        self.assertEqual(oracle.MAX_BOXES, 256)
        self.assertIs(inspect.signature(oracle.validate_certificate_payload).parameters["full_replay"].default, True)
        self.assertIs(inspect.signature(audit_cli.verify_certificate_file).parameters["full_replay"].default, True)

    def test_require_int_and_bool_zone_b(self) -> None:
        self.assertEqual(oracle.require_int(0, lo=0, hi=0), 0)
        max_ok = (1 << oracle.MAX_INT_BITS) - 1
        self.assertEqual(oracle.require_int(max_ok), max_ok)
        self.assertEqual(oracle.require_int(-max_ok), -max_ok)

        for bad in (True, False, 1.0, "1", _IntSub(1), None):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.require_int(bad)
            self.assertEqual(ctx.exception.code, "INVALID_INT_TYPE")

        for huge in (1 << oracle.MAX_INT_BITS, -(1 << oracle.MAX_INT_BITS), 10**5000):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.require_int(huge)
            self.assertEqual(ctx.exception.code, "INT_BIT_LIMIT")

        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.require_int(4, lo=5)
        self.assertEqual(ctx.exception.code, "INT_BELOW_MIN")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.require_int(6, hi=5)
        self.assertEqual(ctx.exception.code, "INT_ABOVE_MAX")
        for bad_bound in (True, _IntSub(1), "1"):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.require_int(5, lo=bad_bound)  # type: ignore[arg-type]
            self.assertEqual(ctx.exception.code, "INVALID_INT_BOUND")
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.require_int(5, hi=bad_bound)  # type: ignore[arg-type]
            self.assertEqual(ctx.exception.code, "INVALID_INT_BOUND")

        self.assertIs(oracle.require_bool(True), True)
        self.assertIs(oracle.require_bool(False), False)
        for bad in (1, 0, "true", None):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.require_bool(bad)
            self.assertEqual(ctx.exception.code, "INVALID_BOOL_TYPE")

    def test_require_str_seq_mapping_hex_and_fractions(self) -> None:
        self.assertEqual(oracle.require_str("a", min_len=1, max_len=1), "a")
        for bad in (1, _StrSub("a"), b"a"):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.require_str(bad)
            self.assertEqual(ctx.exception.code, "INVALID_STR_TYPE")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.require_str("\ud800")
        self.assertEqual(ctx.exception.code, "INVALID_UTF8")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.require_str("", min_len=1)
        self.assertEqual(ctx.exception.code, "STR_TOO_SHORT")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.require_str("ab", max_len=1)
        self.assertEqual(ctx.exception.code, "STR_TOO_LONG")

        self.assertEqual(oracle.require_seq([1, 2], min_len=2, max_len=2), (1, 2))
        self.assertEqual(oracle.require_seq((1,), min_len=1, max_len=1), (1,))
        for bad in ("ab", {1, 2}, _ListSub([1])):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.require_seq(bad)
            self.assertEqual(ctx.exception.code, "INVALID_SEQ_TYPE")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.require_seq([], min_len=1)
        self.assertEqual(ctx.exception.code, "SEQ_TOO_SHORT")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.require_seq([1, 2], max_len=1)
        self.assertEqual(ctx.exception.code, "SEQ_TOO_LONG")

        self.assertEqual(oracle.require_mapping({"k": 1}), {"k": 1})
        for bad in ([], _DictSub({"a": 1})):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.require_mapping(bad)
            self.assertEqual(ctx.exception.code, "INVALID_MAPPING_TYPE")
        for bad in ({1: "a"}, {_StrSub("k"): 1}):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.require_mapping(bad)
            self.assertEqual(ctx.exception.code, "NON_STR_KEY")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.require_mapping({"\ud800": 1})
        self.assertEqual(ctx.exception.code, "INVALID_UTF8")

        good_hex = "0123456789abcdef" * 4
        other_hex = "fedcba9876543210" * 4
        self.assertEqual(oracle.require_hex64_ascii(good_hex), good_hex)
        self.assertTrue(oracle.constant_time_hex64_equal(good_hex, good_hex))
        self.assertFalse(oracle.constant_time_hex64_equal(good_hex, other_hex))
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.require_hex64_ascii(123)
        self.assertEqual(ctx.exception.code, "INVALID_DIGEST_TYPE")
        for bad_len in ("0" * 63, "0" * 65):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.require_hex64_ascii(bad_len)
            self.assertEqual(ctx.exception.code, "INVALID_DIGEST_LEN")
        for bad_ch in ("A" + "0" * 63, "é" + "0" * 63, "\ud800" + "0" * 63, "g" + "0" * 63):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.require_hex64_ascii(bad_ch)
            self.assertEqual(ctx.exception.code, "INVALID_DIGEST_HEX")

        max_ok = (1 << oracle.MAX_INT_BITS) - 1
        self.assertEqual(oracle.fmt_q(Fraction(7, 1)), "7")
        self.assertEqual(oracle.fmt_q(Fraction(-5, 7)), "-5/7")
        self.assertEqual(oracle.fmt_q(Fraction(max_ok, 1)), str(max_ok))
        self.assertEqual(oracle.fmt_q(Fraction(1, max_ok)), f"1/{max_ok}")
        for bad in (1, "1/2", _FracSub(1, 2)):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.fmt_q(bad)  # type: ignore[arg-type]
            self.assertEqual(ctx.exception.code, "INVALID_FRACTION_TYPE")
        for huge_q in (Fraction(1 << oracle.MAX_INT_BITS, 1), Fraction(1, 1 << oracle.MAX_INT_BITS)):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.fmt_q(huge_q)
            self.assertEqual(ctx.exception.code, "FRACTION_BIT_LIMIT")

        self.assertEqual(oracle.parse_q_str("42"), Fraction(42, 1))
        self.assertEqual(oracle.parse_q_str("-5/7"), Fraction(-5, 7))
        self.assertEqual(oracle.parse_q_str("0" * 619 + "1"), Fraction(1, 1))
        self.assertEqual(oracle.parse_q_str(str(max_ok)), Fraction(max_ok, 1))
        for bad_den in ("1/0", "1/-3"):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.parse_q_str(bad_den)
            self.assertEqual(ctx.exception.code, "INVALID_FRACTION_DEN")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.parse_q_str("1/2/3")
        self.assertEqual(ctx.exception.code, "INVALID_FRACTION_STR")
        for bad_tok in ("-", "١٢٣", "12a", "/2", "1/"):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.parse_q_str(bad_tok)
            self.assertEqual(ctx.exception.code, "INVALID_INT_TOKEN")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.parse_q_str("0" * 620 + "1")
        self.assertEqual(ctx.exception.code, "INT_DIGIT_LIMIT")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.parse_q_str(str(1 << oracle.MAX_INT_BITS))
        self.assertEqual(ctx.exception.code, "INT_BIT_LIMIT")

    def test_json_and_path_zone_b(self) -> None:
        nested_ok: object = 1
        for _ in range(4):
            nested_ok = [nested_ok]
        oracle.validate_json_value(nested_ok, max_depth=4)
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_json_value([nested_ok], max_depth=4)
        self.assertEqual(ctx.exception.code, "MAX_DEPTH_EXCEEDED")

        cyclic: list[object] = []
        cyclic.append(cyclic)
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_json_value(cyclic)
        self.assertEqual(ctx.exception.code, "MAX_DEPTH_EXCEEDED")

        oracle.validate_json_value({"a": [True, False, None, (1 << oracle.MAX_INT_BITS) - 1, "ok"]})
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_json_value({1: "x"})
        self.assertEqual(ctx.exception.code, "NON_STR_KEY")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_json_value({"\ud800": "x"})
        self.assertEqual(ctx.exception.code, "INVALID_UTF8")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_json_value({"a": "\ud800"})
        self.assertEqual(ctx.exception.code, "INVALID_UTF8")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_json_value({"a": 1 << oracle.MAX_INT_BITS})
        self.assertEqual(ctx.exception.code, "INT_BIT_LIMIT")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_json_value({"a": 1.5})
        self.assertEqual(ctx.exception.code, "FLOAT_FORBIDDEN")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_json_value({"a": {1, 2}})
        self.assertEqual(ctx.exception.code, "UNSUPPORTED_JSON_TYPE")

        self.assertEqual(oracle.parse_strict_json_text("[1,2]", require_object=False), [1, 2])
        self.assertEqual(oracle.parse_strict_json_text('{"b":1,"a":2}'), {"b": 1, "a": 2})
        self.assertEqual(oracle.canonical_json_bytes({"b": 1, "a": 2}), b'{"a":2,"b":1}')
        self.assertEqual(
            oracle.sha256_canonical({"b": 1, "a": 2}),
            hashlib.sha256(b'{"a":2,"b":1}').hexdigest(),
        )
        exact_pad = "{}" + (" " * (oracle.MAX_JSON_BYTES - 2))
        self.assertEqual(oracle.parse_strict_json_text(exact_pad), {})
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.parse_strict_json_text(exact_pad + " ")
        self.assertEqual(ctx.exception.code, "JSON_TOO_LARGE")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.parse_strict_json_text(b"{}")  # type: ignore[arg-type]
        self.assertEqual(ctx.exception.code, "INVALID_JSON_INPUT")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.parse_strict_json_text("   ")
        self.assertEqual(ctx.exception.code, "EMPTY_JSON_INPUT")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.parse_strict_json_text('{"x":1.5}')
        self.assertEqual(ctx.exception.code, "JSON_FLOAT_FORBIDDEN")
        for tok in ("NaN", "Infinity", "-Infinity"):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.parse_strict_json_text(f'{{"x":{tok}}}')
            self.assertEqual(ctx.exception.code, "JSON_CONSTANT_FORBIDDEN")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.parse_strict_json_text('{"a":1,"a":2}')
        self.assertEqual(ctx.exception.code, "JSON_DUPLICATE_KEY")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.parse_strict_json_text("[1,2]", require_object=True)
        self.assertEqual(ctx.exception.code, "JSON_ROOT_NOT_OBJECT")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.parse_strict_json_text("[1,2]", require_object=0)  # type: ignore[arg-type]
        self.assertEqual(ctx.exception.code, "INVALID_BOOL_TYPE")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.parse_strict_json_text("{bad")
        self.assertEqual(ctx.exception.code, "JSON_PARSE_ERROR")
        deep_json = "[" * 2000 + "]" * 2000
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.parse_strict_json_text(deep_json, require_object=False)
        self.assertIn(ctx.exception.code, ("MAX_DEPTH_EXCEEDED", "JSON_PARSE_ERROR"))

        class _EvilPath(type(Path())):
            pass

        for bad_path in ("/tmp/x", _EvilPath("/tmp/x")):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.require_safe_path(bad_path)  # type: ignore[arg-type]
            self.assertEqual(ctx.exception.code, "INVALID_PATH_TYPE")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.require_safe_path(Path("bad\x00name"))
        self.assertEqual(ctx.exception.code, "INVALID_PATH_NUL")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.require_safe_path(Path("bad\ud800name"))
        self.assertEqual(ctx.exception.code, "INVALID_UTF8")

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            real_file = root / "ok.json"
            oracle.safe_write_bytes(real_file, b'{"longer_key":123456}\n')
            self.assertEqual(oracle.load_safe_bytes(real_file), b'{"longer_key":123456}\n')
            # Overwriting an existing regular file with shorter bytes truncates cleanly
            oracle.safe_write_bytes(real_file, b'{"a":2}\n')
            self.assertEqual(oracle.load_safe_bytes(real_file), b'{"a":2}\n')
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.safe_write_bytes(real_file, "not-bytes")  # type: ignore[arg-type]
            self.assertEqual(ctx.exception.code, "INVALID_BYTES_TYPE")

            # Hardlinked target must be rejected BEFORE truncation (original content preserved)
            hl_target = root / "hl_orig.json"
            oracle.safe_write_bytes(hl_target, b'{"keep":1}\n')
            hl_alias = root / "hl_alias.json"
            os.link(hl_target, hl_alias)
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.safe_write_bytes(hl_alias, b'{"forged":9}\n')
            self.assertEqual(ctx.exception.code, "HARDLINK_FORBIDDEN")
            self.assertEqual(hl_target.read_bytes(), b'{"keep":1}\n')
            hl_alias.unlink()

            # Non-regular writable device (/dev/null) hits fstat S_ISREG guard in safe_write_bytes
            dev_null = Path("/dev/null")
            if dev_null.exists() and not dev_null.is_symlink():
                with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                    oracle.safe_write_bytes(dev_null, b"{}\n")
                self.assertEqual(ctx.exception.code, "NON_REGULAR_FILE")
                with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                    oracle.load_safe_bytes(dev_null)
                self.assertEqual(ctx.exception.code, "NON_REGULAR_FILE")

            sym_file = root / "sym.json"
            sym_file.symlink_to(real_file)
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.load_safe_bytes(sym_file)
            self.assertEqual(ctx.exception.code, "SYMLINK_FORBIDDEN")

            real_dir = root / "sub"
            real_dir.mkdir()
            deep_dir = real_dir / "nested"
            deep_dir.mkdir()
            (real_dir / "f.json").write_bytes(b'{"k":1}\n')
            sym_dir = root / "sym_sub"
            sym_dir.symlink_to(real_dir)
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.load_safe_bytes(sym_dir / "f.json")
            self.assertEqual(ctx.exception.code, "SYMLINK_ANCESTOR_FORBIDDEN")

            orig_cwd = os.getcwd()
            orig_pwd = os.environ.get("PWD")
            try:
                os.chdir(real_dir)
                os.environ["PWD"] = str(sym_dir)
                with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                    oracle.require_safe_path(Path("f.json"))
                self.assertEqual(ctx.exception.code, "SYMLINK_ANCESTOR_FORBIDDEN")

                os.chdir(deep_dir)
                os.environ["PWD"] = str(sym_dir / "nested")
                with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                    oracle.require_safe_path(Path("f.json"))
                self.assertEqual(ctx.exception.code, "SYMLINK_ANCESTOR_FORBIDDEN")

                os.chdir(root)
                os.environ["PWD"] = str(sym_dir)
                self.assertEqual(oracle.require_safe_path(Path("ok.json")), Path("ok.json"))
            finally:
                os.chdir(orig_cwd)
                if orig_pwd is None:
                    os.environ.pop("PWD", None)
                else:
                    os.environ["PWD"] = orig_pwd

            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.load_safe_bytes(root / "missing.json")
            self.assertEqual(ctx.exception.code, "FILE_READ_ERROR")
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.load_safe_bytes(real_dir)
            self.assertEqual(ctx.exception.code, "NON_REGULAR_FILE")

            fifo_path = root / "pipe.fifo"
            os.mkfifo(fifo_path)
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.load_safe_bytes(fifo_path)
            self.assertEqual(ctx.exception.code, "NON_REGULAR_FILE")
            for bad_out in (real_dir, fifo_path):
                with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                    oracle.safe_write_bytes(bad_out, b"{}\n")
                self.assertEqual(ctx.exception.code, "FILE_WRITE_ERROR")

            big_file = root / "big.bin"
            big_file.write_bytes(b"a" * oracle.MAX_JSON_BYTES)
            self.assertEqual(len(oracle.load_safe_bytes(big_file)), oracle.MAX_JSON_BYTES)
            big_file.write_bytes(b"a" * (oracle.MAX_JSON_BYTES + 1))
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.load_safe_bytes(big_file)
            self.assertEqual(ctx.exception.code, "FILE_TOO_LARGE")

            if os.geteuid() != 0:
                unreadable = root / "unreadable.json"
                unreadable.write_bytes(b"{}\n")
                unreadable.chmod(0o000)
                ro_dir = root / "ro_dir"
                ro_dir.mkdir()
                ro_dir.chmod(0o500)
                try:
                    with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                        oracle.load_safe_bytes(unreadable)
                    self.assertEqual(ctx.exception.code, "FILE_READ_ERROR")
                    with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                        oracle.safe_write_bytes(ro_dir / "x.json", b"{}\n")
                    self.assertEqual(ctx.exception.code, "FILE_WRITE_ERROR")
                finally:
                    unreadable.chmod(0o600)
                    ro_dir.chmod(0o700)


class TestB_PureZQIntervalAndLinearAlgebra(unittest.TestCase):
    """Independent zero-shared-code exact rational checks for interval math & Sylvester inertia."""

    def test_interval_primitives_and_coulomb_enclosure(self) -> None:
        p = oracle.ONE_SCALED
        lo3, hi3 = oracle.iv_from_q(Fraction(1, 3))
        self.assertLessEqual(Fraction(lo3, p), Fraction(1, 3))
        self.assertGreaterEqual(Fraction(hi3, p), Fraction(1, 3))
        self.assertEqual(hi3 - lo3, 1)
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.iv_from_q(1)  # type: ignore[arg-type]
        self.assertEqual(ctx.exception.code, "INVALID_FRACTION_TYPE")

        # iv_sqr across all 3 sign regimes
        self.assertEqual(oracle.iv_sqr((0, 2 * p)), (0, 4 * p))
        self.assertEqual(oracle.iv_sqr((-3 * p, 0)), (0, 9 * p))
        self.assertEqual(oracle.iv_sqr((-2 * p, 3 * p)), (0, 9 * p))
        self.assertEqual(oracle.iv_sqr((-4 * p, 3 * p)), (0, 16 * p))

        # iv_mul_int positive, zero, negative
        self.assertEqual(oracle.iv_mul_int((2 * p, 5 * p), 3), (6 * p, 15 * p))
        self.assertEqual(oracle.iv_mul_int((2 * p, 5 * p), 0), (0, 0))
        self.assertEqual(oracle.iv_mul_int((2 * p, 5 * p), -2), (-10 * p, -4 * p))

        # iv_inv_pos & iv_sqrt_pos exactness and non-positive rejection
        inv2 = oracle.iv_inv_pos((2 * p, 2 * p))
        self.assertEqual(inv2, (p // 2, p // 2))
        for bad_iv in ((0, p), (-p, p)):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.iv_inv_pos(bad_iv)
            self.assertEqual(ctx.exception.code, "NON_POSITIVE_INTERVAL")
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.iv_sqrt_pos(bad_iv)
            self.assertEqual(ctx.exception.code, "NON_POSITIVE_INTERVAL")
        self.assertEqual(oracle.iv_sqrt_pos((4 * p, 4 * p)), (2 * p, 2 * p))
        sq2_lo, sq2_hi = oracle.iv_sqrt_pos((2 * p, 2 * p))
        self.assertLess(Fraction(sq2_lo, p) ** 2, 2)
        self.assertGreater(Fraction(sq2_hi, p) ** 2, 2)

        # Independent Pythagorean Coulomb gradient & Hessian check (exact rational distances 5 and 13)
        x = (Fraction(0), Fraction(0), Fraction(0))
        sites = [
            (Fraction(3), Fraction(4), Fraction(0)),
            (Fraction(-5), Fraction(0), Fraction(12)),
        ]
        charges = [Fraction(7, 2), Fraction(11, 3)]
        g_ref, h_ref = _brute_force_coulomb_exact_pythagorean(x, sites, charges, [Fraction(5), Fraction(13)])
        x_iv = [oracle.iv_from_q(c) for c in x]
        sites_iv = [[oracle.iv_from_q(c) for c in s] for s in sites]
        charges_iv = [oracle.iv_from_q(q) for q in charges]
        g_iv, h_iv = oracle.eval_grad_hess_iv(x_iv, sites_iv, charges_iv)
        for i in range(3):
            self.assertLessEqual(Fraction(g_iv[i][0], p), g_ref[i])
            self.assertGreaterEqual(Fraction(g_iv[i][1], p), g_ref[i])
            for j in range(3):
                self.assertLessEqual(Fraction(h_iv[i][j][0], p), h_ref[i][j])
                self.assertGreaterEqual(Fraction(h_iv[i][j][1], p), h_ref[i][j])
        # Trace of Coulomb Hessian must enclose 0
        tr_iv = oracle.iv_add(oracle.iv_add(h_iv[0][0], h_iv[1][1]), h_iv[2][2])
        self.assertLessEqual(tr_iv[0], 0)
        self.assertGreaterEqual(tr_iv[1], 0)

    def test_bareiss_inverse_and_lagrange_inertia_all_regimes(self) -> None:
        self.assertEqual(oracle.bareiss_det_3x3_int([[2, 1, 0], [1, 3, 1], [0, 1, 2]]), 8)
        self.assertEqual(oracle.bareiss_det_3x3_int([[1, 2, 3], [2, 4, 6], [7, 8, 9]]), 0)
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.bareiss_det_3x3_int([[True, 0, 0], [0, 1, 0], [0, 0, 1]])  # type: ignore[list-item]
        self.assertEqual(ctx.exception.code, "INVALID_INT_TYPE")

        # Pivoting inverse check: M * M^-1 == I_3
        m_piv = [
            [Fraction(0), Fraction(2), Fraction(1)],
            [Fraction(1), Fraction(0), Fraction(3)],
            [Fraction(4), Fraction(5), Fraction(0)],
        ]
        m_inv = oracle.exact_inv_3x3_q(m_piv)
        for i in range(3):
            for j in range(3):
                entry = sum(m_piv[i][k] * m_inv[k][j] for k in range(3))
                self.assertEqual(entry, Fraction(1 if i == j else 0))

        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.exact_inv_3x3_q([[Fraction(1), Fraction(2), Fraction(3)], [Fraction(2), Fraction(4), Fraction(6)], [Fraction(0), Fraction(0), Fraction(1)]])
        self.assertEqual(ctx.exception.code, "SINGULAR_MATRIX")

        # Lagrange Sylvester inertia across ALL regimes (incl. swap pivot, hyperbolic 2x2, degenerate/zero)
        cases = [
            ([[Fraction(3), Fraction(1), Fraction(0)], [Fraction(1), Fraction(3), Fraction(1)], [Fraction(0), Fraction(1), Fraction(3)]], (3, 0, 0)),
            ([[Fraction(-4), Fraction(0), Fraction(0)], [Fraction(0), Fraction(-2), Fraction(0)], [Fraction(0), Fraction(0), Fraction(-1)]], (0, 3, 0)),
            ([[Fraction(0), Fraction(0), Fraction(0)], [Fraction(0), Fraction(2), Fraction(0)], [Fraction(0), Fraction(0), Fraction(-5)]], (1, 1, 1)),
            ([[Fraction(0), Fraction(2), Fraction(0)], [Fraction(2), Fraction(0), Fraction(0)], [Fraction(0), Fraction(0), Fraction(-5)]], (1, 2, 0)),
            ([[Fraction(0), Fraction(0), Fraction(3)], [Fraction(0), Fraction(0), Fraction(0)], [Fraction(3), Fraction(0), Fraction(0)]], (1, 1, 1)),
            ([[Fraction(0), Fraction(0), Fraction(0)], [Fraction(0), Fraction(0), Fraction(0)], [Fraction(0), Fraction(0), Fraction(0)]], (0, 0, 3)),
        ]
        for mat, expected in cases:
            self.assertEqual(oracle.exact_lagrange_inertia_3x3_q(mat), expected)

        # Verify public verify_congruence_diagonalization rejects singular S & non-diagonalizing S
        id_m = [[Fraction(1 if i == j else 0) for j in range(3)] for i in range(3)]
        sing_s = [[Fraction(1), Fraction(0), Fraction(0)], [Fraction(1), Fraction(0), Fraction(0)], [Fraction(0), Fraction(0), Fraction(1)]]
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.verify_congruence_diagonalization(id_m, sing_s)
        self.assertEqual(ctx.exception.code, "SINGULAR_CONGRUENCE")
        non_diag_s = [[Fraction(1), Fraction(1), Fraction(0)], [Fraction(0), Fraction(1), Fraction(0)], [Fraction(0), Fraction(0), Fraction(1)]]
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.verify_congruence_diagonalization(id_m, non_diag_s)
        self.assertEqual(ctx.exception.code, "CONGRUENCE_NOT_DIAGONAL")

        # Reject non-Fraction & non-symmetric on all 3 off-diagonal pairs
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.exact_lagrange_inertia_3x3_q([[1, Fraction(0), Fraction(0)], [Fraction(0), Fraction(1), Fraction(0)], [Fraction(0), Fraction(0), Fraction(1)]])  # type: ignore[list-item]
        self.assertEqual(ctx.exception.code, "INVALID_FRACTION_TYPE")
        for r, c in ((0, 1), (0, 2), (1, 2)):
            bad_m = [[Fraction(1 if i == j else 0) for j in range(3)] for i in range(3)]
            bad_m[r][c] = Fraction(5)
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.exact_lagrange_inertia_3x3_q(bad_m)
            self.assertEqual(ctx.exception.code, "NON_SYMMETRIC_MATRIX")

        # _round_q_dyadic positive & negative
        self.assertEqual(oracle._round_q_dyadic(Fraction(3, 8), 2), Fraction(1, 2))
        self.assertEqual(oracle._round_q_dyadic(Fraction(-3, 8), 2), Fraction(-1, 2))


class TestC_CertifierGuardsAndDefensiveControls(unittest.TestCase):
    """Negative controls covering every Krawczyk, geometry, and defensive invariant guard."""

    def test_certify_configuration_negative_controls(self) -> None:
        s2 = [(Fraction(-1), Fraction(0), Fraction(0)), (Fraction(1), Fraction(0), Fraction(0))]
        q2 = [Fraction(1), Fraction(1)]
        c1 = [(Fraction(0), Fraction(0), Fraction(0))]
        hw1 = [Fraction(1, 1 << 20)]

        # Valid baseline passes
        ok = oracle.certify_configuration("n2-ok", s2, q2, c1, hw1)
        self.assertEqual(ok["morse_counts"], [0, 0, 1, 0])
        self.assertEqual(ok["m2_minus_m1"], 1)

        # Length mismatches
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.certify_configuration("bad", s2, [Fraction(1), Fraction(1), Fraction(1)], c1, hw1)
        self.assertEqual(ctx.exception.code, "SITE_CHARGE_LEN_MISMATCH")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.certify_configuration("bad", s2, q2, c1, hw1 + hw1)
        self.assertEqual(ctx.exception.code, "CENTER_HW_LEN_MISMATCH")

        # Type checks on site, charge, center, hw (including non-comparable str to kill guard_if_reject)
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.certify_configuration("bad", [("0", Fraction(0), Fraction(0)), s2[1]], q2, c1, hw1)  # type: ignore[list-item]
        self.assertEqual(ctx.exception.code, "INVALID_FRACTION_TYPE")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.certify_configuration("bad", s2, ["1", Fraction(1)], c1, hw1)  # type: ignore[list-item]
        self.assertEqual(ctx.exception.code, "INVALID_FRACTION_TYPE")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.certify_configuration("bad", s2, q2, [("0", Fraction(0), Fraction(0))], hw1)  # type: ignore[list-item]
        self.assertEqual(ctx.exception.code, "INVALID_FRACTION_TYPE")
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.certify_configuration("bad", s2, q2, c1, ["1"])  # type: ignore[list-item]
        self.assertEqual(ctx.exception.code, "INVALID_FRACTION_TYPE")

        # Duplicate site, non-positive charge, non-positive half-width
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.certify_configuration("bad", [s2[0], s2[0]], q2, c1, hw1)
        self.assertEqual(ctx.exception.code, "DUPLICATE_CHARGE_SITE")
        for bad_q in (Fraction(0), Fraction(-1)):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.certify_configuration("bad", s2, [bad_q, Fraction(1)], c1, hw1)
            self.assertEqual(ctx.exception.code, "NON_POSITIVE_CHARGE")
        for bad_hw in (Fraction(0), Fraction(-1, 100)):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.certify_configuration("bad", s2, q2, c1, [bad_hw])
            self.assertEqual(ctx.exception.code, "NON_POSITIVE_HALFWIDTH")

        # Touching boxes (u == l) must be rejected as OVERLAPPING_BOXES
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.certify_configuration(
                "touch",
                s2,
                q2,
                [(Fraction(0), Fraction(0), Fraction(0)), (Fraction(2, 100), Fraction(0), Fraction(0))],
                [Fraction(1, 100), Fraction(1, 100)],
            )
        self.assertEqual(ctx.exception.code, "OVERLAPPING_BOXES")

        # Charge inside or on boundary of box rejected
        for ctr in ((Fraction(1), Fraction(0), Fraction(0)), (Fraction(1) - Fraction(1, 100), Fraction(0), Fraction(0))):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.certify_configuration("charge-in", s2, q2, [ctr], [Fraction(1, 100)])
            self.assertEqual(ctx.exception.code, "CHARGE_INSIDE_BOX")

        # Wide box -> KRAWCZYK_NOT_CONTRACTIVE
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.certify_configuration("wide", s2, q2, c1, [Fraction(1, 2)])
        self.assertEqual(ctx.exception.code, "KRAWCZYK_NOT_CONTRACTIVE")

        # Contractive tiny box shifted off-root -> KRAWCZYK_NO_INCLUSION (both + and - shifts!)
        for shift in (Fraction(1, 10), Fraction(-1, 10)):
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.certify_configuration("offroot", s2, q2, [(shift, Fraction(0), Fraction(0))], [Fraction(1, 10000)])
            self.assertEqual(ctx.exception.code, "KRAWCZYK_NO_INCLUSION")

    def test_defensive_invariant_guards_in_certifier(self) -> None:
        s2 = [(Fraction(-1), Fraction(0), Fraction(0)), (Fraction(1), Fraction(0), Fraction(0))]
        q2 = [Fraction(1), Fraction(1)]
        c1 = [(Fraction(0), Fraction(0), Fraction(0))]
        hw1 = [Fraction(1, 1 << 20)]
        p = oracle.ONE_SCALED

        # 1. MIDPOINT_OUTSIDE_BOX_HESS (both lo and hi violations)
        orig_eval = oracle.eval_grad_hess_iv
        for delta in (p, -p):
            calls = [0]

            def fake_eval_mid(box, siv, qiv, _delta=delta, _calls=calls):
                g, h = orig_eval(box, siv, qiv)
                _calls[0] += 1
                if _calls[0] == 2:
                    h[0][0] = (h[0][0][0] + _delta, h[0][0][1] + _delta)
                return g, h

            oracle.eval_grad_hess_iv = fake_eval_mid
            try:
                with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                    oracle.certify_configuration("bad-mid", s2, q2, c1, hw1)
                self.assertEqual(ctx.exception.code, "MIDPOINT_OUTSIDE_BOX_HESS")
            finally:
                oracle.eval_grad_hess_iv = orig_eval

        # 2. HARMONIC_TRACE_VIOLATION (both tr > 0 and tr < 0)
        for tr_shift in (p // 2, -p // 2):
            def fake_eval_tr(box, siv, qiv, _sh=tr_shift):
                g, h = orig_eval(box, siv, qiv)
                h[0][0] = (h[0][0][0] + _sh, h[0][0][1] + _sh)
                return g, h

            oracle.eval_grad_hess_iv = fake_eval_tr
            try:
                with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                    oracle.certify_configuration("bad-tr", s2, q2, c1, hw1)
                self.assertEqual(ctx.exception.code, "HARMONIC_TRACE_VIOLATION")
            finally:
                oracle.eval_grad_hess_iv = orig_eval

        # 3. DEGENERATE_BOX_HESSIAN (det_iv straddling 0, incl. exact 0 endpoints)
        orig_det = oracle.det_3x3_iv
        for fake_d in ((-1, 1), (0, 1), (-1, 0)):
            oracle.det_3x3_iv = lambda _m, _d=fake_d: _d  # type: ignore[assignment]
            try:
                with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                    oracle.certify_configuration("bad-det", s2, q2, c1, hw1)
                self.assertEqual(ctx.exception.code, "DEGENERATE_BOX_HESSIAN")
            finally:
                oracle.det_3x3_iv = orig_det

        # 4. Exact det(Y) == 0 preconditioner guard must reject SINGULAR_PRECONDITIONER
        orig_round = oracle._round_q_dyadic
        oracle._round_q_dyadic = lambda _q, _b: Fraction(0)  # type: ignore[assignment]
        try:
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.certify_configuration("sing-y", s2, q2, c1, hw1)
            self.assertEqual(ctx.exception.code, "SINGULAR_PRECONDITIONER")
        finally:
            oracle._round_q_dyadic = orig_round

        # 5. Exact boundary rho_int == ONE_SCALED must reject KRAWCZYK_NOT_CONTRACTIVE
        orig_abs = oracle.iv_abs_max
        abs_calls = [0]

        def fake_abs_exact_one(_iv, _c=abs_calls):
            _c[0] += 1
            return p if _c[0] == 1 else 0

        oracle.iv_abs_max = fake_abs_exact_one
        try:
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                oracle.certify_configuration("rho-eq-1", s2, q2, c1, hw1)
            self.assertEqual(ctx.exception.code, "KRAWCZYK_NOT_CONTRACTIVE")
        finally:
            oracle.iv_abs_max = orig_abs

        # 5. INVALID_HARMONIC_MORSE_INDEX (n_zero != 0 and n_neg not in (1, 2))
        orig_in = oracle.exact_lagrange_inertia_3x3_q
        for bad_in in ((2, 0, 1), (3, 0, 0), (0, 3, 0)):
            oracle.exact_lagrange_inertia_3x3_q = lambda _m, _bi=bad_in: _bi  # type: ignore[assignment]
            try:
                with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                    oracle.certify_configuration("bad-in", s2, q2, c1, hw1)
                self.assertEqual(ctx.exception.code, "INVALID_HARMONIC_MORSE_INDEX")
            finally:
                oracle.exact_lagrange_inertia_3x3_q = orig_in

        # 6. _config_passes_frontier helper checks all 3 clauses independently
        good_cfg = {"n_certified": 51, "exceeds_maxwell_bound": True, "morse_relation_matches_n_minus_1": True}
        self.assertTrue(oracle._config_passes_frontier(good_cfg, 51))
        self.assertFalse(oracle._config_passes_frontier(good_cfg, 44))
        self.assertFalse(oracle._config_passes_frontier({**good_cfg, "exceeds_maxwell_bound": False}, 51))
        self.assertFalse(oracle._config_passes_frontier({**good_cfg, "morse_relation_matches_n_minus_1": False}, 51))


class TestD_CertificateValidationAndCLI(unittest.TestCase):
    """Full 127-equilibria certificate verification, tamper matrix, and CLI tests."""

    def test_real_certificate_schema_and_all_tamper_rejections(self) -> None:
        cert_path = _HERE / "point_charges_n7_n8_certificate.json"
        raw = cert_path.read_bytes()
        cert = oracle.parse_strict_json_text(raw.decode("utf-8"), require_object=True)

        # Fast schema/consistency verification (no replay)
        res_fast = audit_cli.verify_certificate_file(cert_path, full_replay=False)
        self.assertTrue(res_fast["valid"])
        self.assertEqual(res_fast["total_equilibria_certified"], 127)
        self.assertFalse(res_fast["full_replay"])

        verdicts = cert["payload"]["frontier_verdicts"]
        self.assertTrue(verdicts["issue12_c7_44_certified"])
        self.assertTrue(verdicts["n6_31_new_record_certified"])
        self.assertTrue(verdicts["issue5_n8_51_maxwell_counterexample_certified"])
        self.assertEqual(verdicts["total_equilibria_certified"], 127)

        # Verify exact Morse counts for all 4 configurations
        cfgs = cert["payload"]["configurations"]
        self.assertEqual(cfgs["n2-analytic-control"]["morse_counts"], [0, 0, 1, 0])
        self.assertEqual(cfgs["n7-abk-second-pair"]["morse_counts"], [0, 19, 25, 0])
        self.assertEqual(cfgs["n6-abk-milnor8-flat"]["morse_counts"], [0, 13, 18, 0])
        self.assertEqual(cfgs["n8-abk-milnor8-flat"]["morse_counts"], [0, 22, 29, 0])
        self.assertTrue(cfgs["n8-abk-milnor8-flat"]["exceeds_maxwell_bound"])

        # Non-bool full_replay rejected
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_certificate_payload(cert, full_replay=1)  # type: ignore[arg-type]
        self.assertEqual(ctx.exception.code, "INVALID_BOOL_TYPE")

        # Envelope tampers
        c_extra = copy.deepcopy(cert)
        c_extra["extra"] = 1
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_certificate_payload(c_extra, full_replay=False)
        self.assertEqual(ctx.exception.code, "INVALID_CERT_KEYS")

        c_schema = copy.deepcopy(cert)
        c_schema["schema_version"] = "wrong/v0"
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_certificate_payload(c_schema, full_replay=False)
        self.assertEqual(ctx.exception.code, "INVALID_SCHEMA_VERSION")

        c_track = copy.deepcopy(cert)
        c_track["track_id"] = "WRONG"
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_certificate_payload(c_track, full_replay=False)
        self.assertEqual(ctx.exception.code, "INVALID_TRACK_ID")

        c_psha = copy.deepcopy(cert)
        c_psha["payload_sha256"] = "0" * 64
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_certificate_payload(c_psha, full_replay=False)
        self.assertEqual(ctx.exception.code, "PAYLOAD_SHA_MISMATCH")

        c_csha = copy.deepcopy(cert)
        c_csha["certificate_sha256"] = "0" * 64
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_certificate_payload(c_csha, full_replay=False)
        self.assertEqual(ctx.exception.code, "CERT_SHA_MISMATCH")

        # Payload / config / box / verdict tampers (re-hashed so semantic guards fire!)
        tampers = [
            ("INVALID_PAYLOAD_KEYS", lambda p: p.update({"extra": 1})),
            ("INVALID_PAYLOAD_SCHEMA", lambda p: p.update({"schema_version": "bad"})),
            ("INVALID_PAYLOAD_TRACK", lambda p: p.update({"track_id": "BAD"})),
            ("INVALID_ARITHMETIC_MODE", lambda p: p.update({"arithmetic_mode": "float64"})),
            ("INT_ABOVE_MAX", lambda p: p.update({"prec_bits": 385})),
            ("CONTROLS_FAILED", lambda p: p["known_answer_controls"].update({"all_controls_passed": False})),
            ("CONFIG_SET_MISMATCH", lambda p: p["configurations"].pop("n2-analytic-control")),
            ("INVALID_CONFIG_KEYS", lambda p: p["configurations"]["n2-analytic-control"].update({"extra": 1})),
            ("CONFIG_NAME_MISMATCH", lambda p: p["configurations"]["n2-analytic-control"].update({"name": "other"})),
            ("MAXWELL_BOUND_MISMATCH", lambda p: p["configurations"]["n2-analytic-control"].update({"maxwell_bound": 2})),
            ("EXCEEDS_BOUND_MISMATCH", lambda p: p["configurations"]["n2-analytic-control"].update({"exceeds_maxwell_bound": True})),
            ("MORSE_COUNTS_INCONSISTENT", lambda p: p["configurations"]["n2-analytic-control"].update({"morse_counts": [1, 0, 0, 0]})),
            ("MORSE_COUNTS_INCONSISTENT", lambda p: p["configurations"]["n2-analytic-control"].update({"morse_counts": [0, 0, 0, 1]})),
            ("MORSE_COUNTS_INCONSISTENT", lambda p: p["configurations"]["n2-analytic-control"].update({"morse_counts": [0, 1, 1, 0]})),
            ("MORSE_DIFF_MISMATCH", lambda p: p["configurations"]["n2-analytic-control"].update({"m2_minus_m1": 0})),
            ("MORSE_REL_MISMATCH", lambda p: p["configurations"]["n2-analytic-control"].update({"morse_relation_matches_n_minus_1": False})),
            ("INVALID_MAX_RHO", lambda p: p["configurations"]["n2-analytic-control"].update({"max_rho_q": "-1/10"})),
            ("INVALID_MAX_RHO", lambda p: p["configurations"]["n2-analytic-control"].update({"max_rho_q": "1"})),
            ("NON_POSITIVE_CHARGE", lambda p: p["configurations"]["n2-analytic-control"]["charges"].__setitem__(0, "0")),
            ("INVALID_BOX_KEYS", lambda p: p["configurations"]["n2-analytic-control"]["boxes"][0].update({"x": 1})),
            ("NON_POSITIVE_HALFWIDTH", lambda p: p["configurations"]["n2-analytic-control"]["boxes"][0].update({"half_width": "0"})),
            ("INVALID_BOX_RHO", lambda p: p["configurations"]["n2-analytic-control"]["boxes"][0].update({"rho_upper_q": "-1/10"})),
            ("INVALID_BOX_RHO", lambda p: p["configurations"]["n2-analytic-control"]["boxes"][0].update({"rho_upper_q": "1"})),
            ("BOX_SUMMARY_MISMATCH", lambda p: p["configurations"]["n2-analytic-control"]["boxes"][0].update({"morse_index": 1})),
            ("BOX_SUMMARY_MISMATCH", lambda p: p["configurations"]["n2-analytic-control"]["boxes"][0].update({"rho_upper_q": "0"})),
            ("VERDICT_KEYS_MISMATCH", lambda p: p["frontier_verdicts"].pop("total_equilibria_certified")),
            ("TOTAL_EQ_MISMATCH", lambda p: p["frontier_verdicts"].update({"total_equilibria_certified": 126})),
            ("VERDICT_MISMATCH", lambda p: p["frontier_verdicts"].update({"issue12_c7_44_certified": False})),
            ("VERDICT_MISMATCH", lambda p: p["frontier_verdicts"].update({"n6_31_new_record_certified": False})),
            ("VERDICT_MISMATCH", lambda p: p["frontier_verdicts"].update({"issue5_n8_51_maxwell_counterexample_certified": False})),
        ]
        for expected_code, mutate_fn in tampers:
            bad = copy.deepcopy(cert)
            mutate_fn(bad["payload"])
            bad = _rehash_cert(bad)
            with self.assertRaises(oracle.PointChargeFrontierError, msg=expected_code) as ctx:
                oracle.validate_certificate_payload(bad, full_replay=False)
            self.assertEqual(ctx.exception.code, expected_code)

    def test_cli_modes_errors_and_isolated_subprocess(self) -> None:
        cert_path = _HERE / "point_charges_n7_n8_certificate.json"
        cert = oracle.parse_strict_json_text(cert_path.read_text("utf-8"), require_object=True)

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            non_canon = root / "non_canon.json"
            non_canon.write_text(json.dumps(cert, indent=2) + "\n", encoding="utf-8")
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                audit_cli.verify_certificate_file(non_canon, full_replay=False)
            self.assertEqual(ctx.exception.code, "NON_CANONICAL_CERT_BYTES")

            crlf = root / "crlf.json"
            crlf.write_bytes(oracle.canonical_json_bytes(cert) + b"\r\n")
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                audit_cli.verify_certificate_file(crlf, full_replay=False)
            self.assertEqual(ctx.exception.code, "CR_FORBIDDEN")

            bad_u8 = root / "bad_u8.json"
            bad_u8.write_bytes(b"\xff\xfe")
            with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
                audit_cli.verify_certificate_file(bad_u8, full_replay=False)
            self.assertEqual(ctx.exception.code, "INVALID_UTF8")

            # CLI --check --no-replay success & failure
            out_buf = io.StringIO()
            err_buf = io.StringIO()
            with redirect_stdout(out_buf), redirect_stderr(err_buf):
                rc = audit_cli.main(["--check", str(cert_path), "--no-replay"])
            self.assertEqual(rc, 0)
            self.assertEqual(json.loads(out_buf.getvalue())["file_sha256"], hashlib.sha256(cert_path.read_bytes()).hexdigest())

            out_buf = io.StringIO()
            err_buf = io.StringIO()
            with redirect_stdout(out_buf), redirect_stderr(err_buf):
                rc = audit_cli.main(["--check", str(non_canon), "--no-replay"])
            self.assertEqual(rc, 1)
            self.assertEqual(json.loads(err_buf.getvalue())["error_code"], "NON_CANONICAL_CERT_BYTES")

            # CLI invalid flag combinations and invalid argv types/elements
            for bad_argv, exp_code in (
                (["--standalone", "--no-replay"], "INVALID_CLI_ARGS"),
                (["--write", str(root / "w.json"), "--no-replay"], "INVALID_CLI_ARGS"),
                ("bad-argv", "INVALID_CLI_ARGS"),
                ([123], "INVALID_STR_TYPE"),
                ([_StrSub("--no-replay")], "INVALID_STR_TYPE"),
                (["\ud800"], "INVALID_UTF8"),
            ):
                out_buf = io.StringIO()
                err_buf = io.StringIO()
                with redirect_stdout(out_buf), redirect_stderr(err_buf):
                    rc = audit_cli.main(bad_argv)  # type: ignore[arg-type]
                self.assertEqual(rc, 1)
                self.assertEqual(json.loads(err_buf.getvalue())["error_code"], exp_code)

            # CLI --standalone and --write (fast via temporary stub) + dir/FIFO write rejection
            orig_builder = audit_cli.build_certificate
            audit_cli.build_certificate = lambda: cert  # type: ignore[assignment]
            try:
                out_buf = io.StringIO()
                err_buf = io.StringIO()
                with redirect_stdout(out_buf), redirect_stderr(err_buf):
                    rc = audit_cli.main(["--standalone", "--json"])
                self.assertEqual(rc, 0)
                self.assertEqual(out_buf.getvalue().encode("ascii"), cert_path.read_bytes())

                w_out = root / "written.json"
                out_buf = io.StringIO()
                err_buf = io.StringIO()
                with redirect_stdout(out_buf), redirect_stderr(err_buf):
                    rc = audit_cli.main(["--write", str(w_out)])
                self.assertEqual(rc, 0)
                self.assertTrue(json.loads(out_buf.getvalue())["written"])
                self.assertEqual(w_out.read_bytes(), cert_path.read_bytes())

                fifo = root / "out.fifo"
                os.mkfifo(fifo)
                for bad_out in (root, fifo):
                    out_buf = io.StringIO()
                    err_buf = io.StringIO()
                    with redirect_stdout(out_buf), redirect_stderr(err_buf):
                        rc = audit_cli.main(["--write", str(bad_out)])
                    self.assertEqual(rc, 1)
                    self.assertEqual(json.loads(err_buf.getvalue())["error_code"], "FILE_WRITE_ERROR")
            finally:
                audit_cli.build_certificate = orig_builder

            # Isolated -I -O -B subprocess from foreign cwd kills sys.path & __main__ mutants
            proc = subprocess.run(
                [sys.executable, "-I", "-O", "-B", str(_HERE / "point_charges_n7_n8_audit.py"), "--check", str(cert_path), "--no-replay"],
                cwd=td,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue(json.loads(proc.stdout)["valid"])


class TestZ_FullKrawczykReplay127(unittest.TestCase):
    """End-to-end pure-Z/Q 384-bit Krawczyk + Bareiss replay of all 127 equilibria & replay tampers."""

    def test_z_full_replay_127_and_replay_tampers(self) -> None:
        cert_path = _HERE / "point_charges_n7_n8_certificate.json"
        cert = oracle.parse_strict_json_text(cert_path.read_text("utf-8"), require_object=True)

        # Replay mismatch on box rho forgery (caught by KRAWCZYK_REPLAY_MISMATCH when full_replay=True)
        bad_rho = copy.deepcopy(cert)
        bad_rho["payload"]["configurations"]["n2-analytic-control"]["boxes"][0]["rho_upper_q"] = "1/100"
        bad_rho["payload"]["configurations"]["n2-analytic-control"]["max_rho_q"] = "1/100"
        bad_rho = _rehash_cert(bad_rho)
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_certificate_payload(bad_rho, full_replay=True)
        self.assertEqual(ctx.exception.code, "KRAWCZYK_REPLAY_MISMATCH")

        # Canonical replay mismatch on description tamper (passes per-config Krawczyk, fails CANONICAL_REPLAY_MISMATCH)
        bad_desc = copy.deepcopy(cert)
        bad_desc["payload"]["configurations"]["n2-analytic-control"]["description"] = "forged description"
        bad_desc = _rehash_cert(bad_desc)
        with self.assertRaises(oracle.PointChargeFrontierError) as ctx:
            oracle.validate_certificate_payload(bad_desc, full_replay=True)
        self.assertEqual(ctx.exception.code, "CANONICAL_REPLAY_MISMATCH")

        # Full replay on canonical certificate passes end-to-end (all 127 equilibria re-certified)
        res_full = audit_cli.verify_certificate_file(cert_path, full_replay=True)
        self.assertTrue(res_full["valid"])
        self.assertTrue(res_full["full_replay"])
        self.assertEqual(res_full["total_equilibria_certified"], 127)


if __name__ == "__main__":
    unittest.main()
