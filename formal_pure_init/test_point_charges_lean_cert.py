#!/usr/bin/env python3
"""
Adversarial, tamper-resistance, AST-purity, and mutation-killing test suite
for `formal_pure_init/compile_point_charges_lean_cert.py`.

Pure Python 3 standard library only. Zero bare `assert`s. Zero IEEE-754 floats. Zero skips.
"""
from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import copy
from fractions import Fraction
import io
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

BUNDLE_DIR = Path(__file__).resolve().parent
if str(BUNDLE_DIR) not in sys.path:
    sys.path.insert(0, str(BUNDLE_DIR))

import compile_point_charges_lean_cert as mod  # noqa: E402


def _req(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def _expect_err(code: str, fn, *args, **kwargs) -> mod.CertError:
    try:
        fn(*args, **kwargs)
    except mod.CertError as exc:
        _req(exc.code == code, f"expected CertError({code!r}), got {exc.code!r}: {exc.detail}")
        return exc
    raise AssertionError(f"expected CertError({code!r}) but call succeeded")


def _load_real_pack() -> dict:
    raw = mod.read_regular_bytes(BUNDLE_DIR.joinpath(mod.WITNESS_PACK_FILENAME), mod.WITNESS_PACK_FILENAME)
    return mod.parse_strict_json_bytes(raw, mod.WITNESS_PACK_FILENAME)


def _clean_axiom_stdout() -> str:
    lines: list[str] = []
    for idx, nm in enumerate(mod.EXPECTED_THEOREM_NAMES):
        if idx < 11:
            lines.append(f"'{nm}' depends on axioms: [propext, Quot.sound]")
        elif idx < 14:
            lines.append(f"'{nm}' depends on axioms: [propext,")
            lines.append(" Quot.sound]")
        elif idx == 20:
            lines.append(f"'{nm}' depends on axioms: [propext]")
        else:
            lines.append(f"'{nm}' does not depend on any axioms")
    return "\n".join(lines) + "\n"


def _run_cli(argv: list[str]) -> tuple[int, str, str]:
    out_buf = io.StringIO()
    err_buf = io.StringIO()
    with redirect_stdout(out_buf), redirect_stderr(err_buf):
        rc = mod.main(argv)
    return rc, out_buf.getvalue(), err_buf.getvalue()


class TestPureInitPointChargesLeanCert(unittest.TestCase):
    def test_01_ast_purity_self_check_and_negative_controls(self) -> None:
        """Verify compiler & test suite pass hardened AST guard AND every negative control is caught."""
        for fname in ("compile_point_charges_lean_cert.py", "test_point_charges_lean_cert.py"):
            src = BUNDLE_DIR.joinpath(fname).read_text(encoding="utf-8")
            viols = mod.verify_pure_q_ast(src, label=fname)
            _req(viols == [], f"AST purity violations in {fname}: {viols}")
            _req("/" + "dev/" + "null" not in src, f"null device literal in {fname}")
            _req("/" + "home/" not in src, f"home path leak in {fname}")

        # Positive controls allowed by verify_pure_q_ast
        clean_samples = [
            "import math\nfrom fractions import Fraction\nx = Fraction(1, 2)\ny = math.isqrt(9)\nz = 2 ** 3\ns = f'{z:d}'\n",
            "import fractions as fr\ndef h(x: fr.Fraction) -> fr.Fraction:\n    a: fr.Fraction = fr.Fraction(3, 5)\n    return a\ny = fr.gcd if False else 1\n",
            "class Obj:\n    def log(self, x): return x\n    def Fraction(self, x): return x\nobj = Obj()\nobj.log(1)\nobj.Fraction(1)\n",
        ]
        for idx, sample in enumerate(clean_samples):
            viols = mod.verify_pure_q_ast(sample, label=f"clean_{idx}")
            _req(viols == [], f"false positive on clean sample {idx}: {viols}")

        # Negative controls that MUST trigger verify_pure_q_ast
        bad_samples = [
            ("syntax_error", "def (\n"),
            ("bare_assert", "assert True\n"),
            ("float_lit", "x = 1.5\n"),
            ("complex_lit", "x = 1j\n"),
            ("true_div", "x = 1 / 2\n"),
            ("aug_div", "x = 4\nx /= 2\n"),
            ("floor_div", "x = 1 // 2\n"),
            ("aug_floor_div", "x = 4\nx //= 2\n"),
            ("fstring_float", "x = f'{1:.3f}'\n"),
            ("fstring_sci", "x = f'{1:.2e}'\n"),
            ("neg_pow", "x = 2 ** -1\n"),
            ("expr_pow", "x = 2 ** (0 - 1)\n"),
            ("bool_pow", "x = 2 ** True\n"),
            ("var_pow", "k = -1\nx = 2 ** k\n"),
            ("exp_var_pow", "exp = -1\nx = 2 ** exp\n"),
            ("aug_pow", "x = 2\nx **= 2\n"),
            ("builtin_float", "x = float('1.5')\n"),
            ("builtin_float_alias", "f = float\nx = f(1)\n"),
            ("float_fromhex", "x = float.fromhex('0x1.0p0')\n"),
            ("attr_fromhex", "x = obj.fromhex('0x1.0p0')\n"),
            ("builtin_complex", "x = complex(1, 2)\n"),
            ("builtin_pow", "x = pow(2, -1)\n"),
            ("builtin_pow_float", "x = pow(2, 0.5)\n"),
            ("builtin_pow_3arg", "x = pow(2, -1, None)\n"),
            ("builtin_getattr", "import math\nx = getattr(math, 'sqrt')(2)\n"),
            ("builtin_vars", "import math\nx = vars(math)['sqrt'](2)\n"),
            ("builtin_globals", "x = globals()['float'](1)\n"),
            ("builtin_eval", "x = eval('1/2')\n"),
            ("builtins_ref", "x = __builtins__['float'](1)\n"),
            ("import_builtins", "import builtins\nx = builtins.pow(2, -1)\n"),
            ("truediv_dunder", "x = (1).__truediv__(2)\n"),
            ("pow_dunder", "x = (2).__pow__(-1)\n"),
            ("rpow_dunder", "x = (-1).__rpow__(2)\n"),
            ("float_dunder", "x = (1).__float__()\n"),
            ("complex_dunder", "x = (1).__complex__()\n"),
            ("dict_dunder", "import math\nx = math.__dict__['sqrt'](2)\n"),
            ("sys_modules", "import sys\nx = sys.modules['math']\n"),
            ("relative_import", "from . import helper\n"),
            ("import_decimal", "import decimal\n"),
            ("from_cmath", "from cmath import sqrt\n"),
            ("from_math_log", "from math import log\n"),
            ("import_math_alias_log", "import math as m\nx = m.log(2)\n"),
            ("math_pi_attr", "import math\nx = math.pi\n"),
            ("math_sqrt_attr_alias", "import math\nf = math.sqrt\n"),
            ("from_frac_star", "from fractions import *\n"),
            ("frac_str", "from fractions import Fraction\nx = Fraction('1/2')\n"),
            ("frac_1arg", "from fractions import Fraction\nx = Fraction(1)\n"),
            ("frac_kw", "from fractions import Fraction\nx = Fraction(1, 2, _normalize=False)\n"),
            ("frac_starred", "from fractions import Fraction\nx = Fraction(*a, *b)\n"),
            ("frac_mod_starred", "import fractions\nx = fractions.Fraction(*a, *b)\n"),
            ("frac_alias_1arg", "from fractions import Fraction as F\nx = F(1)\n"),
            ("frac_assign_alias", "from fractions import Fraction\nF = Fraction\nx = F('0.5')\n"),
            ("frac_annassign_alias", "from fractions import Fraction\nF: object = Fraction\nx = F('0.5')\n"),
            ("frac_walrus_alias", "from fractions import Fraction\nx = (F := Fraction)\n"),
            ("frac_ann_walrus_alias", "from fractions import Fraction\nx: (F := Fraction) = 1\n"),
            ("frac_attr_alias", "import fractions\nf = fractions.Fraction\n"),
            ("tuple_unpack_alias", "import fractions, math\nF, G = fractions.Fraction, math\n"),
            ("bfs_order_late_import", "def g():\n    return fr.Fraction(1)\nimport fractions as fr\n"),
            ("math_assign_alias", "import math\nm = math\nx = m.sqrt(2)\n"),
            ("frac_mod_assign_alias", "import fractions\nfr = fractions\nx = fr.Fraction(1)\n"),
            ("frac_mod_1arg", "import fractions\nx = fractions.Fraction(1)\n"),
            ("frac_mod_kw", "import fractions as fr\nx = fr.Fraction(1, 2, _normalize=False)\n"),
            ("frac_from_float", "from fractions import Fraction\nx = Fraction.from_float(1)\n"),
            ("frac_from_decimal", "from fractions import Fraction\nx = Fraction.from_decimal(1)\n"),
        ]
        for label, sample in bad_samples:
            viols = mod.verify_pure_q_ast(sample, label=label)
            _req(len(viols) >= 1, f"negative AST control {label!r} missed")
        syn_v = mod.verify_pure_q_ast("def (\n", label="bad_syn")
        _req(len(syn_v) == 1 and "syntax error" in syn_v[0], f"bad syntax violation format: {syn_v}")

    def test_02_exact_primitives_parse_q_int_hex_strict_json_and_io(self) -> None:
        """Exhaustive positive & refusal coverage for exact parsers and strict file IO."""
        _req(mod.parse_q("0", "q") == Fraction(0, 1), "zero")
        _req(mod.parse_q("-42", "q") == Fraction(-42, 1), "neg int")
        _req(mod.parse_q("3/7", "q") == Fraction(3, 7), "reduced frac")
        _req(mod.parse_q("-5/12", "q") == Fraction(-5, 12), "neg frac")
        _req(mod.fmt_q(Fraction(7, 1)) == "7", "fmt int")
        _req(mod.fmt_q(Fraction(-3, 5)) == "-3/5", "fmt frac")
        _req(mod._q_div(Fraction(3, 4), Fraction(2, 5)) == Fraction(15, 8), "q_div")
        _expect_err("ZERO_DIVISOR", mod._q_div, Fraction(1, 1), Fraction(0, 1))
        _req(mod._q_pow_pos(Fraction(2, 3), 1) == Fraction(2, 3), "q_pow_pos 1")
        _expect_err("INVALID_EXP", mod._q_pow_pos, Fraction(2, 3), 0)

        big_ok = "9" * mod.MAX_INT_DIGITS
        big_bad = "9" * (mod.MAX_INT_DIGITS + 1)
        _req(mod.parse_q(big_ok, "b") == Fraction(int(big_ok), 1), "max digits int")
        _req(mod.parse_q(f"-{big_ok}", "b") == Fraction(-int(big_ok), 1), "max digits neg int")
        _req(mod.parse_q(f"{big_ok}/2", "b") == Fraction(int(big_ok), 2), "max digits num")
        _req(mod.parse_q(f"1/{big_ok}", "b") == Fraction(1, int(big_ok)), "max digits den")
        _expect_err("INT_TOO_LARGE", mod.parse_q, big_bad, "b")
        _expect_err("INT_TOO_LARGE", mod.parse_q, f"{big_bad}/2", "b")
        _expect_err("INT_TOO_LARGE", mod.parse_q, f"1/{big_bad}", "b")

        for bad_val in (1, True, None, ["1"]):
            _expect_err("INVALID_Q_TYPE", mod.parse_q, bad_val, "q")
        for bad_fmt in ("", " 1", "01", "-0", "+1", "1.5", "1/0", "1/-2", "1/02", "abc"):
            _expect_err("INVALID_Q_FORMAT", mod.parse_q, bad_fmt, "q")
        _expect_err("NON_CANONICAL_Q", mod.parse_q, "3/1", "q")
        _expect_err("NON_REDUCED_Q", mod.parse_q, "2/4", "q")
        _expect_err("NON_REDUCED_Q", mod.parse_q, "0/5", "q")

        _req(mod.require_exact_int(5, "x", min_val=5, max_val=5) == 5, "exact boundary int")
        _expect_err("INVALID_INT", mod.require_exact_int, True, "x")
        _expect_err("INVALID_INT", mod.require_exact_int, "5", "x")
        _expect_err("INT_OUT_OF_RANGE", mod.require_exact_int, 4, "x", min_val=5)
        _expect_err("INT_OUT_OF_RANGE", mod.require_exact_int, 6, "x", max_val=5)
        h64 = "0123456789abcdef" * 4
        _req(mod.require_hex64(h64, "h") == h64, "valid hex64")
        for bad_h in (123, h64[:-1], h64 + "0", h64.upper(), "g" * 64):
            _expect_err("INVALID_HEX64", mod.require_hex64, bad_h, "h")

        # parse_strict_json_bytes
        _req(mod.parse_strict_json_bytes(b'{"a":1}\n', "j") == {"a": 1}, "valid json")
        exact_json = b'{"a":"' + (b"x" * (mod.MAX_FILE_BYTES - 8)) + b'"}'
        _req(len(mod.parse_strict_json_bytes(exact_json, "j")["a"]) == mod.MAX_FILE_BYTES - 8, "exact max json")
        _expect_err("BAD_FILE_SIZE", mod.parse_strict_json_bytes, b"", "j")
        _expect_err("BAD_FILE_SIZE", mod.parse_strict_json_bytes, b" " * (mod.MAX_FILE_BYTES + 1), "j")
        _expect_err("BAD_BYTES", mod.parse_strict_json_bytes, b'{"a":1}\x00', "j")
        _expect_err("BAD_BYTES", mod.parse_strict_json_bytes, b'{"a":1}\r\n', "j")
        _expect_err("NON_ASCII_JSON", mod.parse_strict_json_bytes, '{"a":"\u00e9"}'.encode("utf-8"), "j")
        _expect_err("FLOAT_FORBIDDEN", mod.parse_strict_json_bytes, b'{"a":1.5}\n', "j")
        _expect_err("JSON_CONSTANT_FORBIDDEN", mod.parse_strict_json_bytes, b'{"a":NaN}\n', "j")
        _expect_err("DUPLICATE_JSON_KEY", mod.parse_strict_json_bytes, b'{"a":1,"a":2}\n', "j")
        _expect_err("INVALID_JSON", mod.parse_strict_json_bytes, b"{bad", "j")
        _expect_err("JSON_NOT_OBJECT", mod.parse_strict_json_bytes, b"[1,2]\n", "j")

        # read_regular_bytes & write_regular_atomic_0644
        with tempfile.TemporaryDirectory(prefix="pc_io_") as td:
            root = Path(td)
            f_ok = root.joinpath("sample.bin")
            mod.write_regular_atomic_0644(f_ok, b"hello")
            _req(mod.read_regular_bytes(f_ok, "s") == b"hello", "roundtrip")
            _req((os.lstat(f_ok).st_mode & 0o777) == 0o644, "0644 mode")

            exact_f = root.joinpath("exact.bin")
            exact_f.write_bytes(b"y" * mod.MAX_FILE_BYTES)
            _req(len(mod.read_regular_bytes(exact_f, "ex")) == mod.MAX_FILE_BYTES, "exact max file")

            _expect_err("PARENT_NOT_DIR", mod.write_regular_atomic_0644, root.joinpath("no_dir", "x"), b"a")
            _expect_err("IO_ERROR", mod.read_regular_bytes, root.joinpath("missing.bin"), "m")
            d_bad = root.joinpath("subdir")
            d_bad.mkdir()
            _expect_err("NON_REGULAR_FILE", mod.read_regular_bytes, d_bad, "d")
            s_bad = root.joinpath("link.bin")
            os.symlink(f_ok, s_bad)
            _expect_err("NON_REGULAR_FILE", mod.read_regular_bytes, s_bad, "sym")
            if hasattr(os, "mkfifo"):
                p_bad = root.joinpath("pipe.bin")
                os.mkfifo(p_bad)
                _expect_err("NON_REGULAR_FILE", mod.read_regular_bytes, p_bad, "fifo")
            big_f = root.joinpath("big.bin")
            big_f.write_bytes(b"x" * (mod.MAX_FILE_BYTES + 1))
            _expect_err("BAD_FILE_SIZE", mod.read_regular_bytes, big_f, "big")
            real_st = os.lstat(f_ok)
            fake_big_st = os.stat_result((real_st.st_mode, 0, 0, 1, 0, 0, mod.MAX_FILE_BYTES + 1, 0, 0, 0))
            with mock.patch.object(os, "lstat", return_value=fake_big_st):
                _expect_err("BAD_FILE_SIZE", mod.read_regular_bytes, f_ok, "lstat_cap")
            fake_fifo_st = os.stat_result((stat.S_IFIFO | 0o600, 0, 0, 1, 0, 0, 5, 0, 0, 0))
            with mock.patch.object(os, "fstat", return_value=fake_fifo_st):
                _expect_err("NON_REGULAR_FILE", mod.read_regular_bytes, f_ok, "fstat_fifo")
            fake_small_st = os.stat_result((real_st.st_mode, 0, 0, 1, 0, 0, 1, 0, 0, 0))
            with mock.patch.object(os, "lstat", return_value=fake_small_st):
                _expect_err("BAD_FILE_SIZE", mod.read_regular_bytes, big_f, "stream_cap")
            with mock.patch.object(os, "replace", side_effect=OSError("boom")):
                try:
                    mod.write_regular_atomic_0644(root.joinpath("fail.bin"), b"abc")
                    raise AssertionError("expected OSError")
                except OSError:
                    pass
            _req(not any(p.name.startswith(".pc_lean_") for p in root.iterdir()), "temp cleaned")

    def test_03_equilibrium_and_off_plane_witness_validation_and_refusals(self) -> None:
        """Positive + every refusal branch on equilibrium & off-plane confinement validators."""
        pack = _load_real_pack()
        eq0 = copy.deepcopy(pack["equilibrium_witnesses"][0])
        res0 = mod.validate_equilibrium_witness(eq0)
        _req(res0["hessian_diag_q"] == ["4", "-2", "-2"], f"bad n2 diag: {res0}")
        _req(res0["morse_index"] == 2, f"bad n2 morse: {res0}")

        _expect_err("INVALID_WITNESS", mod.validate_equilibrium_witness, "not_dict")
        bad = copy.deepcopy(eq0); bad["name"] = 123
        _expect_err("INVALID_WITNESS", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["name"] = "  "
        _expect_err("INVALID_WITNESS", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["sites"] = tuple(bad["sites"])
        _expect_err("INVALID_WITNESS", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["charges"] = tuple(bad["charges"])
        _expect_err("INVALID_WITNESS", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["radii"] = tuple(bad["radii"])
        _expect_err("INVALID_WITNESS", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["sites"] = [["1", "0", "0"]]; bad["charges"] = ["1"]; bad["radii"] = ["1"]
        _expect_err("INVALID_WITNESS", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["charges"] = ["1"]
        _expect_err("INVALID_WITNESS", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["radii"] = ["1"]
        _expect_err("INVALID_WITNESS", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["sites"][0] = ("1", "0", "0")
        _expect_err("INVALID_VEC3", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["sites"][0] = ["1", "0"]
        _expect_err("INVALID_VEC3", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["sites"][1] = list(bad["sites"][0])
        _expect_err("DUPLICATE_SITE", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["charges"][0] = "0"
        _expect_err("NONPOSITIVE_CHARGE", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["radii"][0] = "-1"
        _expect_err("NONPOSITIVE_RADIUS", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["expected_morse_index"] = 3
        _expect_err("INT_OUT_OF_RANGE", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["radii"][0] = "2"
        _expect_err("RADIUS_SQ_MISMATCH", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["charges"][0] = "2"
        _expect_err("NONZERO_GRADIENT", mod.validate_equilibrium_witness, bad)
        bad = copy.deepcopy(eq0); bad["expected_morse_index"] = 1
        _expect_err("MORSE_MISMATCH", mod.validate_equilibrium_witness, bad)

        oct_deg = {
            "name": "oct_degenerate",
            "sites": [["1", "0", "0"], ["-1", "0", "0"], ["0", "1", "0"], ["0", "-1", "0"], ["0", "0", "1"], ["0", "0", "-1"]],
            "charges": ["1", "1", "1", "1", "1", "1"],
            "point": ["0", "0", "0"],
            "radii": ["1", "1", "1", "1", "1", "1"],
            "expected_morse_index": 2,
        }
        _expect_err("DEGENERATE_HESSIAN", mod.validate_equilibrium_witness, oct_deg)

        # Rotated n=2 equilibria exercising each off-diagonal H_xy, H_xz, H_yz refusal individually
        for label, s0, s1 in (
            ("rot_xy", ["3/5", "4/5", "0"], ["-3/5", "-4/5", "0"]),
            ("rot_xz", ["3/5", "0", "4/5"], ["-3/5", "0", "-4/5"]),
            ("rot_yz", ["0", "3/5", "4/5"], ["0", "-3/5", "-4/5"]),
        ):
            rot = {
                "name": label,
                "sites": [s0, s1],
                "charges": ["1", "1"],
                "point": ["0", "0", "0"],
                "radii": ["1", "1"],
                "expected_morse_index": 2,
            }
            _expect_err("NONDIAGONAL_WITNESS", mod.validate_equilibrium_witness, rot)

        # Off-plane confinement witness
        off0 = copy.deepcopy(pack["off_plane_confinement_witness"])
        ro = mod.validate_off_plane_witness(off0)
        _req(ro["z_q"] == "4", f"bad off z: {ro}")
        _req(ro["vz_q"] == "-12/125", f"bad off vz: {ro}")
        off_neg = copy.deepcopy(off0); off_neg["point"][2] = "-4"
        _req(mod.validate_off_plane_witness(off_neg)["vz_q"] == "12/125", "bad neg off_plane")
        # Exact 1-site boundary for off_plane
        off_one = {"name": "one_site", "sites": [["0", "0", "0"]], "charges": ["1"], "point": ["0", "0", "4"], "radii": ["4"]}
        _req(mod.validate_off_plane_witness(off_one)["vz_q"] == "-1/16", "1-site off_plane")

        _expect_err("INVALID_WITNESS", mod.validate_off_plane_witness, [])
        bad = copy.deepcopy(off0); bad["sites"] = tuple(bad["sites"])
        _expect_err("INVALID_WITNESS", mod.validate_off_plane_witness, bad)
        bad = copy.deepcopy(off0); bad["charges"] = tuple(bad["charges"])
        _expect_err("INVALID_WITNESS", mod.validate_off_plane_witness, bad)
        bad = copy.deepcopy(off0); bad["radii"] = tuple(bad["radii"])
        _expect_err("INVALID_WITNESS", mod.validate_off_plane_witness, bad)
        bad = copy.deepcopy(off0); bad["sites"] = []; bad["charges"] = []; bad["radii"] = []
        _expect_err("INVALID_WITNESS", mod.validate_off_plane_witness, bad)
        bad = copy.deepcopy(off0); bad["charges"] = ["1"]
        _expect_err("INVALID_WITNESS", mod.validate_off_plane_witness, bad)
        bad = copy.deepcopy(off0); bad["radii"] = ["5"]
        _expect_err("INVALID_WITNESS", mod.validate_off_plane_witness, bad)
        bad = copy.deepcopy(off0); bad["sites"][0][2] = "1"
        _expect_err("SITE_NOT_PLANAR", mod.validate_off_plane_witness, bad)
        bad = copy.deepcopy(off0); bad["charges"][0] = "0"
        _expect_err("NONPOSITIVE_CHARGE", mod.validate_off_plane_witness, bad)
        bad = copy.deepcopy(off0); bad["point"][2] = "0"
        _expect_err("POINT_ON_PLANE", mod.validate_off_plane_witness, bad)
        bad = copy.deepcopy(off0); bad["radii"][0] = "0"
        _expect_err("NONPOSITIVE_RADIUS", mod.validate_off_plane_witness, bad)
        bad = copy.deepcopy(off0); bad["radii"][0] = "6"
        _expect_err("RADIUS_SQ_MISMATCH", mod.validate_off_plane_witness, bad)

    def test_04_box_collections_thomson_s2q_and_sphere_bridges(self) -> None:
        """Positive + every refusal branch on box collections, Thomson S^2(Q), and sphere bridges."""
        pack = _load_real_pack()
        box0 = copy.deepcopy(pack["box_collection_witnesses"][0])
        rb = mod.validate_box_collection_witness(box0)
        _req(rb["n_boxes"] == 11, f"bad n_boxes: {rb}")
        _req(rb["exceeds_maxwell_bound"], "expected exceeds_maxwell_bound")

        _expect_err("INVALID_WITNESS", mod.validate_box_collection_witness, None)
        bad = copy.deepcopy(box0); bad["name"] = 1
        _expect_err("INVALID_WITNESS", mod.validate_box_collection_witness, bad)
        bad = copy.deepcopy(box0); bad["name"] = ""
        _expect_err("INVALID_WITNESS", mod.validate_box_collection_witness, bad)
        bad = copy.deepcopy(box0); bad["n_charges"] = 1
        _expect_err("INT_OUT_OF_RANGE", mod.validate_box_collection_witness, bad)
        bad = copy.deepcopy(box0); bad["grid_shift_bits"] = 4
        _expect_err("INT_OUT_OF_RANGE", mod.validate_box_collection_witness, bad)
        bad = copy.deepcopy(box0); bad["sites"] = tuple(bad["sites"])
        _expect_err("INVALID_WITNESS", mod.validate_box_collection_witness, bad)
        bad = copy.deepcopy(box0); bad["sites"] = bad["sites"][:-1]
        _expect_err("INVALID_WITNESS", mod.validate_box_collection_witness, bad)
        bad = copy.deepcopy(box0); bad["charges"] = tuple(bad["charges"])
        _expect_err("INVALID_WITNESS", mod.validate_box_collection_witness, bad)
        bad = copy.deepcopy(box0); bad["charges"] = bad["charges"][:-1]
        _expect_err("INVALID_WITNESS", mod.validate_box_collection_witness, bad)
        bad = copy.deepcopy(box0); bad["boxes"] = tuple(bad["boxes"])
        _expect_err("INVALID_WITNESS", mod.validate_box_collection_witness, bad)
        bad = copy.deepcopy(box0); bad["boxes"] = []
        _expect_err("INVALID_WITNESS", mod.validate_box_collection_witness, bad)
        bad = copy.deepcopy(box0); bad["sites"][1] = list(bad["sites"][0])
        _expect_err("DUPLICATE_SITE", mod.validate_box_collection_witness, bad)
        bad = copy.deepcopy(box0); bad["charges"][0] = "-1"
        _expect_err("NONPOSITIVE_CHARGE", mod.validate_box_collection_witness, bad)
        bad = copy.deepcopy(box0); bad["boxes"][0] = "not_dict"
        _expect_err("INVALID_BOX", mod.validate_box_collection_witness, bad)
        bad = copy.deepcopy(box0); bad["boxes"][0]["half_width"] = "0"
        _expect_err("NONPOSITIVE_HALF_WIDTH", mod.validate_box_collection_witness, bad)
        bad = copy.deepcopy(box0); bad["boxes"][0]["morse_index"] = 3
        _expect_err("INT_OUT_OF_RANGE", mod.validate_box_collection_witness, bad)

        tiny_box = {
            "name": "touch_test",
            "n_charges": 2,
            "grid_shift_bits": 8,
            "sites": [["10", "0", "0"], ["-10", "0", "0"]],
            "charges": ["1", "1"],
            "boxes": [
                {"center": ["0", "0", "0"], "half_width": "1", "morse_index": 2},
                {"center": ["2", "0", "0"], "half_width": "1", "morse_index": 1},
            ],
        }
        _expect_err("OVERLAPPING_BOXES", mod.validate_box_collection_witness, tiny_box)
        site_in = copy.deepcopy(tiny_box)
        site_in["boxes"] = [{"center": ["10", "0", "0"], "half_width": "1", "morse_index": 2}]
        _expect_err("SITE_IN_BOX", mod.validate_box_collection_witness, site_in)
        ph_bad = copy.deepcopy(box0)
        ph_bad["boxes"][0]["morse_index"] = 1 if ph_bad["boxes"][0]["morse_index"] == 2 else 2
        _expect_err("POINCARE_HOPF_MISMATCH", mod.validate_box_collection_witness, ph_bad)

        # Thomson S^2(Q)
        _expect_err("BRACKET_UNDERFLOW", mod._sqrt_bracket_q, Fraction(1, 10000), 4)
        _expect_err("INVALID_THOMSON_SITES", mod._eval_thomson_s2_config, [(Fraction(1, 1), Fraction(0, 1), Fraction(0, 1))], 8, "one")
        ev2 = mod._eval_thomson_s2_config(
            [(Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)), (Fraction(-1, 1), Fraction(0, 1), Fraction(0, 1))],
            8,
            "two",
        )
        _req(ev2["e2"] == Fraction(1, 4), "2-site antipodal E2=1/4")
        th0 = copy.deepcopy(pack["thomson_pair_witnesses"][0])
        rt = mod.validate_thomson_pair_witness(th0)
        _req(rt["lower"]["e2_q"] == "83825/36864", f"bad thomson lo: {rt}")
        _req(rt["higher"]["e2_q"] == "5/2", f"bad thomson hi: {rt}")

        _expect_err("INVALID_WITNESS", mod.validate_thomson_pair_witness, 123)
        bad = copy.deepcopy(th0); bad["name"] = 1
        _expect_err("INVALID_WITNESS", mod.validate_thomson_pair_witness, bad)
        bad = copy.deepcopy(th0); bad["name"] = ""
        _expect_err("INVALID_WITNESS", mod.validate_thomson_pair_witness, bad)
        bad = copy.deepcopy(th0); bad["n_charges"] = 1
        _expect_err("INT_OUT_OF_RANGE", mod.validate_thomson_pair_witness, bad)
        bad = copy.deepcopy(th0); bad["sqrt_bracket_bits"] = 2
        _expect_err("INT_OUT_OF_RANGE", mod.validate_thomson_pair_witness, bad)
        bad = copy.deepcopy(th0); bad["lower_energy_sites"] = tuple(bad["lower_energy_sites"])
        _expect_err("INVALID_WITNESS", mod.validate_thomson_pair_witness, bad)
        bad = copy.deepcopy(th0); bad["higher_energy_sites"] = tuple(bad["higher_energy_sites"])
        _expect_err("INVALID_WITNESS", mod.validate_thomson_pair_witness, bad)
        bad = copy.deepcopy(th0); bad["lower_energy_sites"] = bad["lower_energy_sites"][:-1]
        _expect_err("INVALID_WITNESS", mod.validate_thomson_pair_witness, bad)
        bad = copy.deepcopy(th0); bad["higher_energy_sites"] = bad["higher_energy_sites"][:-1]
        _expect_err("INVALID_WITNESS", mod.validate_thomson_pair_witness, bad)
        bad = copy.deepcopy(th0); bad["lower_energy_sites"][1] = list(bad["lower_energy_sites"][0])
        _expect_err("INVALID_THOMSON_SITES", mod.validate_thomson_pair_witness, bad)
        bad = copy.deepcopy(th0); bad["lower_energy_sites"][0] = ["2", "0", "0"]
        _expect_err("NOT_ON_UNIT_SPHERE", mod.validate_thomson_pair_witness, bad)
        bad = copy.deepcopy(th0)
        bad["lower_energy_sites"], bad["higher_energy_sites"] = bad["higher_energy_sites"], bad["lower_energy_sites"]
        _expect_err("THOMSON_E2_NOT_STRICT", mod.validate_thomson_pair_witness, bad)
        bad = copy.deepcopy(th0)
        bad["higher_energy_sites"] = copy.deepcopy(bad["lower_energy_sites"])
        _expect_err("THOMSON_E2_NOT_STRICT", mod.validate_thomson_pair_witness, bad)
        bad = copy.deepcopy(th0); bad["sqrt_bracket_bits"] = 4
        _expect_err("THOMSON_E1_NOT_STRICT", mod.validate_thomson_pair_witness, bad)
        eq_lo = {"e2": Fraction(1, 1), "e1_hi": Fraction(5, 1), "e1_lo": Fraction(4, 1)}
        eq_hi = {"e2": Fraction(2, 1), "e1_hi": Fraction(6, 1), "e1_lo": Fraction(5, 1)}
        with mock.patch.object(mod, "_eval_thomson_s2_config", side_effect=[eq_lo, eq_hi]):
            _expect_err("THOMSON_E1_NOT_STRICT", mod.validate_thomson_pair_witness, th0)

        # Sphere bridges
        br0 = copy.deepcopy(pack["sphere_bridge_witnesses"])
        rbr = mod.validate_sphere_bridges(br0)
        _req(len(rbr["stereographic_pairs"]) == 2, "stereo count")
        _req(len(rbr["cayley_s3_points"]) == 2, "cayley count")
        _req(not rbr["cayley_s3_points"][0]["equatorial_matches_stereographic"], "non-eq flag")
        _req(rbr["cayley_s3_points"][1]["equatorial_matches_stereographic"], "eq flag")

        _expect_err("INVALID_WITNESS", mod.validate_sphere_bridges, "x")
        bad = copy.deepcopy(br0); bad["stereographic_s2q_pairs"] = tuple(bad["stereographic_s2q_pairs"])
        _expect_err("INVALID_WITNESS", mod.validate_sphere_bridges, bad)
        bad = copy.deepcopy(br0); bad["stereographic_s2q_pairs"] = bad["stereographic_s2q_pairs"][:1]
        _expect_err("INVALID_WITNESS", mod.validate_sphere_bridges, bad)
        bad = copy.deepcopy(br0); bad["cayley_s3q_points"] = tuple(bad["cayley_s3q_points"])
        _expect_err("INVALID_WITNESS", mod.validate_sphere_bridges, bad)
        bad = copy.deepcopy(br0); bad["cayley_s3q_points"] = bad["cayley_s3q_points"][:1]
        _expect_err("INVALID_WITNESS", mod.validate_sphere_bridges, bad)
        bad = copy.deepcopy(br0); bad["stereographic_s2q_pairs"][0] = "not_dict"
        _expect_err("INVALID_WITNESS", mod.validate_sphere_bridges, bad)
        bad = copy.deepcopy(br0); bad["cayley_s3q_points"][0] = "not_dict"
        _expect_err("INVALID_WITNESS", mod.validate_sphere_bridges, bad)

        # validate_witness_pack
        _expect_err("INVALID_PACK", mod.validate_witness_pack, [])
        bad_p = copy.deepcopy(pack); bad_p["schema"] = "wrong"
        _expect_err("SCHEMA_MISMATCH", mod.validate_witness_pack, bad_p)
        bad_p = copy.deepcopy(pack); bad_p["upstream_issue"] = "wrong"
        _expect_err("ISSUE_MISMATCH", mod.validate_witness_pack, bad_p)
        bad_p = copy.deepcopy(pack); bad_p["equilibrium_witnesses"] = tuple(bad_p["equilibrium_witnesses"])
        _expect_err("INVALID_PACK", mod.validate_witness_pack, bad_p)
        bad_p = copy.deepcopy(pack); bad_p["equilibrium_witnesses"] = bad_p["equilibrium_witnesses"][:2]
        _expect_err("INVALID_PACK", mod.validate_witness_pack, bad_p)
        bad_p = copy.deepcopy(pack); bad_p["box_collection_witnesses"] = tuple(bad_p["box_collection_witnesses"])
        _expect_err("INVALID_PACK", mod.validate_witness_pack, bad_p)
        bad_p = copy.deepcopy(pack); bad_p["box_collection_witnesses"] = bad_p["box_collection_witnesses"][:1]
        _expect_err("INVALID_PACK", mod.validate_witness_pack, bad_p)
        bad_p = copy.deepcopy(pack); bad_p["thomson_pair_witnesses"] = tuple(bad_p["thomson_pair_witnesses"])
        _expect_err("INVALID_PACK", mod.validate_witness_pack, bad_p)
        bad_p = copy.deepcopy(pack); bad_p["thomson_pair_witnesses"] = bad_p["thomson_pair_witnesses"][:1]
        _expect_err("INVALID_PACK", mod.validate_witness_pack, bad_p)
        bad_p = copy.deepcopy(pack); bad_p["equilibrium_witnesses"][0]["name"] = "renamed"
        _expect_err("WITNESS_SET_MISMATCH", mod.validate_witness_pack, bad_p)
        bad_p = copy.deepcopy(pack); bad_p["box_collection_witnesses"][0]["name"] = "renamed"
        _expect_err("WITNESS_SET_MISMATCH", mod.validate_witness_pack, bad_p)
        bad_p = copy.deepcopy(pack); bad_p["thomson_pair_witnesses"][0]["name"] = "renamed"
        _expect_err("WITNESS_SET_MISMATCH", mod.validate_witness_pack, bad_p)

    def test_05_static_lean4_gate_char_literal_safety_and_escapes(self) -> None:
        """Verify char-literal-safe '\"' stripping, comment/string handling, and trust-escape refusals."""
        validated = mod.validate_witness_pack(_load_real_pack())
        clean_lib = mod.compile_lean4_source(validated)
        clean_aud = mod.compile_axiom_audit_source()

        # Benign nested comments, line comments, char literals, strings with escapes stay clean
        annotated = (
            "/- outer /- sorry bv_decide -/ still comment -/\n"
            "-- native_decide\n"
            'def s1 : String := "sorry \\" \\n\\\nline2"\n'
            "def c1 : Char := '\"'\n"
            "def c2 : Char := '\\''\n"
            "def c3 : Char := '\\x22'\n"
            "def c4 : Char := '\\u0022'\n"
            "def id' : Nat := 1\n"
            "def id_2 : Nat := 2\n"
            + clean_lib
        )
        mod.static_gate_lean4_source(annotated, mode="library")
        mod.static_gate_lean4_source(clean_aud, mode="audit")

        # Direct boundary tests on _match_lean_char_literal_end & strip_lean_comments_and_strings
        _req(mod._match_lean_char_literal_end("'a'", 0, 3) == 3, "simple char")
        _req(mod._match_lean_char_literal_end("'\\n'", 0, 4) == 4, "escaped char")
        _req(mod._match_lean_char_literal_end("'\\x22'", 0, 6) == 6, "hex char")
        _req(mod._match_lean_char_literal_end("'\\u0022'", 0, 8) == 8, "uni char")
        _req(mod._match_lean_char_literal_end("'bad", 0, 4) is None, "unclosed char")
        _req(mod._match_lean_char_literal_end("'\\\n'", 0, 4) is None, "newline esc char")
        _req(mod._match_lean_char_literal_end("'\\xzz'", 0, 6) is None, "bad hex char")
        _req(mod._match_lean_char_literal_end("'\\uzzzz'", 0, 8) is None, "bad uni char")
        # Single-quote identifier prime `a'b'` and `_'b'` must stay intact (not stripped as char literals!)
        _req(mod.strip_lean_comments_and_strings("a'b'") == "a'b'", "alnum ident prime")
        _req(mod.strip_lean_comments_and_strings("_'b'") == "_'b'", "underscore ident prime")
        _req("sorry" in mod.strip_lean_comments_and_strings("def x' := 1\ntheorem t : False := by sorry\n"), "ident prime")
        # Critical char-literal '\"' & 3-char doc-comment `/--/` / `/-!/` attacks MUST NOT desync and hide `sorry`!
        char_attack = clean_lib + "\ndef q : Char := '\"'\ntheorem evil : False := by sorry\n-- \"\n"
        _expect_err("LEAN4_TRUST_ESCAPE", mod.static_gate_lean4_source, char_attack, mode="library")
        doc_attack = clean_lib + '\n/--/ " -/\ntheorem evil : False := by sorry\n/- " -/\n'
        _expect_err("LEAN4_TRUST_ESCAPE", mod.static_gate_lean4_source, doc_attack, mode="library")
        mod_doc_attack = clean_lib + '\n/-!/ " -/\ntheorem evil : False := by sorry\n/- " -/\n'
        _expect_err("LEAN4_TRUST_ESCAPE", mod.static_gate_lean4_source, mod_doc_attack, mode="library")
        _expect_err("LEAN4_UNTERMINATED_COMMENT", mod.strip_lean_comments_and_strings, "/-")

        _expect_err("LEAN4_BAD_CHARS", mod.static_gate_lean4_source, clean_lib + "\r", mode="library")
        _expect_err("LEAN4_BAD_CHARS", mod.static_gate_lean4_source, clean_lib + "\x00", mode="library")
        _expect_err("LEAN4_UNTERMINATED_COMMENT", mod.static_gate_lean4_source, "/- unclosed\n" + clean_lib, mode="library")
        _expect_err("LEAN4_UNTERMINATED_STRING", mod.static_gate_lean4_source, 'def s := "unclosed\n' + clean_lib, mode="library")
        _expect_err("LEAN4_UNTERMINATED_STRING", mod.static_gate_lean4_source, clean_lib + '\ndef s := "abc\\', mode="library")

        escapes = [
            "def «x» := 1\n",
            "theorem bad : False := by sorry\n",
            "theorem bad : False := by admit\n",
            "axiom bad : False\n",
            "theorem bad : True := by native_decide\n",
            "theorem bad : True := by bv_decide\n",
            "theorem bad : True := by bv_check\n",
            "opaque bad : Nat\n",
            "partial def bad : Nat := bad\n",
            "unsafe def bad : Nat := 0\n",
            "def bad := Classical.choice\n",
            "macro \"x\" : term => `(1)\n",
            "set_option debug.skipKernelTC true\n",
            "def bad := dbg_trace 1; 0\n",
            "def bad := panic! 0\n",
            "@[init] def bad := 0\n",
            "@[export foo] def bad := 0\n",
            "#eval 1\n",
            "#print Nat\n",
            "theorem bad : True := by decide +native\n",
            "theorem bad : True := by decide (config := { native := true })\n",
            "def bad := DecideConfig\n",
            "def bad := warn.sorry\n",
            "attribute [simp] t01_coulomb_single_hessian_trace_zero\n",
            "instance : Inhabited Vec3Z := ⟨vecZero⟩\n",
        ]
        for payload in escapes:
            _expect_err("LEAN4_TRUST_ESCAPE", mod.static_gate_lean4_source, clean_lib + "\n" + payload, mode="library")

        # Forbidden imports in library and audit modes
        _expect_err("LEAN4_FORBIDDEN_IMPORT", mod.static_gate_lean4_source, "import Init\n" + clean_lib, mode="library")
        _expect_err("LEAN4_FORBIDDEN_IMPORT", mod.static_gate_lean4_source, "import Std\n" + clean_aud, mode="audit")

        # Missing theorem / missing #print axioms
        no_thm = clean_lib.replace("theorem t01_coulomb_single_hessian_trace_zero", "def renamed_t01", 1)
        _expect_err("LEAN4_MISSING_THEOREM", mod.static_gate_lean4_source, no_thm, mode="library")
        no_ax = clean_aud.replace("#print axioms PointChargesExactCert.t01_coulomb_single_hessian_trace_zero", "")
        _expect_err("LEAN4_MISSING_AXIOMS", mod.static_gate_lean4_source, no_ax, mode="audit")

    def test_06_live_lean_and_leanchecker_verification_mocked_and_refusals(self) -> None:
        """Cover 100% of resolve_lean_binaries and verify_lean4_live_and_fresh paths (works under env -i)."""
        with tempfile.TemporaryDirectory(prefix="pc_mock_") as td:
            root = Path(td)
            prefix = root.joinpath("lean_prefix")
            bin_dir = prefix.joinpath("bin")
            bin_dir.mkdir(parents=True)
            fake_lean = bin_dir.joinpath("lean")
            fake_checker = bin_dir.joinpath("leanchecker")
            fake_lean.write_text("#!/bin/sh\nexit 0\n", encoding="ascii")
            fake_checker.write_text("#!/bin/sh\nexit 0\n", encoding="ascii")
            os.chmod(fake_lean, 0o755)
            os.chmod(fake_checker, 0o755)

            pr_ok = subprocess.CompletedProcess([], 0, stdout=str(prefix) + "\n", stderr="")
            with mock.patch.dict(os.environ, {"POINT_CHARGES_LEAN_BIN": str(fake_lean)}, clear=True):
                with mock.patch.object(subprocess, "run", return_value=pr_ok):
                    lb, cb, pf = mod.resolve_lean_binaries()
                    _req(lb == str(fake_lean), "lean bin ok")
                    _req(cb == str(fake_checker), "checker bin ok")
                    _req(pf == str(prefix), "prefix ok")
                with mock.patch.object(subprocess, "run", side_effect=OSError("noexec")):
                    _expect_err("LEAN_EXEC_ERROR", mod.resolve_lean_binaries)
                with mock.patch.object(subprocess, "run", return_value=subprocess.CompletedProcess([], 1, str(prefix) + "\n", "err")):
                    _expect_err("LEAN_EXEC_ERROR", mod.resolve_lean_binaries)
                with mock.patch.object(subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "   \n", "")):
                    _expect_err("LEAN_EXEC_ERROR", mod.resolve_lean_binaries)

            # Non-executable env binary and executable directory (is_file False, X_OK True)
            noexec_lean = root.joinpath("noexec_lean")
            noexec_lean.write_text("#!/bin/sh\n", encoding="ascii")
            os.chmod(noexec_lean, 0o644)
            with mock.patch.dict(os.environ, {"POINT_CHARGES_LEAN_BIN": str(noexec_lean)}, clear=True):
                _expect_err("LEAN_UNAVAILABLE", mod.resolve_lean_binaries)
            with mock.patch.dict(os.environ, {"POINT_CHARGES_LEAN_BIN": str(bin_dir)}, clear=True):
                _expect_err("LEAN_UNAVAILABLE", mod.resolve_lean_binaries)

            # Empty env (no HOME, no which) -> LEAN_UNAVAILABLE (kills cand is not None and cand.is_file())
            with mock.patch.dict(os.environ, {}, clear=True):
                with mock.patch.object(mod.shutil, "which", return_value=None):
                    _expect_err("LEAN_UNAVAILABLE", mod.resolve_lean_binaries)

            # HOME/.elan fallback (absent, non-executable, then executable)
            with mock.patch.dict(os.environ, {"HOME": str(root)}, clear=True):
                with mock.patch.object(mod.shutil, "which", return_value=None):
                    _expect_err("LEAN_UNAVAILABLE", mod.resolve_lean_binaries)
                    elan_lean = root.joinpath(".elan", "bin", "lean")
                    elan_lean.parent.mkdir(parents=True)
                    elan_lean.write_text("#!/bin/sh\n", encoding="ascii")
                    os.chmod(elan_lean, 0o644)
                    _expect_err("LEAN_UNAVAILABLE", mod.resolve_lean_binaries)
                    os.chmod(elan_lean, 0o755)
                    with mock.patch.object(subprocess, "run", return_value=pr_ok):
                        lb2, _, _ = mod.resolve_lean_binaries()
                        _req(lb2 == str(fake_lean), "HOME/.elan fallback ok")

            # Missing / directory / non-executable real_lean and leanchecker
            wrapper_lean = root.joinpath("shim_lean")
            wrapper_lean.write_text("#!/bin/sh\nexit 0\n", encoding="ascii")
            os.chmod(wrapper_lean, 0o755)
            os.chmod(fake_lean, 0o644)
            with mock.patch.dict(os.environ, {"POINT_CHARGES_LEAN_BIN": str(wrapper_lean)}, clear=True):
                with mock.patch.object(subprocess, "run", return_value=pr_ok):
                    _expect_err("LEAN_UNAVAILABLE", mod.resolve_lean_binaries)
            os.unlink(fake_lean)
            fake_lean.mkdir(mode=0o755)
            with mock.patch.dict(os.environ, {"POINT_CHARGES_LEAN_BIN": str(wrapper_lean)}, clear=True):
                with mock.patch.object(subprocess, "run", return_value=pr_ok):
                    _expect_err("LEAN_UNAVAILABLE", mod.resolve_lean_binaries)
            fake_lean.rmdir()
            fake_lean.write_text("#!/bin/sh\nexit 0\n", encoding="ascii")
            os.chmod(fake_lean, 0o755)

            os.chmod(fake_checker, 0o644)
            with mock.patch.dict(os.environ, {"POINT_CHARGES_LEAN_BIN": str(fake_lean)}, clear=True):
                with mock.patch.object(subprocess, "run", return_value=pr_ok):
                    _expect_err("LEANCHECKER_UNAVAILABLE", mod.resolve_lean_binaries)
            os.unlink(fake_checker)
            fake_checker.mkdir(mode=0o755)
            with mock.patch.dict(os.environ, {"POINT_CHARGES_LEAN_BIN": str(fake_lean)}, clear=True):
                with mock.patch.object(subprocess, "run", return_value=pr_ok):
                    _expect_err("LEANCHECKER_UNAVAILABLE", mod.resolve_lean_binaries)
            fake_checker.rmdir()
            fake_checker.write_text("#!/bin/sh\nexit 0\n", encoding="ascii")
            os.chmod(fake_checker, 0o755)

            # verify_lean4_live_and_fresh
            validated = mod.validate_witness_pack(_load_real_pack())
            proj_files = {
                mod.LEAN_FILENAME: mod.compile_lean4_source(validated).encode("utf-8"),
                mod.AXIOM_AUDIT_FILENAME: mod.compile_axiom_audit_source().encode("utf-8"),
                mod.LAKEFILE_FILENAME: mod.EXPECTED_LAKEFILE_TEXT.encode("ascii"),
                mod.MANIFEST_FILENAME: mod.EXPECTED_MANIFEST_TEXT.encode("ascii"),
                mod.TOOLCHAIN_FILENAME: mod.EXPECTED_TOOLCHAIN_TEXT.encode("ascii"),
            }
            bad_utf = dict(proj_files); bad_utf[mod.LEAN_FILENAME] = b"\xff\xfe"
            _expect_err("INVALID_UTF8", mod.verify_lean4_live_and_fresh, bad_utf)

            ver_ok = subprocess.CompletedProcess([], 0, stdout="Lean (version 4.36.0, x86_64-linux)\n", stderr="")
            comp_ok = subprocess.CompletedProcess([], 0, stdout="", stderr="")
            aud_ok = subprocess.CompletedProcess([], 0, stdout=_clean_axiom_stdout(), stderr="")
            chk_ok = subprocess.CompletedProcess([], 0, stdout="", stderr="")

            with mock.patch.object(mod, "resolve_lean_binaries", return_value=(str(fake_lean), str(fake_checker), str(prefix))):
                # Verify cold tempdir receives ALL 5 files including lake-manifest.json
                seen_files: list[set[str]] = []
                def _spy_run(cmd, *a, **kw):
                    if len(cmd) >= 3 and cmd[1] == "-R":
                        seen_files.append({p.name for p in Path(cmd[2]).iterdir()})
                    return [ver_ok, comp_ok, aud_ok, chk_ok][len(seen_files) if len(cmd) < 2 else (1 if "-o" in cmd else 2)]

                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, comp_ok, aud_ok, chk_ok]):
                    vr = mod.verify_lean4_live_and_fresh(proj_files)
                    _req(vr["theorems_verified"] == 22, f"bad thm count: {vr}")
                    _req(vr["leanchecker_fresh_ec"] == 0, f"bad checker ec: {vr}")

                # Confirm cold tmpdir includes *.json (lake-manifest.json)
                def _capture_run(cmd, *a, **kw):
                    if len(cmd) >= 3 and cmd[1] == "-R" and "-o" in cmd:
                        seen_files.append({p.name for p in Path(cmd[2]).iterdir()})
                        return comp_ok
                    if "--version" in cmd:
                        return ver_ok
                    if "--fresh" in cmd:
                        return chk_ok
                    return aud_ok
                with mock.patch.object(subprocess, "run", side_effect=_capture_run):
                    mod.verify_lean4_live_and_fresh(proj_files)
                _req(seen_files and mod.MANIFEST_FILENAME in seen_files[0], f"cold tmpdir missing lake-manifest.json: {seen_files}")

                # Every failure path in verify_lean4_live_and_fresh
                with mock.patch.object(subprocess, "run", side_effect=OSError("ver")):
                    _expect_err("LEAN_EXEC_ERROR", mod.verify_lean4_live_and_fresh, proj_files)
                with mock.patch.object(subprocess, "run", return_value=subprocess.CompletedProcess([], 1, "Lean 4.36.0", "")):
                    _expect_err("LEAN_VERSION_MISMATCH", mod.verify_lean4_live_and_fresh, proj_files)
                with mock.patch.object(subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "Lean 4.18.0", "")):
                    _expect_err("LEAN_VERSION_MISMATCH", mod.verify_lean4_live_and_fresh, proj_files)
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, OSError("comp")]):
                    _expect_err("LEAN_EXEC_ERROR", mod.verify_lean4_live_and_fresh, proj_files)
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, subprocess.CompletedProcess([], 1, "", "err")]):
                    _expect_err("LEAN_COMPILE_FAIL", mod.verify_lean4_live_and_fresh, proj_files)
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, subprocess.CompletedProcess([], 0, "", "stderr!")]):
                    _expect_err("LEAN_STDERR_NONEMPTY", mod.verify_lean4_live_and_fresh, proj_files)
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, subprocess.CompletedProcess([], 0, "warning: unused\n", "")]):
                    _expect_err("LEAN_WARNING", mod.verify_lean4_live_and_fresh, proj_files)
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, comp_ok, OSError("aud")]):
                    _expect_err("LEAN_EXEC_ERROR", mod.verify_lean4_live_and_fresh, proj_files)
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, comp_ok, subprocess.CompletedProcess([], 1, "", "err")]):
                    _expect_err("LEAN_AUDIT_FAIL", mod.verify_lean4_live_and_fresh, proj_files)
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, comp_ok, subprocess.CompletedProcess([], 0, _clean_axiom_stdout(), "stderr!")]):
                    _expect_err("LEAN_STDERR_NONEMPTY", mod.verify_lean4_live_and_fresh, proj_files)
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, comp_ok, subprocess.CompletedProcess([], 0, "warning: x\n", "")]):
                    _expect_err("LEAN_WARNING", mod.verify_lean4_live_and_fresh, proj_files)
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, comp_ok, subprocess.CompletedProcess([], 0, "  leading space\n", "")]):
                    _expect_err("LEAN_UNEXPECTED_OUT", mod.verify_lean4_live_and_fresh, proj_files)

                dup_none = f"'{mod.EXPECTED_THEOREM_NAMES[14]}' does not depend on any axioms\n" * 2
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, comp_ok, subprocess.CompletedProcess([], 0, dup_none, "")]):
                    _expect_err("LEAN_DUP_THM", mod.verify_lean4_live_and_fresh, proj_files)
                dup_dep = f"'{mod.EXPECTED_THEOREM_NAMES[0]}' depends on axioms: [propext]\n" * 2
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, comp_ok, subprocess.CompletedProcess([], 0, dup_dep, "")]):
                    _expect_err("LEAN_DUP_THM", mod.verify_lean4_live_and_fresh, proj_files)
                bad_ax = f"'{mod.EXPECTED_THEOREM_NAMES[0]}' depends on axioms: [propext, Classical.choice]\n"
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, comp_ok, subprocess.CompletedProcess([], 0, bad_ax, "")]):
                    _expect_err("LEAN_FORBIDDEN_AXIOM", mod.verify_lean4_live_and_fresh, proj_files)
                short_out = f"'{mod.EXPECTED_THEOREM_NAMES[0]}' depends on axioms: [propext]\n"
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, comp_ok, subprocess.CompletedProcess([], 0, short_out, "")]):
                    _expect_err("LEAN_THM_SET_MISMATCH", mod.verify_lean4_live_and_fresh, proj_files)
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, comp_ok, aud_ok, OSError("chk")]):
                    _expect_err("LEANCHECKER_EXEC_ERROR", mod.verify_lean4_live_and_fresh, proj_files)
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, comp_ok, aud_ok, subprocess.CompletedProcess([], 1, "", "")]):
                    _expect_err("LEANCHECKER_FAIL", mod.verify_lean4_live_and_fresh, proj_files)
                with mock.patch.object(subprocess, "run", side_effect=[ver_ok, comp_ok, aud_ok, subprocess.CompletedProcess([], 0, "", "stderr!")]):
                    _expect_err("LEANCHECKER_FAIL", mod.verify_lean4_live_and_fresh, proj_files)

    def test_07_end_to_end_compile_check_tamper_all_pinned_files_and_cli(self) -> None:
        """End-to-end compile, self-bound check, 1-byte tamper battery on ALL 6 pinned files, and CLI."""
        with tempfile.TemporaryDirectory(prefix="pc_e2e_") as td:
            root = Path(td)
            wp_bytes = BUNDLE_DIR.joinpath(mod.WITNESS_PACK_FILENAME).read_bytes()
            root.joinpath(mod.WITNESS_PACK_FILENAME).write_bytes(wp_bytes)

            comp = mod.compile_and_verify_bundle(root, live_lean=False)
            _req(comp["compiled"], f"compile failed: {comp}")
            _req(comp["theorems_compiled"] == 22, "expected 22 theorems")
            _req(not comp["live_lean_checked"], "offline compile must mark live_lean_checked=False")
            _req(comp["lean_live"] is None, "offline compile lean_live must be None")
            cert_path = root.joinpath(mod.DEFAULT_CERT_FILENAME)
            _req((os.lstat(cert_path).st_mode & 0o777) == 0o644, "cert must be 0644")

            chk = mod.check_certificate_bundle(cert_path, live_lean=False)
            _req(chk["valid"], f"check failed: {chk}")
            _req(not chk["live_lean_checked"], "offline check live_lean_checked")

            # Tamper envelope keys, schema, payload_sha256, certificate_sha256
            orig_cert = json.loads(cert_path.read_text(encoding="ascii"))
            t_extra = copy.deepcopy(orig_cert); t_extra["extra"] = 1
            root.joinpath("t_extra.json").write_bytes(mod.canonical_json_bytes(t_extra))
            _expect_err("ENVELOPE_KEYS", mod.check_certificate_bundle, root.joinpath("t_extra.json"), bundle_dir=root, live_lean=False)

            t_sch = copy.deepcopy(orig_cert); t_sch["schema_version"] = "bad"
            root.joinpath("t_sch.json").write_bytes(mod.canonical_json_bytes(t_sch))
            _expect_err("SCHEMA_MISMATCH", mod.check_certificate_bundle, root.joinpath("t_sch.json"), bundle_dir=root, live_lean=False)

            t_psha = copy.deepcopy(orig_cert); t_psha["payload_sha256"] = "0" * 64
            root.joinpath("t_psha.json").write_bytes(mod.canonical_json_bytes(t_psha))
            _expect_err("PAYLOAD_SHA_MISMATCH", mod.check_certificate_bundle, root.joinpath("t_psha.json"), bundle_dir=root, live_lean=False)

            t_csha = copy.deepcopy(orig_cert); t_csha["certificate_sha256"] = "0" * 64
            root.joinpath("t_csha.json").write_bytes(mod.canonical_json_bytes(t_csha))
            _expect_err("CERT_SHA_MISMATCH", mod.check_certificate_bundle, root.joinpath("t_csha.json"), bundle_dir=root, live_lean=False)

            # Re-signed semantic forgery AND pretty-printed non-canonical formatting both fail FULL_REPLAY_MISMATCH
            t_forge = copy.deepcopy(orig_cert)
            t_forge["payload"]["oracle_class"] = "forged"
            t_forge["payload_sha256"] = mod.sha256_hex(mod.canonical_json_bytes(t_forge["payload"]))
            eb = {k: t_forge[k] for k in ("schema_version", "payload_sha256", "payload")}
            t_forge["certificate_sha256"] = mod.sha256_hex(mod.canonical_json_bytes(eb))
            root.joinpath("t_forge.json").write_bytes(mod.canonical_json_bytes(t_forge))
            _expect_err("FULL_REPLAY_MISMATCH", mod.check_certificate_bundle, root.joinpath("t_forge.json"), bundle_dir=root, live_lean=False)

            # Duplicate-key JSON injection on certificate fails DUPLICATE_JSON_KEY
            dup_raw = b'{"schema_version":"bad",' + cert_path.read_bytes()[1:]
            root.joinpath("t_dup.json").write_bytes(dup_raw)
            _expect_err("DUPLICATE_JSON_KEY", mod.check_certificate_bundle, root.joinpath("t_dup.json"), bundle_dir=root, live_lean=False)

            # 1-byte tamper on EACH of the 5 Lean/Lake/manifest/toolchain files fails closed
            for fname, exp_code in (
                (mod.LEAN_FILENAME, "LEAN_DRIFT"),
                (mod.AXIOM_AUDIT_FILENAME, "AUDIT_DRIFT"),
                (mod.LAKEFILE_FILENAME, "LAKEFILE_DRIFT"),
                (mod.MANIFEST_FILENAME, "MANIFEST_DRIFT"),
                (mod.TOOLCHAIN_FILENAME, "TOOLCHAIN_DRIFT"),
            ):
                fp = root.joinpath(fname)
                orig_b = fp.read_bytes()
                fp.write_bytes(orig_b + b"\n")
                _expect_err(exp_code, mod.check_certificate_bundle, cert_path, live_lean=False)
                fp.write_bytes(orig_b)

            # 1-byte whitespace tamper on witness_pack.json fails FULL_REPLAY_MISMATCH (raw-bytes SHA-256 bound!)
            wp_p = root.joinpath(mod.WITNESS_PACK_FILENAME)
            wp_p.write_bytes(wp_bytes + b" ")
            _expect_err("FULL_REPLAY_MISMATCH", mod.check_certificate_bundle, cert_path, live_lean=False)
            wp_p.write_bytes(wp_bytes)

            # Real committed on-disk bundle CLI checks (--standalone, --check --offline, --json, error exit)
            rc_s, out_s, err_s = _run_cli(["--standalone"])
            _req(rc_s == 0, f"--standalone rc={rc_s} err={err_s}")
            _req(err_s == "", f"--standalone stderr={err_s}")
            cert_disk = BUNDLE_DIR.joinpath(mod.DEFAULT_CERT_FILENAME).read_bytes()
            _req(out_s.encode("ascii") == cert_disk, "--standalone output disagrees with committed certificate")

            rc_c1, out_c1, _ = _run_cli(["--check", "--offline"])
            _req(rc_c1 == 0, f"--check --offline rc={rc_c1}")
            _req(out_c1.startswith("OK "), f"--check text={out_c1}")
            rc_c2, out_c2, _ = _run_cli(["--check", str(BUNDLE_DIR.joinpath(mod.DEFAULT_CERT_FILENAME)), "--offline", "--json"])
            _req(rc_c2 == 0, f"--check --json rc={rc_c2}")
            _req(json.loads(out_c2)["valid"] is True, f"--check --json out={out_c2}")
            rc_c3, out_c3, _ = _run_cli(["--offline", "--check"])
            _req(rc_c3 == 0, f"--offline --check rc={rc_c3}")
            _req(out_c3.startswith("OK "), f"--offline --check out={out_c3}")

            with mock.patch.object(mod, "compile_and_verify_bundle", return_value=comp):
                rc_m1, out_m1, _ = _run_cli(["--compile", "--offline"])
                _req(rc_m1 == 0, f"--compile text rc={rc_m1}")
                _req(out_m1.startswith("COMPILED "), f"--compile text out={out_m1}")
                rc_m2, out_m2, _ = _run_cli(["--compile", "--offline", "--json"])
                _req(rc_m2 == 0, f"--compile json rc={rc_m2}")
                _req(json.loads(out_m2)["compiled"] is True, f"--compile json out={out_m2}")

            rc_bad, _, err_bad = _run_cli(["--bogus-flag"])
            _req(rc_bad == 2, f"unknown flag rc={rc_bad}")
            _req(json.loads(err_bad)["code"] == "UNKNOWN_FLAG", f"unknown flag err={err_bad}")


if __name__ == "__main__":
    unittest.main()
