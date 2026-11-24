"""R5 regression suite for Issue #21: verify_arb.py must fail closed under python3 -O.

Zero bare assert statements (100% R4 parity under -I -B, -I -O -B, and -I -B -W error).
Uses real python-flint whenever installed (e.g. in CI), and falls back cleanly inside
the test harness to a self-contained outward-rounded interval test double when running
on zero-flint hosts so every proof guard in verify_arb.py executes live under -O.
"""
from __future__ import annotations

import ast
import copy
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest

ROOT = Path(__file__).resolve().parents[1]
VERIFIER_PATH = ROOT / "verify_arb.py"
CERT_N3 = ROOT / "certificates" / "n3-triangle.json"
CERT_N4 = ROOT / "certificates" / "n4-triangle-plus-centre.json"
HAS_FLINT = importlib.util.find_spec("flint") is not None

# Exact minimal 2-charge certificate (1 non-degenerate zero at origin, Morse index 2).
CERT_N2_PAIR = {
    "name": "n2-pair",
    "sites": [[-1, 0, 0], [1, 0, 0]],
    "charges": ["1", "1"],
    "count": 1,
    "zeros": [{"box": [[-0.015625, -0.015625, -0.015625], [0.015625, 0.015625, 0.015625]], "index": 2}],
}

# Exact 3-charge collinear certificate (2 symmetric non-degenerate zeros at x = +/-1.02724, Morse index 2).
CERT_N3_LINE = {
    "name": "n3-line",
    "sites": [[-2, 0, 0], [0, 0, 0], [2, 0, 0]],
    "charges": ["1", "1", "1"],
    "count": 2,
    "zeros": [
        {"box": [[-1.03724, -0.01, -0.01], [-1.01724, 0.01, 0.01]], "index": 2},
        {"box": [[1.01724, -0.01, -0.01], [1.03724, 0.01, 0.01]], "index": 2},
    ],
}


def _check(cond: bool, msg: str) -> None:
    if not cond:
        raise RuntimeError(f"ISSUE21_R5_FAILED: {msg}")


# Self-contained outward-rounded interval test double for python-flint (test-only fallback).
_FLINT_SHIM_SOURCE = r'''
from fractions import Fraction
import math


class _Ctx:
    prec = 160


ctx = _Ctx()


class fmpq(Fraction):
    pass


def _dn(x):
    return math.nextafter(x, -math.inf)


def _up(x):
    return math.nextafter(x, math.inf)


def _bounds_from_exact(val):
    if isinstance(val, arb):
        return val._lo, val._hi
    if isinstance(val, Fraction):
        f = float(val)
        if not math.isfinite(f):
            return f, f
        qf = Fraction(f)
        lo = f if qf <= val else _dn(f)
        hi = f if qf >= val else _up(f)
        return lo, hi
    f = float(val)
    return f, f


class arb:
    __slots__ = ("_lo", "_hi")

    def __init__(self, x, rad=None):
        if rad is None:
            self._lo, self._hi = _bounds_from_exact(x)
        else:
            if isinstance(x, Fraction) and isinstance(rad, Fraction):
                lo, _ = _bounds_from_exact(x - rad)
                _, hi = _bounds_from_exact(x + rad)
                self._lo, self._hi = lo, hi
            else:
                mlo, mhi = _bounds_from_exact(x)
                rlo, rhi = _bounds_from_exact(rad)
                r = max(abs(rlo), abs(rhi))
                self._lo = _dn(mlo - r)
                self._hi = _up(mhi + r)

    @classmethod
    def _raw(cls, lo, hi):
        obj = cls.__new__(cls)
        obj._lo = lo
        obj._hi = hi
        return obj

    def lower(self):
        return arb._raw(self._lo, self._lo)

    def upper(self):
        return arb._raw(self._hi, self._hi)

    def mid(self):
        m = 0.5 * (self._lo + self._hi)
        return arb._raw(m, m)

    def union(self, other):
        o = other if isinstance(other, arb) else arb(other)
        return arb._raw(min(self._lo, o._lo), max(self._hi, o._hi))

    def contains(self, val):
        vlo, vhi = _bounds_from_exact(val)
        return self._lo <= vlo and vhi <= self._hi

    def sqrt(self):
        if self._lo < 0.0:
            return arb._raw(float("nan"), float("nan"))
        lo = 0.0 if self._lo == 0.0 else _dn(math.sqrt(self._lo))
        hi = _up(math.sqrt(self._hi))
        return arb._raw(lo, hi)

    def __float__(self):
        return 0.5 * (self._lo + self._hi)

    def __neg__(self):
        return arb._raw(-self._hi, -self._lo)

    def __abs__(self):
        if self._lo >= 0.0:
            return arb._raw(self._lo, self._hi)
        if self._hi <= 0.0:
            return arb._raw(-self._hi, -self._lo)
        return arb._raw(0.0, max(-self._lo, self._hi))

    def __add__(self, other):
        o = other if isinstance(other, arb) else arb(other)
        return arb._raw(_dn(self._lo + o._lo), _up(self._hi + o._hi))

    def __radd__(self, other):
        return arb(other).__add__(self)

    def __sub__(self, other):
        o = other if isinstance(other, arb) else arb(other)
        return arb._raw(_dn(self._lo - o._hi), _up(self._hi - o._lo))

    def __rsub__(self, other):
        return arb(other).__sub__(self)

    def __mul__(self, other):
        o = other if isinstance(other, arb) else arb(other)
        ps = (self._lo * o._lo, self._lo * o._hi, self._hi * o._lo, self._hi * o._hi)
        return arb._raw(_dn(min(ps)), _up(max(ps)))

    def __rmul__(self, other):
        return arb(other).__mul__(self)

    def __truediv__(self, other):
        o = other if isinstance(other, arb) else arb(other)
        if o._lo <= 0.0 <= o._hi:
            return arb._raw(float("-inf"), float("inf"))
        inv = arb._raw(_dn(1.0 / o._hi), _up(1.0 / o._lo))
        return self.__mul__(inv)

    def __rtruediv__(self, other):
        return arb(other).__truediv__(self)

    def __lt__(self, other):
        o = other if isinstance(other, arb) else arb(other)
        return self._hi < o._lo

    def __le__(self, other):
        o = other if isinstance(other, arb) else arb(other)
        return self._hi <= o._lo

    def __gt__(self, other):
        o = other if isinstance(other, arb) else arb(other)
        return self._lo > o._hi

    def __ge__(self, other):
        o = other if isinstance(other, arb) else arb(other)
        return self._lo >= o._hi
'''


def _load_verify_arb_module():
    if not HAS_FLINT and "flint" not in sys.modules:
        mod = types.ModuleType("flint")
        exec(_FLINT_SHIM_SOURCE, mod.__dict__)
        sys.modules["flint"] = mod
    root_s = str(ROOT)
    if root_s not in sys.path:
        sys.path.insert(0, root_s)
    import verify_arb
    return verify_arb


def _run_verifier_subprocess(cert_path: Path, flags: list[str]) -> subprocess.CompletedProcess[str]:
    if HAS_FLINT:
        cmd = [sys.executable, *flags, str(VERIFIER_PATH), str(cert_path)]
        return subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    with tempfile.TemporaryDirectory(prefix="pc21_shim_") as shim_dir:
        (Path(shim_dir) / "flint.py").write_text(_FLINT_SHIM_SOURCE, encoding="utf-8")
        runner = (
            "import runpy, sys\n"
            f"sys.path.insert(0, {shim_dir!r})\n"
            f"sys.path.insert(0, {str(ROOT)!r})\n"
            f"sys.argv = [{str(VERIFIER_PATH)!r}, {str(cert_path)!r}]\n"
            f"runpy.run_path({str(VERIFIER_PATH)!r}, run_name='__main__')\n"
        )
        cmd = [sys.executable, *flags, "-c", runner]
        return subprocess.run(cmd, capture_output=True, text=True, timeout=120)


def _write_json(path: Path, obj: object) -> None:
    path.write_text(json.dumps(obj), encoding="utf-8")


class TestVerifyArbOSafe(unittest.TestCase):
    def test_01_zero_bare_asserts_in_verifier_and_suite(self) -> None:
        for path in (VERIFIER_PATH, Path(__file__).resolve()):
            tree = ast.parse(path.read_text(encoding="utf-8"))
            asserts = [n.lineno for n in ast.walk(tree) if isinstance(n, ast.Assert)]
            _check(len(asserts) == 0, f"bare asserts found in {path.name}: {asserts}")

    def test_02_clean_and_corrupted_certificates_across_all_3_isolated_modes(self) -> None:
        base_n4 = json.loads(CERT_N4.read_text(encoding="utf-8"))

        corrupt_cases: list[tuple[str, dict]] = []
        # Full-pipeline corruptions on exact minimal 2-charge certificate (covers all proof steps + completeness)
        c = copy.deepcopy(CERT_N2_PAIR)
        c["zeros"].pop(0)
        c["count"] = 0
        corrupt_cases.append(("n2 one zero removed", c))

        c = copy.deepcopy(CERT_N2_PAIR)
        c["zeros"].append(copy.deepcopy(c["zeros"][0]))
        c["count"] += 1
        corrupt_cases.append(("n2 one zero duplicated", c))

        c = copy.deepcopy(CERT_N2_PAIR)
        c["zeros"][0]["index"] = 1
        corrupt_cases.append(("n2 wrong Morse index", c))

        c = copy.deepcopy(CERT_N2_PAIR)
        c["charges"][0] = "2"
        corrupt_cases.append(("n2 charge perturbed", c))

        c = copy.deepcopy(CERT_N2_PAIR)
        c["count"] += 1
        corrupt_cases.append(("n2 claimed count too high", c))

        c = copy.deepcopy(CERT_N2_PAIR)
        c["zeros"][0]["box"][0][0] = c["zeros"][0]["box"][1][0]
        corrupt_cases.append(("n2 zero box squeezed to a plane", c))

        # Exact Issue #21 four-charge corruptions on committed certificates/n4-triangle-plus-centre.json
        c4_dup = copy.deepcopy(base_n4)
        c4_dup["zeros"].append(copy.deepcopy(c4_dup["zeros"][0]))
        c4_dup["count"] += 1
        corrupt_cases.append(("n4 duplicated zero box (Issue #21)", c4_dup))

        c4_idx = copy.deepcopy(base_n4)
        c4_idx["zeros"][0]["index"] = 3 - c4_idx["zeros"][0]["index"]
        corrupt_cases.append(("n4 wrong Morse index (Issue #21)", c4_idx))

        c4_chg = copy.deepcopy(base_n4)
        c4_chg["charges"][3] = "1/32"
        corrupt_cases.append(("n4 centre charge changed to 1/32", c4_chg))

        c4_sqz = copy.deepcopy(base_n4)
        c4_sqz["zeros"][0]["box"][0][0] = c4_sqz["zeros"][0]["box"][1][0]
        corrupt_cases.append(("n4 zero box squeezed to a plane", c4_sqz))

        with tempfile.TemporaryDirectory(prefix="pc21_certs_") as td:
            root = Path(td)
            p_n2 = root / "clean_n2.json"
            _write_json(p_n2, CERT_N2_PAIR)
            p_n3_line = root / "clean_n3_line.json"
            _write_json(p_n3_line, CERT_N3_LINE)

            for flags in (["-I", "-B"], ["-I", "-O", "-B"], ["-I", "-B", "-W", "error"]):
                r_n2 = _run_verifier_subprocess(p_n2, flags)
                _check(
                    r_n2.returncode == 0 and "VERIFIED n2-pair: exactly 1 zeros" in r_n2.stdout,
                    f"clean n2-pair failed under {flags}: rc={r_n2.returncode} out={r_n2.stdout!r} err={r_n2.stderr!r}",
                )
                r_l3 = _run_verifier_subprocess(p_n3_line, flags)
                _check(
                    r_l3.returncode == 0 and "VERIFIED n3-line: exactly 2 zeros" in r_l3.stdout,
                    f"clean n3-line failed under {flags}: rc={r_l3.returncode} out={r_l3.stdout!r} err={r_l3.stderr!r}",
                )
                if HAS_FLINT:
                    r_ok = _run_verifier_subprocess(CERT_N3, flags)
                    _check(
                        r_ok.returncode == 0 and "VERIFIED n3-triangle: exactly 4 zeros" in r_ok.stdout,
                        f"clean n3-triangle failed under {flags}: rc={r_ok.returncode} err={r_ok.stderr!r}",
                    )
                for idx, (label, payload) in enumerate(corrupt_cases):
                    p = root / f"bad_{idx}.json"
                    _write_json(p, payload)
                    r_bad = _run_verifier_subprocess(p, flags)
                    _check(
                        r_bad.returncode != 0 and "VERIFIED" not in r_bad.stdout,
                        f"corrupt '{label}' accepted under {flags}: rc={r_bad.returncode} out={r_bad.stdout!r}",
                    )

    def test_03_in_process_guards_empty_intersection_and_zone_b_schema(self) -> None:
        va = _load_verify_arb_module()
        base = copy.deepcopy(CERT_N2_PAIR)

        def run_cert(cert_obj: dict) -> tuple[int, int, int]:
            with tempfile.TemporaryDirectory(prefix="pc21_ok_") as td:
                p = Path(td) / "ok.json"
                _write_json(p, cert_obj)
                return va.verify(str(p))

        # Positive controls: 1-zero cert + 2-zero certs in BOTH orders (kills L215 or->and & directional checks)
        _check(run_cert(CERT_N2_PAIR) == (1, 0, 1), "in-process n2-pair verify failed")
        _check(run_cert(CERT_N3_LINE) == (2, 0, 2), "in-process n3-line LR verify failed")
        rev_line = copy.deepcopy(CERT_N3_LINE)
        rev_line["zeros"].reverse()
        _check(run_cert(rev_line) == (2, 0, 2), "in-process n3-line RL verify failed")

        # Single-charge Config boundary (kills len(charges) >= 1 -> len(charges) > 1)
        cfg1 = va.Config([[0, 0, 0]], ["1"])
        _check(cfg1.n == 1, "single-charge Config failed")

        def expect_fail(
            cert_obj: object,
            label: str,
            *,
            raw_text: str | None = None,
            substr: str | None = None,
        ) -> None:
            with tempfile.TemporaryDirectory(prefix="pc21_zb_") as td:
                p = Path(td) / "c.json"
                if raw_text is not None:
                    p.write_text(raw_text, encoding="utf-8")
                else:
                    _write_json(p, cert_obj)
                try:
                    va.verify(str(p))
                except va.VerificationError as exc:
                    if substr is not None:
                        _check(substr in str(exc), f"{label}: expected {substr!r} in {str(exc)!r}")
                    return
                raise RuntimeError(f"Expected VerificationError for {label}, but verify() passed")

        # Missing file & raw JSON syntax / duplicate key / NaN constant / non-dict root
        try:
            va.verify("/nonexistent_pc21_cert_path.json")
            raise RuntimeError("missing file accepted")
        except va.VerificationError:
            pass
        expect_fail(None, "syntax error", raw_text="{bad")
        expect_fail(None, "duplicate key", raw_text='{"count": 1, "count": 1}', substr="duplicate JSON key")
        expect_fail(None, "NaN constant", raw_text='{"name": "x", "count": NaN}', substr="non-finite JSON constant")
        expect_fail([1, 2, 3], "list root", substr="JSON object")

        # Missing required keys
        for key in ("name", "sites", "charges", "zeros", "count"):
            c = copy.deepcopy(base)
            del c[key]
            expect_fail(c, f"missing key {key}", substr=f"missing required key: {key}")

        # Bad scalar/container types, exact boundaries, & bool-as-int
        c = copy.deepcopy(base); c["name"] = 123
        expect_fail(c, "non-str name", substr="name must be str")
        c = copy.deepcopy(base); c["count"] = True
        expect_fail(c, "bool count", substr="non-negative int")
        c = copy.deepcopy(base); c["count"] = -1
        expect_fail(c, "negative count", substr="non-negative int")
        c = copy.deepcopy(base); c["zeros"] = {}
        expect_fail(c, "non-list zeros", substr="zeros must be a list")
        c = copy.deepcopy(base); c["zeros"][0] = "bad"
        expect_fail(c, "non-dict zero entry", substr="dict with box and index")
        c = copy.deepcopy(base); del c["zeros"][0]["box"]
        expect_fail(c, "zero missing box", substr="dict with box and index")
        c = copy.deepcopy(base); del c["zeros"][0]["index"]
        expect_fail(c, "zero missing index", substr="dict with box and index")
        c = copy.deepcopy(base); c["zeros"][0]["index"] = True
        expect_fail(c, "bool index == 1", substr="exact int")
        c = copy.deepcopy(base); c["zeros"][0]["box"] = [[0, 0, 0]]
        expect_fail(c, "box len == 1", substr="[lo, hi]")
        c = copy.deepcopy(base); c["zeros"][0]["box"] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        expect_fail(c, "box len == 3", substr="[lo, hi]")
        c = copy.deepcopy(base); c["zeros"][0]["box"] = "ab"
        expect_fail(c, "non-list box len 2", substr="[lo, hi]")
        c = copy.deepcopy(base); c["zeros"][0]["box"][0] = [0, 0]
        expect_fail(c, "box bound len == 2", substr="3-element list")
        c = copy.deepcopy(base); c["zeros"][0]["box"][0] = [0, 0, 0, 0]
        expect_fail(c, "box bound len == 4", substr="3-element list")
        c = copy.deepcopy(base); c["zeros"][0]["box"][0] = "bad"
        expect_fail(c, "box bound non-list len 3", substr="3-element list")
        c = copy.deepcopy(base); c["zeros"][0]["box"][0][0] = True
        expect_fail(c, "bool in box", substr="bool forbidden")
        c = copy.deepcopy(base); c["zeros"][0]["box"][0][0] = "not-a-number"
        expect_fail(c, "invalid numeric box bound", substr="invalid numeric")
        c = copy.deepcopy(base); c["zeros"][0]["box"][0][0], c["zeros"][0]["box"][1][0] = 2.0, 1.0
        expect_fail(c, "inverted box bound", substr="inverted zero box bound")

        # Config validation: sites/charges types, exact length boundaries, non-empty, bool, non-positive charges
        c = copy.deepcopy(base); c["sites"] = "bad"
        expect_fail(c, "non-list sites", substr="must be lists")
        c = copy.deepcopy(base); c["charges"] = "bad"
        expect_fail(c, "non-list charges", substr="must be lists")
        c = copy.deepcopy(base); c["charges"] = ["1"]
        expect_fail(c, "sites > charges", substr="same length")
        c = copy.deepcopy(base); c["charges"] = ["1", "1", "1"]
        expect_fail(c, "sites < charges", substr="same length")
        c = copy.deepcopy(base); c["sites"] = []; c["charges"] = []
        expect_fail(c, "empty sites and charges", substr="non-empty")
        c = copy.deepcopy(base); c["sites"][0] = [1, 0]
        expect_fail(c, "site len == 2", substr="3-element list")
        c = copy.deepcopy(base); c["sites"][0] = [1, 0, 0, 0]
        expect_fail(c, "site len == 4", substr="3-element list")
        c = copy.deepcopy(base); c["sites"][0] = "bad"
        expect_fail(c, "site non-list len 3", substr="3-element list")
        c = copy.deepcopy(base); c["sites"][0][0] = True
        expect_fail(c, "bool site coord", substr="bool forbidden")
        c = copy.deepcopy(base); c["charges"][0] = True
        expect_fail(c, "bool charge", substr="bool forbidden")
        c = copy.deepcopy(base); c["charges"][0] = "0"
        expect_fail(c, "zero charge", substr="strictly positive")
        c = copy.deepcopy(base); c["charges"][0] = "-1"
        expect_fail(c, "negative charge", substr="strictly positive")

        # Direct in-process hits for the 5 replaced proof guards & exact boundary comparisons
        # 1. Krawczyk guard (lo == hi passes _load_cert_json <= check, then fails Krawczyk -> kills L315 <= -> < too!)
        c = copy.deepcopy(base); c["zeros"][0]["box"][0][0] = c["zeros"][0]["box"][1][0]
        expect_fail(c, "krawczyk failure on squeezed box", substr="Krawczyk failed")
        # 2. Morse index mismatch in BOTH directions (kills L209 == -> <= and == -> >=)
        c = copy.deepcopy(base); c["zeros"][0]["index"] = 1
        expect_fail(c, "index mismatch low", substr="index mismatch")
        c = copy.deepcopy(base); c["zeros"][0]["index"] = 3
        expect_fail(c, "index mismatch high", substr="index mismatch")
        # 3. Overlapping enclosures
        c = copy.deepcopy(base); c["zeros"].append(copy.deepcopy(c["zeros"][0])); c["count"] += 1
        expect_fail(c, "overlapping enclosures", substr="overlap")
        # 4. Count mismatch in BOTH directions (count=0 also kills L305 >= 0 -> > 0!)
        c = copy.deepcopy(base); c["count"] = 2
        expect_fail(c, "count mismatch high", substr="certificate claims 2 zeros")
        c = copy.deepcopy(base); c["count"] = 0
        expect_fail(c, "count mismatch zero", substr="certificate claims 0 zeros")

        # 5a. Empty-intersection guard (disjoint K on shrinking step -> kills L204 require)
        orig_krawczyk = va.krawczyk
        orig_morse = va.morse_index
        try:
            def disjoint_shrink_krawczyk(cfg, lo, hi):
                return True, [va.arb(100), va.arb(100), va.arb(100)]

            va.krawczyk = disjoint_shrink_krawczyk
            expect_fail(base, "empty intersection on shrinking step", substr="empty intersection")

            # 5b. Boundary-touching shrinking step (nlo[0] == nhi[0] == -0.015625 -> kills L204 <= -> <)
            calls = {"n": 0}
            touch_val = math.nextafter(-0.015625, -math.inf)

            def boundary_touch_krawczyk(cfg, lo, hi):
                calls["n"] += 1
                if calls["n"] == 1:
                    return True, [va.arb(0), va.arb(0), va.arb(0)]
                if calls["n"] == 2:
                    return True, [va.arb(touch_val), va.arb(0), va.arb(0)]
                return True, None

            va.krawczyk = boundary_touch_krawczyk
            va.morse_index = lambda cfg, lo, hi: 2
            _check(run_cert(base) == (1, 0, 1), "boundary-touching shrink (nlo == nhi) should pass L204")

            # 5c. Exact boundary-touching enclosures in BOTH orders (e1[3][0] == e2[2][0] -> kills BOTH L215 < -> <=)
            z_left = {"box": [[-0.015625, -0.015625, -0.015625], [0.0, 0.015625, 0.015625]], "index": 2}
            z_right = {"box": [[0.0, -0.015625, -0.015625], [0.015625, 0.015625, 0.015625]], "index": 2}
            va.krawczyk = lambda cfg, lo, hi: (True, None)
            va.morse_index = lambda cfg, lo, hi: 2
            for lbl, zs in (("touching LR", [z_left, z_right]), ("touching RL", [z_right, z_left])):
                c_touch = copy.deepcopy(base)
                c_touch["zeros"] = zs
                c_touch["count"] = 2
                expect_fail(c_touch, lbl, substr="overlap")
        finally:
            va.krawczyk = orig_krawczyk
            va.morse_index = orig_morse


if __name__ == "__main__":
    unittest.main()
