#!/usr/bin/env python3
"""Assertion-free R5 verification, mutation-killing & Zone-B adversarial test suite
for Foundation-0/point-charges Issue #4 (Three-Charge Maxwell Bound Exact Study).

100% pure Python 3 stdlib unittest, zero bare asserts, zero AST floats, zero float(),
zero .isdigit(), zero mathematical RNG.
"""
from __future__ import annotations

import ast
import copy
import ctypes
from fractions import Fraction
import io
import json
import os
from pathlib import Path
import stat
import sys
import tempfile
import threading
import unittest
from unittest import mock

_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

import point_charges_n3_audit as audit  # noqa: E402
import point_charges_n3_oracle as oracle  # noqa: E402


def _resign_payload(payload: dict[str, object]) -> dict[str, object]:
    p_sha = audit._sha256_canonical(payload)
    outer = {
        "schema_version": oracle.SCHEMA_VERSION,
        "track_id": oracle.TRACK_ID,
        "payload_sha256": p_sha,
        "payload": payload,
    }
    c_sha = audit._sha256_canonical(outer)
    return {
        "schema_version": oracle.SCHEMA_VERSION,
        "track_id": oracle.TRACK_ID,
        "payload_sha256": p_sha,
        "certificate_sha256": c_sha,
        "payload": payload,
    }


def _call_main(argv: list[str]) -> tuple[int, str, str]:
    out_buf = io.StringIO()
    err_buf = io.StringIO()
    with mock.patch.object(sys, "stdout", out_buf), mock.patch.object(sys, "stderr", err_buf):
        rc = audit.main(argv)
    return rc, out_buf.getvalue(), err_buf.getvalue()


class TestPointChargesN3R5(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.cert = audit.build_standalone_certificate()
        cls.payload = cls.cert["payload"]

    def _expect_err(self, code: str, fn, *args, **kwargs) -> oracle.PointChargesN3Error:
        with self.assertRaises(oracle.PointChargesN3Error) as ctx:
            fn(*args, **kwargs)
        self.assertEqual(ctx.exception.code, code)
        return ctx.exception

    # ------------------------------------------------------------------
    # 1. Static AST Hygiene (Zero Floats, Zero Bare Asserts, Zero .isdigit)
    # ------------------------------------------------------------------
    def test_01_static_ast_hygiene(self) -> None:
        files = (
            _THIS_DIR / "point_charges_n3_oracle.py",
            _THIS_DIR / "point_charges_n3_audit.py",
            _THIS_DIR / "test_point_charges_n3_r5.py",
        )
        for fp in files:
            src = fp.read_text(encoding="utf-8")
            self.assertNotIn("\r", src)
            tree = ast.parse(src, filename=str(fp))
            for node in ast.walk(tree):
                self.assertNotIsInstance(node, ast.Assert, f"bare assert in {fp.name}")
                if isinstance(node, ast.Constant):
                    self.assertIsNot(type(node.value), float, f"float literal in {fp.name}")
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                    self.assertNotEqual(node.func.id, "float", f"float() call in {fp.name}")
                if isinstance(node, ast.Attribute):
                    self.assertNotEqual(node.attr, "isdigit", f".isdigit() in {fp.name}")

    # ------------------------------------------------------------------
    # 2. Strict Primitives, Bit/Digit Caps, Formatters
    # ------------------------------------------------------------------
    def test_02_strict_primitives_and_formatters(self) -> None:
        oracle.require_ok(True, "OK", "ok")
        self._expect_err("BOOM", oracle.require_ok, False, "BOOM", "fail")
        self._expect_err("INVALID_CONDITION_TYPE", oracle.require_ok, 1, "X", "bad")
        self.assertIs(oracle.require_exact_bool(True), True)
        self.assertIs(oracle.require_exact_bool(False), False)
        self._expect_err("INVALID_BOOL_TYPE", oracle.require_exact_bool, 1)
        self._expect_err("INVALID_BOOL_TYPE", oracle.require_exact_bool, "true")

        max_ok = (1 << oracle.MAX_INT_BITS) - 1
        self.assertEqual(oracle.require_exact_int(max_ok, "m", lo=0, hi=max_ok), max_ok)
        self.assertEqual(oracle.require_exact_int(0, "z", lo=0, hi=0), 0)
        self._expect_err("INVALID_INT_TYPE", oracle.require_exact_int, True)
        self._expect_err("INT_BIT_LIMIT", oracle.require_exact_int, 1 << oracle.MAX_INT_BITS)
        self._expect_err("INVALID_BOUND_TYPE", oracle.require_exact_int, 5, lo="0")
        self._expect_err("INVALID_BOUND_TYPE", oracle.require_exact_int, 5, hi="9")
        self._expect_err("INT_OUT_OF_RANGE", oracle.require_exact_int, -1, lo=0)
        self._expect_err("INT_OUT_OF_RANGE", oracle.require_exact_int, 10, hi=9)

        self.assertEqual(oracle.require_exact_fraction(Fraction(max_ok, 1)), Fraction(max_ok, 1))
        self.assertEqual(oracle.require_exact_fraction(Fraction(1, max_ok)), Fraction(1, max_ok))
        self._expect_err("INVALID_FRACTION_TYPE", oracle.require_exact_fraction, 3)
        self._expect_err("RATIONAL_BIT_LIMIT", oracle.require_exact_fraction, Fraction(1 << oracle.MAX_RATIONAL_BITS, 1))
        self._expect_err("RATIONAL_BIT_LIMIT", oracle.require_exact_fraction, Fraction(1, 1 << oracle.MAX_RATIONAL_BITS))

        self.assertEqual(oracle.format_rational(Fraction(5, 1)), "5")
        self.assertEqual(oracle.format_rational(Fraction(-3, 4)), "-3/4")
        self.assertEqual(oracle.format_decimal(Fraction(1, 200), 2), "0.01")
        self.assertEqual(oracle.format_decimal(Fraction(1, 201), 2), "0.00")
        self.assertEqual(oracle.format_decimal(Fraction(-1, 10000), 2), "0.00")
        self.assertEqual(oracle.format_decimal(Fraction(-7, 2), 0), "-4")
        self.assertEqual(oracle.format_decimal(Fraction(7, 2), 0), "4")
        self.assertEqual(oracle.format_decimal(Fraction(0), 0), "0")
        self.assertEqual(oracle.format_decimal(Fraction(-1, 4), 4), "-0.2500")
        self._expect_err("INT_OUT_OF_RANGE", oracle.format_decimal, Fraction(1), -1)
        self._expect_err("INT_OUT_OF_RANGE", oracle.format_decimal, Fraction(1), 13)

    # ------------------------------------------------------------------
    # 3. Exact QInterval & 192-bit Outward-Rounded DInterval Arithmetic
    # ------------------------------------------------------------------
    def test_03_qinterval_and_dinterval(self) -> None:
        q_pt = oracle.QInterval(Fraction(3))
        self.assertEqual((q_pt.lo, q_pt.hi), (Fraction(3), Fraction(3)))
        self._expect_err("INVALID_INTERVAL", oracle.QInterval, Fraction(2), Fraction(1))

        a = oracle.QInterval(Fraction(2), Fraction(4))
        b = oracle.QInterval(Fraction(-4), Fraction(-2))
        c = oracle.QInterval(Fraction(-2), Fraction(3))
        d = oracle.QInterval(Fraction(-3), Fraction(2))
        self.assertEqual((a.add(b).lo, a.add(b).hi), (Fraction(-2), Fraction(2)))
        self.assertEqual((a.sub(b).lo, a.sub(b).hi), (Fraction(4), Fraction(8)))
        self.assertEqual((a.neg().lo, a.neg().hi), (Fraction(-4), Fraction(-2)))
        self.assertEqual((a.mul(b).lo, a.mul(b).hi), (Fraction(-16), Fraction(-4)))
        self.assertEqual((a.scale(Fraction(2)).lo, a.scale(Fraction(2)).hi), (Fraction(4), Fraction(8)))
        self.assertEqual((a.scale(Fraction(-2)).lo, a.scale(Fraction(-2)).hi), (Fraction(-8), Fraction(-4)))
        self.assertEqual((a.sq().lo, a.sq().hi), (Fraction(4), Fraction(16)))
        self.assertEqual((b.sq().lo, b.sq().hi), (Fraction(4), Fraction(16)))
        self.assertEqual((c.sq().lo, c.sq().hi), (Fraction(0), Fraction(9)))
        self.assertEqual((d.sq().lo, d.sq().hi), (Fraction(0), Fraction(9)))
        self.assertEqual((a.inv().lo, a.inv().hi), (Fraction(1, 4), Fraction(1, 2)))
        self.assertEqual((b.inv().lo, b.inv().hi), (Fraction(-1, 2), Fraction(-1, 4)))
        self._expect_err("INTERVAL_DIV_ZERO", c.inv)
        self._expect_err("INTERVAL_DIV_ZERO", oracle.QInterval(Fraction(0), Fraction(2)).inv)
        self._expect_err("INTERVAL_DIV_ZERO", oracle.QInterval(Fraction(-2), Fraction(0)).inv)

        # 192-bit DInterval exact outward rounding
        scale = 1 << oracle.DYADIC_PREC_BITS
        d_third = oracle.DInterval.from_fraction(Fraction(1, 3))
        self.assertLess(Fraction(d_third.lo, scale), Fraction(1, 3))
        self.assertGreater(Fraction(d_third.hi, scale), Fraction(1, 3))
        self.assertEqual(d_third.hi - d_third.lo, 1)
        d_eq = oracle.DInterval.from_box(Fraction(2), Fraction(2))
        self.assertEqual((d_eq.lo, d_eq.hi), (2 << oracle.DYADIC_PREC_BITS, 2 << oracle.DYADIC_PREC_BITS))
        self._expect_err("INVALID_INTERVAL", oracle.DInterval.from_box, Fraction(2), Fraction(1))

        d_pos = oracle.DInterval.from_box(Fraction(2), Fraction(3))
        d_neg = oracle.DInterval.from_box(Fraction(-3), Fraction(-2))
        d_str1 = oracle.DInterval.from_box(Fraction(-2), Fraction(3))
        d_str2 = oracle.DInterval.from_box(Fraction(-3), Fraction(2))
        self.assertEqual((d_pos.sq().lo, d_pos.sq().hi), (4 << oracle.DYADIC_PREC_BITS, 9 << oracle.DYADIC_PREC_BITS))
        self.assertEqual((d_neg.sq().lo, d_neg.sq().hi), (4 << oracle.DYADIC_PREC_BITS, 9 << oracle.DYADIC_PREC_BITS))
        self.assertEqual((d_str1.sq().lo, d_str1.sq().hi), (0, 9 << oracle.DYADIC_PREC_BITS))
        self.assertEqual((d_str2.sq().lo, d_str2.sq().hi), (0, 9 << oracle.DYADIC_PREC_BITS))
        self.assertEqual(d_str1.abs_hi(), 3 << oracle.DYADIC_PREC_BITS)
        self.assertEqual(d_str2.abs_hi(), 3 << oracle.DYADIC_PREC_BITS)
        self.assertEqual(d_pos.midpoint_fraction(), Fraction(5, 2))
        self.assertEqual(d_pos.neg().lo, -(3 << oracle.DYADIC_PREC_BITS))
        self.assertEqual(d_pos.neg().hi, -(2 << oracle.DYADIC_PREC_BITS))

        sq4 = oracle.DInterval.from_fraction(Fraction(4)).sqrt_pos()
        self.assertEqual((sq4.lo, sq4.hi), (2 << oracle.DYADIC_PREC_BITS, 2 << oracle.DYADIC_PREC_BITS))
        sq2 = oracle.DInterval.from_fraction(Fraction(2)).sqrt_pos()
        self.assertLess(Fraction(sq2.lo, scale) ** 2, Fraction(2))
        self.assertGreater(Fraction(sq2.hi, scale) ** 2, Fraction(2))
        self.assertEqual(sq2.hi - sq2.lo, 1)
        self._expect_err("SQRT_NONPOSITIVE", oracle.DInterval.from_box(Fraction(0), Fraction(1)).sqrt_pos)

        inv3 = oracle.DInterval.from_fraction(Fraction(3)).inv_pos()
        self.assertLess(Fraction(inv3.lo, scale), Fraction(1, 3))
        self.assertGreater(Fraction(inv3.hi, scale), Fraction(1, 3))
        self.assertEqual(inv3.hi - inv3.lo, 1)
        self._expect_err("INV_NONPOSITIVE", oracle.DInterval.from_box(Fraction(0), Fraction(1)).inv_pos)

        # Multiplication outward rounding check
        prod = d_third.mul(d_third)
        self.assertLess(Fraction(prod.lo, scale), Fraction(1, 9))
        self.assertGreater(Fraction(prod.hi, scale), Fraction(1, 9))

    # ------------------------------------------------------------------
    # 4. Mathematical Operators & Negative Controls (All 5 Pillars)
    # ------------------------------------------------------------------
    def test_04_oracle_operators_and_negative_controls(self) -> None:
        # eval_exact_pmin_and_dl positive & all 3 barycentric failure axes + degenerate triangle
        ev = oracle.eval_exact_pmin_and_dl(Fraction(1), Fraction(1), Fraction(1), Fraction(1, 3), Fraction(1, 3))
        self.assertEqual(ev["pmin"], Fraction(1, 216))
        self.assertEqual(ev["det_h_poly"], Fraction(1, 4))
        self.assertEqual(ev["det_h_phase"], Fraction(1, 4))
        self.assertEqual(ev["det_dl"], Fraction(27, 4))
        self.assertEqual(ev["det_dl_pred"], Fraction(27, 4))
        self.assertEqual(ev["dl_matrix"], ((Fraction(-3), Fraction(-3, 2)), (Fraction(-3, 2), Fraction(-3))))

        self._expect_err("OUTSIDE_TRIANGLE", oracle.eval_exact_pmin_and_dl, Fraction(1), Fraction(1), Fraction(1), Fraction(0), Fraction(1, 3))
        self._expect_err("OUTSIDE_TRIANGLE", oracle.eval_exact_pmin_and_dl, Fraction(1), Fraction(1), Fraction(1), Fraction(1, 3), Fraction(0))
        self._expect_err("OUTSIDE_TRIANGLE", oracle.eval_exact_pmin_and_dl, Fraction(1), Fraction(1), Fraction(1), Fraction(2, 3), Fraction(1, 3))
        self._expect_err("DEGENERATE_TRIANGLE", oracle.eval_exact_pmin_and_dl, Fraction(1), Fraction(1), Fraction(4), Fraction(1, 3), Fraction(1, 3))
        self._expect_err("DEGENERATE_TRIANGLE", oracle.eval_exact_pmin_and_dl, Fraction(1), Fraction(1), Fraction(9), Fraction(1, 3), Fraction(1, 3))

        # weierstrass_rational_triangle positive & negative controls
        wt = oracle.weierstrass_rational_triangle(Fraction(1, 3), Fraction(1, 3))
        self.assertEqual(wt["v3"], Fraction(1, 7))
        self.assertEqual(wt["grad_incenter"], (Fraction(0), Fraction(0)))
        self.assertEqual(wt["det_h_tilde_incenter"], Fraction(1, 4))
        self._expect_err("INVALID_WEIERSTRASS_PARAM", oracle.weierstrass_rational_triangle, Fraction(0), Fraction(1, 3))
        self._expect_err("INVALID_WEIERSTRASS_PARAM", oracle.weierstrass_rational_triangle, Fraction(1, 3), Fraction(0))
        self._expect_err("INVALID_WEIERSTRASS_PARAM", oracle.weierstrass_rational_triangle, Fraction(1, 3), Fraction(1, 2))
        self._expect_err("INVALID_WEIERSTRASS_PARAM", oracle.weierstrass_rational_triangle, Fraction(1, 2), Fraction(1, 2))

        # certify_2d_krawczyk_box negative controls (all 5 error paths)
        sites = wt["sites"]
        charges = wt["charges"]
        self._expect_err("INVALID_HALF_WIDTH", oracle.certify_2d_krawczyk_box, sites, charges, Fraction(0), Fraction(0), Fraction(0))
        self._expect_err(
            "SINGULAR_PRECONDITIONER",
            oracle.certify_2d_krawczyk_box,
            sites,
            (Fraction(0), Fraction(0), Fraction(0)),
            Fraction(0),
            Fraction(0),
            Fraction(1, 1 << 20),
        )
        # Off-equilibrium tiny box fails KRAWCZYK_NOT_INTERIOR
        self._expect_err(
            "KRAWCZYK_NOT_INTERIOR",
            oracle.certify_2d_krawczyk_box,
            sites,
            charges,
            Fraction(1, 10),
            Fraction(0),
            Fraction(1, 1 << 20),
        )
        # Wide box around (0,0) fails NONPOSITIVE_TRACE
        self._expect_err(
            "NONPOSITIVE_TRACE",
            oracle.certify_2d_krawczyk_box,
            sites,
            charges,
            Fraction(0),
            Fraction(0),
            Fraction(1, 5),
        )
        # Intermediate box fails DEGENERATE_HESSIAN_BOX
        self._expect_err(
            "DEGENERATE_HESSIAN_BOX",
            oracle.certify_2d_krawczyk_box,
            sites,
            charges,
            Fraction(0),
            Fraction(0),
            Fraction(1, 32),
        )
        # Tight-but-too-wide box fails KRAWCZYK_NO_CONTRACTION
        self._expect_err(
            "KRAWCZYK_NO_CONTRACTION",
            oracle.certify_2d_krawczyk_box,
            sites,
            charges,
            Fraction(0),
            Fraction(0),
            Fraction(1, 36),
        )
        # Exact boundary tr_h.lo == 0 and rho_scaled == 1 << DYADIC_PREC_BITS
        d0 = oracle.DInterval(0, 0)
        d1 = oracle.DInterval(1 << oracle.DYADIC_PREC_BITS, 1 << oracle.DYADIC_PREC_BITS)
        dm1 = oracle.DInterval(-(1 << oracle.DYADIC_PREC_BITS), -(1 << oracle.DYADIC_PREC_BITS))
        d2 = oracle.DInterval(2 << oracle.DYADIC_PREC_BITS, 2 << oracle.DYADIC_PREC_BITS)
        with mock.patch.object(oracle, "_eval_2d_coulomb_dyadic", side_effect=[(d0, d0, d1, d0, dm1), (d0, d0, d1, d0, dm1)]):
            self._expect_err("NONPOSITIVE_TRACE", oracle.certify_2d_krawczyk_box, sites, charges, Fraction(0), Fraction(0), Fraction(1, 1 << 20))
        with mock.patch.object(oracle, "_eval_2d_coulomb_dyadic", side_effect=[(d0, d0, d1, d0, d1), (d0, d0, d2, d0, d1)]):
            self._expect_err("KRAWCZYK_NO_CONTRACTION", oracle.certify_2d_krawczyk_box, sites, charges, Fraction(0), Fraction(0), Fraction(1, 1 << 20))

        # _boxes_strictly_disjoint: all 4 separating inequalities + touching/overlapping false
        b0 = ((Fraction(0), Fraction(0)), (Fraction(1), Fraction(1)))
        self.assertTrue(oracle._boxes_strictly_disjoint(b0[0], b0[1], (Fraction(2), Fraction(0)), (Fraction(3), Fraction(1))))
        self.assertTrue(oracle._boxes_strictly_disjoint((Fraction(2), Fraction(0)), (Fraction(3), Fraction(1)), b0[0], b0[1]))
        self.assertTrue(oracle._boxes_strictly_disjoint(b0[0], b0[1], (Fraction(0), Fraction(2)), (Fraction(1), Fraction(3))))
        self.assertTrue(oracle._boxes_strictly_disjoint((Fraction(0), Fraction(2)), (Fraction(1), Fraction(3)), b0[0], b0[1]))
        self.assertFalse(oracle._boxes_strictly_disjoint(b0[0], b0[1], (Fraction(1), Fraction(0)), (Fraction(2), Fraction(1))))
        self.assertFalse(oracle._boxes_strictly_disjoint(b0[0], b0[1], (Fraction(1, 2), Fraction(1, 2)), (Fraction(2), Fraction(2))))

        # _check_4d_gale_nikaido_box negative control outside R_min
        ok_bad, pm_bad, _s00_bad, _ds_bad = oracle._check_4d_gale_nikaido_box(
            oracle.QInterval(Fraction(1)),
            oracle.QInterval(Fraction(1)),
            oracle.QInterval(Fraction(1)),
            oracle.QInterval(Fraction(1, 10), Fraction(1, 9)),
            oracle.QInterval(Fraction(1, 10), Fraction(1, 9)),
            ((Fraction(4, 9), Fraction(-2, 9)), (Fraction(-2, 9), Fraction(4, 9))),
        )
        self.assertFalse(ok_bad)
        self.assertLess(pm_bad, 0)

    # ------------------------------------------------------------------
    # 5. Exact Five-Pillar Payload Verification
    # ------------------------------------------------------------------
    def test_05_five_pillars_payload_exact_values(self) -> None:
        pl = self.payload
        self.assertEqual(pl["track_id"], oracle.TRACK_ID)
        self.assertEqual(pl["issue_ref"], "Foundation-0/point-charges#4")
        self.assertEqual(pl["novelty_class"], "N1")

        p1 = pl["pillar1_algebraic_identities"]
        self.assertTrue(p1["all_verified"])
        self.assertEqual(p1["out_of_plane_vz_probe_q"], "-12/125")
        self.assertEqual(len(p1["algebraic_identity_cases"]), 7)
        self.assertEqual([r["in_R_min"] for r in p1["algebraic_identity_cases"]], [True, True, True, True, True, True, False])

        p2 = pl["pillar2_incenter_universal_theorem"]
        self.assertTrue(p2["all_verified"])
        self.assertTrue(p2["closes_487_flat_triangle_nonemptiness_gap"])
        self.assertEqual(len(p2["weierstrass_incenter_witnesses"]), 6)
        self.assertEqual(p2["weierstrass_incenter_witnesses"][4]["aspect_h_over_max_side_dec6"], "0.000020")

        p3 = pl["pillar3_isosceles_axis_uniqueness"]
        self.assertTrue(p3["all_verified"])
        self.assertEqual([r["k"] for r in p3["isosceles_axis_cases"]], [2, 3, 10, 100])
        self.assertEqual(p3["isosceles_axis_cases"][0]["height_h_q"], "3")
        self.assertEqual(p3["isosceles_axis_cases"][1]["height_h_q"], "9/7")

        p4 = pl["pillar4_gale_nikaido_compact_univalence"]
        self.assertTrue(p4["all_verified"])
        self.assertEqual(p4["total_4d_subboxes_certified"], 128)
        self.assertEqual(p4["compact_4d_boxes"][0]["min_pmin_lo_dec6"], "0.000576")
        self.assertEqual(p4["compact_4d_boxes"][0]["min_det_sym_lo_dec6"], "0.594345")
        self.assertEqual(p4["compact_4d_boxes"][1]["min_pmin_lo_dec6"], "2.782278")
        self.assertEqual(p4["compact_4d_boxes"][1]["min_det_sym_lo_dec6"], "0.670078")

        p5 = pl["pillar5_krawczyk_configurations"]
        self.assertTrue(p5["all_verified"])
        self.assertEqual(p5["total_equilibria_certified"], 14)
        self.assertEqual([c["equilibria_certified"] for c in p5["configurations"]], [4, 4, 4, 2])
        self.assertEqual([c["m1_local_minima"] for c in p5["configurations"]], [1, 1, 1, 0])
        self.assertEqual([c["m2_saddles"] for c in p5["configurations"]], [3, 3, 3, 2])
        self.assertEqual([c["pairwise_disjoint_pairs"] for c in p5["configurations"]], [6, 6, 6, 1])

        summ = pl["summary"]
        self.assertTrue(summ["all_five_pillars_verified"])
        self.assertTrue(summ["global_n3_conjecture_open_outside_compact_boxes"])
        self.assertEqual(summ["float_count_on_trust_path"], 0)
        self.assertEqual(summ["bare_assert_count"], 0)
        self.assertEqual(summ["rng_calls"], 0)

    # ------------------------------------------------------------------
    # 6. Strict Parsers, Aligned 3,321-Bit / 1,000-Digit Cap & Canonical JSON
    # ------------------------------------------------------------------
    def test_06_parsers_aligned_caps_and_canonical_json(self) -> None:
        self.assertEqual(audit.require_hex64_lower_ascii("a" * 64), "a" * 64)
        self._expect_err("INVALID_HEX64", audit.require_hex64_lower_ascii, "a" * 63)
        self._expect_err("INVALID_HEX64", audit.require_hex64_lower_ascii, "A" * 64)
        self._expect_err("INVALID_HEX64", audit.require_hex64_lower_ascii, "g" * 64)
        self._expect_err("INVALID_STR", audit.require_hex64_lower_ascii, 123)

        self.assertEqual(audit.parse_exact_rational("0"), Fraction(0))
        self.assertEqual(audit.parse_exact_rational("-5"), Fraction(-5))
        self.assertEqual(audit.parse_exact_rational("-7/11"), Fraction(-7, 11))
        for bad in ("-0", "01", "", "-", "1/0", "1/1", "2/4", "1/2/3", "\u0661", "1/02", "9" * 1001):
            with self.assertRaises(oracle.PointChargesN3Error):
                audit.parse_exact_rational(bad)

        # Aligned 3,321-bit boundary round-trip (1,000 digits)
        max_ok = (1 << oracle.MAX_INT_BITS) - 1
        self.assertEqual(len(str(max_ok)), oracle.MAX_INT_DIGITS)
        rt = audit.parse_strict_json_text(f'{{"x":{max_ok}}}')
        self.assertEqual(rt["x"], max_ok)
        self.assertEqual(audit.canonical_json_bytes(rt), f'{{"x":{max_ok}}}'.encode("ascii"))

        self._expect_err("INT_BIT_LIMIT", audit._parse_bounded_int, str(1 << oracle.MAX_INT_BITS))
        self._expect_err("INT_BIT_LIMIT", audit.canonical_json_bytes, {"x": 1 << oracle.MAX_INT_BITS})
        self._expect_err("INT_BIT_LIMIT", audit.parse_strict_json_text, f'{{"x":{1 << oracle.MAX_INT_BITS}}}')
        self._expect_err("INT_DIGIT_LIMIT", audit.parse_strict_json_text, '{"x":1' + "0" * 4300 + "}")
        self._expect_err("FLOAT_LITERAL_FORBIDDEN", audit.parse_strict_json_text, '{"x":1.5}')
        self._expect_err("FLOAT_LITERAL_FORBIDDEN", audit.parse_strict_json_text, '{"x":1e2}')
        self._expect_err("NONFINITE_LITERAL_FORBIDDEN", audit.parse_strict_json_text, '{"x":NaN}')
        self._expect_err("DUPLICATE_JSON_KEY", audit.parse_strict_json_text, '{"a":1,"a":2}')
        self._expect_err("JSON_PARSE_ERROR", audit.parse_strict_json_text, '{"a":')
        self._expect_err("JSON_ROOT_NOT_OBJECT", audit.parse_strict_json_text, "[1,2]")
        self.assertEqual(audit.parse_strict_json_text("[1,null,true,false]", require_object=False), [1, None, True, False])
        self.assertEqual(audit.parse_strict_json_text('{"a":1}', max_bytes=7), {"a": 1})
        self._expect_err("JSON_TOO_LARGE", audit.parse_strict_json_text, '{"a":"0123456789"}', max_bytes=10)
        self._expect_err("INT_OUT_OF_RANGE", audit.parse_strict_json_text, "{}", max_bytes=1 << 100)
        self._expect_err("NUL_BYTE_FORBIDDEN", audit.parse_strict_json_text, '{"a":"\x00"}')
        self._expect_err("INVALID_UTF8_SURROGATE", audit.parse_strict_json_text, '{"a":"\ud800"}')
        self._expect_err("FORBIDDEN_JSON_TYPE", audit.canonical_json_bytes, {"a": {1, 2}})

        # Exact depth limit 64 passes, 65 fails
        d64 = "[" * 64 + "0" + "]" * 64
        d65 = "[" * 65 + "0" + "]" * 65
        self.assertIsInstance(audit.parse_strict_json_text(d64, require_object=False), list)
        self._expect_err("JSON_DEPTH_LIMIT", audit.parse_strict_json_text, d65, require_object=False)

        # Canonical ASCII & key sorting
        self.assertEqual(audit.canonical_json_bytes({"b": "\u00e9", "a": 1}), b'{"a":1,"b":"\\u00e9"}')

    # ------------------------------------------------------------------
    # 7. Deep Secondary Auditor & Re-Signed full_replay=False Tamper Matrix
    # ------------------------------------------------------------------
    def test_07_secondary_audit_and_resigned_tamper_matrix(self) -> None:
        v_full = audit.verify_certificate_dict(self.cert, full_replay=True)
        self.assertTrue(v_full["valid"])
        self.assertTrue(v_full["full_replay"])
        v_fast = audit.verify_certificate_dict(self.cert, full_replay=False)
        self.assertTrue(v_fast["valid"])
        self.assertFalse(v_fast["full_replay"])
        self._expect_err("INVALID_BOOL_TYPE", audit.verify_certificate_dict, self.cert, full_replay=1)

        # Envelope tampers
        self._expect_err("INVALID_CERT_TYPE", audit.verify_certificate_dict, [])
        c_extra = copy.deepcopy(self.cert)
        c_extra["rogue"] = 1
        self._expect_err("INVALID_CERT_KEYS", audit.verify_certificate_dict, c_extra, full_replay=False)
        c_sch = copy.deepcopy(self.cert)
        c_sch["schema_version"] = "bad"
        self._expect_err("INVALID_SCHEMA", audit.verify_certificate_dict, c_sch, full_replay=False)
        c_tid = copy.deepcopy(self.cert)
        c_tid["track_id"] = "bad"
        self._expect_err("INVALID_TRACK_ID", audit.verify_certificate_dict, c_tid, full_replay=False)
        c_pt = copy.deepcopy(self.cert)
        c_pt["payload"] = []
        self._expect_err("INVALID_PAYLOAD_TYPE", audit.verify_certificate_dict, c_pt, full_replay=False)
        c_psha = copy.deepcopy(self.cert)
        c_psha["payload_sha256"] = "0" * 64
        self._expect_err("PAYLOAD_SHA256_MISMATCH", audit.verify_certificate_dict, c_psha, full_replay=False)
        c_csha = copy.deepcopy(self.cert)
        c_csha["certificate_sha256"] = "0" * 64
        self._expect_err("CERT_SHA256_MISMATCH", audit.verify_certificate_dict, c_csha, full_replay=False)

        # Re-signed semantic tamper matrix under full_replay=False
        def mutate_and_expect(mutator, exp_code: str) -> None:
            pl = copy.deepcopy(self.payload)
            mutator(pl)
            resigned = _resign_payload(pl)
            self._expect_err(exp_code, audit.verify_certificate_dict, resigned, full_replay=False)

        mutate_and_expect(lambda p: p.__setitem__("track_id", "WRONG"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p.__setitem__("issue_ref", "WRONG"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p.__setitem__("novelty_class", "N0"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p.__setitem__("extra", 1), "INVALID_SCHEMA_KEYS")
        mutate_and_expect(lambda p: p["pillar1_algebraic_identities"].__setitem__("out_of_plane_vz_probe_q", "-1/10"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar1_algebraic_identities"]["algebraic_identity_cases"].pop(), "INVALID_SCHEMA_LIST")
        mutate_and_expect(lambda p: p["pillar1_algebraic_identities"].__setitem__("algebraic_identity_cases", {str(k): k for k in range(7)}), "INVALID_SCHEMA_LIST")
        self._expect_err("INVALID_SCHEMA_DICT", audit._req_dict, [], "row")
        self._expect_err("MISSING_FIELD", audit._req_int_f, {}, "k", "r")
        self._expect_err("MISSING_FIELD", audit._req_bool_f, {}, "b", "r")
        mutate_and_expect(lambda p: p["pillar1_algebraic_identities"]["algebraic_identity_cases"][0].__setitem__("verified", False), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar1_algebraic_identities"]["algebraic_identity_cases"][0].__setitem__("in_R_min", 1), "INVALID_BOOL_TYPE")
        mutate_and_expect(lambda p: p["pillar1_algebraic_identities"]["algebraic_identity_cases"][0].__setitem__("in_R_min", False), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar1_algebraic_identities"]["algebraic_identity_cases"][0].pop("verified"), "INVALID_SCHEMA_KEYS")
        mutate_and_expect(lambda p: p["pillar2_incenter_universal_theorem"]["weierstrass_incenter_witnesses"][0].__setitem__("det_h_tilde_q", "1/5"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar3_isosceles_axis_uniqueness"]["isosceles_axis_cases"][0].__setitem__("k", 99), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar3_isosceles_axis_uniqueness"]["isosceles_axis_cases"][0].__setitem__("y_pitchfork_q", "3/2"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar3_isosceles_axis_uniqueness"]["isosceles_axis_cases"][0].__setitem__("y_pitchfork_eq_2_y_fold", 1), "INVALID_BOOL_TYPE")
        mutate_and_expect(lambda p: p["pillar3_isosceles_axis_uniqueness"]["isosceles_axis_cases"][0].pop("k"), "INVALID_SCHEMA_KEYS")
        mutate_and_expect(lambda p: p["pillar4_gale_nikaido_compact_univalence"]["compact_4d_boxes"][0].__setitem__("subboxes_passed", 63), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar4_gale_nikaido_compact_univalence"].__setitem__("total_4d_subboxes_certified", 127), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar5_krawczyk_configurations"]["configurations"][0]["equilibria"][0].__setitem__("det_h_sign", True), "INVALID_INT_TYPE")
        mutate_and_expect(lambda p: p["pillar5_krawczyk_configurations"]["configurations"][0]["equilibria"][0].__setitem__("det_h_sign", -1), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar5_krawczyk_configurations"]["configurations"][0].__setitem__("pairwise_disjoint_pairs", True), "INVALID_INT_TYPE")
        mutate_and_expect(lambda p: p["pillar5_krawczyk_configurations"]["configurations"][0]["equilibria"][0].__setitem__("morse_index_3d", 2), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar5_krawczyk_configurations"]["configurations"][0]["equilibria"][0].__setitem__("morse_index_3d", 3), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar5_krawczyk_configurations"]["configurations"][0]["equilibria"][0].__setitem__("krawczyk_contraction_lt_1", False), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar5_krawczyk_configurations"]["configurations"][0].__setitem__("m2_minus_m1", 1), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["summary"].__setitem__("float_count_on_trust_path", 1), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["summary"].__setitem__("pillar_count", True), "INVALID_INT_TYPE")
        mutate_and_expect(lambda p: p["summary"].__setitem__("all_five_pillars_verified", 1), "INVALID_BOOL_TYPE")
        mutate_and_expect(lambda p: p["summary"].__setitem__("global_n3_conjecture_open_outside_compact_boxes", False), "SECONDARY_AUDIT_FAILED")

        # Shadow semantic tamper matrix (16 re-signed value/display/geometry/Krawczyk mutations caught under full_replay=False)
        mutate_and_expect(lambda p: p["pillar1_algebraic_identities"]["algebraic_identity_cases"][0].__setitem__("name", "tampered_name"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar1_algebraic_identities"]["algebraic_identity_cases"][0].__setitem__("det_h_tilde_q", "1/2"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar1_algebraic_identities"]["algebraic_identity_cases"][0].__setitem__("det_h_tilde_dec6", "0.250001"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar1_algebraic_identities"]["algebraic_identity_cases"][0].__setitem__("det_dl_q", "27/5"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar2_incenter_universal_theorem"]["weierstrass_incenter_witnesses"][0].__setitem__("name", "wrong_witness"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar2_incenter_universal_theorem"]["weierstrass_incenter_witnesses"][0].__setitem__("v1_q", "1/4"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar2_incenter_universal_theorem"]["weierstrass_incenter_witnesses"][0].__setitem__("aspect_h_over_max_side_dec6", "0.960001"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar3_isosceles_axis_uniqueness"]["isosceles_axis_cases"][0].__setitem__("height_h_q", "4"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar3_isosceles_axis_uniqueness"]["isosceles_axis_cases"][0].__setitem__("y_fold_q", "1/3"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar4_gale_nikaido_compact_univalence"]["compact_4d_boxes"][0].__setitem__("a2_range_q", ["1", "2"]), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar4_gale_nikaido_compact_univalence"]["compact_4d_boxes"][0].__setitem__("min_pmin_lo_dec6", "0.000500"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar5_krawczyk_configurations"]["configurations"][0].__setitem__("charges_q", ["1", "1", "1"]), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar5_krawczyk_configurations"]["configurations"][0]["equilibria"][0].__setitem__("center_q", ["1/10", "0"]), "KRAWCZYK_NOT_INTERIOR")
        mutate_and_expect(lambda p: p["pillar5_krawczyk_configurations"]["configurations"][0]["equilibria"][0].__setitem__("box_lo_q", ["0", "0"]), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar5_krawczyk_configurations"]["configurations"][0]["equilibria"][0].__setitem__("box_hi_q", ["0", "0"]), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar5_krawczyk_configurations"]["configurations"][0]["equilibria"][0].__setitem__("krawczyk_rho_ub_dec6", "0.000099"), "SECONDARY_AUDIT_FAILED")
        mutate_and_expect(lambda p: p["pillar5_krawczyk_configurations"]["configurations"][0]["equilibria"][0].__setitem__("half_width_q", "1/2097152"), "SECONDARY_AUDIT_FAILED")

        # Full replay guard independently rejects digest mismatch against recomputed certificate
        with mock.patch.object(audit, "build_standalone_certificate", return_value={"payload_sha256": "0" * 64}):
            self._expect_err("REPLAY_MISMATCH", audit.verify_certificate_dict, self.cert, full_replay=True)

    # ------------------------------------------------------------------
    # 8. Zone-B File I/O, Symlink/PWD/Hardlink/FIFO, Streaming Cap & CLI
    # ------------------------------------------------------------------
    def test_08_zone_b_io_and_cli(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td).resolve()
            cert_path = root / "sub" / "cert.json"
            wres = audit.write_certificate_file(cert_path)
            self.assertTrue(wres["written"])
            self.assertEqual(stat.S_IMODE(os.lstat(cert_path).st_mode), 0o600)
            cres = audit.check_certificate_file(cert_path, full_replay=True)
            self.assertTrue(cres["valid"])
            self.assertEqual(cres["file_sha256"], wres["file_sha256"])

            # Missing parent read MUST fail without creating directory
            missing_dir = root / "never_created"
            self._expect_err("FILE_READ_ERROR", audit.read_regular_bytes_capped, missing_dir / "x.json")
            self.assertFalse(missing_dir.exists())

            # Non-canonical, CRLF, missing LF, bad UTF-8
            crlf_p = root / "crlf.json"
            crlf_p.write_bytes(cert_path.read_bytes().replace(b"\n", b"\r\n"))
            self._expect_err("CRLF_FORBIDDEN", audit.check_certificate_file, crlf_p)

            nolf_p = root / "nolf.json"
            nolf_p.write_bytes(cert_path.read_bytes()[:-1])
            self._expect_err("MISSING_TRAILING_LF", audit.check_certificate_file, nolf_p)

            badu_p = root / "badu.json"
            badu_p.write_bytes(b"\xff\xfe\n")
            self._expect_err("INVALID_UTF8", audit.check_certificate_file, badu_p)

            pretty_p = root / "pretty.json"
            pretty_p.write_bytes((json.dumps(self.cert, indent=2) + "\n").encode("utf-8"))
            self._expect_err("NON_CANONICAL_BYTES", audit.check_certificate_file, pretty_p)

            # Symlink leaf, symlink ancestor, dir target, root, "..", FIFO, /dev/null, hardlink, 70KB streaming cap
            sym_leaf = root / "sym_leaf.json"
            os.symlink(cert_path, sym_leaf)
            self._expect_err("SYMLINK_FORBIDDEN", audit.read_regular_bytes_capped, sym_leaf)
            self._expect_err("SYMLINK_FORBIDDEN", audit.write_regular_bytes_safe, sym_leaf, b"x")

            sym_dir = root / "sym_dir"
            os.symlink(root / "sub", sym_dir)
            self._expect_err("SYMLINK_FORBIDDEN", audit.read_regular_bytes_capped, sym_dir / "cert.json")
            self._expect_err("SYMLINK_FORBIDDEN", audit.write_regular_bytes_safe, sym_dir / "new.json", b"x")

            self._expect_err("NOT_A_REGULAR_FILE", audit.read_regular_bytes_capped, root / "sub")
            self._expect_err("NOT_A_REGULAR_FILE", audit.write_regular_bytes_safe, root / "sub", b"x")
            self._expect_err("NOT_A_REGULAR_FILE", audit.read_regular_bytes_capped, Path("/"))
            self._expect_err("NOT_A_REGULAR_FILE", audit.read_regular_bytes_capped, missing_dir / "..")
            self._expect_err("NOT_A_REGULAR_FILE", audit.read_regular_bytes_capped, Path("/dev/null"))

            fifo_p = root / "pipe.fifo"
            os.mkfifo(fifo_p)
            self._expect_err("NOT_A_REGULAR_FILE", audit.read_regular_bytes_capped, fifo_p)

            hl_p = root / "hardlink.json"
            os.link(cert_path, hl_p)
            self._expect_err("HARDLINK_FORBIDDEN", audit.write_regular_bytes_safe, hl_p, b"x")
            os.unlink(hl_p)

            big_p = root / "big70k.bin"
            big_p.write_bytes(b"A" * 70_000)
            self.assertEqual(len(audit.read_regular_bytes_capped(big_p, "exact", 70_000)), 70_000)
            self._expect_err("FILE_TOO_LARGE", audit.read_regular_bytes_capped, big_p, "big", 65_536)
            if Path("/proc/kallsyms").exists():
                self._expect_err("FILE_TOO_LARGE", audit.read_regular_bytes_capped, Path("/proc/kallsyms"), "ks", 1024)

            # Pre-replace fstat st_nlink != 1 guard
            fake_st = os.stat_result((stat.S_IFREG | 0o600, 1, 1, 2, 0, 0, 1, 0, 0, 0))
            with mock.patch.object(os, "fstat", return_value=fake_st):
                self._expect_err("HARDLINK_FORBIDDEN", audit.write_regular_bytes_safe, root / "race_nlink.bin", b"x")

            self._expect_err("INVALID_PATH_TYPE", audit.read_regular_bytes_capped, str(cert_path))  # type: ignore[arg-type]
            self._expect_err("INVALID_PAYLOAD_BYTES", audit.write_regular_bytes_safe, root / "t.bin", "str")  # type: ignore[arg-type]

            # Relative path + symlinked PWD + deleted CWD checks
            orig_cwd = os.getcwd()
            orig_pwd = os.environ.get("PWD")
            try:
                os.chdir(root)
                os.environ["PWD"] = str(root)
                rel_p = Path("rel_out.json")
                audit.write_regular_bytes_safe(rel_p, b"{}\n")
                self.assertEqual(audit.read_regular_bytes_capped(rel_p), b"{}\n")

                os.environ["PWD"] = str(sym_dir)
                self._expect_err("SYMLINK_ANCESTOR_FORBIDDEN", audit.read_regular_bytes_capped, Path("rel_out.json"))

                ghost = root / "ghost_dir"
                ghost.mkdir()
                os.chdir(ghost)
                os.environ["PWD"] = str(ghost)
                ghost.rmdir()
                self._expect_err("CWD_UNAVAILABLE", audit.read_regular_bytes_capped, Path("rel.json"))
            finally:
                os.chdir(orig_cwd)
                if orig_pwd is None:
                    os.environ.pop("PWD", None)
                else:
                    os.environ["PWD"] = orig_pwd

            # CLI tests (space & equals syntax + errors)
            rc, out, err = _call_main([])
            self.assertEqual((rc, err), (0, ""))
            self.assertIn("OK ", out)

            rc, out, err = _call_main(["--standalone", "--json"])
            self.assertEqual((rc, err), (0, ""))
            self.assertEqual(json.loads(out)["certificate_sha256"], self.cert["certificate_sha256"])

            cli_w = root / "cli_w.json"
            rc, out, err = _call_main(["--write", str(cli_w), "--json"])
            self.assertEqual((rc, err), (0, ""))
            self.assertTrue(json.loads(out)["written"])

            rc, out, err = _call_main(["--check", str(cli_w), "--json"])
            self.assertEqual((rc, err), (0, ""))
            self.assertTrue(json.loads(out)["valid"])

            cli_w2 = root / "cli_w2.json"
            rc, out, err = _call_main([f"--write={cli_w2}"])
            self.assertEqual((rc, err), (0, ""))
            rc, out, err = _call_main([f"--check={cli_w2}", "--json"])
            self.assertEqual((rc, err), (0, ""))
            self.assertTrue(json.loads(out)["valid"])

            for bad_argv, exp_code in (
                (["--unknown"], "UNKNOWN_FLAG"),
                (["--json", "--json"], "DUPLICATE_FLAG"),
                (["--standalone", "--check", str(cli_w)], "CONFLICTING_MODE"),
                (["--check", str(cli_w), "--standalone"], "CONFLICTING_MODE"),
                (["--standalone", f"--check={cli_w}"], "CONFLICTING_MODE"),
                (["--check"], "MISSING_PATH_ARG"),
                (["--check", ""], "INVALID_PATH_ARG"),
                (["--check", "-bad"], "INVALID_PATH_ARG"),
                (["--check="], "INVALID_PATH_ARG"),
                (["--write=-bad"], "INVALID_PATH_ARG"),
                (["--check", "a\x00b"], "NUL_BYTE_FORBIDDEN"),
            ):
                rc, _o, e = _call_main(bad_argv)
                self.assertEqual(rc, 2)
                self.assertEqual(json.loads(e)["error_code"], exp_code)

            rc, _o, e = _call_main("not-a-list")  # type: ignore[arg-type]
            self.assertEqual(rc, 2)
            self.assertEqual(json.loads(e)["error_code"], "INVALID_ARGV")

    # ------------------------------------------------------------------
    # 9. Concurrency & Kernel Races (16T/80T + 550 renameat2 + 550 hardlink)
    # ------------------------------------------------------------------
    def test_09_concurrency_and_kernel_races(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td).resolve()

            # 16T distinct 0o600 writes
            errs: list[Exception] = []
            def w_distinct(idx: int) -> None:
                try:
                    p = root / f"d_{idx}.bin"
                    audit.write_regular_bytes_safe(p, f"payload-{idx}\n".encode("ascii"))
                    if stat.S_IMODE(os.lstat(p).st_mode) != 0o600:
                        raise RuntimeError("bad mode")
                except Exception as exc:
                    errs.append(exc)

            ts = [threading.Thread(target=w_distinct, args=(i,)) for i in range(16)]
            for t in ts:
                t.start()
            for t in ts:
                t.join()
            self.assertEqual(errs, [])

            # 16T same target + 80T burst
            shared = root / "shared.bin"
            def w_same(idx: int) -> None:
                try:
                    audit.write_regular_bytes_safe(shared, f"v-{idx}\n".encode("ascii"))
                except Exception as exc:
                    errs.append(exc)

            for count in (16, 80):
                ts = [threading.Thread(target=w_same, args=(i,)) for i in range(count)]
                for t in ts:
                    t.start()
                for t in ts:
                    t.join()
            self.assertEqual(errs, [])
            self.assertEqual([p.name for p in root.iterdir() if p.name.endswith(".tmp")], [])

            # 550 renameat2(RENAME_EXCHANGE) symlink-vs-dir races
            libc = ctypes.CDLL("libc.so.6", use_errno=True)
            sys_renameat2 = 316
            rename_exchange = 2
            at_fdcwd = -100
            victim_dir = root / "victim_outside"
            victim_dir.mkdir()
            real_sub = root / "race_slot_a"
            real_sub.mkdir()
            sym_sub = root / "race_slot_b"
            os.symlink(victim_dir, sym_sub)

            stop_flag = threading.Event()
            b_a = os.fsencode(str(real_sub))
            b_b = os.fsencode(str(sym_sub))

            def flipper() -> None:
                while not stop_flag.is_set():
                    libc.syscall(sys_renameat2, at_fdcwd, b_a, at_fdcwd, b_b, rename_exchange)

            fl = threading.Thread(target=flipper)
            fl.start()
            try:
                for i in range(550):
                    try:
                        audit.write_regular_bytes_safe(real_sub / "x.bin", f"{i}\n".encode("ascii"))
                    except oracle.PointChargesN3Error:
                        pass
            finally:
                stop_flag.set()
                fl.join()
            # ZERO escape writes into victim_dir!
            self.assertEqual(list(victim_dir.iterdir()), [])

            # 550 hardlink races: victim file must never be mutated
            hl_dir = root / "hl_race"
            hl_dir.mkdir()
            hl_target = hl_dir / "target.bin"
            hl_victim = hl_dir / "victim.bin"
            hl_victim.write_bytes(b"ORIGINAL_VICTIM\n")

            stop_hl = threading.Event()
            def hl_attacker() -> None:
                while not stop_hl.is_set():
                    try:
                        os.unlink(hl_target)
                    except OSError:
                        pass
                    try:
                        os.link(hl_victim, hl_target)
                    except OSError:
                        pass

            ha = threading.Thread(target=hl_attacker)
            ha.start()
            try:
                for i in range(550):
                    try:
                        audit.write_regular_bytes_safe(hl_target, f"safe-{i}\n".encode("ascii"))
                    except oracle.PointChargesN3Error:
                        pass
            finally:
                stop_hl.set()
                ha.join()
            self.assertEqual(hl_victim.read_bytes(), b"ORIGINAL_VICTIM\n")


if __name__ == "__main__":
    unittest.main()
