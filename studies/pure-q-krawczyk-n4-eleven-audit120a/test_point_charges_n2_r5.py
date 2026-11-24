"""Assertion-free R5 verification & adversarial test suite for Track S20a.

Covers positive pure-Z/Q Krawczyk + Bareiss Sylvester-inertia proofs,
Law-2b 11 -> 9 saddle-node fold wall & zero-free exclusion, boundary tests
killing every arithmetic/comparison/BoolOp mutant, and all 16 Zone-B classes.
"""
from __future__ import annotations

import copy
from fractions import Fraction
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import point_charges_n2_audit as audit
import point_charges_n2_oracle as oracle
from point_charges_n2_oracle import PointChargeCertError


class EvilDict(dict):
    pass


class EvilInt(int):
    pass


class EvilStr(str):
    pass


class EvilList(list):
    pass


class TestPureZQPrimitives(unittest.TestCase):
    def test_parse_and_format_exact_rational_positive_and_boundaries(self) -> None:
        self.assertEqual(oracle.parse_exact_rational(0), Fraction(0, 1))
        self.assertEqual(oracle.parse_exact_rational(-7), Fraction(-7, 1))
        self.assertEqual(oracle.parse_exact_rational(Fraction(3, 4)), Fraction(3, 4))
        self.assertEqual(oracle.parse_exact_rational("-5/8"), Fraction(-5, 8))
        self.assertEqual(oracle.parse_exact_rational("42"), Fraction(42, 1))
        exact_1024 = (1 << 1024) - 1
        over_1024 = 1 << 1024
        self.assertEqual(oracle.parse_exact_rational(exact_1024), Fraction(exact_1024, 1))
        self.assertEqual(oracle.parse_exact_rational(Fraction(exact_1024, 1)), Fraction(exact_1024, 1))
        self.assertEqual(oracle.parse_exact_rational(Fraction(1, exact_1024)), Fraction(1, exact_1024))
        with self.assertRaises(PointChargeCertError):
            oracle.parse_exact_rational(over_1024)
        with self.assertRaises(PointChargeCertError):
            oracle.parse_exact_rational(Fraction(over_1024, 1))
        with self.assertRaises(PointChargeCertError):
            oracle.parse_exact_rational(Fraction(1, over_1024))
        s_256 = "0" * 255 + "1"
        self.assertEqual(oracle.parse_exact_rational(s_256), Fraction(1, 1))
        with self.assertRaises(PointChargeCertError):
            oracle.parse_exact_rational("0" * 256 + "1")
        self.assertEqual(oracle.format_rational(Fraction(5, 1)), "5")
        self.assertEqual(oracle.format_rational(Fraction(-5, 3)), "-5/3")
        exact_fmt = (1 << oracle.MAX_FORMAT_BITS) - 1
        over_fmt = 1 << oracle.MAX_FORMAT_BITS
        self.assertEqual(oracle.format_rational(Fraction(exact_fmt, 1)), str(exact_fmt))
        self.assertEqual(oracle.format_rational(Fraction(1, exact_fmt)), f"1/{exact_fmt}")
        with self.assertRaises(PointChargeCertError):
            oracle.format_rational(5)  # type: ignore[arg-type]
        with self.assertRaises(PointChargeCertError):
            oracle.format_rational(Fraction(over_fmt, 1))
        with self.assertRaises(PointChargeCertError):
            oracle.format_rational(Fraction(1, over_fmt))

    def test_parse_exact_rational_rejections(self) -> None:
        bad_inputs = [
            True,
            False,
            1.5,
            None,
            EvilInt(1),
            EvilStr("1"),
            "",
            "1/2/3",
            "1/0",
            "1/-2",
            "+1/2",
            "1.5",
            " 1",
            "1 ",
            "1/\u0661",
            "a/2",
            "1/b",
            "-/2",
            "-",
        ]
        for bad in bad_inputs:
            with self.assertRaises(PointChargeCertError):
                oracle.parse_exact_rational(bad)

    def test_dyadic_interval_ops_exact_enclosure(self) -> None:
        for val in (Fraction(0, 1), Fraction(1, 3), Fraction(-1, 3), Fraction(7, 8), Fraction(-7, 8)):
            lo, hi = oracle.frac_to_iv(val)
            self.assertLessEqual(Fraction(lo, oracle.SCALE), val)
            self.assertGreaterEqual(Fraction(hi, oracle.SCALE), val)
            self.assertLessEqual(hi - lo, 1)
        self.assertEqual(oracle.iv_add((10, 20), (-5, 7)), (5, 27))
        self.assertEqual(oracle.iv_sub((10, 20), (-5, 7)), (3, 25))
        S = oracle.SCALE
        self.assertEqual(oracle.iv_mul((2 * S, 2 * S), (-3 * S, -S)), (-6 * S, -2 * S))
        self.assertEqual(oracle.iv_mul((S, S), (-S, 2 * S)), (-S, 2 * S))
        self.assertEqual(oracle.iv_sq((-2 * S, 3 * S)), (0, 9 * S))
        self.assertEqual(oracle.iv_sq((-3 * S, 2 * S)), (0, 9 * S))
        self.assertEqual(oracle.iv_sq((0, 3 * S)), (0, 9 * S))
        self.assertEqual(oracle.iv_sq((-3 * S, 0)), (0, 9 * S))
        self.assertEqual(oracle.iv_sq((2 * S, 3 * S)), (4 * S, 9 * S))
        self.assertEqual(oracle.iv_sq((-3 * S, -2 * S)), (4 * S, 9 * S))
        self.assertEqual(oracle.iv_sq((2 * S, 2 * S)), (4 * S, 4 * S))
        self.assertEqual(oracle.iv_abs_max((-10, 3)), 10)
        self.assertEqual(oracle.iv_abs_max((-3, 10)), 10)
        inv3_1, inv5_1 = oracle.iv_inv_r3_r5((S, S))
        self.assertEqual(inv3_1, (S, S))
        self.assertEqual(inv5_1, (S, S))
        inv3_2, inv5_2 = oracle.iv_inv_r3_r5((2 * S, 2 * S))
        self.assertLessEqual(Fraction(inv3_2[0], S) ** 2, Fraction(1, 8))
        self.assertGreaterEqual(Fraction(inv3_2[1], S) ** 2, Fraction(1, 8))
        self.assertLessEqual(Fraction(inv5_2[0], S) ** 2, Fraction(1, 32))
        self.assertGreaterEqual(Fraction(inv5_2[1], S) ** 2, Fraction(1, 32))
        with self.assertRaises(PointChargeCertError) as cm0:
            oracle.iv_inv_r3_r5((0, S))
        self.assertEqual(cm0.exception.code, "SITE_OVERLAP")
        with self.assertRaises(PointChargeCertError) as cm1:
            oracle.iv_inv_r3_r5((-1, S))
        self.assertEqual(cm1.exception.code, "SITE_OVERLAP")
        h_diag = [[(4 * S, 4 * S), (0, 0), (0, 0)], [(0, 0), (-2 * S, -2 * S), (0, 0)], [(0, 0), (0, 0), (-2 * S, -2 * S)]]
        self.assertEqual(oracle.iv_det3(h_diag), (16 * S, 16 * S))

    def test_bareiss_det_mat3_inv_and_sylvester_inertia(self) -> None:
        F = Fraction
        self.assertEqual(oracle.bareiss_det_exact([[F(7, 3)]]), F(7, 3))
        mat_swap = [[F(0), F(1), F(2)], [F(1, 2), F(3), F(4)], [F(5), F(6), F(7, 2)]]
        self.assertEqual(oracle.bareiss_det_exact(mat_swap), F(-23, 4))
        inv_m, det_m = oracle.mat3_inv_exact(mat_swap)
        self.assertEqual(det_m, F(-23, 4))
        for r in range(3):
            for c in range(3):
                self.assertEqual(sum(mat_swap[r][k] * inv_m[k][c] for k in range(3)), F(1 if r == c else 0))
        sing = [[F(1), F(2), F(3)], [F(2), F(4), F(6)], [F(0), F(1), F(1)]]
        self.assertEqual(oracle.bareiss_det_exact(sing), F(0))
        with self.assertRaises(PointChargeCertError) as cm_s:
            oracle.mat3_inv_exact(sing)
        self.assertEqual(cm_s.exception.code, "SINGULAR_MATRIX")
        with self.assertRaises(PointChargeCertError):
            oracle.bareiss_det_exact([])
        with self.assertRaises(PointChargeCertError):
            oracle.bareiss_det_exact([[F(1), F(2)]])
        self.assertEqual(oracle.sylvester_inertia_exact([[F(4), F(0), F(0)], [F(0), F(-2), F(0)], [F(0), F(0), F(-2)]])[:3], (1, 2, 0))
        self.assertEqual(oracle.sylvester_inertia_exact([[F(0), F(1), F(0)], [F(1), F(0), F(0)], [F(0), F(0), F(5)]])[:3], (2, 1, 0))
        self.assertEqual(oracle.sylvester_inertia_exact([[F(0), F(0), F(0)], [F(0), F(0), F(1)], [F(0), F(1), F(0)]])[:3], (1, 1, 1))
        self.assertEqual(oracle.sylvester_inertia_exact([[F(0)] * 3 for _ in range(3)])[:3], (0, 0, 3))
        # Shape & symmetry negative controls (2x2, 3x2, asymmetric)
        with self.assertRaises(PointChargeCertError) as cm_sh1:
            oracle.sylvester_inertia_exact([[F(1), F(0)], [F(0), F(1)]])
        self.assertEqual(cm_sh1.exception.code, "INVALID_MATRIX")
        with self.assertRaises(PointChargeCertError) as cm_sh2:
            oracle.sylvester_inertia_exact([[F(1), F(0)], [F(0), F(1)], [F(0), F(0)]])
        self.assertEqual(cm_sh2.exception.code, "INVALID_MATRIX")
        with self.assertRaises(PointChargeCertError) as cm_asym:
            oracle.sylvester_inertia_exact([[F(1), F(2), F(0)], [F(0), F(1), F(0)], [F(0), F(0), F(1)]])
        self.assertEqual(cm_asym.exception.code, "ASYMMETRIC_MATRIX")
        # Direct verify_congruence_exact negative controls
        eye = [[F(1 if r == c else 0) for c in range(3)] for r in range(3)]
        zero_mat = [[F(0) for _ in range(3)] for _ in range(3)]
        with self.assertRaises(PointChargeCertError) as cm_cs:
            oracle.verify_congruence_exact(eye, zero_mat, [F(1), F(1), F(1)])
        self.assertEqual(cm_cs.exception.code, "CONGRUENCE_SINGULAR")
        with self.assertRaises(PointChargeCertError) as cm_cm:
            oracle.verify_congruence_exact(eye, eye, [F(2), F(1), F(1)])
        self.assertEqual(cm_cm.exception.code, "CONGRUENCE_MISMATCH")


class TestKrawczykAndBifurcationControls(unittest.TestCase):
    def test_n2_analytic_control_and_negative_krawczyk_controls(self) -> None:
        F = Fraction
        sites = [[F(-1), F(0), F(0)], [F(1), F(0), F(0)]]
        charges = [F(1), F(1)]
        hw = F(1, 1 << 20)
        good_box = ((-hw, hw), (-hw, hw), (-hw, hw))
        res = oracle.certify_equilibrium_box(sites, charges, good_box)
        self.assertEqual(res["morse_index"], 2)
        self.assertEqual(res["inertia"], [1, 2, 0])
        self.assertEqual(res["det_sign"], 1)
        self.assertTrue(res["harmonic_trace_contains_zero"])
        shifted = ((F(1, 10), F(1, 10) + hw), (-hw, hw), (-hw, hw))
        with self.assertRaises(PointChargeCertError) as cm:
            oracle.certify_equilibrium_box(sites, charges, shifted)
        self.assertEqual(cm.exception.code, "KRAWCZYK_NOT_INTERIOR")
        huge = ((F(-1, 2), F(1, 2)), (F(-1, 2), F(1, 2)), (F(-1, 2), F(1, 2)))
        with self.assertRaises(PointChargeCertError) as cm2:
            oracle.certify_equilibrium_box(sites, charges, huge)
        self.assertEqual(cm2.exception.code, "KRAWCZYK_NOT_CONTRACTIVE")
        site_box = ((F(9, 10), F(11, 10)), (-hw, hw), (-hw, hw))
        with self.assertRaises(PointChargeCertError) as cm3:
            oracle.certify_equilibrium_box(sites, charges, site_box)
        self.assertEqual(cm3.exception.code, "SITE_OVERLAP")
        # Exact error-code validation on system, charges, and box guards
        with self.assertRaises(PointChargeCertError) as cm_sys1:
            oracle.certify_equilibrium_box([[F(-1), F(0), F(0)]], [F(1)], good_box)
        self.assertEqual(cm_sys1.exception.code, "INVALID_SYSTEM")
        with self.assertRaises(PointChargeCertError) as cm_sys2:
            oracle.certify_equilibrium_box(sites, [F(1)], good_box)
        self.assertEqual(cm_sys2.exception.code, "INVALID_SYSTEM")
        with self.assertRaises(PointChargeCertError) as cm_sys3:
            oracle.certify_equilibrium_box([[F(-1), F(0)], [F(1), F(0)]], charges, good_box)
        self.assertEqual(cm_sys3.exception.code, "INVALID_SYSTEM")
        for bad_q in (F(0), F(-1)):
            with self.assertRaises(PointChargeCertError) as cm_q:
                oracle.certify_equilibrium_box(sites, [F(1), bad_q], good_box)
            self.assertEqual(cm_q.exception.code, "NONPOSITIVE_CHARGE")
        for bad_box in (
            ((F(0), F(0)), (-hw, hw), (-hw, hw)),
            ((F(1), F(0)), (-hw, hw), (-hw, hw)),
            ((-hw, hw), (-hw, hw)),
        ):
            with self.assertRaises(PointChargeCertError) as cm_b:
                oracle.certify_equilibrium_box(sites, charges, bad_box)  # type: ignore[arg-type]
            self.assertEqual(cm_b.exception.code, "INVALID_BOX")
        # Pairwise disjointness & configuration guards
        with self.assertRaises(PointChargeCertError) as cm_ov:
            oracle.verify_pairwise_disjoint([good_box, good_box])
        self.assertEqual(cm_ov.exception.code, "OVERLAPPING_BOXES")
        b_left = ((F(0), F(1)), (F(0), F(1)), (F(0), F(1)))
        b_touch = ((F(1), F(2)), (F(0), F(1)), (F(0), F(1)))
        with self.assertRaises(PointChargeCertError) as cm_t:
            oracle.verify_pairwise_disjoint([b_left, b_touch])
        self.assertEqual(cm_t.exception.code, "OVERLAPPING_BOXES")
        oracle.verify_pairwise_disjoint([b_left, ((F(0), F(1)), (F(0), F(1)), (F(2), F(3)))])
        for bad_hw in (F(0), F(-1)):
            with self.assertRaises(PointChargeCertError) as cm_hw:
                oracle.certify_configuration("bad", sites, charges, [[F(0), F(0), F(0)]], bad_hw, 1)
            self.assertEqual(cm_hw.exception.code, "INVALID_HALF_WIDTH")
        with self.assertRaises(PointChargeCertError) as cm_cnt:
            oracle.certify_configuration("bad", sites, charges, [[F(0), F(0), F(0)]], hw, 2)
        self.assertEqual(cm_cnt.exception.code, "COUNT_MISMATCH")
        with self.assertRaises(PointChargeCertError) as cm_ctr:
            oracle.certify_configuration("bad", sites, charges, [[F(0), F(0)]], hw, 1)
        self.assertEqual(cm_ctr.exception.code, "INVALID_CENTER")
        # Real physical Morse-sum mismatch: only 1 of 11 equilibria on n4-eleven-short -> m2-m1 = 1 != 3
        s_sites, s_charges, s_centers, s_hw = oracle._load_spec_system(oracle.CONFIG_SPECS["n4-eleven-short"])
        with self.assertRaises(PointChargeCertError) as cm_ms:
            oracle.certify_configuration("subset", s_sites, s_charges, [s_centers[0]], s_hw, 1)
        self.assertEqual(cm_ms.exception.code, "MORSE_SUM_MISMATCH")

    def test_hessian_krawczyk_and_bifurcation_helper_guards(self) -> None:
        F = Fraction
        S = oracle.SCALE
        # _verify_hessian_in_box boundaries
        box3 = [[(10, 20) for _ in range(3)] for _ in range(3)]
        oracle._verify_hessian_in_box([[10, 15, 20]] * 3, box3)
        with self.assertRaises(PointChargeCertError) as cm_hlo:
            oracle._verify_hessian_in_box([[9, 15, 15], [15, 15, 15], [15, 15, 15]], box3)
        self.assertEqual(cm_hlo.exception.code, "HESSIAN_MIDPOINT_OUTSIDE")
        with self.assertRaises(PointChargeCertError) as cm_hhi:
            oracle._verify_hessian_in_box([[21, 15, 15], [15, 15, 15], [15, 15, 15]], box3)
        self.assertEqual(cm_hhi.exception.code, "HESSIAN_MIDPOINT_OUTSIDE")
        # _check_krawczyk_contraction exact boundary
        oracle._check_krawczyk_contraction(S - 1)
        with self.assertRaises(PointChargeCertError) as cm_keq:
            oracle._check_krawczyk_contraction(S)
        self.assertEqual(cm_keq.exception.code, "KRAWCZYK_NOT_CONTRACTIVE")
        with self.assertRaises(PointChargeCertError) as cm_kgt:
            oracle._check_krawczyk_contraction(S + 1)
        self.assertEqual(cm_kgt.exception.code, "KRAWCZYK_NOT_CONTRACTIVE")
        # _verify_harmonic_hessian positive & all negative branches
        h_idx2 = [[F(4), F(0), F(0)], [F(0), F(-2), F(0)], [F(0), F(0), F(-2)]]
        hb_idx2 = [[(4 * S, 4 * S), (0, 0), (0, 0)], [(0, 0), (-2 * S, -2 * S), (0, 0)], [(0, 0), (0, 0), (-2 * S, -2 * S)]]
        self.assertEqual(oracle._verify_harmonic_hessian(h_idx2, F(16), hb_idx2), (1, 2, 0, 1, 1))
        h_idx1 = [[F(-4), F(0), F(0)], [F(0), F(2), F(0)], [F(0), F(0), F(2)]]
        hb_idx1 = [[(-4 * S, -4 * S), (0, 0), (0, 0)], [(0, 0), (2 * S, 2 * S), (0, 0)], [(0, 0), (0, 0), (2 * S, 2 * S)]]
        self.assertEqual(oracle._verify_harmonic_hessian(h_idx1, F(-16), hb_idx1), (2, 1, 0, -1, -1))
        # Earnshaw violation: positive definite (+I) and singular diag(1,-1,0)
        eye3 = [[F(1 if r == c else 0) for c in range(3)] for r in range(3)]
        sing3 = [[F(1), F(0), F(0)], [F(0), F(-1), F(0)], [F(0), F(0), F(0)]]
        for bad_h in (eye3, sing3):
            with self.assertRaises(PointChargeCertError) as cm_e:
                oracle._verify_harmonic_hessian(bad_h, F(1), hb_idx2)
            self.assertEqual(cm_e.exception.code, "EARNSHAW_VIOLATION")
        # det_box encloses zero: [0, 4], [-4, 0], [-4, 4]
        for z_entry in ((0, S), (-S, 0), (-S, S)):
            hb_zero = [[z_entry, (0, 0), (0, 0)], [(0, 0), (-2 * S, -2 * S), (0, 0)], [(0, 0), (0, 0), (-2 * S, -2 * S)]]
            with self.assertRaises(PointChargeCertError) as cm_dz:
                oracle._verify_harmonic_hessian(h_idx2, F(16), hb_zero)
            self.assertEqual(cm_dz.exception.code, "DET_INTERVAL_CONTAINS_ZERO")
        # Inertia-det mismatch: det_tilde wrong sign alone, det_box wrong sign alone
        with self.assertRaises(PointChargeCertError) as cm_dm1:
            oracle._verify_harmonic_hessian(h_idx2, F(-16), hb_idx2)
        self.assertEqual(cm_dm1.exception.code, "INERTIA_DET_MISMATCH")
        with self.assertRaises(PointChargeCertError) as cm_dm2:
            oracle._verify_harmonic_hessian(h_idx2, F(16), hb_idx1)
        self.assertEqual(cm_dm2.exception.code, "INERTIA_DET_MISMATCH")
        # Harmonic trace mismatch: tr > 0 and tr < 0
        hb_tr_pos = [[(5 * S, 5 * S), (0, 0), (0, 0)], [(0, 0), (-2 * S, -2 * S), (0, 0)], [(0, 0), (0, 0), (-2 * S, -2 * S)]]
        hb_tr_neg = [[(3 * S, 3 * S), (0, 0), (0, 0)], [(0, 0), (-2 * S, -2 * S), (0, 0)], [(0, 0), (0, 0), (-2 * S, -2 * S)]]
        for bad_hb in (hb_tr_pos, hb_tr_neg):
            with self.assertRaises(PointChargeCertError) as cm_tr:
                oracle._verify_harmonic_hessian(h_idx2, F(16), bad_hb)
            self.assertEqual(cm_tr.exception.code, "HARMONIC_TRACE_MISMATCH")
        # _verify_bifurcation_checks negative paths
        b_ann = ((F(0), F(10)), (F(0), F(10)), (F(0), F(10)))
        b_in1 = ((F(1), F(2)), (F(1), F(2)), (F(1), F(2)))
        b_in2 = ((F(3), F(4)), (F(3), F(4)), (F(3), F(4)))
        b_out = ((F(20), F(21)), (F(20), F(21)), (F(20), F(21)))
        self.assertEqual(oracle._verify_bifurcation_checks([b_in1, b_in2], [b_out], [1, 2], 0, 1, b_ann, {"excluded": True}), (1, 2))
        with self.assertRaises(PointChargeCertError) as cm_fi:
            oracle._verify_bifurcation_checks([b_in1, b_in2], [b_out], [2, 2], 0, 1, b_ann, {"excluded": True})
        self.assertEqual(cm_fi.exception.code, "FOLD_INDEX_MISMATCH")
        with self.assertRaises(PointChargeCertError) as cm_ab1:
            oracle._verify_bifurcation_checks([b_out, b_in2], [b_out], [1, 2], 0, 1, b_ann, {"excluded": True})
        self.assertEqual(cm_ab1.exception.code, "ANNIHILATION_BOX_MISMATCH")
        with self.assertRaises(PointChargeCertError) as cm_ab2:
            oracle._verify_bifurcation_checks([b_in1, b_out], [b_out], [1, 2], 0, 1, b_ann, {"excluded": True})
        self.assertEqual(cm_ab2.exception.code, "ANNIHILATION_BOX_MISMATCH")
        with self.assertRaises(PointChargeCertError) as cm_so:
            oracle._verify_bifurcation_checks([b_in1, b_in2], [b_in1], [1, 2], 0, 1, b_ann, {"excluded": True})
        self.assertEqual(cm_so.exception.code, "SURVIVING_BOX_OVERLAP")
        with self.assertRaises(PointChargeCertError) as cm_ex:
            oracle._verify_bifurcation_checks([b_in1, b_in2], [b_out], [1, 2], 0, 1, b_ann, {"excluded": False})
        self.assertEqual(cm_ex.exception.code, "EXCLUSION_FAILED")
        # _check_expected_morse_counts negative paths
        sum_ok = {"index_counts": {"m1": 4, "m2": 7}}
        oracle._check_expected_morse_counts(sum_ok, 4, 7)
        with self.assertRaises(PointChargeCertError) as cm_m1:
            oracle._check_expected_morse_counts(sum_ok, 3, 7)
        self.assertEqual(cm_m1.exception.code, "EXPECTED_INERTIA_MISMATCH")
        with self.assertRaises(PointChargeCertError) as cm_m2:
            oracle._check_expected_morse_counts(sum_ok, 4, 6)
        self.assertEqual(cm_m2.exception.code, "EXPECTED_INERTIA_MISMATCH")

    def test_zero_free_exclusion_positive_and_negative_controls(self) -> None:
        F = Fraction
        sites = [[F(-1), F(0), F(0)], [F(1), F(0), F(0)]]
        charges = [F(1), F(1)]
        far_box = ((F(1, 5), F(3, 10)), (F(-1, 20), F(1, 20)), (F(-1, 20), F(1, 20)))
        # Exact boundary max_leaves=1 kills '>' -> '>=' immediately!
        res_ok = oracle.certify_zero_free_box(sites, charges, far_box, [F(1), F(0), F(0)], max_leaves=1)
        self.assertTrue(res_ok["excluded"])
        self.assertEqual(res_ok["leaves"], 1)
        self.assertEqual(res_ok["unresolved"], 0)
        # Box containing root (0,0,0) halts cleanly at leaves + unresolved == max_leaves + 1 (kills 'if False:' cleanly!)
        root_box = ((F(-1, 10), F(1, 10)), (F(-1, 10), F(1, 10)), (F(-1, 10), F(1, 10)))
        res_no = oracle.certify_zero_free_box(sites, charges, root_box, [F(1), F(0), F(0)], max_leaves=4)
        self.assertFalse(res_no["excluded"])
        self.assertEqual(res_no["leaves"] + res_no["unresolved"], 5)
        with self.assertRaises(PointChargeCertError):
            oracle.certify_zero_free_box(sites, charges, far_box, [F(0), F(0), F(0)], max_leaves=8)
        with self.assertRaises(PointChargeCertError):
            oracle.certify_zero_free_box(sites, charges, far_box, [F(1), F(0)], max_leaves=8)
        with self.assertRaises(PointChargeCertError):
            oracle.certify_zero_free_box(sites, charges, far_box, [F(1), F(0), F(0)], max_leaves=0)
        with self.assertRaises(PointChargeCertError):
            oracle.certify_zero_free_box(sites, charges, far_box, [F(1), F(0), F(0)], max_leaves=True)  # type: ignore[arg-type]

    def test_box_containment_and_disjoint_boundaries(self) -> None:
        F = Fraction
        outer = ((F(0), F(10)), (F(0), F(10)), (F(0), F(10)))
        inner = ((F(1), F(9)), (F(1), F(9)), (F(1), F(9)))
        self.assertTrue(oracle._box_strictly_inside(inner, outer))
        self.assertFalse(oracle._box_strictly_inside(((F(0), F(9)), (F(1), F(9)), (F(1), F(9))), outer))
        self.assertFalse(oracle._box_strictly_inside(((F(1), F(10)), (F(1), F(9)), (F(1), F(9))), outer))
        self.assertTrue(oracle._boxes_disjoint(outer, ((F(11), F(12)), (F(0), F(10)), (F(0), F(10)))))
        self.assertFalse(oracle._boxes_disjoint(outer, ((F(10), F(12)), (F(0), F(10)), (F(0), F(10)))))


class TestFullCertificationAndZoneBHardening(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.cert = audit.build_certificate()

    def test_full_certificate_scientific_invariants(self) -> None:
        payload = self.cert["payload"]
        self.assertEqual(payload["track"], "TRACK-S20A-POINT-CHARGES-N2-EQUILIBRIA")
        self.assertTrue(payload["arithmetic"]["zero_floats"])
        cfgs = payload["configurations"]
        self.assertEqual(cfgs["n2-analytic-control"]["certified_equilibria_lower_bound"], 1)
        self.assertEqual(cfgs["n2-analytic-control"]["morse_alternating_sum"], 1)
        for name in ("n4-eleven-short", "n4-eleven-flat", "n4-eleven-jet"):
            c = cfgs[name]
            self.assertEqual(c["n_charges"], 4)
            self.assertEqual(c["maxwell_bound"], 9)
            self.assertEqual(c["certified_equilibria_lower_bound"], 11)
            self.assertTrue(c["exceeds_maxwell_bound"])
            self.assertTrue(c["pairwise_disjoint"])
            self.assertTrue(c["all_nondegenerate"])
            self.assertEqual(c["index_counts"], {"m0": 0, "m1": 4, "m2": 7, "m3": 0})
            self.assertEqual(c["morse_alternating_sum"], 3)
            self.assertLess(oracle.parse_exact_rational(c["max_krawczyk_rho"]), Fraction(1, 1000))
        bif = payload["bifurcation_wall_11_to_9"]
        self.assertEqual(bif["bracket_width"], "1/1000000000")
        self.assertEqual(bif["pre_fold"]["certified_equilibria_lower_bound"], 11)
        self.assertEqual(bif["post_fold"]["certified_equilibria_lower_bound"], 9)
        self.assertEqual(bif["pre_fold"]["index_counts"], {"m0": 0, "m1": 4, "m2": 7, "m3": 0})
        self.assertEqual(bif["post_fold"]["index_counts"], {"m0": 0, "m1": 3, "m2": 6, "m3": 0})
        self.assertEqual(sorted(bif["colliding_morse_indices"]), [1, 2])
        self.assertTrue(bif["post_fold_zero_free_cover"]["excluded"])
        self.assertTrue(bif["poincare_hopf_invariant_across_fold"])
        ver = audit.verify_certificate_dict(self.cert)
        self.assertTrue(ver["valid"])

    def test_zone_b_in_memory_adversarial_vectors(self) -> None:
        # ZB01: non-ASCII digest & hex validation
        c_zb01 = copy.deepcopy(self.cert)
        c_zb01["payload_sha256"] = "\u00e9" * 64
        with self.assertRaises(PointChargeCertError) as cm_z1:
            audit.verify_certificate_dict(c_zb01)
        self.assertEqual(cm_z1.exception.code, "INVALID_SHA256")
        for bad_hex in ("A" * 64, "a" * 63, "a" * 65, "g" * 64, 123):
            c_h = copy.deepcopy(self.cert)
            c_h["certificate_sha256"] = bad_hex
            with self.assertRaises(PointChargeCertError):
                audit.verify_certificate_dict(c_h)
        # ZB02: lone UTF-16 surrogates
        c_sur = copy.deepcopy(self.cert)
        c_sur["certificate_sha256"] = "\ud800" * 64
        with self.assertRaises(PointChargeCertError):
            audit.verify_certificate_dict(c_sur)
        with self.assertRaises(PointChargeCertError):
            audit.validate_json_tree_iterative({"\ud800": 1})
        with self.assertRaises(PointChargeCertError):
            audit.validate_json_tree_iterative({"k": "\ud800"})
        # ZB03: huge int & MAX_INT_BITS exact boundary
        audit.validate_json_tree_iterative({"k": (1 << audit.MAX_INT_BITS) - 1})
        with self.assertRaises(PointChargeCertError):
            audit.validate_json_tree_iterative({"k": 1 << audit.MAX_INT_BITS})
        with self.assertRaises(PointChargeCertError):
            audit.validate_json_tree_iterative({"k": 10 ** 5000})
        # ZB04: exact depth & exact node count boundaries
        ok_nest: object = 0
        for _ in range(audit.MAX_JSON_DEPTH):
            ok_nest = [ok_nest]
        audit.validate_json_tree_iterative(ok_nest)
        with self.assertRaises(PointChargeCertError) as cm_dp:
            audit.validate_json_tree_iterative([ok_nest])
        self.assertEqual(cm_dp.exception.code, "JSON_TOO_DEEP")
        audit.validate_json_tree_iterative([0] * (audit.MAX_JSON_NODES - 1))
        with self.assertRaises(PointChargeCertError) as cm_nd:
            audit.validate_json_tree_iterative([0] * audit.MAX_JSON_NODES)
        self.assertEqual(cm_nd.exception.code, "JSON_TOO_LARGE")
        # ZB05 / ZB06: bool-as-int re-signed SHA forgery
        c_bool = copy.deepcopy(self.cert)
        c_bool["payload"]["configurations"]["n2-analytic-control"]["certified_equilibria_lower_bound"] = True
        p_sha = audit.sha256_hex(audit.canonical_json_bytes(c_bool["payload"]))
        c_bool["payload_sha256"] = p_sha
        c_bool["certificate_sha256"] = audit.sha256_hex(
            audit.canonical_json_bytes({"schema": c_bool["schema"], "payload": c_bool["payload"], "payload_sha256": p_sha})
        )
        with self.assertRaises(PointChargeCertError) as cm_bl:
            audit.verify_certificate_dict(c_bool)
        self.assertEqual(cm_bl.exception.code, "ORACLE_PAYLOAD_MISMATCH")
        # Direct exact_json_equal checks
        self.assertFalse(audit.exact_json_equal(1, True))
        self.assertFalse(audit.exact_json_equal({"a": 1}, {"a": 1, "b": 2}))
        self.assertFalse(audit.exact_json_equal({"a": 1}, {"b": 1}))
        self.assertFalse(audit.exact_json_equal([1], [1, 2]))
        self.assertFalse(audit.exact_json_equal([1], [2]))
        self.assertTrue(audit.exact_json_equal({"a": [1, None, True]}, {"a": [1, None, True]}))
        # ZB08 / ZB09 / ZB12
        with self.assertRaises(PointChargeCertError):
            audit.validate_json_tree_iterative({"x": 1.5})
        for bad_root in ([], None, "x", EvilDict(self.cert)):
            with self.assertRaises(PointChargeCertError):
                audit.verify_certificate_dict(bad_root)
        for bad_leaf in (EvilInt(1), EvilStr("a"), EvilList([1]), (1, 2)):
            with self.assertRaises(PointChargeCertError):
                audit.validate_json_tree_iterative({"x": bad_leaf})
        with self.assertRaises(PointChargeCertError):
            audit.validate_json_tree_iterative({1: "a"})
        # ZB13 + distinct error codes for schema, keys, payload type, payload sha, cert sha
        c_hol = copy.deepcopy(self.cert)
        c_hol["payload"] = {}
        h_sha = audit.sha256_hex(audit.canonical_json_bytes(c_hol["payload"]))
        c_hol["payload_sha256"] = h_sha
        c_hol["certificate_sha256"] = audit.sha256_hex(
            audit.canonical_json_bytes({"schema": c_hol["schema"], "payload": c_hol["payload"], "payload_sha256": h_sha})
        )
        with self.assertRaises(PointChargeCertError) as cm_hol:
            audit.verify_certificate_dict(c_hol)
        self.assertEqual(cm_hol.exception.code, "ORACLE_PAYLOAD_MISMATCH")
        c_bad_schema = copy.deepcopy(self.cert)
        c_bad_schema["schema"] = "wrong"
        c_bad_schema["certificate_sha256"] = audit.sha256_hex(
            audit.canonical_json_bytes({"schema": "wrong", "payload": c_bad_schema["payload"], "payload_sha256": c_bad_schema["payload_sha256"]})
        )
        with self.assertRaises(PointChargeCertError) as cm_sch:
            audit.verify_certificate_dict(c_bad_schema)
        self.assertEqual(cm_sch.exception.code, "SCHEMA_MISMATCH")
        c_extra_key = copy.deepcopy(self.cert)
        c_extra_key["extra"] = 1
        with self.assertRaises(PointChargeCertError) as cm_ek:
            audit.verify_certificate_dict(c_extra_key)
        self.assertEqual(cm_ek.exception.code, "CERT_KEYS_MISMATCH")
        c_list_pay = copy.deepcopy(self.cert)
        c_list_pay["payload"] = []
        with self.assertRaises(PointChargeCertError) as cm_lp:
            audit.verify_certificate_dict(c_list_pay)
        self.assertEqual(cm_lp.exception.code, "PAYLOAD_NOT_DICT")
        c_bad_psha = copy.deepcopy(self.cert)
        c_bad_psha["payload_sha256"] = "0" * 64
        c_bad_psha["certificate_sha256"] = audit.sha256_hex(
            audit.canonical_json_bytes({"schema": c_bad_psha["schema"], "payload": c_bad_psha["payload"], "payload_sha256": "0" * 64})
        )
        with self.assertRaises(PointChargeCertError) as cm_ps:
            audit.verify_certificate_dict(c_bad_psha)
        self.assertEqual(cm_ps.exception.code, "PAYLOAD_SHA_MISMATCH")
        c_bad_csha = copy.deepcopy(self.cert)
        c_bad_csha["certificate_sha256"] = "0" * 64
        with self.assertRaises(PointChargeCertError) as cm_cs:
            audit.verify_certificate_dict(c_bad_csha)
        self.assertEqual(cm_cs.exception.code, "CERT_SHA_MISMATCH")

    def test_zone_b_parser_file_and_cli_controls(self) -> None:
        for empty_or_bad in ("", "   ", 123):
            with self.assertRaises(PointChargeCertError) as cm_em:
                audit.parse_strict_json_text(empty_or_bad)  # type: ignore[arg-type]
            self.assertEqual(cm_em.exception.code, "EMPTY_JSON")
        for bad_json in (
            "[]",
            "null",
            '{"a": 1.5}',
            '{"a": NaN}',
            '{"a": Infinity}',
            '{"a": -Infinity}',
            '{"a": 1, "a": 2}',
            '{"a": 1, "\\u0061": 2}',
            '{"a": ' + "9" * 5000 + "}",
            "[" * 1200 + "]" * 1200,
            '{"a":',
            '{"a": "\\ud800"}',
            "\ud800",
        ):
            with self.assertRaises(PointChargeCertError):
                audit.parse_strict_json_text(bad_json)

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            cert_path = tdir / "cert.json"
            out_w = io.StringIO()
            with patch("sys.stdout", out_w):
                rc_w = audit.main(["--write", str(cert_path)])
            self.assertEqual(rc_w, 0)
            out_c = io.StringIO()
            with patch("sys.stdout", out_c):
                rc_c = audit.main(["--check", str(cert_path)])
            self.assertEqual(rc_c, 0)
            self.assertTrue(json.loads(out_c.getvalue())["valid"])
            out_s = io.StringIO()
            with patch("sys.stdout", out_s):
                rc_s = audit.main(["--standalone", "--json"])
            self.assertEqual(rc_s, 0)
            self.assertEqual(out_s.getvalue().encode("ascii"), cert_path.read_bytes())
            out_bare = io.StringIO()
            with patch("sys.stdout", out_bare):
                rc_bare = audit.main([])
            self.assertEqual(rc_bare, 0)
            self.assertEqual(out_bare.getvalue().encode("ascii"), cert_path.read_bytes())
            # ZB10: CRLF (explicit code check!) and pretty-printed non-canonical JSON
            crlf_path = tdir / "crlf.json"
            crlf_path.write_bytes(cert_path.read_bytes().replace(b"\n", b"\r\n"))
            with self.assertRaises(PointChargeCertError) as cm_cr:
                audit.check_certificate_file(crlf_path)
            self.assertEqual(cm_cr.exception.code, "CRLF_FORBIDDEN")
            pretty_path = tdir / "pretty.json"
            pretty_path.write_text(json.dumps(self.cert, indent=2) + "\n", encoding="utf-8")
            with self.assertRaises(PointChargeCertError) as cm_pr:
                audit.check_certificate_file(pretty_path)
            self.assertEqual(cm_pr.exception.code, "NON_CANONICAL_JSON")
            # ZB14: symlink file and symlink ancestor
            sym_file = tdir / "sym.json"
            os.symlink(cert_path, sym_file)
            with self.assertRaises(PointChargeCertError) as cm_sf:
                audit.check_certificate_file(sym_file)
            self.assertEqual(cm_sf.exception.code, "SYMLINK_FORBIDDEN")
            real_sub = tdir / "real_sub"
            real_sub.mkdir()
            (real_sub / "cert.json").write_bytes(cert_path.read_bytes())
            sym_dir = tdir / "sym_sub"
            os.symlink(real_sub, sym_dir)
            with self.assertRaises(PointChargeCertError) as cm_sd:
                audit.check_certificate_file(sym_dir / "cert.json")
            self.assertEqual(cm_sd.exception.code, "SYMLINK_FORBIDDEN")
            # ZB15: missing file and directory path (explicit code checks!)
            with self.assertRaises(PointChargeCertError) as cm_nf:
                audit.check_certificate_file(tdir / "missing.json")
            self.assertEqual(cm_nf.exception.code, "FILE_NOT_FOUND")
            with self.assertRaises(PointChargeCertError) as cm_dr:
                audit.check_certificate_file(tdir)
            self.assertEqual(cm_dr.exception.code, "NOT_A_FILE")
            # ZB16: invalid UTF-8
            bad_utf8 = tdir / "bad.json"
            bad_utf8.write_bytes(b"\xff\xfe\x00\x80")
            with self.assertRaises(PointChargeCertError) as cm_u8:
                audit.check_certificate_file(bad_utf8)
            self.assertEqual(cm_u8.exception.code, "INVALID_UTF8")
            # CLI invalid args (covers --check/--write arity and flag mismatches)
            for bad_argv in (["--bogus"], ["--check"], ["--write"], ["--foo", "bar"]):
                err_buf = io.StringIO()
                with patch("sys.stderr", err_buf):
                    rc_bad = audit.main(bad_argv)
                self.assertEqual(rc_bad, 2)
                self.assertIn("INVALID_CLI_ARGS", err_buf.getvalue())
            err_buf2 = io.StringIO()
            with patch("sys.stderr", err_buf2):
                rc_bad2 = audit.main(["--check", str(tdir / "missing.json")])
            self.assertEqual(rc_bad2, 2)
            self.assertIn("FILE_NOT_FOUND", err_buf2.getvalue())


            # R2 S20A-BUG-02 regressions: string path acceptance, invalid path types, --write dir/OSError/empty
            res_str = audit.check_certificate_file(str(cert_path))
            self.assertTrue(res_str["valid"])
            for bad_path in ("", None, 123, EvilStr(str(cert_path))):
                with self.assertRaises(PointChargeCertError) as cm_bp:
                    audit.check_certificate_file(bad_path)  # type: ignore[arg-type]
                self.assertEqual(cm_bp.exception.code, "INVALID_PATH")
            err_w_dir = io.StringIO()
            with patch("sys.stderr", err_w_dir):
                rc_w_dir = audit.main(["--write", str(tdir)])
            self.assertEqual(rc_w_dir, 2)
            self.assertIn("NOT_A_FILE", err_w_dir.getvalue())
            err_w_os = io.StringIO()
            with patch("sys.stderr", err_w_os):
                rc_w_os = audit.main(["--write", str(tdir / "missing_parent" / "cert.json")])
            self.assertEqual(rc_w_os, 2)
            self.assertIn("FILE_WRITE_ERROR", err_w_os.getvalue())
            err_w_em = io.StringIO()
            with patch("sys.stderr", err_w_em):
                rc_w_em = audit.main(["--write", ""])
            self.assertEqual(rc_w_em, 2)
            self.assertIn("INVALID_PATH", err_w_em.getvalue())

    def test_r2_public_api_zone_b_type_hardening(self) -> None:
        # S20A-BUG-01: bareiss_det_exact, mat3_inv_exact, verify_congruence_exact, sylvester_inertia_exact
        f0, f1 = Fraction(0, 1), Fraction(1, 1)
        id3 = [[f1, f0, f0], [f0, f1, f0], [f0, f0, f1]]
        for bad_mat in (None, [], EvilList([[f1]]), [[True]], [[1.0]], [[1]], [None], [[f1, f0]]):
            with self.assertRaises(PointChargeCertError) as cm_bm:
                oracle.bareiss_det_exact(bad_mat)  # type: ignore[arg-type]
            self.assertEqual(cm_bm.exception.code, "INVALID_MATRIX")
        for bad_3x3 in (
            None,
            [[f1]],
            [[True, False, False], [False, True, False], [False, False, True]],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            EvilList(id3),
        ):
            with self.assertRaises(PointChargeCertError) as cm_sy:
                oracle.sylvester_inertia_exact(bad_3x3)  # type: ignore[arg-type]
            self.assertEqual(cm_sy.exception.code, "INVALID_MATRIX")
            with self.assertRaises(PointChargeCertError) as cm_inv:
                oracle.mat3_inv_exact(bad_3x3)  # type: ignore[arg-type]
            self.assertEqual(cm_inv.exception.code, "INVALID_MATRIX")
        for bad_diag in (None, [f1, f1], [True, f1, f1], EvilList([f1, f1, f1])):
            with self.assertRaises(PointChargeCertError) as cm_dg:
                oracle.verify_congruence_exact(id3, id3, bad_diag)  # type: ignore[arg-type]
            self.assertEqual(cm_dg.exception.code, "INVALID_MATRIX")
        # certify_equilibrium_box, certify_zero_free_box, verify_pairwise_disjoint, certify_configuration
        valid_sites = [[Fraction(-1, 1), f0, f0], [f1, f0, f0]]
        valid_charges = [f1, f1]
        hw = Fraction(1, 1024)
        valid_box = ((-hw, hw), (-hw, hw), (-hw, hw))
        for bad_sites, bad_charges in (
            (None, valid_charges),
            (valid_sites, None),
            (EvilList(valid_sites), valid_charges),
            (valid_sites, EvilList(valid_charges)),
            ([[True, False, False], [f1, f0, f0]], valid_charges),
            (valid_sites, [True, True]),
            ([[1.0, 0.0, 0.0], [f1, f0, f0]], valid_charges),
            (valid_sites, [1.0, 1.0]),
            ([None, [f1, f0, f0]], valid_charges),
        ):
            with self.assertRaises(PointChargeCertError) as cm_sc:
                oracle.certify_equilibrium_box(bad_sites, bad_charges, valid_box)  # type: ignore[arg-type]
            self.assertEqual(cm_sc.exception.code, "INVALID_SYSTEM")
        for bad_box in (
            None,
            EvilList(list(valid_box)),
            ((-hw, hw), (-hw, hw)),
            ((hw,), (-hw, hw), (-hw, hw)),
            (None, (-hw, hw), (-hw, hw)),
            ((False, True), (-hw, hw), (-hw, hw)),
            ((-0.1, 0.1), (-hw, hw), (-hw, hw)),
        ):
            with self.assertRaises(PointChargeCertError) as cm_bx:
                oracle.certify_equilibrium_box(valid_sites, valid_charges, bad_box)  # type: ignore[arg-type]
            self.assertEqual(cm_bx.exception.code, "INVALID_BOX")
        with self.assertRaises(PointChargeCertError) as cm_vpd:
            oracle.verify_pairwise_disjoint(None)  # type: ignore[arg-type]
        self.assertEqual(cm_vpd.exception.code, "INVALID_BOX")
        for bad_norm in (None, EvilList([f1, f0, f0]), [f1, f0], [True, False, False], [1.0, 0.0, 0.0], [f0, f0, f0]):
            with self.assertRaises(PointChargeCertError) as cm_fn:
                oracle.certify_zero_free_box(valid_sites, valid_charges, valid_box, bad_norm)  # type: ignore[arg-type]
            self.assertEqual(cm_fn.exception.code, "INVALID_NORMAL")
        with self.assertRaises(PointChargeCertError) as cm_nm:
            oracle.certify_configuration(123, valid_sites, valid_charges, [[f0, f0, f0]], hw, 1)  # type: ignore[arg-type]
        self.assertEqual(cm_nm.exception.code, "INVALID_NAME")
        for bad_hw in (0.01, True, None):
            with self.assertRaises(PointChargeCertError) as cm_hw:
                oracle.certify_configuration("x", valid_sites, valid_charges, [[f0, f0, f0]], bad_hw, 1)  # type: ignore[arg-type]
            self.assertEqual(cm_hw.exception.code, "INVALID_HALF_WIDTH")
        for bad_ec in (True, 0, -1, 1.0, None):
            with self.assertRaises(PointChargeCertError) as cm_ec:
                oracle.certify_configuration("x", valid_sites, valid_charges, [[f0, f0, f0]], hw, bad_ec)  # type: ignore[arg-type]
            self.assertEqual(cm_ec.exception.code, "COUNT_MISMATCH")
        for bad_ctrs in (None, EvilList([[f0, f0, f0]]), [None], [[True, False, False]], [[0.0, 0.0, 0.0]]):
            with self.assertRaises(PointChargeCertError) as cm_ct:
                oracle.certify_configuration("x", valid_sites, valid_charges, bad_ctrs, hw, 1)  # type: ignore[arg-type]
            self.assertEqual(cm_ct.exception.code, "INVALID_CENTER")


if __name__ == "__main__":
    unittest.main()
