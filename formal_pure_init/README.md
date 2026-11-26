# Pure-`Init` Lean `4.36.0` Exact-Rational Certificate Compiler (`formal_pure_init/`)

Self-contained, zero-external-dependency (`pure Init`, `lake-manifest.json` `packages: []`)
Lean `4.36.0` certificate module and deterministic Python 3 standard-library compiler
for exact-rational electrostatic point-charge witnesses (Refs `#6`).

## Overview

Whereas `formal/PointCharges.lean` sketches real-analysis problem statements against Mathlib `v4.18.0`,
`formal_pure_init/` compiles exact cleared-denominator integer and rational witnesses from
`witness_pack.json` directly into a standalone **pure-`Init` Lean `4.36.0`** library
(`PointChargesExactCert.lean`) paired with a separate axiom-audit module (`AxiomAudit.lean`)
that builds in a few seconds with **zero external packages**, **zero warnings**,
**zero `Classical.choice`**, **zero `bv_decide`**, **zero `native_decide`**, **zero `sorryAx`**,
and passes cold kernel replay via `leanchecker --fresh` (exit code `0`).

## Contents

| File | Purpose |
|---|---|
| `witness_pack.json` | Canonical exact-rational witness bundle (equilibria, off-plane confinement, outward-rounded $2^{40}$ integer-grid box collections for $n=4$ and $n=8$, unit-sphere $S^2(\mathbb{Q})$ Thomson energy comparisons for $N=4, 6$, and rational stereographic $S^2(\mathbb{Q})$ / Cayley $S^3(\mathbb{Q})$ sphere bridges) |
| `compile_point_charges_lean_cert.py` | Pure-stdlib (`int`, `fractions.Fraction`, `math.isqrt`, `math.lcm`) exact validator, deterministic Lean `4.36.0` code generator, hardened pure-$\mathbb{Q}$ AST checker, char-literal-safe static Lean trust-escape scanner, cold `lean` + `leanchecker --fresh` runner, and self-bound canonical JSON certificate verifier |
| `PointChargesExactCert.lean` | Deterministically compiled pure-`Init` Lean `4.36.0` module (`0` imports, `22` theorems: `14` universal `∀` + `8` compiled exact-rational witness theorems) |
| `AxiomAudit.lean` | Standalone `#print axioms` audit module importing only `PointChargesExactCert` |
| `lakefile.toml`, `lake-manifest.json`, `lean-toolchain` | Pinned pure-`Init` Lake `4.36.0` package configuration (`packages: []`) |
| `point_charges_lean_cert.json` | Self-bound canonical ASCII JSON certificate pinning raw-bytes SHA-256 digests of all 6 source/configuration files |
| `test_point_charges_lean_cert.py` | Zero-bare-assert, zero-skip unit, AST-purity, char-literal-escape, tamper-battery, and mutation-killing test suite |

## Verified Theorem Inventory (`22 / 22`, Axioms $\subseteq \{\text{propext}, \text{Quot.sound}\}$)

### Universal (`∀`) Theorems (`t01`–`t14`, axioms `[propext, Quot.sound]`)

| Theorem | Exact Statement Verified in Lean `4.36.0` |
|---|---|
| `t01_coulomb_single_hessian_trace_zero` | `∀ w dx dy dz : Int`, single-charge cleared-denominator Coulomb Hessian has trace `0` |
| `t02_earnshaw_total_hessian_trace_zero` | `∀ terms : List CoulombTermZ`, arbitrary finite-sum Coulomb Hessian has trace `0` |
| `t03_earnshaw_no_strict_definite_quadform` | `∀ h : Sym3Z`, `symTrace h = 0` rules out both strict positive and strict negative definiteness on nonzero `Vec3Z` |
| `t04_planar_offplane_force_and_z_curvature` | `∀ z ws`, strictly positive planar weights imply restoring off-plane force sign (`z > 0 → F_z < 0`, `z < 0 → F_z > 0`), `H_zz < 0`, and positive in-plane trace |
| `t05_lagrange_congruence_quadform_identity` | `∀ s h u`, exact $3\times 3$ congruence quadratic-form identity $Q_h(S u) = Q_{S^\top h S}(u)$ |
| `t06_lagrange_diagonal_morse_index2_witness` | `∀ s h d1 d2 d3`, diagonal congruence with `d1 > 0, d2 < 0, d3 < 0` yields `(+, -, -)` quadratic-form signs on `S e1, S e2, S e3` |
| `t07_lagrange_diagonal_morse_index1_witness` | `∀ s h d1 d2 d3`, diagonal congruence with `d1 > 0, d2 > 0, d3 < 0` yields `(+, +, -)` quadratic-form signs on `S e1, S e2, S e3` |
| `t08_disjoint_boxes_no_common_point` | `∀ b1 b2`, `boxDisjoint b1 b2 = true` implies no integer point lies in both boxes |
| `t09_inv_sqrt_bracket_order_preservation` | `∀ ln ld un ud r2n r2d : Nat`, `0 < r2d` and valid rational square-root enclosures preserve `ln * ud ≤ un * ld` |
| `t10_stereographic_s2q_unit_sphere` | `∀ un ud vn vd : Int`, cleared-denominator rational stereographic projection lands on $S^2$ ($a^2+b^2+c^2=d^2$) |
| `t11_stereographic_s2q_chord_distance_sq` | `∀ u1 v1 u2 v2 : Int`, cleared-denominator stereographic chord-squared identity $dx^2+dy^2+dz^2 = 4((u_1-u_2)^2+(v_1-v_2)^2) d_1 d_2$ |
| `t12_cayley_s3q_unit_norm_and_equatorial_bridge` | `∀ w1 w2 w3 : Int`, cleared-denominator Cayley $S^3$ unit-norm identity `(1 - s)^2 + (2*w1)^2 + (2*w2)^2 + (2*w3)^2 = (1 + s)^2` (with `s = w1^2 + w2^2 + w3^2`) and two-variable sign identity `(1 - s2) = -(s2 - 1)` for `s2 = w1^2 + w2^2` (the 4-coordinate equatorial Cayley-to-stereographic bridge `equatorialCayleyMatchesStereo` is verified in `t22`) |
| `t13_n3_isosceles_cubic_cancellation_and_pythagorean_roots` | `∀ h y k : Int`, exact algebraic factorization for the $n=3$ isosceles geometry and its rational root family $h = 3k / (k^2 - 2)$ |
| `t14_regular_tetrahedron_s2q_sqrt3_obstruction` | `∀ a b : Nat`, `0 < b → 3 * (b * b) ≠ a * a` (strong-induction descent showing $\sqrt{3} \notin \mathbb{Q}$, obstructing an exact regular tetrahedron on $S^2(\mathbb{Q})$) |

### Compiled Exact-Rational Witness Theorems (`t15`–`t22`)

| Theorem | Exact Witness Checked in Lean Kernel | Axioms |
|---|---|---|
| `t15_compiled_n2_analytic_equilibrium_and_hessian` | `n = 2` analytic equilibrium: `totalCoulombGrad w_n2_grad = vecZero`, `totalCoulombHessian w_n2_hess = diagSym 4 (-2) (-2)`, `tr = 0`, `det ≠ 0`, quadratic-form signs `(+, -, -)` (Morse index `2`) | none (`∅`) |
| `t16_compiled_n3_weierstrass_incenter_iso_acute` | `n = 3` acute isosceles Weierstrass incenter (cleared integer witness): `totalCoulombGrad w_n3_iso_grad = vecZero`, `totalCoulombHessian w_n3_iso_hess = diagSym 882000000 882000000 (-1764000000)`, `tr = 0`, `det ≠ 0`, quadratic-form signs `(+, +, -)` (Morse index `1`; reduced rational Hessian `diag(384/625, 384/625, -768/625)` is validated in Python and recorded in `point_charges_lean_cert.json`) | none (`∅`) |
| `t17_compiled_n3_weierstrass_incenter_scalene` | `n = 3` scalene Weierstrass incenter (cleared integer witness): `totalCoulombGrad w_n3_sca_grad = vecZero`, `totalCoulombHessian w_n3_sca_hess = diagSym 1703087936222400 1703087936222400 (-3406175872444800)`, `tr = 0`, `det ≠ 0`, quadratic-form signs `(+, +, -)` (Morse index `1`; reduced rational Hessian `diag(30780/48841, 30780/48841, -61560/48841)` is validated in Python and recorded in `point_charges_lean_cert.json`) | none (`∅`) |
| `t18_compiled_n3_pythagorean_off_plane_confinement` | Cleared integer weight-list inequalities `planarForceZ 4 w_offplane_weights < 0` (`-12 < 0`), `planarHessZZ w_offplane_weights < 0`, and `0 < planarInPlaneTrace w_offplane_weights` for `w_offplane_weights = [1, 1, 1]` at `z = 4` (the full rational `r_k = 5`, `V_z = -12/125` Coulomb evaluation is validated by `validate_off_plane_witness` and recorded in `point_charges_lean_cert.json`, not re-encoded inside `t18`) | none (`∅`) |
| `t19_compiled_n4_eleven_disjointness_and_morse` | `n = 4` (`11` outward-rounded $2^{40}$ boxes): `55/55` pairwise disjoint, `44/44` site-avoiding, `m1=4, m2=7`, $m_2 - m_1 = 3$, `11 > 9` | none (`∅`) |
| `t20_compiled_n8_fiftyone_disjointness_and_morse` | `n = 8` (`51` outward-rounded $2^{40}$ boxes): `1,275/1,275` pairwise disjoint, `408/408` site-avoiding, `m1=22, m2=29`, $m_2 - m_1 = 7$, `51 > 49` | none (`∅`) |
| `t21_compiled_thomson_s2q_energy_dominance` | `N=4` (`w_th_n4_lo_*` vs `w_th_n4_hi_*`) and `N=6` (`w_th_n6_lo_*` vs `w_th_n6_hi_*`): unit-sphere (`siteOnUnitSphere`), exact pairwise squared-distance matching (`pairDistsMatch`), valid rational square-root enclosures (`validSqrtBracket` with reduced dyadic bracket denominators dividing $2^{12}$: `{1, 128, 256, 1024, 2048, 4096}` from 12-bit floor/ceil bracketing), strict `sumE2` & `sumE1Hi < sumE1Lo` separation via `qLt`, and `qLt ⟨9, 4⟩ (sumE2 w_th_n4_lo_br)` (reduced summary values `E_2 = 83825/36864, 5/2, 27/4, 2641/288` are recorded in `point_charges_lean_cert.json`) | `[propext]` |
| `t22_compiled_stereographic_s2q_and_cayley_s3q_bridges` | Exact rational stereographic $S^2(\mathbb{Q})$ map, unit-sphere, & chord-squared checks (`w_stereo_p0a/b`, `w_stereo_p1a/b`), Cayley $S^3(\mathbb{Q})$ map & unit-norm checks (`w_cayley_q0`, `w_cayley_q1`), and full 4-coordinate equatorial bridge `equatorialCayleyMatchesStereo w_cayley_q1 w_stereo_p0a` | none (`∅`) |

## Scope Boundaries

1. **Box Geometry vs. Interval Krawczyk Root Existence (`t19`, `t20`):** `PointChargesExactCert.lean`
   formally checks in Lean's kernel that the outward-rounded $2^{40}$ integer enclosures of the `11`
   (`n=4`) and `51` (`n=8`) boxes are strictly pairwise disjoint (`t08` + `t19`/`t20`), strictly avoid
   all charge sites, and satisfy the Poincaré–Hopf Morse-index count ($m_2 - m_1 = n - 1$) exceeding
   Maxwell's $(n-1)^2$ bound. For exact dyadic charge sites whose denominators divide $2^{40}$
   (0-based indices `3`–`6` of `w_n8_sites`: `5/4`, `3/4`, `2049/2048`, `2047/2048`), outward rounding produces
   exact zero-width point intervals (`lo == hi` on all axes) rather than width-1 boxes; non-dyadic
   sites round outward to width-1 intervals (`hi == lo + 1`), and `t08` proves strict closed-interval
   avoidance in both cases. Analytic Krawczyk contraction proving root existence and uniqueness
   inside each box over $\mathbb{R}^3$ is verified by the repository's interval-arithmetic verifiers
   (`verify_arb.py` / `studies/`) and is not re-proved in real analysis inside pure-`Init` Lean.
2. **Regular Tetrahedron on $S^2(\mathbb{Q})$ (`t14`, `t21`):** Because a regular tetrahedron inscribed
   in the unit sphere has irrational volume $8/(9\sqrt{3})$ while any four rational unit-sphere points
   span a rational volume, `t14` proves $\sqrt{3} \notin \mathbb{Q}$ and `t21` certifies strict $E_2$
   and bracketed $E_1$ dominance of an exact $S^2(\mathbb{Q})$ near-tetrahedron over the equatorial
   square (reduced summary values $E_2 = 83825/36864$ vs. $5/2$ recorded in `point_charges_lean_cert.json`).

## Reproduction & Verification Commands

```bash
# Offline deterministic replay & self-bound certificate verification (zero Lean install required)
python3 -I -B formal_pure_init/compile_point_charges_lean_cert.py --check --offline --json

# Standalone canonical JSON certificate emission (matches point_charges_lean_cert.json byte-for-byte)
python3 -I -B formal_pure_init/compile_point_charges_lean_cert.py --standalone

# Full cold Lean 4.36.0 compilation + AxiomAudit.lean + leanchecker --fresh replay
python3 -I -B formal_pure_init/compile_point_charges_lean_cert.py --check --json

# Unit, AST-purity, char-literal-escape, and tamper-resistance test suite
python3 -I -B formal_pure_init/test_point_charges_lean_cert.py -v
```
