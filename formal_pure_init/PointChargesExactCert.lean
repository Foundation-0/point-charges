-- Pure-Init Lean 4.36.0 Point-Charge Exact Certificate Module
-- Deterministically compiled from exact rational witnesses (Foundation-0/point-charges#6).
-- Zero external imports, zero Mathlib, zero floats, zero sorry, zero Classical.choice,
-- zero bv_decide, zero native_decide, zero warnings.

namespace PointChargesExactCert

-- ============================================================================
-- 0. Core Exact Integer / Rational Structures (Pure Init)
-- ============================================================================

structure Vec3Z where
  x : Int
  y : Int
  z : Int
  deriving DecidableEq

def vecZero : Vec3Z := ⟨0, 0, 0⟩
def e1 : Vec3Z := ⟨1, 0, 0⟩
def e2 : Vec3Z := ⟨0, 1, 0⟩
def e3 : Vec3Z := ⟨0, 0, 1⟩

def vecAdd (a b : Vec3Z) : Vec3Z := ⟨a.x + b.x, a.y + b.y, a.z + b.z⟩
def vecScale (s : Int) (v : Vec3Z) : Vec3Z := ⟨s * v.x, s * v.y, s * v.z⟩
def vecNormSq (v : Vec3Z) : Int := v.x * v.x + v.y * v.y + v.z * v.z

structure Sym3Z where
  h11 : Int
  h22 : Int
  h33 : Int
  h12 : Int
  h13 : Int
  h23 : Int
  deriving DecidableEq

def symZero : Sym3Z := ⟨0, 0, 0, 0, 0, 0⟩
def diagSym (d1 d2 d3 : Int) : Sym3Z := ⟨d1, d2, d3, 0, 0, 0⟩
def symAdd (a b : Sym3Z) : Sym3Z :=
  ⟨a.h11 + b.h11, a.h22 + b.h22, a.h33 + b.h33, a.h12 + b.h12, a.h13 + b.h13, a.h23 + b.h23⟩
def symTrace (h : Sym3Z) : Int := h.h11 + h.h22 + h.h33
def symDet (h : Sym3Z) : Int :=
  h.h11 * (h.h22 * h.h33 - h.h23 * h.h23)
  - h.h12 * (h.h12 * h.h33 - h.h23 * h.h13)
  + h.h13 * (h.h12 * h.h23 - h.h22 * h.h13)

def quadFormSym (h : Sym3Z) (v : Vec3Z) : Int :=
  h.h11 * v.x * v.x + h.h22 * v.y * v.y + h.h33 * v.z * v.z
  + 2 * h.h12 * v.x * v.y + 2 * h.h13 * v.x * v.z + 2 * h.h23 * v.y * v.z

structure Mat3Z where
  a11 : Int
  a12 : Int
  a13 : Int
  a21 : Int
  a22 : Int
  a23 : Int
  a31 : Int
  a32 : Int
  a33 : Int
  deriving DecidableEq

def matVec (s : Mat3Z) (u : Vec3Z) : Vec3Z :=
  ⟨s.a11 * u.x + s.a12 * u.y + s.a13 * u.z,
   s.a21 * u.x + s.a22 * u.y + s.a23 * u.z,
   s.a31 * u.x + s.a32 * u.y + s.a33 * u.z⟩

def bilSym (h : Sym3Z) (u v : Vec3Z) : Int :=
  h.h11 * u.x * v.x + h.h22 * u.y * v.y + h.h33 * u.z * v.z
  + h.h12 * (u.x * v.y + u.y * v.x)
  + h.h13 * (u.x * v.z + u.z * v.x)
  + h.h23 * (u.y * v.z + u.z * v.y)

def congruenceSym (s : Mat3Z) (h : Sym3Z) : Sym3Z :=
  let c1 : Vec3Z := ⟨s.a11, s.a21, s.a31⟩
  let c2 : Vec3Z := ⟨s.a12, s.a22, s.a32⟩
  let c3 : Vec3Z := ⟨s.a13, s.a23, s.a33⟩
  ⟨quadFormSym h c1, quadFormSym h c2, quadFormSym h c3,
   bilSym h c1 c2, bilSym h c1 c3, bilSym h c2 c3⟩

structure CoulombTermZ where
  w : Int
  dx : Int
  dy : Int
  dz : Int
  deriving DecidableEq

structure GradTermZ where
  w : Int
  d : Vec3Z
  deriving DecidableEq

def coulombSingleHessian (w dx dy dz : Int) : Sym3Z :=
  let r2 := dx * dx + dy * dy + dz * dz
  ⟨w * (3 * dx * dx - r2), w * (3 * dy * dy - r2), w * (3 * dz * dz - r2),
   w * (3 * dx * dy), w * (3 * dx * dz), w * (3 * dy * dz)⟩

def totalCoulombHessian : List CoulombTermZ → Sym3Z
  | [] => symZero
  | t :: ts => symAdd (coulombSingleHessian t.w t.dx t.dy t.dz) (totalCoulombHessian ts)

def totalCoulombGrad : List GradTermZ → Vec3Z
  | [] => vecZero
  | t :: ts => vecAdd (vecScale (-t.w) t.d) (totalCoulombGrad ts)

def sumPosList : List Int → Int
  | [] => 0
  | w :: ws => w + sumPosList ws

def allPosNonEmpty : List Int → Prop
  | [] => False
  | [w] => 0 < w
  | w :: w2 :: ws => 0 < w ∧ allPosNonEmpty (w2 :: ws)

def planarForceZ (z : Int) (ws : List Int) : Int := -z * sumPosList ws
def planarHessZZ (ws : List Int) : Int := -sumPosList ws
def planarInPlaneTrace (ws : List Int) : Int := sumPosList ws

structure Box3Z where
  lx : Int
  ux : Int
  ly : Int
  uy : Int
  lz : Int
  uz : Int
  deriving DecidableEq

def inBox (x : Vec3Z) (b : Box3Z) : Prop :=
  b.lx ≤ x.x ∧ x.x ≤ b.ux ∧ b.ly ≤ x.y ∧ x.y ≤ b.uy ∧ b.lz ≤ x.z ∧ x.z ≤ b.uz

def boxDisjoint (a b : Box3Z) : Bool :=
  decide (a.ux < b.lx) || decide (b.ux < a.lx) ||
  decide (a.uy < b.ly) || decide (b.uy < a.ly) ||
  decide (a.uz < b.lz) || decide (b.uz < a.lz)

def allPairsDisjoint : List Box3Z → Bool
  | [] => true
  | b :: bs => bs.all (fun c => boxDisjoint b c) && allPairsDisjoint bs

def allBoxesAvoidSites (boxes sites : List Box3Z) : Bool :=
  boxes.all (fun b => sites.all (fun s => boxDisjoint b s))

def countVal (k : Nat) (xs : List Nat) : Nat :=
  xs.foldl (fun acc x => if x == k then acc + 1 else acc) 0

structure RationalSite3 where
  xn : Int
  xd : Nat
  yn : Int
  yd : Nat
  zn : Int
  zd : Nat
  deriving DecidableEq

def siteOnUnitSphere (p : RationalSite3) : Bool :=
  let xd2 : Int := (p.xd * p.xd : Nat)
  let yd2 : Int := (p.yd * p.yd : Nat)
  let zd2 : Int := (p.zd * p.zd : Nat)
  decide (0 < p.xd ∧ 0 < p.yd ∧ 0 < p.zd) &&
  decide (p.xn * p.xn * yd2 * zd2 + p.yn * p.yn * xd2 * zd2 + p.zn * p.zn * xd2 * yd2 = xd2 * yd2 * zd2)

structure SqrtBracket where
  r2n : Nat
  r2d : Nat
  ln : Nat
  ld : Nat
  un : Nat
  ud : Nat
  deriving DecidableEq

def validSqrtBracket (b : SqrtBracket) : Bool :=
  decide (0 < b.r2n ∧ 0 < b.r2d ∧ 0 < b.ln ∧ 0 < b.ld ∧ 0 < b.un ∧ 0 < b.ud) &&
  decide (b.ln * b.ln * b.r2d ≤ b.r2n * b.ld * b.ld) &&
  decide (b.r2n * b.ud * b.ud ≤ b.un * b.un * b.r2d)

structure QNat where
  n : Nat
  d : Nat
  deriving DecidableEq

def qAdd (a b : QNat) : QNat := ⟨a.n * b.d + b.n * a.d, a.d * b.d⟩
def qLt (a b : QNat) : Bool := decide (a.n * b.d < b.n * a.d)
def sumE2 (bs : List SqrtBracket) : QNat :=
  bs.foldl (fun acc b => qAdd acc ⟨b.r2d, b.r2n⟩) ⟨0, 1⟩
def sumE1Lo (bs : List SqrtBracket) : QNat :=
  bs.foldl (fun acc b => qAdd acc ⟨b.ud, b.un⟩) ⟨0, 1⟩
def sumE1Hi (bs : List SqrtBracket) : QNat :=
  bs.foldl (fun acc b => qAdd acc ⟨b.ld, b.ln⟩) ⟨0, 1⟩

def checkPairDistSq (p q : RationalSite3) (r2n r2d : Nat) : Bool :=
  let dxn := p.xn * (q.xd : Int) - q.xn * (p.xd : Int)
  let dxd : Int := (p.xd * q.xd : Nat)
  let dyn := p.yn * (q.yd : Int) - q.yn * (p.yd : Int)
  let dyd : Int := (p.yd * q.yd : Nat)
  let dzn := p.zn * (q.zd : Int) - q.zn * (p.zd : Int)
  let dzd : Int := (p.zd * q.zd : Nat)
  let num := dxn * dxn * dyd * dyd * dzd * dzd + dyn * dyn * dxd * dxd * dzd * dzd + dzn * dzn * dxd * dxd * dyd * dyd
  let den := dxd * dxd * dyd * dyd * dzd * dzd
  decide (num * (r2d : Int) = (r2n : Int) * den)

def pairDistsMatch : List RationalSite3 → List SqrtBracket → Bool
  | [], [] => true
  | p :: ps, bs =>
    let k := ps.length
    let head := bs.take k
    let tail := bs.drop k
    (head.length == k) &&
    ((ps.zip head).all (fun (q, b) => checkPairDistSq p q b.r2n b.r2d)) &&
    pairDistsMatch ps tail
  | _, _ => false

/-- Rational point on the 3-sphere S^3(Q) stored as four (numerator, positive denominator) pairs. -/
structure RationalSite4 where
  wn : Int
  wd : Nat
  xn : Int
  xd : Nat
  yn : Int
  yd : Nat
  zn : Int
  zd : Nat
  deriving DecidableEq

/-- Exact rational equality with positive-denominator guards: aN / aD = bN / bD. -/
def ratEq (aN : Int) (aD : Nat) (bN : Int) (bD : Nat) : Bool :=
  decide (0 < aD ∧ 0 < bD ∧ aN * (bD : Int) = bN * (aD : Int))

/-- Cleared-denominator S^3(Q) unit-norm check: w^2 + x^2 + y^2 + z^2 = 1. -/
def siteOnUnitSphere4 (q : RationalSite4) : Bool :=
  let wd2 : Int := (q.wd * q.wd : Nat)
  let xd2 : Int := (q.xd * q.xd : Nat)
  let yd2 : Int := (q.yd * q.yd : Nat)
  let zd2 : Int := (q.zd * q.zd : Nat)
  decide (0 < q.wd ∧ 0 < q.xd ∧ 0 < q.yd ∧ 0 < q.zd) &&
  decide (
    q.wn * q.wn * xd2 * yd2 * zd2 +
    q.xn * q.xn * wd2 * yd2 * zd2 +
    q.yn * q.yn * wd2 * xd2 * zd2 +
    q.zn * q.zn * wd2 * xd2 * yd2 = wd2 * xd2 * yd2 * zd2
  )

/-- Verify that `p` equals the exact inverse stereographic image of `(un/ud, vn/vd)` on S^2(Q). -/
def stereographicMatches (un : Int) (ud : Nat) (vn : Int) (vd : Nat) (p : RationalSite3) : Bool :=
  let ud2 : Int := (ud * ud : Nat)
  let vd2 : Int := (vd * vd : Nat)
  let sNum : Int := un * un * vd2 + vn * vn * ud2
  let dDen : Int := ud2 * vd2
  let qDen : Int := sNum + dDen
  decide (0 < ud ∧ 0 < vd ∧ 0 < p.xd ∧ 0 < p.yd ∧ 0 < p.zd ∧ 0 < qDen) &&
  decide (p.xn * qDen = 2 * un * (ud : Int) * vd2 * (p.xd : Int)) &&
  decide (p.yn * qDen = 2 * vn * (vd : Int) * ud2 * (p.yd : Int)) &&
  decide (p.zn * qDen = (sNum - dDen) * (p.zd : Int))

/-- Verify that `q` equals the exact Cayley S^3(Q) map of `(w1n/w1d, w2n/w2d, w3n/w3d)`. -/
def cayleyS3Matches (w1n : Int) (w1d : Nat) (w2n : Int) (w2d : Nat) (w3n : Int) (w3d : Nat) (q : RationalSite4) : Bool :=
  let d1 : Int := (w1d * w1d : Nat)
  let d2 : Int := (w2d * w2d : Nat)
  let d3 : Int := (w3d * w3d : Nat)
  let sNum : Int := w1n * w1n * d2 * d3 + w2n * w2n * d1 * d3 + w3n * w3n * d1 * d2
  let dDen : Int := d1 * d2 * d3
  let qDen : Int := dDen + sNum
  decide (0 < w1d ∧ 0 < w2d ∧ 0 < w3d ∧ 0 < q.wd ∧ 0 < q.xd ∧ 0 < q.yd ∧ 0 < q.zd ∧ 0 < qDen) &&
  decide (q.wn * qDen = (dDen - sNum) * (q.wd : Int)) &&
  decide (q.xn * qDen = 2 * w1n * (w1d : Int) * d2 * d3 * (q.xd : Int)) &&
  decide (q.yn * qDen = 2 * w2n * (w2d : Int) * d1 * d3 * (q.yd : Int)) &&
  decide (q.zn * qDen = 2 * w3n * (w3d : Int) * d1 * d2 * (q.zd : Int))

/-- Verify the 4-coordinate equatorial S^3(Q) <-> S^2(Q) bridge: (q_w, q_x, q_y, q_z) = (-p_z, p_x, p_y, 0). -/
def equatorialCayleyMatchesStereo (q : RationalSite4) (p : RationalSite3) : Bool :=
  ratEq q.wn q.wd (-p.zn) p.zd &&
  ratEq q.xn q.xd p.xn p.xd &&
  ratEq q.yn q.yd p.yn p.yd &&
  ratEq q.zn q.zd 0 1

-- ============================================================================
-- 1. Fourteen Universal (forall) Pure-Init Electrostatic, Box, & Sphere Theorems
-- ============================================================================

private theorem two_mul_int (x : Int) : 2 * x = x + x := by omega
private theorem three_mul_int (x : Int) : 3 * x = x + x + x := by omega
private theorem four_mul_int (x : Int) : 4 * x = x + x + x + x := by omega
private theorem mul_four_int (x : Int) : x * 4 = x + x + x + x := by omega

private theorem symTrace_add (a b : Sym3Z) :
    symTrace (symAdd a b) = symTrace a + symTrace b := by
  simp only [symTrace, symAdd]
  omega

theorem t01_coulomb_single_hessian_trace_zero (w dx dy dz : Int) :
    symTrace (coulombSingleHessian w dx dy dz) = 0 := by
  simp only [symTrace, coulombSingleHessian]
  rw [← Int.mul_add, ← Int.mul_add]
  have h0 : 3 * dx * dx - (dx * dx + dy * dy + dz * dz)
      + (3 * dy * dy - (dx * dx + dy * dy + dz * dz))
      + (3 * dz * dz - (dx * dx + dy * dy + dz * dz)) = 0 := by
    rw [Int.mul_assoc 3 dx dx, Int.mul_assoc 3 dy dy, Int.mul_assoc 3 dz dz]
    omega
  rw [h0, Int.mul_zero]

theorem t02_earnshaw_total_hessian_trace_zero (terms : List CoulombTermZ) :
    symTrace (totalCoulombHessian terms) = 0 := by
  induction terms with
  | nil => rfl
  | cons t ts ih =>
    simp only [totalCoulombHessian]
    rw [symTrace_add, t01_coulomb_single_hessian_trace_zero, ih]
    omega

private theorem quadFormSym_e1 (h : Sym3Z) : quadFormSym h e1 = h.h11 := by
  simp only [quadFormSym, e1, Int.mul_one, Int.mul_zero, Int.add_zero]

private theorem quadFormSym_e2 (h : Sym3Z) : quadFormSym h e2 = h.h22 := by
  simp only [quadFormSym, e2, Int.mul_one, Int.mul_zero, Int.zero_add, Int.add_zero]

private theorem quadFormSym_e3 (h : Sym3Z) : quadFormSym h e3 = h.h33 := by
  simp only [quadFormSym, e3, Int.mul_one, Int.mul_zero, Int.zero_add, Int.add_zero]

theorem t03_earnshaw_no_strict_definite_quadform (h : Sym3Z) (htr : symTrace h = 0) :
    (¬ ∀ v : Vec3Z, v ≠ vecZero → 0 < quadFormSym h v) ∧
    (¬ ∀ v : Vec3Z, v ≠ vecZero → quadFormSym h v < 0) := by
  have he1 : e1 ≠ vecZero := by decide
  have he2 : e2 ≠ vecZero := by decide
  have he3 : e3 ≠ vecZero := by decide
  constructor
  · intro hpd
    have h1 := hpd e1 he1; rw [quadFormSym_e1] at h1
    have h2 := hpd e2 he2; rw [quadFormSym_e2] at h2
    have h3 := hpd e3 he3; rw [quadFormSym_e3] at h3
    simp only [symTrace] at htr
    omega
  · intro hnd
    have h1 := hnd e1 he1; rw [quadFormSym_e1] at h1
    have h2 := hnd e2 he2; rw [quadFormSym_e2] at h2
    have h3 := hnd e3 he3; rw [quadFormSym_e3] at h3
    simp only [symTrace] at htr
    omega

private theorem sumPosList_pos (ws : List Int) (h : allPosNonEmpty ws) : 0 < sumPosList ws := by
  match ws with
  | [] => exact False.elim h
  | [w] => simp only [allPosNonEmpty, sumPosList] at h ⊢; omega
  | w :: w2 :: rest =>
    simp only [allPosNonEmpty] at h
    have ih := sumPosList_pos (w2 :: rest) h.2
    change 0 < w + sumPosList (w2 :: rest)
    omega

theorem t04_planar_offplane_force_and_z_curvature (z : Int) (ws : List Int) (hpos : allPosNonEmpty ws) :
    (0 < z → planarForceZ z ws < 0) ∧
    (z < 0 → 0 < planarForceZ z ws) ∧
    planarHessZZ ws < 0 ∧
    0 < planarInPlaneTrace ws := by
  have hs : 0 < sumPosList ws := sumPosList_pos ws hpos
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro hz
    simp only [planarForceZ]
    have hp : 0 < z * sumPosList ws := Int.mul_pos hz hs
    rw [Int.neg_mul]; omega
  · intro hz
    simp only [planarForceZ]
    have hnz : 0 < -z := by omega
    exact Int.mul_pos hnz hs
  · simp only [planarHessZZ]; omega
  · simp only [planarInPlaneTrace]; exact hs

theorem t05_lagrange_congruence_quadform_identity (s : Mat3Z) (h : Sym3Z) (u : Vec3Z) :
    quadFormSym h (matVec s u) = quadFormSym (congruenceSym s h) u := by
  simp only [quadFormSym, matVec, congruenceSym, bilSym]
  simp only [two_mul_int, Int.add_mul, Int.mul_add]
  simp only [Int.mul_comm, Int.mul_left_comm]
  omega

private theorem quadFormSym_diag (d1 d2 d3 : Int) (u : Vec3Z) :
    quadFormSym (diagSym d1 d2 d3) u = d1 * u.x * u.x + d2 * u.y * u.y + d3 * u.z * u.z := by
  simp only [quadFormSym, diagSym, Int.mul_zero, Int.zero_mul, Int.add_zero]

theorem t06_lagrange_diagonal_morse_index2_witness (s : Mat3Z) (h : Sym3Z) (d1 d2 d3 : Int)
    (hcong : congruenceSym s h = diagSym d1 d2 d3) (hd1 : 0 < d1) (hd2 : d2 < 0) (hd3 : d3 < 0) :
    0 < quadFormSym h (matVec s e1) ∧ quadFormSym h (matVec s e2) < 0 ∧ quadFormSym h (matVec s e3) < 0 := by
  rw [t05_lagrange_congruence_quadform_identity, t05_lagrange_congruence_quadform_identity,
      t05_lagrange_congruence_quadform_identity, hcong]
  rw [quadFormSym_diag, quadFormSym_diag, quadFormSym_diag]
  simp only [e1, e2, e3, Int.mul_one, Int.mul_zero, Int.add_zero, Int.zero_add]
  exact ⟨hd1, hd2, hd3⟩

theorem t07_lagrange_diagonal_morse_index1_witness (s : Mat3Z) (h : Sym3Z) (d1 d2 d3 : Int)
    (hcong : congruenceSym s h = diagSym d1 d2 d3) (hd1 : 0 < d1) (hd2 : 0 < d2) (hd3 : d3 < 0) :
    0 < quadFormSym h (matVec s e1) ∧ 0 < quadFormSym h (matVec s e2) ∧ quadFormSym h (matVec s e3) < 0 := by
  rw [t05_lagrange_congruence_quadform_identity, t05_lagrange_congruence_quadform_identity,
      t05_lagrange_congruence_quadform_identity, hcong]
  rw [quadFormSym_diag, quadFormSym_diag, quadFormSym_diag]
  simp only [e1, e2, e3, Int.mul_one, Int.mul_zero, Int.add_zero, Int.zero_add]
  exact ⟨hd1, hd2, hd3⟩

theorem t08_disjoint_boxes_no_common_point (b1 b2 : Box3Z) (hd : boxDisjoint b1 b2 = true)
    (x : Vec3Z) : ¬(inBox x b1 ∧ inBox x b2) := by
  intro ⟨h1, h2⟩
  simp only [inBox] at h1 h2
  simp only [boxDisjoint, Bool.or_eq_true, decide_eq_true_eq] at hd
  omega

theorem t09_inv_sqrt_bracket_order_preservation (ln ld un ud r2n r2d : Nat)
    (hr2d : 0 < r2d)
    (hlo : ln * ln * r2d ≤ r2n * ld * ld) (hhi : r2n * ud * ud ≤ un * un * r2d) :
    ln * ud ≤ un * ld := by
  have hmul1 : (ln * ln * r2d) * (ud * ud) ≤ (r2n * ld * ld) * (ud * ud) := Nat.mul_le_mul_right _ hlo
  have hmul2 : (r2n * ud * ud) * (ld * ld) ≤ (un * un * r2d) * (ld * ld) := Nat.mul_le_mul_right _ hhi
  have hmid : (r2n * ld * ld) * (ud * ud) = (r2n * ud * ud) * (ld * ld) := by
    simp only [Nat.mul_left_comm, Nat.mul_comm]
  have hchain : (ln * ln * r2d) * (ud * ud) ≤ (un * un * r2d) * (ld * ld) := by omega
  have hlhs : (ln * ln * r2d) * (ud * ud) = r2d * ((ln * ud) * (ln * ud)) := by
    simp only [Nat.mul_left_comm, Nat.mul_comm]
  have hrhs : (un * un * r2d) * (ld * ld) = r2d * ((un * ld) * (un * ld)) := by
    simp only [Nat.mul_left_comm, Nat.mul_comm]
  rw [hlhs, hrhs] at hchain
  have hsq : (ln * ud) * (ln * ud) ≤ (un * ld) * (un * ld) := Nat.le_of_mul_le_mul_left hchain hr2d
  if hle : ln * ud ≤ un * ld then
    exact hle
  else
    have hgt : un * ld < ln * ud := by omega
    have h1 : (un * ld) * (un * ld) < (ln * ud) * (ln * ud) := Nat.mul_lt_mul_of_lt_of_lt hgt hgt
    omega

theorem t10_stereographic_s2q_unit_sphere (un ud vn vd : Int) :
    let a := 2 * un * ud * vd * vd
    let b := 2 * vn * vd * ud * ud
    let c := (un * un * vd * vd + vn * vn * ud * ud) - ud * ud * vd * vd
    let d := (un * un * vd * vd + vn * vn * ud * ud) + ud * ud * vd * vd
    a * a + b * b + c * c = d * d := by
  dsimp only
  simp only [two_mul_int, Int.sub_mul, Int.mul_sub, Int.add_mul, Int.mul_add]
  simp only [Int.mul_comm, Int.mul_left_comm]
  omega

theorem t11_stereographic_s2q_chord_distance_sq (u1 v1 u2 v2 : Int) :
    let d1 := 1 + u1 * u1 + v1 * v1
    let d2 := 1 + u2 * u2 + v2 * v2
    let dx := 2 * u1 * d2 - 2 * u2 * d1
    let dy := 2 * v1 * d2 - 2 * v2 * d1
    let dz := (u1 * u1 + v1 * v1 - 1) * d2 - (u2 * u2 + v2 * v2 - 1) * d1
    dx * dx + dy * dy + dz * dz = 4 * ((u1 - u2) * (u1 - u2) + (v1 - v2) * (v1 - v2)) * d1 * d2 := by
  dsimp only
  simp only [two_mul_int, four_mul_int, Int.sub_mul, Int.mul_sub, Int.add_mul, Int.mul_add, Int.mul_one, Int.one_mul]
  simp only [Int.mul_comm, Int.mul_left_comm]
  omega

theorem t12_cayley_s3q_unit_norm_and_equatorial_bridge (w1 w2 w3 : Int) :
    let s := w1 * w1 + w2 * w2 + w3 * w3
    ((1 - s) * (1 - s) + (2 * w1) * (2 * w1) + (2 * w2) * (2 * w2) + (2 * w3) * (2 * w3)
      = (1 + s) * (1 + s)) ∧
    (let s2 := w1 * w1 + w2 * w2; (1 - s2) = -(s2 - 1)) := by
  dsimp only
  constructor
  · simp only [two_mul_int, Int.sub_mul, Int.mul_sub, Int.add_mul, Int.mul_add, Int.mul_one, Int.one_mul]
    simp only [Int.mul_comm, Int.mul_left_comm]
    omega
  · omega

theorem t13_n3_isosceles_cubic_cancellation_and_pythagorean_roots (h y k : Int) :
    ((h - y) * (h - 3 * y) * (1 + y * y) - 3 * y * (h - y) * (h - y) * y
      = (h - y) * (h - 3 * y - 2 * h * y * y)) ∧
    (let hNum := 3 * k; let hDen := k * k - 2;
      hNum * k * k - 3 * hDen * k - 2 * hNum = 0 ∧
      2 * hNum * k * k - 3 * hDen * (2 * k) - hNum * 4 = 0) := by
  dsimp only
  refine ⟨?_, ?_, ?_⟩
  · simp only [two_mul_int, three_mul_int, Int.sub_mul, Int.mul_sub, Int.add_mul, Int.mul_add, Int.mul_one]
    simp only [Int.mul_comm, Int.mul_left_comm]
    omega
  · simp only [two_mul_int, three_mul_int, Int.sub_mul, Int.mul_sub, Int.add_mul, Int.mul_add]
    simp only [Int.mul_comm]
    omega
  · simp only [two_mul_int, three_mul_int, mul_four_int, Int.sub_mul, Int.mul_sub, Int.add_mul, Int.mul_add]
    simp only [Int.mul_comm]
    omega

theorem t14_regular_tetrahedron_s2q_sqrt3_obstruction : ∀ (a b : Nat), 0 < b → 3 * (b * b) ≠ a * a := by
  intro a
  induction a using Nat.strongRecOn with
  | ind a ih =>
    intro b hb hEq
    have ha3 : a % 3 = 0 := by
      have hmod : (a * a) % 3 = 0 := by omega
      have hm := Nat.mul_mod a a 3
      have hlt : a % 3 < 3 := Nat.mod_lt a (by decide)
      rcases (show a % 3 = 0 ∨ a % 3 = 1 ∨ a % 3 = 2 by omega) with h0 | h1 | h2
      · exact h0
      · rw [h1] at hm; omega
      · rw [h2] at hm; omega
    have ha_eq : a = 3 * (a / 3) := by omega
    have hb_eq : b * b = 3 * ((a / 3) * (a / 3)) := by
      have h1 : 3 * (b * b) = 9 * ((a / 3) * (a / 3)) := by
        rw [hEq]
        conv => lhs; rw [ha_eq]
        rw [Nat.mul_assoc, Nat.mul_left_comm (a / 3) 3 (a / 3), ← Nat.mul_assoc]
      omega
    have hb3 : b % 3 = 0 := by
      have hmod : (b * b) % 3 = 0 := by omega
      have hm := Nat.mul_mod b b 3
      have hlt : b % 3 < 3 := Nat.mod_lt b (by decide)
      rcases (show b % 3 = 0 ∨ b % 3 = 1 ∨ b % 3 = 2 by omega) with h0 | h1 | h2
      · exact h0
      · rw [h1] at hm; omega
      · rw [h2] at hm; omega
    have hba : a / 3 < a := by
      have ha_pos : 0 < a := by
        rcases Nat.eq_zero_or_pos a with rfl | hpos
        · have : 0 < b * b := Nat.mul_pos hb hb
          omega
        · exact hpos
      omega
    have hb_div_pos : 0 < b / 3 := by omega
    have hb_sub : 3 * ((b / 3) * (b / 3)) = (a / 3) * (a / 3) := by
      have hb_eq2 : b = 3 * (b / 3) := by omega
      have h2 : b * b = 9 * ((b / 3) * (b / 3)) := by
        conv => lhs; rw [hb_eq2]
        rw [Nat.mul_assoc, Nat.mul_left_comm (b / 3) 3 (b / 3), ← Nat.mul_assoc]
      omega
    exact ih (a / 3) hba (b / 3) hb_div_pos hb_sub

-- ============================================================================
-- 2. Eight Compiled Exact-Rational Witness Theorems (Emitted From Witness Pack)
-- ============================================================================

def w_n2_grad : List GradTermZ := [⟨1, ⟨1, 0, 0⟩⟩, ⟨1, ⟨(-1 : Int), 0, 0⟩⟩]
def w_n2_hess : List CoulombTermZ := [⟨1, 1, 0, 0⟩, ⟨1, (-1 : Int), 0, 0⟩]

theorem t15_compiled_n2_analytic_equilibrium_and_hessian :
    totalCoulombGrad w_n2_grad = vecZero ∧
    totalCoulombHessian w_n2_hess = diagSym 4 (-2 : Int) (-2 : Int) ∧
    symTrace (totalCoulombHessian w_n2_hess) = 0 ∧
    symDet (totalCoulombHessian w_n2_hess) ≠ 0 ∧
    0 < quadFormSym (totalCoulombHessian w_n2_hess) e1 ∧
    quadFormSym (totalCoulombHessian w_n2_hess) e2 < 0 ∧
    quadFormSym (totalCoulombHessian w_n2_hess) e3 < 0 := by
  decide

def w_n3_iso_grad : List GradTermZ := [⟨100, ⟨175, 0, 0⟩⟩, ⟨100, ⟨(-49 : Int), 168, 0⟩⟩, ⟨56, ⟨(-225 : Int), (-300 : Int), 0⟩⟩]
def w_n3_iso_hess : List CoulombTermZ := [⟨22500, 175, 0, 0⟩, ⟨22500, (-49 : Int), 168, 0⟩, ⟨2744, (-225 : Int), (-300 : Int), 0⟩]

theorem t16_compiled_n3_weierstrass_incenter_iso_acute :
    totalCoulombGrad w_n3_iso_grad = vecZero ∧
    totalCoulombHessian w_n3_iso_hess = diagSym 882000000 882000000 (-1764000000 : Int) ∧
    symTrace (totalCoulombHessian w_n3_iso_hess) = 0 ∧
    symDet (totalCoulombHessian w_n3_iso_hess) ≠ 0 ∧
    0 < quadFormSym (totalCoulombHessian w_n3_iso_hess) e1 ∧
    0 < quadFormSym (totalCoulombHessian w_n3_iso_hess) e2 ∧
    quadFormSym (totalCoulombHessian w_n3_iso_hess) e3 < 0 := by
  decide

def w_n3_sca_grad : List GradTermZ := [⟨1014, ⟨10115, 0, 0⟩⟩, ⟨867, ⟨(-7840 : Int), 9576, 0⟩⟩, ⟨1197, ⟨(-2890 : Int), (-6936 : Int), 0⟩⟩]
def w_n3_sca_hess : List CoulombTermZ := [⟨10967424, 10115, 0, 0⟩, ⟨6264075, (-7840 : Int), 9576, 0⟩, ⟨23461200, (-2890 : Int), (-6936 : Int), 0⟩]

theorem t17_compiled_n3_weierstrass_incenter_scalene :
    totalCoulombGrad w_n3_sca_grad = vecZero ∧
    totalCoulombHessian w_n3_sca_hess = diagSym 1703087936222400 1703087936222400 (-3406175872444800 : Int) ∧
    symTrace (totalCoulombHessian w_n3_sca_hess) = 0 ∧
    symDet (totalCoulombHessian w_n3_sca_hess) ≠ 0 ∧
    0 < quadFormSym (totalCoulombHessian w_n3_sca_hess) e1 ∧
    0 < quadFormSym (totalCoulombHessian w_n3_sca_hess) e2 ∧
    quadFormSym (totalCoulombHessian w_n3_sca_hess) e3 < 0 := by
  decide

def w_offplane_weights : List Int := [1, 1, 1]

theorem t18_compiled_n3_pythagorean_off_plane_confinement :
    planarForceZ 4 w_offplane_weights < 0 ∧
    planarHessZZ w_offplane_weights < 0 ∧
    0 < planarInPlaneTrace w_offplane_weights := by
  decide

def w_n4_sites : List Box3Z := [
  ⟨(-386808190652 : Int), (-386808190651 : Int), (-300969317872 : Int), (-300969317871 : Int), (-221716519742 : Int), (-221716519741 : Int)⟩,
  ⟨(-62034446040 : Int), (-62034446039 : Int), 490470146918, 490470146919, (-235691312531 : Int), (-235691312530 : Int)⟩,
  ⟨(-3012661861 : Int), (-3012661860 : Int), 7751556975, 7751556976, 930021910354, 930021910355⟩,
  ⟨449282441341, 449282441342, (-191359003699 : Int), (-191359003698 : Int), (-217142551370 : Int), (-217142551369 : Int)⟩
]
def w_n4_boxes : List Box3Z := [
  ⟨(-19207522901 : Int), (-19207522898 : Int), (-14424789389 : Int), (-14424789386 : Int), 17158819990, 17158819993⟩,
  ⟨(-5939541728 : Int), (-5939541725 : Int), 2326302587, 2326302590, (-5647147432 : Int), (-5647147429 : Int)⟩,
  ⟨(-5352873505 : Int), (-5352873502 : Int), (-3957788058 : Int), (-3957788055 : Int), 4691547086, 4691547089⟩,
  ⟨(-3376824206 : Int), (-3376824203 : Int), 25330256727, 25330256730, 17657870133, 17657870136⟩,
  ⟨(-971817672 : Int), (-971817669 : Int), 6765432912, 6765432915, 4739135691, 4739135694⟩,
  ⟨(-233280007 : Int), (-233280004 : Int), 236621840, 236621843, 1700045808, 1700045811⟩,
  ⟨156214821, 156214824, (-599745119 : Int), (-599745116 : Int), (-2757398988 : Int), (-2757398985 : Int)⟩,
  ⟨454961082, 454961085, (-3608502087 : Int), (-3608502084 : Int), (-3760315553 : Int), (-3760315550 : Int)⟩,
  ⟨4340929498, 4340929501, 3160236795, 3160236798, (-4995974421 : Int), (-4995974418 : Int)⟩,
  ⟨6415395664, 6415395667, (-2676840557 : Int), (-2676840554 : Int), 5008269348, 5008269351⟩,
  ⟨21383785193, 21383785196, (-8852553856 : Int), (-8852553853 : Int), 16810242638, 16810242641⟩
]
def w_n4_morse : List Nat := [2, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2]

theorem t19_compiled_n4_eleven_disjointness_and_morse :
    w_n4_boxes.length = 11 ∧
    allPairsDisjoint w_n4_boxes = true ∧
    allBoxesAvoidSites w_n4_boxes w_n4_sites = true ∧
    countVal 1 w_n4_morse = 4 ∧ countVal 2 w_n4_morse = 7 ∧
    countVal 2 w_n4_morse - countVal 1 w_n4_morse = 4 - 1 ∧
    (4 - 1) * (4 - 1) < w_n4_boxes.length := by
  decide

def w_n8_sites : List Box3Z := [
  ⟨3298534884614, 3298534884615, (-713 : Int), (-712 : Int), (-1209 : Int), (-1208 : Int)⟩,
  ⟨(-420 : Int), (-419 : Int), 3298534885337, 3298534885338, (-528 : Int), (-527 : Int)⟩,
  ⟨(-1064 : Int), (-1063 : Int), (-676 : Int), (-675 : Int), 3298534884639, 3298534884640⟩,
  ⟨1374389534720, 1374389534720, 1374389534720, 1374389534720, 1374389534720, 1374389534720⟩,
  ⟨824633720832, 824633720832, 824633720832, 824633720832, 824633720832, 824633720832⟩,
  ⟨1100048498688, 1100048498688, 1100048498688, 1100048498688, 1100048498688, 1100048498688⟩,
  ⟨1098974756864, 1098974756864, 1098974756864, 1098974756864, 1098974756864, 1098974756864⟩,
  ⟨1099518311421, 1099518311422, 1099518311421, 1099518311422, 1099518311421, 1099518311422⟩
]
def w_n8_boxes : List Box3Z := [
  ⟨674789642115, 674789642116, 1311872620797, 1311872620798, 1311872620563, 1311872620564⟩,
  ⟨831986169413, 831986169414, 1233274357140, 1233274357141, 1233274356910, 1233274356911⟩,
  ⟨1059217960396, 1059217960397, 1059217960288, 1059217960289, 1100151021191, 1100151021192⟩,
  ⟨1059217963993, 1059217963994, 1100151020345, 1100151020346, 1059217963909, 1059217963910⟩,
  ⟨1084164486308, 1084164486309, 1084164485767, 1084164485768, 1092116227293, 1092116227294⟩,
  ⟨1084164482510, 1084164482511, 1092116229940, 1092116229941, 1084164482373, 1084164482374⟩,
  ⟨1087677942742, 1087677942743, 1087677940343, 1087677940344, 1087677942448, 1087677942449⟩,
  ⟨1092116227166, 1092116227167, 1084164486049, 1084164486050, 1084164486452, 1084164486453⟩,
  ⟨1097341313407, 1097341313408, 1100596785189, 1100596785190, 1100596784963, 1100596784964⟩,
  ⟨1099328762277, 1099328762278, 1099603060688, 1099603060689, 1099603060460, 1099603060461⟩,
  ⟨1099470353546, 1099470353547, 1099470353401, 1099470353402, 1099511943348, 1099511943349⟩,
  ⟨1099470357594, 1099470357595, 1099511942313, 1099511942314, 1099470357505, 1099470357506⟩,
  ⟨1099493966987, 1099493966988, 1099493966415, 1099493966416, 1099503935041, 1099503935042⟩,
  ⟨1099493962578, 1099493962579, 1099503937920, 1099503937921, 1099493962438, 1099493962439⟩,
  ⟨1099498484881, 1099498484882, 1099498482314, 1099498482315, 1099498484568, 1099498484569⟩,
  ⟨1099503934896, 1099503934897, 1099493966738, 1099493966739, 1099493967169, 1099493967170⟩,
  ⟨1099513762055, 1099513762056, 1099513771133, 1099513771134, 1099513770806, 1099513770807⟩,
  ⟨1099511943402, 1099511943403, 1099470353368, 1099470353369, 1099470353425, 1099470353426⟩,
  ⟨1099513770797, 1099513770798, 1099513761886, 1099513761887, 1099513770684, 1099513770685⟩,
  ⟨1099513770835, 1099513770836, 1099513771045, 1099513771046, 1099513762055, 1099513762056⟩,
  ⟨1099522671571, 1099522671572, 1099522672265, 1099522672266, 1099533662936, 1099533662937⟩,
  ⟨1099522671923, 1099522671924, 1099533660802, 1099533660803, 1099522671937, 1099522671938⟩,
  ⟨1099525155067, 1099525155068, 1099525156642, 1099525156643, 1099525155259, 1099525155260⟩,
  ⟨1099525172305, 1099525172306, 1099525172756, 1099525172757, 1099566733549, 1099566733550⟩,
  ⟨1099525172580, 1099525172581, 1099566735170, 1099566735171, 1099525172564, 1099525172565⟩,
  ⟨1099533663717, 1099533663718, 1099522672200, 1099522672201, 1099522671520, 1099522671521⟩,
  ⟨1099566732639, 1099566732640, 1099525172655, 1099525172656, 1099525172188, 1099525172189⟩,
  ⟨1099603060790, 1099603060791, 1099328761539, 1099328761540, 1099603060741, 1099603060742⟩,
  ⟨1099603060381, 1099603060382, 1099603060560, 1099603060561, 1099328762440, 1099328762441⟩,
  ⟨1100151021242, 1100151021243, 1059217960265, 1059217960266, 1059217960289, 1059217960290⟩,
  ⟨1100596784587, 1100596784588, 1097341313797, 1097341313798, 1100596784538, 1100596784539⟩,
  ⟨1100596785046, 1100596785047, 1100596785223, 1100596785224, 1097341313212, 1097341313213⟩,
  ⟨1109557608010, 1109557608011, 1109557608906, 1109557608907, 1117509352454, 1117509352455⟩,
  ⟨1109557608638, 1109557608639, 1117509349365, 1117509349366, 1109557608676, 1109557608677⟩,
  ⟨1111345312809, 1111345312810, 1111345315215, 1111345315216, 1111345313101, 1111345313102⟩,
  ⟨1112516587540, 1112516587541, 1112516588002, 1112516588003, 1153449646869, 1153449646870⟩,
  ⟨1112516587854, 1112516587855, 1153449648582, 1153449648583, 1112516587839, 1112516587840⟩,
  ⟨1117509353345, 1117509353346, 1109557608794, 1109557608795, 1109557607936, 1109557607937⟩,
  ⟨1153449646028, 1153449646029, 1112516587899, 1112516587900, 1112516587421, 1112516587422⟩,
  ⟨1233274358051, 1233274358052, 831986167045, 831986167046, 1233274358001, 1233274358002⟩,
  ⟨1233274356607, 1233274356608, 1233274356787, 1233274356788, 831986170025, 831986170026⟩,
  ⟨1311872619473, 1311872619474, 674789644184, 674789644185, 1311872619421, 1311872619422⟩,
  ⟨1311872620880, 1311872620881, 1311872621063, 1311872621064, 674789641484, 674789641485⟩,
  ⟨1099513771475, 1099513771476, 1099513777592, 1099513777593, 1099513771577, 1099513771578⟩,
  ⟨1099513784014, 1099513784015, 1099513771734, 1099513771735, 1099513771331, 1099513771332⟩,
  ⟨1099513772026, 1099513772027, 1099513773481, 1099513773482, 1099513772440, 1099513772441⟩,
  ⟨1099513770036, 1099513770037, 1099513769662, 1099513769663, 1099513769592, 1099513769593⟩,
  ⟨1099513771365, 1099513771366, 1099513771845, 1099513771846, 1099513779892, 1099513779893⟩,
  ⟨1099514376325, 1099514376326, 1099513783116, 1099513783117, 1099513782749, 1099513782750⟩,
  ⟨1099513782972, 1099513782973, 1099514382774, 1099514382775, 1099513782943, 1099513782944⟩,
  ⟨1099513782888, 1099513782889, 1099513783226, 1099513783227, 1099514380229, 1099514380230⟩
]
def w_n8_morse : List Nat := [2, 1, 2, 2, 1, 1, 2, 1, 2, 1, 2, 2, 1, 1, 2, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 1, 2, 1, 1, 2, 2, 2, 1, 1, 2, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 2, 2]

theorem t20_compiled_n8_fiftyone_disjointness_and_morse :
    w_n8_boxes.length = 51 ∧
    allPairsDisjoint w_n8_boxes = true ∧
    allBoxesAvoidSites w_n8_boxes w_n8_sites = true ∧
    countVal 1 w_n8_morse = 22 ∧ countVal 2 w_n8_morse = 29 ∧
    countVal 2 w_n8_morse - countVal 1 w_n8_morse = 8 - 1 ∧
    (8 - 1) * (8 - 1) < w_n8_boxes.length := by
  decide

def w_th_n4_lo_sites : List RationalSite3 := [
  ⟨0, 1, 0, 1, 1, 1⟩,
  ⟨0, 1, 24, 25, (-7 : Int), 25⟩,
  ⟨96, 125, (-72 : Int), 125, (-7 : Int), 25⟩,
  ⟨(-96 : Int), 125, (-72 : Int), 125, (-7 : Int), 25⟩
]
def w_th_n4_hi_sites : List RationalSite3 := [
  ⟨1, 1, 0, 1, 0, 1⟩,
  ⟨0, 1, 1, 1, 0, 1⟩,
  ⟨(-1 : Int), 1, 0, 1, 0, 1⟩,
  ⟨0, 1, (-1 : Int), 1, 0, 1⟩
]
def w_th_n4_lo_br : List SqrtBracket := [
  ⟨64, 25, 6553, 4096, 3277, 2048⟩,
  ⟨64, 25, 6553, 4096, 3277, 2048⟩,
  ⟨64, 25, 6553, 4096, 3277, 2048⟩,
  ⟨9216, 3125, 3517, 2048, 7035, 4096⟩,
  ⟨9216, 3125, 3517, 2048, 7035, 4096⟩,
  ⟨36864, 15625, 6291, 4096, 1573, 1024⟩
]
def w_th_n4_hi_br : List SqrtBracket := [
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨4, 1, 2, 1, 8193, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨4, 1, 2, 1, 8193, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩
]
def w_th_n6_lo_sites : List RationalSite3 := [
  ⟨1, 1, 0, 1, 0, 1⟩,
  ⟨(-1 : Int), 1, 0, 1, 0, 1⟩,
  ⟨0, 1, 1, 1, 0, 1⟩,
  ⟨0, 1, (-1 : Int), 1, 0, 1⟩,
  ⟨0, 1, 0, 1, 1, 1⟩,
  ⟨0, 1, 0, 1, (-1 : Int), 1⟩
]
def w_th_n6_hi_sites : List RationalSite3 := [
  ⟨1, 1, 0, 1, 0, 1⟩,
  ⟨3, 5, 4, 5, 0, 1⟩,
  ⟨(-3 : Int), 5, 4, 5, 0, 1⟩,
  ⟨(-1 : Int), 1, 0, 1, 0, 1⟩,
  ⟨(-3 : Int), 5, (-4 : Int), 5, 0, 1⟩,
  ⟨3, 5, (-4 : Int), 5, 0, 1⟩
]
def w_th_n6_lo_br : List SqrtBracket := [
  ⟨4, 1, 2, 1, 8193, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨4, 1, 2, 1, 8193, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨2, 1, 181, 128, 5793, 4096⟩,
  ⟨4, 1, 2, 1, 8193, 4096⟩
]
def w_th_n6_hi_br : List SqrtBracket := [
  ⟨4, 5, 3663, 4096, 229, 256⟩,
  ⟨16, 5, 7327, 4096, 229, 128⟩,
  ⟨4, 1, 2, 1, 8193, 4096⟩,
  ⟨16, 5, 7327, 4096, 229, 128⟩,
  ⟨4, 5, 3663, 4096, 229, 256⟩,
  ⟨36, 25, 4915, 4096, 1229, 1024⟩,
  ⟨16, 5, 7327, 4096, 229, 128⟩,
  ⟨4, 1, 2, 1, 8193, 4096⟩,
  ⟨64, 25, 6553, 4096, 3277, 2048⟩,
  ⟨4, 5, 3663, 4096, 229, 256⟩,
  ⟨64, 25, 6553, 4096, 3277, 2048⟩,
  ⟨4, 1, 2, 1, 8193, 4096⟩,
  ⟨4, 5, 3663, 4096, 229, 256⟩,
  ⟨16, 5, 7327, 4096, 229, 128⟩,
  ⟨36, 25, 4915, 4096, 1229, 1024⟩
]

theorem t21_compiled_thomson_s2q_energy_dominance :
    w_th_n4_lo_sites.all siteOnUnitSphere = true ∧
    w_th_n4_hi_sites.all siteOnUnitSphere = true ∧
    pairDistsMatch w_th_n4_lo_sites w_th_n4_lo_br = true ∧
    pairDistsMatch w_th_n4_hi_sites w_th_n4_hi_br = true ∧
    w_th_n4_lo_br.all validSqrtBracket = true ∧
    w_th_n4_hi_br.all validSqrtBracket = true ∧
    qLt (sumE2 w_th_n4_lo_br) (sumE2 w_th_n4_hi_br) = true ∧
    qLt (sumE1Hi w_th_n4_lo_br) (sumE1Lo w_th_n4_hi_br) = true ∧
    qLt ⟨9, 4⟩ (sumE2 w_th_n4_lo_br) = true ∧
    w_th_n6_lo_sites.all siteOnUnitSphere = true ∧
    w_th_n6_hi_sites.all siteOnUnitSphere = true ∧
    pairDistsMatch w_th_n6_lo_sites w_th_n6_lo_br = true ∧
    pairDistsMatch w_th_n6_hi_sites w_th_n6_hi_br = true ∧
    w_th_n6_lo_br.all validSqrtBracket = true ∧
    w_th_n6_hi_br.all validSqrtBracket = true ∧
    qLt (sumE2 w_th_n6_lo_br) (sumE2 w_th_n6_hi_br) = true ∧
    qLt (sumE1Hi w_th_n6_lo_br) (sumE1Lo w_th_n6_hi_br) = true := by
  decide

def w_stereo_p0a : RationalSite3 := ⟨36, 61, 48, 61, (-11 : Int), 61⟩
def w_stereo_p0b : RationalSite3 := ⟨(-96 : Int), 197, 504, 985, 697, 985⟩
def w_stereo_p1a : RationalSite3 := ⟨24, 25, 0, 1, (-7 : Int), 25⟩
def w_stereo_p1b : RationalSite3 := ⟨0, 1, 1, 1, 0, 1⟩
def w_cayley_q0 : RationalSite4 := ⟨(-37 : Int), 325, 144, 325, (-192 : Int), 325, 216, 325⟩
def w_cayley_q1 : RationalSite4 := ⟨11, 61, 36, 61, 48, 61, 0, 1⟩

/-- Compiled rational stereographic S^2(Q), Cayley S^3(Q), and 4-coordinate equatorial bridge theorem. -/
theorem t22_compiled_stereographic_s2q_and_cayley_s3q_bridges :
    stereographicMatches 1 2 2 3 w_stereo_p0a = true ∧
    stereographicMatches (-5 : Int) 3 7 4 w_stereo_p0b = true ∧
    siteOnUnitSphere w_stereo_p0a = true ∧
    siteOnUnitSphere w_stereo_p0b = true ∧
    checkPairDistSq w_stereo_p0a w_stereo_p0b 24336 12017 = true ∧
    stereographicMatches 3 4 0 1 w_stereo_p1a = true ∧
    stereographicMatches 0 1 1 1 w_stereo_p1b = true ∧
    siteOnUnitSphere w_stereo_p1a = true ∧
    siteOnUnitSphere w_stereo_p1b = true ∧
    checkPairDistSq w_stereo_p1a w_stereo_p1b 2 1 = true ∧
    cayleyS3Matches 1 2 (-2 : Int) 3 3 4 w_cayley_q0 = true ∧
    siteOnUnitSphere4 w_cayley_q0 = true ∧
    cayleyS3Matches 1 2 2 3 0 1 w_cayley_q1 = true ∧
    siteOnUnitSphere4 w_cayley_q1 = true ∧
    equatorialCayleyMatchesStereo w_cayley_q1 w_stereo_p0a = true := by
  decide

end PointChargesExactCert
