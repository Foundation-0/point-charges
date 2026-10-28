import Mathlib

/- Statement definitions only. No Maxwell bound or certificate count is proved. -/
namespace PointCharges

abbrev Space := EuclideanSpace ℝ (Fin 3)

def potential {n : ℕ} (sites : Fin n → Space) (charges : Fin n → ℝ) (x : Space) : ℝ :=
  ∑ i, charges i / dist x (sites i)

def equilibrium {n : ℕ} (sites : Fin n → Space) (charges : Fin n → ℝ) (x : Space) : Prop :=
  (∀ i, x ≠ sites i) ∧ fderiv ℝ (potential sites charges) x = 0

def positiveCharges {n : ℕ} (charges : Fin n → ℝ) : Prop :=
  ∀ i, 0 < charges i

def exactlyEquilibria {n : ℕ} (sites : Fin n → Space) (charges : Fin n → ℝ) (count : ℕ) : Prop :=
  {x | equilibrium sites charges x}.Finite ∧
    {x | equilibrium sites charges x}.ncard = count

-- This names the historical proposed bound. It is not asserted as a theorem.
def maxwellBoundFor {n : ℕ} (sites : Fin n → Space) (charges : Fin n → ℝ) : Prop :=
  ∀ count, exactlyEquilibria sites charges count → count ≤ (n - 1) ^ 2

end PointCharges
