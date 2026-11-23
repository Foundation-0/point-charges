# Track S20a — Zero-Float Pure-$\mathbb{Z}/\mathbb{Q}$ Krawczyk & Bareiss Sylvester-Inertia Certification of Four-Charge `11`-Equilibria Configurations (`Foundation-0/point-charges#22`)

- **Track ID:** `TRACK-S20A-POINT-CHARGES-N2-EQUILIBRIA`
- **Oracle Class:** `exact-dyadic-interval+bareiss-sylvester` (Python 3 standard library only: `int`, `fractions.Fraction`, `math.isqrt`; **zero IEEE-754 floats**, **zero external dependencies**, **zero bare `assert`s**, **zero RNG**)
- **Law-4 Novelty Accounting:** `N1` independent zero-shared-code pure-$\mathbb{Z}/\mathbb{Q}$ certification & Law-2b bifurcation-wall proof over existing `N2` four-charge configurations (`Foundation-0/point-charges#22` by `@zero-aska` at commit `da0aaa4f44780e27eef1ca44d80f4d743592095a`, following the classical Maxwell 1873 conjecture, Gabrielov–Novikov–Shapiro 2007 formalization, Arathoon et al. Jul 2026 `arXiv:2607.27197` $n=5$ counterexample, `arXiv:2607.28785`, and Shiva Kintali Aug 2026 proving $M(4) \ge 11$).

---

## 1. OWN_TERMS_FIRST — Problem Statement & Scope Boundary

For $n$ positive point charges $q_1, \dots, q_n \in \mathbb{Q}_{>0}$ at distinct sites $a_1, \dots, a_n \in \mathbb{Q}^3$, the electrostatic Coulomb potential on $\mathbb{R}^3 \setminus \{a_1, \dots, a_n\}$ is
$$V(x) = \sum_{j=1}^n \frac{q_j}{\|x - a_j\|}, \qquad \nabla V(x) = -\sum_{j=1}^n q_j \frac{x - a_j}{\|x - a_j\|^3}, \qquad H(x) = \nabla^2 V(x) = \sum_{j=1}^n q_j \left(\frac{3(x - a_j)(x - a_j)^\top}{\|x - a_j\|^5} - \frac{I_3}{\|x - a_j\|^3}\right).$$
Maxwell's 1873 conjecture (formalized by Gabrielov–Novikov–Shapiro 2007) asserted at most $(n - 1)^2$ non-degenerate critical points ($\nabla V(x) = 0$, $\det H(x) \ne 0$), which equals $9$ for $n = 4$.

In `Foundation-0/point-charges#22` (commit `da0aaa4f44780e27eef1ca44d80f4d743592095a`, tree `92121f912928589a0e9930631c09d52a6d22d01a`), `@zero-aska` constructed three explicit rational four-charge configurations (`n4-eleven-short`, `n4-eleven-flat`, `n4-eleven-jet`) and checked them with C-library ball arithmetic (`python-flint`/`Arb`, `NumPy`/`SciPy`/`mpmath`).

**What Track S20a independently proves (zero shared code, zero floats, stdlib only):**
1. **Air-tight lower bound $N_{\mathrm{eq}} \ge 11 > 9$ on ALL three `#22` configurations (`33/33` equilibria):** For each configuration, `11` pairwise-disjoint closed rational boxes $X_1, \dots, X_{11} \subset \mathbb{Q}^3$ avoid all charge sites, strictly contract under the exact Krawczyk operator ($K(X_k) \subset \mathrm{int}(X_k)$ and $\rho_k = \|I_3 - Y_k H(X_k)\|_\infty < 1$), and carry certified non-degenerate Sylvester inertia (`7` of Morse index $\lambda = 2$, inertia $(1,2,0)$; `4` of Morse index $\lambda = 1$, inertia $(2,1,0)$; Earnshaw $m_0 = m_3 = 0$; Poincaré–Hopf alternating sum $\sum_{k=1}^{11} (-1)^{\lambda_k} = m_2 - m_1 = 3 = n - 1$).
2. **Exact Law-2b `11 -> 9` Saddle-Node Bifurcation Wall (`bracket width = 1/10^9`):** Perturbing a single charge $q_2(\tau) = 1 + \tau$ on `n4-eleven-short`, we certify **both sides** of the fold singularity ($\gamma = 0$): `11` equilibria at $\tau_{\mathrm{lo}} = 4239/125000000$ (`3.3912e-5`) with colliding opposite-index roots `(6, 7)` strictly inside a rational box $B_{\mathrm{annih}}$, and `9` surviving equilibria (`m_1 = 3, m_2 = 6, m_2 - m_1 = 3`) strictly outside $B_{\mathrm{annih}}$ at $\tau_{\mathrm{hi}} = 33913/1000000000$ (`3.3913e-5`) paired with a certified `48`-leaf centered-form zero-free exclusion cover proving $0 \notin \nabla V_{\tau_{\mathrm{hi}}}(B_{\mathrm{annih}})$.

**Explicit Scope Boundary (Law 4 Honesty):**
- Local Krawczyk certification proves **at least `11`** pairwise-disjoint non-degenerate equilibria per configuration — which is **mathematically sufficient by itself to refute Maxwell's upper bound $\le 9$**.
- It does **not** perform global completeness exclusion over all of $\mathbb{R}^3 \setminus \{a_1, \dots, a_4\}$ (i.e., ruling out a 12th equilibrium elsewhere in $\mathbb{R}^3$; that global upper bound was handled separately in `#22` via `Arb` subdivision). Likewise, the post-fold zero-free cover proves annihilation specifically on $B_{\mathrm{annih}}$, not globally on $\mathbb{R}^3$.

---

## 2. Pure-$\mathbb{Z}/\mathbb{Q}$ Mathematical Proof Architecture

1. **Outward-Rounded Exact Dyadic Interval Arithmetic (`PREC_BITS = 160`, $S = 2^{160}$):**
   - Every interval is stored as an exact integer pair `(lo, hi)` denoting $[lo / 2^{160}, hi / 2^{160}] \subset \mathbb{Q}$.
   - Rational conversion uses exact floor/ceiling shifts (`(p << P) // q` and `-(((-p) << P) // q)`), multiplication uses exact four-corner `>> P` / `-((-max) >> P)`, and $r^{-3} = u^{-3/2}$, $r^{-5} = u^{-5/2}$ use exact integer `math.isqrt` with explicit site-avoidance guard $u_{\mathrm{lo}} > 0$ (`SITE_OVERLAP`).
2. **Exact Krawczyk Existence, Uniqueness & Box-Wide Invertibility:**
   - For each box $X \subset \mathbb{Q}^3$ with midpoint $m$, preconditioner $Y \in \mathbb{Q}^{3\times 3}$ (quantized to $2^{-80}$ from the exact inverse of midpoint symmetric matrix $\tilde{H} \in H(X)$), and $C = I_3 - Y H(X)$, we verify both $\rho = \|C\|_\infty < 1$ and strict interior inclusion $X_{\mathrm{lo}} < K(X)_{\mathrm{lo}} \le K(X)_{\mathrm{hi}} < X_{\mathrm{hi}}$ in exact `Fraction` arithmetic.
   - By Banach's Fixed-Point Theorem and Neumann series invertibility ($\|I_3 - Y A\|_\infty \le \rho < 1$ for every $A \in H(X)$), $X$ contains a **unique** zero $x^*$ of $\nabla V$, and **every** matrix $A \in H(X)$ is invertible.
3. **Exact Bareiss + Lagrange Sylvester-Inertia Transfer (`S^T \tilde{H} S = D`):**
   - Because both $\tilde{H}$ and $H(x^*)$ lie in the convex symmetric interval box $H(X)$, the straight segment $(1-t)\tilde{H} + t H(x^*) \subset H(X)$ consists exclusively of invertible symmetric matrices; hence $\mathrm{inertia}(H(x^*)) = \mathrm{inertia}(\tilde{H})$.
   - `sylvester_inertia_exact` performs exact rational Lagrange congruence diagonalization over $\mathbb{Q}$, constructs the explicit transform $S \in \mathrm{GL}_3(\mathbb{Q})$, **self-verifies** $\det(S) \ne 0$ via fraction-free integer Bareiss elimination (`bareiss_det_exact`) and exact identity $S^\top \tilde{H} S = \mathrm{diag}(d_1, d_2, d_3)$, and cross-checks against the direct $3\times 3$ interval determinant `iv_det3(H(X))` strictly excluding `0` plus harmonic trace enclosure $0 \in \mathrm{tr}(H(X))$.
4. **Centered-Form Mean-Value Zero-Free Exclusion (`certify_zero_free_box`):**
   - For any subbox $B$ with midpoint $m$ and fold normal $v \in \mathbb{Q}^3$, MVT yields $\nabla V(x) \in \nabla V(m) + H(B)(B - m)$ and $v^\top \nabla V(x) \in v^\top \nabla V(m) + (v^\top H(B))(B - m)$. Excluding `0` in any coordinate or along $v$ (which suppresses linear overestimation near $\ker H(x_*)$) certifies $0 \notin \nabla V(B)$.

---

## 3. Certified Results Summary (`point_charges_n2_certificate.json`)

| Configuration | $n$ | Maxwell $(n-1)^2$ | Certified Equilibria ($\ge$) | Exceeds Maxwell? | $(m_0, m_1, m_2, m_3)$ | $m_2 - m_1$ | Box Half-Width | $\max \rho_{\mathrm{Kraw}}$ |
|---|---:|---:|---:|:---:|---|---:|---|---|
| `n2-analytic-control` | `2` | `1` | `1` | No (exact control) | `(0, 0, 1, 0)` | `1` | `1/2^20` | `~1.431e-05` |
| `n4-eleven-short` (`#22`) | `4` | `9` | **`11`** | **Yes (`11 > 9`)** | `(0, 4, 7, 0)` | `3` | `1/2^40` | `~5.126e-08` |
| `n4-eleven-flat` (`#22`) | `4` | `9` | **`11`** | **Yes (`11 > 9`)** | `(0, 4, 7, 0)` | `3` | `1/2^40` | `~1.851e-07` |
| `n4-eleven-jet` (`#22`) | `4` | `9` | **`11`** | **Yes (`11 > 9`)** | `(0, 4, 7, 0)` | `3` | `1/2^40` | `~5.412e-06` |

### Law-2b Exact Saddle-Node Fold Wall (`n4-eleven-short`, $q_2(\tau) = 1 + \tau$)

| State | $\tau$ | Certified Equilibria | $(m_1, m_2)$ | $m_2 - m_1$ | $\max \rho_{\mathrm{Kraw}}$ | Annihilation Region $B_{\mathrm{annih}}$ Proof |
|---|---|---:|---|---:|---|---|
| Pre-fold | `4239/125000000` (`3.3912e-5`) | **`11`** | `(4, 7)` | `3` | `~2.885e-06` | Colliding pair `(6, 7)` (indices `1` & `2`, $\ell_\infty$ gap `~6.570e-06`) **strictly inside** $B_{\mathrm{annih}}$ |
| Post-fold | `33913/1000000000` (`3.3913e-5`) | **`9`** | `(3, 6)` | `3` | `~7.154e-08` | Surviving `9` boxes **strictly disjoint** from $B_{\mathrm{annih}}$; **`48`-leaf** centered-form cover proves $0 \notin \nabla V_{\tau_{\mathrm{hi}}}(B_{\mathrm{annih}})$ |

- **Certified fold bracket width:** $\Delta\tau = \tau_{\mathrm{hi}} - \tau_{\mathrm{lo}} = 1 / 10^9$ (`10^-9`), preserving the Poincaré–Hopf invariant $m_2 - m_1 = 3$ identically across the annihilation of the `(index 1, index 2)` pair.
- **Source provenance hashes (`Foundation-0/point-charges#22 @ da0aaa4f`):**
  - `certificates/n4-eleven-short.json`: `921842862facd0a837e3400cf8dd1fffe6e57efdafa9ff1c043fd013ba33efcc`
  - `certificates/n4-eleven-flat.json`: `a17b38e01f23a8535b076ce06d604cba1554f80713ab0bc5beec8e4c66cf6ebe`
  - `certificates/n4-eleven-jet.json`: `0d9e335eeeda62d4302da95793cb09c77c0bd7cffc9e4ac359953ca5aeeb32e8`

---

## 4. Verification, Mutation, Zone-B Hardening & Frozen Manifest

- **4-Mode Isolated `-I` Parity (`-I -B`, `-I -O -B`, `-I -B -W error`, `PYTHONHASHSEED=42 -I -B`):**
  - `point_charges_n2_audit.py --standalone --json` matches `point_charges_n2_certificate.json` byte-for-byte across all `4` modes (`payload_sha256: 11be56ab97dc4b09de1ab52a23a6d800c3f49f7feba8f1fc508a9a9abf34d326`, `certificate_sha256: df784478d34ea9b3e29b0af97843ed75f2963598f674c8320e96b9d7ee71416e`).
  - `point_charges_n2_audit.py --check point_charges_n2_certificate.json`: `valid: true` × 4.
  - `test_point_charges_n2_r5.py`: **`12/12 OK`** × 4 (`0` bare `assert`s, `0` skips; includes R2 public-API Zone-B regressions closing `S20A-BUG-01` & `S20A-BUG-02`).
- **Live `aletheia v2.0.0 mutate --max-mutants 0 --work-root . --add-cwd` (`-I -O -B -W error`):**
  - `point_charges_n2_oracle.py`: **`72/72` (`100.0%`, `survived=0`, `truncated=False`) KILLED**.
  - `point_charges_n2_audit.py`: **`36/36` (`100.0%`, `survived=0`, `truncated=False`) KILLED**.
  - **Combined:** **`108/108` (`100.0%`) KILLED**.
- **Exhaustive AST `BoolOp` (`and` <-> `or`) Sweep (`-I -O -B -W error`):**
  - `point_charges_n2_oracle.py`: **`23/23` (`100.0%`, `survived=[]`) KILLED**.
  - `point_charges_n2_audit.py`: **`5/5` (`100.0%`, `survived=[]`) KILLED**.
  - **Combined:** **`28/28` (`100.0%`) KILLED**.
- **Live `Track S19c` Zone-B Adversarial Fuzzer (`zone_b_fuzzer.py`, all 16 classes `ZB01`..`ZB16`):**
  - CLI surface (`--target-cli point_charges_n2_audit.py`): **`23/23` `REJECTED_CLEAN`** (`bypasses=0, crashes=0, passed=true`).
  - In-memory surface (`verify_certificate_dict`): **`18/18` `REJECTED_CLEAN`** (`bypasses=0, crashes=0, passed=true`).
  - **Combined:** **`41/41` (`100.0%`) `REJECTED_CLEAN`**.
- **`aletheia v2.0.0 audit-check . --strict-ancestors --safety-scan --scripts point_charges_n2_audit.py test_point_charges_n2_r5.py`:** `PASS` (`5` files, `2/2` R4, `0` bare asserts, `0` CRLF, `0` pycache, `0` path/secret leaks, `0` safety findings).

### Frozen File Manifest (R2; excluding self-hash of this report)

| File | Lines | Bytes | SHA-256 |
|---|---:|---:|---|
| `point_charges_n2_oracle.py` | `951` | `47,216` | `87b55d1b692b4b3ab64d7e923beb51f405dd42519a5c771345e33c3150b3bf96` |
| `point_charges_n2_audit.py` | `289` | `11,957` | `d8d3b9e7b0f26ac6df46040ed2bd5eb79b82d6b99c6351e33e20070ad445d0ce` |
| `test_point_charges_n2_r5.py` | `723` | `41,803` | `379f144a969005761dad3b61e843461957116e0cfd3bc8970792df49f01a6def` |
| `point_charges_n2_certificate.json` | `1` | `4,871` | `2f204756a7d5ff8ee99de05f69430c30031d95a9783544985841a1c5d52f41a1` |
| `POINT_CHARGES_N2_REPORT.md` | — | — | Verified in handoff receipt |
