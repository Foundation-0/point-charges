# Track S21d: Pure-$\mathbb{Z}/\mathbb{Q}$ Krawczyk & Bareiss Sylvester-Inertia Certification of Positive Point-Charge Equilibria at the $n=7$ (`44`, `#12`) and $n=8$ (`51 > 49`, `#5`) Maxwell Frontier

**Track ID:** `TRACK-S21D` | **Oracle Class:** `exact` (pure stdlib $\mathbb{Z}/\mathbb{Q}$, zero floats, zero external dependencies) | **Novelty:** `N1` (first zero-float, stdlib-only pure-$\mathbb{Z}/\mathbb{Q}$ 384-bit Krawczyk + self-verified Lagrange/Bareiss Sylvester-inertia machine-verified certificate resolving `Foundation-0/point-charges#12` (`n=7 -> 44`) and `#5` (`n=6 -> 31 > 29`, `n=8 -> 51 > 49`), building on the $n \ge 4$ Maxwell counterexamples of Philip Arathoon, Gavin Ball, and Matthew D. Kvalheim (ABK, Jul 2026, `arXiv:2607.27197`) and Shiva Kintali (Aug 2026))

---

## 1. Executive Summary (Own Terms First)

For $n$ positive point charges $q_1, \dots, q_n > 0$ located at distinct sites $a_1, \dots, a_n \in \mathbb{R}^3$, the electrostatic potential is
$$V(x) = \sum_{j=1}^{n} \frac{q_j}{\|x - a_j\|_2}, \qquad x \in \mathbb{R}^3 \setminus \{a_1, \dots, a_n\},$$
with gradient $\nabla V(x) = -\sum_{j=1}^n q_j \frac{x - a_j}{\|x - a_j\|_2^3}$ and harmonic Hessian $H(x) = \nabla^2 V(x)$ ($\operatorname{tr} H(x) \equiv 0$). Maxwell (1873, *A Treatise on Electricity and Magnetism*, Vol. 1, §113) conjectured that $n$ point charges can produce at most $(n-1)^2$ isolated equilibrium points ($\nabla V(x) = 0$). Gabrielov, Novikov, and Shapiro (2007) established finiteness bounds, **Philip Arathoon, Gavin Ball, and Matthew D. Kvalheim** (ABK, Jul 2026, `arXiv:2607.27197`, *"The Maxwell Conjecture is False"*) refuted Maxwell's bound via nested axial pairs ($n=5 \to 24 > 16$, $n=7 \to 44 > 36$, etc.), and **Shiva Kintali** (Aug 2026) constructed counterexamples including $n=4 \to 11 > 9$ and $n=8 \to N \ge 51 > 49$. Within the open `Foundation-0/point-charges` repository (`main @ a2fd159b`, Draft `#22 @ da0aaa4f`), however, two concrete open issues remained unresolved:
- **Issue `#12` ($n=7$, `C_7`, 44 equilibria):** The repository's float64 certifier `certify.py` crashed at 53 GB / 60M boxes without certifying `certificates/n7-abk-second-pair.json`.
- **Issue `#5` ($n=8$, beat $(8-1)^2 = 49$):** The best certified configuration in the repository tied $49$ (`C_7` plus one charge at the symmetric origin, `44 + 5 = 49`), and the best certified $n=6$ configuration had `29` equilibria.

This study delivers a self-contained, zero-float, stdlib-only **384-bit outward-rounded pure-$\mathbb{Z}/\mathbb{Q}$ dyadic Krawczyk certifier** paired with **explicit non-singular preconditioner checking ($\det(Y) \ne 0$)** and **self-verified Lagrange/Bareiss Sylvester-inertia diagonalization** ($S^\top \tilde{H} S = \operatorname{diag}(d)$ with exact fraction-free $\det(S) \ne 0$), resolving both issues:

1. **Resolution of `Foundation-0/point-charges#12` ($n=7$, `C_7`, 44 Equilibria):** Certifies all **`44/44`** pairwise-disjoint non-degenerate equilibria of `certificates/n7-abk-second-pair.json` (`m_1 = 19, m_2 = 25, m_2 - m_1 = 6 = n - 1`, `44 > 36`, $\max_k \rho_k \approx 9.195 \times 10^{-8}$) in **~1.2 seconds** with zero floating-point operations, and isolates why float64 `certify.py` stalls.
2. **Resolution of `Foundation-0/point-charges#5` ($n=8$, `51 > 49`) & $n=6$ (`31 > 29`):** Using an **ABK + Axial Milnor-8 Flat-Point Unfolding** construction that improves the ABK odd-to-even charge addition from `+5` (`9 + 20m`) to **`+7` (`11 + 20m`)** for every $m \ge 0$, we construct explicit exact-rational positive-charge configurations and certify:
   - **`n = 6` (`n6-abk-milnor8-flat`):** **`31` pairwise-disjoint non-degenerate equilibria** (`m_1 = 13, m_2 = 18, m_2 - m_1 = 5`, beats the repo's `29` by `+2` and Maxwell `25` by `+6`, $\max_k \rho_k \approx 6.250 \times 10^{-9}$, 64-bit dyadic sites/charges).
   - **`n = 8` (`n8-abk-milnor8-flat`):** **`51` pairwise-disjoint non-degenerate equilibria** (`m_1 = 22, m_2 = 29, m_2 - m_1 = 7`, **strictly beats Maxwell's bound `(8-1)^2 = 49` by `+2`**, $\max_k \rho_k \approx 9.895 \times 10^{-12}$, 96-bit dyadic sites/charges).

---

## 2. Certified Configurations Summary (`127` Total Pairwise-Disjoint Non-Degenerate Equilibria)

All four exact-rational configurations below are certified deterministically from scratch by `point_charges_n7_n8_oracle.py` and sealed in `point_charges_n7_n8_certificate.json`:

| Configuration | $n$ | Maxwell $(n-1)^2$ | Certified Equilibria | Beats Maxwell? | Morse $(m_0, m_1, m_2, m_3)$ | $m_2 - m_1$ | Box Half-Width | $\max_k \rho_k$ (Krawczyk) |
|---|---:|---:|---:|---|---|---:|---:|---:|
| `n2-analytic-control` | 2 | 1 | **1** | Control (`1 = 1`) | `(0, 0, 1, 0)` | `1 = n-1` | $2^{-20}$ | $\approx 1.431 \times 10^{-5}$ |
| `n7-abk-second-pair` (`#12`) | 7 | 36 | **44** | **YES (`44 > 36`)** | `(0, 19, 25, 0)` | `6 = n-1` | $2^{-50}$ | $\approx 9.195 \times 10^{-8}$ |
| `n6-abk-milnor8-flat` (`#5` family) | 6 | 25 | **31** | **YES (`31 > 29 > 25`)** | `(0, 13, 18, 0)` | `5 = n-1` | $2^{-55}$ | $\approx 6.250 \times 10^{-9}$ |
| `n8-abk-milnor8-flat` (`#5` target) | 8 | 49 | **51** | **YES (`51 > 49`)** | `(0, 22, 29, 0)` | `7 = n-1` | $2^{-75}$ | $\approx 9.895 \times 10^{-12}$ |

Every certified box $X_k = c_k + [-h_k, h_k]^3 \subset \mathbb{R}^3$ satisfies all seven exact pure-$\mathbb{Z}/\mathbb{Q}$ proofs simultaneously:
1. All charges strictly positive ($q_j > 0$, all $q_j \in \mathbb{Q}$) and all charge sites pairwise distinct and strictly outside every closed box $X_k$.
2. Exact pairwise strict disjointness over $\mathbb{Q}$ across all $\binom{N}{2}$ box pairs ($\exists d \in \{x,y,z\}: u_{i,d} < \ell_{j,d} \lor u_{j,d} < \ell_{i,d}$).
3. Exact non-singular rational preconditioner $\det(Y_k) \ne 0$ (verified directly in $\mathbb{Q}$ before forming the Krawczyk operator).
4. Strict Krawczyk contraction $\rho_k = \|I_3 - Y_k H(X_k)\|_\infty < 1$ over 384-bit outward-rounded dyadic enclosures.
5. Strict interior inclusion $K(X_k) = c_k - Y_k \nabla V(c_k) + (I_3 - Y_k H(X_k))(X_k - c_k) \subset \operatorname{int}(X_k)$ on all three axes (proving existence and uniqueness of a single equilibrium $x_k^* \in \operatorname{int}(X_k)$ by Banach's fixed-point / Krawczyk theorem).
6. Harmonic trace consistency $0 \in \operatorname{tr}(H(X_k))$ and box-wide non-degeneracy $0 \notin \det(H(X_k))$ (guaranteeing constant Sylvester inertia throughout the connected box $X_k$).
7. Exact midpoint Sylvester inertia $(n_+, n_-, 0)$ via Lagrange congruence $S_k^\top \tilde{H}_k S_k = \operatorname{diag}(d)$ with exact integer Bareiss $\det(S_k) \ne 0$ and $\tilde{H}_k \in H(X_k)$ entrywise.

---

## 3. Mathematical Construction: The ABK + Axial Milnor-8 Flat-Point Unfolding (`11 + 20m`)

### Why Repository Even-$n$ Constructions Stalled at $9 + 20m$ (`29` for $n=6$, `49` for $n=8$)

The ABK nested-pair family $C_{3+2m}$ consists of a regular equilateral triangle of three unit charges at $(3,0,0), (0,3,0), (0,0,3)$ centered at $c = (1,1,1)$ along the $C_{3v}$ symmetry axis $e_z = \frac{1}{\sqrt{3}}(1,1,1)$, plus $m$ symmetric axial pairs $c \pm \eta_r e_z$ with charges $q^{(r)}$. Each pair splits the central $D_{3h}$-symmetric index-1 equilibrium into a new central index-1 point surrounded by 20 new equilibria (giving $4 + 20m$: `4, 24, 44, ...` for $n = 3, 5, 7, \dots$).

Prior `Foundation-0/point-charges` constructions (`NOTE.md` Prop. 2(b) and Draft `#22`) added one more positive charge directly *at* the central index-1 point $c$ of $C_{3+2m}$. Because the central Hessian has signature $(+, +, -)$ with negative axial eigenvalue $-\lambda_z < 0$ and positive transverse eigenvalues $+\frac{\lambda_z}{2}, +\frac{\lambda_z}{2}$, placing a tiny positive charge at $c$ creates only **6** new inner equilibria in place of the destroyed central equilibrium — a net gain of **`+5`**, yielding $M(4 + 2m) \ge 9 + 20m$ (`9` for $n=4$, `29` for $n=6$, `49 = (8-1)^2` tie for $n=8$).

### Displaced Axial Milnor-8 Flat Point (`+7` Net Gain $\implies M(4+2m) \ge 11 + 20m$)

Near the central index-1 equilibrium $c$ of $C_{3+2m}$, the potential has $D_{3h}$ Taylor expansion in local cylindrical coordinates $(u, v, z)$ aligned with the triangle ($e_z = \frac{1}{\sqrt{3}}(1,1,1)$, $e_u = \frac{1}{\sqrt{6}}(2,-1,-1)$, $e_v = \frac{1}{\sqrt{2}}(0,1,-1)$):
$$V_{C_{3+2m}}(u, v, z) = V_0 + \alpha\left(\tfrac{1}{2}(u^2 + v^2) - z^2\right) + \beta (u^3 - 3 u v^2) + O(r^4), \qquad \alpha, \beta > 0.$$
Placing the $(4 + 2m)$-th positive charge $q_*$ **displaced along the $C_{3v}$ axis** at $a_{\text{new}} = c + \delta e_z$ ($\delta > 0$) gives along the axis ($u = v = 0$) $f(z) = -\alpha z^2 + \frac{q_*}{\delta - z} + O(z^4)$, which satisfies $f'(z_*) = f''(z_*) = 0$ simultaneously at:
$$z_* \approx \frac{\delta}{3}, \qquad q_* \approx \frac{16}{27}\alpha \delta^3 > 0.$$
By harmonicity ($\Delta V = 0$) and $C_{3v}$ rotational symmetry around $e_z$, vanishing of the axial second derivative $V_{zz}(0,0,z_*) = 0$ forces the transverse second derivatives $V_{uu} = V_{vv} = -\frac{1}{2} V_{zz} = 0$ simultaneously. Thus $x_{\text{flat}} = c + z_* e_z$ is an **exact degenerate flat point** ($\nabla V = 0$, $H = 0$, Milnor number $\mu = 8$) whose leading cubic harmonic Taylor jet is:
$$V(x_{\text{flat}} + (u, v, w)) = \text{const} + \frac{T_{zzz}}{6}\left(w^3 - \tfrac{3}{2} w (u^2 + v^2)\right) + \beta (u^3 - 3 u v^2) + O(r^4),$$
where $T_{zzz} \approx -\frac{9\alpha}{\delta}$ while $\beta$ is fixed by the equilateral triangle. Hence the dimensionless cubic harmonic shape ratio
$$\left|\frac{T_{zzz}}{6\beta}\right| \approx \frac{3\alpha}{2\beta \delta}$$
is continuously tunable across $(0, \infty)$ via $\delta$. Tuning $\delta$ to match the all-real Milnor-8 unfolding ratio of the `n=4` 11-equilibrium flat point (`n4-eleven-flat`, $|T_{zzz} / T_{xxx}| \approx 1.27761707$), and applying tiny $O(10^{-6}\text{--}10^{-9})$ unfolding perturbations $(Q, \ell)$ to the outer triangle and the displaced axial charge, splits the flat point into **`8` distinct non-degenerate inner equilibria (`4` index-1, `4` index-2)** while preserving all `3 + 20m` outer equilibria of $C_{3+2m}$ (net `+7` instead of `+5`, adding `+3` to $m_1$ and `+4` to $m_2$, preserving $m_2 - m_1 = n - 1$).

- For $m = 1$ ($n = 6$, built on `C_5`): optimal $\delta \approx 8.82 \times 10^{-3}$, $q_6 \approx 7.674 \times 10^{-11}$ (`707809575 / 2^63`), yielding **`31` certified equilibria** (`m_1 = 13, m_2 = 18`).
- For $m = 2$ ($n = 8$, built on `C_7`): optimal $\delta \approx 1.053 \times 10^{-5}$, $q_8 \approx 1.548 \times 10^{-22}$ (`6134079 / 2^95`), yielding **`51` certified equilibria** (`m_1 = 22, m_2 = 29`) — **strictly exceeding Maxwell's $(8-1)^2 = 49$.**

---

## 4. Root-Cause Diagnosis of `Foundation-0/point-charges#12` (Why Float64 `certify.py` Stalls on `C_7`)

In `Foundation-0/point-charges` (`main @ a2fd159b`), `certificates/n7-abk-second-pair.json` sits at scale $\eta_2 = 1/(2048\sqrt{2}) \approx 3.45 \times 10^{-4}$ with second-pair charge $q_6 = q_7 = 134713 / 2^{60} \approx 1.168 \times 10^{-13}$. At the inner cluster around $c = (1,1,1)$, the outer equilateral triangle produces individual gradient and Hessian terms of magnitude $O(1)$ that cancel down to $O(10^{-7}\text{--}10^{-13})$ near the inner equilibria (where the smallest Hessian eigenvalue is $\sim 5 \times 10^{-7}$).

1. **Float64 Cancellation Collapse:** `certify.py` evaluates interval enclosures in IEEE-754 float64 (`53`-bit significand, $\varepsilon_{\text{mach}} \approx 2.22 \times 10^{-16}$). Subtracting $O(1)$ terms to isolate $O(10^{-13})$ inner residuals leaves at most 10 bits of relative precision — float64 outward rounding blows up the interval width of $I - Y H(X)$ past $1$ before Krawczyk contraction can trigger.
2. **Breadth-First Queue Memory Explosion:** When Krawczyk fails to contract, `certify.py` bisects boxes breadth-first while materializing `(m, n, 3)` float64 arrays across the active frontier, blowing past 60M boxes and 53 GB RAM without converging.
3. **Pure-$\mathbb{Z}/\mathbb{Q}$ 384-Bit Resolution:** Our pure-$\mathbb{Z}/\mathbb{Q}$ 384-bit ($\sim 115$ decimal digits) outward-rounded dyadic certifier completely eliminates cancellation — certifying all `44/44` equilibria of `C_7` in **~1.2 seconds** with $\max_k \rho_k \approx 9.195 \times 10^{-8}$, and likewise resolving `n = 8 -> 51` where the unfolding gradient target $\sim 10^{-20}$ lies four orders of magnitude below float64 machine epsilon.

---

## 5. Honest Scope Boundary (Law 4) & Attribution

- **What Is Rigorously Proved:** Local pure-$\mathbb{Z}/\mathbb{Q}$ Krawczyk certification proves exact **lower bounds** on the number of isolated, pairwise-disjoint, non-degenerate equilibria in $\mathbb{R}^3 \setminus \{a_j\}$:
  $$N_{\text{eq}}(C_7) \ge 44 > 36, \qquad N_{\text{eq}}(C_6^{\text{flat}}) \ge 31 > 29 > 25, \qquad N_{\text{eq}}(C_8^{\text{flat}}) \ge 51 > 49 = (8-1)^2.$$
  Because Maxwell's conjecture asserts an *upper bound* $N_{\text{eq}} \le (n-1)^2$, proving the certified lower bound $N_{\text{eq}}(C_8^{\text{flat}}) \ge 51 > 49$ is **mathematically complete and sufficient to refute Maxwell's bound for $n = 8$**.
- **What Is Explicitly Out of Scope:** Local Krawczyk verification on 51 disjoint boxes does *not* perform a global branch-and-bound completeness exclusion over all of $\mathbb{R}^3 \setminus \{a_j\}$ (it does not rule out a hypothetical 52nd/53rd equilibrium outside the 51 boxes, though $m_2 - m_1 = 29 - 22 = 7 = n - 1$ matches the exact Morse index sum for a complete generic set of equilibria).
- **Prior-Art Attribution & `N1` Classification:**
  - J. C. Maxwell (1873), *A Treatise on Electricity and Magnetism*, Vol. 1, §113.
  - A. Gabrielov, D. Novikov, and B. Shapiro (2007), finiteness bounds for electrostatic equilibria.
  - **Philip Arathoon, Gavin Ball, and Matthew D. Kvalheim** (ABK, Jul 2026, `arXiv:2607.27197`, *"The Maxwell Conjecture is False"*) — refutation of Maxwell's conjecture via nested axial pairs ($C_{3+2m}$).
  - **Shiva Kintali** (Aug 2026) — earlier counterexamples for $n = 4$ (`11 > 9`) and $n = 8$ ($N \ge 51 > 49$).
  - `Foundation-0/point-charges` (`main @ a2fd159b` `certificates/n7-abk-second-pair.json`, Draft `#22 @ da0aaa4f` `n4-eleven-flat` & `NOTE.md`).
  - Accordingly, this study is classified **`N1`** (first zero-float, stdlib-only pure-$\mathbb{Z}/\mathbb{Q}$ 384-bit Krawczyk + Bareiss Sylvester-inertia machine-verified certificate suite resolving `Foundation-0/point-charges#12` and `#5`).

---

## 6. Verification Receipts

- **Frozen Canonical Certificate (`point_charges_n7_n8_certificate.json`):**
  - `file_sha256`: `e3042b534a6fb4ed94bf491857e3005a10b9d011af542b6710ca5ba5835d4689`
  - `payload_sha256`: `8501ea4867788279e5b4b6c3e8070309ced321e6159162a64d3002cdf9e82752`
  - `certificate_sha256`: `1350e5f5493a5b27158e24b2a2a40cfeecf9280912b15702fe20bb964644260b`
- **5-Mode Isolated `-I` Parity (`-I -B`, `-I -O -B`, `-I -B -W error`, `-I -O -B -W error`, `PYTHONHASHSEED=42 -I -B`):**
  - `point_charges_n7_n8_audit.py --standalone --json` x5: bit-identical `e3042b534a6fb4ed94bf491857e3005a10b9d011af542b6710ca5ba5835d4689`.
  - `point_charges_n7_n8_audit.py --check point_charges_n7_n8_certificate.json` (default full pure-$\mathbb{Z}/\mathbb{Q}$ 384-bit Krawczyk + Bareiss replay across all 127 equilibria) x5: `valid: true`, `full_replay: true`, `total_equilibria_certified: 127`.
  - `test_point_charges_n7_n8_r5.py` x5: `11/11 OK` (zero bare `assert`s; independent Pythagorean exact-rational Coulomb cross-checks, full Lagrange/Bareiss inertia regimes, explicit `det(Y) != 0` preconditioner guard, atomic `O_NOFOLLOW | O_NONBLOCK` fd I/O + streaming cap on `/proc/kallsyms` + `st_nlink > 1` hardlink refusal, and full tamper matrix).
- **Live Mutation & AST `BoolOp` Testing:**
  - `aletheia v2.0.0 mutate --max-mutants 0 --work-root . --add-cwd` (`-s -O -B -W error`): `point_charges_n7_n8_oracle.py` **`154/154`** + `point_charges_n7_n8_audit.py` **`5/5`** = **`159/159` (`100.0%`, `survived=0`, `truncated=False`) KILLED**.
  - Exhaustive AST `BoolOp` (`and` <-> `or`) sweep: `point_charges_n7_n8_oracle.py` **`19/19`** + `point_charges_n7_n8_audit.py` **`2/2`** = **`21/21` (`100.0%`, `surv=[]`) KILLED**.
- **Reusable 16-Vector Zone-B Adversarial Fuzzer (`zone_b_fuzzer.py`):**
  - CLI surface (`point_charges_n7_n8_audit.py --check`, default full replay): **`23/23` `REJECTED_CLEAN`** (`baseline_ok=true`, `bypasses=0`, `crashes=0`, `PASS`).
  - In-memory validator surface (`validate_certificate_payload`): **`18/18` `REJECTED_CLEAN`** (`bypasses=0`, `crashes=0`, `PASS`).
  - Combined: **`41/41` (`100.0%`) `REJECTED_CLEAN`**, plus live `/proc/kallsyms` (`FILE_TOO_LARGE`), hardlink (`HARDLINK_FORBIDDEN`), `/dev/null` (`NON_REGULAR_FILE`), symlinked-`PWD`, FIFO, and CLI `argv` rejection.
- **Static & Hygiene Gate (`aletheia v2.0.0 audit-check . --strict-ancestors --safety-scan --scripts point_charges_n7_n8_audit.py test_point_charges_n7_n8_r5.py`):**
  - `PASS` (`5` files scanned, `3` `.py`, `1` `.json`, `2/2` R4 3-mode parity, `0` bare asserts, `0` AST float literals / `float()` calls on oracle & audit CLI, `0` CRLF, `0` `__pycache__`, `0` path leaks, `0` ancestor permission findings, `0` safety findings).
