# Pure-$\mathbb{Z}/\mathbb{Q}$ Krawczyk & Bareiss Sylvester-Inertia Certification of Four-Charge `11`-Equilibria Configurations (`AUDIT-120AR`)

Independent, zero-dependency, zero-float Python 3 standard-library (`int`, `fractions.Fraction`, `math.isqrt`) certification study for the three four-charge `11`-equilibria configurations introduced in [`Foundation-0/point-charges#22`](https://github.com/Foundation-0/point-charges/pull/22) (`@zero-aska`, commit `da0aaa4f44780e27eef1ca44d80f4d743592095a`).

## Novelty Accounting & Prior Art

- **Historical lineage:** Maxwell (1873) conjectured at most $(n - 1)^2$ non-degenerate electrostatic equilibria for $n$ positive point charges; Gabrielov–Novikov–Shapiro (2007) formalized the problem; Arathoon–Ball–Kvalheim (Jul 2026, `arXiv:2607.27197`, `arXiv:2607.28785`) disproved it for $n = 5$; Shiva Kintali (Aug 2026) proved $M(4) \ge 11$; `@zero-aska` (`Foundation-0/point-charges#22`) constructed three explicit rational four-charge configurations (`n4-eleven-short`, `n4-eleven-flat`, `n4-eleven-jet`) and certified them with C-library `python-flint`/`Arb` ball arithmetic.
- **Contribution of this study (Law-4 `N1` certification over `#22`'s `N2` configurations):**
  1. Independently certifies all **`33/33`** equilibria across `n4-eleven-short`, `n4-eleven-flat`, and `n4-eleven-jet` (`11` pairwise-disjoint closed rational boxes per configuration, plus the analytic `n2-analytic-control`) using **zero shared code**, **zero IEEE-754 floats**, and **zero third-party libraries** — proving $N_{\mathrm{eq}} \ge 11 > 9$, Morse-index profile `(m_0, m_1, m_2, m_3) = (0, 4, 7, 0)`, and Poincaré–Hopf sum $m_2 - m_1 = 3 = n - 1$ via `160`-bit outward-rounded exact dyadic intervals, exact Krawczyk contraction, and exact Lagrange congruence ($S^\top \tilde{H} S = D$, self-verified with fraction-free integer Bareiss determinant $\det S \ne 0$).
  2. Certifies the exact **`11 -> 9` saddle-node bifurcation wall** on `n4-eleven-short` under $q_2(\tau) = 1 + \tau$ within bracket width $\Delta\tau = 1/10^9$: `11` equilibria at $\tau_{\mathrm{lo}} = 4239/125000000$ with colliding opposite-index roots `(6, 7)` strictly inside rational box $B_{\mathrm{annih}}$, and `9` surviving equilibria (`m_1 = 3, m_2 = 6, m_2 - m_1 = 3`) strictly outside $B_{\mathrm{annih}}$ at $\tau_{\mathrm{hi}} = 33913/1000000000$ paired with a certified `48`-leaf centered-form mean-value cover proving $0 \notin \nabla V_{\tau_{\mathrm{hi}}}(B_{\mathrm{annih}})$.

## Scope Boundary

- **Lower bound vs. global completeness:** Local Krawczyk certification proves **at least `11`** distinct non-degenerate equilibria per configuration, which is mathematically sufficient on its own to refute Maxwell's $(4 - 1)^2 = 9$ upper bound. It does **not** perform global completeness exclusion over all of $\mathbb{R}^3 \setminus \{a_1, \dots, a_4\}$ (handled separately in `#22` via `Arb` global subdivision). Likewise, the post-fold zero-free exclusion proof certifies annihilation specifically on $B_{\mathrm{annih}}$.
- **CI isolation:** This directory is intentionally placed under `studies/` rather than `certificates/` so the repository's `.github/workflows/verify.yml` (`python verify_arb.py certificates/*.json`) remains completely untouched.

## Reproduction (Zero External Dependencies)

Run from inside `studies/pure-q-krawczyk-n4-eleven-audit120a/` with any standard Python 3.11+ interpreter (verified across all four isolated `-I` modes: `-I -B`, `-I -O -B`, `-I -B -W error`, and `PYTHONHASHSEED=42 -I -B`):

```bash
# 1. Regenerate canonical certificate from scratch and compare byte-for-byte
python3 -I -B point_charges_n2_audit.py --standalone --json | diff -u point_charges_n2_certificate.json -

# 2. Verify frozen certificate independently (exits 0 with valid: true)
python3 -I -B point_charges_n2_audit.py --check point_charges_n2_certificate.json

# 3. Run the 12-test R5 regression, negative-control & Zone-B hardening suite
python3 -I -O -B -W error test_point_charges_n2_r5.py
```

## Frozen Study File Manifest

All five core study files are bit-identical copies of the frozen `Track S20a` R2 (`AUDIT-120AR`) verification artifacts; full mathematical details and mutation/fuzzing results are documented in [POINT_CHARGES_N2_REPORT.md](POINT_CHARGES_N2_REPORT.md).

| File | Lines | Bytes | SHA-256 |
|---|---:|---:|---|
| `point_charges_n2_oracle.py` | `951` | `47,216` | `87b55d1b692b4b3ab64d7e923beb51f405dd42519a5c771345e33c3150b3bf96` |
| `point_charges_n2_audit.py` | `289` | `11,957` | `d8d3b9e7b0f26ac6df46040ed2bd5eb79b82d6b99c6351e33e20070ad445d0ce` |
| `test_point_charges_n2_r5.py` | `723` | `41,803` | `379f144a969005761dad3b61e843461957116e0cfd3bc8970792df49f01a6def` |
| `point_charges_n2_certificate.json` | `1` | `4,871` | `2f204756a7d5ff8ee99de05f69430c30031d95a9783544985841a1c5d52f41a1` |
| `POINT_CHARGES_N2_REPORT.md` | `95` | `11,718` | `9d3cbd1463738b29ca285792bf7eaeccbd9b033f0be0958924afb17eec38e8be` |
