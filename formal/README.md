# Formal statement draft

The definitions distinguish the finite equilibrium set from an exact count. The Maxwell-bound definition names a historical proposed statement rather than asserting it. It does not supply a theorem for any charge configuration.

Next steps are a nondegenerate-Hessian definition, exact rational configuration encoding, and a proved interval certificate replay or the analytic lemmas it relies on. Python certificates remain outside Lean trust until that replay exists.

## Build (from repository root)

Requires [elan](https://github.com/leanprover/elan). The repo pins Lean and Mathlib via `lean-toolchain`, `lakefile.lean`, and `lake-manifest.json`.

```bash
lake exe cache get
lake env lean formal/PointCharges.lean
```

A successful check prints nothing and exits with status 0. The first `cache get` downloads prebuilt Mathlib artifacts (several hundred MB); later runs are fast when the cache is warm.
