# Stage 2.1.1 zeta shape parameterisation test plan

Stage 2.1 added target-behaviour tests directly under `tests/testthat`.
Those tests correctly exposed that the current implementation still treats zeta
`shape` as VGAM's shifted parameter in several paths, but they also prevented
the strict staged workflow from completing before the Stage 2.2 implementation
refactor.

Stage 2.1.1 therefore moves those target expectations out of the active test
suite and keeps them as the Stage 2 implementation plan.

## Target convention for Stage 2.2

- User-facing and package-facing `shape` means the standard zeta parameter
  alpha.
- `shape` must be greater than 1.
- VGAM receives `vgamShape = shape - 1` at explicit VGAM boundaries.
- `fitDist()`, `fitZIDist()`, `probfun()`, `predict.psFit()`,
  `fitted.psFit()`, `rzeta()`, and `rZIzeta()` should agree on that convention.

## Test expectations to reinstate during Stage 2.2

When the implementation refactor is ready, reinstate active tests that verify:

- `fitDist(..., start = 1)` and `fitDist(..., shape = 1)` fail clearly because
  user-facing zeta shape must be greater than 1.
- `rzeta(n = 1, shape = 1)` and `rZIzeta(n = 1, pi = 0.5, shape = 1)` fail
  clearly with the same convention before calling VGAM.
- Zeta probability calculations call VGAM with `shape - 1`, not `shape`.
- Zero-inflated zeta probability calculations use the same conversion for the
  non-zero component.
- `predict.psFit()` and `fitted.psFit()` return values consistent with the
  standard-shape convention.

## Why these tests are not active in Stage 2.1.1

The current code has not yet been refactored. Active target-behaviour tests are
therefore expected to fail under strict validation. Keeping the expectations in
`dev/` preserves the audit result without making the package unbuildable during
Stage 2.1.1.
