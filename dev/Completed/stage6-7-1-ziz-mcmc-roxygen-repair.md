# Stage 6.7.1 - ZIZ MCMC roxygen and fitted-shape repair

## Purpose

Repair the attempted Stage 6.7 migration without changing its statistical algorithm.

## Changes

- Register the internal Stage 6 model and posterior-engine S3 methods with roxygen `@exportS3Method` tags so `devtools::document()` can generate valid S3 registrations without exposing the helpers as user-facing functions.
- Preserve the historical `psFit$fitted` contract for ZIZ MCMC fits as a named numeric vector rather than a one-row matrix.
- Keep the Stage 6.7 test-suite cleanup: test filenames and test descriptions use stable behavioural names rather than development stage numbers.
- Carry forward the complete intended Stage 6.7 ZIZ MCMC migration because Stage 6.7 failed before commit.

## Validation

Run the full fitPS validation workflow. In particular, roxygen regeneration must complete without S3-registration diagnostics, and the deterministic ZIZ MCMC regression tests must pass for both P and S data.
