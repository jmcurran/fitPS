# Stage 6.8 - Laplace and importance ZIZ engine migration

## Purpose

Stage 6.8 completes migration of the existing zero-inflated zeta posterior approximation methods onto the model/engine architecture introduced in Stage 6.3. It covers the ZIZ-only Laplace and importance-sampling methods. It does not add Laplace or importance fitting for the plain zeta model.

## Changes

- Add model-specific Laplace and importance fitting protocols for `zizModel`.
- Add `fitPosterior()`, `summarisePosterior()`, `posteriorDiagnostics()`, and `posteriorPointEstimate()` methods for the Laplace and importance engines.
- Preserve the existing transformed-coordinate Laplace optimisation, Hessian and delta-method covariance calculations.
- Preserve the existing Laplace-based Gaussian importance proposal, seeded sampling, normalized weights, weighted moments, and diagnostics.
- Route `fitZIDistBayesLaplace()` and `fitZIDistBayesImportance()` through the shared model/engine protocol while retaining their current `psFit` fields and `psPosterior` representation payloads.
- Keep fitted probabilities as named numeric vectors at the `psFit` boundary.
- Add deterministic regression tests for P and S data using stable behaviour-based test names and descriptions.

## Deliberate non-changes

- Plain-zeta Laplace and importance engines remain unsupported.
- `psPosterior` representation harmonisation is deferred to Stage 6.9.
- Inflation-probability dispatch and common `psFit` finalisation remain deferred to Stage 6.9.
- The underlying Laplace and importance statistical algorithms are unchanged.

## Validation expectation

This is a full package stage because R code and tests change. The runner must regenerate roxygen documentation, run strict tests, run strict package check, commit the validated state, create `stage6_8_completed.zip`, build the package, and install the exact built source archive.
