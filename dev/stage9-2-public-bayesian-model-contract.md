# Stage 9.2 public Bayesian model contract

## Purpose

Stage 9.2 implements the public model-level mathematics contract identified by the Stage 9.1 Bayesian extensibility audit. It does not yet migrate any posterior engine to generic external-model fitting.

## Public extension points

The stage adds five exported S3 generics:

- `modelLogPrior()` evaluates the complete natural-scale model prior and leaves the prior representation under model control.
- `modelBayesControl()` supplies model-owned Bayesian controls such as named natural-scale starting values.
- `modelToWorking()` maps natural parameters to engine working coordinates.
- `modelFromWorking()` maps working coordinates back to natural parameters.
- `modelWorkingLogJacobian()` evaluates the inverse-transformation log absolute Jacobian.

The base `psModel` transformation methods use an identity transform and zero Jacobian. This keeps unconstrained external models free from transformation boilerplate. The base prior/control methods fail explicitly because those quantities are model mathematics.

## Built-in compatibility methods

Zeta, logarithmic, and zero-inflated zeta expose their existing prior and transformation mathematics through the new public methods without changing their current fitting paths.

The ZIZ prior method combines the beta prior on `pi` with the existing shape `psPrior`. The existing ZIZ logit and `log(shape - 1)` transformations remain the underlying implementation.

## Stage boundary

Stage 9.2 deliberately does not change `fitPosterior()` or `fitModel.psModel()`. External Bayesian fitting therefore remains unavailable until a posterior engine is migrated to consume this contract. The next implementation stage can use the new methods to construct a model-neutral working-scale log posterior, with MCMC as the preferred first proof.

## Validation

The stage adds deterministic tests for external prior/control methods, identity defaults, built-in prior evaluation, transform round trips, Jacobians, and named built-in Bayesian starts. The stage runner performs the full package validation workflow because package code and tests change.
