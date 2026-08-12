# fitPS Stage 6.4 numerical zeta engine migration

## Purpose

Stage 6.4 is the first behavioural migration onto the Stage 6 model/posterior-engine architecture. It routes the plain-zeta one-dimensional numerical Bayesian path through `zetaModel` and `numericalPosteriorEngine` while preserving the existing statistical calculation and `psFit` behaviour.

## Scope

The stage implements `modelLogLikelihood.zetaModel()` and the numerical engine operations `fitPosterior()`, `summarisePosterior()`, `posteriorDiagnostics()`, and `posteriorPointEstimate()`. A secondary `fitNumericalPosteriorModel()` dispatch keeps the numerical engine independent of concrete model names; Stage 6.4 provides only the `zetaModel` implementation.

The numerical representation retains the normalized posterior density together with posterior mean, variance, bounds, mode, and numerical integration error metadata. `fitDistBayesIntegrate()` becomes a thin orchestration layer that constructs the model and engine, obtains the representation and summary, and builds the existing `psFit` result.

## Behaviour preservation

The migrated implementation intentionally preserves the previous one-dimensional algorithm: optimize the negative log posterior over the prior range, shift by the posterior mode before exponentiation, normalize with `integrate()`, and obtain posterior first and second moments by numerical integration.

The fitted P/S probabilities continue to be evaluated at the posterior mean shape. Stage 6.4 does not redefine them as posterior-mean probabilities. The legacy `pdf` component remains available and now points to the density function stored in the numerical posterior representation.

## Architecture boundary

ZIZ numerical fitting remains explicitly unmigrated. Stage 6.4 does not add zeta Laplace or importance support and does not alter the MCMC path. This keeps the first engine proof one-dimensional and provides a regression-tested reference before wider migrations.

## Validation

Tests compare the migrated numerical zeta fit directly with an in-test implementation of the pre-Stage-6 algorithm for both P and S data. Additional tests cover the zeta model likelihood, numerical representation, summary, point estimate, diagnostics, normalization, public `fitDist()` routing, and explicit failure of the still-unmigrated numerical ZIZ protocol path.
