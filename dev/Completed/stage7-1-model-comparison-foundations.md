# fitPS Stage 7.1 - model-comparison foundations

## Purpose

Stage 7.1 establishes model-comparison infrastructure before the logarithmic distribution is migrated into the `psModel` architecture.

## Contract

- The model owns likelihood and deviance mathematics through `modelLogLikelihood()` and `modelDeviance()`.
- `logLik.psFit()` is restricted to MLE fits and supplies the metadata required by `stats::AIC()` and `stats::BIC()`.
- `DIC()` is restricted to Bayesian fits and obtains posterior expected deviance from the common posterior representation.
- Numerical, MCMC, and importance posterior representations provide DIC directly from their integration weights or draws.
- Laplace DIC currently fails explicitly because the existing representation stores a local approximation but not a posterior deviance expectation. A future extension should improve that representation rather than silently substitute an unrelated criterion.

## API direction

`fitDist()` and `fitZIDist()` remain established entry points. A later stage may add a model-oriented `fit()`-style interface in which the user specifies the model explicitly. Stage 7 should avoid decisions that would prevent that unification.

## Logarithmic distribution

The existing logarithmic-distribution fitting code is provisional and has not been exposed in the literature. It is not a backwards-compatibility target. Stage 7.2 may replace or remove it so that the logarithmic distribution becomes a clean `psModel` implementation and a genuine test of extensibility.

## Deferred Bayesian bootstrap context

The former `dev/stage7-1-bayesian-bootstrap-context.md` has been retained as `dev/future-bayesian-bootstrap-context.md`. It remains future work and no longer claims Stage 7.1.
