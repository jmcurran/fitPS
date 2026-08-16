# Stage 9.9 final Bayesian architecture audit

## Purpose

Stage 9 set out to make Bayesian fitting model-extensible in the same sense that
Stage 8 made maximum-likelihood fitting model-extensible. The required design
principle was:

> fitPS owns orchestration; the model owns the mathematics.

This audit records the final Stage 9 architecture and checks the original
completion criteria against the package state after Stage 9.8.

## Final public Bayesian model contract

External `psModel` classes can participate in fitPS-owned Bayesian engines by
providing model mathematics and controls rather than implementing a sampler or
integrator. The public contract now includes:

- `modelParameterNames()` for the ordered model parameter vector;
- `modelObservationData()` and `modelProbabilities()` for the established public
  model-data and probability interfaces;
- `modelLogLikelihood()` for likelihood evaluation;
- `modelLogPrior()` for complete model-specific prior evaluation;
- `modelBayesControl()` for starts and engine-specific bounds or controls;
- `modelToWorking()`, `modelFromWorking()`, and
  `modelWorkingLogJacobian()` for constrained MCMC parameterisations;
- `supportedPosteriorEngines()` for declaring supported fitPS-owned posterior
  engines.

External models do not need to implement numerical integration or an MCMC
sampler.

## Engine policy

The core Bayesian engine policy at the end of Stage 9 is:

- one-parameter models may use deterministic numerical integration through the
  generic numerical engine;
- two-parameter models may use deterministic adaptive integration through
  `cubature::hcubature()` and model-supplied finite natural-scale bounds;
- models with three or more parameters use MCMC rather than deterministic
  numerical integration;
- generic MCMC uses the public model likelihood, prior, controls, transformations,
  and Jacobian contract;
- Laplace and importance sampling remain available existing capabilities, but
  they are secondary and are not part of the core external-model extension
  promise.

When an external Bayesian fit does not explicitly request an engine, numerical
fitting is preferred for supported one- and two-parameter models and MCMC is
preferred for higher-dimensional models. An explicit unsupported numerical
request for a model above two dimensions fails clearly rather than silently
changing engines.

## Proof models

Stage 9 retains the Stage 8 external proof models entirely outside `R/`:

1. Poisson, as a one-parameter external model;
2. Poisson-normal, as a two-parameter external model whose marginal likelihood
   includes numerical integration over a latent normal log-rate.

The Poisson-normal probability calculation uses a standard-normal latent
variable transformation so that small values of `sigma` do not create a narrow
normal density that base `integrate()` can miss numerically.

## Completion criteria

### Built-in Bayesian fitting remains supported

Zeta, logarithmic, and zero-inflated zeta retain their established Bayesian
interfaces and regression coverage.

All three built-in models also satisfy the generic MCMC mathematics contract.
Specialised built-in samplers remain only where the regression suite treats
historical stochastic behaviour or exact RNG ordering as a compatibility
contract. They are compatibility overrides, not architectural requirements for
new models.

### External Poisson Bayesian fitting

Covered through both generic numerical fitting and generic MCMC without
`fitPS:::` in downstream model definitions.

### External Poisson-normal Bayesian fitting

Covered through generic two-dimensional numerical fitting and generic MCMC.
The fitted object retains the external model descriptor and supports the public
posterior workflow.

### Unsupported engines fail deliberately

Regression coverage checks deliberate failure for unsupported external posterior
engines and for numerical fitting above the two-parameter ceiling.

### Posterior summaries and comparison criteria

External Bayesian fits participate in posterior parameter summaries, derived
probability summaries, fitted probabilities, and DIC where the posterior
representation supports it. MLE fits retain AIC and BIC behaviour.

### P/S survey support mapping

External Poisson and Poisson-normal regression tests verify that matched P and S
surveys preserve the same underlying zero-based support sequence. No truncation
or renormalisation is introduced.

### Serialization

External MLE and Bayesian fitted objects are regression-tested through
`saveRDS()` / `readRDS()`, including retention of the originating external model
descriptor and posterior representation.

### Documentation

The package vignettes now teach the generic `fit()` model-descriptor interface,
the public Bayesian extension contract, the one-/two-dimensional numerical
policy, MCMC as the core simulation route, and Laplace/importance as secondary
capabilities. Deprecated distribution-specific fitting wrappers are described as
compatibility paths rather than the primary interface.

## Architectural decisions retained after Stage 9

The following are deliberate decisions, not incomplete Stage 9 work:

- specialised built-in MCMC methods may remain where exact historical stochastic
  behaviour is explicitly protected by tests;
- Laplace and importance sampling remain secondary existing engines rather than
  being expanded into the principal public external-model extension contract;
- deterministic numerical integration is intentionally capped at two model
  parameters;
- external proof distributions remain in tests and vignettes and are not added
  to package implementation code.

## Deferred work

The following work is outside the Stage 9 completion boundary:

- Bayesian Bootstrap as a later model-neutral uncertainty engine;
- adding the logarithmic distribution as an external extension tutorial only if
  a future documentation goal calls for another extension example;
- any future decision to relax exact-RNG compatibility and remove specialised
  built-in MCMC overrides;
- further model-file restructuring that would move still-needed compatibility
  implementations without changing behaviour.

## Conclusion

Stage 9 meets its original architectural goal. A downstream model author can now
provide the public `psModel` likelihood/prior/control/transform contract and use
fitPS-owned numerical or MCMC Bayesian machinery without rebuilding, forking, or
editing fitPS and without implementing a sampler. Numerical fitting is generic
for one- and two-parameter models, MCMC is generic for arbitrary supported model
dimensions, and the package retains compatibility-specific built-in overrides
only where they protect established behaviour.
