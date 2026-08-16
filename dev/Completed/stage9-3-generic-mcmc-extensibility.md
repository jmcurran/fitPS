# Stage 9.3 generic MCMC extensibility

## Purpose

Stage 9.3 makes MCMC the primary generic Bayesian extension path for external
`psModel` classes. The model owns likelihood, prior, starting values, and
natural/working-scale transformations; fitPS owns random-walk Metropolis
sampling and posterior finalisation.

## Scope

- Add a model-neutral random-walk Metropolis implementation to the `psModel`
  fallback of the MCMC posterior engine.
- Allow external `fit(..., method = "bayes")` calls to use model-specific prior
  objects rather than requiring the legacy scalar `psPrior` structure.
- Default generic external Bayesian fitting to MCMC and reject other generic
  posterior engines deliberately for now.
- Prove the path with the external Poisson and Poisson-normal models used by the
  Stage 8 extensibility tests.
- Preserve built-in zeta, ZIZ, and logarithmic specialised MCMC methods so their
  established behaviour and random-number sequences are not changed by this
  stage.

## Secondary posterior engines

Laplace and importance sampling remain available where already implemented,
but Stage 9 does not make them a central part of the public extension contract.
They are retained as secondary/legacy capabilities rather than expanded or
promoted. Numerical posterior support likewise remains separate until there is
a clear model-generic use case.

## Generic MCMC contract

An external model advertising `supportedEngines = "mcmc"` supplies:

- `modelLogLikelihood()`;
- `modelLogPrior()`;
- `modelBayesControl()` with a natural-scale `start` vector;
- `modelToWorking()`;
- `modelFromWorking()`;
- `modelWorkingLogJacobian()`; and
- `modelProbabilities()` for posterior fitted/probability summaries.

The generic sampler proposes a symmetric Gaussian random walk on the working
scale, evaluates likelihood + prior + inverse-transform log-Jacobian, stores
retained draws on the natural scale, and lets the existing common posterior
orchestration build the `psFit` and `psPosterior` objects.

## Prior policy

The generic external path requires an explicit prior. fitPS does not invent a
default prior and does not validate the prior against the legacy `psPrior`
shape. The external model's `modelLogPrior()` method owns that validation.
This permits multi-parameter and correlated prior representations in later
extensions without changing the engine.

## Validation expectations

The Stage 9.3 full validation path should prove that existing package tests
remain green and that deterministic external Poisson and Poisson-normal MCMC
fits work without `fitPS:::` access. The Poisson-normal test is intentionally
small because its likelihood performs numerical integration at every posterior
evaluation.
