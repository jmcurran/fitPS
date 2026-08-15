# Stage 6.6 - numerical ZIZ engine migration

## Purpose

Stage 6.6 migrates the zero-inflated zeta numerical Bayesian fitting path onto
the shared Stage 6 model/posterior-engine architecture. It is the first ZIZ
engine migration after both plain-zeta engines established the protocol.

## Scope

The stage:

- implements `modelLogLikelihood.zizModel()` through the characterised ZIZ
  likelihood;
- implements numerical posterior fitting for `zizModel` through
  `numericalPosteriorEngine`;
- preserves the existing two-dimensional rectangular posterior-grid algorithm,
  including beta prior, shape prior, mode rescaling, grid integration, moments,
  and marginal densities;
- wraps the existing posterior grid in the common numerical posterior
  representation without forcing it into the one-dimensional zeta storage
  layout;
- routes `fitZIDistBayesNumerical()` through `zizModel()` and
  `numericalPosteriorEngine()`;
- uses `modelProbabilities()` for fitted P/S probabilities while preserving the
  previous posterior-mean parameter calculation;
- retains the existing raw posterior grid inside `psPosterior` so current
  posterior probability and inflation methods continue to behave unchanged;
- adds deterministic regression tests for both P and S data;
- leaves ZIZ MCMC, Laplace, and importance fitting unchanged.

## Architecture

The numerical engine now supports both current models through model dispatch:

```text
numericalPosteriorEngine
  +-- zetaModel -> one-dimensional adaptive integration
  +-- zizModel  -> two-dimensional rectangular posterior grid
```

The engine-level summary and point-estimate methods are shared because both
representations expose named posterior means and variance-covariance matrices.
The physical numerical representation remains model-appropriate.

## Compatibility decision

Stage 6.6 preserves the numerical ZIZ posterior grid, posterior probability
summaries, marginal density functions, fitted probabilities, parameter moments,
and current `psPosterior` representation. The new engine representation is
added alongside those existing components so architectural migration can be
validated before Stage 6.9 consolidates the final posterior interface.

## Deferred work

Stage 6.6 does not migrate ZIZ MCMC, Laplace, or importance sampling and does
not yet refactor inflation-probability dispatch. Those remain Stages 6.7-6.9.
