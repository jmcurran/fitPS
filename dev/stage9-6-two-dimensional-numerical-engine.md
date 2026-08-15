# Stage 9.6: two-dimensional numerical posterior engine

## Purpose

Stage 9.6 extends the generic Bayesian numerical engine from one to two model
parameters while putting a deliberate ceiling on deterministic integration.

The fitPS policy is now:

- one parameter: deterministic numerical integration with base `integrate()`;
- two parameters: adaptive deterministic integration with `cubature::hcubature()`;
- three or more parameters: numerical fitting is refused and MCMC is the
  supported scalable route.

When an external model supports both numerical and MCMC fitting and the user
does not choose an engine explicitly, fitPS selects numerical fitting for one-
and two-parameter models and MCMC for models with three or more parameters.
An explicit request for the numerical engine above two dimensions fails rather
than silently changing the requested method.

## Architecture

Two-dimensional integration uses the public Stage 9 Bayesian model contract.
The model owns likelihood, prior, and finite natural-scale numerical bounds.
`modelBayesControl()` supplies named `start`, `lower`, and `upper` vectors, and
the numerical engine integrates directly over that model-supplied rectangle.
Working-scale transformations and Jacobians remain available to MCMC and other
algorithms that benefit from unconstrained coordinates, but deterministic
cubature does not impose an additional whole-real-line transformation.

The cubature call integrates normalization, first moments, second moments, and
expected deviance together as a vector-valued integrand. A small deterministic
natural-scale grid is retained for derived posterior probability summaries.

ZIZ now obtains its posterior normalization, moments, and DIC from the generic
cubature engine. Its established rectangular `posteriorGrid` field is retained
as a compatibility view for plotting and inflation-probability helpers, but it
is no longer the source of the fitted posterior moments.

## Proofs

- Built-in ZIZ exercises generic two-dimensional cubature while retaining the
  existing user-facing numerical Bayesian workflow.
- The external Poisson-normal model exercises the same two-dimensional engine
  without package-specific sampler or integration code.
- The external Poisson and built-in zeta/logarithmic models continue to use the
  generic one-dimensional `integrate()` path.
- A three-parameter model is rejected by the numerical engine with an explicit
  recommendation to use MCMC.

## Dependency

`cubature` is now an imported package because two-dimensional deterministic
posterior integration is core fitPS functionality rather than an optional
example.
