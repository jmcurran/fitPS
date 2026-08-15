# fitPS Stage 4.1 audit and design: Bayesian zero-inflated posterior approximations

## Purpose

Stage 4 investigates deterministic and semi-deterministic posterior approximations for Bayesian zero-inflated zeta fits in `fitZIDist(..., method = "bayes")`. The immediate goal is not to change user-facing defaults. The goal is to define the posterior calculation boundary clearly enough that later stages can add tested helpers, a deterministic grid approximation, and optional experimental approximations without breaking the existing MCMC path.

## Current Bayesian zero-inflated zeta behavior

The current zero-inflated Bayesian implementation is concentrated in `R/fitZIDistBayes.R` and is called from `fitZIDist(..., method = "bayes")`.

Observed behavior from the current code base:

- The fitted object has class `psFit`, `model = "ziz"`, and `method = "bayes"`.
- The posterior output is stored as an MCMC chain in `fit$chain` with columns `pi` and `shape`.
- Point estimates are posterior chain means stored in `fit$par`, `fit$pi`, and `fit$shape`.
- `fit$var.cov` is the covariance matrix of the retained chain.
- The fitted probabilities are computed from the posterior mean using a local zero-inflated zeta probability helper.
- For P-type data, the observed counts are shifted by one before zeta evaluation; S-type data use the stored counts directly.
- The zeta shape parameter is on the standard fitPS scale and must satisfy `shape > 1`.
- The default prior for `pi` is `Beta(shape1, shape2)`, with default `shape1 = 1` and `shape2 = 1`.
- The shape prior is documented as uniform on `log(shape - 1)` over `[a, b]`, with default `a = -2` and `b = 2`.
- The sampler stores only retained draws after burn-in and does not currently store a reusable log-posterior object, marginal density functions, or deterministic posterior summaries.

The Stage 3 posterior plotting work already makes `plotPosterior()` useful for Bayesian zero-inflated fits because it can draw marginal densities from `fit$chain$pi` and `fit$chain$shape`. It can also plot stored one-dimensional density functions for non-zero-inflated integration fits. That means Stage 4 should preserve `fit$chain` for backward compatibility while adding new approximation structures that can be consumed by `plotPosterior()` without rerunning an approximation.

## Current reusable pieces

The following pieces are useful for Stage 4 but should be rationalized before implementation:

- `validateZetaShape()` already enforces the standard `shape > 1` convention.
- `dzetaStandard()` is the common zeta probability function on the current standard shape scale.
- `makePrior()` defines one-dimensional priors for ordinary zeta fits and documents the standard `shape > 1` scale.
- `fitDistBayesIntegrate()` shows a useful pattern for storing a normalized one-dimensional posterior density function as `fit$pdf`.
- `plotPosterior()` already prefers posterior samples where available and falls back to stored posterior density functions for supported cases.

However, zero-inflated Bayesian posterior evaluation is currently embedded inside `fitZIDistBayes()`. Stage 4.2 should extract this into small non-exported helpers so grid, Laplace, importance sampling, and MCMC can share the same likelihood and prior definitions.

## Posterior target and transformed scale

The recommended internal working scale is:

```text
eta = logit(pi)
tau = log(shape - 1)
```

with inverse transformations:

```text
pi = exp(eta) / (1 + exp(eta))
shape = 1 + exp(tau)
```

The standard-scale posterior is proportional to:

```text
likelihood(data | pi, shape) * priorPi(pi) * priorShape(shape)
```

When evaluating on the unconstrained transformed scale, the log posterior density for `eta, tau` must include Jacobian terms if the target is meant to be a density with respect to `eta, tau`:

```text
log posterior(eta, tau | data)
  = log likelihood(data | pi, shape)
  + log priorPi(pi)
  + log priorShape(shape)
  + log(pi) + log(1 - pi)
  + log(shape - 1)
```

For the default shape prior that is uniform on `tau = log(shape - 1)`, this simplifies because the standard-scale density contains `-log(shape - 1)` and the transformation adds `+log(shape - 1)`. In that case the prior contribution on the `tau` scale is constant within `[a, b]` and `-Inf` outside.

Stage 4.2 should make the density scale explicit in helper names or documentation. A suggested API is:

```r
zizTransformParams(etaTau)
zizInverseTransformParams(piShape)
zizLogLik(params, x)
zizLogPriorStandard(params, priorSpec)
zizLogPosteriorStandard(params, x, priorSpec)
zizLogPosteriorTransformed(etaTau, x, priorSpec)
```

The exact names can change, but the distinction between standard-scale and transformed-scale posterior density should not be ambiguous.

## MCMC code issues to verify before reuse

Before new approximations are compared against the existing sampler, Stage 4.2 should review the current Metropolis-Hastings algebra carefully. Points to verify include:

- The current proposal mechanism alternates between a beta proposal for `pi` and a uniform draw for `shape` over the transformed-prior support on the standard shape scale.
- The acceptance-ratio code combines likelihood, prior, and proposal terms inline, which makes it difficult to audit.
- The `pi` prior/proposal terms should be checked for symmetry and sign consistency between the current and proposed states.
- The normalizing constant used for the shape proposal or prior should be checked against the documented prior on `log(shape - 1)`.
- Any correction should be made in a targeted stage with regression tests, because the existing chain behavior is user visible.

This audit does not change the sampler. It recommends extracting the posterior target first, then using tests to compare old and new calculations on small deterministic examples.

## Candidate approximation methods

### Two-dimensional grid quadrature

Grid quadrature should be the first implemented approximation and the reference for later approximations.

Recommended approach:

- Build a rectangular grid on `eta, tau` rather than directly on `pi, shape`.
- Evaluate the transformed-scale log posterior at every grid point.
- Stabilize exponentiation by subtracting the maximum finite log posterior before exponentiating.
- Normalize with cell widths on the transformed scale.
- Store the normalized grid, marginal densities, marginal CDF helpers or quantile summaries, posterior means, variances, covariance, and credible intervals.
- Convert marginal density functions to standard parameter scales for `plotPosterior()` compatibility.

Expected storage shape:

```r
fit$posteriorApprox = list(
  method = "grid",
  scale = "etaTau",
  grid = ...,
  marginalPdf = list(pi = ..., shape = ...),
  posteriorSummary = ...,
  diagnostics = ...
)
```

The exact structure can be refined in Stage 4.3, but it should avoid overloading `fit$chain` with artificial samples unless compatibility requires a transitional path.

### Laplace approximation

Laplace approximation should initially be internal only.

Recommended approach:

- Optimize the transformed-scale log posterior over `eta, tau`.
- Use the negative inverse Hessian as the approximate covariance matrix.
- Guard against non-finite modes, non-positive-definite Hessians, and boundary modes.
- Store diagnostics, not just estimates.
- Use it to help choose grid bounds or as a proposal for importance sampling.

Laplace should not be exposed as a public approximation method until comparisons against grid results show acceptable behavior for representative small examples.

### Importance sampling

Importance sampling should also be internal or experimental at first.

Recommended approach:

- Use a bivariate normal proposal on `eta, tau`, initialized from the Laplace mode and an inflated covariance matrix.
- Compute log weights as transformed log posterior minus log proposal density.
- Normalize weights with log-sum-exp stabilization.
- Store weighted transformed samples, standard-scale samples, weights, effective sample size, and maximum normalized weight.
- Warn or mark diagnostics as failed when weights degenerate.

Importance sampling should not be treated as reliable without diagnostics.

## Recommended Stage 4 path

### Stage 4.2: posterior core helpers

Implement non-exported helpers for transformed parameters, zero-inflated zeta log likelihood, priors, and posterior evaluation. Add deterministic offline tests that do not rely on long MCMC runs. This stage should not change the default user API.

### Stage 4.3: deterministic grid posterior

Implement grid approximation as the first non-MCMC reference method. Store enough approximation data to compute marginal densities and credible intervals without rerunning the grid.

### Stage 4.4: Laplace helper

Add internal Laplace mode and covariance helper. Use it for diagnostics and possibly grid-bound initialization.

### Stage 4.5: importance sampling prototype

Add weighted posterior sampling with effective sample size and maximum-weight diagnostics.

### Stage 4.6: API and plotting integration

Decide whether to expose `posteriorMethod` through `fitZIDist(..., method = "bayes")` or to use a separate helper. Update `plotPosterior()` to consume `posteriorApprox$marginalPdf` for `pi` and `shape`.

## API recommendation

The safest user-facing design is to keep the existing default unchanged:

```r
fitZIDist(x, method = "bayes")
```

should continue to mean the current MCMC behavior unless the user explicitly asks for another posterior method. A future explicit API could be:

```r
fitZIDist(x, method = "bayes", posteriorMethod = "mcmc")
fitZIDist(x, method = "bayes", posteriorMethod = "grid")
```

Laplace and importance sampling should remain internal until their diagnostics and accuracy have been tested.

## Testing principles

Future implementation stages should use deterministic offline tests that cover:

- parameter transformation round trips;
- rejection of invalid `pi` and `shape` values;
- zero-inflated likelihood calculations for tiny hand-checkable examples;
- finite posterior values inside the prior support and `-Inf` values outside it;
- grid weights summing to one after normalization;
- marginal posterior summaries lying inside valid parameter ranges;
- unchanged default MCMC API behavior;
- `plotPosterior()` compatibility for both `pi` and `shape` where posterior approximations are stored.

Tests should be small and fast. Long MCMC comparisons should not be required for normal package checks.

## Stage 4.1 conclusion

Stage 4 should first make posterior evaluation explicit and shared. The recommended implementation path is:

1. Extract and test zero-inflated posterior core helpers.
2. Implement transformed-scale grid quadrature as the deterministic reference approximation.
3. Add Laplace as an internal diagnostic and proposal helper.
4. Add importance sampling only with diagnostics.
5. Integrate new marginal posterior storage with `plotPosterior()` while preserving existing `fit$chain` behavior.
