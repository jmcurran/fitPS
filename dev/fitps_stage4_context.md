# fitPS Stage 4 context: Bayesian zero-inflated posterior approximations

I want to continue work on the `jmcurran/fitPS` R package in a new stage series.

## Current state

- Work is on the `bayes` branch.
- Stage 3 is complete through Stage 3.4.
- Stage 3 added and stabilized posterior plotting support via `plotPosterior()`.
- `plotPosterior()` supports posterior densities for Bayesian `psFit` objects using:
  - stored MCMC samples in `fit$chain`, including zero-inflated zeta fits with `chain$pi` and `chain$shape`;
  - stored one-dimensional posterior density functions, where available, using `approxfun()` / `integrate()` / `uniroot()` rather than hand-written trapezoid helpers.
- Stage 3 also fixed stage-runner workflow issues, including the controlled git staging bug and the no-longer-used ChatGPT bundle step.

## Version convention for Stage 4

Use fitPS versions of the form:

```text
1.0.6.xxx
```

where `xxx` starts at `001` and is incremented for every build attempt, regardless of success or failure.

Stage 4 should start at `4.1`.

## Development conventions

Use James Curran's R development style:

- `=` for assignment;
- camelCase identifiers;
- braces for all control structures;
- roxygen2 documentation;
- deterministic offline `testthat` tests;
- targeted changes only;
- do not manually edit `NAMESPACE`; regenerate with `devtools::document()`;
- prefer `@importFrom` roxygen imports over `pkg::fun()` in package code where sensible.

Use the corrected fitPS stage-script-generator workflow:

- deliver a separate `stage4_X_changes.zip` and `run_stage4_X.sh`;
- support `--install-files` and `--changes-zip`;
- if a stage has no new or modified R files, tests, generated docs, or package-impacting files, use the lightweight version/news/commit workflow only;
- for package-impacting stages, run document, tests, strict check, update NEWS, commit, archive, build, install;
- no ChatGPT bundle step.

## Stage 4 main goal

Investigate and implement non-MCMC Bayesian posterior approximations for `fitZIDist(..., method = "bayes")` or a related explicit method, with focus on the two-parameter zero-inflated zeta posterior over:

```text
pi
shape
```

The current `fitZIDist(..., method = "bayes")` appears to use a simple Metropolis-Hastings sampler and stores posterior samples in `fit$chain`, including columns `pi` and `shape`.

Stage 4 should explore whether fitPS can also support deterministic or semi-deterministic posterior approximation for the zero-inflated zeta model.

## Candidate posterior methods to evaluate

### 1. Two-dimensional grid quadrature

This is the baseline deterministic approach.

For each grid point `(pi, shape)`, compute the unnormalized posterior:

```text
posterior(pi, shape | data) proportional to
  likelihood(data | pi, shape) * priorPi(pi) * priorShape(shape)
```

Then normalize over the two-dimensional grid.

Expected advantages:

- deterministic and reproducible;
- straightforward for two parameters;
- allows marginal posterior densities for `pi` and `shape` by summing/integrating over the other dimension;
- integrates naturally with `plotPosterior()` if the fitted object stores marginal posterior density functions or grids.

Expected disadvantages:

- can be expensive if the grid is too fine;
- accuracy depends on sensible parameter bounds and grid spacing;
- inefficient for narrow or highly correlated posteriors unless adaptive grid refinement is added.

Implementation notes:

- Work on stable transformed scales where helpful:
  - logit scale for `pi`;
  - `log(shape - 1)` for `shape`, preserving `shape > 1`.
- Be explicit about Jacobian terms if evaluating/posterior-normalizing on transformed scales.
- Store enough posterior approximation information to allow marginal plotting and credible intervals without rerunning the approximation.
- Prefer marginal densities/functions for compatibility with `plotPosterior()`.

### 2. Laplace approximation

Evaluate a Laplace approximation around the posterior mode.

This would approximate the posterior on a transformed unconstrained scale, probably:

```text
eta = logit(pi)
tau = log(shape - 1)
```

Then:

1. optimize the log posterior over `(eta, tau)`;
2. compute or approximate the Hessian at the posterior mode;
3. use the inverse negative Hessian as an approximate covariance matrix;
4. transform samples or marginal summaries back to `pi` and `shape`.

Expected advantages:

- fast;
- useful as a starting approximation or diagnostic;
- can provide proposal distributions for importance sampling.

Expected disadvantages:

- can be inaccurate when the posterior is skewed, bounded, multimodal, or strongly non-Gaussian;
- credible intervals after back-transformation need care;
- should probably not be the only approximation without diagnostics or comparison to MCMC/grid results.

Potential Stage 4 use:

- add a non-exported helper to compute a Laplace approximation;
- use it to initialize grid bounds or as an importance-sampling proposal;
- optionally expose it later as `method = "laplace"` only if the approximation is good enough.

### 3. Importance sampling

Use an importance sampler to approximate the posterior with weighted draws.

A natural proposal is a bivariate normal approximation on the transformed scale from the Laplace approximation:

```text
proposal(eta, tau) = N(mode, inflated covariance)
```

Then:

1. draw proposal samples on `(eta, tau)`;
2. transform to `(pi, shape)`;
3. compute log posterior and log proposal density;
4. compute normalized importance weights;
5. estimate marginal posterior summaries and densities using weighted samples.

Expected advantages:

- more flexible than pure Laplace;
- can be deterministic-ish if seeded;
- can reuse the Laplace approximation as a proposal;
- can approximate posterior quantities without a rectangular grid.

Expected disadvantages:

- can be unstable if the proposal misses posterior tails;
- requires diagnostics such as effective sample size and maximum weight share;
- weighted density estimation and credible intervals need careful implementation.

Potential Stage 4 use:

- implement as an internal experimental approximation first;
- compare against MCMC on small deterministic examples;
- store weighted samples in the fit object, perhaps as `fit$weightedSamples` or similar, without breaking existing `fit$chain` behavior.

## Suggested Stage 4 plan

### Stage 4.1: audit and design

Create a dev audit/design note that:

- inventories current `fitZIDistBayes()` MCMC behavior;
- identifies reusable likelihood/prior code;
- specifies a stable transformed-parameter posterior API;
- compares grid quadrature, Laplace approximation, and importance sampling;
- recommends the first implementation path.

This stage may be documentation-only, so it can use the lightweight runner if it only adds a `dev/` note.

### Stage 4.2: posterior core helpers

Add internal helpers for:

- transformed parameters `(eta, tau)` to `(pi, shape)` and back;
- log likelihood for zero-inflated zeta fits;
- log prior density on the chosen working scale;
- log posterior evaluation;
- deterministic validation of parameter bounds and inputs.

Add offline tests for helper behavior.

### Stage 4.3: two-dimensional grid posterior

Implement a deterministic grid posterior approximation for zero-inflated zeta fits.

Possible object storage:

```r
fit$posteriorGrid
fit$marginalPdf
fit$posteriorSummary
```

or a similarly named structure that does not break existing `psFit` behavior.

Add tests using tiny fake or small real `psData` objects that run quickly.

### Stage 4.4: Laplace approximation

Implement an internal Laplace approximation helper and compare it to grid/MCMC on small examples.

Do not expose as public API until diagnostics and accuracy are acceptable.

### Stage 4.5: importance sampling prototype

Implement weighted posterior sampling using the Laplace approximation as the proposal.

Include diagnostics such as:

- effective sample size;
- maximum normalized weight;
- warnings for weight degeneracy.

### Stage 4.6: API and plotting integration

Decide whether users should request approximations with:

```r
fitZIDist(..., method = "bayes", posteriorMethod = "mcmc")
fitZIDist(..., method = "bayes", posteriorMethod = "grid")
fitZIDist(..., method = "bayes", posteriorMethod = "laplace")
fitZIDist(..., method = "bayes", posteriorMethod = "importance")
```

or a simpler separate helper/API.

Ensure `plotPosterior()` works for the new approximation objects, especially marginal posterior densities for `pi` and `shape`.

## Accuracy and safety principles

- Prefer deterministic grid quadrature as the reference method for two-parameter examples.
- Treat Laplace as a fast approximation and proposal generator, not automatically as the gold standard.
- Treat importance sampling as useful only when diagnostics indicate stable weights.
- Keep existing MCMC behavior backward compatible.
- Avoid changing `fitZIDist()` defaults until the new methods are tested.
- Add NEWS entries through the stage runner, version-numbered under `1.0.6.xxx`.

## Merge and next branch workflow

Before starting Stage 4, the current `bayes` branch should be merged back into `master`. After the merge, create and push a new branch called:

```text
bayes-ii
```

A merge helper script has been generated for this handoff.
