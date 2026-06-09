# Stage 3.1 Bayesian posterior plotting audit

## Context

Stage 2 completed the public shape-parameterisation change: throughout the public API, `shape` now means the standard zeta parameter alpha and must satisfy `shape > 1`. VGAM's shifted parameterisation is an internal implementation detail only. Stage 3 should preserve that boundary while improving the Bayesian user experience.

The immediate Stage 3 goal is to add a posterior-density plotting function for Bayesian fit objects. The function should display posterior density for at least `shape`, mark the posterior point estimate, optionally show a credible interval, and work with both MCMC chains and stored numerical posterior-density approximations where practical.

## Current implementation findings

### Bayesian zeta fit by MCMC

`fitDist(..., method = "bayes")` delegates to `fitDistBayes()` and returns a `psFit` object with:

- `model = "zeta"`;
- `method = "bayes"`;
- `shape` as the posterior mean of the sampled chain;
- `var.shape` as `var(chain)`;
- `chain` as a numeric vector of shape samples;
- `pdf` as a spline function built from `density(chain, from = a)`.

This is the simplest target for `plotPosterior()`: the parameter is `shape`, the chain is one-dimensional, and both a sample and a density approximation are already stored.

### Bayesian zeta fit by numerical integration

`fitDist(..., method = "integrate")` delegates to `fitDistBayesIntegrate()` and returns a `psFit` object with:

- `model = "zeta"`;
- `method = "integrate"`;
- `shape` as the posterior mean computed by numerical integration;
- `var.shape` as the integrated posterior variance;
- `pdf` as a vectorised posterior-density function;
- no `chain` element.

Although this object is not currently marked as `method = "bayes"`, it is a Bayesian posterior approximation and should be supported by the new plotting function. The user-facing documentation should explicitly say that `plotPosterior()` supports Bayesian MCMC fits and integrated posterior-density fits.

### Bayesian zero-inflated zeta fit

`fitZIDist(..., method = "bayes")` delegates to `fitZIDistBayes()` and returns a `psFit` object with:

- `model = "ziz"`;
- `method = "bayes"`;
- `pi` as the posterior mean of the chain's `pi` column;
- `shape` as the posterior mean of the chain's `shape` column;
- `var.cov` as `cov(chain)`;
- `chain` as a data frame with columns `pi` and `shape`;
- no stored marginal `pdf` element.

This can support one-dimensional marginal posterior-density plots for `parameter = "shape"` and `parameter = "pi"` using `density(object$chain[[parameter]])`.

### Existing plot method

`plot.psFit()` currently plots fitted probabilities against observed frequencies. It is a useful fitted-distribution diagnostic, not a posterior diagnostic. Extending it with a posterior mode would increase argument complexity and risks breaking existing plotting expectations. A separate exported `plotPosterior()` function is cleaner.

### Existing credible interval function

`credint()` supports Bayesian fits only and currently checks `psFit$method != "bayes"`, which excludes integrated zeta posterior objects even though they contain a posterior-density function. For zeta MCMC fits, it computes intervals using `ks::kcde()`. For ZIZ fits, it computes two-dimensional credible regions using `ks::kde()` and contour extraction.

The new plotting function should not initially depend on the heavy two-dimensional `ks` machinery. For the Stage 3.2 implementation, a fast and deterministic one-dimensional interval based on sample quantiles is sufficient for MCMC-backed marginal plots. For integrated zeta fits, a numerical CDF built from a grid over the density support is the most direct route.

## Recommended API

Add an exported generic and an S3 method:

```r
plotPosterior = function(object, ...){
  UseMethod("plotPosterior")
}

plotPosterior.psFit = function(object,
                               parameter = "shape",
                               level = 0.95,
                               showEstimate = TRUE,
                               showInterval = TRUE,
                               ...){
  ...
}
```

Recommended behaviour:

- Accept only `psFit` objects.
- Accept `method = "bayes"` and `method = "integrate"`.
- Default to `parameter = "shape"`.
- For `model = "zeta"`, support only `parameter = "shape"`.
- For `model = "ziz"`, support `parameter = "shape"` and `parameter = "pi"`.
- Prefer MCMC samples when the requested parameter exists in `object$chain`.
- Fall back to `object$pdf` for zeta/integrate objects and zeta/bayes objects if needed.
- Draw a base R density plot to avoid adding a new plotting dependency.
- Mark `object[[parameter]]` with a vertical line when `showEstimate = TRUE`.
- Mark the credible interval with vertical lines or a lightly shaded interval when `showInterval = TRUE`.
- Return an invisible structured list with at least `parameter`, `estimate`, `interval`, and the plotted density grid. This makes deterministic tests possible without inspecting a device.

## Internal helper structure

Keep the public function thin and implement small internal helpers in a new file such as `R/plotPosterior.R`:

- `posteriorSamples(object, parameter)` extracts numeric samples from numeric-vector chains or data-frame chains.
- `posteriorDensityFromSamples(samples, ...)` wraps `stats::density()` and returns a standard list with `x` and `y`.
- `posteriorDensityFromPdf(object, parameter, n = 512)` evaluates `object$pdf` over a sensible grid.
- `posteriorIntervalFromSamples(samples, level)` returns equal-tailed quantiles.
- `posteriorIntervalFromGrid(x, y, level)` normalises the grid density and returns equal-tailed grid-quantile bounds.
- `posteriorParameterSupport(object, parameter)` gives support limits: `shape` uses lower bound just above 1 and a data-driven upper bound; `pi` uses `[0, 1]`.

The helpers should be unexported, deterministic, and directly unit tested where useful.

## Density-source decisions

### MCMC samples

Use `density(samples, ...)`. Remove non-finite samples before plotting and fail clearly if fewer than two finite samples remain. For ZIZ fits, use the requested chain column.

### Stored `pdf` function

For integrated zeta fits, evaluate `object$pdf` over a grid. The current integrated fit does not store the prior range used by the integration, so Stage 3.2 has two options:

1. infer an evaluation range from the estimate and variance, bounded below by `1 + .Machine$double.eps`; or
2. update `fitDistBayesIntegrate()` to store posterior support metadata, for example `posteriorSupport = c(lower = a, upper = b)` transformed onto the public shape scale.

The second option is more reliable and should be preferred if it can be added without changing user-facing behaviour. It also improves future numerical credible interval work.

## Shape-parameterisation safeguards

- All public plots and labels should use standard zeta `shape`, not VGAM shifted shape.
- Shape-density grids must not evaluate at or below 1.
- Tests should include at least one shape sample vector with values just above 1 to prevent accidental shifted-parameter assumptions.
- Labels should use `shape` and not `s - 1`, `alpha - 1`, or VGAM terminology.

## Tests to add in Stage 3.2

Add deterministic offline tests using small fake `psFit` objects rather than running long MCMC chains:

1. A zeta MCMC object with numeric `chain`, `shape`, `model = "zeta"`, and `method = "bayes"` returns an invisible result with a `shape` estimate and interval.
2. A zeta integrated object with `pdf`, `shape`, `model = "zeta"`, and `method = "integrate"` can plot and return a grid-based interval without `chain`.
3. A ZIZ MCMC object with data-frame `chain` supports `parameter = "shape"`.
4. The same ZIZ object supports `parameter = "pi"`.
5. Unsupported parameters fail clearly.
6. Non-Bayesian methods fail clearly.
7. Invalid `level` values fail clearly.
8. The function returns invisibly so scripts can call it without unwanted console output.

Use a temporary graphics device in tests, for example `pdf(tempfile(fileext = ".pdf"))`, and always close it with `on.exit(dev.off(), add = TRUE)`.

## Documentation to add in Stage 3.2

- New roxygen page for `plotPosterior()`.
- Fast examples based on fake objects or very small stored chains, not full MCMC runs.
- Cross-reference `fitDist()`, `fitZIDist()`, `plot.psFit()`, and `credint()`.
- Explain that `plot.psFit()` remains the fitted-probability plot while `plotPosterior()` is for posterior parameter uncertainty.

## NEWS guidance

The Stage 3 runner should update `NEWS.md` after successful validation and before commit using a version-numbered heading such as `## fitPS 1.0.5.001`. The NEWS entry should describe the audit in Stage 3.1 and, in later stages, the posterior plotting API and tests. It should not create a durable `## Stage 3.1` heading.

## Stage plan

### Stage 3.1: audit and workflow reset

- Add this audit document under `dev/`.
- Update the stage runner convention to use `1.0.5.xxx` versions, beginning at `1.0.5.001`.
- Ensure the stage runner removes old completed stage archives and old ChatGPT bundles while keeping only the newest relevant files.

### Stage 3.2: posterior plotting implementation

- Add `plotPosterior()` and `plotPosterior.psFit()`.
- Add internal density, interval, and parameter-extraction helpers.
- Add deterministic tests using fake `psFit` objects.
- Add roxygen documentation and regenerate `NAMESPACE`/`man` with `devtools::document()`.
- Add a version-numbered NEWS entry.

### Stage 3.3: posterior metadata hardening

- Store posterior support metadata from `fitDistBayes()` and `fitDistBayesIntegrate()` if Stage 3.2 reveals that grid inference is too fragile.
- Consider allowing `credint()` to support integrated zeta posterior objects.
- Add tests for integrated posterior support metadata.

### Stage 3.4: optional ZIZ posterior diagnostics

- Consider two-dimensional posterior contour plotting for ZIZ fits separately from the one-dimensional marginal `plotPosterior()` function.
- Reuse or improve `credint()` contour logic if this remains useful.
- Keep the one-dimensional `plotPosterior()` API stable.

## Open questions for implementation

- Should `fitDist(..., method = "integrate")` remain `method = "integrate"`, or should the object include an additional flag such as `posterior = TRUE`? The audit recommends leaving `method` unchanged and making `plotPosterior()` explicitly accept both `bayes` and `integrate`.
- Should credible intervals for MCMC plots use equal-tailed quantiles or highest posterior-density intervals? The audit recommends equal-tailed intervals for Stage 3.2 because they are deterministic, fast, dependency-free, and easy to explain. HPD intervals can be added later if needed.
- Should the plotting result include the raw samples? The audit recommends not returning raw samples by default, to avoid unnecessarily large invisible return values.
