# Stage 5.7 audit: bootstrap distributions of fitted P and S probabilities

## Purpose

Stage 5.7 reviews the existing frequentist bootstrap implementation in fitPS and designs a staged path for obtaining uncertainty summaries for the fitted P and S probabilities that parallels the Bayesian `psPosterior` work.

The key frequentist quantity is not a Bayesian conditional expectation. For bootstrap replicate estimates `thetaHatStar_b`, the natural derived quantity is the bootstrap distribution

```text
P_k(thetaHatStar_b), b = 1, ..., B
```

with bootstrap mean

```text
E_*[P_k(thetaHatStar)] ~= mean_b P_k(thetaHatStar_b).
```

This generally differs from the plug-in probability `P_k(thetaHat)` because the probability mapping is nonlinear.

## Existing bootstrap architecture

The current bootstrap implementation is concentrated in `R/bootCI.R`.

### `bootFit()`

`bootFit()` performs a nonparametric case bootstrap by expanding the grouped survey counts to the observation vector, sampling `n * B` observations with replacement, converting each bootstrap sample back to `psData`, and refitting the selected model.

For the zeta model it returns a numeric vector of bootstrap shape estimates.

For the zero-inflated zeta model it returns a data frame with columns:

```text
pi
shape
```

These retained parameter replicates are exactly the information needed to construct bootstrap distributions of P and S probabilities.

### `bootCI()`

`bootCI()` calls `bootFit()` and then converts the bootstrap parameter replicates into parameter inference:

- for the zeta model, a smoothed bootstrap confidence interval for shape using `ks::hscv()` and `ks::kcde()`;
- for the zero-inflated zeta model, a smoothed joint confidence region for `(pi, shape)` using `ks::Hscv()`, `ks::kde()`, and contour levels.

With `returnBootValues = TRUE`, the underlying bootstrap parameter replicates are exposed as `bootVals`.

This means fitPS already generates the required bootstrap sampling distribution, but it does not currently represent that distribution as an object attached to a fit or transform the replicates into derived probabilities.

## Relationship to the Bayesian implementation

The Stage 5.2-5.6 Bayesian work established a useful pattern:

```text
psFit
  posterior -> psPosterior
      parameters
      probabilities
      representation
      level
      diagnostics
```

The bootstrap analogue should be structurally parallel without pretending that a bootstrap distribution is a posterior:

```text
psFit
  bootstrap -> psBootstrap
      parameters
      probabilities
      replicates
      level
      diagnostics
```

The interpretation remains different:

- `psPosterior$probabilities$estimate` is a posterior mean probability;
- `psBootstrap$probabilities$estimate` is a bootstrap mean probability.

The common table shape is nevertheless valuable because printing, plotting, prediction, and fitted-value methods can share internal presentation helpers.

## Recommended `psBootstrap` class

Introduce one internal constructor, provisionally:

```r
newPsBootstrap(
  method,
  parameters,
  probabilities,
  replicates,
  level = 0.95,
  diagnostics = NULL
)
```

The class should be:

```r
"psBootstrap"
```

Recommended components are:

### `method`

A short description of the bootstrap method. Initially this should identify the existing nonparametric observation bootstrap rather than overloading `model`.

### `parameters`

A data frame summarising the bootstrap distribution of fitted parameters. For ZIZ fits this can contain:

```text
parameter
estimate
sd
lower
upper
level
```

where `estimate` is the bootstrap mean of the parameter replicate, not the original MLE. The original point estimate remains available on the parent `psFit` object.

For the zeta model the same structure can contain only `shape`.

### `probabilities`

A data frame parallel to the posterior-probability summary:

```text
term
estimate
sd
lower
upper
level
bootstrapMethod
```

The `estimate` column is the bootstrap mean probability.

### `replicates`

Retain the raw fitted parameter replicates because they are the primary bootstrap representation and permit later recalculation of derived quantities without repeating model fits.

For ZIZ this should be a data frame with `pi` and `shape` columns. For zeta this should preferably also be normalised to a one-column data frame with `shape`, rather than retaining two different representation types.

### `diagnostics`

At minimum record:

```text
B
nSuccessful
nFailed
failureRate
```

This becomes important because some bootstrap resamples may not contain enough support to fit a ZIZ model, or optimisation may fail for individual resamples.

## Reuse the existing probability transformation

The zero-inflated bootstrap implementation should call the existing internal `zizProbabilities()` helper for every retained `(pi, shape)` replicate.

Do not implement a second P/S mapping inside the bootstrap code.

For the ordinary zeta model, either add a small matching internal helper or generalise the existing probability transformation so that both zeta and ZIZ fits use one probability-mapping boundary. A broad rewrite of `probfun()` is not required for the first bootstrap stage.

For ZIZ bootstrap replicate `b` and requested term `k`:

```text
T_bk = P_k(piHatStar_b, shapeHatStar_b)
```

or the corresponding S term.

Summaries should then be computed column-wise over the `B` transformed replicates.

## Bootstrap probability summaries

For each P or S term calculate at least:

- bootstrap mean;
- bootstrap standard deviation;
- percentile confidence interval.

The initial implementation should use percentile intervals because they follow directly from the empirical bootstrap distribution and do not require additional jackknife calculations.

Do not initially label the current `bootCI()` smoothed parameter contours as confidence intervals for the derived probabilities. They are parameter-space procedures and are not equivalent to transforming the bootstrap probability distribution.

### Possible later intervals

Basic bootstrap and BCa intervals could be added later, but BCa would require an acceleration estimate, normally from jackknife calculations. That is a separate computational feature and should not be added merely to mirror `bootCI()`.

## Should the bootstrap be attached automatically to every frequentist fit?

No.

Unlike the posterior object, a bootstrap requires an explicit and potentially expensive resampling computation. `fitDist()` and `fitZIDist()` should remain fast by default.

The recommended design is therefore:

```r
fit = fitZIDist(x)
fit = bootstrap(fit, B = 2000, seed = 1234)
```

or an equivalently clear fitPS-specific verb if `bootstrap()` would create an undesirable generic conflict.

After bootstrap computation:

```r
fit$bootstrap
```

would contain the `psBootstrap` object.

An alternative is for `bootCI(fit, ...)` to attach or return a bootstrap object, but mutating an existing R object passed by value is not natural, and replacing the established return value of `bootCI()` would be a compatibility risk. Therefore `bootCI()` should not become the primary constructor for `psBootstrap`.

## Recommended public API

The cleanest parallel to `posteriorProbs()` is:

```r
bootstrapProbs(fit)
bootstrapProbs(fit$bootstrap)
```

The generic should dispatch on both `psFit` and `psBootstrap`.

For the parent fit, it should fail clearly if no bootstrap has been computed.

Selection by `n` should match `posteriorProbs()` exactly so users do not have to learn two indexing conventions.

## `fitted()` semantics

Preserve the existing default:

```r
fitted(fit)
fitted(fit, type = "plugIn")
```

For a frequentist fit with an attached bootstrap object, add:

```r
fitted(fit, type = "bootstrapMean")
```

For Bayesian fits retain:

```r
fitted(fit, type = "posteriorMean")
```

Do not make bootstrap means the default frequentist fitted values. The ordinary MLE plug-in probabilities remain the conventional fitted values and changing the default would be a backward-compatibility break.

## `predict()` semantics

Once `psBootstrap` exists, `predict.psFit()` can eventually support the same distinction for requested P or S indices:

```r
predict(fit, newdata = ..., type = "plugIn")
predict(fit, newdata = ..., type = "bootstrapMean")
predict(fitBayes, newdata = ..., type = "posteriorMean")
```

For bootstrap predictions with intervals, percentile bounds can come directly from transformed bootstrap replicates.

This should be implemented only after the bootstrap object and `bootstrapProbs()` API are stable.

## `bootCI()` compatibility strategy

`bootCI()` is an established public function and should remain backward compatible.

Recommended approach:

1. factor the resampling/refitting code out of `bootFit()` into a reusable internal function that returns a stable bootstrap-replicate representation;
2. keep `bootFit()` as a compatibility wrapper if it is used internally or externally despite being unexported;
3. make `bootCI()` consume the new internal representation but preserve its existing return values;
4. make the new bootstrap-object constructor consume the same representation.

This prevents two independent bootstrap engines from developing.

## Bootstrap failures and invalid resamples

The existing `bootFit()` assumes every resample can be fitted. This is a significant issue for the new derived-probability use case, especially for ZIZ surveys.

A bootstrap resample can lose values needed for model identification. The current `fitZIDist()` checks can therefore reject some resamples, and optimisation may fail on others.

The new bootstrap engine should explicitly handle replicate failure:

- wrap each replicate fit safely;
- record failed replicates;
- summarise only successful fits;
- report `nSuccessful`, `nFailed`, and `failureRate`;
- stop if too few successful replicates remain for reliable summaries;
- never silently replace failed fits or reuse the original fitted parameters.

Tests should deliberately include a small survey capable of producing some invalid ZIZ resamples so this behaviour is characterised.

## Reproducibility and parallel execution

The existing bootstrap code calls `sample()` before parallel fitting. This is helpful because the resamples themselves are generated in the main R process, but a new public bootstrap API should still expose an explicit `seed` argument.

Recommended behaviour:

- `seed = NULL` preserves the caller's normal RNG behaviour;
- a supplied seed makes generated resamples reproducible;
- parallel and serial modes should use the same pre-generated resamples so results match for a given seed apart from any optimisation nondeterminism;
- tests should run serially by default for deterministic offline validation.

Avoid relying on the number of detected machine cores in unit tests.

## Parallel implementation concerns

The current code uses a mixture of `foreach`, `doParallel`, `pbapply`, and commented-out `parallel` alternatives. Stage 5.8 should not perform a broad parallelism refactor unless required for correctness.

However, cluster cleanup should use `on.exit()` or an equivalent safeguard so a fitting error cannot leak a parallel cluster.

The new internal replicate engine should separate:

1. creation of bootstrap samples;
2. conversion to `psData`;
3. model fitting;
4. summary/transformation.

This will make serial tests much easier and preserve optional parallel fitting for users.

## Plotting

A `plot.psBootstrap()` method should eventually parallel `plot.psPosterior()`:

- P or S term on the horizontal axis;
- bootstrap mean probability as a point;
- percentile confidence interval as a vertical interval.

Presentation code should ideally be shared between `psPosterior` and `psBootstrap` via an internal plotting helper that accepts a probability-summary data frame and labels. Do not create a new graphics dependency.

The existing `plot.psFit()` should remain the plot of ordinary fitted/observed probabilities unless an explicit uncertainty mode is added later.

## Printing and summaries

Recommended S3 methods are:

```r
print.psBootstrap()
summary.psBootstrap()
print.summary.psBootstrap()
fitted.psBootstrap()
plot.psBootstrap()
```

`fitted.psBootstrap()` should return the bootstrap mean probabilities, mirroring `fitted.psPosterior()` returning posterior mean probabilities.

The parent `print.psFit()` and `summary.psFit()` can report the presence of an attached bootstrap object without dumping all bootstrap replicates.

## Relationship between bootstrap mean and point estimate

The package documentation should distinguish:

```text
plug-in frequentist probability:
P_k(thetaHat)

bootstrap mean probability:
E_*[P_k(thetaHatStar)]

Bayesian posterior mean probability:
E[P_k(theta) | x]
```

These are three different quantities.

The bootstrap mean should not be described as `E[P_k | x]` because there is no posterior distribution in the frequentist analysis.

The difference between the plug-in and bootstrap mean is useful information about the effect of estimator variability and nonlinearity, but the bootstrap mean should not automatically replace the MLE plug-in estimate as the primary frequentist point estimate.

## Testing requirements

Stage 5.8 implementation tests should cover at least:

1. seeded serial reproducibility;
2. stable bootstrap replicate shape for zeta and ZIZ fits;
3. correct P and S term names;
4. bootstrap mean probabilities equal direct column means of transformed successful replicates;
5. bootstrap standard deviations equal direct transformed-replicate standard deviations;
6. percentile interval bounds in `[0, 1]`;
7. bootstrap mean probabilities differ from plug-in estimates for a deliberately nonlinear example;
8. probabilities from each successful replicate sum appropriately over a sufficiently long truncation;
9. failed ZIZ bootstrap fits are counted and not silently substituted;
10. existing `bootCI()` output remains compatible;
11. existing frequentist `fitted()`, `predict()`, `print()`, and `summary()` behaviour remains unchanged when no bootstrap is attached;
12. serial tests remain offline and modest in size.

## Compatibility risks

### `bootCI()` return values

High compatibility risk if changed. Preserve the existing interval/region return shapes and `returnBootValues` behaviour.

### Automatic bootstrap during fitting

High performance and behavioural risk. Do not bootstrap automatically in `fitDist()` or `fitZIDist()`.

### Meaning of `fitted()`

Keep `plugIn` as the default. Add `bootstrapMean` explicitly later.

### Failed bootstrap replicates

Current code effectively assumes all fits work. Handling failures will improve robustness but can change the number of retained bootstrap values. Report this explicitly in diagnostics.

### Parallel behaviour

Do not promise bit-for-bit equality between arbitrary historical parallel runs. Seeded serial reproducibility should be authoritative in tests.

### Existing raw `bootVals`

If the internal representation is normalised, ensure `bootCI(..., returnBootValues = TRUE)` still returns values in the legacy form expected by existing users.

## Proposed stages

### Stage 5.8 - Bootstrap object and replicate engine

- add `psBootstrap` internal class and constructors;
- factor bootstrap sample generation and model refitting into a reusable internal engine;
- add explicit seed support;
- retain successful parameter replicates and failure diagnostics;
- calculate parameter and P/S probability summaries;
- preserve `bootCI()` behaviour by routing it through the same engine;
- keep all new public probability accessors internal at this stage if API details still need settling.

### Stage 5.9 - Public bootstrap probability API and plotting

- export `bootstrapProbs()` for `psFit` and `psBootstrap`;
- add `fitted.psBootstrap()`;
- add `plot.psBootstrap()` using shared plotting infrastructure where practical;
- allow `fitted.psFit(type = "bootstrapMean")` only when a bootstrap object is attached;
- add user-facing help and examples for bootstrap mean probabilities.

### Stage 5.10 - Harmonise prediction and fit presentation

- extend `predict.psFit()` with explicit `type = "plugIn"`, `"bootstrapMean"`, and `"posteriorMean"` where available;
- support percentile intervals for bootstrap probability predictions;
- refine `print.psFit()` and `summary.psFit()` so posterior/bootstrap components are concise and consistent;
- repair remaining historical fitted/printing paths rather than duplicating probability calculations.

### Stage 5.11 - Documentation, vignette, and final Bayesian/bootstrap audit

- document the three probability estimands clearly;
- document all Bayesian posterior engines and bootstrap probability inference;
- add built-in-survey reproducible examples;
- explain credible versus confidence intervals;
- explain importance weights and Laplace approximation;
- update README and vignette material;
- perform a final backward-compatibility and regression audit.

## Recommendation

Proceed with the bootstrap extension before finalising the Bayesian documentation.

The existing `bootFit()` already supplies the essential frequentist sampling distribution of the fitted parameters, so the new functionality can be built without inventing a second inferential framework. The main architectural task is to make that distribution a first-class `psBootstrap` object, transform its retained replicates through the same P/S probability machinery, and preserve `bootCI()` as a backward-compatible parameter-inference interface.

This gives fitPS a coherent distinction between plug-in fitted probabilities, bootstrap mean probabilities, and Bayesian posterior mean probabilities while allowing the two uncertainty frameworks to share presentation and derived-probability infrastructure where their mathematics genuinely overlaps.
