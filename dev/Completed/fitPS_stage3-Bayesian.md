I want to continue work on the `jmcurran/fitPS` R package.

Current project state:
- Stage 2 shape reparameterisation is complete.
- `shape` now means the standard zeta parameter alpha, with `shape > 1`.
- VGAM's shifted parameterisation is now an internal implementation detail only.
- Frequentist compatibility against CRAN was sanity-checked on the Roux example:
  - `fitDist()` local shape equals CRAN shape + 1.
  - `fitDist()` `var.shape` is unchanged.
  - `fitZIDist()` `pi`, translated shape, `var.cov`, and fitted probabilities are unchanged.
- The package now has R Markdown vignettes after Stage 2.4.
- Use James Curran's R style:
  - `=` for assignment.
  - camelCase identifiers.
  - braces for all control structures.
  - roxygen2 documentation.
  - deterministic offline `testthat` tests.
  - keep changes targeted.
  - do not manually edit `NAMESPACE`; regenerate with `devtools::document()`.

New task:
I want to improve the Bayesian side of fitPS.

Main goal:
Add a plotting function that can display the posterior density for Bayesian fits, regardless of whether the Bayesian approximation came from:
- MCMC samples, or
- numerical integration / stored posterior density approximation.

The plot should show:
- posterior density of the relevant parameter, especially `shape`;
- the posterior point estimate, probably `fit$shape`;
- optionally a credible interval, e.g. 95%;
- support both `fitDist(..., method = "bayes")` and `fitZIDist(..., method = "bayes")` where practical.

Existing behaviour to inspect:
- `fitDist(..., method = "bayes")` appears to store `fit$chain` and `fit$pdf`.
- `fitZIDist(..., method = "bayes")` appears to store posterior samples in `fit$chain`, likely with `fit$chain$shape`.
- Existing `plot.psFit()` plots fitted probabilities, not posterior densities.
- Existing `credint()` uses KDE or posterior summaries internally and may contain reusable logic.
- Need to inspect current object structures before implementing.

Preferred implementation direction:
- Add an exported function such as `plotPosterior()` or an S3 method such as `plotPosterior.psFit()`.
- Avoid overloading existing `plot.psFit()` unless there is a clean `type = "posterior"` extension that preserves backward compatibility.
- The function should work when posterior samples are available and when only a posterior density function/grid is available.
- It should probably default to `parameter = "shape"`.
- For zero-inflated Bayesian fits, consider allowing `parameter = "shape"` and `parameter = "pi"`.
- For numerical-density objects, use the stored density function where available.
- For MCMC objects, use `density()` on the stored samples.
- Add deterministic tests that do not require long MCMC runs or external services. Use small/fake psFit objects if needed.
- Add roxygen documentation and examples that are fast and safe for `R CMD check`.
- Update NEWS with a version-numbered entry through the stage runner.

Possible API:
```r
plotPosterior(
  object,
  parameter = "shape",
  level = 0.95,
  showEstimate = TRUE,
  showInterval = TRUE,
  ...
)
