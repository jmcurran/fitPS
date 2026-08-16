# fitPS

`fitPS` fits probability models to forensic P- and S-survey data, especially clothing-survey data describing the background presence of glass and similar trace material. The package provides a common model-fitting interface for maximum likelihood, parametric Bayesian inference, the ordinary nonparametric bootstrap, and Rubin's Bayesian Bootstrap.

The two survey types are:

- `P` data: counts of the number of groups or sources found on an item;
- `S` data: counts of the sizes of those groups.

These probability terms arise in activity-level interpretation of trace evidence. The package is intended to make the fitted models, uncertainty calculations, and their assumptions explicit rather than treating estimated P/S probabilities as known constants.

## Installation

Install from the repository root during development:

```r
devtools::install()
```

or from GitHub:

```r
remotes::install_github("jmcurran/fitPS")
```

## A common fitting interface

The canonical public entry point is `fit()`:

```r
library(fitPS)
data("Psurveys")
roux = Psurveys$roux

mleFit = fit(
  roux,
  model = zetaModel(),
  method = "mle"
)
```

The model and inferential method are chosen separately. Built-in model descriptors include:

- `zetaModel()`;
- `zizModel()` for the zero-inflated zeta model;
- `logarithmicModel()`.

The main inferential methods are:

```r
fit(roux, model = zetaModel(), method = "mle")
fit(roux, model = zetaModel(), method = "bayes")
fit(roux, model = zetaModel(), method = "bootstrap")
fit(roux, model = zetaModel(), method = "bayesianBootstrap")
```

They answer related but distinct questions. Maximum likelihood gives a fitted parametric model. Parametric Bayesian inference propagates uncertainty in model parameters under an explicit prior. The ordinary bootstrap resamples surveyed observational units and refits the model. Rubin's Bayesian Bootstrap keeps the observed empirical support fixed and places Dirichlet uncertainty on its weights before weighted refitting.

## Data format

Input files for `readData()` contain two columns:

- one column named `P` or `S`;
- one column named `count`.

For `P` data, the `P` column contains values such as `0`, `1`, `2`, and so on. For `S` data, the `S` column contains group sizes such as `1`, `2`, `3`, and so on.

Example CSV:

```csv
P,count
0,98
1,1
2,1
```

Data can also be constructed directly:

```r
pData = makePSData(
  n = c(0, 1, 2),
  count = c(98, 1, 1),
  type = "P"
)

sData = makePSData(
  n = 1:3,
  count = c(1, 1, 1),
  type = "S"
)
```

or read from a file:

```r
pData = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
sData = readData(system.file("extdata", "s.xlsx", package = "fitPS"))
```

## Fitted probabilities

A fitted `psFit` can be printed, plotted, predicted from, or converted to a probability function:

```r
mleFit = fit(roux, model = zizModel())

mleFit
fitted(mleFit)
predict(mleFit)
pFun = probfun(mleFit)
pFun(0:5)
```

For zeta-based models, `shape` is the zeta shape parameter on the standard `shape > 1` scale.

## Parametric Bayesian inference

Bayesian fitting is selected with `method = "bayes"`. The posterior engine is selected separately through `bayesOptions`.

```r
bayesFit = fit(
  roux,
  model = zizModel(),
  nterms = 6,
  method = "bayes",
  bayesOptions = list(posteriorMethod = "numerical")
)

posteriorProbs(bayesFit, n = 6)
fitted(bayesFit, n = 6, type = "posteriorMean")
plot(bayesFit$posterior, n = 6)
```

Posterior probability summaries average the model-implied probabilities over the posterior distribution of the parameters. They are not, in general, equal to probabilities evaluated at posterior mean parameters.

## Ordinary nonparametric bootstrap

The frequentist bootstrap is also requested through `fit()`:

```r
bootFit = fit(
  roux,
  model = zizModel(),
  method = "bootstrap",
  nterms = 6,
  B = 2000,
  seed = 1234,
  silent = TRUE
)

bootstrapProbs(bootFit, n = 6)
fitted(bootFit, n = 6, type = "bootstrapMean")
plot(bootFit$bootstrap, n = 6)
```

The ordinary bootstrap resamples observational units with replacement. A replicate can therefore lose occupied support values and, in sparse data, may occasionally produce a sample for which the requested MLE does not exist. Such failures are part of the uncertainty behaviour of the estimator rather than being repaired by smoothing or redraws.

`bootstrapFit()` is retained only as a deprecated compatibility wrapper. New code should use `fit(..., method = "bootstrap")`.

## Rubin's Bayesian Bootstrap

Rubin's Bayesian Bootstrap is requested with `method = "bayesianBootstrap"`:

```r
bayesBoot = fit(
  roux,
  model = zizModel(),
  method = "bayesianBootstrap",
  nterms = 6,
  B = 2000,
  seed = 1234
)

summary(bayesBoot, nterms = 6)
```

At the individual-observation level, Rubin's Bayesian Bootstrap assigns `Dirichlet(1, ..., 1)` weights. For aggregated fitPS survey data, equivalent category weights are drawn from `Dirichlet(n1, ..., nk)`, where the `nj` are the observed category counts. All originally occupied categories therefore retain positive weight with probability one. This differs from both the ordinary bootstrap and a parametric Bayesian posterior.

## Plotting parameter uncertainty

`plotUncertainty()` gives the four inferential methods a common plotting interface while retaining their different statistical interpretations:

```r
plotUncertainty(mleFit, level = c(0.80, 0.95))
plotUncertainty(bootFit, level = c(0.80, 0.95))
plotUncertainty(bayesFit, level = c(0.80, 0.95))
plotUncertainty(bayesBoot, level = c(0.80, 0.95))
```

For a two-parameter model such as ZIZ, these are respectively profile-likelihood confidence regions, smoothed bootstrap confidence regions, posterior credible regions, and Rubin Bayesian-Bootstrap weighted-fit regions. Bootstrap-based displays show their stored parameter realizations beneath the KDE contours by default. Sample-based methods reuse stored parameter draws, importance sampling retains its weights, and Laplace fitting uses its stored Gaussian covariance approximation. Numerical Bayesian fits reuse the posterior representation retained during fitting; in one dimension the stored cumulative representation uses Simpson quadrature where the grid permits it, and in two dimensions stored posterior density and quadrature mass are used to determine credible-region contour thresholds without rerunning adaptive cubature.

When numerical results rather than a plot are required, `confint()` remains the frequentist interval/region extractor and `credint()` is the corresponding Bayesian extractor. `credint()` reuses the stored parametric-posterior or Rubin Bayesian-Bootstrap representation; `plotUncertainty()` is the common visual layer.

## Model comparison

Fitted models can be compared with AIC and BIC, and Bayesian fits can also provide DIC where the posterior representation supports it. These criteria address model comparison; they are conceptually separate from bootstrap or posterior probability uncertainty.

See the model-comparison vignette for worked examples.

## Extending fitPS

New probability models can be supplied through the public `psModel` contract rather than by adding another specialised fitting function. The extension vignette demonstrates both a simple external Poisson model and a two-parameter Poisson-normal model using the common likelihood, probability, parameter-transformation, and Bayesian-engine interfaces.

## Vignettes

The package contains four main vignettes:

- **fitPS basics**: survey data, zeta/ZIZ fitting, P/S probabilities, and survey comparison;
- **Probability uncertainty with fitPS**: plug-in estimates, parametric Bayes, ordinary bootstrap, and Rubin's Bayesian Bootstrap;
- **Model comparison with fitPS**: AIC, BIC, DIC, and fitted-probability comparisons;
- **Extending fitPS with new distributions**: the public model-extension contract.

Use:

```r
browseVignettes("fitPS")
```

## Backward compatibility

Older functions such as `fitDist()`, `fitZIDist()`, and `bootstrapFit()` remain only where needed for compatibility and issue deprecation guidance. New analyses should use `fit()` with an explicit model descriptor and inferential method.

## References

Coulson, S. A., Buckleton, J. S., Gummer, A. B., and Triggs, C. M. (2001). Glass on clothing and shoes of members of the general population and people suspected of breaking crimes. *Science & Justice*, 41(1), 39-48. https://doi.org/10.1016/S1355-0306(01)71847-3

Curran, J. M., Buzzini, P., and Trejos, T. (2024). Estimating probability terms for the background presence of glass when considering activity in forensic casework. *Forensic Science International*, 364, 112221. https://doi.org/10.1016/j.forsciint.2024.112221

Evett, I. W., and Buckleton, J. S. (1990). The interpretation of glass evidence. A practical approach. *Journal of the Forensic Science Society*, 30(4), 215-223.

Efron, B. (1979). Bootstrap methods: Another look at the jackknife. *The Annals of Statistics*, 7(1), 1-26.

Rubin, D. B. (1981). The Bayesian bootstrap. *The Annals of Statistics*, 9(1), 130-134.

## License

GPL (>= 2).
