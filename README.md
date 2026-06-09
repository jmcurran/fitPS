# fitPS

`fitPS` fits zeta-distribution models to forensic survey data, especially count data from clothing surveys involving glass or paint transfer.

The package supports two related survey data types:

- `P` data: counts of the number of groups or sources found on clothing.
- `S` data: counts of group sizes.

The main workflow is:

1. Create or read a `psData` object.
2. Fit a model with `fitDist()` or `fitZIDist()`.
3. Inspect, predict, or compare fitted `psFit` objects.

## Installation

Install from the repository root during development:

```r
devtools::install()
```

or from GitHub when the repository is available:

```r
remotes::install_github("jmcurran/fitPS")
```

## Data format

Input files for `readData()` must contain exactly two columns:

- one column named `P` or `S`;
- one column named `count`.

For `P` data, the `P` column contains counts such as `0`, `1`, `2`, and so on. For `S` data, the `S` column contains group sizes such as `1`, `2`, `3`, and so on.

Example CSV:

```csv
P,count
0,98
1,1
2,1
```

## Creating data manually

```r
library(fitPS)

pData = makePSData(
  n = c(0, 1, 2),
  count = c(98, 1, 1),
  type = "P"
)

pData
```

For `S` data:

```r
sData = makePSData(
  n = 1:3,
  count = c(1, 1, 1),
  type = "S"
)

sData
```

## Reading data from a file

```r
pData = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
sData = readData(system.file("extdata", "s.xlsx", package = "fitPS"))
```

CSV files with the same two-column layout can also be read:

```r
pData = readData("survey.csv")
```


## Zeta shape parameterisation

fitPS uses `shape` for the zeta distribution shape parameter, with `shape > 1`.
Users should supply, inspect, compare, and report `shape` on this scale.

## Fitting a zeta model

```r
fit = fitDist(pData)
fit
```

`fitDist()` returns a `psFit` object. Standard methods include printing, plotting, prediction, fitted values, confidence intervals, and log-likelihood extraction.

```r
logLik(fit)
predict(fit)
```

## Fitting a zero-inflated zeta model

```r
ziFit = fitZIDist(pData)
ziFit
predict(ziFit)
```

## Probability functions

Use `probfun()` to create a function for computing fitted `P` or `S` probabilities.

```r
pFun = probfun(fit)
pFun(0:5)
```

## Comparing surveys

The package includes survey-comparison helpers:

```r
compareSurveys(fit1, fit2)
compareSurveysLRT(data1, data2)
```

See the function documentation and vignettes for details.

## Vignettes

The `vignettes/` directory contains worked examples for simple fitting and confidence-region workflows. During staged package checks, built vignette artifacts are excluded so strict checks remain focused on package code, documentation, and tests.

## Development checks

The stage 1 stabilization work added baseline `testthat` coverage and strict package checks. During development, use:

```r
devtools::document()
devtools::test(stop_on_failure = TRUE, stop_on_warning = TRUE)
devtools::check(
  build_args = c("--no-build-vignettes"),
  args = c("--no-manual", "--ignore-vignettes", "--no-tests"),
  error_on = "note"
)
```

## License

GPL (>= 2).
