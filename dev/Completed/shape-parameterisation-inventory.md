# fitPS Stage 2.1 shape parameterisation inventory

## Purpose

This inventory records where the zeta shape parameterisation appears before the Stage 2 refactor. It is intended to guide the follow-up implementation stage and to make the current VGAM-shifted assumptions visible.

## Package convention to enforce

- `shape` is the user-facing and package-facing zeta parameter, equal to the standard Riemann zeta exponent alpha.
- `shape` must be greater than 1.
- VGAM uses a shifted parameter internally, equal to `shape - 1`.
- Calls to VGAM zeta probability or random-generation helpers should pass a clearly named shifted value such as `vgamShape`.

## Inventory summary

### User-facing `shape = alpha`

These locations expose or describe `shape` as part of the package API and should use the standard zeta parameter after the refactor:

- `R/fitDist.R`: return value `shape`, documented fitting parameter, start/legacy `shape` argument handling.
- `R/fitZIDist.R`: return value `shape`, documented fitting parameter, start/legacy `shape` argument handling.
- `R/probfun.R`: probability helper for fitted zeta and zero-inflated zeta models.
- `R/predict.psFit.R`: S3 prediction method for fitted models.
- `R/fitted.psFit.R`: S3 fitted method through `probfun()` when `n` is supplied.
- `R/rzeta.R`: exported random generation helper.
- `R/rZIzeta.R`: exported zero-inflated random generation helper.
- `R/print.psFit.R`, `R/summary.psFit.R`, `R/confint.psFit.R`, `R/bootCI.R`, `R/credint.R`, `R/compareSurveys.R`, and `R/compareSurveysLRT.R`: reported, compared, intervalled, or resampled shape values.

### VGAM shifted-shape calculations

These locations currently pass `shape` directly to VGAM and therefore need review in Stage 2.2:

- `R/fitDist.R`: likelihood, fitted values, variance calculation comments, and stored `fit$par` currently treat `shape` as VGAM's shifted parameter.
- `R/fitDistBayes.R` and `R/fitDistBayesIntegrate.R`: Bayesian likelihood and fitted values pass `shape` directly to VGAM.
- `R/fitZIDist.R`: likelihood, fitted values, helper density, and stored `fit$par[2]` currently treat `shape` as VGAM's shifted parameter.
- `R/fitZIDistBayes.R`: Bayesian zero-inflated likelihood and fitted values need the same review.
- `R/probfun.R`, `R/predict.psFit.R`, `R/rzeta.R`, and `R/rZIzeta.R`: exported helpers call VGAM directly using package `shape`.
- `R/confint.psFit.R`, `R/profileLikelihoodZIZ.R`, and any interval/profile calculations that call `VGAM::dzeta()` should convert at the VGAM boundary.

### Documentation-only references

These references should be updated after the code refactor so they describe `shape = alpha > 1` and mention VGAM's shifted parameter only as an implementation detail:

- `R/fitDist.R` and `R/fitZIDist.R` roxygen notes currently describe stored values as VGAM parameterisation.
- `R/rzeta.R` and `R/rZIzeta.R` roxygen parameter text currently permits shifted values greater than zero.
- `R/confint.psFit.R` roxygen currently says intervals are for VGAM parameterisation.
- Generated `man/*.Rd` files mirror the same roxygen text and should be regenerated with `devtools::document()`.
- `README.md` and `NEWS.md` should be reviewed for any user-visible shape wording.

## Stage 2.1 tests added

`tests/testthat/test-zeta-shape-parameterisation.R` adds desired-behavior tests for the standard-shape convention. These tests are intentionally written against the target behaviour rather than the current implementation. They are expected to expose failures until Stage 2.2 converts package code to use `vgamShape = shape - 1` at VGAM boundaries.

## Follow-up implementation notes

- Start by adding small helpers or local variables named `vgamShape` where VGAM is called.
- Preserve the public `shape` name and returned list element names.
- Treat `shape <= 1` as invalid for user-facing arguments.
- Be cautious with intervals, bootstraps, comparison methods, saved fit objects, and Bayesian methods because they may store or resample the shifted value today.
