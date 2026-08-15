# Stage 9.7 source-layout consolidation

Stage 9.7 is a behaviour-preserving source-layout cleanup following the Stage 9 Bayesian engine work.

## Purpose

The package had accumulated several very small R files that each contained one method or a small wrapper from the same conceptual API. This made the source tree harder to scan without creating useful architectural boundaries.

Stage 9.7 consolidates only those clearly fragmented areas. Function bodies and roxygen documentation are moved without changing their public interfaces or statistical behaviour.

## Consolidated files

The following psData fragments are moved into `R/psData-methods.R`:

- `R/as.data.frame.psData.R`
- `R/mean.psData.R`
- `R/var.fitPS.R`
- `R/print.psData.R`
- `R/operators.psData.R`
- `R/add.R`

The following small psFit methods are moved into `R/psFit-methods.R`:

- `R/fitted.psFit.R`
- `R/logLik.psFit.R`
- `R/print.psFit.R`
- `R/summary.psFit.R`

The survey-comparison family is moved into `R/survey-comparison.R`:

- `R/compareSurveys.R`
- `R/compareSurveysLRT.R`
- `R/fitCompare.R`

The small bootstrap user-facing wrappers are moved into `R/bootstrap-inference.R`:

- `R/bootstrapFit.R`
- `R/bootstrapProbs.R`

## Deliberately unchanged

This stage does not split or rewrite `R/model-zeta.R` or `R/model-ziz.R`. Those files still contain specialised Bayesian code whose final ownership depends on a later decision about whether the built-in MCMC implementations should remain for exact compatibility or be replaced by the generic MCMC engine.

The shared posterior-engine files, Bayesian model contract, posterior orchestration, prediction, confidence-interval, plotting, and substantive bootstrap implementation files remain separate because they represent useful architectural or user-facing boundaries.

## Validation intent

The full existing test suite is the regression contract for this stage. Roxygen documentation is regenerated to verify that moving the source blocks does not alter exports, S3 registration, aliases, or generated Rd documentation.
