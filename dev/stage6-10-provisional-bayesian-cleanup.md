# Stage 6.10 - provisional Bayesian cleanup

Stage 6.10 removes the duplicated Bayesian fitting paths that remained after the
model/engine migration. The public fitting entry points are now `fitDist()` and
`fitZIDist()`; both construct a model, select a posterior engine, and call the
same internal `fitBayesianModel()` orchestration.

The method-specific development functions `fitDistBayes()`,
`fitDistBayesIntegrate()`, `fitZIDistBayes()`, `fitZIDistBayesNumerical()`,
`fitZIDistBayesLaplace()`, and `fitZIDistBayesImportance()` have been removed.
They were provisional 1.0.7 implementation details rather than part of the
established CRAN 1.0.6 interface.

`psPosterior` is now the source of truth for Bayesian posterior summaries,
representations, and diagnostics. ZIZ fits no longer duplicate grids, chains,
Laplace objects, importance samples, diagnostics, posterior summaries, or
posterior representations at the top level of `psFit`. Public posterior methods
read from `psPosterior` instead.

The older plain-zeta Bayesian fields needed for the established pre-1.0.7 API
(`var.shape`, MCMC `chain`, and marginal `pdf`) remain available. Their
construction is isolated in `establishedBayesianFitFields.zetaModel()` so that
compatibility does not leak back into the engine implementations.

The four separate inflation-probability calculation helpers have also been
removed. Their small engine-specific calculations now live directly in the S3
methods for `posteriorInflationProbability()`, which is the polymorphic interface
introduced during Stage 6.

The unused `attachZizPosterior()`, `makeZizParameterSummary()`,
`makeZizMarginalPdf()`, and the old proof-of-concept `fitdistLaplace()` helper are
removed as dead migration-era code.

Stage 6.11 should perform the final regression and architecture audit, including
strict package validation and a search for remaining dead compatibility code or
avoidable model/engine duplication.
