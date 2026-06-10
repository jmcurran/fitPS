# fitPS development news

## Purpose

This file records user-facing and developer-facing changes for fitPS. It is a release-note summary, not a commit-by-commit history.










## fitPS 1.0.6.013

- Rationalise the Bayesian fitting API around method = "bayes" and bayesOptions$posteriorMethod.
- Add deprecated legacy aliases for numerical integration, MCMC, Laplace, and importance posterior methods.
- Preserve compatibility by translating legacy Bayesian method values into canonical bayesOptions.
- Add tests for Bayesian method alias translation and default option handling.
- Validated by the Stage 4.6.1.1 full package runner.

## fitPS 1.0.6.011

- Added an internal importance-sampling helper for zero-inflated zeta Bayesian posteriors using the Laplace approximation as a proposal.
- Returned weighted samples, posterior means, covariance estimates, and diagnostics including effective sample size and maximum normalized weight.
- Added deterministic offline tests for weighted samples, seeded reproducibility, and fitZIDist dispatch through bayesOptions$posteriorMethod = "importance".
- Added ignore/build-ignore coverage for generated built-package path files so package checks do not inspect runner artifacts.
- Ran the full package validation workflow for this package-impacting stage.

## fitPS 1.0.6.010

- Repaired the Stage 4.4 numerical posterior grid file by restoring the missing bracket in the returned pi value.
- Kept the Stage 4.4 Laplace approximation work unchanged while making the installed R code parsable for roxygen documentation.
- Preserved ignore/build-ignore coverage for generated built-package path files so package checks do not inspect runner artifacts.
- Ran the full package validation workflow for this package-impacting repair stage.

## fitPS 1.0.6.006

- Added Bayesian options plumbing so method = "bayes" can select posterior engines through bayesOptions$posteriorMethod.
- Stripped inherited names from scalar transform inputs and outputs before constructing working-scale vectors.
- Routed fitDist() numerical integration through method = "bayes" with posteriorMethod = "numerical" while retaining method = "integrate" as a legacy alias.
- Kept zero-inflated zeta MCMC available through posteriorMethod = "mcmc" until the numerical Stage 4 posterior is implemented.
- Kept deterministic offline test coverage unchanged while repairing the helper implementation.

## fitPS 1.0.6.001

- Added the Stage 4.1 design audit for deterministic and semi-deterministic Bayesian posterior approximations for zero-inflated zeta fits.
- Recorded the recommended path for transformed posterior helpers, grid quadrature, Laplace diagnostics, importance sampling diagnostics, and plotting integration.
- Preserved existing Bayesian MCMC behavior and user-facing defaults; this lightweight documentation stage did not run package validation.

## fitPS 1.0.5.008

- Stabilised the generated stage runner so controlled paths are staged one at a time.
- Prevented absent optional directories from masking `git add` failures and leaving new files untracked.
- Preserved the controlled commit scope while keeping the workflow compatible with package roots that do not contain every optional directory.
- Kept the no-ChatGPT-bundle workflow and completed-stage archive cleanup behaviour.
- Validated by the strict fitPS stage workflow.

## fitPS 1.0.5.007

- Added `plotPosterior()` for posterior-density plots from Bayesian `psFit` objects.
- Supported zeta MCMC chains, zeta numerical posterior density functions, and zero-inflated MCMC chains for `shape` and `pi`.
- Kept existing `plot.psFit()` fitted-probability plots unchanged for backward compatibility.
- Added deterministic offline tests using small fake `psFit` objects rather than long MCMC runs.
- Validated by the strict fitPS stage workflow.

## fitPS 1.0.5.003

- Added Stage 3.2 workflow documentation for the Windows Ghostscript `R CMD check` NOTE under `dev/`.
- Kept strict package validation while allowing only the known missing-Ghostscript PDF-size reduction NOTE in the stage runner.
- Documented why the Ghostscript NOTE is a local check-environment issue rather than a package failure.
- Preserved the Stage 3 `1.0.5.xxx` build convention, including one consumed build number per attempt.
- Validated by the strict fitPS stage workflow.

## fitPS 1.0.5.002

- Added the Stage 3.1 Bayesian posterior plotting audit under `dev/`.
- Recommended a separate `plotPosterior()` API for posterior parameter-density plots while preserving `plot.psFit()` for fitted probabilities.
- Documented support requirements for MCMC-backed zeta and zero-inflated zeta fits, integrated zeta posterior densities, credible intervals, deterministic tests, and shape-parameterisation safeguards.
- Reset the Stage 3 development series to the `1.0.5.xxx` build convention.
- Validated by the strict fitPS stage workflow.

## fitPS 1.0.4.018

- Migrated the user-facing vignette sources from Sweave `.Rnw` files to R Markdown `.Rmd` files.
- Removed generated vignette build artifacts from the staged source tree so they can be regenerated by the package toolchain.
- Added `rmarkdown` to `Suggests` for the R Markdown vignette builder.
- Preserved the Stage 2 zeta shape wording while making the vignettes easier to maintain.

## fitPS 1.0.4.017

- Removed user-facing discussion of VGAM's shifted zeta parameterisation from help and vignette source text.
- Kept documentation focused on the fitPS shape parameter, where shape > 1.
- Corrected spelling, grammar, and punctuation issues in the vignette source using NZ English conventions.
- Regenerated help with devtools::document() during the strict stage workflow.

## fitPS 1.0.4.016

- Cleaned roxygen comments and README wording for the Stage 2 zeta shape parameterisation documentation.
- Fixed spelling, punctuation, and grammar in the fitDist(), fitZIDist(), and predict.psFit() help source.
- Kept documentation on NZ English wording, including parameterisation and optimisation.
- Left generated help files to be refreshed by devtools::document() during the stage workflow.

## fitPS 1.0.4.015

- Updated documentation so `shape` consistently refers to the standard zeta parameter alpha with `shape > 1`.
- Removed obsolete VGAM-shifted wording from `fitDist()` and `fitZIDist()` help text.
- Documented that the default Bayesian prior is placed on `log(shape - 1)`, giving support only on valid standard-shape values.
- Added README guidance for users who compare fitPS results with VGAM zeta functions.
- Validated by the strict fitPS stage workflow.
## fitPS 1.0.4.014

- Refactored zeta and zero-inflated zeta workflows so fitPS `shape` means the standard zeta parameter alpha with `shape > 1`.
- Converted to VGAM's shifted parameter only at VGAM boundaries using `shape - 1`.
- Added active regression tests for fitting, probability functions, prediction, fitted values, and random generation under the standard-shape convention.
- Validated by the strict fitPS stage workflow.

## fitPS 1.0.4.012

- Added Stage 2.1.1 corrective notes for the zeta shape parameterisation workflow.
- Removed active target-behaviour tests from `tests/testthat` because they intentionally fail before the Stage 2.2 refactor.
- Preserved the standard-shape expectations in `dev/` so Stage 2.2 can reinstate them with the implementation repair.
- Validated by the strict fitPS stage workflow.

## fitPS 1.0.4.010

- Added README.md with installation, data-format, fitting, prediction, comparison, and development-check examples.
- Updated `.gitignore` to exclude the fitPS built-package path marker produced by stage runners.
- Kept the stage focused on documentation and local-build hygiene without changing package APIs.
- Validated by the strict stage 1.4 package workflow before commit.


## fitPS 1.0.4.009

- Added NEWS.md as durable package release notes rather than a stage transcript.
- Recorded the stage 1 stabilization series, including baseline tests, build hygiene, start-value handling, and zero-inflated prediction fixes.
- Kept entries version-numbered and newest-first for future fitPS maintenance.
- Validated by the strict stage 1.3 package workflow before commit.


## fitPS 1.0.4 development series

- Added baseline `testthat` infrastructure for core `psData` and `psFit` workflows.
- Added deterministic tests for data construction, CSV import, model fitting, prediction, and `logLik.psFit()`.
- Repaired stage check issues by documenting the `fitDist()` `prior` argument and adding required `stats` imports.
- Cleaned package-build inputs by excluding local development files, generated archives, stage runners, scripts, and generated vignette artifacts from R package builds.
- Repaired start-value handling for `fitDist()` and `fitZIDist()` while preserving backward-compatible `shape` aliases.
- Repaired zero-inflated zeta prediction logic so fitted predictions use the model object's `pi` component.
- Validated through the staged fitPS workflow before completed stage commits.
