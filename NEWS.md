# fitPS development news

## Purpose

This file records user-facing and developer-facing changes for fitPS. It is a release-note summary, not a commit-by-commit history.



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
