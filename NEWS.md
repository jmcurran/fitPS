# fitPS development news

## Purpose

This file records user-facing and developer-facing changes for fitPS. It is a release-note summary, not a commit-by-commit history.






















## fitPS 1.0.9.007

- Completed the Stage 7 model-comparison and distribution-extensibility audit across zeta, ZIZ, and logarithmic models.
- Corrected the cross-model audit test to accept the numeric nobs attribute returned by the shared logLik contract.
- Added durable cross-model regression tests for the shared MLE comparison contract and intentional posterior-engine support matrix.
- Confirmed that logarithmic integration reuses the model/engine abstractions without distribution-name branches in generic posterior engines.
- Documented that the current extension contract remains internal and recommended a later deliberate third-party model API rather than exporting evolving internals piecemeal.
- Recorded a future thin fit(data, model = ...) front end as the natural place to accept external model objects while retaining existing fitters.
- Validated with vignette rendering, roxygen, strict tests, and strict package check.

## fitPS 1.0.9.005

- Added the developer-facing Extending fitPS with a New Distribution vignette using the logarithmic-series model as the worked example.
- Documented the psModel extension contract for parameters, probabilities, likelihoods, supported posterior engines, fitting, posterior summaries, and model comparison.
- Documented the expected behavioural tests for new distributions and the current boundary between internal package extensibility and a future public model-registration API.
- Recorded the possible future fit(data, model = ...) front end while retaining fitDist() and fitZIDist() as compatibility and convenience interfaces.
- Validated the vignette by rendering it in addition to the strict fitPS test and package-check workflow.

## fitPS 1.0.9.004

- Added the logarithmic-series distribution as a psModel with shared P/S support, probability, and likelihood dispatch.
- Added numerical and MCMC Bayesian logarithmic fits through the existing posterior-engine extension points and common psPosterior finalisation path.
- Enabled logarithmic MLE fits to participate in logLik, deviance, AIC, BIC, fitted probabilities, prediction, and standard psFit summaries.
- Enabled DIC for one-parameter numerical and MCMC models by removing residual zeta-specific parameter-name assumptions.
- Generalised uniform/custom prior ranges while preserving zeta-domain safeguards for zeta and ZIZ posterior fits.
- Replaced the provisional logarithmic fitter internals rather than treating them as a compatibility contract.
- Validated by the strict fitPS stage workflow.

## fitPS 1.0.9.002

- Added a shared model-comparison contract for likelihood, deviance, and parameter counting.
- Enabled AIC and BIC for maximum-likelihood psFit objects through the common logLik interface.
- Added DIC for numerical, MCMC, and importance Bayesian posterior representations; Laplace fits fail explicitly until a posterior deviance expectation is available.
- Preserved the Bayesian-bootstrap proposal as future work under a non-stage-specific context and recorded the future direction toward a model-oriented fit interface.
- Corrected the DIC contract test to compare the numeric criterion separately from its Dbar, Dhat, and pD attributes.
- Validated by the strict fitPS stage workflow.

## fitPS 1.0.8.014

- Completed the final Bayesian architecture and regression audit for the model/engine migration.
- Added durable tests for the supported model/engine matrix, removal of provisional fitter paths, and model-descriptor extensibility.
- Confirmed all named R functions remain documented and removed the last migration-stage wording from internal roxygen text.
- Recorded logarithmic-distribution extensibility and AIC/DIC model comparison as post-migration follow-on work.
- Re-ran the full fitPS validation workflow before closing Stage 6.

## fitPS 1.0.8.013

- Completed the provisional Bayesian cleanup and repaired controlled Git staging for deleted generated documentation.
- Kept psPosterior as the source of truth for ZIZ posterior representations, probability summaries, and diagnostics.
- Preserved established pre-1.0.7 plain-zeta Bayesian compatibility fields behind the model-specific compatibility method.
- Carried forward the full Stage 6.10 cleanup because the failed attempt did not commit.
- Re-ran the full fitPS validation workflow before committing the repaired cleanup.

## fitPS 1.0.8.012

- Removed superseded provisional Bayesian fitter wrappers and routed fitDist() and fitZIDist() through the shared model/engine orchestration.
- Made psPosterior the source of truth for ZIZ posterior representations, probability summaries, and diagnostics instead of duplicating native engine payloads at psFit top level.
- Preserved established pre-1.0.7 plain-zeta Bayesian compatibility fields behind a model-specific compatibility method.
- Removed dead migration helpers and consolidated inflation calculations directly into representation-specific S3 methods.
- Validated with the full fitPS stage workflow because package R source, S3 registration metadata, generated documentation, and tests changed.

## fitPS 1.0.8.011

- Consolidated Bayesian psFit finalisation and psPosterior construction across plain zeta and zero-inflated zeta models.
- Added the common psPosterior contract and posterior mean probability summaries to plain-zeta numerical and MCMC fits.
- Standardised psPosterior representation storage on typed engine wrappers while retaining legacy top-level representation fields during migration.
- Replaced inflation method-name switching with S3 dispatch on posterior representation classes and centralised derived probability summarisation.
- Validated with the full fitPS stage workflow because package R source, S3 registration metadata, and tests changed.

## fitPS 1.0.8.010

- Migrated the existing zero-inflated zeta Laplace and importance-sampling approximations onto the shared model and posterior-engine architecture.
- Preserved the transformed-coordinate Laplace optimisation, Hessian/delta-method covariance, importance proposal, seeded weighted sampling, posterior summaries, and fitted probability behaviour.
- Added engine-level fitting, summary, diagnostic, and point-estimate methods while keeping plain-zeta Laplace and importance combinations explicitly unsupported.
- Added deterministic P and S regression coverage using stable behaviour-based test names and descriptions.
- Validated with the full fitPS stage workflow because package R source, S3 registration metadata, and tests changed.

## fitPS 1.0.8.009

- Completed the zero-inflated zeta MCMC migration onto the shared model and posterior-engine architecture while preserving the established Metropolis-Hastings sampler and seeded RNG sequence.
- Registered the internal Stage 6 S3 methods explicitly so roxygen documentation generation succeeds without exposing those helpers as user-facing API.
- Preserved the historical ZIZ MCMC fitted-probability component as a named numeric vector rather than a one-row matrix.
- Kept the behaviour-based test filename cleanup and strengthened deterministic P and S regression coverage for the migrated ZIZ MCMC path.
- Validated with the full fitPS stage workflow because package R source, S3 registration metadata, and tests changed.

## fitPS 1.0.8.007

- Migrated numerical zero-inflated zeta Bayesian fitting onto the shared Stage 6 model and posterior-engine protocols while preserving the existing two-dimensional posterior grid calculation.
- Added the ZIZ model log-likelihood method and numerical-engine model dispatch while reusing the shared numerical summary, diagnostics, and point-estimate methods.
- Retained the legacy posterior grid, marginal densities, posterior probability summaries, parameter moments, fitted probabilities, and psPosterior representation while also storing the engine wrapper.
- Added deterministic P and S regression tests against the pre-migration numerical ZIZ orchestration; ZIZ MCMC, Laplace, and importance fitting remain unmigrated.
- Validated with the full fitPS stage workflow because package R source and tests changed.

## fitPS 1.0.8.006

- Migrated the plain-zeta MCMC Bayesian fit onto the Stage 6 model and posterior-engine protocols while preserving the existing sampler calculation and RNG ordering.
- Added MCMC-engine fit, summary, diagnostics, point-estimate, and representation methods for zetaModel.
- Retained the legacy chain, posterior-mean shape, variance, fitted probabilities, and spline density while also storing the engine representation.
- Added deterministic P and S regression tests against the pre-refactor MCMC algorithm; ZIZ MCMC fitting remains explicitly unmigrated.
- Validated with the full fitPS stage workflow because package R source and tests changed.

## fitPS 1.0.8.005

- Migrated the plain-zeta numerical Bayesian fit onto the Stage 6 model and posterior-engine protocols while preserving the existing integration calculation.
- Added the zeta model log-likelihood method and numerical-engine fit, summary, diagnostics, point-estimate, and representation methods.
- Kept fitted P/S probabilities evaluated at the posterior mean shape and retained the legacy numerical posterior density component for compatibility.
- Added deterministic regression tests against the pre-refactor numerical algorithm for both P and S data; ZIZ numerical fitting remains explicitly unmigrated.
- Validated with the full fitPS stage workflow because package R source and tests changed.

## fitPS 1.0.8.004

- Added internal S3 psModel descriptors for plain zeta and zero-inflated zeta models.
- Added posterior-engine descriptors and common fit, summary, diagnostics, point-estimate, and representation protocols without migrating existing fitting paths yet.
- Declared the supported model/engine matrix explicitly so unsupported combinations fail clearly rather than being inferred from method switches.
- Consolidated duplicated P/S probability-index validation, latent-support conversion, term naming, and observation mapping shared by zeta and ZIZ code.
- Added deterministic offline protocol tests; validated with the full fitPS stage workflow because package R source and tests changed.

## fitPS 1.0.8.003

- Completed the Stage 6.2 internal roxygen2 documentation baseline and Bayesian rationale comments.
- Repaired shared summary help topics by documenting the x arguments contributed by print.summary.psBootstrap() and print.summary.psPosterior().
- Preserved executable R statements while carrying forward the complete Stage 6.2 documentation change set after its pre-commit check failure.
- Added the Stage 6.2.1 repair note explaining the roxygen @describeIn argument-documentation issue.
- Validated with the full fitPS stage workflow because package R source files changed.

## fitPS 1.0.8.001

- Began Stage 6 as the fitPS Bayesian architecture, documentation, and consolidation series.
- Added the Stage 6.1 audit defining orthogonal S3 model and posterior-engine abstractions for plain zeta, zero-inflated zeta, and future distributions.
- Recorded consolidation targets for shared observation handling, probability naming, prior validation, posterior summaries, and Bayesian psFit construction.
- Reserved the 1.0.8.xxx version series for Stage 6, with every runner attempt consuming a new build number.
- Renumbered the proposed Bayesian-bootstrap context as future Stage 7.1 work so it remains separate from the architecture refactor.
- Stage 6.1.1 repairs the Stage 6.1 delivery workflow and installs the same planning documentation; no R package validation is claimed for this planning stage.

## fitPS 1.0.7.018

- Added posteriorInflation() for Pr(pi < epsilon | data) under Bayesian zero-inflated zeta fits.
- Implemented the practical-inflation probability consistently for numerical, MCMC, importance, and Laplace posterior engines.
- Added inflationEpsilon to summary.psPosterior(), defaulting to 0.01, and report the resulting practical-negligibility probability.
- Replaced the vignette inflation discussion with an epsilon-based interpretation using the Roux footwear example and epsilon = 0.01.
- Documented that epsilon is application-specific and that the diagnostic differs from assigning posterior mass to the exact no-inflation model.
- Added deterministic engine-specific tests and a Roux numerical-posterior regression test.
  - Repaired summary compatibility for posterior objects without a model tag, corrected Roux grid-option placement, and updated psPosterior structure expectations.
- Validated through the full fitPS documentation, test, check, vignette, build, and installation workflow.
## fitPS 1.0.7.016

- Corrected the Roux vignette description to identify the response as the number of different glass sources found on each surveyed pair of shoes.
- Added posterior interpretation of the zero/one-inflation parameter, including the distinction between concentration near zero and formal posterior probability of a no-inflation model.
- Added the Efron bootstrap/Bayesian connection to the bootstrap discussion while preserving the distinct inferential interpretations.
- Added future-development context for Rubin's Bayesian bootstrap, centred on weighted-likelihood support and a possible psBayesianBootstrap object.
- Added long-term architecture context for distribution-pluggable fitPS models and evaluation of the distributional package as a possible distribution representation layer.
- Validated through the full fitPS documentation, test, check, vignette, build, and installation workflow.
## fitPS 1.0.7.015

- Added a dedicated Bayesian and bootstrap probability vignette using the real Roux et al. footwear survey supplied with fitPS.
- Documented plug-in, bootstrap-mean, and posterior-mean P probabilities and their distinct frequentist and Bayesian interpretations.
- Documented the numerical, MCMC, Laplace, and importance posterior engines, including priors, importance weighting, and Laplace simulation for derived probabilities.
- Added Roux-based README examples for posteriorProbs(), bootstrapProbs(), fitted(), plotting, and bootstrapFit().
- Expanded Bayesian probability help and fitZIDist() cross-references to the public uncertainty APIs.
- Added an end-to-end real-data regression test covering MLE, numerical Bayesian, bootstrap, fitted probabilities, interval bounds, and legacy probfun() behaviour.
- Added the Stage 5.11 final regression audit documenting completed interfaces and optional future cleanup only.
- Validated through roxygen, strict offline tests, strict package check, vignette build, source build, and installation.
## fitPS 1.0.7.014

- Harmonised predict.psFit() across plug-in, posterior-mean, and bootstrap-mean probability definitions while preserving plug-in predictions as the default.
- Added equal-tailed credible prediction intervals for stored posterior summaries and percentile confidence intervals for stored bootstrap summaries.
- Required posterior and bootstrap interval levels to match the level stored in the corresponding uncertainty object, avoiding unsupported implicit recomputation.
- Updated print.psFit() to route fitted probabilities through fitted(), fixing the legacy long-output path that could ignore zero inflation.
- Added concise posterior and bootstrap printing with optional nterms selection and bootstrap-aware psFit summaries.
- Clarified that probfun() evaluates plug-in probabilities and directs users to posteriorProbs() or bootstrapProbs() for uncertainty-distribution means.
- Added deterministic offline regression tests for prediction semantics, intervals, term selection, printing, and summary presentation.
- Validated through roxygen, strict offline tests, and strict package check.
## fitPS 1.0.7.012

- Added bootstrapFit() as the public entry point for attaching a psBootstrap distribution to an MLE fit.
- Added bootstrapProbs() methods for psFit and psBootstrap objects with the same P/S term selection rules as posteriorProbs().
- Added fitted(..., type = "bootstrapMean") while preserving plug-in fitted probabilities as the default.
- Added plot.psBootstrap() for bootstrap mean probabilities and percentile confidence intervals.
- Kept Bayesian posteriorMean fitted values and all existing frequentist defaults unchanged.
- Added deterministic offline API, indexing, fitted-value, plotting, and validation tests.
- Validated through roxygen, strict offline tests, and strict package check.
## fitPS 1.0.7.011

- Added an internal psBootstrap S3 object for frequentist bootstrap parameter and P/S probability distributions.
- Added a reusable nonparametric bootstrap replicate engine with seeded reproducibility and failed-fit diagnostics.
- Added bootstrap mean, standard deviation, and percentile confidence intervals for parameter and probability replicates.
- Reused the existing zero-inflated probability transformation and added the matching standard-zeta P/S transformation.
- Kept bootCI() return structures backward compatible while routing its replicate generation through the shared engine.
- Added deterministic offline tests for zeta and zero-inflated bootstrap summaries, reproducibility, failures, and legacy bootFit() shapes.
- Validated through roxygen, strict offline tests, and strict package check.
## fitPS 1.0.7.009

- Added the Stage 5.7 planning audit for bootstrap distributions of fitted P and S probabilities.
- Recommended a psBootstrap S3 object that parallels psPosterior while preserving the distinct frequentist interpretation.
- Documented reuse of retained bootstrap parameter replicates, zizProbabilities(), percentile probability intervals, failure diagnostics, and seeded serial reproducibility.
- Preserved bootCI() return values and ordinary plug-in fitted probabilities as compatibility constraints for the implementation stages.
- Planned Stages 5.8-5.11 for bootstrap objects, public bootstrap probability APIs, harmonised prediction, and final documentation.

## fitPS 1.0.7.008

- Added plot.psPosterior() for posterior mean P or S probabilities with equal-tailed credible intervals.
- Supported plotting all stored terms, leading terms, or explicitly selected P and S indices through posteriorProbs().
- Kept plot.psFit() and plotPosterior() behaviour unchanged for backward compatibility.
- Added deterministic offline plotting tests for P and S posterior objects and interval suppression.
- Validated through roxygen, strict offline tests, and strict package check.
## fitPS 1.0.7.007

- Exported posteriorProbs() for Bayesian psFit and psPosterior objects.
- Added consistent selection of leading or explicitly indexed P and S posterior probability summaries.
- Extended fitted.psFit() with type = "posteriorMean" while preserving plug-in probabilities as the default.
- Kept posterior probability access unavailable for frequentist fits with clear errors.
- Added deterministic offline tests for S3 dispatch, subset selection, fitted-value semantics, and compatibility.
- Validated through roxygen, strict offline tests, and strict package check.
## fitPS 1.0.7.006

- Added the psPosterior S3 class as the coherent posterior component of Bayesian psFit objects.
- Stored parameter summaries, posterior probability summaries, engine-specific representations, interval level, and diagnostics under fit$posterior.
- Preserved chain, posteriorGrid, weightedSamples, laplace, importance, and posteriorProbs fields for backward compatibility.
- Added print, summary, and fitted methods for psPosterior objects and delegated Bayesian psFit summaries to the posterior object.
- Corrected the Stage 5.3 tests by comparing probability estimates without incidental names and by using valid compact S-survey data.
- Validated through roxygen, strict offline tests, and strict package check.
## fitPS 1.0.7.004

- Repaired the existing parallel import documentation in bootCI so every @importFrom tag occupies one physical line.
- Removed the roxygen warning reported during Stage 5.2 documentation regeneration without changing package behaviour.
- Scanned the remaining R sources and found no other wrapped @importFrom continuation lines.
- Validated through roxygen, strict offline tests, and strict package check.
## fitPS 1.0.7.003

- Added one shared internal transformation from zero-inflated zeta parameters to requested P or S probabilities.
- Used the package standard zeta shape parameterisation and supported vectorised posterior parameter values.
- Added deterministic tests for P and S formulas, names, validation, scalar recycling, and probability-mass truncation.
- Removed and ignored generated Rplots.pdf test artifacts.
- Corrected stage-runner source-package path capture so build console output cannot be mistaken for a filename.
- Validated through roxygen, strict offline tests, and strict package check.
## fitPS 1.0.7.002

- Repaired the zero-inflated zeta Metropolis-Hastings sampler so each component update evaluates only the current state and its proposed component.
- Implemented the beta-prior and independence-proposal cancellation explicitly for pi updates and used the selected psPrior for shape updates.
- Preserved legacy a/b shape-prior bounds when no explicit prior is supplied and added seeded reproducibility support.
- Added non-uniform Beta(2, 5) tests comparing MCMC posterior means and variances with deterministic numerical integration.
- Updated the stage runner delivery to use a single outer bundle with an exact sibling change-set and optional ChatGPT bundle creation.
- Validated through roxygen, strict offline tests, and strict package check.
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
