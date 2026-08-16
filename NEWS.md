# fitPS development news

## Purpose

This file records user-facing and developer-facing changes for fitPS. It is a release-note summary, not a commit-by-commit history.









## fitPS 1.1.3.019

- Corrected the four-panel uncertainty vignette so shared axes include the final bootstrap KDE contour coordinates as well as the realization clouds, preventing contour clipping.
- Kept the four ZIZ uncertainty panels on exactly the same padded pi and shape scales for direct visual comparison.
- Hid panel-layout, KDE-boundary, and shared-axis bookkeeping from the rendered vignette while leaving only the four meaningful plotUncertainty() calls visible to readers.
- Preserved the established profile-likelihood, ordinary-bootstrap, parametric-Bayesian, and Rubin Bayesian-Bootstrap interpretations and plotting behavior.
- Validated through the strict fitPS documentation, test, package-check, vignette-build, build, and installation workflow.

## fitPS 1.1.3.018

- Made the uncertainty vignette compare MLE, ordinary bootstrap, parametric Bayes, and Rubin Bayesian Bootstrap ZIZ regions in a single 2 by 2 panel with shared parameter axes.
- Showed ordinary- and Bayesian-Bootstrap parameter realizations beneath their KDE contours by default while leaving parametric Bayesian point displays opt-in.
- Preserved the original ZIZ smoothed-bootstrap construction by using Hscv() KDE smoothing with positivity enforced for the pi/shape parameterization.
- Standardized user-facing terminology on profile-likelihood confidence intervals and regions while retaining the likelihood-ratio cutoff as their mathematical construction.
- Added regression coverage for method-specific bootstrap realization defaults and validated through the strict fitPS package workflow.

## fitPS 1.1.3.017

- Made plotUncertainty() the documented parameter-uncertainty plotting interface across MLE, ordinary bootstrap, parametric Bayes, and Rubin Bayesian Bootstrap fits.
- Kept credint()/credInt() as the canonical Bayesian credible-interval and credible-region extractor, extending it across stored posterior representations and Rubin Bayesian Bootstrap results.
- Kept confint() as the frequentist numerical confidence-interval and likelihood-ratio contour extractor, with credint() as its Bayesian counterpart and plotUncertainty() as the common visual layer.
- Updated the basics, uncertainty, and README documentation so plotUncertainty() examples are rendered and the numerical-versus-visual uncertainty interfaces are explicit.
- Deprecated only the redundant bootCI() compatibility path and added regression coverage that credint() remains an active Bayesian extractor.
- Validated with the strict fitPS documentation, test, package-check, build, and installation workflow.

## fitPS 1.1.3.016

- Repaired the cumulative posterior integration helper to use direct vector indexing instead of unqualified head() and tail() calls.
- Removed the resulting R CMD check NOTE without adding unnecessary package imports or changing numerical results.
- Carried forward the Stage 11.7 uncertainty-plotting interface, Simpson cumulative integration, and weighted KDE repair unchanged.
- Validated with the strict fitPS documentation, test, package-check, build, and installation workflow.

## fitPS 1.1.3.013

- Corrected plotPosterior() credible-region hatching so the shaded polygon follows the posterior density curve rather than joining only the interval endpoint heights.
- Added deterministic regression coverage for the posterior interval polygon geometry.
- Modernized plotPosterior() documentation and examples around the fit() Bayesian interface.
- Clarified why the uncertainty vignette uses the zero-inflated zeta model as a rich worked example and stated that model choice is a separate task.
- Validated with the strict fitPS documentation, test, and package-check workflow.

## fitPS 1.1.3.012

- Rewrote README.md around fit() as the canonical interface for MLE, parametric Bayes, the ordinary bootstrap, and Rubin's Bayesian Bootstrap.
- Audited the user-facing vignettes for stale API language and added foundational references for forensic P/S modelling, model-comparison criteria, and uncertainty methods.
- Added Curran, Buzzini and Trejos (2024), Estimating probability terms for the background presence of glass when considering activity in forensic casework, as the forensic sampling-uncertainty reference and retained Rubin (1981) for the Bayesian Bootstrap.
- Modernized the package title, description, and CITATION entry and made the cited package version follow DESCRIPTION automatically.
- Corrected the Coulson et al. bibliography year typo and normalized stored DOI fields.
- Kept the extension vignette focused on the fitPS public model contract rather than adding unnecessary external citations.
- Validated with the strict fitPS documentation, test, and package-check workflow.

## fitPS 1.1.3.011

- Made fit() the canonical public entry point for maximum likelihood, parametric Bayesian inference, the ordinary nonparametric bootstrap, and Rubin's Bayesian Bootstrap.
- Added method = "bootstrap" to return the maximum-likelihood psFit with psBootstrap uncertainty attached, and method = "bayesianBootstrap" to return psBayesianBootstrap inference.
- Removed the standalone bayesianBootstrap() public wrapper while retaining the psBayesianBootstrap class and its summary/print methods.
- Deprecated bootstrapFit() with guidance to use fit(..., method = "bootstrap") while preserving its historical compatibility behaviour.
- Updated bootstrap diagnostics, examples, README material, and uncertainty/model-comparison vignettes to use the unified fit interface and explain all four inferential routes.
- Added deterministic regression coverage for fit-based bootstrap dispatch, deprecation behaviour, Bayesian Bootstrap dispatch, and unsupported ordinary-bootstrap model coverage.
- Validated with the strict fitPS documentation, test, and package-check workflow.

## fitPS 1.1.3.007

- Added the public bayesianBootstrap() interface for Rubin Bayesian Bootstrap inference with built-in zeta, zero-inflated zeta, and logarithmic models.
- Added print() and summary() methods for psBayesianBootstrap while keeping raw parameter/probability replicates and failed-fit diagnostics available on the object.
- Documented that Rubin Bayesian Bootstrap uncertainty is over empirical observation weights and is distinct from the parametric Bayesian psPosterior architecture.
- Added deterministic public-API tests covering object identity, summaries, printing, validation, and separation from psPosterior and psBootstrap.
- Recorded future work on likelihood-matched conjugate zeta priors, including effective-sample-size and prior-geometric-mean interpretations, without expanding Stage 11 scope.
- Validated with the strict fitPS documentation, test, and package-check workflow.

## fitPS 1.1.3.006

- Added an internal Rubin Bayesian Bootstrap replicate engine using grouped Dirichlet category weights derived directly from observed survey counts.
- Scaled grouped weight draws to the observed sample size and refitted built-in zeta, zero-inflated zeta, and logarithmic models through the Stage 11 weighted-MLE boundary.
- Added the provisional psBayesianBootstrap object with parameter/probability replicates, equal-tail summaries, seed metadata, and explicit failed-fit diagnostics.
- Kept failed weighted fits as missing replicate rows rather than redrawing or smoothing them, including genuinely singleton-support surveys.
- Added deterministic tests for positive normalized weights, seeded reproducibility, grouped-Dirichlet mean behaviour, replicate dimensions, failure retention, and built-in model coverage.
- No exported Bayesian Bootstrap user-facing function is introduced yet; public API integration is deferred to the next stage.
- Validated with the strict fitPS documentation, test, and package-check workflow.

## fitPS 1.1.3.004

- Added an internal weighted-MLE fitting boundary for zeta, zero-inflated zeta, and logarithmic models using arbitrary positive real category weights.
- Kept observed psData counts unchanged by substituting weights only in a private likelihood copy, avoiding makePSData() rounding and rep()-based expansion.
- Reused the established model log-likelihood and probability contracts so weighted fits are ready for grouped Dirichlet Bayesian-Bootstrap replicates.
- Added deterministic regression coverage for ordinary-count equivalence, fractional weights, common-scale invariance, invalid weights, and sparse-support handling.
- Corrected the psData rounding regression test so warning detection and returned-object inspection are tested separately across testthat versions.
- No public Bayesian Bootstrap API or uncertainty class is introduced in this stage.
- Validated with the strict fitPS documentation, test, and package-check workflow.

## fitPS 1.1.3.001

- Audited Bayesian Bootstrap weighted fitting before introducing any public Bayesian Bootstrap API.
- Confirmed that zeta, zero-inflated zeta, and logarithmic core log likelihoods already accept arbitrary positive real weights mathematically.
- Identified public psData construction and legacy rep()-based expansion paths as integer-count boundaries that should remain separate from Bayesian-Bootstrap weights.
- Established grouped Dirichlet draws with parameters equal to observed category counts as the efficient observation-level-equivalent representation.
- Recommended a small internal weighted-MLE layer for the next stage while keeping psBootstrap, psPosterior, and the future psBayesianBootstrap statistically distinct.
- Package validation was intentionally skipped because this stage changes only development audit documentation.

## fitPS 1.1.2.008

- Repaired Stage 10.4 source-package completion after verbose `devtools::build()` output obscured the returned archive path.
- Preserved the completed unconstrained-parameter API rename and extension-vignette documentation from Stage 10.4.1.
- Recorded the exact source archive path directly from R before installing that verified build.
- Retained portable cleanup of older completed repository archives after the new archive is verified.
- Revalidated documentation, the extending vignette, tests, package checks, source build, and installation.

## fitPS 1.1.2.007

- Renamed the public MCMC transformation contract to `modelToUnconstrained()`, `modelFromUnconstrained()`, and `modelLogJacobian()` so the API states its mathematical purpose directly.
- Removed the unreleased `modelToWorking()`, `modelFromWorking()`, and `modelWorkingLogJacobian()` names rather than retaining compatibility aliases.
- Expanded the extending vignette to explain natural and unconstrained parameter scales, inverse transformations, Jacobian direction, and one- and multi-parameter examples.
- Updated built-in models, generic MCMC, external-model fixtures, and public-contract regression tests to use the new terminology.
- Replaced non-portable completed-archive cleanup with a macOS-safe Bash glob loop in the Stage 10.4.1 runner.
- Loaded the changed source namespace before direct vignette rendering so extension-method registration validates the renamed generics before package installation.
- Validated by roxygen regeneration, direct rendering of the extending vignette, strict tests, package check, source build, and installation.

## fitPS 1.1.2.005

- Simplified vignette filenames to `basics`, `model-comparison`, `uncertainty`, and `extending`, while removing the redundant `simple_P_fit` vignette.
- Kept Roux as the worked dataset inside the model-comparison vignette rather than encoding the dataset name in the filename.
- Added `.gitignore` and `.Rbuildignore` rules for transient LaTeX, knitr, and rmarkdown files created while rendering vignettes.
- Cleaned ignored vignette rendering products after direct rendering so validation leaves only deliberate source changes.
- Validated by direct vignette rendering followed by the full roxygen, strict test, package check, build, and install workflow.

## fitPS 1.1.2.003

- Allowed Bayesian fitting to use surveys with a single occupied support value instead of rejecting them during shared observation conversion.
- Preserved the historical MLE singleton-support rejection while moving that restriction into the MLE fitting layer.
- Added an explicit Bayesian warning that singleton-support posterior inference may be strongly influenced by the prior.
- Added regression coverage using the Pettard S survey, whose six observations are all S1, across zeta, logarithmic, and zero-inflated zeta Bayesian paths.
- Stage 10.2.1 repaired the warning-capture and numeric-type assertions exposed by the initial validation attempt.
- Validated with roxygen regeneration, strict tests, and package check under the full fitPS stage workflow.

## fitPS 1.1.2.001

- Added a Roux-data model-comparison vignette covering the built-in zeta, zero-inflated zeta, and logarithmic models through the public `fit()` API.
- Compared maximum-likelihood fits with log-likelihood, AIC, and BIC, and Bayesian numerical fits with DIC without treating the criteria as interchangeable.
- Added observed-versus-fitted probability diagnostics and clarified plug-in probabilities, posterior parameter uncertainty, and posterior predictive probabilities.
- Kept Stage 10.1 documentation-only; Bayesian Bootstrap implementation remains deferred to later Stage 10 work.
- Validated by directly rendering the executable vignette with the package vignette toolchain; full package-wide validation was intentionally skipped for this vignette-only stage.

## fitPS 1.1.1.021

- Completed the Stage 9 Bayesian architecture audit and confirmed that external models can use fitPS-owned numerical and MCMC posterior engines through the public model mathematics contract.
- Recorded the numerical-engine policy of deterministic integration for one- and two-parameter models and MCMC for models with three or more parameters.
- Recorded specialised built-in MCMC methods as compatibility overrides rather than extension requirements, and kept Laplace and importance as secondary existing capabilities.
- Confirmed external Poisson and Poisson-normal regression coverage for posterior summaries, DIC, P/S support mapping, serialization, and deliberate unsupported-engine failures.
- No installed package behaviour changed; package validation was intentionally skipped because this stage changes only development audit documentation.

## fitPS 1.1.1.020

- Confirmed that zeta, logarithmic, and zero-inflated zeta models all satisfy the generic MCMC model contract.
- Retained specialised built-in MCMC samplers as compatibility overrides where established stochastic behaviour remains part of the regression contract.
- Clarified that external and future models do not need model-specific sampler implementations because the generic MCMC fallback is complete.
- Added deterministic regression coverage for the built-in generic MCMC fallback and documented the Stage 9.8 consolidation decision.

## fitPS 1.1.1.019

- Updated all executable vignettes to teach the generic `fit()` model-descriptor API rather than deprecated distribution-specific fitting wrappers.
- Updated the Bayesian vignette to present deterministic numerical integration for one- and two-parameter models and MCMC as the core simulation route, while keeping Laplace and importance sampling secondary.
- Expanded the external-model vignette to document the public Bayesian prior, control, transform, numerical, and MCMC extension contract with Poisson and Poisson-normal examples.
- Documented stable standard-Normal evaluation of the Poisson-normal marginal probability and the numerical-engine policy that models with three or more parameters use MCMC.
- Validated the documentation-only change by rendering every changed vignette directly with its established R Markdown output format.

## fitPS 1.1.1.018

- Consolidated fragmented `psData` methods and operations into `R/psData-methods.R` without changing their interfaces or calculations.
- Consolidated small `psFit` fitted, log-likelihood, print, and summary methods into `R/psFit-methods.R`.
- Grouped survey-comparison helpers and bootstrap inference wrappers into cohesive source files while preserving public behaviour.
- Left the large zeta and ZIZ model files unchanged pending a separate decision about specialised versus generic MCMC implementations.
- Regenerated documentation and validated the behaviour-preserving source-layout cleanup with the strict fitPS test and check workflow.

## fitPS 1.1.1.017

- Extended generic numerical posterior fitting to two-parameter models with adaptive `cubature::hcubature()` integration over model-supplied natural-scale bounds.
- Adopted a dimensional policy of deterministic numerical fitting for one and two parameters and MCMC for models with three or more parameters.
- Migrated ZIZ posterior normalization, moments, and DIC to the generic two-dimensional engine while retaining a compatibility grid for established plotting helpers.
- Stabilized the external Poisson-normal proof likelihood by integrating over a standard-Normal latent variable rather than a collapsing narrow Normal density.
- Added targeted numerical regression validation through `devtools::test()` so the package and test helpers are loaded before the full suite runs.
- Added `cubature` as a package import for core two-dimensional numerical integration and kept Laplace and importance as secondary existing capabilities.
- Validated through targeted numerical tests, roxygen regeneration, strict full tests, package check, source build, and installation.

## fitPS 1.1.1.011

- Made one-dimensional numerical posterior fitting generic through the public model likelihood, prior, and Bayesian-control contract.
- Migrated zeta and logarithmic numerical fitting and probability summaries onto the shared numerical engine while preserving their model mathematics.
- Added external Poisson numerical fitting with infinite natural-scale support and conjugate posterior regression checks.
- Updated the external-model regression contract to recognise both numerical and MCMC support for the Poisson proof model.
- Added the explicit `stats::setNames` roxygen import required by the generic numerical engine.
- Kept ZIZ numerical fitting model-specific because its established posterior representation is two-dimensional.
- Validated through roxygen regeneration, strict tests, package check, source build, and installation.

## fitPS 1.1.1.008

- Completed the multi-parameter external Bayesian proof using the Poisson-normal test model through the generic MCMC contract.
- Verified posterior parameter summaries, generic derived probability summaries, fitted probabilities, and DIC for an external two-parameter model.
- Verified matched P- and S-survey data retain identical zero-based support behaviour through Bayesian fitting.
- Verified serialization preserves the external model descriptor, MCMC representation, posterior summaries, and fitted values.
- Kept Laplace and importance as secondary existing capabilities without expanding their public extension role.
- Validated through roxygen regeneration, strict tests, package check, source build, and installation.

## fitPS 1.1.1.007

- Added a model-neutral random-walk Metropolis fallback so external psModel classes can use fitPS-owned MCMC posterior orchestration.
- Let external Bayesian models supply model-specific prior objects through modelLogPrior() rather than requiring the legacy one-dimensional psPrior representation.
- Proved the generic Bayesian path with deterministic external Poisson and two-parameter Poisson-normal models, including transformed positive parameters.
- Updated the public-model regression expectation so missing external priors now fail with the new explicit model-specific-prior contract rather than the retired MLE-only error.
- Preserved the specialised zeta, ZIZ, and logarithmic MCMC implementations and kept Laplace and importance as secondary existing capabilities.
- Validated through roxygen regeneration, strict tests, package check, source build, and installation.

## fitPS 1.1.1.005

- Added a public Bayesian model contract for complete prior evaluation, model-owned Bayesian controls, parameter transformations, and working-scale Jacobians.
- Added identity transform defaults so unconstrained external models can use the contract without transformation boilerplate.
- Exposed existing zeta, logarithmic, and zero-inflated zeta prior and transformation mathematics through the new public S3 methods without changing their current Bayesian fitting paths.
- Added deterministic contract tests for external methods, built-in prior evaluation, transform round trips, Jacobians, and Bayesian starting values.
- Documented method-specific Bayesian prior and starting-value arguments so generated Rd usage and argument sections remain synchronized.
- Generic external-model Bayesian fitting remains intentionally deferred until posterior engines are migrated to consume the new contract.

## fitPS 1.1.1.002
- Repaired the Stage 9.1 version formatter, which had recycled the build component and produced an invalid six-part version; build 001 remains consumed by that attempted stage.

- Audited the Bayesian architecture for external model extensibility and defined the minimum public model contract needed for fitPS-owned posterior engines.
- Identified model-specific posterior fitting dispatch, scalar prior assumptions, and built-in parameter transformations as the main barriers to third-party Bayesian models.
- Recorded a staged implementation path beginning with public prior, Bayesian-control, and parameter-transformation methods, followed by generic MCMC proof models.
- No installed package behaviour changed; package validation was intentionally skipped because this stage changes only development audit documentation.

## fitPS 1.1.0.018

- Completed the Stage 8 public model API audit and repaired the multi-parameter external-summary regression assertion.
- Verified external Poisson and integral-based Poisson-normal models use only the public extension contract, including P/S support shifting without truncation.
- Preserved deprecated compatibility fitters while keeping zeta, ZIZ, and logarithmic implementation consolidated in model-specific source files.
- Documented that generic third-party fitting is MLE-only at this boundary; unsupported Bayesian external fitting fails explicitly.
- Preserved generic print/summary presentation for external MLE models and validated with roxygen regeneration, strict tests, package check, source build, and installation.

## fitPS 1.1.0.015

- Deprecated fitDist(), fitZIDist(), and fitlogDist() as compatibility entry points and direct users to the generic fit() interface.
- Kept legacy fitters functional through common model-oriented dispatch and explicitly regression-tested their retained model metadata.
- Consolidated distribution-specific zeta, zero/one-inflated zeta, and logarithmic implementations into model-zeta.R, model-ziz.R, and model-logarithmic.R.
- Renamed shared numerical, MCMC, Laplace, and importance posterior code as model-neutral posterior-engine source files.
- Corrected the external Poisson and Poisson-normal examples so S surveys shift the zero-based probability sequence from P0, P1, ... to S1, S2, ... without truncation or renormalisation.
- Expanded extension documentation and tests to make P/S support mapping an explicit part of the public model contract.
- Validated with roxygen regeneration, strict tests, a strict package check, source build, and installation.

## fitPS 1.1.0.010

- Corrected the external Poisson-normal example to use a latent normal log-rate integrated out of the Poisson likelihood rather than the Hermite distribution.
- Added the Poisson probability mass function, weighted likelihood, Poisson-normal marginal integral, marginal likelihood, and moment relationships to the public extension vignette.
- Updated the external Poisson-normal tests to evaluate marginal probabilities by numerical integration while keeping all distribution-specific code outside R/.
- Preserved the public mu and sigma parameterization for the latent normal distribution and the generic fitPS MLE, prediction, serialization, and model-comparison interfaces.
- Validated with roxygen regeneration, strict tests, a strict package check, source build, and vignette build during package construction.

## fitPS 1.1.0.009

- Corrected the external Poisson-normal example to use a latent normal log-rate integrated out of the Poisson likelihood rather than the Hermite distribution.
- Added the Poisson probability mass function, weighted likelihood, Poisson-normal marginal integral, marginal likelihood, and moment relationships to the public extension vignette.
- Updated the external Poisson-normal tests to evaluate marginal probabilities by numerical integration while keeping all distribution-specific code outside R/.
- Preserved the public mu and sigma parameterization for the latent normal distribution and the generic fitPS MLE, prediction, serialization, and model-comparison interfaces.
- Validated with roxygen regeneration, strict tests, a strict package check, source build, and vignette build during package construction.

## fitPS 1.1.0.008

- Replaced the earlier maintainer-facing extension vignette with a public guide for defining models outside fitPS.
- Added worked external Poisson and Poisson-normal examples using only the exported psModel contract and generic fit() interface.
- Documented downstream S3 registration, model-specific MLE controls, retained model objects, prediction, serialization, and likelihood-based model comparison.
- Kept all Poisson and Poisson-normal distribution mathematics outside fitPS package implementation code.
- Validated with roxygen regeneration, strict tests, a strict package check, source build, and vignette build during package construction.

## fitPS 1.1.0.007

- Demonstrated the public model-extension contract with a second distribution defined entirely outside fitPS package code.
- Added an external two-parameter Poisson-normal (Hermite) test model parameterized by `mu` and `sigma`, with model-specific mathematics and S3 registration confined to test helper code.
- Confirmed the unchanged generic MLE path handles coupled parameter constraints, fitted probabilities, prediction, log likelihood, deviance, AIC, BIC, and serialization for a multi-parameter external model.
- Added no distribution-specific implementation under `R/`; Stage 8.4 therefore tests the Stage 8.3 public API without extending fitPS internals.
- Validated with roxygen regeneration, strict tests, and a strict package check.

## fitPS 1.1.0.006

- Added a supported public `psModel()` constructor and exported the minimal S3 generics needed for third-party model extensions.
- Added fitPS-owned generic maximum-likelihood optimisation using model-declared starting values and bounds, with optional `modelMleControl()` overrides for more demanding models.
- Demonstrated true external extensibility with a Poisson model defined entirely in test helper code using public APIs only.
- External Poisson fits now participate in fitted probabilities, point prediction, log likelihood, deviance, AIC, BIC, and fitted-object serialization without modifying fitPS source code.
- Kept Bayesian support for arbitrary external models deliberately out of the initial contract until posterior-engine compatibility is defined explicitly.
- Validated with roxygen regeneration, strict tests, and a strict package check.

## fitPS 1.1.0.004

- Added the model-oriented `fit()` entry point for built-in zeta, zero-inflated zeta, and logarithmic model objects.
- Fixed generic `fit()` method dispatch so its default method vector is resolved to the scalar `mle` choice before delegation to model-specific fitters.
- Corrected the built-in model interface regression test so list names retained by `vapply()` do not cause a false failure when the fitted model identifiers are otherwise identical.
- New `psFit` objects retain their originating model descriptor while preserving the established character model identifier.
- Model-aware operations now prefer the retained descriptor and fall back to built-in name reconstruction for legacy fitted objects.
- Exported the built-in model constructors needed by the new fitting interface; third-party model construction remains deferred until the external Poisson contract is validated.
- Validated with roxygen regeneration, strict tests, and a strict package check.

## fitPS 1.1.0.001

- Audited the internal model protocol and defined the Stage 8 public-extension direction.
- Established external Poisson and Poisson-normal models as proof-of-concept targets that must work without rebuilding or modifying fitPS.
- Chose retained model descriptors with legacy name fallback as the preferred route away from hard-coded model reconstruction.
- Started the Stage 8 development series at 1.1.0.xxx; package validation was intentionally skipped because this stage changes only a development audit document.

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
