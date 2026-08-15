# Stage 6.11 final regression and architecture audit

## Purpose

Stage 6.11 closes the Bayesian architecture migration by auditing the committed
Stage 6.10.1 source, adding durable architecture-contract tests, and removing
remaining migration-era wording from internal documentation. It deliberately
does not add another distribution or model-comparison criterion.

## Architecture status

The Bayesian implementation now separates model and posterior-engine concerns.
The model layer defines parameter names, observation handling, probabilities,
log likelihoods, and supported engines. The engine layer defines posterior
fitting, summaries, diagnostics, point estimates, and representation-specific
probability operations. Bayesian fitting converges through `fitBayesianModel()`
and `finaliseBayesianPsFit()`, and fitted objects expose the common
`psPosterior` contract.

The supported model/engine matrix is intentionally sparse:

| Model | numerical | mcmc | laplace | importance |
| --- | --- | --- | --- | --- |
| zeta | yes | yes | no | no |
| ZIZ | yes | yes | yes | yes |

Unsupported combinations fail through the shared model/engine validation rather
than acquiring placeholder implementations merely for symmetry.

## Regression coverage

The existing behaviour-based test suite characterises each supported migration:
plain-zeta numerical and MCMC fits; ZIZ numerical, MCMC, Laplace, and importance
fits; posterior probability summaries; inflation probabilities; posterior
presentation; credible intervals; and public Bayesian orchestration.

Stage 6.11 adds architecture-contract tests that verify the supported-engine
matrix, confirm the unpublished provisional Bayesian fitter functions are absent
from the package namespace, and demonstrate that a new model descriptor can
declare supported engines without changing the engine abstraction.

All test filenames and `test_that()` descriptions use stable behavioural names,
not development-stage numbers.

## Documentation audit

A static scan of `R/` found 218 named functions and no function lacking an
immediately preceding roxygen block. Remaining internal roxygen text was also
reviewed for migration-era language. Two references to the development stage
were replaced with architecture-neutral descriptions.

## Duplication and dispatch audit

No Bayesian fitting or posterior operation retains a switch on posterior method
names. Engine construction maps a validated method name to an engine object once;
after that point operations use S3 dispatch. The older `method = "integrate"`
translation in Bayesian options is retained because it supports established
plain-zeta behaviour predating the unpublished migration work.

The method-specific Bayesian fitter wrappers removed in Stage 6.10 remain absent.
The common finalisation path is the single source for `psPosterior` construction,
posterior summaries, point-estimate fitted probabilities, and core Bayesian
`psFit` assembly. Established pre-1.0.7 plain-zeta compatibility fields remain
isolated behind the model-specific compatibility method.

## Extensibility assessment

The architecture is sufficiently small to test with another distribution without
changing the engine abstraction. A future distribution should be able to define
a `psModel` subclass, its natural parameters, likelihood and probability methods,
and the engines it supports. Engine-specific mathematics is needed only where the
chosen approximation genuinely requires model-specific work.

The logarithmic distribution is therefore a suitable post-migration proof of
extensibility and a strong candidate for a developer-facing vignette showing how
to add a distribution. That work is deferred until after Stage 6 is closed.

## Model comparison follow-on

AIC and DIC should also be added after the migration. AIC should use a common
likelihood/parameter-count contract. DIC should use the model deviance together
with posterior summaries supplied through the posterior representation/engine
contract. Keeping these criteria outside Stage 6 avoids mixing new statistical
features with the architecture migration.

## Stage 6 conclusion

Subject to the authoritative Stage 6.11 roxygen, strict test, and package-check
workflow, the migration can be considered complete. The package has a shared
model/engine architecture, a consistent posterior contract, no provisional
method-specific Bayesian fitter paths, comprehensive behavioural regression
coverage, and a clear route for testing extensibility with a new distribution.
