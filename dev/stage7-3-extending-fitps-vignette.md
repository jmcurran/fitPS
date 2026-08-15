# Stage 7.3 extending fitPS vignette

Stage 7.3 documents the model-extension contract established by the Stage 6
architecture and exercised by the Stage 7.2 logarithmic-series implementation.

## Scope

The new developer-facing vignette, `Extending fitPS with a New Distribution`,
uses the logarithmic model as a real example of the extension process. It covers:

- defining a `psModel` descriptor;
- declaring and validating natural parameters;
- implementing P/S probabilities through shared support helpers;
- implementing the model log likelihood;
- declaring supported posterior engines;
- supplying model-specific numerical/MCMC methods at the existing secondary
  dispatch boundary;
- fitting maximum-likelihood and Bayesian models;
- reusing posterior summaries, fitted probabilities, and diagnostics;
- obtaining deviance, AIC, BIC, and DIC through the shared comparison contract;
- writing stable behavioural tests for a new distribution.

## Architectural position

The vignette explicitly describes the current extension surface as an internal
package-developer contract, not a public third-party plug-in API. This is an
important limitation: a distribution can now be added coherently without
copying unrelated engine machinery, but its S3 model methods still live inside
the fitPS source tree. In practical terms, additions currently require changes
by the maintainer or by someone working from a fork or modified copy of fitPS.
A downstream package cannot yet register its own model with an installed fitPS
package.

The logarithmic example also records the residual duplication between the zeta
and logarithmic one-parameter numerical/MCMC model methods. Stage 7.3 does not
refactor that duplication. Stage 7.4 should assess whether there is enough
evidence for a small reusable one-parameter posterior helper.


## Stage 7.4 audit question

Stage 7.4 should explicitly assess whether the internal model contract should be
made into a supported third-party extension API. The audit should determine:

- whether an external package could define a `psModel` subclass and methods
  without modifying fitPS;
- which constructors, generics, validators, or registration hooks would need to
  become public and stable;
- whether a future `fit(data, model = ...)` entry point should accept externally
  constructed model objects; and
- whether the maintenance burden of supporting that public contract is justified.

The audit should distinguish a technically possible S3 extension from a
deliberately supported and documented plug-in API.

## Future fitting API

The vignette notes the longer-term possibility of a model-oriented entry point
such as `fit(data, model = ...)`, while retaining `fitDist()` and `fitZIDist()`
as compatibility/convenience functions. Any future generic fitter should be a
thin front end to the existing model/engine contracts rather than a parallel
implementation.
