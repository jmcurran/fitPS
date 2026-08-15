# Stage 7.2 logarithmic-distribution extensibility

Stage 7.2 adds the logarithmic-series distribution as a `psModel` and uses it
as a practical test of the Stage 6 model/engine abstraction.

## Extension surface

The logarithmic implementation adds a `logarithmicModel` descriptor with one
natural parameter, `pi`, and declares only the posterior engines that have a
clear one-dimensional implementation: numerical integration and MCMC.
Laplace and importance sampling remain unsupported rather than being added for
symmetry.

The model owns its probability and log-likelihood mathematics. P/S observation
mapping and probability naming reuse the shared fitPS helpers. Bayesian fitting
reuses the existing `fitPosterior()` -> `psPosterior` -> finalisation path.
The engine implementations themselves do not contain logarithmic-model branches;
logarithmic numerical and MCMC behaviour is supplied through the existing
secondary model-dispatch extension points.

## Architectural findings

Adding a new one-parameter model exposed two remaining zeta assumptions in
otherwise generic code: scalar MCMC chains were labelled `shape`, and numerical
DIC explicitly recognised `zetaModel`. Stage 7.2 removes both assumptions by
using `modelParameterNames()` for one-parameter models.

The current numerical and MCMC model methods still duplicate some mechanics
between zeta and logarithmic models. That duplication is intentionally visible
rather than hidden by model-name conditionals. Stage 7.4 should decide whether a
small reusable one-parameter posterior helper is warranted after the logarithmic
extension has been exercised in practice.

The prior constructor also previously enforced the zeta `shape > 1` domain for
all priors. Stage 7.2 separates general finite ordered prior support from the
zeta-specific domain check, while retaining the original restriction for the
loguniform zeta prior.

## Legacy logarithmic implementation

The former `fitlogDist()` implementation is not treated as a compatibility
contract. It had not been part of the published fitPS methodology. Stage 7.2
replaces its internals with a wrapper around the common model architecture while
retaining the function name and aliases as convenient user entry points.

## Model comparison

Maximum-likelihood logarithmic fits participate in the Stage 7.1 `logLik()`,
AIC, BIC, and deviance contract. Numerical and MCMC Bayesian logarithmic fits
participate in DIC through the same model-level deviance calculation used by the
existing models.
