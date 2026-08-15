# Stage 8.7 public API and regression audit

## Purpose

This audit closes Stage 8 by checking the public model-extension architecture after the generic `fit()` interface, external Poisson and Poisson-normal proofs, compatibility-wrapper deprecations, P/S support mapping, and model-source consolidation.

## Public fitting contract

The public fitting entry point is `fit(x, model = ...)`. Built-in constructors `zetaModel()`, `zizModel()`, and `logarithmicModel()` are exported, and external code can construct a `psModel` subclass with `psModel()`.

The supported external MLE contract is deliberately small. External models can provide methods for `modelObservationData()`, `modelProbabilities()`, `modelLogLikelihood()`, and, when needed, `modelMleControl()`. Parameter names and declared posterior-engine support are available through `modelParameterNames()` and `supportedPosteriorEngines()`.

Object-based fitting does not require a model registry or a built-in name switch. New fitted objects retain their originating descriptor in `modelObject`; `modelFromFit()` uses the retained descriptor first and falls back to built-in character-name recovery only for legacy fits.

## External-model proof

Two external test models demonstrate the contract without adding their distribution mathematics under `R/` and without using `fitPS:::`:

- Poisson is a one-parameter closed-form MLE example.
- Poisson-normal is a two-parameter MLE example whose probabilities require numerical integration over a latent normal log-rate.

Both examples support P and S surveys by shifting the same natural probability sequence between P0/P1/... and S1/S2/... labels. No truncation or renormalisation is introduced.

The examples participate in generic fitting, fitted probabilities, prediction, `probfun()`, `logLik()`, deviance, AIC, BIC, printing, summaries, and serialization. The Stage 8.7 presentation fallback ensures an unknown external model no longer needs a fitPS-specific branch merely to print or summarize its natural-scale MLE parameters.

## Built-in compatibility

The built-in zeta, zero/one-inflated zeta, and logarithmic models retain their established implementations and model-specific presentation where needed. Distribution-specific implementation is consolidated in:

- `R/model-zeta.R`
- `R/model-ziz.R`
- `R/model-logarithmic.R`

Shared posterior algorithms remain in model-neutral posterior-engine files.

The older `fitDist()`, `fitZIDist()`, and `fitlogDist()` entry points remain functional compatibility wrappers but emit deprecation warnings directing users to `fit()`.

## Model comparison and posterior support

MLE fits use the common model contract for `logLik()`, deviance, AIC, and BIC. Built-in Bayesian models retain the sparse engine matrix established in Stages 6 and 7, and DIC remains available where the posterior representation supports it.

External model fitting is intentionally MLE-only at the end of Stage 8. Passing a generic external `psModel` to `fit(..., method = "bayes")` fails explicitly. Stage 8 therefore establishes third-party extensibility along the MLE/model axis; it does not claim that arbitrary external models can yet opt into fitPS Bayesian posterior engines through the public `fit()` fallback.

That limitation is preferable to implying unsupported Bayesian compatibility and can be revisited after the posterior-engine architecture is exercised further.

## Stage 8 conclusion

Stage 8 achieves the principal user-facing goal: a user or downstream package can define a new probability model outside fitPS, fit it through the common public interface, and use the resulting ordinary `psFit` with the main likelihood-based fitted-object and model-comparison machinery without rebuilding or forking fitPS.

The external Poisson and Poisson-normal examples provide regression proofs of both a simple closed-form model and a more demanding numerically integrated two-parameter model. The remaining known boundary is public Bayesian-engine participation for third-party models, which is explicitly outside the MLE contract stabilized here.
