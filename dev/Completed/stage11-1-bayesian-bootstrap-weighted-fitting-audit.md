# Stage 11.1 Bayesian Bootstrap weighted-fitting architecture audit

## Purpose

Stage 11.1 audits the fitPS fitting architecture before any public Rubin Bayesian Bootstrap API is implemented. The immediate question is whether built-in models can be fitted from arbitrary positive real observation weights, especially when an aggregated `psData` object represents sums of individual Bayesian-Bootstrap weights within occupied support categories.

This stage changes no installed package behaviour.

## Statistical representation

Rubin's Bayesian Bootstrap assigns observation-level weights

`(w1, ..., wn) ~ Dirichlet(1, ..., 1)`.

If occupied support category `j` contains `n_j` identical observations, Dirichlet aggregation gives

`(W1, ..., Wk) ~ Dirichlet(n_1, ..., n_k)`.

Therefore fitPS does not need to expand a survey into individual observations merely to generate Bayesian-Bootstrap weights. A draw over occupied categories with Dirichlet parameters equal to the observed integer counts is mathematically equivalent to drawing observation-level Dirichlet weights and summing them within categories.

The scale of the weights is immaterial for the MLE because multiplying every likelihood weight by one common positive constant multiplies the log likelihood by that constant without changing its maximiser. For compatibility with existing count-based calculations, a convenient internal representation is likely to scale each grouped Dirichlet draw to total sample size `N`, so the replicate weights sum to `N` rather than one. This preserves the ordinary fit when weights equal the original counts and avoids introducing fractional effective sample-size conventions accidentally.

## Current ordinary bootstrap

`bootstrapParameterReplicates()` in `R/psBootstrap.R` expands the aggregated survey using

`rep(x$data$n, x$data$rn)`

and resamples the resulting individual observations with replacement. Each replicate is then re-aggregated into integer frequencies before fitting.

This is the correct observation-level nonparametric bootstrap. It must remain distinct from Rubin's Bayesian Bootstrap.

Failed ordinary-bootstrap MLEs are already retained as missing replicate rows and summarised in diagnostics through `nSuccessful`, `nFailed`, and `failureRate`. Stage 11 should not silently redraw failed replicates or apply Laplace/add-one smoothing to make them fit. Such smoothing would change the estimator and, for infinite-support models, would require an additional support/prior assumption.

## Built-in likelihood audit

### Zeta

Both the established MLE path in `fitDistImpl()` and `modelLogLikelihood.zetaModel()` evaluate a weighted sum of per-support log probabilities:

`sum(rn * logProbability)`.

Nothing in that likelihood expression requires `rn` to be integer. The zeta MLE variance helper uses `N = sum(rn)`, which also remains meaningful when Bayesian-Bootstrap weights are scaled to sum to the observed sample size.

Conclusion: the core zeta likelihood is already compatible with arbitrary positive real weights.

### Zero-inflated zeta

The established MLE objective and `modelLogLikelihood.zizModel()` / `zizLogLikelihood()` also use

`sum(rn * logProbability)`.

The generic likelihood mathematics therefore accepts arbitrary positive real weights. Several legacy helpers still create expanded vectors using `rep(obsData, x$data$rn)`. The principal current MLE objective does not need that expanded vector, but these expansions are integer-count assumptions and should not be part of the Bayesian-Bootstrap fitting path.

The specialised Bayesian engines also read `x$data$rn` as counts and use weighted likelihood expressions. Stage 11 Bayesian Bootstrap does not require running a parametric Bayesian engine inside each replicate, so those paths need not be changed merely to implement Bayesian Bootstrap MLE refits.

Conclusion: the core ZIZ likelihood is compatible with arbitrary positive real weights, but the weighted fitting path should bypass legacy count expansion.

### Logarithmic

`fitlogDistImpl()` delegates its MLE objective to `modelLogLikelihood.logarithmicModel()`, which uses

`sum(rn * logProbability)`.

Nothing in this expression requires integer frequencies.

Conclusion: the logarithmic core likelihood is already compatible with arbitrary positive real weights.

### External model contract

The generic `fitModel.psModel()` MLE path evaluates `modelLogLikelihood(model, parameters, data = x)`. External model methods commonly receive `x$data$rn` through the existing `psData` structure, but the public contract currently describes survey data rather than arbitrary real fitting weights.

Stage 11 should not silently redefine the public `psData` contract merely to support Bayesian Bootstrap. Bayesian-Bootstrap weighted fitting should first be an internal facility for built-in models. A later stage can decide whether arbitrary weighted external-model fitting deserves a separate public contract.

## Main architectural blocker: psData construction

`makePSData()` currently treats non-integer counts as invalid survey counts. When either support labels or `rn` values are non-integer it rounds both and warns:

`Non-integer input detected. Both input vectors have been rounded.`

That behaviour is appropriate for the historical meaning of `psData`: observed survey counts are integer frequencies. It is incompatible with Bayesian-Bootstrap category weights.

Changing `makePSData()` globally to accept fractional counts would blur the distinction between observed frequencies and inferential weights and could alter existing user expectations.

Recommendation: preserve ordinary `psData` construction semantics and introduce a narrow internal weighted-data representation or constructor used only during weighted fitting. The representation may retain the same support/type information while keeping a separate positive numeric weight vector rather than pretending the weights are observed integer frequencies.

## Other integer-count assumptions

The following areas expand or otherwise assume integer frequencies and should not be used directly with fractional Bayesian-Bootstrap weights:

- `R/psBootstrap.R`: ordinary bootstrap expansion with `rep()`; correct for the existing frequentist bootstrap and should remain unchanged.
- `R/psData-methods.R`: conversion to expanded data frames uses `rep()`; this describes observed survey units and should remain count-based.
- `R/combineSurveys.R`: expands surveys with `rep()`; not required for Bayesian Bootstrap replicate fitting.
- legacy/profile ZIZ helpers contain `rep(obsData, x$data$rn)`; weighted fitting should use the model likelihood rather than these expansion paths.
- some confidence-region code expands observations using `rep()`; Bayesian Bootstrap need not use that code for replicate fitting.

Functions that compute weighted means or empirical proportions using multiplication and `sum(rn)` are mathematically capable of accepting real weights, but they should not automatically be interpreted as Bayesian-Bootstrap summaries without an explicit interface.

## Starting values and optimisation

Current built-in starting values are fixed defaults or user-supplied values rather than functions requiring integer counts:

- zeta: default shape 2;
- ZIZ: default `(pi, shape) = (0.5, 2)`;
- logarithmic: default `pi = 0.5`.

The corresponding `optim()` objectives are weighted log likelihoods. No built-in starting-value calculation found in the principal MLE paths requires integer frequencies.

Because a Bayesian-Bootstrap draw assigns strictly positive weight to every originally occupied category with probability one, it cannot lose occupied support categories. Thus a non-degenerate original survey will not become singleton-support merely because of a Bayesian-Bootstrap draw. This removes the ordinary-bootstrap mechanism that can create all-zero or all-one resamples.

If the original survey itself has only one occupied support category, Bayesian Bootstrap does not repair that limitation. The existing MLE singleton-support rule should continue to apply to weighted MLE refits. Parametric Bayesian fitting under a proper prior remains a distinct approach that may still be meaningful for such surveys.

## Probability calculations

`modelProbabilities()` evaluates model-implied P/S probabilities from fitted model parameters and does not depend on observation frequencies. It is therefore directly reusable after each weighted MLE fit.

Stage 11 should keep two possible Bayesian-Bootstrap quantities conceptually distinct:

1. direct probabilities of the random empirical distribution represented by a Dirichlet weight draw;
2. model-implied probabilities after refitting a parametric model to that weighted empirical distribution.

The planned model-based Bayesian Bootstrap should compute the second quantity. A direct empirical Bayesian-Bootstrap summary, if ever added, should be a separate explicitly named operation.

## What can be reused from psBootstrap

Useful design and implementation ideas that can be shared without sharing the statistical class include:

- validation of `B`, `level`, and seed inputs;
- seeded reproducibility conventions;
- replicate orchestration and optional parallel execution patterns;
- parameter and probability summary shapes;
- retention of failed fits as missing replicate rows;
- diagnostics with requested, successful, and failed replicate counts;
- probability transformation through model probability functions;
- plotting conventions where they remain semantically appropriate.

The `psBootstrap` class itself should not be reused because it denotes frequentist observation resampling. Bayesian Bootstrap should remain a separate class, provisionally `psBayesianBootstrap`.

## Proposed minimal internal weighted-fitting interface

Stage 11.2 should add a deliberately small internal abstraction rather than changing the public data API. A suitable direction is:

- preserve the original `psData` object as the source of support labels, survey type, notes, and observed integer counts;
- carry a separate positive finite numeric vector of fitting weights aligned with the occupied support rows;
- provide an internal model log-likelihood evaluator that substitutes these weights for `x$data$rn` without mutating or reclassifying the original survey;
- provide an internal weighted-MLE helper that uses the model descriptor, `modelMleControl()`, `modelLogLikelihood()`-compatible mathematics, `optim()`, and `modelProbabilities()`;
- add built-in regression tests proving equal/original weights reproduce ordinary MLE parameter estimates and fitted probabilities within numerical tolerance;
- add tests proving arbitrary positive non-integer weights work for zeta, ZIZ, and logarithmic models.

A possible implementation is an internal lightweight copy of the survey with a separate weight field, but the preferred design is to avoid making arbitrary external model methods accidentally interpret non-integer `rn` as observed counts. The exact helper signature should be decided during Stage 11.2 implementation after checking the least invasive route through model likelihood dispatch.

## Bayesian-Bootstrap draw representation

For an observed aggregated survey with counts `n_j`, draw grouped weights directly as

`G_j ~ Gamma(n_j, 1)` independently,

then normalise

`W_j = G_j / sum(G)`.

This produces `Dirichlet(n_1, ..., n_k)` without adding a new package dependency. For MLE fitting, scale to

`N * W_j`, where `N = sum(n_j)`.

This gives positive real weights summing to the original sample size. The draw remains mathematically equivalent to observation-level `Dirichlet(1, ..., 1)` aggregation.

The implementation should not reconstruct thousands of individual observations for surveys such as Roux.

## Proposed object and API direction

Do not finalise the public API in Stage 11.1. The current best direction remains a distinct object:

`psBayesianBootstrap`

with likely components for:

- raw parameter replicates;
- model-implied probability replicates or summaries;
- interval level;
- replicate diagnostics;
- seed/reproducibility information;
- method/weighting metadata.

There is no need to store every full Dirichlet draw by default if parameter/probability replicates and reproducibility metadata are retained. Optional diagnostic retention can be considered later.

The public operation should probably attach this object to an existing MLE `psFit`, paralleling the ergonomics of `bootstrapFit()` while using statistically explicit naming. Final names should wait until the internal weighted-MLE path and object are implemented.

## Failure handling

Bayesian-Bootstrap weighted fits should follow the useful existing bootstrap principle that every requested replicate remains accounted for. A failed optimisation should produce a retained failed replicate and diagnostic count, not a replacement redraw.

No Laplace/add-one smoothing should be introduced to rescue failed weighted MLEs. For infinite-support distributions it has no natural finite-category interpretation, and any pseudo-count scheme would introduce an additional modelling assumption.

Because grouped Dirichlet weights remain strictly positive on all observed categories, Stage 11 should expect materially fewer support-collapse failures than the ordinary bootstrap for non-degenerate surveys. Remaining failures would indicate optimisation/numerical issues or an original-data/model identifiability problem rather than bootstrap resampling deleting support.

## Stage 11.1 decision

The architecture is suitable for Bayesian Bootstrap with a small internal weighted-fitting layer.

The decisive findings are:

1. zeta, ZIZ, and logarithmic core log likelihoods already support arbitrary positive real weights mathematically;
2. `makePSData()` deliberately rounds fractional counts, so Bayesian-Bootstrap weights must not be routed through ordinary public `psData` construction;
3. several legacy helpers expand integer frequencies with `rep()`, but the principal likelihood paths do not require expansion;
4. grouped `Dirichlet(n_j)` draws are exactly equivalent to aggregating observation-level `Dirichlet(1, ..., 1)` weights and should be used directly;
5. model-implied probability calculations are already reusable after each weighted fit;
6. `psBootstrap` orchestration and diagnostic ideas can be reused, but the class and interpretation must remain separate;
7. no public Bayesian Bootstrap API should be added until the internal weighted-MLE path has deterministic regression coverage.

## Recommended Stage 11.2 objective

Implement and test the minimal internal weighted-MLE infrastructure for built-in models.

Stage 11.2 should not yet implement the full Bayesian-Bootstrap public feature. It should establish deterministic weighted fitting for zeta, ZIZ, and logarithmic models, verify equality with ordinary fitting under original/equal weights, verify arbitrary positive fractional weights, and avoid any `rep()` expansion or change to public `makePSData()` count semantics.
