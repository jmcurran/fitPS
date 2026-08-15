# Stage 6.9 - common posterior orchestration and derived quantities

## Purpose

Stage 6.9 consolidates the Bayesian fit finalisation that remained duplicated after the individual posterior engines migrated in Stages 6.4 through 6.8.

The stage is intentionally architectural. It does not change the likelihoods, samplers, numerical grids, Laplace approximation, or importance-sampling mathematics.

## Common fit finalisation

All supported Bayesian model/engine combinations now finish through `finaliseBayesianPsFit()`.

That function is responsible for:

- obtaining the common parameter summary from the engine;
- obtaining the engine point estimate used for plug-in fitted probabilities;
- constructing the named fitted P/S probability vector;
- obtaining posterior mean probability summaries;
- constructing the common `psPosterior` object;
- storing the typed `psPosteriorRepresentation` consistently;
- setting common `psFit` fields such as `model`, `method`, and `posteriorMethod`.

Model- and engine-specific fitters continue to retain established top-level fields such as `chain`, `pdf`, `posteriorGrid`, `laplace`, `importance`, and `weightedSamples`. Those compatibility fields are deliberately left in place until Stage 6.10.

## Plain zeta posterior contract

Plain-zeta numerical and MCMC fits now attach the same `psPosterior` class used by ZIZ fits.

MCMC posterior probability summaries are calculated directly from retained shape draws. Numerical zeta probability summaries are obtained deterministically from the stored one-dimensional posterior density over the prior support.

This means `summary()`, `posteriorProbs()`, `fitted(..., type = "posteriorMean")`, and `plot()` on the attached posterior can use the same interfaces for zeta and ZIZ fits.

## Typed posterior representations

`psPosterior$representation` now stores the typed representation wrapper for every migrated Bayesian fit. The raw engine payload remains available inside that wrapper and, where required for compatibility, in the historical top-level fit field.

This establishes one invariant:

```text
fit$posterior$representation is fit$posteriorRepresentation
```

for all supported Bayesian fits.

## Derived probability dispatch

The new `summarisePosteriorProbabilities()` protocol dispatches on the posterior engine. Numerical fitting uses a secondary model dispatch because zeta uses one-dimensional integration while ZIZ uses a two-dimensional grid.

MCMC and importance summaries transform posterior parameter samples through `modelProbabilities()`. Laplace summaries preserve the existing seeded Gaussian working-scale draw calculation.

## Inflation diagnostic dispatch

`posteriorInflation()` no longer switches on character method names. It delegates to `posteriorInflationProbability()`, whose methods dispatch on:

- `numericalPosteriorRepresentation`;
- `mcmcPosteriorRepresentation`;
- `laplacePosteriorRepresentation`;
- `importancePosteriorRepresentation`.

The four established calculation helpers remain temporarily as implementation details so Stage 6.9 does not mix orchestration consolidation with deletion. Stage 6.10 can remove superseded compatibility code after the common pathway has been validated.

## Compatibility

The plug-in fitted probability vector remains evaluated at the same engine point estimate as before. ZIZ posterior probability summaries retain their existing algorithms and random-number behaviour. Existing top-level compatibility fields remain present.

The deliberate user-visible addition is that plain-zeta Bayesian fits now expose the common `psPosterior` interface and posterior mean P/S probability summaries.

## Validation expectations

The full package workflow must verify:

- existing model/engine characterisation tests still pass;
- ZIZ posterior probability summaries are unchanged;
- plug-in fitted values remain named vectors;
- zeta numerical and MCMC fits attach valid `psPosterior` objects;
- inflation probabilities remain unchanged while using S3 representation dispatch;
- roxygen regenerates all internal S3 registrations cleanly;
- no stage-numbered test filenames or test descriptions are introduced.
