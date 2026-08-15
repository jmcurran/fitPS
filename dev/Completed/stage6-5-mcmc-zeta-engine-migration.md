# Stage 6.5 - MCMC zeta engine migration

## Purpose

Stage 6.5 migrates the plain-zeta MCMC Bayesian fitting path onto the shared
Stage 6 model/posterior-engine architecture introduced in Stage 6.3. It is the
second proof of the architecture, following the numerical zeta migration in
Stage 6.4.

## Scope

The stage:

- implements the MCMC posterior-engine protocol for `zetaModel`;
- keeps the existing independent-uniform-proposal Metropolis calculation;
- preserves the legacy random-number generation order for reproducibility under
  the same `set.seed()`;
- preserves posterior-mean shape, chain variance, fitted probabilities, stored
  chain, and spline posterior density in the returned `psFit`;
- adds the internal MCMC posterior representation alongside legacy components;
- routes `fitDistBayes()` through `zetaModel()` and `mcmcPosteriorEngine()`;
- adds deterministic P- and S-data regression tests against an independent copy
  of the pre-migration algorithm;
- leaves all ZIZ fitting paths unchanged.

## Architecture

`fitDistBayes()` becomes orchestration rather than an algorithm container:

1. construct `zetaModel()`;
2. construct `mcmcPosteriorEngine()`;
3. obtain a posterior representation via `fitPosterior()`;
4. obtain posterior summaries and the point estimate through engine methods;
5. obtain fitted probabilities through `modelProbabilities()`;
6. assemble the existing `psFit` shape while retaining the engine
   representation for later consolidation.

The MCMC engine stores the retained chain, posterior mean, variance matrix and
run metadata. Model-specific likelihood evaluation remains on `zetaModel`.

## Compatibility decision

The 1.0.7 Bayesian API is not treated as a published compatibility constraint,
but Stage 6.5 deliberately preserves current fit values and legacy `psFit`
components so this architecture migration can be validated independently of
later API consolidation.

## Deferred work

Stage 6.5 does not migrate ZIZ MCMC fitting and does not redesign final
`psPosterior` construction. Those remain later Stage 6 tasks after both plain
zeta engines have demonstrated the shared protocol.
