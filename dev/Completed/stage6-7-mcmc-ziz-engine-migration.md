# Stage 6.7 - MCMC ZIZ engine migration and test naming cleanup

## Purpose

Stage 6.7 migrates the zero-inflated zeta Metropolis-Hastings fitting path onto the shared model/posterior-engine architecture introduced in Stage 6.3. It also removes development-stage numbering from the test suite so tests are named by stable behaviour rather than by implementation history.

## Scope

- Add `fitMcmcPosteriorModel.zizModel()` as the ZIZ implementation of the shared MCMC engine contract.
- Preserve the established ZIZ sampler, including its RNG ordering, beta proposal, shape proposal, burn-in handling, acceptance diagnostics, posterior summaries, and probability summaries.
- Reduce `fitZIDistBayes()` to model/engine orchestration and `psFit` assembly.
- Store the common `mcmcPosteriorRepresentation` on the returned fit while retaining the raw MCMC chain as the representation used by `psPosterior` for current probability and inflation workflows.
- Keep numerical, Laplace, and importance ZIZ mathematics unchanged.
- Rename stage-numbered test files to stable component/behaviour names and remove stage numbers from test descriptions.
- Remove the now-obsolete test asserting that ZIZ MCMC is not yet migrated.

## Behaviour-preservation target

For fixed inputs and seed, the migrated ZIZ MCMC path must reproduce the pre-migration implementation for:

- the complete retained MCMC chain;
- posterior means of `pi` and `shape`;
- the posterior variance-covariance matrix;
- per-parameter acceptance fractions;
- fitted P/S probabilities;
- posterior probability summaries; and
- the raw chain and acceptance diagnostics exposed through `psPosterior`.

The migration intentionally preserves the existing RNG sequence. The update-type vector and acceptance uniforms are drawn before the iteration-specific beta/uniform proposals, matching the previous implementation.

## Test naming policy

Tests are repository behaviour specifications, not stage records. After this stage:

- no test filename may contain a development stage number;
- no `test_that()` description may refer to a development stage number;
- redundant transition-only assertions should be removed once the transition is complete.

Development history remains recorded under `dev/` and in `NEWS.md`, not in test names.

## Deferred work

Later Stage 6 work will migrate the remaining ZIZ Laplace and importance engines and then consolidate shared result construction, posterior summaries, and inflation/probability dispatch where that reduces duplication without obscuring method-specific mathematics.
