# fitPS Stage 6.3 shared model and posterior-engine protocols

## Purpose

Stage 6.3 introduces the internal S3 protocols that separate statistical-model responsibilities from posterior-approximation responsibilities. It deliberately does not migrate any existing Bayesian fitting path yet; Stage 6.4 will use plain-zeta numerical integration as the first proof of the architecture.

## Model protocol

The stage adds internal `psModel` descriptors with concrete `zetaModel` and `zizModel` subclasses. Models declare their natural parameter names and the posterior engines they support. The current support matrix remains intentionally sparse: zeta declares numerical integration and MCMC, while ZIZ additionally declares Laplace and importance sampling.

The initial model protocol includes parameter-name discovery, supported-engine discovery, common P/S observation conversion, probability dispatch, and an explicit log-likelihood generic whose base method fails until each model path is migrated. Existing `zetaProbabilities()` and `zizProbabilities()` are used by the model probability methods rather than duplicated.

## Posterior-engine protocol

The stage adds internal `psPosteriorEngine` descriptors for numerical integration, MCMC, Laplace approximation, and importance sampling. It defines the intended strategy operations `fitPosterior()`, `summarisePosterior()`, `posteriorDiagnostics()`, and `posteriorPointEstimate()`.

No concrete fitting implementation is attached to these generics in Stage 6.3. The base methods validate the declared model/engine pairing and then fail explicitly with a not-yet-implemented message. This prevents the new architecture from silently changing existing Bayesian behaviour before migration tests are in place.

## Representation contract

`newPsPosteriorRepresentation()` introduces a lightweight common wrapper around an engine-specific payload. The payload remains unconstrained, so a numerical grid, MCMC chain, Gaussian approximation, and weighted sample retain their natural structures. The wrapper records the owning engine method and optional metadata and receives both an engine-specific representation subclass and the common `psPosteriorRepresentation` class.

## Consolidated P/S helpers

The stage centralises three pieces of behaviour already duplicated between zeta and ZIZ code:

- P/S term-index validation;
- mapping P indices onto the positive-integer latent support;
- standard P/S probability term naming.

It also introduces shared `psObservationData()` and makes the existing `zizObservationData()` delegate to it. Existing probability functions now use the shared helpers while retaining their distribution-specific mathematics.

## Compatibility and scope

This stage does not alter the established CRAN 1.0.6 API and does not route existing `fitDistBayes*()` or `fitZIDistBayes*()` functions through the new protocols. The new functions are internal and documented with roxygen2. Stage 6.4 will characterise and migrate the one-dimensional numerical zeta path before any wider engine migration.

## Validation

The change set adds deterministic offline tests covering model descriptors, the sparse engine support matrix, model probability dispatch, shared observation and term helpers, representation validation, and explicit failure for unsupported or not-yet-migrated model/engine combinations.
