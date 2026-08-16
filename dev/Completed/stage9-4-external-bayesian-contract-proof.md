# Stage 9.4 external Bayesian contract proof

## Purpose

Stage 9.4 completes the multi-parameter external-model proof for the generic
MCMC architecture introduced in Stage 9.3. The Poisson-normal model remains
entirely outside package implementation code and supplies only the public model
mathematics contract.

## Scope

- Exercise posterior parameter summaries for a two-parameter external model.
- Exercise generic posterior probability summaries derived from retained MCMC
  draws through `modelProbabilities()`.
- Verify that fitted probabilities are derived from the posterior point estimate
  through the external model contract.
- Verify DIC for a multi-parameter external MCMC posterior.
- Verify serialization preserves the originating external model descriptor,
  posterior representation, summaries, and fitted values.
- Verify matched P- and S-survey data produce the same posterior chain and
  probability sequence when they represent the same zero-based support.

## Engine policy

MCMC remains the primary public extension engine for external Bayesian models.
Laplace and importance remain available where already implemented but are not
expanded or promoted by this stage. Numerical posterior fitting is likewise
left unchanged because its current one-dimensional integration design is not a
necessary part of the Stage 9 extensibility proof.

## Stage 9 status after this stage

After Stage 9.4, the public external-model architecture has been exercised for:

- one- and two-parameter Bayesian models;
- constrained parameter transformations and Jacobians;
- model-specific multivariate prior specifications;
- fitPS-owned MCMC sampling;
- posterior parameter and probability summaries;
- DIC;
- P/S support mapping; and
- serialization and restoration.

The remaining Stage 9 work should therefore focus on consolidation and
regression auditing rather than adding more posterior engines unless a concrete
user-facing need emerges.
