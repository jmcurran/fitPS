# Stage 9.8 built-in MCMC consolidation audit

## Purpose

Stage 9.8 checks whether the built-in zeta, logarithmic, and zero-inflated zeta models still need model-specific MCMC implementations now that fitPS has a generic model-neutral MCMC engine.

## Finding

All three built-in model descriptors satisfy the public Bayesian mathematics contract and can be sampled successfully by `fitMcmcPosteriorModel.psModel()`. The generic sampler therefore does not require model-name switches or built-in-only mathematics.

The specialised samplers are nevertheless retained for compatibility in Stage 9.8.

- The zeta regression suite explicitly requires the specialised sampler to reproduce the historical chain, fitted values, posterior density, and RNG sequence exactly for a fixed seed.
- The ZIZ regression suite explicitly requires exact reproduction of the historical chain, posterior probabilities, fitted values, covariance, acceptance rates, and RNG sequence.
- The logarithmic sampler follows the same long-standing compatibility pattern and is retained rather than changing its stochastic output solely to reduce source size.

This means the specialised methods are compatibility adapters, not architectural requirements. External models and future built-in models can use the generic MCMC engine without implementing their own sampler.

## Architectural decision

Stage 9.8 does not remove the specialised built-in samplers. Their presence does not prevent Bayesian model extensibility because the generic `psModel` fallback is complete and independently tested.

A future major or explicitly compatibility-breaking release may remove these methods after announcing the change. Until then, exact historical stochastic behaviour is considered sufficient value to justify keeping them.

## Consequence for Stage 9 completion

The MCMC architecture is considered consolidated at the contract level:

- fitPS owns the generic MCMC algorithm;
- models own likelihoods, priors, controls, constraints, and transformations;
- built-in specialised samplers are optional compatibility overrides;
- downstream models do not need sampler implementations.
