# Stage 9.6.4 Poisson-normal likelihood repair

Stage 9.6.4 repairs the external Poisson-normal proof model used to validate the generic two-dimensional numerical posterior engine.

The previous probability helper integrated directly over the latent Normal variable `z` using `dnorm(z, mean = mu, sd = sigma)`. When `sigma` was very small, the Normal density became so concentrated that `integrate()` could miss essentially all mass and return zero probabilities. That made the model log-likelihood `-Inf` before `cubature` was reached.

The repair uses the equivalent standard-Normal representation

`z = mu + sigma * u`, with `u ~ N(0, 1)`,

and integrates the Poisson probability against `dnorm(u)`. This removes the collapsing numerical spike while preserving the exact Poisson-normal mixture probability.

A targeted diagnostic on the failing regression case showed that the old helper returned zero probabilities and `-Inf` log-likelihood at the supplied start, while the standard-Normal formulation returned finite probabilities and allowed the generic two-dimensional `hcubature()` posterior fit to complete successfully.

The Stage 9.6.4 runner executes `test-generic-numerical-posterior.R` before the full test suite. Full package validation runs only after that targeted numerical regression passes.
