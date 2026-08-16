# Stage 9.5.1 external engine regression repair

Stage 9.5 intentionally extended the external Poisson model from MCMC-only Bayesian support to both numerical and MCMC posterior engines.

The existing public-model regression test still expected `supportedPosteriorEngines(model)` to return only `"mcmc"`, so the Stage 9.5 strict test suite failed after the version bump.

This repair updates that stale expectation to `c("numerical", "mcmc")`. No numerical-posterior mathematics or engine behaviour is changed by the repair.
