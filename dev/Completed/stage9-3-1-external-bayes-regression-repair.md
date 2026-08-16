# Stage 9.3.1 external Bayesian regression repair

Stage 9.3 intentionally changed the public external-model Bayesian path. An
external model that advertises the MCMC engine is no longer rejected merely
because it is external. Instead, generic Bayesian fitting proceeds through the
public model contract and requires the downstream model author to supply an
explicit model-specific prior.

The existing public-model regression test still asserted the pre-Stage 9.3
error message. Stage 9.3.1 updates that stable behavioural test to assert the
new contract: external Bayesian fitting without a prior fails clearly because
an explicit model-specific prior is required.

No posterior algorithms or model mathematics are changed by this repair.
