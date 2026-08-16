# Stage 10.2.1 sparse-support Bayesian fitting audit

## Motivation

The Pettard S survey contains six observations and every observed group contains
one fragment. Its aggregated support is therefore a single occupied value,
`S1`, with count six.

The previous shared observation conversion rejected every survey having fewer
than two occupied support values. That made the data unavailable to both MLE
and Bayesian engines before the estimation method could decide whether the
support pattern was estimable.

## Statistical distinction

A singleton occupied support value does not imply that the observations contain
no information. For Pettard S data, the likelihood strongly favours models that
place high mass on one fragment. For zeta and logarithmic models the likelihood
is driven toward a parameter boundary, so the historical MLE rejection remains
appropriate. A proper finite-support prior, however, makes the corresponding
Bayesian posterior proper and allows the data to update the prior.

For the zero/one-inflated zeta model, the six S1 observations similarly favour
parameter combinations assigning high probability to one. The two parameters
are weakly separated by these data alone, making prior sensitivity especially
important, but this is an inferential limitation rather than invalid input.

## Stage 10.2.1 decision

- `psObservationData()` now performs observation mapping only.
- The historical singleton-support error is applied by MLE fitting rather than
  by shared observation conversion.
- Bayesian fitting proceeds for singleton occupied support with an explicit
  warning that posterior inference may be strongly influenced by the prior.
- Pettard S data are used as the real regression case.
- No Bayesian Bootstrap functionality is introduced in this stage.

This keeps data representation neutral while allowing estimation engines to
make method-specific decisions.

## Repair note

The initial Stage 10.2 validation attempt exposed test-harness issues rather than an implementation failure: numeric survey support was compared with integer identity, and `expect_warning()` return values were assigned in place of the fitted objects. Stage 10.2.1 corrects those regression tests while preserving the intended sparse-support implementation.
