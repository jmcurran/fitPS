# fitPS future development context: Bayesian bootstrap

## Purpose

This note records a possible future extension of fitPS based on Rubin's Bayesian bootstrap. It follows the Stage 5 work that introduced parallel `psPosterior` and `psBootstrap` uncertainty objects and posterior/bootstrap summaries for derived P and S probabilities.

This is a development context only. It does not commit fitPS to a Bayesian-bootstrap API or implementation.

## Statistical motivation

The ordinary nonparametric bootstrap repeatedly resamples observations and refits the model. Equivalently, each resample assigns multinomial counts, normalised to empirical weights, to the observed sample points. Rubin's Bayesian bootstrap instead keeps the observed support points fixed and draws continuous random weights

```text
(w1, ..., wn) ~ Dirichlet(1, ..., 1).
```

Each draw defines a random discrete distribution supported on the observed sample. A statistic or fitted-model functional can then be recomputed under those weights. For fitPS the corresponding derived probability calculation would be

```text
weighted data -> weighted zeta/ZIZ fit -> P or S probabilities
```

repeated over Bayesian-bootstrap weight draws.

The result is Bayesian in interpretation, but it differs from the current parametric Bayesian analysis in fitPS. The current `psPosterior` describes uncertainty in explicit zeta/ZIZ model parameters under stated priors. A Bayesian-bootstrap object would instead describe uncertainty induced by a nonparametric posterior on the empirical support weights.

## Relationship to existing fitPS architecture

The Stage 5 architecture makes this extension substantially easier than it would previously have been. The following components should be reusable:

- `zetaProbabilities()` and `zizProbabilities()` for derived P/S probabilities;
- the structured probability summary format used by `psPosterior` and `psBootstrap`;
- percentile/quantile summary helpers where mathematically appropriate;
- the print, summary, fitted, prediction and plotting conventions established for uncertainty objects;
- seeded replicate workflows and failed-fit diagnostics.

The main new requirement is weighted fitting.

## First technical question: weighted likelihoods

Before implementing a Bayesian bootstrap, audit whether every relevant fitPS likelihood can accept arbitrary positive observation weights without altering its statistical meaning. This should include at least:

1. ordinary zeta fitting for P and S surveys;
2. zero/one-inflated zeta fitting;
3. starting-value calculations;
4. likelihood and log-likelihood helpers;
5. any optimisation or profile-likelihood helpers that currently assume integer frequency counts;
6. handling of survey data already stored in aggregated frequency form.

For aggregated survey data, Dirichlet weights conceptually belong to the individual surveyed units, not merely to distinct count categories. An efficient implementation may aggregate those individual Dirichlet weights by observed count after drawing them, but the observational unit must remain clear.

## Suggested object design

Do not overload `psPosterior` or `psBootstrap`. If implemented, use a distinct class such as

```r
class(fit$bayesianBootstrap)
# "psBayesianBootstrap"
```

A possible structure is

```text
psBayesianBootstrap
  parameters
  probabilities
  weights / replicate representation
  level
  diagnostics
  seed
```

This preserves the statistical distinction among:

```text
psPosterior          parametric Bayesian posterior
psBootstrap          frequentist nonparametric bootstrap
psBayesianBootstrap  Rubin Bayesian bootstrap
```

## Public API possibilities

Do not choose the final public names until weighted fitting has been audited. A likely interface would parallel existing functions, for example:

```r
fit = bayesianBootstrapFit(fit, B = 2000, seed = 1234)
bayesianBootstrapProbs(fit)
plot(fit$bayesianBootstrap)
```

`fitted.psFit()` could eventually support an explicit type such as `"bayesianBootstrapMean"`, but only if that does not make the `type` interface unwieldy.

## Interpretation

Keep the Bayesian-bootstrap interpretation separate from both existing frameworks. In particular:

- ordinary bootstrap means estimate the centre of a bootstrap sampling distribution;
- parametric Bayesian posterior means average over `p(theta | x)`;
- Bayesian-bootstrap means average over the posterior distribution on empirical support weights.

Do not describe the Bayesian bootstrap as simply another way to obtain the current zeta/ZIZ parameter posterior.

## Suggested future stages

### Stage 6.1 - weighted-likelihood audit

- inspect all fitting paths for arbitrary positive observation-weight support;
- identify where integer counts are assumed;
- design representation of individual versus aggregated survey weights;
- produce no public API yet.

### Stage 6.2 - weighted fitting infrastructure

- add internal weighted likelihood support;
- characterise equivalence to current behaviour when all weights are equal;
- add deterministic tests for weighted fits.

### Stage 6.3 - Bayesian-bootstrap object

- generate Dirichlet weight draws;
- refit under each draw;
- construct parameter and P/S probability summaries;
- retain diagnostics and seeded reproducibility.

### Stage 6.4 - public API and documentation

- expose a coherent user-facing interface;
- add plotting and prediction integration only after the statistical object is stable;
- compare Roux results across plug-in, frequentist bootstrap, parametric Bayes and Bayesian bootstrap.

## References

Rubin, D. B. (1981). The Bayesian Bootstrap. *The Annals of Statistics*, 9(1), 130-134. DOI: 10.1214/aos/1176345338.

Efron, B. (2012). Bayesian inference and the parametric bootstrap. *The Annals of Applied Statistics*, 6(4), 1971-1997.
