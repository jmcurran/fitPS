# Stage 5.11 final probability-inference audit

## Purpose

Stage 5.11 closes the Bayesian and bootstrap probability-inference development series. The final user-facing design distinguishes three quantities:

- plug-in probability: `P_k(thetaHat)`;
- bootstrap mean probability: `E*[P_k(thetaHat*)]`;
- Bayesian posterior mean probability: `E[P_k(theta) | x]`.

The worked documentation uses the real Roux et al. (2001) footwear survey supplied as `Psurveys$roux`.

## Public interfaces reviewed

- `fitDist()` and `fitZIDist()` retain existing frequentist defaults.
- Bayesian ZIZ fitting uses `method = "bayes"` and `bayesOptions$posteriorMethod`.
- Bayesian engines are `numerical`, `mcmc`, `laplace`, and `importance`.
- `posteriorProbs()` exposes posterior means, posterior SDs, and equal-tailed credible intervals.
- `bootstrapFit()` attaches a `psBootstrap` object to an MLE fit.
- `bootstrapProbs()` exposes bootstrap means, bootstrap SDs, and percentile confidence intervals.
- `fitted()` retains plug-in output by default and supports explicit `posteriorMean` and `bootstrapMean` types.
- `predict()` retains plug-in output by default and supports posterior credible and bootstrap percentile intervals when the corresponding uncertainty object is available.
- `probfun()` remains explicitly a plug-in probability function.
- `plot.psPosterior()` and `plot.psBootstrap()` display derived probability estimates with uncertainty intervals.

## Regression coverage

Existing staged tests cover:

- corrected component-wise Metropolis-Hastings sampling with non-uniform beta priors;
- seeded MCMC reproducibility and comparison with numerical integration;
- numerical-grid weighting for derived probabilities;
- weighted importance-sampling summaries;
- seeded Laplace probability summaries;
- P and S term naming and interval bounds;
- the distinction between plug-in and posterior-mean probabilities;
- bootstrap replicate transformation, failed-resample diagnostics, and reproducibility;
- `posteriorProbs()` and `bootstrapProbs()` S3 dispatch;
- `fitted()` and `predict()` semantics across plug-in, posterior, and bootstrap definitions;
- preservation of legacy plug-in and frequentist behaviour.

Stage 5.11 adds an end-to-end regression test using the real Roux survey for MLE, numerical Bayesian, and bootstrap probability workflows.

## Documentation coverage

Stage 5.11 adds a dedicated vignette, `bayesian-bootstrap-probabilities.Rmd`, using the Roux survey to demonstrate:

- Bayesian zero-inflated fitting;
- the `psPosterior` object;
- posterior P probabilities and credible intervals;
- the difference between plug-in and posterior-mean probabilities;
- all four Bayesian posterior engines;
- the beta prior on `pi` and shape prior;
- importance weights and the Laplace simulation approximation;
- the `psBootstrap` object;
- bootstrap mean probabilities and percentile intervals;
- side-by-side plug-in, bootstrap-mean, and posterior-mean interpretations;
- harmonised `predict()` behaviour.

The README contains a shorter Roux-based introduction and links users to the vignette for the complete workflow.

## Optional future work

The Bayesian/bootstrap probability programme is functionally complete after Stage 5.11 if the strict stage workflow passes. A later stage is only warranted if the final validation exposes a defect or if there is a separate decision to modernise legacy facilities such as `credint()` or consolidate older plotting helpers. Those are not required for the probability-inference architecture completed here.
