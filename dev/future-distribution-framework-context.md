# fitPS future development context: distribution-pluggable models

## Purpose

fitPS is currently organised around zeta and zero/one-inflated zeta distributions, but the surrounding infrastructure has become more general than those models. The long-term opportunity is to separate the package's survey/inference machinery from the specific probability distribution so users can switch among suitable discrete distributions.

This is a strategic context note, not a proposal to refactor the package immediately. The Bayesian and bootstrap work should remain stable before a distribution abstraction is introduced.

## Long-term goal

A future fitPS model could be described by three largely independent layers:

```text
survey semantics    P versus S indexing and interpretation
distribution model  zeta, alternative count distributions, inflated variants
inference engine    MLE, bootstrap, Bayesian posterior, Bayesian bootstrap
```

The public workflow should ideally allow the distribution to change without rewriting fitted-value, prediction, bootstrap, posterior, plotting and derived P/S probability code.

## Why the current Stage 5 work helps

Recent work already separates several concepts that were previously entangled:

- fitted models are represented by `psFit`;
- posterior uncertainty is represented by `psPosterior`;
- bootstrap uncertainty is represented by `psBootstrap`;
- derived P/S probabilities are calculated from parameter representations rather than hard-coded only in presentation methods;
- `fitted()` and `predict()` distinguish plug-in, posterior-mean and bootstrap-mean quantities.

These are useful foundations for a distribution-pluggable architecture.

## Candidate role for the distributional package

The CRAN package `distributional` provides vectorised distribution objects and generic operations for density/mass, CDF, quantiles, random generation, means, variances, intervals and related summaries. Its design goal is to let modelling packages return distribution objects rather than exposing only raw parameter vectors. That makes it worth evaluating as part of a future fitPS abstraction.

Potential advantages include:

- a common distribution object interface;
- vectorised density/CDF/quantile/random-generation behaviour;
- an ecosystem-oriented representation rather than a fitPS-only distribution class;
- easier downstream interaction with predictive distributions.

However, `distributional` does not by itself solve the fitPS modelling problem. fitPS still needs to know:

- which parameters are free and their constraints;
- how to optimise or integrate a likelihood;
- how P and S survey indices map onto distribution support;
- how inflation or hurdle components are represented;
- how priors are assigned to model parameters;
- how starting values are obtained;
- how each model behaves under bootstrap and Bayesian inference.

Therefore `distributional` should initially be evaluated as a possible **distribution representation layer**, not assumed to be the fitting framework.

## Possible internal abstraction

A future internal model specification might expose operations conceptually like

```text
validateParameters(parameters)
logLikelihood(data, parameters, weights = NULL)
probabilityTerms(parameters, surveyType, n)
randomSample(n, parameters)
startingValues(data)
parameterConstraints()
makeDistribution(parameters)
```

The final function could return a `distributional` object where an appropriate distribution class exists. Custom fitPS distribution classes or adapters may still be required for zeta and inflated-zeta models.

The key design rule should be that inference engines depend on this model interface rather than directly on zeta-specific functions.

## Mixture and inflation models

Inflated models should not be implemented as one-off special cases for every base distribution. A longer-term design should consider a compositional model representation, for example:

```text
base distribution + inflation specification
```

rather than separate monolithic implementations such as `zeta` and `zizeta`. This would make it easier to ask whether another base distribution should also have a zero/one-inflated form.

## Compatibility requirements

Any future refactor must preserve:

- current zeta and ZIZ numerical results;
- P and S survey semantics;
- existing public calls where practical;
- current `psFit`, `psPosterior`, and `psBootstrap` user workflows;
- package-level reproducibility and offline tests;
- ability to fit models without forcing users to understand the internal distribution adapter.

The first implementation should wrap the existing zeta model behind the new interface and prove exact behaviour preservation before adding another distribution.

## Suggested exploratory programme

### Phase A - inventory

- list every zeta-specific assumption in fitting, prediction, simulation, bootstrap, posterior and plotting code;
- identify the minimum model/distribution interface;
- evaluate whether `distributional` can represent the required zeta and inflated-zeta objects directly or through extension classes.

### Phase B - adapter around existing zeta

- introduce the internal abstraction without adding a new public distribution choice;
- require regression equivalence with the current package.

### Phase C - second distribution

- choose one scientifically plausible alternative distribution;
- implement it through the same interface;
- test whether fitted, bootstrap and Bayesian workflows genuinely remain distribution-agnostic.

### Phase D - public model selection

Only after the abstraction works for at least two distributions should fitPS expose a public `distribution = ...` or model-specification argument.

## Questions to answer before implementation

1. What alternative distributions are scientifically plausible for the forensic P/S applications?
2. Must candidate distributions have support on all positive integers, or are truncated/support-shifted models useful?
3. How should P versus S indexing be represented independently of a particular distribution?
4. Should inflation be a generic wrapper around a base distribution?
5. Can priors be specified generically across models with different parameter dimensions?
6. Should fitPS depend on `distributional`, suggest it, or merely provide adapters compatible with it?
7. How much of `psFit` should remain stable if model parameters cease to be simply `shape` and optional `pi`?

## Current recommendation

Do not add `distributional` as a dependency yet. First perform an architecture audit and prototype an adapter around the existing zeta implementation. The package is a promising candidate because its vectorised distribution-object philosophy aligns with the direction fitPS is moving, but fitPS requires a model-fitting contract above that representation layer.

Useful current documentation: https://cran.r-project.org/package=distributional and https://pkg.mitchelloharawild.com/distributional/.
