# Stage 9.1 Bayesian model-extensibility audit

## Purpose

Stage 9.1 audits the Bayesian fitting architecture left by Stage 8 and defines the minimum public model contract needed for external `psModel` objects to use fitPS-owned posterior engines without `fitPS:::` access or package modification.

Stage 8 achieved this separation for maximum-likelihood fitting: fitPS owns the generic optimisation while an external model owns its likelihood and probability mathematics. Stage 9 should establish the same boundary for Bayesian fitting.

## Current architecture

The current architecture already contains several useful pieces of the intended design:

- `psModel()` is public and external model subclasses can be created outside fitPS.
- `modelParameterNames()`, `modelObservationData()`, `modelProbabilities()`, `modelLogLikelihood()`, `modelMleControl()`, and `supportedPosteriorEngines()` are public model-level extension points.
- Posterior engines are represented independently by `psPosteriorEngine` subclasses.
- `validateEngineModelPair()` already rejects model/engine combinations that the model does not declare as supported.
- Each engine has common representation, summary, diagnostic, and point-estimate machinery.
- MCMC and importance probability summaries already call `modelProbabilities()` generically once posterior draws are available.

These are strong foundations. The missing layer is generic construction of the posterior representation itself.

## Primary blocker: the engines delegate fitting back to the model

Each generic engine currently performs only validation and then dispatches to a model-specific fitting method:

- numerical -> `fitNumericalPosteriorModel()`
- MCMC -> `fitMcmcPosteriorModel()`
- Laplace -> `fitLaplacePosteriorModel()`
- importance -> `fitImportancePosteriorModel()`

The default `psModel` methods stop with "not implemented" errors. Built-in zeta, ZIZ, and logarithmic models therefore implement the actual integration, MCMC, Laplace, or importance algorithm themselves.

This is model-dispatch rather than model-generic Bayesian fitting. An external model author can declare an engine in `supportedPosteriorEngines()`, but that declaration is insufficient: the author must also reproduce the corresponding sampler or approximation algorithm. That violates the Stage 9 principle that fitPS should own posterior orchestration while the external model owns only its mathematics.

## Secondary blocker: the current prior object is not a general model contract

`makePrior()` and `validateBayesPrior()` assume a single finite two-element `range` and one `logd` function. This representation is useful for the legacy one-parameter zeta/logarithmic paths, but it is not a sufficient universal prior protocol.

In particular:

- a two-parameter external Poisson-normal model cannot naturally describe its complete prior using one two-element range;
- parameter-specific or correlated priors are awkward or impossible to express through the current validation rule;
- a model may need to interpret a legacy `psPrior` differently from a general external prior specification;
- engine code should not need to know the physical representation chosen for a model's prior.

The generic engine therefore needs a model-level operation that evaluates the complete log prior for a named natural-scale parameter vector.

## Tertiary blocker: transformations and Jacobians are built-in mathematics

Laplace and importance code for ZIZ uses private helpers such as natural-to-working transformations, inverse transformations, and the log absolute Jacobian. Those are model mathematics, not engine mathematics.

A generic engine needs public model-level access to:

- natural -> working parameter transformation;
- working -> natural parameter transformation;
- log absolute Jacobian of the inverse transformation.

The default transformation should be identity with zero log-Jacobian so unconstrained models such as a simple Poisson log-rate parameter do not need boilerplate.

## Additional blocker: Bayesian starts and bounds are model-specific

Built-in posterior methods currently obtain starting values, finite integration ranges, proposal ranges, and related controls inside distribution-specific methods. Generic engines need a uniform way to ask a model for Bayesian fitting controls without knowing the model class or parameter names.

These controls should be separate from the likelihood and prior mathematics and should be allowed to vary by engine.

## Minimum public Bayesian model contract

Stage 9 should add the smallest set of public model generics needed to make the existing engines model-neutral. The recommended contract is below.

### Existing public methods to retain

1. `modelParameterNames(model)`
2. `modelObservationData(model, x)`
3. `modelLogLikelihood(model, parameters, data, ...)`
4. `modelProbabilities(model, parameters, n, type, ...)`
5. `supportedPosteriorEngines(model)`

No replacement is needed for these Stage 8 methods.

### New method: `modelLogPrior()`

Recommended shape:

```r
modelLogPrior(model, parameters, prior, ...)
```

Responsibilities:

- receive a complete named natural-scale parameter vector;
- return one finite log-density value or `-Inf` outside prior support;
- interpret whatever prior object/specification the model supports;
- permit multivariate and correlated priors;
- provide backward-compatible built-in methods that delegate to existing `psPrior` mathematics where appropriate.

This method should be the only prior-density operation required by generic engines.

### New method: `modelBayesControl()`

Recommended shape:

```r
modelBayesControl(model, x, engine, prior, ...)
```

It should return a small validated named list containing only controls that an engine needs, for example:

- `start`: named natural-scale starting vector;
- optional natural or working bounds where mathematically required;
- engine-specific tuning defaults that genuinely depend on the model.

Engine-owned controls such as requested iteration count, burn-in, number of importance draws, numerical tolerances, and random seed should remain engine arguments rather than becoming part of the model contract.

### New transformation methods

Recommended generics:

```r
modelToWorking(model, parameters, ...)
modelFromWorking(model, working, ...)
modelWorkingLogJacobian(model, working, ...)
```

Defaults for `psModel` should implement the identity transformation and a zero log-Jacobian. Models with constrained natural parameters override them.

The working representation should be a named numeric vector with length equal to `modelParameterNames(model)`. Generic engines can then operate in unconstrained Euclidean space when needed.

### Optional method: `modelPosteriorSupport()`

A separate public support/bounds method should only be added if implementation experience shows that `modelBayesControl()` becomes overloaded. Stage 9 should not add it pre-emptively.

## Generic posterior kernel owned by fitPS

With the methods above, fitPS can own one internal model-neutral posterior evaluator:

```text
working parameters
    -> modelFromWorking()
    -> modelLogLikelihood()
    +  modelLogPrior()
    +  modelWorkingLogJacobian()
    -> log posterior on working scale
```

This kernel should become the mathematical boundary consumed by MCMC, Laplace, and importance engines. A natural-scale variant can be used by simple numerical integration where appropriate.

External model authors should not have to construct `psPosteriorRepresentation` objects directly.

## Engine-by-engine implications

### MCMC

MCMC is the best first proof of genuine extensibility because a generic random-walk Metropolis implementation needs only:

- parameter names;
- starting values;
- the generic log-posterior kernel;
- proposal/tuning controls owned by the engine.

The current zeta MCMC implementation is effectively already a model-specific version of this pattern. Generalising this engine should be the first implementation target after the public contract exists.

### Laplace

Laplace can become generic once the engine can optimise the working-scale log posterior and map the resulting mode/covariance back to natural parameters. The ZIZ-specific transform helpers currently embedded in this path should move behind the public model transformation methods.

### Importance

Importance sampling can use the same working-scale posterior kernel and generic Gaussian proposal machinery used by Laplace. It should therefore follow naturally after generic Laplace support.

### Numerical

The current numerical engine is based on direct one-dimensional integration and is not naturally a universal multi-parameter integrator. It should remain explicitly restricted to suitable models rather than forcing a misleading generic implementation.

A model can declare numerical support when its requirements are satisfied. Initially, a generic numerical implementation may reasonably support only one-parameter models with finite natural-scale bounds supplied by `modelBayesControl()`.

## Probability-summary implications

MCMC and importance posterior probability summaries are already close to model-neutral because they evaluate `modelProbabilities()` on posterior samples.

Two remaining areas need attention:

- Laplace probability summaries currently require `zizModel` and call ZIZ-specific proposal/transform helpers.
- numerical probability summaries dispatch back to `summariseNumericalPosteriorProbabilities.<model>()`.

These should eventually use generic posterior draws or a generic one-dimensional quadrature representation so external models do not need model-specific summary methods.

This is part of Bayesian extensibility, but it can be implemented after generic fitting itself is proven.

## `fit()` orchestration gap

`fitModel.psModel()` currently rejects every method except `"mle"`. Stage 9 must add the external-model Bayesian path there (or in a thin helper called from it):

1. normalise `method` / `bayesOptions`;
2. resolve and validate the requested engine against `supportedPosteriorEngines()`;
3. resolve the prior without assuming the legacy scalar `psPrior` is universal;
4. ask the model for Bayesian controls;
5. call the fitPS-owned posterior engine;
6. build the established `psPosterior` / `psFit` result while retaining `modelObject`.

Built-in compatibility wrappers can continue to delegate through their established paths while they are progressively migrated onto the same generic contract.

## Recommended compatibility strategy

Do not remove the existing `fit*PosteriorModel.<built-in>` methods immediately. Instead:

1. add the new public model contract and generic engine implementations;
2. prove the generic path using external Poisson;
3. migrate built-in models engine by engine to the generic path while retaining regression tests;
4. remove or reduce model-specific engine methods only after behavioural equivalence is demonstrated.

This limits Stage 9 risk and avoids a large simultaneous rewrite.

Legacy `psPrior` objects should remain accepted by built-in models. `makePrior()` does not need to become a universal multivariate prior constructor in the first Stage 9 implementation.

## Proposed implementation sequence after Stage 9.1

### Stage 9.2 - public Bayesian mathematics contract

- add/export `modelLogPrior()`;
- add/export `modelBayesControl()`;
- add/export identity transformation generics and defaults;
- add contract validation and deterministic tests;
- retain current built-in fitting behaviour.

### Stage 9.3 - generic MCMC engine and external Poisson proof

- implement a fitPS-owned generic MCMC path using the public contract;
- connect external-model Bayesian orchestration through `fit()`;
- demonstrate external Poisson fitting without `fitPS:::`;
- prove clear unsupported-engine errors and model retention/serialization.

### Stage 9.4 - multi-parameter proof with Poisson-normal

- exercise transformations, vector parameters, starts, and multivariate prior evaluation;
- fit the external Poisson-normal model through generic MCMC;
- retain numerical integration only inside the model likelihood, not inside the posterior engine.

### Later Stage 9 substages

- generalise Laplace and importance using the shared working-scale posterior kernel;
- generalise posterior probability summaries;
- decide whether a restricted generic numerical engine adds sufficient value;
- migrate built-in Bayesian methods onto the generic path with regression coverage.

The exact substage boundaries can be adjusted after each completed archive is reviewed.

## Proof criteria for Stage 9 completion

Stage 9 is complete when an external model defined entirely outside `R/` can, without `fitPS:::`:

1. construct a public `psModel`;
2. provide likelihood, prior, transform, and Bayesian-control methods;
3. declare a supported fitPS posterior engine;
4. call `fit(..., method = "bayes", ...)` or the equivalent supported engine selection;
5. receive a normal `psFit`/`psPosterior` result retaining the external model descriptor;
6. use posterior summaries and derived probabilities where the engine supports them;
7. serialize and restore the fit successfully.

Poisson and Poisson-normal remain test/vignette models only and must not be added under `R/`.

## Stage 9.1 decision

The recommended Stage 9 architecture is therefore:

**fitPS owns a generic posterior kernel and the posterior algorithms; models expose likelihood, prior, transformation, and minimal control mathematics through public S3 methods.**

The existing engine-specific model-fitting generics are the principal architectural seam to replace for the external-model path. The existing engine descriptors, engine-support declaration, posterior representations, and Stage 8 likelihood/probability contract should be retained and built upon rather than redesigned.
