# fitPS Stage 6.1 code review and Bayesian architecture audit

## Purpose

Stage 6 is the fitPS maintainability, documentation, and Bayesian architecture series. It starts from the completed Stage 5.12.1 codebase (`Version: 1.0.7.018`) and uses the version series `1.0.8.xxx`.

The current CRAN release is 1.0.6. The Bayesian work developed during 1.0.7 has not been exposed in the literature and should be treated as provisional. Stage 6 therefore preserves the established pre-1.0.7 public API, but it does not preserve provisional 1.0.7 Bayesian function names, internal structures, or dispatch patterns merely for compatibility.

Stage 6.1 is a design and audit stage. It makes no intentional R behaviour changes. Its purpose is to define the architecture before code is reorganised.

## Executive summary

The code review identifies two related maintainability problems.

First, internal documentation is incomplete. A source scan found 138 named functions in `R/*.R`, with 71 lacking an immediately preceding roxygen2 block. The gaps are concentrated in newer Bayesian, plotting, bootstrap, validation, and numerical helper code. Stage 6 should document every function, exported or not, and add rationale comments where the mathematics or numerical choices are not obvious from the code.

Second, the Bayesian implementation has grown along two dimensions at once: statistical model and posterior approximation method. Plain zeta and zero-inflated zeta (ZIZ) duplicate validation, observation handling, fitted-probability construction, posterior summarisation, and `psFit` assembly, while the ZIZ implementation additionally duplicates method-specific dispatch across numerical integration, MCMC, Laplace, and importance sampling.

The central Stage 6 recommendation is therefore not merely to consolidate ZIZ functions. It is to separate **model responsibilities** from **posterior-engine responsibilities** and make both explicit S3 abstractions.

The target design is:

```text
psModel
    |
    +-- zetaModel
    +-- zizModel
    +-- futureDistributionModel

psPosteriorEngine
    |
    +-- numericalPosteriorEngine
    +-- mcmcPosteriorEngine
    +-- laplacePosteriorEngine
    +-- importancePosteriorEngine
```

A Bayesian fit becomes a composition of a model and an engine. The model determines what the parameters mean and how probabilities and likelihoods are computed. The engine determines how the posterior is approximated and represented.

This is deliberately close to a C++ strategy/object-oriented design while remaining idiomatic within fitPS's existing S3 style.

## 1. Documentation audit

### 1.1 Document every function

All functions should receive roxygen2 documentation, including internal helpers. Internal documentation need not create a separate public help topic for every helper. Related helpers should use `@keywords internal` and shared `@rdname` topics where that improves navigation.

The current source contains 138 named functions, with 71 lacking an immediately preceding roxygen2 block. Important concentrations include:

- posterior plotting helpers;
- Bayesian option and validation helpers;
- bootstrap helpers;
- posterior inflation helpers;
- numerical, Laplace, and importance-sampling machinery;
- internal construction and conversion helpers.

Stage 6.2 should establish a complete documentation baseline before structural refactoring begins.

### 1.2 Internal comments

Comments should explain decisions, not restate syntax. The following areas deserve explicit rationale comments:

- rescaling posterior densities before numerical integration;
- construction and use of integration weights;
- working versus natural parameterisations, especially `(eta, tau)` transformations;
- Jacobian terms introduced by parameter transformations;
- Laplace Hessian/covariance calculations;
- importance proposal construction, log weights, and normalisation;
- why `Pr(pi < epsilon | data)` is the relevant practical inflation diagnostic rather than `Pr(pi = 0 | data)` under a continuous prior;
- any numerical clipping used only to avoid floating-point boundary failures.

## 2. The architectural problem

The current code mixes model identity and posterior method in function names and orchestration. Examples include plain-zeta functions such as `fitDistBayes()` and `fitDistBayesIntegrate()` alongside ZIZ functions for MCMC, numerical integration, Laplace approximation, and importance sampling.

This produces a matrix that tends to grow by multiplication:

```text
                      posterior engine
                 numerical  mcmc  laplace  importance
zeta                 yes     yes     -        -
ZIZ                  yes     yes    yes      yes
future model          ?       ?      ?        ?
```

A procedural design encourages one fitting function for every populated cell. That is the duplication Stage 6 should avoid.

The abstraction should also avoid pretending that every model must support every engine. The matrix is sparse by design. Adding a common architecture is about sharing contracts and plumbing, not forcing unsupported statistical combinations into existence.

## 3. Model abstraction

### 3.1 Model responsibilities

A `psModel` object should own statistical-model behaviour that is independent of the posterior approximation algorithm. At minimum, the design should provide model-level operations for:

- validating model parameters and prior compatibility;
- converting `psData` into the observation representation required by the likelihood;
- evaluating the log likelihood, and where useful the log posterior contribution supplied to an engine;
- defining natural parameter names and constraints;
- defining working-parameter transformations where an engine needs unconstrained coordinates;
- converting parameter values into P/S probability vectors;
- defining parameter summaries needed by the common posterior object;
- declaring which posterior engines are supported.

The exact generic names can be settled during implementation. The architectural boundary matters more than the final spelling.

### 3.2 Zeta model

The plain zeta implementation currently duplicates work across `fitDistBayes()` and `fitDistBayesIntegrate()`, including:

- repeated `psData` validation;
- repeated P/S observation conversion;
- repeated prior-range validation;
- repeated zeta likelihood construction;
- repeated fitted probability construction and term naming;
- repeated `psFit` assembly;
- separate storage conventions for MCMC chains and numerical posterior density functions.

A `zetaModel` should make these common operations single-source.

### 3.3 ZIZ model

The ZIZ implementation has already factored out some model-level helpers such as `zizLogLikelihood()` and probability helpers, but engine implementations still repeat:

- fitted ZIZ probability construction;
- beta-prior parameter validation;
- parameter-summary assembly;
- common `psFit` fields;
- translation between engine-specific representation and common posterior summaries.

These should move behind the common model contract where appropriate.

### 3.4 Future distributions

The model abstraction should intentionally leave the door open for alternative P/S distributions without attempting to design them in Stage 6.

The design goal is that adding a new distribution should primarily require:

1. a new model constructor/class;
2. model-specific likelihood and probability methods;
3. model parameter validation and transforms;
4. tests establishing which existing engines it supports.

It should not require copying a complete Bayesian fitting stack.

This reduces the marginal cost of future model development while avoiding an over-general framework built for distributions that do not yet exist.

## 4. Posterior-engine abstraction

### 4.1 Engine responsibilities

A `psPosteriorEngine` object should own algorithm-dependent behaviour. Candidate engine classes are:

```text
numericalPosteriorEngine
mcmcPosteriorEngine
laplacePosteriorEngine
importancePosteriorEngine
```

The engine should not be named after zeta or ZIZ. Numerical integration and MCMC are posterior approximation strategies, not model identities.

### 4.2 Common protocol

A strong initial protocol is conceptually:

```text
fitPosterior(engine, model, data, prior, ...)
summarisePosterior(engine, model, representation, ...)
posteriorDiagnostics(engine, representation, ...)
posteriorPointEstimate(engine, model, representation, ...)
```

The exact methods may change during Stage 6.3 after the simplest zeta numerical implementation is used as a proof of design.

`fitPosterior()` should return an engine-specific representation, not a complete `psFit` object.

`summarisePosterior()` should translate that representation into the common `psPosterior` contract.

`posteriorDiagnostics()` should return algorithm-specific diagnostics or `NULL`.

`posteriorPointEstimate()` is useful because engines may naturally differ in which representative parameter value is most appropriate for the top-level fit. For example, an MCMC implementation naturally provides posterior means, while a Laplace approximation naturally centres on a posterior mode.

### 4.3 Representation objects

The representations should remain different when the mathematics is different. A numerical grid, MCMC sample, Laplace approximation, and weighted importance sample should not be forced into one physical data structure.

They can receive lightweight internal classes so their contracts are explicit, for example:

```text
numericalPosteriorRepresentation
mcmcPosteriorRepresentation
laplacePosteriorRepresentation
importancePosteriorRepresentation
```

The common `psPosterior` object then presents a uniform user-facing summary while retaining the engine-specific representation internally.

## 5. Model versus engine ownership

The most important Stage 6 rule is to keep the two axes separate.

The **model** owns:

- likelihood mathematics;
- parameter meaning and constraints;
- parameter transformations required by the model;
- conversion from parameters to P/S probabilities;
- model-specific diagnostics such as the meaning of the ZIZ inflation parameter.

The **engine** owns:

- numerical integration;
- MCMC transitions and retained draws;
- Laplace optimisation/Hessian approximation;
- importance proposal sampling and weighting;
- engine-specific convergence or numerical diagnostics;
- representation-specific extraction of moments, intervals, or draws.

The **fit orchestration layer** owns:

- selecting/constructing the model;
- selecting/constructing the posterior engine;
- invoking the common protocol;
- constructing the final `psPosterior` and `psFit` objects;
- preserving the established pre-1.0.7 fitPS public interface where required.

## 6. Inflation probability and other model-specific derived quantities

The original audit proposed consolidating:

- `numericalInflationProbability()`;
- `mcmcInflationProbability()`;
- `importanceInflationProbability()`;
- `laplaceInflationProbability()`.

They should still be consolidated, but the revised architecture changes where that responsibility sits.

Inflation is a ZIZ model concept, while the operation needed to compute it depends on posterior representation. Therefore the final design should avoid both extremes:

- a central `switch()` over engine names; and
- treating inflation as a mandatory method on every posterior engine, including models that have no inflation parameter.

A better design is a generic representation-level posterior probability operation, or a ZIZ model method that delegates representation-specific probability evaluation. Conceptually:

```text
posteriorProbability(representation, parameter = "pi", condition = pi < epsilon)
```

The exact API should be kept simpler than this pseudo-signature if possible. The key requirement is polymorphic dispatch without making ZIZ-specific semantics part of the base engine contract.

## 7. Common posterior object for zeta and ZIZ

Bayesian plain-zeta fits and Bayesian ZIZ fits should both use `psPosterior` consistently.

Currently, plain-zeta Bayesian functions store objects such as `chain` and `pdf` directly on `psFit`, while the newer ZIZ implementation attaches a `psPosterior` object containing parameter summaries, probability summaries, engine representation, diagnostics, method, level, and model.

Stage 6 should remove this asymmetry.

The target is:

```text
psFit
    |
    +-- posterior : psPosterior
                       |
                       +-- method
                       +-- parameters
                       +-- probabilities
                       +-- representation
                       +-- diagnostics
                       +-- model
```

This should make posterior printing, summaries, fitted probabilities, credible intervals, and plotting consistent across zeta and ZIZ.

## 8. Consolidation opportunities across zeta and ZIZ

### 8.1 Observation conversion and support validation

Both model families repeatedly transform P observations by adding one and leave S observations unchanged. This should have one well-documented helper or model method.

The repeated checks ensuring there is information above the minimum support should likewise have one source of truth where the semantic rule is common.

### 8.2 Probability construction and term naming

Both zeta and ZIZ calculate probability vectors and then apply the same P/S term naming convention. The distribution-specific probability calculation should remain model-specific, while naming and output-shape logic should be shared.

Existing `zetaProbabilities()` and `zizProbabilities()` provide a natural basis for this consolidation.

### 8.3 Prior validation

General validation of prior object structure/range should be shared. Model-specific prior requirements should remain with the model. ZIZ beta-prior shape validation should be consolidated into one helper rather than repeated by individual engines.

### 8.4 Parameter summaries

`makeZizParameterSummary()` is currently ZIZ-specific, while plain zeta constructs parameter and variance fields differently. Stage 6 should define a common parameter-summary contract and let each model provide its parameter names and model-specific values.

### 8.5 Common `psFit` construction

Bayesian engines should not construct complete `psFit` objects independently. Once model and engine results conform to common contracts, a shared finaliser should construct the top-level Bayesian fit.

This is one of the largest opportunities to prevent future duplication.

## 9. Things that should remain separate

Consolidation is not an objective by itself. Keep separate code where mathematical distinctions are real.

Examples include:

- zeta versus ZIZ likelihoods;
- numerical grid construction versus MCMC sampling;
- Laplace Hessian calculations versus importance weighting;
- one-dimensional versus multidimensional numerical representations where common code would obscure rather than clarify;
- natural-to-working transforms that genuinely differ between models;
- model-specific diagnostics such as ZIZ inflation.

A good Stage 6 result may contain many small functions. The goal is clear ownership and reuse, not the smallest possible function count.

## 10. API policy

Stage 6 should use two compatibility rules.

1. Established functionality present in the CRAN 1.0.6 API should remain compatible unless separately and explicitly deprecated.
2. Bayesian interfaces introduced during 1.0.7 may be renamed, reorganised, or removed when the Stage 6 design provides a cleaner replacement.

This means provisional functions such as `fitDistBayes()`, `fitDistBayesIntegrate()`, and the method-specific ZIZ Bayesian implementations do not need permanent compatibility wrappers solely because they existed in a development version.

The eventual user-facing API should preferably select a Bayesian method through the established fitting functions and a coherent Bayesian option object rather than expose the internal engine matrix directly.

## 11. Recommended implementation sequence

### Stage 6.1 - audit and architecture definition

- add this audit under `dev/`;
- reserve the `1.0.8.xxx` version series, beginning at `1.0.8.001`;
- define the two-axis model/engine architecture;
- make future-distribution extensibility an explicit design objective;
- renumber the previously proposed Bayesian-bootstrap Stage 6.1 context to future Stage 7.1;
- make no R behaviour changes.

### Stage 6.2 - complete documentation baseline

- document every function with roxygen2;
- use grouped internal topics where appropriate;
- add rationale comments for mathematically or numerically non-obvious implementation choices;
- regenerate `NAMESPACE` and `man/*.Rd` with `devtools::document()`;
- run strict tests and package check;
- make no intentional statistical changes.

### Stage 6.3 - shared model and engine protocols

- introduce lightweight internal S3 constructors/classes for `psModel` and `psPosteriorEngine`;
- define the minimum generics and explicit default failure methods;
- define the representation contract and supported-engine declaration;
- add shared observation/probability naming/validation helpers where they are already demonstrably duplicated;
- add deterministic protocol tests;
- do not migrate all fitting implementations in this stage.

### Stage 6.4 - numerical zeta proof of architecture

- migrate `fitDistBayesIntegrate()` through `zetaModel + numericalPosteriorEngine`;
- return the common `psPosterior` representation;
- characterise numerical results against the current implementation;
- adjust the abstractions if this simplest case feels forced before continuing.

### Stage 6.5 - MCMC zeta migration

- migrate the plain-zeta MCMC implementation to `zetaModel + mcmcPosteriorEngine`;
- consolidate duplicated zeta validation, likelihood, fitted probabilities, and result construction;
- verify agreement with the current implementation under deterministic seeds.

### Stage 6.6 - numerical ZIZ migration

- migrate numerical ZIZ fitting to `zizModel + numericalPosteriorEngine`;
- share numerical-engine infrastructure with zeta where the mathematics permits;
- preserve ZIZ-specific multidimensional integration logic behind model/representation methods.

### Stage 6.7 - MCMC ZIZ migration

- migrate ZIZ MCMC to the common MCMC engine protocol;
- consolidate common MCMC representation, summaries, and diagnostics where appropriate;
- keep ZIZ-specific proposals or block-update behaviour explicit if they are genuinely model dependent.

### Stage 6.8 - Laplace and importance ZIZ migration

- migrate the existing ZIZ-only Laplace and importance engines;
- do not add plain-zeta Laplace or importance fitting merely for symmetry;
- use the supported-engine contract to make unavailable combinations explicit.

### Stage 6.9 - common posterior orchestration and derived quantities

- make both zeta and ZIZ Bayesian fits construct the same `psPosterior` contract;
- consolidate top-level `psFit` finalisation;
- replace method-name switches for posterior operations with S3 dispatch;
- refactor ZIZ inflation probability onto the new representation/model dispatch path;
- harmonise posterior summaries, fitted probabilities, diagnostics, and plotting.

### Stage 6.10 - remove superseded provisional Bayesian paths

- delete duplicated 1.0.7 Bayesian implementations after replacement tests are established;
- remove compatibility branches that exist only for unpublished provisional 1.0.7 interfaces;
- update documentation to describe the final architecture rather than migration history.

### Stage 6.11 - final regression and architecture audit

- run comprehensive deterministic regression tests across supported model/engine combinations;
- compare key posterior estimates and P/S probability summaries with characterised pre-refactor results;
- run `devtools::document()`, strict tests, and strict package check;
- audit all R files for remaining undocumented functions, avoidable method switches, dead compatibility code, and duplicated model/engine plumbing;
- confirm the design remains simple enough to support a future alternative distribution without changing the engine abstraction.

## 12. Testing strategy for the refactor

The architecture refactor must be protected by characterisation tests before old implementations are removed.

For each migrated model/engine combination, tests should cover:

- parameter estimates or posterior moments within appropriate numerical tolerance;
- fitted P/S probabilities;
- posterior probability summary columns and term names;
- credible intervals;
- engine-specific diagnostics;
- deterministic MCMC behaviour under a fixed seed where practical;
- error handling for invalid priors, parameter values, unsupported engines, and malformed representations;
- consistent `psPosterior` structure for zeta and ZIZ.

Tests should remain deterministic and offline.

## 13. Design constraints for Stage 6

Stage 6 should resist three forms of overengineering.

First, do not introduce R6 or another object-system dependency when S3 dispatch is sufficient.

Second, do not build a universal distribution framework in anticipation of unknown future models. Define the smallest `psModel` contract that cleanly supports zeta and ZIZ, and verify that a future distribution could implement the same contract without copying the Bayesian stack.

Third, do not force common code where dimensions or mathematics are genuinely different. Shared interfaces are valuable even when concrete implementations remain separate.

## Conclusion

Stage 6 should be understood as a package-wide Bayesian architecture and maintainability stage, not a ZIZ-only cleanup.

The principal design is the composition of two orthogonal abstractions:

```text
Bayesian fit = psModel + psPosteriorEngine
```

For the current package, `zetaModel` and `zizModel` provide the model dimension, while numerical integration, MCMC, Laplace approximation, and importance sampling provide the engine dimension. A common `psPosterior` object presents their results consistently.

This design removes current duplication, gives plain zeta and ZIZ a coherent Bayesian workflow, and leaves a practical path for future alternative distributions without making Stage 6 depend on hypothetical future requirements.
