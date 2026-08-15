# Stage 8.1 public model API audit

## Purpose

Stage 8 converts the internal `psModel` architecture into a supported public extension API. The target is stronger than adding another distribution inside fitPS: a user or another R package must be able to define a model outside the fitPS namespace and fit it without rebuilding, modifying, or forking fitPS.

Stage 8 uses the `1.1.0.xxx` development series. The first Stage 8 attempt is `1.1.0.001`; every later attempted build consumes the next final component.

The proof-of-concept targets for the completed public API are:

- a Poisson model defined outside fitPS;
- a Poisson-normal model defined outside fitPS.

Stage 8.1 is an audit/design stage. It deliberately does not implement either model and does not guess the precise statistical definition of the Poisson-normal model. That definition should be fixed before its implementation test is written.

## Current internal model protocol

`R/psModel.R` currently provides these model-level operations:

- `newPsModel()` constructs an internal descriptor with a model identifier, natural parameter names, supported posterior engines, and a concrete S3 subclass;
- `modelParameterNames()` obtains the natural parameter names;
- `supportedPosteriorEngines()` declares supported Bayesian engines;
- `supportsPosteriorEngine()` tests one model/engine pairing;
- `modelObservationData()` maps `psData` to model observation support;
- `modelProbabilities()` computes fitted P/S probabilities;
- `modelLogLikelihood()` computes the model log likelihood.

The built-in descriptors are `zetaModel`, `zizModel`, and `logarithmicModel`.

This is a coherent internal protocol, but it is not yet a supported extension API. The constructor and generics are internal and documented with `@noRd`, and external users have no public supported route for constructing a model and passing it into a common fitter.

## Current generic infrastructure that already dispatches through models

The Stage 6/7 architecture already uses model dispatch in important places:

- Bayesian finalisation obtains parameter names through `modelParameterNames()` and probabilities through `modelProbabilities()`;
- likelihood-based comparison uses `modelLogLikelihood()` and `modelParameterNames()`;
- `probfun()` evaluates probabilities through `modelProbabilities()`;
- DIC uses the model likelihood together with posterior representations;
- posterior engines validate declared model/engine support.

This means Stage 8 should expose and complete the existing architecture rather than create a parallel fitting system.

## Remaining built-in assumptions

### 1. `modelFromFit()` is hard coded

`R/modelComparison.R` reconstructs a model from `object$model` with an explicit switch over `zeta`, `ziz`, and `logarithmic`/`log`.

This blocks genuine third-party models. A fitted external model cannot later participate in `logLik()`, deviance, DIC, `probfun()`, or model-aware plotting unless fitPS already knows its character identifier.

### 2. Fitted objects retain only a model name

Current fit construction stores `model = model$model` rather than retaining the originating descriptor. Generic operations therefore have to reconstruct the model after fitting.

The preferred Stage 8 direction is for every newly created `psFit` to retain the originating model object, for example in a dedicated `modelObject` component, while preserving the existing character `model` component for compatibility and display.

`modelFromFit()` can then become compatibility logic:

1. return the retained model object when present;
2. fall back to the legacy built-in-name resolver for older `psFit` objects.

This avoids requiring a registry for object-based fitting and preserves old serialized fits as far as practical.

### 3. Public fitters still own model-specific orchestration

`fitDist()`, `fitZIDist()`, and `fitlogDist()` each perform model-specific orchestration. Bayesian branches are closer to the common architecture than the MLE branches, but there is not yet one public model-oriented entry point.

The long-term API should have a single common fitting path. Existing functions remain public convenience wrappers.

### 4. MLE does not yet have a model-level fitting contract

The current public model protocol sketch has likelihood and probability operations, but no generic declaration of:

- starting values;
- parameter bounds or transformations;
- parameter validation;
- generic MLE optimisation policy;
- model-specific MLE fitting when generic optimisation is insufficient.

A third-party Poisson model is a useful test because it should be possible to fit it with a small general contract. The Poisson-normal target is a stronger test because it is likely to require more than a single scalar parameter and may expose where a generic optimizer needs model-supplied bounds, starts, transformations, or an optional custom fitting method.

### 5. Some posterior implementation remains model-specific

The posterior-engine matrix is intentionally sparse. Exporting `supportedPosteriorEngines()` does not make every existing engine automatically usable by an arbitrary model.

For Stage 8 the public contract should distinguish:

- model methods required for MLE and common fitted-object behaviour;
- model methods required for a particular existing Bayesian engine;
- optional capabilities that may legitimately be absent.

External model support must not imply support for every Bayesian engine.

### 6. DIC still contains one concrete model branch

`posteriorExpectedDeviance.numericalPosteriorRepresentation()` has a special `inherits(model, "zizModel")` branch for the two-dimensional numerical posterior grid. This is an internal posterior-representation limitation rather than a reason to expose ZIZ-specific behaviour publicly.

Stage 8 should avoid promising numerical-engine DIC for arbitrary multidimensional external models until the representation contract is generic enough to support it. MLE AIC/BIC can be part of the initial external-model success criterion independently.

### 7. Prediction has zeta-specific interval logic

Plug-in point predictions already route through `probfun()` and therefore the model protocol. Profile/Wald interval handling in `predict.psFit()` still contains a zeta-specific branch.

The external public contract should initially promise generic point prediction/fitted probabilities. Model-specific interval procedures should remain optional capabilities unless a clean dispatch abstraction is introduced deliberately.

## Proposed public contract

Stage 8 should expose the smallest coherent contract, not every internal helper.

### Public constructor

Introduce a supported constructor, tentatively retaining the name `psModel()` or using a similarly explicit name, that creates a validated `psModel` descriptor from:

- a stable model identifier;
- natural parameter names;
- supported posterior engines;
- a concrete subclass supplied by the extension author.

The constructor should be sufficient for ordinary external use and should remove the need for `fitPS:::`.

Built-in constructors `zetaModel()`, `zizModel()`, and `logarithmicModel()` should become public if they are intended to be passed directly to the generic fitter.

### Required model generics

The minimum likely public S3 extension surface is:

- `modelParameterNames()`;
- `modelObservationData()`;
- `modelProbabilities()`;
- `modelLogLikelihood()`;
- `supportedPosteriorEngines()`.

The base `psModel` methods can satisfy declaration-style methods when values are stored directly in the descriptor. External subclasses should only need methods where model-specific mathematics or data mapping differs.

### MLE fitting contract

Stage 8.2/8.3 should add one deliberately small model-aware MLE contract. Two designs should be evaluated before implementation:

1. **Generic optimiser plus declarations**: models expose start values, parameter bounds/transformations, and likelihood; fitPS owns `optim()` orchestration.
2. **Model MLE generic**: a generic such as `fitModelMle(model, x, ...)` returns a small standard result consumed by common finalisation.

The preferred design should minimize the public surface while still supporting both the simple Poisson proof and the more demanding Poisson-normal proof. A hybrid is acceptable if generic defaults cover simple models and an optional model-specific method handles exceptional optimisation requirements.

### Fitted-object retention

New `psFit` objects should retain their originating model descriptor. The compatibility character field `model` should remain.

Generic code should ask `modelFromFit()` for the retained descriptor and should not depend on registration or hard-coded character resolution for new fits.

### Model registry

A registry is **not required** for the core external model API if users pass model objects directly:

```r
fit(x, model = myModel())
```

A registry may later support convenience character resolution:

```r
fit(x, model = "zeta")
```

Object-based third-party fitting must not depend on registration. This keeps the extension contract simple and avoids global mutable state as a prerequisite for fitting.

## Generic fitting interface

The recommended direction remains a concise public function:

```r
fit(x, model, method = c("mle", "bayes"), ...)
```

with Bayesian controls remaining compatible with the existing `bayesOptions` approach.

Before exporting literal `fit()`, Stage 8.2 should check namespace/import consequences, but the existence of similarly named functions in other packages is not by itself a reason to reject it.

The generic path should:

1. validate `x` and the model descriptor;
2. resolve `method` and Bayesian options;
3. validate model/engine support for Bayesian fitting;
4. delegate model-specific mathematics through the public model contract;
5. run common fit finalisation;
6. store both the compatibility model identifier and originating model object;
7. return an ordinary `psFit`.

`fitDist()`, `fitZIDist()`, and `fitlogDist()` should ultimately become thin wrappers around this path while preserving established arguments and behaviour.

## Poisson proof of concept

The Poisson extension should be defined entirely in test/helper code outside the fitPS namespace and should use public APIs only.

Its purpose is to prove the smallest external contract. At minimum it should demonstrate:

- external model construction;
- parameter declaration, probably a scalar rate parameter;
- external likelihood implementation;
- external P/S probability implementation;
- external observation-support mapping if the standard fitPS mapping is not statistically appropriate;
- MLE through `fit()`;
- returned `psFit` containing/recovering the external model object;
- `fitted()` and `predict(..., interval = "none")`;
- `logLik()`, AIC, BIC, and deviance;
- serialization round-trip of the fitted object when the model class/methods remain available.

The test must not call unexported fitPS functions or use `fitPS:::`.

## Poisson-normal proof of concept

Poisson-normal is the stronger extension test. Before implementation, Stage 8 should record the exact intended probability model and parameterization rather than silently infer one from the name.

Once defined, it should be implemented outside fitPS using the same public API as Poisson. It should test whether the contract remains adequate for a more complicated model, especially:

- multiple parameters;
- parameter domains/bounds or transformations;
- nontrivial likelihood evaluation;
- nontrivial P/S probabilities;
- generic model comparison;
- optional custom MLE behaviour if the general optimiser contract is insufficient.

If Poisson works but Poisson-normal requires fitPS source edits, Stage 8 has not met its success criterion.

## Proposed Stage 8 sequence

### Stage 8.1 - public model-contract audit and design

This document. No package code is changed. The key decisions are:

- use `1.1.0.xxx`, beginning with `1.1.0.001`;
- retain the originating model object in new fits;
- keep legacy character model identifiers for compatibility;
- do not require a registry for object-based fitting;
- expose only the model contract needed by external implementations;
- use external Poisson and Poisson-normal models as the proof of concept;
- do not begin the Bayesian bootstrap during Stage 8.

### Stage 8.2 - fitted-object model retention and generic fit foundation

Introduce the common model-oriented fitting path, model retention/fallback behaviour, and the initial `fit()` API. Migrate built-in paths only as far as needed to establish a stable foundation.

### Stage 8.3 - public extension contract and external Poisson

Export/document the minimal constructor/generics and define a Poisson model entirely outside the fitPS namespace in tests. Require public-only fitting, prediction/fitted probabilities, and MLE model comparison.

### Stage 8.4 - Poisson-normal stress test and MLE contract refinement

Implement the agreed Poisson-normal definition externally. Refine the public MLE contract only if this more demanding model reveals a genuine abstraction gap.

### Stage 8.5 - compatibility-wrapper consolidation and resolution audit

Make `fitDist()`, `fitZIDist()`, and `fitlogDist()` thin wrappers where safe, preserve legacy behaviour, and reduce remaining built-in model reconstruction assumptions.

### Stage 8.6 - public extension vignette

Turn the existing maintainer-oriented extension vignette into a true public third-party guide using public APIs only. Poisson is a likely compact teaching example; Poisson-normal can be discussed as the stronger validation example if its implementation would distract from the guide.

### Stage 8.7 - public API and regression audit

Audit built-in models, external models, comparison criteria, posterior support enforcement, serialization/reconstruction, examples, and absence of `:::` dependencies. Record any intentionally unsupported capabilities explicitly.

## Compatibility policy

- Preserve established pre-1.0.7 public behaviour unless a deliberate change is justified.
- Keep `fitDist()`, `fitZIDist()`, and `fitlogDist()`.
- Preserve the existing character `psFit$model` field while adding model-object retention.
- Continue to support older built-in fits that do not contain the retained model object through a compatibility fallback.
- Do not promise that arbitrary external models support every posterior engine, DIC representation, diagnostic, or interval method.

## Success criterion

Stage 8 is successful only if both proof models can be authored outside fitPS and used without rebuilding or modifying the package.

A user should be able to define a model with the documented public constructor and required S3 methods, pass that object to the generic fitter, receive a normal `psFit`, and use the generic facilities supported by that model. FitPS must not need a new hard-coded model name, source-file edit, or internal-namespace call for either proof model.
