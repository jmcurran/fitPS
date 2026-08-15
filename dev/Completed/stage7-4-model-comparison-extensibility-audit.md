# Stage 7.4 model-comparison and extensibility audit

Stage 7.4 closes the Stage 7 model-comparison and distribution-extensibility
work. It reviews the zeta, zero-inflated zeta (ZIZ), and logarithmic-series
models against the architecture established in Stage 6 and exercised in Stages
7.1-7.3.

## Model-comparison contract

All three built-in maximum-likelihood models now participate in the same
comparison contract through `modelLogLikelihood()` and `logLik.psFit()`.
Consequently, deviance, AIC, and BIC are available without model-specific
criterion implementations. A durable regression test fits zeta, ZIZ, and
logarithmic models to the same survey data and verifies finite likelihood and
criterion values, common sample size, and the expected parameter counts of one,
two, and one respectively.

For Bayesian fits, DIC is calculated from the shared model deviance plus the
posterior representation. Numerical and MCMC representations support DIC for
zeta and logarithmic models. ZIZ additionally supports importance sampling.
Laplace fits deliberately reject DIC because the stored Laplace representation
does not currently provide a posterior deviance expectation. This sparse matrix
is intentional; Stage 7 does not add statistically meaningless engine/model
combinations merely for symmetry.

AIC/BIC and DIC answer different questions and should not be mixed casually.
AIC and BIC use the maximised likelihood and are therefore restricted to MLE
fits. DIC is restricted to Bayesian fits and can depend on the prior as well as
the likelihood. Comparisons should use the same observations and otherwise
comparable model specifications.

## Did logarithmic require inappropriate engine changes?

No model-name branches were added to the generic posterior-engine dispatch.
The logarithmic model supplies its mathematics through the existing model
methods and the secondary numerical/MCMC model dispatch points. Common P/S
support mapping, posterior finalisation, parameter summaries, probability
summaries, diagnostics, fitted values, and comparison criteria are reused.

The extension did expose two zeta assumptions in otherwise generic code:
scalar MCMC chains had been labelled `shape`, and one-parameter numerical DIC
had explicitly recognised `zetaModel`. Those assumptions were generalised in
Stage 7.2 using `modelParameterNames()`. This is evidence that the extension
exercise improved the abstraction rather than hiding incompatibilities behind
logarithmic-specific conditionals.

Some duplication remains between the zeta and logarithmic one-parameter
numerical and MCMC methods. It is currently small and located at the intended
secondary model-dispatch boundary. Stage 7.4 does not factor it out yet: two
examples are not enough evidence that a general one-parameter posterior helper
would remain appropriate for the next distribution. The duplication should be
revisited when another model exercises the same pathway.

## Internal extensibility result

The Stage 7 success criterion is substantially met for models maintained inside
fitPS. A new distribution can be represented by a small `psModel` subclass with
parameter declarations, probability and likelihood methods, supported-engine
declarations, and only the model-specific posterior methods required by those
engines. It can then reuse the common fitting/posterior infrastructure and
participate naturally in prediction, summaries, diagnostics, fitted
probabilities, and model comparison.

The logarithmic implementation therefore demonstrates a coherent internal
package-developer extension contract rather than another copied fitting stack.

## Third-party extension API audit

The current contract is not yet a supported plug-in API for downstream
packages. Several concrete barriers remain:

- `newPsModel()` and the built-in model constructors are internal;
- the key model generics are internal rather than exported as a supported API;
- `modelFromFit()` resolves built-in models through a hard-coded model-name
  switch;
- `fitDist()`, `fitZIDist()`, and `fitlogDist()` construct their own model
  descriptors rather than accepting an externally supplied model object;
- there is no public model registry, model validation contract, or namespace
  policy for third-party model classes; and
- supporting external posterior methods would require a documented stability
  promise for several currently internal representation and engine contracts.

Exporting the existing internals alone would therefore be insufficient and
would prematurely freeze implementation details that are still evolving.

## Recommendation for public extensibility

Do not expose the current internal machinery piecemeal in Stage 7. Instead,
treat third-party distributions as a later, deliberate API project. A sensible
future design should begin with a thin model-oriented fitting entry point such
as `fit(data, model = ..., ...)`, while retaining `fitDist()` and `fitZIDist()`
as compatibility/convenience wrappers.

That future work should decide, as one coherent contract:

1. how external model objects are constructed and validated;
2. which model generics are public and stable;
3. how supported posterior engines are declared;
4. whether models are passed directly, registered by name, or both;
5. which posterior representation interfaces third-party packages may rely on;
6. how model names, parameter names, and fitted-object reconstruction are
   resolved without hard-coded switches; and
7. what compatibility guarantees fitPS is willing to make for external models.

A direct model-object path is likely preferable to mandatory global
registration because it is explicit, testable, and naturally fits a future
`fit(data, model = ...)` interface. Named registration could be added later if
there is a real need for it.

## Remaining architectural limitations

- MLE orchestration is not yet unified behind a model-level fitting generic;
  the three public fitters still contain model-specific MLE orchestration.
- `modelFromFit()` is built-in-name based rather than reconstructing or storing
  a model descriptor on the fitted object.
- One-parameter numerical/MCMC posterior methods retain some duplicated
  mechanics between zeta and logarithmic models.
- Laplace DIC remains unavailable because the current representation does not
  retain the posterior information needed for a defensible expected deviance.
- The public API still exposes distribution-specific fitters rather than a
  single model-oriented front end.

None of these limitations undermines the internal Stage 7 extensibility result,
but they identify the natural next architectural boundary.

## Stage 7 conclusion

Stage 7 establishes common AIC, BIC, deviance, and DIC infrastructure and shows
that the logarithmic-series distribution can be integrated through the model
contract without copying unrelated posterior-engine machinery. The architecture
is now demonstrably extensible inside fitPS.

The next major API step, if desired, should be a deliberately supported public
model-extension layer and thin generic `fit()` front end, not further ad hoc
exports of internal functions.
