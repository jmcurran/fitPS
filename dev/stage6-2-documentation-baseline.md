# fitPS Stage 6.2 documentation baseline

## Purpose

Stage 6.2 establishes the documentation baseline required before the Stage 6 Bayesian architecture refactor. It deliberately changes source documentation and explanatory comments only; it does not intentionally change package behaviour.

The baseline is the completed Stage 6.1.1 repository state at fitPS version 1.0.8.001. The Stage 6 runner will consume the next 1.0.8.xxx build number before validation.

## Function documentation

A source scan of `R/*.R` identifies 139 named function definitions. Before this stage, 71 had no immediately preceding roxygen2 block. Stage 6.2 adds internal roxygen2 documentation for all 71, including purpose, parameters, return value, and internal status.

A second quality pass also identified four functions whose existing roxygen blocks contained tags but no descriptive prose: `makeZizProposalDraws()`, `var.default()`, `print.summary.psBootstrap()`, and `print.summary.psPosterior()`. Stage 6.2 gives those functions explicit source descriptions while retaining their existing import, export, or shared-topic declarations.

Internal helpers use `@keywords internal` and `@noRd` where a separate installed help topic would not improve the public package manual. Existing exported functions and S3 methods retain their established public documentation topics.

The principal documentation concentrations are:

- Bayesian option validation and natural/working parameter transformations;
- plain-zeta Bayesian fitting helpers;
- zero-inflated zeta numerical, Laplace, and importance posterior helpers;
- posterior plotting and inflation-diagnostic helpers;
- bootstrap construction and summarisation helpers;
- zeta parameterisation and probability helpers;
- legacy likelihood, profile-likelihood, and model-comparison helpers.

## Rationale comments

Stage 6.2 adds comments only where the reason for the implementation is not obvious from the syntax. In particular, the source now explains:

- why `pi` and zeta shape are transformed to unconstrained `(eta, tau)` working coordinates;
- the Jacobian term required when the ZIZ posterior is represented on the working scale;
- clipping a numerical `pi` grid just inside its open `(0, 1)` support;
- maximum-log-posterior rescaling before exponentiation for stable numerical integration;
- rectangular grid integration weights and marginal-density scaling;
- inverse-Hessian covariance in the Laplace approximation and its delta-method transformation back to natural parameters;
- covariance inflation for the importance proposal, log-weight construction, normalisation, and effective use of weights;
- why `Pr(pi < epsilon | data)` is the practical zero-inflation diagnostic under a continuous prior rather than `Pr(pi = 0 | data)`;
- the corresponding threshold calculation for numerical-grid, MCMC, importance, and Laplace posterior representations;
- mode-based rescaling in the older one-dimensional zeta numerical-integration path.

Comments are intentionally sparse elsewhere: they explain mathematical or numerical decisions rather than restating code.

## Behaviour preservation

Stage 6.2 does not refactor or consolidate implementation code. In particular, it does not yet introduce `psModel` or `psPosteriorEngine` classes and it does not change the provisional Bayesian dispatch paths. Those changes remain Stage 6.3 and later work.

As a static behaviour-preservation check, all `R/*.R` files changed in this stage are compared with the Stage 6.1.1 baseline after removing full-line comments and blank lines. The executable source is expected to be identical.

Some older files do not follow the current fitPS formatting conventions. They are not reformatted in this stage because broad style-only changes would obscure the documentation-only diff and make behavioural review harder. Touched implementation code will be normalised when it is structurally refactored in later Stage 6 work.

## Validation expectation

Because this stage deliberately changes files under `R/`, the Stage 6.2 runner uses full package validation even though executable statements are unchanged. It should run:

1. change-set installation;
2. R library-path configuration;
3. the next `1.0.8.xxx` version bump;
4. stale build-path cleanup;
5. the R string-escape preflight when available;
6. `devtools::document()`;
7. strict `devtools::test()`;
8. strict `devtools::check()`;
9. version-numbered `NEWS.md` update;
10. controlled commit;
11. creation of `stage6_2_completed.zip` from committed `HEAD`;
12. package build and installation;
13. cleanup of temporary build-path files and older completed-stage archives while preserving `stage6_2_completed.zip`.

The completed archive is the handoff artifact for Stage 6.3.
