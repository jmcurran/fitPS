# Future work: likelihood-matched conjugate priors for fitPS

## Purpose

Record a Bayesian-prior direction that arose during Stage 11 without expanding
Bayesian Bootstrap scope. This note is intentionally future work. It should be
audited mathematically and against the literature before implementation or
publication claims are made.

## Motivation

Bill Bolstad advocated choosing priors whose mathematical shape matches the
likelihood. For the ordinary zeta model this appears to lead to a particularly
simple and interpretable conjugate family.

The zeta probability mass function can be written as

```text
p(x | s) = exp{-s log(x) - log(zeta(s))},  s > 1.
```

For observations `x_1, ..., x_n`, the likelihood is therefore proportional to

```text
exp{-s sum(log(x_i))} zeta(s)^(-n).
```

A natural likelihood-matched prior is

```text
pi(s | m, mu) proportional to exp(-m mu s) zeta(s)^(-m),  s > 1,
```

or, with `g = exp(mu)`,

```text
pi(s | m, g) proportional to g^(-m s) zeta(s)^(-m).
```

The resulting posterior has the same form, with updates

```text
m* = m + n
m* log(g*) = m log(g) + sum(log(x_i)).
```

This suggests an intuitive interpretation:

- `m` behaves like a prior effective sample size;
- `g` behaves like a prior geometric mean;
- posterior updating averages prior and observed log-geometric information on
  the sufficient-statistic scale.

## Relevance to sparse all-one surveys

For an all-one survey,

```text
sum(log(x_i)) = 0,
```

and the ordinary zeta likelihood can drive the MLE toward `s = Inf`. A proper
likelihood-matched prior with `g > 1` supplies explicit prior information on
the sufficient-statistic scale and can give a proper finite posterior without
pretending that the MLE exists.

This could give fitPS a much more interpretable sparse-survey Bayesian story
than a generic prior directly on `s`.

## Propriety questions to establish carefully

Before implementation, establish the exact conditions for prior and posterior
propriety. The working tail argument is:

- as `s -> Inf`, `zeta(s) -> 1`, so `g > 1` gives exponential right-tail decay
  through `g^(-m s)`;
- near `s = 1`, `zeta(s)` diverges, so positive `m` appears to suppress rather
  than create a boundary-integrability problem.

These observations are useful motivation, not a substitute for a formal
proof.

## Relationship to earlier prior work

An earlier Bernardo-and-Smith-style route investigated by James required high
order derivatives of the Riemann zeta function, with derivative order tied to
sample size. Even if theoretically valid, that is unattractive as a practical
fitPS implementation boundary.

The likelihood-matched family above is computationally much simpler and has
hyperparameters with a direct data-scale interpretation. The associated paper
should compare these approaches carefully rather than assuming the simpler
family is preferable under every prior philosophy.

## Extensions to investigate

After establishing the ordinary zeta case:

1. Verify the conjugate-family derivation from exponential-family theory and
   identify the appropriate primary references.
2. Prove prior and posterior propriety conditions.
3. Decide whether `(m, g)` is the best public elicitation parameterisation.
4. Compare the conjugate prior with the priors currently supported by fitPS.
5. Study behaviour for all-one and otherwise sparse surveys.
6. Derive the corresponding structure, if any, for zero-inflated zeta.
7. Derive the corresponding structure, if any, for the logarithmic model.
8. Decide whether this should become an optional prior family, a recommended
   default, or remain primarily a paper-level result.

## Scope boundary

Do not mix this work into Rubin Bayesian Bootstrap implementation. Parametric
Bayesian prior uncertainty and Bayesian-Bootstrap uncertainty over empirical
observation weights answer different inferential questions and should remain
separate in the package and documentation.
