# Stage 6.2.1: roxygen summary-method repair

## Purpose

Stage 6.2.1 repairs the documentation-only Stage 6.2 attempt after `R CMD check`
reported undocumented `x` arguments in the generated `summary.psBootstrap.Rd` and
`summary.psPosterior.Rd` topics.

Stage 6.2 had deliberately attached the corresponding `print.summary.*` methods to
the summary help topics using `@describeIn`. Roxygen therefore added the print
method usages, whose first argument is `x`, to those shared Rd files. The shared
roxygen topics documented the summary methods' `object` arguments but did not yet
document the print methods' `x` arguments.

## Repair

- Add `@param x An object of class summary.psBootstrap` to the
  `print.summary.psBootstrap()` roxygen block.
- Add `@param x An object of class summary.psPosterior` to the
  `print.summary.psPosterior()` roxygen block.
- Preserve the complete Stage 6.2 documentation baseline because Stage 6.2 failed
  before the commit step.
- Make no executable R-code changes.

## Validation

The repair must run the full fitPS validation workflow because the change set
contains files under `R/`. In particular, `devtools::document()` must regenerate
both shared summary Rd topics and the subsequent strict package check must confirm
that their `x` arguments are documented.

The Stage 6.2 attempt consumed version `1.0.8.002`; the first Stage 6.2.1 attempt
therefore advances the package to `1.0.8.003` in the user's current working tree.
