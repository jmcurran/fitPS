# fitPS Stage 2 Context: Shape Parameterisation Consistency

## Current project state

The repository is `jmcurran/fitPS`.

Stage 1 has focused on stabilisation and package hygiene:

- Added baseline `testthat` infrastructure.
- Added tests for core `psData` and `psFit` workflows.
- Added or repaired tests around fitting, prediction, and `logLik.psFit()`.
- Fixed start-value handling for `fitDist()` and `fitZIDist()`.
- Fixed zero-inflated zeta prediction to use the model object's `pi`.
- Added `NEWS.md`.
- Added `README.md`.
- Updated `.gitignore` and `.Rbuildignore` to exclude generated archives, stage runners, development files, build markers, and generated vignette artifacts.
- Stage runners use versions of the form `1.0.4.xxx`, incrementing the final numeric component before validation so every build attempt consumes a build number.

The usual workflow is controlled staged development:

- Work from the current branch, usually `test-infrastructure` or the next controlled stage branch.
- Generate a `fitPS_stageX_changes.zip` patch.
- Generate a matching external `run_stageX.sh` runner.
- Put both files in Downloads.
- Run from the repo root with `scripts/runStage.sh X`.
- Stage runners should remain external and should not be copied into the package root.
- Strict validation should fail on notes:
  ```r
  devtools::check(
    build_args = c("--no-build-vignettes"),
    args = c("--no-manual", "--ignore-vignettes", "--no-tests"),
    error_on = "note"
  )
  ```

## Major aim of the next stage

The next major aim is **parameterisation consistency** across the entire package.

The standard Riemann zeta / zeta distribution uses a parameter usually written as:

```text
alpha
```

with:

```text
alpha > 1
```

VGAM, however, uses a parameter named `s` or `shape` where:

```text
s = alpha - 1
```

There is inconsistency in the current code:

- Some package code appears to use `alpha` as if it were VGAM's `s`.
- Some code uses `s` or `sPrime` to represent `alpha - 1`.
- Some user-facing code or documentation may call the parameter `shape`, but it is not always clear whether that means standard zeta `alpha` or VGAM `shape`.

## Desired package-wide convention

The whole package should use the **original/statistical zeta definition** internally and at the user API level, but call it:

```text
shape
```

That means:

```text
shape = alpha
shape > 1
```

For any call into VGAM, convert explicitly:

```r
vgamShape = shape - 1
```

or another clearly named internal variable.

The user should specify:

```r
shape > 1
```

The user should never need to know VGAM's shifted parameterisation.

VGAM calls should receive:

```r
shape - 1
```

not `shape`.

## Naming rules for the refactor

Use these names consistently:

- `shape`: the user-facing and package-facing zeta parameter, equal to standard `alpha`, constrained to `shape > 1`.
- `vgamShape`: the shifted VGAM parameter, equal to `shape - 1`.
- Avoid using `s`, `sPrime`, `s.prime`, or `alpha` in user-facing code unless documenting the mathematical relationship.
- If a helper truly needs the shifted value, name it `vgamShape` or another explicit name that communicates the conversion.
- If documentation mentions alpha, say clearly:
  ```text
  fitPS uses `shape` for the standard zeta parameter alpha, with shape > 1. VGAM uses a shifted parameter internally, equal to shape - 1.
  ```

## R style requirements

Follow James Curran's R style:

- Use `=` for assignment.
- Use camelCase identifiers.
- Always use braces for control structures, even single-line bodies.
- Use roxygen2 documentation.
- Do not manually edit `NAMESPACE`; regenerate through `devtools::document()`.
- Prefer `@importFrom` rather than `pkg::fun()` in package code where practical.
- Keep changes targeted.
- Preserve backward compatibility where possible.
- Add deterministic offline `testthat` tests for all behavior changes.

## Suggested stage breakdown

### Stage 2.1: Inventory and tests

Goal: identify every place where zeta shape parameterisation appears and add tests that lock the desired convention.

Suggested tasks:

1. Search for:
   - `shape`
   - `alpha`
   - `sPrime`
   - `s.prime`
   - `VGAM`
   - `dzeta`
   - `pzeta`
   - `qzeta`
   - `rzeta`
   - `zeta`
   - `zetaff`

2. Classify each occurrence as one of:
   - user-facing `shape = alpha`;
   - internal standard-shape calculation;
   - VGAM shifted-shape calculation;
   - documentation only.

3. Add tests verifying:
   - user-specified `shape` must be `> 1`;
   - fitted `psFit$shape` is standard `alpha`, not VGAM's shifted shape;
   - VGAM-backed probability calculations use `shape - 1`;
   - `probfun()`, `predict.psFit()`, `fitted.psFit()`, and random generation helpers agree with the standard-shape convention;
   - zero-inflated zeta workflows also report and use standard `shape`.

Do not do the full refactor until the inventory and tests make the current inconsistencies visible.

### Stage 2.2: Refactor core fitting and prediction

Goal: make core fitting and prediction consistently use `shape = alpha`.

Likely files to inspect and update:

- `R/fitDist.R`
- `R/fitZIDist.R`
- `R/fitlogDist.R` if shape-like terminology appears there
- `R/predict.psFit.R`
- `R/fitted.psFit.R`
- `R/probfun.R`
- `R/rzeta.R`
- `R/rZIzeta.R`
- `R/internalfunctions.R`

Expected implementation pattern:

```r
shape = object$shape
vgamShape = shape - 1
```

For VGAM calls, pass `vgamShape`.

### Stage 2.3: Documentation and NEWS

Goal: make documentation unambiguous.

Update roxygen and man pages so:

- `shape` always means standard zeta `alpha`.
- `shape > 1`.
- VGAM's shifted parameterisation is mentioned only where necessary.
- NEWS includes a version-numbered entry summarising the consistency repair.
- README gets a short note if needed.

## Important caveats

- This is a semantic consistency refactor, not just a rename.
- Be careful with saved objects, printed output, comparisons, confidence intervals, and bootstrap/credible intervals.
- Preserve public API names where possible.
- If existing code allowed `shape <= 1` because it was treating `shape` as VGAM's shifted parameter, that should be considered a bug under the new convention.
- Any backward-compatibility handling should warn clearly if legacy shifted values are detected, but avoid ambiguous automatic conversion unless tests define the intended behavior.

## Stage runner expectations

For any generated stage runner:

- Name it `run_stage2_1.sh`, `run_stage2_2.sh`, etc.
- Name the patch `fitPS_stage2_1_changes.zip`, etc.
- Runner stays external in Downloads.
- Do not copy the runner into the package root.
- Patch zips should contain only package files.
- Use strict staged validation.
- Bump `DESCRIPTION` version before validation.
- Keep `1.0.4.xxx` versioning until explicitly changed.
- Update `NEWS.md` after validation and before commit when the stage has user/developer-facing changes.
