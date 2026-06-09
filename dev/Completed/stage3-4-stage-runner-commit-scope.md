# Stage 3.4 stage-runner commit-scope fix

## Purpose

Stage 3.4 stabilises the generated stage runner after Stage 3.3 exposed a controlled staging bug. The Stage 3.3 runner attempted to add all controlled paths in a single `git add` command and then masked failures with `|| true`. On repositories where optional controlled directories such as `inst/` or `vignettes/` were absent, that command could fail before staging newly installed files.

## Behavioural change

The Stage 3.4 runner stages controlled paths one at a time and only when they exist. This keeps the controlled commit scope while ensuring new package files under tracked controlled directories are added correctly.

The corrected pattern is:

```bash
add_existing_path() {
  local path="$1"
  if [[ -e "$path" ]]; then
    git add "$path"
  fi
}
```

The runner then calls `add_existing_path` for each controlled package path. Root-level Markdown and R project files are still added separately.

## Scope

This is a workflow-stabilisation stage. It does not change the public fitPS API.

## Validation intent

The corrected runner should:

- install `stage3_4_changes.zip` when requested;
- bump the `DESCRIPTION` build number before validation;
- run documentation, tests, and strict check;
- update `NEWS.md`;
- stage controlled paths safely, including new files;
- commit staged changes when present;
- create `stage3_4_completed.zip`;
- build and install the package;
- clean older completed stage archives after the current archive exists.
