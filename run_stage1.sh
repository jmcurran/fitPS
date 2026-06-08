#!/usr/bin/env bash
set -euo pipefail

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  echo "Do not source this script. Run it with: bash ${BASH_SOURCE[0]}"
  return 1 2>/dev/null || exit 1
fi

fitps_root_override=""
start_step_number=1
total_steps=10
built_pkg_file=".fitps_built_pkg_path"
current_archive="stage1_completed.zip"

usage() {
  cat <<'USAGE'
Usage: bash run_stage1.sh [--start-step-number N|-sn N]

Runs the fitPS stage 1 package workflow. The version bump is performed before
validation so every build attempt consumes the next 1.0.4.xxx build number.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --start-step-number|-sn)
      if [[ $# -lt 2 ]]; then
        echo "Missing value for $1"
        exit 1
      fi
      start_step_number="$2"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1"
      usage
      exit 1
      ;;
  esac
done

case "$start_step_number" in
  ''|*[!0-9]*)
    echo "--start-step-number must be a positive integer"
    exit 1
    ;;
esac

if (( start_step_number < 1 || start_step_number > total_steps )); then
  echo "--start-step-number must be between 1 and $total_steps"
  exit 1
fi

should_run_step() {
  local step_number="$1"
  (( step_number >= start_step_number ))
}

resolve_fitps_root() {
  if [[ -n "$fitps_root_override" ]]; then
    printf '%s\n' "$fitps_root_override"
    return
  fi

  local script_dir
  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

  if [[ -f "$script_dir/DESCRIPTION" ]] && grep -q '^Package: fitPS$' "$script_dir/DESCRIPTION"; then
    printf '%s\n' "$script_dir"
    return
  fi

  if [[ -f "$(pwd)/DESCRIPTION" ]] && grep -q '^Package: fitPS$' "$(pwd)/DESCRIPTION"; then
    pwd
    return
  fi

  case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*)
      printf '%s\n' "D:/Dropbox/Code/git/fitPS"
      ;;
    Darwin*)
      printf '%s\n' "$HOME/Dropbox/Code/git/fitPS"
      ;;
    *)
      printf '%s\n' "$HOME/Dropbox/Code/git/fitPS"
      ;;
  esac
}

fitps_root="$(resolve_fitps_root)"

if [[ ! -f "$fitps_root/DESCRIPTION" ]] || ! grep -q '^Package: fitPS$' "$fitps_root/DESCRIPTION"; then
  echo "Could not resolve the fitPS package root. Edit fitps_root_override near the top of this script."
  exit 1
fi

if should_run_step 1; then
  echo "[1/$total_steps] Configuring R library paths"
  r_library_paths="$(${R_SCRIPT:-Rscript} -e 'cat(paste(normalizePath(.libPaths(), winslash = "/", mustWork = FALSE), collapse = .Platform$path.sep))')"

  if [[ -n "$r_library_paths" ]]; then
    export R_LIBS="$r_library_paths"
    export R_LIBS_USER="$r_library_paths"
  fi
else
  echo "Skipping step 1/$total_steps: Configuring R library paths"
fi

cd "$fitps_root"

if should_run_step 2; then
  echo "[2/$total_steps] Bumping DESCRIPTION version"
  Rscript -e 'lines = readLines("DESCRIPTION", warn = FALSE); versionIndex = which(startsWith(lines, "Version:")); if (length(versionIndex) != 1) stop("DESCRIPTION must contain exactly one Version line"); versionText = trimws(substring(lines[versionIndex], nchar("Version:") + 1)); parts = strsplit(versionText, ".", fixed = TRUE)[[1]]; if (length(parts) < 3) stop("Version must contain at least major.minor.patch components"); if (length(parts) == 3) { parts = c("1", "0", "4", "000") }; if (length(parts) != 4 || paste(parts[1:3], collapse = ".") != "1.0.4") { parts = c("1", "0", "4", "000") }; buildNumber = suppressWarnings(as.integer(parts[4])); if (is.na(buildNumber)) stop("Build component must be numeric"); parts[4] = sprintf("%03d", buildNumber + 1); lines[versionIndex] = paste0("Version: ", paste(parts, collapse = ".")); writeLines(lines, "DESCRIPTION")'
else
  echo "Skipping step 2/$total_steps: Bumping DESCRIPTION version"
fi

if should_run_step 3; then
  echo "[3/$total_steps] Running R string escape preflight"
  if [[ -f scripts/checkRStringEscapes.R ]]; then
    Rscript scripts/checkRStringEscapes.R
  else
    echo "Warning: scripts/checkRStringEscapes.R was not found; skipping escape preflight."
  fi
else
  echo "Skipping step 3/$total_steps: Running R string escape preflight"
fi

if should_run_step 4; then
  echo "[4/$total_steps] Running devtools::document()"
  Rscript -e 'options(warn = 2); devtools::document()'
else
  echo "Skipping step 4/$total_steps: Running devtools::document()"
fi

if should_run_step 5; then
  echo "[5/$total_steps] Running strict tests"
  Rscript -e 'options(warn = 2); devtools::test(stop_on_failure = TRUE, stop_on_warning = TRUE)'
else
  echo "Skipping step 5/$total_steps: Running strict tests"
fi

if should_run_step 6; then
  echo "[6/$total_steps] Running strict check"
  Rscript -e 'options(warn = 2); devtools::check(args = c("--no-manual", "--ignore-vignettes", "--no-tests"), error_on = "note")'
else
  echo "Skipping step 6/$total_steps: Running strict check"
fi

if should_run_step 7; then
  echo "[7/$total_steps] Committing changes"
  echo "Git status before staging:"
  git status --short

  paths_to_add=(
    DESCRIPTION
    NAMESPACE
    R
    tests
    man
    scripts
    inst
    dev
    vignettes
    .gitignore
    .Rbuildignore
    run_stage1.sh
  )

  for path_to_add in "${paths_to_add[@]}"; do
    if [[ -e "$path_to_add" ]]; then
      git add "$path_to_add"
    fi
  done

  shopt -s nullglob
  for file_to_add in *.md *.Rproj; do
    if [[ -f "$file_to_add" ]]; then
      git add "$file_to_add"
    fi
  done
  shopt -u nullglob

  echo "Git status after staging:"
  git status --short

  git commit -F - <<'COMMIT'
Add baseline test infrastructure for fitPS

- Add testthat infrastructure for core psData and psFit workflows
- Add baseline tests for data construction, CSV import, fitting, prediction, and log-likelihood extraction
- Ignore generated compressed archives in git and R package builds
- Add a reusable stage 1 runner that bumps 1.0.4.xxx build numbers before validation
- Not validated here unless this script reaches the strict test and check steps successfully
COMMIT
else
  echo "Skipping step 7/$total_steps: Committing changes"
fi

if should_run_step 8; then
  echo "[8/$total_steps] Creating stage archive"
  git archive --format=zip --output="$current_archive" HEAD
  shopt -s nullglob
  for old_archive in stage*_completed.zip; do
    if [[ "$old_archive" != "$current_archive" ]]; then
      rm -f "$old_archive"
    fi
  done
  shopt -u nullglob
else
  echo "Skipping step 8/$total_steps: Creating stage archive"
fi

if should_run_step 9; then
  echo "[9/$total_steps] Building package"
  Rscript -e 'options(warn = 2); pkg = devtools::build(); writeLines(pkg, ".fitps_built_pkg_path"); cat(pkg, "\n")'
else
  echo "Skipping step 9/$total_steps: Building package"
fi

if should_run_step 10; then
  echo "[10/$total_steps] Installing built package"
  if [[ ! -f "$built_pkg_file" ]]; then
    echo "Cannot install because $built_pkg_file is missing. Re-run from the build step or earlier."
    exit 1
  fi
  Rscript -e 'options(warn = 2); pkg = readLines(".fitps_built_pkg_path", warn = FALSE); if (length(pkg) != 1 || !nzchar(pkg)) stop("Built package path file did not contain exactly one path"); install.packages(pkg, repos = NULL, type = "source")'
else
  echo "Skipping step 10/$total_steps: Installing built package"
fi

echo "Stage 1 completed successfully."
