#!/usr/bin/env bash
set -euo pipefail

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  echo "Do not source this script. Run it with: bash ${BASH_SOURCE[0]}"
  return 1 2>/dev/null || exit 1
fi

show_usage() {
  cat <<'HELP'
Usage:
  scripts/runStage.sh STAGE [stage-runner-options]
  scripts/runStage.sh PATH_TO_BUNDLE [stage-runner-options]
  scripts/runStage.sh PATH_TO_RUNNER [stage-runner-options]

Examples:
  scripts/runStage.sh 3.1
  scripts/runStage.sh 3_1
  scripts/runStage.sh 3.1 -sn 7
  scripts/runStage.sh ~/Downloads/fitps_stage3_1_bundle.zip
  scripts/runStage.sh ~/Downloads/fitps_run_stage3_1.sh

Default bundle workflow:
  For Stage 3.1, the wrapper looks for:

    ~/Downloads/fitps_stage3_1_bundle.zip

  It extracts the bundle into:

    ~/Downloads/fitps_stage3_1_bundle/

  It then runs fitps_run_stage3_1.sh with --install-files. The bundled
  runner uses its sibling fitps_stage3_1_changes.zip.

Legacy workflow:
  If no matching bundle is found, the wrapper also supports older runner
  names such as run_stage3_1.sh and changes zips such as
  stage3_1_changes.zip.

Override behavior:
  Pass --no-install-files to run without installing files.
  Pass -cz PATH or --changes-zip PATH to choose a specific changes zip.
  Set FITPS_ROOT to override the default fitPS package root.

The wrapper runs the selected stage runner as a child process. It does not
source it.
HELP
}

absolute_path() {
  local input_path="$1"
  local input_dir
  local input_file

  input_dir="$(dirname "$input_path")"
  input_file="$(basename "$input_path")"
  printf '%s/%s\n' "$(cd "$input_dir" && pwd)" "$input_file"
}

resolve_downloads_dir() {
  case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*)
      if [[ -d "/c/Users/james/Downloads" ]]; then
        printf '%s\n' "/c/Users/james/Downloads"
      else
        printf '%s\n' "$HOME/Downloads"
      fi
      ;;
    *)
      printf '%s\n' "$HOME/Downloads"
      ;;
  esac
}

resolve_fitps_root() {
  if [[ -n "${FITPS_ROOT:-}" ]]; then
    printf '%s\n' "$FITPS_ROOT"
    return
  fi

  if [[ -f "$PWD/DESCRIPTION" ]] && grep -q '^Package:[[:space:]]*fitPS[[:space:]]*$' "$PWD/DESCRIPTION"; then
    printf '%s\n' "$PWD"
    return
  fi

  case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*)
      printf '%s\n' 'D:/Dropbox/Code/git/fitPS'
      ;;
    *)
      printf '%s\n' "$HOME/Dropbox/Code/git/fitPS"
      ;;
  esac
}

if [[ $# -lt 1 ]]; then
  show_usage
  exit 1
fi

stage_or_path="$1"
shift

if [[ "$stage_or_path" == "--help" || "$stage_or_path" == "-h" ]]; then
  show_usage
  exit 0
fi

downloads_dir="$(resolve_downloads_dir)"
fitps_root="$(resolve_fitps_root)"

if [[ ! -d "$fitps_root" ]]; then
  echo "Could not find the fitPS package root: $fitps_root"
  echo "Set FITPS_ROOT to the correct path and try again."
  exit 1
fi

if [[ ! -f "$fitps_root/DESCRIPTION" ]] || ! grep -q '^Package:[[:space:]]*fitPS[[:space:]]*$' "$fitps_root/DESCRIPTION"; then
  echo "The resolved directory is not a fitPS package root: $fitps_root"
  echo "Set FITPS_ROOT to the correct path and try again."
  exit 1
fi

stage_id=""
stage_file_id=""
runner_path=""
bundle_path=""
bundle_mode=0
cleanup_download_bundle=0
extraction_dir=""

cleanup_transient_bundle_files() {
  local exit_status=$?

  if [[ -n "$extraction_dir" && -d "$extraction_dir" ]]; then
    rm -rf "$extraction_dir"
  fi

  if [[ "$cleanup_download_bundle" -eq 1 && -n "$bundle_path" && -f "$bundle_path" ]]; then
    rm -f "$bundle_path"
  fi

  if [[ "$cleanup_download_bundle" -eq 1 && -n "$stage_file_id" ]]; then
    rm -f "$downloads_dir/fitps_stage${stage_file_id}_changes.zip"
  fi

  exit "$exit_status"
}

if [[ -f "$stage_or_path" ]]; then
  supplied_path="$(absolute_path "$stage_or_path")"

  case "$supplied_path" in
    *.zip)
      bundle_path="$supplied_path"
      bundle_file="$(basename "$bundle_path")"
      case "$bundle_file" in
        fitps_stage*_bundle.zip)
          stage_file_id="${bundle_file#fitps_stage}"
          stage_file_id="${stage_file_id%_bundle.zip}"
          stage_id="${stage_file_id//_/.}"
          bundle_mode=1
          ;;
        *)
          echo "The supplied zip is not a fitPS stage bundle: $bundle_file"
          exit 1
          ;;
      esac
      ;;
    *.sh)
      runner_path="$supplied_path"
      runner_file="$(basename "$runner_path")"
      case "$runner_file" in
        fitps_run_stage*.sh)
          stage_file_id="${runner_file#fitps_run_stage}"
          stage_file_id="${stage_file_id%.sh}"
          stage_id="${stage_file_id//_/.}"
          ;;
        run_stage*.sh)
          stage_file_id="${runner_file#run_stage}"
          stage_file_id="${stage_file_id%.sh}"
          stage_id="${stage_file_id//_/.}"
          ;;
        *)
          echo "The supplied shell script is not a recognised fitPS stage runner: $runner_file"
          exit 1
          ;;
      esac
      ;;
    *)
      echo "Expected a stage bundle zip or stage runner shell script: $supplied_path"
      exit 1
      ;;
  esac
else
  stage_id="${stage_or_path#stage}"
  stage_id="${stage_id//_/.}"
  stage_file_id="${stage_id//./_}"

  bundle_name="fitps_stage${stage_file_id}_bundle.zip"
  candidate_bundles=(
    "$PWD/$bundle_name"
    "$downloads_dir/$bundle_name"
  )

  for candidate in "${candidate_bundles[@]}"; do
    if [[ -f "$candidate" ]]; then
      bundle_path="$(absolute_path "$candidate")"
      bundle_mode=1
      if [[ "$bundle_path" == "$downloads_dir/"* ]]; then
        cleanup_download_bundle=1
      fi
      break
    fi
  done

  if [[ "$bundle_mode" -eq 0 ]]; then
    new_runner_name="fitps_run_stage${stage_file_id}.sh"
    legacy_runner_name="run_stage${stage_file_id}.sh"
    candidate_runners=(
      "$PWD/$new_runner_name"
      "$PWD/scripts/$new_runner_name"
      "$downloads_dir/$new_runner_name"
      "$PWD/$legacy_runner_name"
      "$PWD/scripts/$legacy_runner_name"
      "$downloads_dir/$legacy_runner_name"
    )

    for candidate in "${candidate_runners[@]}"; do
      if [[ -f "$candidate" ]]; then
        runner_path="$(absolute_path "$candidate")"
        break
      fi
    done

    if [[ -z "$runner_path" ]]; then
      echo "Could not find the Stage $stage_id bundle or runner."
      echo "Looked for bundle:"
      printf '  %s\n' "${candidate_bundles[@]}"
      echo "Looked for runners:"
      printf '  %s\n' "${candidate_runners[@]}"
      exit 1
    fi
  fi
fi

if [[ "$bundle_mode" -eq 1 ]]; then
  extraction_dir="$downloads_dir/fitps_stage${stage_file_id}_bundle"
  expected_runner="$extraction_dir/fitps_run_stage${stage_file_id}.sh"
  expected_changes_zip="$extraction_dir/fitps_stage${stage_file_id}_changes.zip"

  trap cleanup_transient_bundle_files EXIT INT TERM

  rm -rf "$extraction_dir"
  mkdir -p "$extraction_dir"
  unzip -q "$bundle_path" -d "$extraction_dir"

  if [[ ! -f "$expected_runner" ]]; then
    echo "The bundle does not contain the expected runner:"
    echo "  $expected_runner"
    exit 1
  fi

  if [[ ! -f "$expected_changes_zip" ]]; then
    echo "The bundle does not contain the expected change set:"
    echo "  $expected_changes_zip"
    exit 1
  fi

  chmod +x "$expected_runner"
  runner_path="$expected_runner"
fi

install_files=1
has_install_arg=0
has_changes_zip_arg=0
runner_args=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h)
      runner_args+=("$1")
      shift
      ;;
    --no-install-files)
      install_files=0
      shift
      ;;
    --install-files|-if)
      has_install_arg=1
      runner_args+=("$1")
      shift
      ;;
    --changes-zip|-cz)
      if [[ $# -lt 2 ]]; then
        echo "Missing value for $1"
        exit 1
      fi
      has_changes_zip_arg=1
      runner_args+=("$1" "$(absolute_path "$2")")
      shift 2
      ;;
    *)
      runner_args+=("$1")
      shift
      ;;
  esac
done

if [[ "$install_files" -eq 1 && "$has_install_arg" -eq 0 ]]; then
  runner_args=("--install-files" "${runner_args[@]}")
fi

if [[ "$bundle_mode" -eq 0 && "$install_files" -eq 1 && "$has_changes_zip_arg" -eq 0 && -n "$stage_file_id" ]]; then
  new_changes_zip_name="fitps_stage${stage_file_id}_changes.zip"
  legacy_changes_zip_name="stage${stage_file_id}_changes.zip"
  changes_zip_path=""
  candidate_zips=(
    "$(dirname "$runner_path")/$new_changes_zip_name"
    "$PWD/$new_changes_zip_name"
    "$PWD/scripts/$new_changes_zip_name"
    "$downloads_dir/$new_changes_zip_name"
    "$(dirname "$runner_path")/$legacy_changes_zip_name"
    "$PWD/$legacy_changes_zip_name"
    "$PWD/scripts/$legacy_changes_zip_name"
    "$downloads_dir/$legacy_changes_zip_name"
  )

  for candidate in "${candidate_zips[@]}"; do
    if [[ -f "$candidate" ]]; then
      changes_zip_path="$(absolute_path "$candidate")"
      break
    fi
  done

  if [[ -n "$changes_zip_path" ]]; then
    runner_args=("--changes-zip" "$changes_zip_path" "${runner_args[@]}")
  else
    echo "Warning: no matching changes zip was found."
    echo "The runner will use its own default change-set lookup."
  fi
fi

cd "$fitps_root"
bash "$runner_path" "${runner_args[@]}"
