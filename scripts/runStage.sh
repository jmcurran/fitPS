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

Examples:
  scripts/runStage.sh 1
  scripts/runStage.sh 1 -sn 2

Behavior:
  The wrapper detects the platform, locates the matching runner and change-set
  zip in Downloads, and runs:

    <downloads-prefix>/run_stageSTAGE.sh --install-files --changes-zip <downloads-prefix>/fitPS_stageSTAGE_changes.zip [stage-runner-options]

Platform defaults:
  Windows Git Bash/MSYS/MINGW/Cygwin: /c/Users/james/Downloads
  macOS/Linux:                       ~/Downloads

The wrapper runs the selected stage runner as a child process. It does not source it.
HELP
}

if [[ $# -lt 1 ]]; then
  show_usage
  exit 1
fi

stage_id="$1"
shift

if [[ "$stage_id" == "--help" || "$stage_id" == "-h" ]]; then
  show_usage
  exit 0
fi

stage_id="${stage_id#stage}"
stage_id="${stage_id//./_}"

case "$(uname -s)" in
  MINGW*|MSYS*|CYGWIN*)
    downloads_prefix="/c/Users/james/Downloads"
    ;;
  Darwin*)
    downloads_prefix="$HOME/Downloads"
    ;;
  *)
    downloads_prefix="$HOME/Downloads"
    ;;
esac

runner_path="$downloads_prefix/run_stage${stage_id}.sh"
changes_zip_path="$downloads_prefix/fitPS_stage${stage_id}_changes.zip"

if [[ ! -f "$runner_path" ]]; then
  echo "Could not find stage runner:"
  echo "  $runner_path"
  exit 1
fi

if [[ ! -f "$changes_zip_path" ]]; then
  echo "Could not find stage change-set zip:"
  echo "  $changes_zip_path"
  exit 1
fi

bash "$runner_path" --install-files --changes-zip "$changes_zip_path" "$@"
