#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd -P)"

usage() {
  cat <<'EOF'
Usage: pipeline/run_pipeline.sh [options]

Run GO-CAM release pipeline steps in a local development workspace.

Options:
  --models-dir PATH  Directory containing raw Minerva model files.
  --output-root-dir PATH
                     Root directory for pipeline outputs.
  --steps NAMES      Comma-separated step names (default: all steps).
                     Steps: convert, filter, index, index-files,
                     browser-search, stats, summary.
  --clean            Remove the entire output root directory before running.
  --dry-run          Print planned cleanup and commands without changing files.
  -h, --help         Show this help.
EOF
}

models_dir_from_environment="${GOCAM_PIPELINE_MODELS_DIR:-}"
output_root_dir_from_environment="${GOCAM_PIPELINE_OUTPUT_ROOT_DIR:-}"

if [[ -f "$REPO_ROOT/.env" ]]; then
  set -a
  # shellcheck disable=SC1091
  source "$REPO_ROOT/.env"
  set +a
fi

models_dir="${models_dir_from_environment:-${GOCAM_PIPELINE_MODELS_DIR:-}}"
output_root_dir="${output_root_dir_from_environment:-${GOCAM_PIPELINE_OUTPUT_ROOT_DIR:-/tmp/gocam-work}}"
steps="all"
clean=false
dry_run=false

while (($#)); do
  case "$1" in
    --models-dir)
      models_dir="${2:?--models-dir requires a path}"
      shift 2
      ;;
    --output-root-dir)
      output_root_dir="${2:?--output-root-dir requires a path}"
      shift 2
      ;;
    --steps)
      steps="${2:?--steps requires a comma-separated list}"
      shift 2
      ;;
    --clean)
      clean=true
      shift
      ;;
    --dry-run)
      dry_run=true
      shift
      ;;
    -h | --help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$models_dir" ]]; then
  echo "GOCAM_PIPELINE_MODELS_DIR or --models-dir is required." >&2
  exit 2
fi

if [[ ! -d "$models_dir" ]]; then
  echo "Models directory does not exist: $models_dir" >&2
  exit 2
fi

canonical_path() {
  python3 -c 'import os, sys; print(os.path.realpath(os.path.abspath(os.path.expanduser(sys.argv[1]))))' "$1"
}

models_dir="$(canonical_path "$models_dir")"
output_root_dir="$(canonical_path "$output_root_dir")"
home_dir="$(canonical_path "${HOME:-/}")"

path_is_within() {
  [[ "$1" == "$2" || "$1" == "$2/"* ]]
}

if [[ "$clean" == true ]] \
  && { [[ "$output_root_dir" == "/" ]] \
    || path_is_within "$REPO_ROOT" "$output_root_dir" \
    || path_is_within "$home_dir" "$output_root_dir" \
    || path_is_within "$models_dir" "$output_root_dir"; }; then
  echo "Refusing to clean unsafe output root directory: $output_root_dir" >&2
  exit 2
fi

logs_dir="$output_root_dir/logs"
known_steps=(convert filter index index-files browser-search stats summary)

if [[ "$steps" == "all" ]]; then
  selected_steps="$(IFS=,; echo "${known_steps[*]}")"
else
  selected_steps="$steps"
fi

IFS=',' read -r -a requested_steps <<<"$selected_steps"
for requested_step in "${requested_steps[@]}"; do
  known=false
  for known_step in "${known_steps[@]}"; do
    if [[ "$requested_step" == "$known_step" ]]; then
      known=true
      break
    fi
  done
  if [[ "$known" == false ]]; then
    echo "Unknown step: $requested_step" >&2
    exit 2
  fi
done

step_selected() {
  [[ ",$selected_steps," == *",$1,"* ]]
}

require_input_dir() {
  local input_dir="$1"
  local consumer="$2"
  local producer="$3"
  if step_selected "$producer"; then
    return
  fi
  if [[ "$clean" == true || ! -d "$input_dir" ]]; then
    echo "Input directory for step '$2' does not exist: $1" >&2
    exit 2
  fi
}

if step_selected filter; then
  require_input_dir "$output_root_dir/01-gocam-models" filter convert
fi
if step_selected index; then
  require_input_dir "$output_root_dir/02-true-gocams" index filter
fi
if step_selected index-files; then
  require_input_dir "$output_root_dir/03-indexed-true-gocams" index-files index
fi
if step_selected browser-search; then
  require_input_dir "$output_root_dir/03-indexed-true-gocams" browser-search index
fi
if step_selected stats; then
  require_input_dir "$output_root_dir/03-indexed-true-gocams" stats index
fi
if step_selected summary \
  && ! step_selected convert \
  && ! step_selected filter \
  && ! step_selected index \
  && ! step_selected index-files \
  && ! step_selected browser-search; then
  if [[ "$clean" == true || ! -d "$logs_dir" ]]; then
    echo "Logs directory for step 'summary' does not exist: $logs_dir" >&2
    exit 2
  fi
  if ! compgen -G "$logs_dir/*.jsonl" >/dev/null; then
    echo "No JSONL reports are available for step 'summary': $logs_dir" >&2
    exit 2
  fi
fi

printf 'Models directory: %s\n' "$models_dir"
printf 'Output root directory: %s\n' "$output_root_dir"

if [[ "$clean" == true ]]; then
  if [[ "$dry_run" == true ]]; then
    printf 'Would remove: %s\n' "$output_root_dir"
  else
    rm -rf -- "$output_root_dir"
  fi
fi

prepare_output_dir() {
  if [[ "$dry_run" == true ]]; then
    printf 'Would replace output directory: %s\n' "$1"
  else
    rm -rf -- "$1"
    mkdir -p "$1"
  fi
}

prepare_output_file() {
  if [[ "$dry_run" == true ]]; then
    printf 'Would replace output file: %s\n' "$1"
  else
    rm -f -- "$1"
    mkdir -p "$(dirname "$1")"
  fi
}

run_command() {
  printf '+'
  printf ' %q' "$@"
  printf '\n'
  if [[ "$dry_run" == false ]]; then
    "$@"
  fi
}

if [[ "$dry_run" == false ]]; then
  mkdir -p "$logs_dir"
fi

if step_selected convert; then
  output_dir="$output_root_dir/01-gocam-models"
  report_file="$logs_dir/01-convert.jsonl"
  prepare_output_dir "$output_dir"
  prepare_output_file "$report_file"
  run_command uv run --project "$REPO_ROOT" python "$SCRIPT_DIR/convert_minerva_models_to_gocam_models.py" \
    --input-dir "$models_dir" \
    --output-dir "$output_dir" \
    --report-file "$report_file"
fi

if step_selected filter; then
  output_dir="$output_root_dir/02-true-gocams"
  pseudo_output_dir="$output_root_dir/02-pseudo-gocams"
  report_file="$logs_dir/02-filter.jsonl"
  prepare_output_dir "$output_dir"
  prepare_output_dir "$pseudo_output_dir"
  prepare_output_file "$report_file"
  run_command uv run --project "$REPO_ROOT" python "$SCRIPT_DIR/filter_true_gocam_models.py" \
    --input-dir "$output_root_dir/01-gocam-models" \
    --output-dir "$output_dir" \
    --pseudo-gocam-output-dir "$pseudo_output_dir" \
    --report-file "$report_file"
fi

if step_selected index; then
  output_dir="$output_root_dir/03-indexed-true-gocams"
  report_file="$logs_dir/03-index.jsonl"
  prepare_output_dir "$output_dir"
  prepare_output_file "$report_file"
  run_command uv run --project "$REPO_ROOT" python "$SCRIPT_DIR/add_query_index_to_models.py" \
    --input-dir "$output_root_dir/02-true-gocams" \
    --output-dir "$output_dir" \
    --report-file "$report_file"
fi

if step_selected index-files; then
  output_dir="$output_root_dir/04-index-files"
  report_file="$logs_dir/04-index-files.jsonl"
  prepare_output_dir "$output_dir"
  prepare_output_file "$report_file"
  run_command uv run --project "$REPO_ROOT" python "$SCRIPT_DIR/generate_index_files.py" \
    --input-dir "$output_root_dir/03-indexed-true-gocams" \
    --output-dir "$output_dir" \
    --report-file "$report_file"
fi

if step_selected browser-search; then
  output_dir="$output_root_dir/05-browser-search-docs"
  report_file="$logs_dir/05-browser-search.jsonl"
  prepare_output_dir "$output_dir"
  prepare_output_file "$report_file"
  run_command uv run --project "$REPO_ROOT" python "$SCRIPT_DIR/generate_go_cam_browser_search_docs.py" \
    --input-dir "$output_root_dir/03-indexed-true-gocams" \
    --output "$output_dir/go-cam-browser-search-docs.json" \
    --report-file "$report_file"
fi

if step_selected stats; then
  output_dir="$output_root_dir/06-stats"
  prepare_output_dir "$output_dir"
  run_command uv run --project "$REPO_ROOT" python "$SCRIPT_DIR/output_stats_for_gocam_models.py" \
    --input-dir "$output_root_dir/03-indexed-true-gocams" \
    --output-dir "$output_dir"
fi

if step_selected summary; then
  excel_output_file="$logs_dir/summary.xlsx"
  html_output_dir="$logs_dir/summary-html"
  prepare_output_file "$excel_output_file"
  prepare_output_dir "$html_output_dir"
  run_command uv run --project "$REPO_ROOT" python "$SCRIPT_DIR/generate_log_summary.py" \
    --logs-dir "$logs_dir" \
    --excel-output "$excel_output_file" \
    --html-output-dir "$html_output_dir"
fi
