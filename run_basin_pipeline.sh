#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/$(basename "$0")"
CONTROL_ENV_DIR="${CONTROL_ENV_DIR:-${SCRIPT_DIR}/.controller_env}"
ORIGINAL_ARGS=("$@")

PROCESS_ORDER=(
  BIOCHEM_MERGE
  BIOCHEM_DENSITY
  BIOCHEM_STRAT_METRICS
  BIOCHEM_CUSTOM_CLEAN
  BIOCHEM_EIGENVECTORS
  BIOCHEM_SELECTK
  BIOCHEM_GMM
  BIOCHEM_O2_SOFT
  BIOCHEM_HYBRID
  BIOCHEM_COMPARE
  BIOCHEM_SPLIT_O2_BY_GMM
  BIOCHEM_STRAT_ANOMALY
  BIOCHEM_STATE_TRANSITIONS
  BIOCHEM_SUCCESSION_GRAPH
  BIOCHEM_FEATURE_ASSOC
  BIOCHEM_EOF_PIPELINE
  BIOCHEM_EOF_STATE_CLUSTER
  BIOCHEM_EOF_MODE_PLOTS
  BIOCHEM_WITHIN_GMM_HDBSCAN
  MASTER_SUMMARY
)

usage() {
  cat <<'EOF'
Usage:
  run_basin_pipeline.sh [CONFIG_FILE] [--rerun-from PROCESS_NAME] [--resume-run RUN_NAME] [--resume-policy POLICY] [--no-resume] [--list-stages] [-- NEXTFLOW_ARGS...]

Options:
  --rerun-from PROCESS_NAME  Force rerun starting at PROCESS_NAME and all later processes in controller order.
  --resume-run RUN_NAME      Resume from a specific Nextflow run name/id.
  --resume-policy POLICY     Resume baseline selection when --resume-run is not set.
                             Allowed: last-with-tasks (default), latest
  --no-resume                Disable resume for this run.
  --list-stages              Print known process names for --rerun-from and exit.
  --help, -h                 Show this help.

Examples:
  ./run_basin_pipeline.sh basin_pipeline_nextflow.yml
  ./run_basin_pipeline.sh examples/si.local.yml --rerun-from BIOCHEM_EIGENVECTORS
  ./run_basin_pipeline.sh my_run.yml -- -with-report report.html -with-trace trace.tsv
EOF
}

if [[ -z "${IN_CONTROLLER_ENV:-}" ]]; then
  if ! command -v mamba >/dev/null 2>&1; then
    echo "mamba is required to bootstrap the controller environment." >&2
    exit 1
  fi
  if [[ ! -d "$CONTROL_ENV_DIR" ]]; then
    echo "[controller] Creating mamba env at $CONTROL_ENV_DIR"
    mamba env create --yes --prefix "$CONTROL_ENV_DIR" --file "${SCRIPT_DIR}/processes/shared_envs/controller.yml"
  elif [[ ! -x "$CONTROL_ENV_DIR/bin/nextflow" || ! -x "$CONTROL_ENV_DIR/bin/yq" || ! -x "$CONTROL_ENV_DIR/bin/python" ]]; then
    echo "[controller] Repairing incomplete controller env at $CONTROL_ENV_DIR"
    mamba env update --prune --prefix "$CONTROL_ENV_DIR" --file "${SCRIPT_DIR}/processes/shared_envs/controller.yml"
  fi
  exec env IN_CONTROLLER_ENV=1 CONTROL_ENV_DIR="$CONTROL_ENV_DIR" conda run --no-capture-output -p "$CONTROL_ENV_DIR" "$SCRIPT_PATH" "$@"
fi

CONFIG_FILE="basin_pipeline_nextflow.yml"
CONFIG_SET=0
LIST_STAGES=0
RERUN_FROM=""
NEXTFLOW_ARGS=()
RESUME_ENABLED=1
RESUME_POLICY="last-with-tasks"
RESUME_RUN=""
RERUN_CONFIG_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --list-stages)
      LIST_STAGES=1
      shift
      ;;
    --rerun-from)
      if [[ $# -lt 2 ]]; then
        echo "--rerun-from requires a process name" >&2
        exit 1
      fi
      RERUN_FROM="$2"
      shift 2
      ;;
    --rerun-from=*)
      RERUN_FROM="${1#*=}"
      shift
      ;;
    --resume-run)
      if [[ $# -lt 2 ]]; then
        echo "--resume-run requires a run name or id" >&2
        exit 1
      fi
      RESUME_RUN="$2"
      shift 2
      ;;
    --resume-run=*)
      RESUME_RUN="${1#*=}"
      shift
      ;;
    --resume-policy)
      if [[ $# -lt 2 ]]; then
        echo "--resume-policy requires a value: last-with-tasks|latest" >&2
        exit 1
      fi
      RESUME_POLICY="$2"
      shift 2
      ;;
    --resume-policy=*)
      RESUME_POLICY="${1#*=}"
      shift
      ;;
    --no-resume)
      RESUME_ENABLED=0
      shift
      ;;
    --)
      shift
      NEXTFLOW_ARGS+=("$@")
      break
      ;;
    -*)
      NEXTFLOW_ARGS+=("$1")
      shift
      ;;
    *)
      if [[ $CONFIG_SET -eq 0 ]]; then
        CONFIG_FILE="$1"
        CONFIG_SET=1
      else
        NEXTFLOW_ARGS+=("$1")
      fi
      shift
      ;;
  esac
done

if [[ "$LIST_STAGES" -eq 1 ]]; then
  printf '%s\n' "${PROCESS_ORDER[@]}"
  exit 0
fi

if [[ "$RESUME_POLICY" != "last-with-tasks" && "$RESUME_POLICY" != "latest" ]]; then
  echo "Invalid --resume-policy '${RESUME_POLICY}'. Allowed: last-with-tasks, latest" >&2
  exit 1
fi

if [[ "$RESUME_ENABLED" -eq 0 && -n "$RESUME_RUN" ]]; then
  echo "--no-resume cannot be combined with --resume-run" >&2
  exit 1
fi

sanitize_nextflow_args() {
  local sanitized=()
  local skip_next=0
  local i arg next_arg
  for ((i=0; i<${#NEXTFLOW_ARGS[@]}; i++)); do
    if [[ $skip_next -eq 1 ]]; then
      skip_next=0
      continue
    fi
    arg="${NEXTFLOW_ARGS[$i]}"
    case "$arg" in
      -resume)
        if (( i + 1 < ${#NEXTFLOW_ARGS[@]} )); then
          next_arg="${NEXTFLOW_ARGS[$((i + 1))]}"
          if [[ "$next_arg" != -* ]]; then
            skip_next=1
          fi
        fi
        echo "[controller] Ignoring passthrough '-resume' argument; use --resume-run/--resume-policy/--no-resume." >&2
        ;;
      -resume=*)
        echo "[controller] Ignoring passthrough '-resume=...' argument; use --resume-run/--resume-policy/--no-resume." >&2
        ;;
      *)
        sanitized+=("$arg")
        ;;
    esac
  done
  NEXTFLOW_ARGS=("${sanitized[@]}")
}

has_tasks_for_run() {
  local run_name="$1"
  local rows rc
  set +e
  rows="$(nextflow log "$run_name" -f 'process,workdir,status' 2>/dev/null)"
  rc=$?
  set -e
  [[ $rc -eq 0 ]] || return 1
  while IFS=$'\t' read -r process_name workdir _status; do
    [[ -n "$process_name" ]] || continue
    [[ "$process_name" == "process" && "$workdir" == "workdir" ]] && continue
    return 0
  done <<< "$rows"
  return 1
}

select_baseline_run() {
  local selected=""
  local latest=""
  local runs=()
  local i run_name

  if [[ -n "$RESUME_RUN" ]]; then
    if ! nextflow log "$RESUME_RUN" >/dev/null 2>&1; then
      echo "Specified --resume-run not found in Nextflow history: ${RESUME_RUN}" >&2
      exit 1
    fi
    echo "$RESUME_RUN"
    return 0
  fi

  mapfile -t runs < <(nextflow log -q 2>/dev/null | sed '/^[[:space:]]*$/d')
  if [[ ${#runs[@]} -eq 0 ]]; then
    echo ""
    return 0
  fi
  latest="${runs[$(( ${#runs[@]} - 1 ))]}"

  if [[ "$RESUME_POLICY" == "latest" ]]; then
    echo "$latest"
    return 0
  fi

  for ((i=${#runs[@]}-1; i>=0; i--)); do
    run_name="${runs[$i]}"
    if has_tasks_for_run "$run_name"; then
      selected="$run_name"
      break
    fi
  done

  echo "${selected:-$latest}"
}

sanitize_nextflow_args

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Config file not found: $CONFIG_FILE" >&2
  exit 1
fi

CONFIG_FILE="$(realpath "$CONFIG_FILE")"
CONFIG_ROOT="$(dirname "$CONFIG_FILE")"

resolve_config_path() {
  local value="$1"
  if [[ "$value" = /* ]]; then
    realpath -m "$value"
  else
    realpath -m "${CONFIG_ROOT}/${value}"
  fi
}

WORK_DIR=$(yq -r '.paths.work_dir // empty' "$CONFIG_FILE")
OUTPUT_DIR=$(yq -r '.paths.output_dir // empty' "$CONFIG_FILE")
CONDA_CACHE_DIR=$(yq -r '.paths.conda_cache_dir // empty' "$CONFIG_FILE")
RUNTIME_DIR=$(yq -r '.paths.runtime_dir // empty' "$CONFIG_FILE")
KEEP_RUNTIME_DIR=$(yq -r '.paths.keep_runtime_dir // true' "$CONFIG_FILE")

if [[ -z "$OUTPUT_DIR" ]]; then
  echo "paths.output_dir must be set in $CONFIG_FILE" >&2
  exit 1
fi

OUTPUT_DIR="$(resolve_config_path "$OUTPUT_DIR")"
if [[ -z "$RUNTIME_DIR" || "$RUNTIME_DIR" == "null" ]]; then
  RUNTIME_DIR="${OUTPUT_DIR}/.basin"
else
  RUNTIME_DIR="$(resolve_config_path "$RUNTIME_DIR")"
fi
if [[ -z "$WORK_DIR" || "$WORK_DIR" == "null" ]]; then
  WORK_DIR="${RUNTIME_DIR}/nf_work"
else
  WORK_DIR="$(resolve_config_path "$WORK_DIR")"
fi
if [[ -z "$CONDA_CACHE_DIR" || "$CONDA_CACHE_DIR" == "null" ]]; then
  CONDA_CACHE_DIR="${RUNTIME_DIR}/conda_cache"
else
  CONDA_CACHE_DIR="$(resolve_config_path "$CONDA_CACHE_DIR")"
fi
PUBLICATION_STAGING_DIR="${RUNTIME_DIR}/publication_staging"

case "${KEEP_RUNTIME_DIR,,}" in
  true|false) ;;
  *) echo "paths.keep_runtime_dir must be true or false: ${KEEP_RUNTIME_DIR}" >&2; exit 1 ;;
esac

mkdir -p "$WORK_DIR" "$CONDA_CACHE_DIR" "${CONDA_CACHE_DIR}/pkgs"
mkdir -p "${OUTPUT_DIR}/modules" "${OUTPUT_DIR}/intermediates" \
  "${OUTPUT_DIR}/references" "${OUTPUT_DIR}/summary" "${OUTPUT_DIR}/logs"

exec 7>"${CONDA_CACHE_DIR}/.basin-run.lock"
if ! flock -n 7; then
  echo "Another BASIN run is using Conda cache: ${CONDA_CACHE_DIR}" >&2
  echo "Wait for it to finish or configure another paths.conda_cache_dir." >&2
  exit 1
fi
stale_env_locks=()
while IFS= read -r -d '' stale_lock; do
  stale_env_locks+=("$stale_lock")
done < <(find "$CONDA_CACHE_DIR" -maxdepth 1 -type f -name '.env-*.lock' -print0)
if (( ${#stale_env_locks[@]} > 0 )); then
  rm -f "${stale_env_locks[@]}"
  echo "[controller] Removed ${#stale_env_locks[@]} stale Nextflow Conda environment lock marker(s)."
fi
if [[ "$RESUME_ENABLED" -eq 0 && -d "$PUBLICATION_STAGING_DIR" ]]; then
  staging_real="$(realpath -m "$PUBLICATION_STAGING_DIR")"
  runtime_real="$(realpath -m "$RUNTIME_DIR")"
  if [[ "$staging_real" == "$runtime_real"/* ]]; then
    rm -rf "$staging_real"
  else
    echo "Refusing to clear publication staging outside runtime directory: ${staging_real}" >&2
    exit 1
  fi
fi
mkdir -p "${PUBLICATION_STAGING_DIR}/logs"

has_nextflow_arg() {
  local expected="$1" arg
  for arg in "${NEXTFLOW_ARGS[@]}"; do
    [[ "$arg" == "$expected" || "$arg" == "${expected}="* ]] && return 0
  done
  return 1
}

DEFAULT_REPORT_ARGS=()
has_nextflow_arg -with-report || DEFAULT_REPORT_ARGS+=(-with-report "${PUBLICATION_STAGING_DIR}/logs/nextflow_report.html")
has_nextflow_arg -with-timeline || DEFAULT_REPORT_ARGS+=(-with-timeline "${PUBLICATION_STAGING_DIR}/logs/nextflow_timeline.html")
has_nextflow_arg -with-trace || DEFAULT_REPORT_ARGS+=(-with-trace "${PUBLICATION_STAGING_DIR}/logs/nextflow_trace.tsv")
has_nextflow_arg -with-dag || DEFAULT_REPORT_ARGS+=(-with-dag "${PUBLICATION_STAGING_DIR}/logs/nextflow_dag.html")
rm -f "${PUBLICATION_STAGING_DIR}/logs/nextflow_report.html" \
  "${PUBLICATION_STAGING_DIR}/logs/nextflow_timeline.html" \
  "${PUBLICATION_STAGING_DIR}/logs/nextflow_trace.tsv" \
  "${PUBLICATION_STAGING_DIR}/logs/nextflow_dag.html"
printf '%q ' "$SCRIPT_PATH" "${ORIGINAL_ARGS[@]}" > "${PUBLICATION_STAGING_DIR}/logs/launch_command.txt"
printf '\n' >> "${PUBLICATION_STAGING_DIR}/logs/launch_command.txt"
nextflow -version > "${PUBLICATION_STAGING_DIR}/logs/nextflow_version.txt" 2>&1

export NXF_WORK="$WORK_DIR"
export NXF_CONDA_CACHEDIR="$CONDA_CACHE_DIR"
export CONDA_PKGS_DIRS="${CONDA_CACHE_DIR}/pkgs"
export NXF_SYNTAX_PARSER="${NXF_SYNTAX_PARSER:-v1}"

BASELINE_RUN="$(select_baseline_run)"
if [[ -n "$BASELINE_RUN" ]]; then
  echo "[controller] Baseline run for cache/history lookup: ${BASELINE_RUN} (policy: ${RESUME_POLICY})"
else
  echo "[controller] No prior Nextflow run history found."
fi

if [[ -n "$RERUN_FROM" ]]; then
  RERUN_FROM_UPPER="$(printf '%s' "$RERUN_FROM" | tr '[:lower:]' '[:upper:]')"
  start_idx=-1
  for i in "${!PROCESS_ORDER[@]}"; do
    if [[ "${PROCESS_ORDER[$i]}" == "$RERUN_FROM_UPPER" ]]; then
      start_idx=$i
      break
    fi
  done

  if [[ $start_idx -lt 0 ]]; then
    echo "Unknown process for --rerun-from: $RERUN_FROM" >&2
    echo "Use --list-stages to see valid process names." >&2
    exit 1
  fi

  declare -A RERUN_STAGE_SET=()
  for ((i=start_idx; i<${#PROCESS_ORDER[@]}; i++)); do
    RERUN_STAGE_SET["${PROCESS_ORDER[$i]}"]=1
  done

  RERUN_CONFIG_FILE="$(mktemp /tmp/basin-rerun-cache-XXXXXX.config)"
  {
    echo "process {"
    for stage_name in "${PROCESS_ORDER[@]}"; do
      [[ -n "${RERUN_STAGE_SET[$stage_name]:-}" ]] || continue
      printf "  withName: /(^|.*:)%s\$/ { cache = false }\n" "$stage_name"
    done
    echo "}"
  } > "$RERUN_CONFIG_FILE"
  echo "[controller] Generated temporary Nextflow config to disable cache from ${RERUN_FROM_UPPER} onward: ${RERUN_CONFIG_FILE}"
fi

NEXTFLOW_RERUN_ARGS=()
if [[ -n "$RERUN_CONFIG_FILE" ]]; then
  NEXTFLOW_RERUN_ARGS=(-c "$RERUN_CONFIG_FILE")
fi

if [[ "$RESUME_ENABLED" -eq 1 ]]; then
  if [[ -n "$BASELINE_RUN" ]]; then
    RESUME_ARGS=(-resume "$BASELINE_RUN")
  else
    RESUME_ARGS=(-resume)
  fi
else
  RESUME_ARGS=()
  echo "[controller] Resume disabled for this run (--no-resume)."
fi

BASIN_PIPELINE_CONFIG="$CONFIG_FILE" \
nextflow run "${SCRIPT_DIR}/basin_pipeline.nf" \
  --params-file "$CONFIG_FILE" \
  --pipeline_config "$CONFIG_FILE" \
  "${RESUME_ARGS[@]}" \
  "${NEXTFLOW_RERUN_ARGS[@]}" \
  -w "$NXF_WORK" \
  "${DEFAULT_REPORT_ARGS[@]}" \
  "${NEXTFLOW_ARGS[@]}" 2>&1 | tee "${PUBLICATION_STAGING_DIR}/logs/controller.log"

COMPLETED_RUN="$(nextflow log -q 2>/dev/null | tail -n 1 || true)"
if [[ -n "$COMPLETED_RUN" ]]; then
  nextflow log "$COMPLETED_RUN" -f 'process,hash,workdir,status,exit,duration' \
    > "${PUBLICATION_STAGING_DIR}/logs/task_execution.tsv" || true
fi
if [[ -f "${SCRIPT_DIR}/.nextflow.log" ]]; then
  cp "${SCRIPT_DIR}/.nextflow.log" "${PUBLICATION_STAGING_DIR}/logs/nextflow.log"
fi

python "${SCRIPT_DIR}/processes/output_layout/organize_outputs.py" \
  --staging-dir "$PUBLICATION_STAGING_DIR" \
  --output-dir "$OUTPUT_DIR" \
  --config "$CONFIG_FILE"

if [[ "${KEEP_RUNTIME_DIR,,}" == "true" ]]; then
  echo "[controller] Retaining Nextflow resume state: ${RUNTIME_DIR}"
else
  runtime_real="$(realpath -m "$RUNTIME_DIR")"
  output_real="$(realpath -m "$OUTPUT_DIR")"
  script_real="$(realpath -m "$SCRIPT_DIR")"
  if [[ "$runtime_real" == "/" || "$runtime_real" == "$output_real" ||
        "$output_real" == "$runtime_real"/* || "$script_real" == "$runtime_real" ||
        "$script_real" == "$runtime_real"/* ]]; then
    echo "Refusing to remove unsafe runtime directory: ${runtime_real}" >&2
    exit 1
  fi
  rm -rf "$runtime_real"
  echo "[controller] Removed runtime directory after successful publication: ${runtime_real}"
fi
