#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 && $# -ne 10 ]]; then
  echo "Usage: $0 INPUT GENOME_ID TAXONOMY ALIGNER MMSEQS2_BIN THREADS OUTDIR [RETRY_ENABLED RETRY_MIN_GROWTH RETRY_SOLVER]" >&2
  exit 2
fi

input=$1
genome_id=$2
taxonomy=$3
aligner=$4
mmseqs2_bin=$5
threads=$6
outdir=$7
retry_enabled=${8:-true}
retry_min_growth=${9:-0.001}
retry_solver=${10:-glpk}

case "$taxonomy" in Bacteria|Archaea|auto) ;; *) echo "Invalid taxonomy: $taxonomy" >&2; exit 2 ;; esac
case "$aligner" in blast|diamond|mmseqs2) ;; *) echo "Invalid aligner: $aligner" >&2; exit 2 ;; esac
if [[ "$aligner" == "mmseqs2" ]]; then
  [[ -x "$mmseqs2_bin" ]] || { echo "MMseqs2 executable is missing: $mmseqs2_bin" >&2; exit 2; }
  export PATH="$(dirname "$mmseqs2_bin"):$PATH"
fi
case "$threads" in ''|*[!0-9]*) echo "THREADS must be a positive integer" >&2; exit 2 ;; esac
(( threads > 0 )) || { echo "THREADS must be positive" >&2; exit 2; }
case "$retry_enabled" in true|false) ;; *) echo "RETRY_ENABLED must be true or false" >&2; exit 2 ;; esac
case "$retry_solver" in glpk|cplex) ;; *) echo "RETRY_SOLVER must be glpk or cplex" >&2; exit 2 ;; esac
awk -v value="$retry_min_growth" 'BEGIN { exit !(value >= 0.001 && value < 0.01) }' || {
  echo "RETRY_MIN_GROWTH must be at least 0.001 and below the gapseq default of 0.01" >&2
  exit 2
}
[[ -s "$input" ]] || { echo "Genome input is missing or empty: $input" >&2; exit 2; }
mkdir -p "$outdir"
default_log="$outdir/gapseq_doall.log"
set +e
gapseq doall \
  -A "$aligner" \
  -K "$threads" \
  -t "$taxonomy" \
  -f "$outdir" \
  "$input" 2>&1 | tee "$default_log"
default_status=${PIPESTATUS[0]}
set -e

retry_attempted=false
retry_status=""
reconstruction_mode="gapseq_doall_default"
effective_min_growth=0.01
effective_gapfill_solver=auto
if (( default_status != 0 )); then
  if [[ "$retry_enabled" != true ]] || ! grep -Eq 'Final model cannot grow \(status=[0-9]+\)' "$default_log"; then
    echo "gapseq doall failed without an eligible numerical final-model failure; see $default_log" >&2
    exit "$default_status"
  fi
  mapfile -t draft_models < <(
    find "$outdir" -maxdepth 1 -type f -name '*-draft.RDS' -print
  )
  mapfile -t predicted_media < <(
    find "$outdir" -maxdepth 1 -type f -name '*-medium.csv' -print
  )
  if (( ${#draft_models[@]} != 1 || ${#predicted_media[@]} != 1 )); then
    printf 'Numerical retry requires one draft RDS and one predicted medium; found %d and %d\n' \
      "${#draft_models[@]}" "${#predicted_media[@]}" >&2
    exit "$default_status"
  fi
  retry_attempted=true
  reconstruction_mode="gapseq_gapfill_reduced_minimum_growth_retry"
  effective_min_growth=$retry_min_growth
  effective_gapfill_solver=$retry_solver
  retry_log="$outdir/gapseq_gapfill_retry.log"
  set +e
  gapseq fill \
    -m "${draft_models[0]}" \
    -n "${predicted_media[0]}" \
    -f "$outdir" \
    -k "$retry_min_growth" \
    -z "$retry_solver" 2>&1 | tee "$retry_log"
  retry_status=${PIPESTATUS[0]}
  set -e
  {
    printf 'genome_id\tretry_attempted\tdefault_exit_status\tretry_exit_status\tretry_minimum_growth\tretry_solver\tdefault_log\tretry_log\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$genome_id" "$retry_attempted" "$default_status" "$retry_status" \
      "$retry_min_growth" "$retry_solver" "$(basename "$default_log")" "$(basename "$retry_log")"
  } > "$outdir/reconstruction_retry_audit.tsv"
  if (( retry_status != 0 )); then
    echo "Reduced-minimum-growth gapseq retry also failed; see $retry_log" >&2
    exit "$retry_status"
  fi
else
  {
    printf 'genome_id\tretry_attempted\tdefault_exit_status\tretry_exit_status\tretry_minimum_growth\tretry_solver\tdefault_log\tretry_log\n'
    printf '%s\tfalse\t0\t\t\t\t%s\t\n' "$genome_id" "$(basename "$default_log")"
  } > "$outdir/reconstruction_retry_audit.tsv"
fi

mapfile -t models < <(
  find "$outdir" -maxdepth 1 -type f -name '*.xml' ! -iname '*-draft.xml' -print
)
if (( ${#models[@]} != 1 )); then
  printf 'Expected exactly one non-draft final SBML XML, found %d in %s\n' \
    "${#models[@]}" "$outdir" >&2
  if (( ${#models[@]} > 0 )); then
    printf 'Final-model candidates:\n' >&2
    printf '  %s\n' "${models[@]}" >&2
  fi
  exit 1
fi
if [[ "$(readlink -f "${models[0]}")" != "$(readlink -m "$outdir/${genome_id}.xml")" ]]; then
  cp "${models[0]}" "$outdir/${genome_id}.xml"
fi

{
  printf 'genome_id\tinput_file\ttaxonomy\taligner\tthreads\tgapfill_medium_source\tdatabase_source\tmodel_file\treconstruction_mode\tdefault_doall_exit_status\teffective_gapfill_minimum_growth\teffective_gapfill_solver\n'
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$genome_id" "$input" "$taxonomy" "$aligner" "$threads" "gapseq_doall_default" "gapseq_environment_default" "${genome_id}.xml" \
    "$reconstruction_mode" "$default_status" "$effective_min_growth" "$effective_gapfill_solver"
} > "$outdir/reconstruction_provenance.tsv"
