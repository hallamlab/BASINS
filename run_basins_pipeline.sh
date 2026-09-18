#!/usr/bin/env bash
# BASINS entrypoint; preserve the established workflow and resume/cache paths.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${SCRIPT_DIR}/run_basin_pipeline.sh" "$@"
