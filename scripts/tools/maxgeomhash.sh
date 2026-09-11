#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# maxgeomhash — research code (2024/2025), no stable CLI known at plan time.
# W1: install from GitHub, then implement build/query here following mash.sh;
# output pairs.tsv with the standard columns. Until then the wrapper exits 3 so
# stage 04 records a clean failure instead of a fake result.
echo "maxgeomhash: not implemented — see comment" >&2; exit 3
