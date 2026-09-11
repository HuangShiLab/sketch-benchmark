#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# CAUTIONARY ROW: FracMinHash at sylph's default scale applied per read. Expect most reads
# to have no sketch k-mer at all (the estimand-mismatch demonstration). Capped at
# SB_T4_CAUTION_MAXREADS reads. Build: T2T .syldb.
need sylph; need seqkit
K=$(param k 31); C=$(param c 200)
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  sylph sketch -k "$K" -c "$C" -t "$THREADS" -l "$IDX/refs.list" -o "$IDX/ref" >/dev/null 2>&1
else
  W="$6"; IDX="$7"; OUT="$8"
  seqkit head -n "$SB_T4_CAUTION_MAXREADS" "$IN1" -o "$W/cap_1.fq.gz"
  # every read becomes its own "genome" query: split, sketch, query
  python "$REPO_DIR/scripts/tools/py/read_sketch_filter.py" --mode sylph --reads "$W/cap_1.fq.gz" --db "$IDX/ref.syldb" \
      --k "$K" --scaled "$C" --threads "$THREADS" --workdir "$W" --out "$OUT/calls.tsv"
fi
