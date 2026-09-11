#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# sourmash FracMinHash (F2). T1: sketch refset; query = Jaccard, containment, ANI+CI on listed pairs (Python API).
need sourmash
K=$(param k 21); SC=$(param scaled 1000)
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  sourmash sketch dna -p "k=$K,scaled=$SC,abund" --from-file "$IDX/refs.list" -o "$IDX/ref.zip" >/dev/null 2>&1
else
  W="$6"; IDX="$7"; OUT="$8"
  if [ "$KIND" = allvsall ]; then
    sourmash compare "$IDX/ref.zip" -k "$K" -p "$THREADS" --csv "$OUT/cmp.csv" >/dev/null 2>&1
    python "$REPO_DIR/scripts/tools/py/matrix_to_pairs.py" --matrix "$OUT/cmp.csv" --csv --task T1 --k "$K" --out "$OUT/pairs.tsv"
  else
    pairs_list "$W/pairs.list"
    python "$REPO_DIR/scripts/tools/py/sourmash_pairs.py" --sig "$IDX/ref.zip" --k "$K" --pairs "$W/pairs.list" --out "$OUT/pairs.tsv"
  fi
fi
