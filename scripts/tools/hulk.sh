#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# HULK histosketch (F3, weighted Jaccard over the k-mer count histogram). VERIFY(W1): flags of the installed hulk.
need hulk
K=$(param k 31); S=$(param s 512)
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/samples.tsv"
  while IFS=$'\t' read -r sid r1 r2; do
    zcat -f "$r1" ${r2:+"$r2"} | hulk sketch -k "$K" -s "$S" -p "$THREADS" -o "$IDX/$sid" >/dev/null 2>&1
  done < "$IDX/samples.tsv"
else
  W="$6"; IDX="$7"; OUT="$8"
  hulk smash -d "$IDX" -m weightedjaccard -p "$THREADS" -o "$OUT/hulk" >/dev/null 2>&1   # -> hulk.csv matrix (VERIFY)
  python "$REPO_DIR/scripts/tools/py/matrix_to_pairs.py" --matrix "$OUT/hulk.csv" --csv --task T3 --estimand weighted_jaccard --out "$OUT/dist.tsv"
fi
