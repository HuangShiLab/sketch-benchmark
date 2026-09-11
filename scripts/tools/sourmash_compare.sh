#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# sourmash compare (F2). abund=1: angular similarity on abundance-weighted sketches; abund=0: Jaccard.
need sourmash
K=$(param k 31); SC=$(param scaled 1000); AB=$(param abund 1)
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/samples.tsv"
  while IFS=$'\t' read -r sid r1 r2; do
    sourmash sketch dna -p "k=$K,scaled=$SC,abund" --name "$sid" -o "$IDX/$sid.sig" "$r1" ${r2:+"$r2"} >/dev/null 2>&1
  done < "$IDX/samples.tsv"
else
  W="$6"; IDX="$7"; OUT="$8"
  FLAG=""; EST="angular"; [ "$AB" = 0 ] && { FLAG="--ignore-abundance"; EST="jaccard"; }
  sourmash compare "$IDX"/*.sig -k "$K" -p "$THREADS" $FLAG --csv "$OUT/cmp.csv" >/dev/null 2>&1
  python "$REPO_DIR/scripts/tools/py/matrix_to_pairs.py" --matrix "$OUT/cmp.csv" --csv --task T3 --estimand "$EST" --out "$OUT/dist.tsv"
fi
