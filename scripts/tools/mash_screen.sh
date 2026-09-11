#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# Mash Screen (F1, asymmetric: sketch refs, stream all sample k-mers).
need mash
K=$(param k 31); S=$(param s 10000); I=$(param i 0.9)
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  mash sketch -k "$K" -s "$S" -p "$THREADS" -l "$IDX/refs.list" -o "$IDX/ref"
else
  W="$6"; IDX="$7"; OUT="$8"
  mash screen -p "$THREADS" -w "$IDX/ref.msh" "$IN1" ${IN2:+"$IN2"} > "$OUT/screen.tsv" 2>/dev/null
  # identity shared-hashes median-multiplicity p-value query-ID [comment]
  awk -F'\t' -v thr="$I" 'BEGIN{OFS="\t"; print "genome","score","called_present","ani_est","cov_est","abund_est"}
       { g=$5; sub(/.*\//,"",g); sub(/(_genomic)?\.(fna|fa|fasta)(\.gz)?$/,"",g); if (match(g,/^GC[AF]_[0-9]+\.[0-9]+/)) g=substr(g,RSTART,RLENGTH);
         print g, $1, ($1>=thr)?1:0, $1, $3, "NA" }' "$OUT/screen.tsv" > "$OUT/detect.tsv"
fi
