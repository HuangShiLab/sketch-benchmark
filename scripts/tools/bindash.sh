#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# BinDash (F1 one-permutation b-bit). T1 only.
# VERIFY(W1): flag names against `bindash sketch --help` (2.x); sketchsize64 = s/64.
need bindash
K=$(param k 21); S=$(param s 1000); B=$(( S / 64 ))
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  bindash sketch --kmerlen="$K" --sketchsize64="$B" --nthreads="$THREADS" --listfname="$IDX/refs.list" --outfname="$IDX/ref.sketch"
else
  W="$6"; IDX="$7"; OUT="$8"
  bindash dist --nthreads="$THREADS" "$IDX/ref.sketch" > "$OUT/bindash.tsv"
  # columns: query target mutation-distance p-value shared/total  (VERIFY)
  awk -F'\t' 'BEGIN{OFS="\t"; print "genome_a","genome_b","j_est","ani_est","ci_lo","ci_hi"}
       $1<$2 { split($5,sh,"/"); a=$1; b=$2; sub(/.*\//,"",a); sub(/.*\//,"",b); print a, b, sh[1]/sh[2], 1-$3, "NA", "NA" }' "$OUT/bindash.tsv" > "$OUT/pairs.tsv"
fi
