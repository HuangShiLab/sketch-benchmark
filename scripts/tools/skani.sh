#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# skani (sparse-chaining ANI; non-sketch baseline). T1 query only.
need skani
[ "$MODE" = build ] && exit 0
W="$6"; IDX="$7"; OUT="$8"
if [ "$KIND" = allvsall ]; then LIST="$IN1"; else ls "$IN1"/*.fa.gz > "$W/g.list"; LIST="$W/g.list"; fi
skani triangle -l "$LIST" -t "$THREADS" --sparse -o "$OUT/skani.tsv" >/dev/null 2>&1
# Ref_file Query_file ANI Align_fraction_ref Align_fraction_query Ref_name Query_name
awk -F'\t' 'BEGIN{OFS="\t"; print "genome_a","genome_b","j_est","ani_est","ci_lo","ci_hi"}
     NR>1 { a=$1; b=$2; sub(/.*\//,"",a); sub(/.*\//,"",b); print a, b, "NA", $3/100, "NA", "NA" }' "$OUT/skani.tsv" > "$OUT/pairs.tsv"
