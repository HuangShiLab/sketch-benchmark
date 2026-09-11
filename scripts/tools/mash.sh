#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# Mash (F1 bottom-k). T1: sketch refset, all-vs-all dist -> pairs.tsv.
# T3: sketch each sample's reads (-r, -m), all-vs-all -> dist.tsv (estimand mash_dist).
need mash
K=$(param k 21); S=$(param s 1000); M=$(param m 2)
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  case "$TASK" in
    T1) mash sketch -k "$K" -s "$S" -p "$THREADS" -l "$IDX/refs.list" -o "$IDX/ref" ;;
    T3) while IFS=$'\t' read -r sid r1 r2; do
          mash sketch -r -m "$M" -k "$K" -s "$S" -p "$THREADS" -o "$IDX/$sid" -I "$sid" "$r1" ${r2:+"$r2"}
        done < "$IDX/refs.list"
        mash paste "$IDX/all" "$IDX"/g*.msh "$IDX"/[!ag]*.msh 2>/dev/null || mash paste "$IDX/all" "$IDX"/*.msh ;;
    *) die "mash: not a contestant for $TASK" ;;
  esac
else
  W="$6"; IDX="$7"; OUT="$8"
  case "$TASK" in
    T1) # all-vs-all; sketch entries are named by file path, map to genome id
        mash dist -p "$THREADS" "$IDX/ref.msh" "$IDX/ref.msh" > "$OUT/mash.tsv"
        awk -F'\t' 'BEGIN{OFS="\t"; print "genome_a","genome_b","j_est","ani_est","ci_lo","ci_hi"}
             $1<$2 { split($5,sh,"/"); j=sh[1]/sh[2]; a=$1; b=$2; sub(/.*\//,"",a); sub(/.*\//,"",b);
                     print a, b, j, 1-$3, "NA", "NA" }' "$OUT/mash.tsv" > "$OUT/pairs.tsv" ;;
    T3) mash dist -p "$THREADS" "$IDX/all.msh" "$IDX/all.msh" > "$OUT/mash.tsv"
        awk -F'\t' 'BEGIN{OFS="\t"; print "sample_a","sample_b","d_est","estimand"} $1<$2 {print $1,$2,$3,"mash_dist"}' "$OUT/mash.tsv" > "$OUT/dist.tsv" ;;
    *) die "mash: not a contestant for $TASK" ;;
  esac
fi
