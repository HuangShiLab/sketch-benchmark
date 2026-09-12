#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# Generic in-silico 2b-RAD tag sampler (F2, motif-defined context-free sampling; the
# *mechanism* row — Fast2bRAD-M / Syn2bANI are the tool rows). Tool-free Python.
#   T1: tag-set containment -> ANI per pair          params enz=BcgI[+AlfI]
#   T2: reference tag db, containment per genome    params enz=BcgI,min_hits=5,min_frac=0.001
#   T3: per-sample tag counters, Bray-Curtis/Jaccard params enz=BcgI,abund=1
PYI="$REPO_DIR/scripts/tools/py/i2brad.py"
ENZ=$(param enz BcgI)
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  case "$TASK" in
    T1) python "$PYI" build-db --enz "$ENZ" --threads "$THREADS" --genomes "$IDX/refs.list" -o "$IDX/db" ;;   # bytes of the reference sketch
    T2) python "$PYI" build-db --enz "$ENZ" --threads "$THREADS" --genomes "$IDX/refs.list" -o "$IDX/db" ;;
    T3) while IFS=$'\t' read -r sid r1 r2; do
          python "$PYI" sketch --enz "$ENZ" --reads "$r1" ${r2:+"$r2"} -o "$IDX/$sid.tags.tsv.gz"
        done < "$IDX/refs.list" ;;
    *) die "i2brad: not a contestant for $TASK" ;;
  esac
else
  W="$6"; IDX="$7"; OUT="$8"
  case "$TASK" in
    T1) # $IN1 is a directory (pairs kind) or a list file (allvsall kind); both accepted by --genomes
        python "$PYI" pairs --enz "$ENZ" --threads "$THREADS" --genomes "$IN1" -o "$OUT/pairs.tsv" ;;
    T2) MH=$(param min_hits 5); MF=$(param min_frac 0.001)
        python "$PYI" detect --enz "$ENZ" --db "$IDX/db" --reads "$IN1" ${IN2:+"$IN2"} --min-hits "$MH" --min-frac "$MF" --read-len "$READ_LENGTH" -o "$OUT/detect.tsv" ;;
    T3) AB=$(param abund 1)
        ls "$IDX"/*.tags.tsv.gz | awk '{n=$0; sub(/.*\//,"",n); sub(/\.tags\.tsv\.gz$/,"",n); print n"\t"$0}' > "$W/samples.tsv"
        python "$PYI" bc --samples "$W/samples.tsv" --abund "$AB" -o "$OUT/dist.tsv" ;;
    *) die "i2brad: not a contestant for $TASK" ;;
  esac
fi
