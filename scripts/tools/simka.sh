#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# Simka (count-based, all k-mers). Build does the whole computation; query formats the matrix.
need simka
K=$(param k 31)
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/samples.tsv"
  awk -F'\t' '{printf "%s: %s", $1, $2; if ($3!="") printf " ; %s", $3; print ""}' "$IDX/samples.tsv" > "$IDX/simka_input.txt"
  simka -in "$IDX/simka_input.txt" -out "$IDX/out" -out-tmp "$IDX/tmp" -kmer-size "$K" -abundance-min 2 -nb-cores "$THREADS" -max-memory 48000 >/dev/null 2>&1
  rm -rf "$IDX/tmp"
else
  W="$6"; IDX="$7"; OUT="$8"
  zcat -f "$IDX/out/mat_abundance_braycurtis.csv.gz" > "$OUT/bc.csv"
  python - "$OUT/bc.csv" "$OUT/dist.tsv" <<'PY'
import sys, csv
rows = list(csv.reader(open(sys.argv[1]), delimiter=";")); names = rows[0][1:]
with open(sys.argv[2], "w") as out:
    out.write("sample_a\tsample_b\td_est\testimand\n")
    for i, r in enumerate(rows[1:]):
        for j in range(i + 1, len(names)):
            out.write(f"{names[i]}\t{names[j]}\t{float(r[j+1]):.6f}\tbray_curtis\n")
PY
fi
