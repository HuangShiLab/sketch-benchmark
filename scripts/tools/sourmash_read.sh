#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# CAUTIONARY ROW: FracMinHash forced to minimizer density (scaled=8) per read against T2T.
# Same index bytes as the minimizer index, no window guarantee. Python; capped reads.
need sourmash
K=$(param k 31); SC=$(param scaled 8); A=$(param a 2)
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  sourmash sketch dna -p "k=$K,scaled=$SC" --from-file "$IDX/refs.list" -o "$IDX/ref.zip" >/dev/null 2>&1
  python - "$IDX/ref.zip" "$K" "$IDX/hashes.npy" <<'PY'
import sys, numpy as np, sourmash
hs = set()
for s in sourmash.load_file_as_signatures(sys.argv[1], ksize=int(sys.argv[2])): hs.update(s.minhash.hashes)
np.save(sys.argv[3], np.array(sorted(hs), dtype=np.uint64))
PY
else
  W="$6"; IDX="$7"; OUT="$8"
  python "$REPO_DIR/scripts/tools/py/read_sketch_filter.py" --mode sourmash --reads "$IN1" ${IN2:+--reads2 "$IN2"} --db "$IDX/hashes.npy" \
      --k "$K" --scaled "$SC" --abs "$A" --max-reads "$SB_T4_CAUTION_MAXREADS" --workdir "$W" --out "$OUT/calls.tsv"
fi
