#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# sylph profile -> Bray-Curtis on species abundance (the "what people actually compute" baseline).
need sylph
K=$(param k 31); C=$(param c 200)
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/samples.tsv"
  cat "$T3_POOL_LIST" > "$IDX/refs.list"   # the T3 pool is the reference set for profiling
  sylph sketch -k "$K" -c "$C" -t "$THREADS" -l "$IDX/refs.list" -o "$IDX/ref" >/dev/null 2>&1
  while IFS=$'\t' read -r sid r1 r2; do
    if [ -n "$r2" ]; then sylph sketch -k "$K" -c "$C" -t "$THREADS" -1 "$r1" -2 "$r2" -d "$IDX/sk_$sid" >/dev/null 2>&1
    else sylph sketch -k "$K" -c "$C" -t "$THREADS" -r "$r1" -d "$IDX/sk_$sid" >/dev/null 2>&1; fi
  done < "$IDX/samples.tsv"
else
  W="$6"; IDX="$7"; OUT="$8"
  sylph profile -t "$THREADS" "$IDX/ref.syldb" "$IDX"/sk_*/*.sylsp -o "$OUT/profile.tsv" >/dev/null 2>&1
  python - "$OUT/profile.tsv" "$IDX/samples.tsv" "$OUT/dist.tsv" <<'PY'
import sys, csv, itertools, collections
from pathlib import Path
prof = collections.defaultdict(dict)
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    s = Path(r["Sample_file"]).name.split(".")[0]; prof[s][r["Genome_file"]] = float(r["Sequence_abundance"])
# map sample file stem -> sample id via samples.tsv (r1 basename stem)
sid_of = {}
for l in open(sys.argv[2]):
    sid, r1, *_ = l.rstrip("\n").split("\t"); sid_of[Path(r1).name.split(".")[0]] = sid
names = sorted(prof)
with open(sys.argv[3], "w") as out:
    out.write("sample_a\tsample_b\td_est\testimand\n")
    for a, b in itertools.combinations(names, 2):
        keys = set(prof[a]) | set(prof[b]); x = [prof[a].get(k, 0) for k in keys]; y = [prof[b].get(k, 0) for k in keys]
        bc = sum(abs(p - q) for p, q in zip(x, y)) / max(1e-9, sum(p + q for p, q in zip(x, y)))
        out.write(f"{sid_of.get(a,a)}\t{sid_of.get(b,b)}\t{bc:.6f}\tbray_curtis_profile\n")
PY
fi
