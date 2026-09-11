#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# sourmash gather (F2 minimum set cover). Build: reference zip of GTDB-5k. Query: sketch reads, gather.
need sourmash
K=$(param k 31); SC=$(param scaled 1000); TBP=$(param threshold_bp 50000)
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  sourmash sketch dna -p "k=$K,scaled=$SC" --from-file "$IDX/refs.list" -o "$IDX/ref.zip" >/dev/null 2>&1
  sourmash index -k "$K" "$IDX/ref.sbt.zip" "$IDX/ref.zip" >/dev/null 2>&1 || true   # SBT speeds gather; optional
else
  W="$6"; IDX="$7"; OUT="$8"; mkdir -p "$OUT/sketch"
  sourmash sketch dna -p "k=$K,scaled=$SC,abund" --name "$NAME" -o "$OUT/sketch/sample.sig" "$IN1" ${IN2:+"$IN2"} >/dev/null 2>&1
  DB="$IDX/ref.zip"; [ -f "$IDX/ref.sbt.zip" ] && DB="$IDX/ref.sbt.zip"
  sourmash gather "$OUT/sketch/sample.sig" "$DB" -k "$K" --scaled "$SC" --threshold-bp "$TBP" -o "$OUT/gather.csv" >/dev/null 2>&1 || true
  python - "$OUT/gather.csv" "$OUT/detect.tsv" <<'PY'
import csv, sys, re, os
from pathlib import Path
def gid(p):
    n = Path(p).name
    for e in (".fa.gz",".fna.gz",".fasta.gz",".fa",".fna",".fasta"):
        if n.endswith(e): n = n[:-len(e)]; break
    m = re.match(r"^(GC[AF]_\d+\.\d+)", n); return m.group(1) if m else n.split("_genomic")[0]
out = open(sys.argv[2], "w"); out.write("genome\tscore\tcalled_present\tani_est\tcov_est\tabund_est\n")
if os.path.exists(sys.argv[1]):
    for r in csv.DictReader(open(sys.argv[1])):
        g = gid(r.get("filename") or r.get("name") or r.get("match_name", ""))
        if not g: g = gid(r.get("name", ""))
        ani = r.get("match_containment_ani") or r.get("query_containment_ani") or "NA"
        out.write(f"{g}\t{r.get('intersect_bp','NA')}\t1\t{ani}\t{r.get('average_abund','NA')}\t{r.get('f_unique_weighted','NA')}\n")
PY
fi
