#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# sylph profile (F2 + zero-inflated Poisson coverage). Build: .syldb of GTDB-5k. Query: sketch sample, profile.
need sylph
K=$(param k 31); C=$(param c 200); MNK=$(param min_kmers 50)
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  sylph sketch -k "$K" -c "$C" -t "$THREADS" -l "$IDX/refs.list" -o "$IDX/ref" >/dev/null 2>&1
else
  W="$6"; IDX="$7"; OUT="$8"; mkdir -p "$OUT/sketch"
  if [ -n "$IN2" ]; then sylph sketch -k "$K" -c "$C" -t "$THREADS" -1 "$IN1" -2 "$IN2" -d "$OUT/sketch" >/dev/null 2>&1
  else sylph sketch -k "$K" -c "$C" -t "$THREADS" -r "$IN1" -d "$OUT/sketch" >/dev/null 2>&1; fi
  sylph profile -t "$THREADS" --min-number-kmers "$MNK" "$IDX/ref.syldb" "$OUT"/sketch/*.sylsp -o "$OUT/profile.tsv" >/dev/null 2>&1
  python - "$OUT/profile.tsv" "$OUT/detect.tsv" <<'PY'
import csv, sys, re
from pathlib import Path
def gid(p):
    n = Path(p).name
    for e in (".fa.gz",".fna.gz",".fasta.gz",".fa",".fna",".fasta"):
        if n.endswith(e): n = n[:-len(e)]; break
    m = re.match(r"^(GC[AF]_\d+\.\d+)", n); return m.group(1) if m else n.split("_genomic")[0]
out = open(sys.argv[2], "w"); out.write("genome\tscore\tcalled_present\tani_est\tcov_est\tabund_est\n")
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    out.write(f"{gid(r['Genome_file'])}\t{r.get('Eff_cov','NA')}\t1\t{float(r['Adjusted_ANI'])/100:.6f}\t{r.get('Eff_cov','NA')}\t{r.get('Sequence_abundance','NA')}\n")
PY
fi
