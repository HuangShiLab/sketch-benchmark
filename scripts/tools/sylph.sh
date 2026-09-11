#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# sylph (F2 + coverage model). T1: db of refset; each genome queried as a "sample" -> ANI with 5-95 percentile.
# VERIFY(W1): `sylph query` accepts a genome FASTA as reads via -r; column names of the TSV.
need sylph
K=$(param k 21); C=$(param c 200)
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  sylph sketch -k "$K" -c "$C" -t "$THREADS" -l "$IDX/refs.list" -o "$IDX/ref" >/dev/null 2>&1   # -> ref.syldb
else
  W="$6"; IDX="$7"; OUT="$8"
  ls "$IN1"/*.fa.gz 2>/dev/null > "$W/q.list" || cat "$IN1" > "$W/q.list"
  sylph sketch -k "$K" -c "$C" -t "$THREADS" -r $(cat "$W/q.list") -d "$W/qs" >/dev/null 2>&1
  sylph query -t "$THREADS" "$W"/qs/*.sylsp -d "$IDX/ref.syldb" -o "$OUT/sylph.tsv" >/dev/null 2>&1
  python - "$OUT/sylph.tsv" "$OUT/pairs.tsv" <<'PY'
import sys, csv, re
from pathlib import Path
def gid(p):
    n = Path(p).name
    for e in (".fa.gz",".fna.gz",".fasta.gz",".sylsp",".fa",".fna",".fasta"):
        if n.endswith(e): n = n[:-len(e)]; break
    n = n.replace(".fastq.gz","")
    m = re.match(r"^(GC[AF]_\d+\.\d+)", n); return m.group(1) if m else n.split("_genomic")[0]
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
with open(sys.argv[2], "w") as out:
    out.write("genome_a\tgenome_b\tj_est\tani_est\tci_lo\tci_hi\n")
    for r in rows:
        a, b = gid(r["Sample_file"]), gid(r["Genome_file"])
        if a == b: continue
        lo, hi = (r.get("ANI_5-95_percentile", "NA-NA") or "NA-NA").split("-")[:2]
        out.write(f"{a}\t{b}\tNA\t{float(r['Adjusted_ANI'])/100:.6f}\t{float(lo)/100 if lo!='NA' else 'NA'}\t{float(hi)/100 if hi!='NA' else 'NA'}\n")
PY
fi
