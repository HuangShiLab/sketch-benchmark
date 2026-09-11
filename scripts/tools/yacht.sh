#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# YACHT (F2 hypothesis test with stated FPR). VERIFY(W1): CLI of the installed yacht version.
need yacht
K=$(param k 31); SC=$(param scaled 1000); FPR=$(param fpr 0.05)
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  sourmash sketch dna -p "k=$K,scaled=$SC" --from-file "$IDX/refs.list" -o "$IDX/ref.zip" >/dev/null 2>&1
  yacht train --ref_file "$IDX/ref.zip" --ksize "$K" --num_threads "$THREADS" --ani_thresh 0.95 --outdir "$IDX/train" --prefix ref >/dev/null 2>&1
else
  W="$6"; IDX="$7"; OUT="$8"; mkdir -p "$OUT/sketch"
  sourmash sketch dna -p "k=$K,scaled=$SC,abund" --name "$NAME" -o "$OUT/sketch/sample.sig.zip" "$IN1" ${IN2:+"$IN2"} >/dev/null 2>&1
  yacht run --json "$IDX/train/ref_config.json" --sample_file "$OUT/sketch/sample.sig.zip" --significance "$(python -c "print(1-$FPR)")" \
      --num_threads "$THREADS" --min_coverage_list 1 0.5 0.1 0.05 --out "$OUT/yacht.xlsx" >/dev/null 2>&1
  python - "$OUT/yacht.xlsx" "$OUT/detect.tsv" <<'PY'
import sys, re
import pandas as pd
from pathlib import Path
def gid(p):
    n = Path(str(p)).name
    for e in (".fa.gz",".fna.gz",".fasta.gz",".fa",".fna",".fasta"):
        if n.endswith(e): n = n[:-len(e)]; break
    m = re.match(r"^(GC[AF]_\d+\.\d+)", n); return m.group(1) if m else n.split("_genomic")[0]
xl = pd.read_excel(sys.argv[1], sheet_name=None)
# lowest min_coverage sheet is the most permissive; report its in_sample_est with the coverage as score
sheet = sorted(xl, key=lambda s: float(re.sub(r"[^0-9.]", "", s) or 1))[0]
df = xl[sheet]
with open(sys.argv[2], "w") as out:
    out.write("genome\tscore\tcalled_present\tani_est\tcov_est\tabund_est\n")
    for _, r in df.iterrows():
        g = gid(r.get("organism_name", r.get("genome", "")))
        out.write(f"{g}\t{r.get('num_matches', 'NA')}\t{int(bool(r.get('in_sample_est', False)))}\tNA\tNA\tNA\n")
PY
fi
