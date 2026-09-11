#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# minimap2 -x sr mapping -> breadth/depth per reference genome (alignment gold for T2, 10 M sets only).
need minimap2; need samtools
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  python "$REPO_DIR/scripts/tools/py/kraken_library.py" --refs "$IDX/refs.list" --out "$IDX/refs.fa" --no-taxid   # concatenates, contig->genome map
  minimap2 -x sr -t "$THREADS" -d "$IDX/refs.mmi" "$IDX/refs.fa" >/dev/null 2>&1; samtools faidx "$IDX/refs.fa"
else
  W="$6"; IDX="$7"; OUT="$8"
  minimap2 -ax sr -t "$THREADS" --secondary=no "$IDX/refs.mmi" "$IN1" ${IN2:+"$IN2"} 2>/dev/null \
    | samtools view -b -F 0x904 -q 20 - | samtools sort -@ 4 -o "$W/aln.bam" - 2>/dev/null
  samtools index "$W/aln.bam"; samtools coverage "$W/aln.bam" > "$OUT/coverage.tsv"
  python - "$OUT/coverage.tsv" "$IDX/refs.fa.map" "$OUT/detect.tsv" <<'PY'
import sys, csv, collections
cmap = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(sys.argv[2])}   # contig -> genome
L = collections.defaultdict(int); cov = collections.defaultdict(float); reads = collections.defaultdict(int)
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    g = cmap.get(r["#rname"], r["#rname"]); n = int(r["endpos"]) - int(r["startpos"]) + 1
    L[g] += n; cov[g] += float(r["meandepth"]) * n; reads[g] += int(r["numreads"])
with open(sys.argv[3], "w") as out:
    out.write("genome\tscore\tcalled_present\tani_est\tcov_est\tabund_est\n")
    for g in L:
        d = cov[g] / L[g]
        out.write(f"{g}\t{reads[g]}\t{int(d >= 0.01)}\tNA\t{d:.6g}\tNA\n")
PY
fi
