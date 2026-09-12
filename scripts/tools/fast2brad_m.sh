#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# Fast2bRAD-M (F2 motif-defined tags; HuangShiLab/Fast2bRAD-M, Rust port of 2bRAD-M / MAP2B).
#   T2 build: qualitative tag database of the GTDB-5k pool  (build-qual-db, species level)
#   T2 query: extract tags from reads (-t 2) -> quantify -> per-genome detection
#   T3 build: the same database (reused from T2 when present)
#   T3 query: per-sample quantify -> species profile -> Bray-Curtis (estimand bc_profile)
# params enz=BcgI,g=5   (G-score threshold)
# VERIFY(W1): column names of {sample}.{enz}.GCF_detected.xls and .xls (README documents the
#   abundance table; the per-genome detail file is parsed by header sniffing); `quantify -l` list
#   format `sample<TAB>path.iibsp`; GTDB-format taxonomy file is auto-detected.
need fast2bRAD-M
ENZ=$(param enz BcgI); G=$(param g 5)
META="$SB_REFS/bac120_metadata_r202.tsv"
taxonomy_files() {   # genome_list.tsv (gid path) + taxonomy.tsv (gid gtdb_taxonomy) for a list of FASTA paths
  local list="$1" out="$2"
  python - "$list" "$META" "$out" "$REPO_DIR/scripts/tools/py" <<'PY'
import csv, sys
from pathlib import Path
sys.path.insert(0, sys.argv[4]); from names import genome_name
tax = {}
try:
    for r in csv.DictReader(open(sys.argv[2]), delimiter="\t"):
        tax[r["accession"].replace("RS_", "").replace("GB_", "")] = r["gtdb_taxonomy"]
except FileNotFoundError: pass
out = Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=True)
with open(out / "genome_list.tsv", "w") as gl, open(out / "taxonomy.tsv", "w") as tx:
    for l in open(sys.argv[1]):
        p = l.strip()
        if not p: continue
        g = genome_name(p); gl.write(f"{g}\t{p}\n")
        tx.write(f"{g}\t{tax.get(g, 'd__Unknown;p__Unknown;c__Unknown;o__Unknown;f__Unknown;g__Unknown;s__' + g)}\n")
PY
}
if [ "$MODE" = build ]; then
  IDX="$6"
  case "$TASK" in
    T2|T3)
      T2IDX="$SB_IDX/T2/fast2brad_m/gtdb5k/${PARAMS//[,=]/_}"
      if [ "$TASK" = T3 ] && [ -f "$T2IDX/DONE" ]; then ln -sfn "$T2IDX/db" "$IDX/db"; cp "$T2IDX"/genome_list.tsv "$T2IDX"/taxonomy.tsv "$IDX/"; exit 0; fi
      refset_list T2 gtdb5k "$IDX/refs.list"; taxonomy_files "$IDX/refs.list" "$IDX"
      fast2bRAD-M build-qual-db -l "$IDX/genome_list.tsv" --taxonomy "$IDX/taxonomy.tsv" -s "$ENZ" -t species -o "$IDX/db" -r yes -j "$THREADS" >/dev/null 2>&1 ;;
    *) die "fast2brad_m: not a contestant for $TASK" ;;
  esac
else
  W="$6"; IDX="$7"; OUT="$8"; mkdir -p "$OUT/sketch"
  quantify_one() {   # <sid> <r1> [r2] -> $OUT/q/<sid>/<sid>.<ENZ>.xls ; the .iibsp tag file is the query-side sketch
    local sid="$1" r1="$2" r2="${3:-}"
    fast2bRAD-M extract -i "$r1" ${r2:+"$r2"} -t 2 -s "$ENZ" --od "$OUT/sketch" --op "$sid" -j "$THREADS" --qc yes >/dev/null 2>&1
    printf "%s\t%s\n" "$sid" "$OUT/sketch/$sid.$ENZ.iibsp" > "$W/$sid.list"
    fast2bRAD-M quantify -l "$W/$sid.list" -d "$IDX/db" -t species -s "$ENZ" -o "$OUT/q" -g "$G" -v yes -j "$THREADS" >/dev/null 2>&1
  }
  case "$TASK" in
    T2) quantify_one sample "$IN1" ${IN2:+"$IN2"}
        python - "$OUT/q/sample" "$ENZ" "$IDX/taxonomy.tsv" "$G" "$OUT/detect.tsv" <<'PY'
import csv, glob, sys, os
d, enz, taxf, gthr, out = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
species = {}                                                    # species name -> [genome ids]
for l in open(taxf):
    g, t = l.rstrip("\n").split("\t"); species.setdefault(t.split(";")[-1].replace("s__", ""), []).append(g)
rows = []
det = glob.glob(os.path.join(d, f"*.{enz}.GCF_detected.xls"))
if det:                                                         # per-genome detail: sniff columns
    r = list(csv.DictReader(open(det[0]), delimiter="\t"))
    if r:
        cols = {c.lower(): c for c in r[0]}
        gcol = next((cols[c] for c in cols if "gcf" in c or "genome" in c or "assembly" in c), None)
        tcol = next((cols[c] for c in cols if "sequenced_tag" in c or c.startswith("tag")), None)
        scol = next((cols[c] for c in cols if "g_score" in c or "gscore" in c), None)
        pcol = next((cols[c] for c in cols if "percent" in c or "abund" in c), None)
        for x in r:
            if not gcol: break
            g = x[gcol]; sc = float(x[scol]) if scol and x[scol] not in ("", "NA") else float("nan")
            rows.append((g, sc, int(sc >= gthr) if sc == sc else 1, "NA", x.get(tcol, "NA") if tcol else "NA", (float(x[pcol]) / 100 if pcol and x[pcol] not in ("", "NA") else "NA")))
if not rows:                                                    # fall back to the species table mapped onto pool genomes
    for f in glob.glob(os.path.join(d, f"*.{enz}.xls")):
        if f.endswith("GCF_detected.xls"): continue
        for x in csv.DictReader(open(f), delimiter="\t"):
            sp = x.get("Species", ""); sc = float(x.get("G_Score", "nan") or "nan"); pc = x.get("Percent", "NA")
            for g in species.get(sp, []):
                rows.append((g, sc, int(sc >= gthr) if sc == sc else 1, "NA", x.get("Sequenced_Tag_Num", "NA"), float(pc) / 100 if pc not in ("", "NA") else "NA"))
with open(out, "w") as fh:
    fh.write("genome\tscore\tcalled_present\tani_est\tcov_est\tabund_est\n")
    for g, sc, c, a, cov, ab in rows: fh.write(f"{g}\t{sc}\t{c}\t{a}\t{cov}\t{ab}\n")
PY
        ;;
    T3) awk -F'\t' 'NR>1{print $1"\t"$3"\t"$4}' "$SB_DATA/T3/datasets.tsv" > "$W/samples.tsv"
        while IFS=$'\t' read -r sid r1 r2; do quantify_one "$sid" "$r1" ${r2:+"$r2"}; done < "$W/samples.tsv"
        python - "$OUT/q" "$ENZ" "$W/samples.tsv" "$OUT/dist.tsv" <<'PY'
import csv, sys, os, glob
import numpy as np
qd, enz, sl, out = sys.argv[1:5]
sids = [l.split("\t")[0] for l in open(sl) if l.strip()]
prof = {}
for sid in sids:
    v = {}
    for f in glob.glob(os.path.join(qd, sid, f"*.{enz}.xls")):
        if f.endswith("GCF_detected.xls"): continue
        for x in csv.DictReader(open(f), delimiter="\t"):
            try: v[x["Species"]] = float(x["Percent"])
            except (KeyError, ValueError): pass
    prof[sid] = v
sp = sorted({s for v in prof.values() for s in v})
M = np.array([[prof[s].get(x, 0.0) for x in sp] for s in sids], float)
with open(out, "w") as fh:
    fh.write("sample_a\tsample_b\td_est\testimand\n")
    for i in range(len(sids)):
        for j in range(i + 1, len(sids)):
            den = M[i].sum() + M[j].sum(); d = np.abs(M[i] - M[j]).sum() / den if den else float("nan")
            fh.write(f"{sids[i]}\t{sids[j]}\t{d:.6f}\tbc_profile\n")
PY
        ;;
    *) die "fast2brad_m: not a contestant for $TASK" ;;
  esac
fi
