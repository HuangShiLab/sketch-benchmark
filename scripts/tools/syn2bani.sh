#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# Syn2bANI (F2 motif-defined tags + chain-restricted MLE ANI; HuangShiLab/Syn2bANI). T1 only.
# params panel=BcgI+AlfI+AloI+FalI  ('+' separated; the default 4-enzyme panel)
# Build: `syn2bani sketch` of the refset so sketch bytes are recorded. Query: `syn2bani dist` on FASTA.
# VERIFY(W1): ani columns are fractions (README: --min-ani is a fraction); std_err on the same scale;
#             `sketch` accepts --enzymes; `dist --ql/--rl` list files hold one FASTA path per line.
need syn2bani
PANEL=$(param panel "BcgI+AlfI+AloI+FalI"); PANEL="${PANEL//+/,}"
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  [ "$TASK" = T1 ] || die "syn2bani: not a contestant for $TASK"
  mkdir -p "$IDX/db"
  xargs -a "$IDX/refs.list" syn2bani sketch --enzymes "$PANEL" -o "$IDX/db" >/dev/null 2>&1
else
  W="$6"; IDX="$7"; OUT="$8"
  [ "$TASK" = T1 ] || die "syn2bani: not a contestant for $TASK"
  if [ "$KIND" = allvsall ]; then cp "$IN1" "$W/g.list"; else ls "$IN1"/*.fa.gz > "$W/g.list"; fi
  syn2bani dist --ql "$W/g.list" --rl "$W/g.list" --enzymes "$PANEL" --threads "$THREADS" -o "$OUT/syn2bani.tsv" >/dev/null 2>&1
  python - "$OUT/syn2bani.tsv" "$OUT/pairs.tsv" <<'PY'
import csv, sys, math
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
def f(x):
    try: v = float(x)
    except (TypeError, ValueError): return None
    return None if math.isnan(v) else v
with open(sys.argv[2], "w") as out:
    out.write("genome_a\tgenome_b\tj_est\tani_est\tci_lo\tci_hi\tflag\n")
    seen = set()
    for r in rows:
        a, b = r["query"], r["reference"]
        if a == b: continue
        key = tuple(sorted((a, b)))
        if key in seen: continue
        seen.add(key)
        ani = f(r.get("ani_gated")) or f(r.get("ani")); se = f(r.get("std_err"))
        if ani is None: continue
        if ani > 1.5: ani /= 100; se = se / 100 if se is not None else None     # VERIFY(W1): percent vs fraction
        lo = f"{max(0, ani - 1.96 * se):.6f}" if se is not None else "NA"; hi = f"{min(1, ani + 1.96 * se):.6f}" if se is not None else "NA"
        out.write(f"{a}\t{b}\tNA\t{ani:.6f}\t{lo}\t{hi}\t{r.get('flag', '')}\n")
PY
fi
