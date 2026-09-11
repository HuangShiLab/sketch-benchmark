#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# Kraken2. T4: existing T2T-only (or kraken16 mixed) DB; host = reads assigned to 9606 or a descendant.
# T2: GTDB-5k DB built here (headers rewritten with kraken:taxid from GTDB metadata ncbi_taxid) + Bracken.
need kraken2
CONF=$(param conf 0); DBSEL=$(param db t2t)
if [ "$MODE" = build ]; then
  [ "$TASK" = T2 ] || exit 0
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  python "$REPO_DIR/scripts/tools/py/kraken_library.py" --refs "$IDX/refs.list" --metadata "$SB_REFS/bac120_metadata_r202.tsv" --out "$IDX/library.fa"
  mkdir -p "$IDX/db/taxonomy"; cp "$SOURCE_TAXONOMY"/{names,nodes}.dmp "$IDX/db/taxonomy/"
  kraken2-build --add-to-library "$IDX/library.fa" --db "$IDX/db" --threads "$THREADS" >/dev/null 2>&1
  kraken2-build --build --db "$IDX/db" --threads "$THREADS" --kmer-len 35 --minimizer-len 31 >/dev/null 2>&1
  bracken-build -d "$IDX/db" -t "$THREADS" -k 35 -l 150 >/dev/null 2>&1 || log "bracken-build failed (optional)"
  kraken2-build --clean --db "$IDX/db" >/dev/null 2>&1; rm -f "$IDX/library.fa"
  # genome -> taxid map for detect.tsv
  cut -f1,2 "$IDX/library.map" > "$IDX/genome_taxid.tsv" 2>/dev/null || true
else
  W="$6"; IDX="$7"; OUT="$8"
  DB="$IDX"; [ -d "$IDX/db" ] && DB="$IDX/db"
  if [ -n "$IN2" ]; then PAIRED=(--paired "$IN1" "$IN2"); else PAIRED=("$IN1"); fi
  kraken2 --db "$DB" --threads "$THREADS" --confidence "$CONF" --gzip-compressed --report "$OUT/report.txt" \
      --output "$OUT/kraken.txt" "${PAIRED[@]}" >/dev/null 2>&1
  if [ "$TASK" = T4 ]; then
    # host taxid set: 9606 and descendants (taxonkit against the DB's taxonomy)
    if command -v taxonkit >/dev/null && [ -d "$DB/taxonomy" ]; then taxonkit list --data-dir "$DB/taxonomy" --ids 9606 --indent "" 2>/dev/null | awk 'NF' > "$W/host_taxids"; fi
    [ -s "$W/host_taxids" ] || echo 9606 > "$W/host_taxids"
    awk 'NR==FNR{h[$1]=1; next} $1=="C" && ($3 in h) {print $2}' "$W/host_taxids" "$OUT/kraken.txt" | norm_ids | sort -u > "$W/host_ids"
    all_ids "$W/all_ids"; write_calls "$W/all_ids" "$W/host_ids" "$OUT/calls.tsv"
    rm -f "$OUT/kraken.txt"
  else
    # T2: per-genome read counts from taxid (species-level via Bracken when available)
    python "$REPO_DIR/scripts/tools/py/kraken_detect.py" --report "$OUT/report.txt" --map "$IDX/genome_taxid.tsv" --out "$OUT/detect.tsv"
    rm -f "$OUT/kraken.txt"
  fi
fi
