#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# T4: build a host index (T2T-only), filter in keep-matching mode, matched read ids = host.
# deacon_panhuman uses the prebuilt panhuman-1 index (no build). deacon_syncmer uses the
# patched binary with --scheme syncmer (k = key length l, w = s-mer length).
BIN="${DEACON_BIN:-deacon}"; need "$BIN"
K=$(param k 31); W_=$(param w 15); S_=$(param s ""); A=$(param a 2); RR=$(param r 0.01)
if [ "$MODE" = build ]; then
  [ "$BIN" = deacon ] && [ "${DEACON_PREBUILT:-0}" = 1 ] && exit 0
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  if [ -n "$S_" ]; then "$BIN" index build "$(cat "$IDX/refs.list")" -k "$K" -w "$S_" --scheme syncmer -t "$THREADS" -o "$IDX/host.idx" -q
  else "$BIN" index build "$(cat "$IDX/refs.list")" -k "$K" -w "$W_" -t "$THREADS" -o "$IDX/host.idx" -q; fi
else
  W="$6"; IDX="$7"; OUT="$8"
  IDXFILE="$IDX/host.idx"; [ -f "$IDX" ] && IDXFILE="$IDX"     # prebuilt index is a file
  "$BIN" filter "$IDXFILE" "$IN1" ${IN2:+"$IN2"} -a "$A" -r "$RR" -t "$THREADS" -o "$W/hit.fq" --summary "$OUT/summary.json" >/dev/null 2>&1
  awk 'NR%4==1' "$W/hit.fq" | norm_ids | sort -u > "$W/host_ids"; rm -f "$W/hit.fq"
  all_ids "$W/all_ids"; write_calls "$W/all_ids" "$W/host_ids" "$OUT/calls.tsv"
fi
