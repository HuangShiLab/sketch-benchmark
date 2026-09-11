#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# minimap2 -x sr against the existing T2T .mmi; mapped = host.
need minimap2; need samtools
Q=$(param mapq 0)
[ "$MODE" = build ] && exit 0
W="$6"; IDX="$7"; OUT="$8"
MMI="$IDX"; [ -d "$IDX" ] && MMI="$IDX/$(basename "$MINIMAP2_INDEX")"
minimap2 -ax sr -t "$THREADS" --secondary=no "$MMI" "$IN1" ${IN2:+"$IN2"} 2>/dev/null \
  | samtools view -F 4 -q "$Q" - | cut -f1 | norm_ids | sort -u > "$W/host_ids"
all_ids "$W/all_ids"; write_calls "$W/all_ids" "$W/host_ids" "$OUT/calls.tsv"
