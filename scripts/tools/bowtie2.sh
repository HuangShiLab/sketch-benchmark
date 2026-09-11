#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# Bowtie2 against the existing T2T index (Hostile-style local alignment); mapped (MAPQ >= q) = host.
need bowtie2; need samtools
Q=$(param mapq 0)
[ "$MODE" = build ] && exit 0
W="$6"; IDX="$7"; OUT="$8"
PREFIX="$IDX/$(basename "$BOWTIE2_INDEX")"
if [ -n "$IN2" ]; then IN=(-1 "$IN1" -2 "$IN2"); else IN=(-U "$IN1"); fi
bowtie2 -x "$PREFIX" --very-sensitive-local -k 1 -p "$THREADS" "${IN[@]}" 2>"$OUT/bowtie2.log" \
  | samtools view -F 4 -q "$Q" - | cut -f1 | norm_ids | sort -u > "$W/host_ids"
all_ids "$W/all_ids"; write_calls "$W/all_ids" "$W/host_ids" "$OUT/calls.tsv"
