#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# Hostile (its own defaults, T2T+HLA index); host = input ids minus cleaned-output ids.
need hostile
[ "$MODE" = build ] && exit 0
W="$6"; IDX="$7"; OUT="$8"
if [ -n "$IN2" ]; then hostile clean --fastq1 "$IN1" --fastq2 "$IN2" --index "$IDX" --threads "$THREADS" --out-dir "$W/clean" > "$OUT/hostile.json" 2>"$OUT/hostile.log"
else hostile clean --fastq1 "$IN1" --index "$IDX" --threads "$THREADS" --out-dir "$W/clean" > "$OUT/hostile.json" 2>"$OUT/hostile.log"; fi
all_ids "$W/all_ids"
fq_ids "$(ls "$W"/clean/*1.fastq.gz "$W"/clean/*.clean.fastq.gz 2>/dev/null | head -1)" | sort -u > "$W/kept_ids"
comm -23 "$W/all_ids" "$W/kept_ids" > "$W/host_ids"
write_calls "$W/all_ids" "$W/host_ids" "$OUT/calls.tsv"
