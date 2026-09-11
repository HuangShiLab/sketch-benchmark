#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/../hpc/config.sh"; source "$(dirname "$0")/common.sh"
MODE="$1"; TASK="$2"; NAME="$3"; export PARAMS="$4"; THREADS="$5"
# Dashing 2 (F4 SetSketch; --bottom-k and --pminhash modes give F1/F3 rows).
# VERIFY(W1): dashing2 CLI flags (`dashing2 sketch --help`): -k, -S (log2 sketch size),
# --cmpout, -F (file list), --bottom-k / --pminhash / --fullsetsketch, --fastq for reads.
need dashing2
K=$(param k 21); L2=$(param log2m 12); MODE_=$(param mode setsketch)
case "$MODE_" in setsketch) MF="--fullsetsketch";; bottomk) MF="--bottom-k";; pminhash) MF="--pminhash";; *) die "mode $MODE_";; esac
if [ "$MODE" = build ]; then
  IDX="$6"; refset_list "$TASK" "$NAME" "$IDX/refs.list"
  case "$TASK" in
    T1) dashing2 sketch -k "$K" -S "$L2" $MF -p "$THREADS" --cache -F "$IDX/refs.list" -o "$IDX/ref" ;;
    T3) awk -F'\t' '{print $2}' "$IDX/refs.list" > "$IDX/r1.list"   # mate 1 only; dashing2 has no paired mode
        dashing2 sketch --fastq -k "$K" -S "$L2" $MF -p "$THREADS" --cache -F "$IDX/r1.list" -o "$IDX/samples" ;;
  esac
else
  W="$6"; IDX="$7"; OUT="$8"
  case "$TASK" in
    T1) dashing2 sketch -k "$K" -S "$L2" $MF -p "$THREADS" --cache -F "$IDX/refs.list" --cmpout "$OUT/cmp.tsv" --similarity-measure jaccard >/dev/null 2>&1
        python "$REPO_DIR/scripts/tools/py/matrix_to_pairs.py" --matrix "$OUT/cmp.tsv" --names "$IDX/refs.list" --task T1 --k "$K" --out "$OUT/pairs.tsv" ;;
    T3) dashing2 sketch --fastq -k "$K" -S "$L2" $MF -p "$THREADS" --cache -F "$IDX/r1.list" --cmpout "$OUT/cmp.tsv" --similarity-measure jaccard >/dev/null 2>&1
        python "$REPO_DIR/scripts/tools/py/matrix_to_pairs.py" --matrix "$OUT/cmp.tsv" --names "$IDX/refs.list" --task T3 --estimand "$([ "$MODE_" = pminhash ] && echo prob_jaccard || echo jaccard)" --out "$OUT/dist.tsv" ;;
  esac
fi
