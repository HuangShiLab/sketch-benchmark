#!/bin/bash
# Blocking sbatch: retries while the per-user submit limit refuses the job,
# prints the job id, appends (stage, jobid) to logs/jobs.tsv.
#   scripts/hpc/submit.sh <stage-label> [sbatch args...] script.sh
set -euo pipefail
source "$(dirname "$0")/config.sh"
LABEL="$1"; shift
if [ "${DRY_RUN:-0}" = "1" ]; then
  echo "DRY-RUN sbatch --output=$LOG_DIR/%x-%j.out --error=$LOG_DIR/%x-%j.err $*" >&2
  echo "dry-$LABEL"; exit 0
fi
while :; do
  out=$(sbatch --parsable --partition="$SB_PARTITION" --qos="$SB_QOS" \
        --output="$LOG_DIR/%x-%A_%a.out" --error="$LOG_DIR/%x-%A_%a.err" "$@" 2>&1) && break
  if echo "$out" | grep -qiE "MaxSubmitJob|QOSMaxSubmitJob|Job violates|limit"; then
    echo "[$LABEL] submit limit; retry in 120s" >&2; sleep 120; continue
  fi
  echo "[$LABEL] sbatch failed: $out" >&2; exit 1
done
jid="${out%%;*}"
printf '%s\t%s\t%s\n' "$LABEL" "$jid" "$(date -Iseconds)" >> "$LOG_DIR/jobs.tsv"
echo "$jid"
