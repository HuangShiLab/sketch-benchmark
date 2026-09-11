#!/bin/bash
# =============================================================================
# Sketch benchmark — driver
# =============================================================================
# Two phases, because array sizes for stages 03–05 depend on what stage 02
# produced:
#   phase sim : 00 fetch -> 02 simulate -> 02b manifests   (02b re-invokes
#               `run_all.sh --phase run` when it finishes, unless --no-continue)
#   phase run : 03 build -> 04 query -> 05 metrics -> 06 figures, arrays sized
#               from data/T*/manifest.tsv
#
#   bash scripts/run_all.sh --minimal            # stage 01: everything at 10%
#   bash scripts/run_all.sh                      # full: phase sim, then run
#   bash scripts/run_all.sh --phase run          # after sim finished / to redo
#   bash scripts/run_all.sh --phase run --tasks T1,T4
#   bash scripts/run_all.sh --dry-run ...
#
# Submission blocks on the per-user job limit; run detached:
#   nohup bash scripts/run_all.sh > logs/submit.log 2>&1 &
# =============================================================================
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
export REPO_DIR="$REPO"

PHASE=all; TASKS="T1,T2,T3,T4"; export DRY_RUN=0; CONTINUE=1; AFTER=""
while [ $# -gt 0 ]; do
  case "$1" in
    --minimal) export SB_MINIMAL=1 ;;
    --dry-run) export DRY_RUN=1 ;;
    --phase) PHASE="$2"; shift ;;
    --tasks) TASKS="$2"; shift ;;
    --after) AFTER="$2"; shift ;;
    --no-continue) CONTINUE=0 ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac; shift
done
source "$REPO/scripts/hpc/config.sh"
export SB_TASKS="$TASKS"
SUB="$REPO/scripts/hpc/submit.sh"
S="$REPO/scripts"

lock() {
  [ "$DRY_RUN" = "1" ] && return 0
  local lf="$LOG_DIR/.run_all.lock"
  exec 9>"$lf"
  if ! flock -n 9; then echo "ERROR: another run_all.sh holds $lf" >&2; exit 1; fi
}
dep() { [ -n "$1" ] && echo "--dependency=afterok:$1" || true; }
count() { [ -f "$1" ] && wc -l < "$1" || echo 0; }
has() { echo ",$TASKS," | grep -q ",$1,"; }

lock
echo "phase=$PHASE tasks=$TASKS minimal=$SB_MINIMAL dry=$DRY_RUN"

if [ "$PHASE" = "all" ] || [ "$PHASE" = "sim" ]; then
  j00=$($SUB 00-fetch $(dep "$AFTER") --job-name=sb-00-fetch "$S/stage00_fetch.sh")
  echo "00 fetch      $j00"
  sims=""
  has T1 && { j=$($SUB 02-t1 $(dep "$j00") --job-name=sb-02-t1 "$S/stage02_simulate/t1.sh"); sims="$sims:$j"; echo "02 T1 mutate  $j"; }
  has T2 && { j=$($SUB 02-t2 $(dep "$j00") --job-name=sb-02-t2 "$S/stage02_simulate/t2.sh"); sims="$sims:$j"; echo "02 T2 sim     $j"; }
  has T3 && { j=$($SUB 02-t3 $(dep "$j00") --job-name=sb-02-t3 "$S/stage02_simulate/t3.sh"); sims="$sims:$j"; echo "02 T3 sim     $j"; }
  has T4 && { j=$($SUB 02-t4 $(dep "$j00") --job-name=sb-02-t4 "$S/stage02_simulate/t4.sh"); sims="$sims:$j"; echo "02 T4 sweeps  $j"; }
  sims="${sims#:}"
  export SB_CONTINUE="$CONTINUE"
  j02b=$($SUB 02b-manifests $(dep "$sims") --job-name=sb-02b-manifests "$S/stage02b_manifests.sh")
  echo "02b manifests $j02b  (continues into phase run: $CONTINUE)"
  [ "$PHASE" = "sim" ] || [ "$CONTINUE" = "1" ] && exit 0
fi

if [ "$PHASE" = "run" ]; then
  fig_deps=""
  for T in T1 T2 T3 T4; do
    has "$T" || continue
    nb=$(count "$SB_DATA/$T/build_manifest.tsv"); nq=$(count "$SB_DATA/$T/manifest.tsv")
    if [ "$nq" -eq 0 ]; then echo "$T: no manifest — run phase sim first" >&2; continue; fi
    j03=""
    if [ "$nb" -gt 0 ]; then
      j03=$($SUB 03-$T $(dep "$AFTER") --job-name=sb-03-$T --export=ALL,SB_TASK=$T --array=0-$((nb-1))%20 "$S/stage03_build.sh")
      echo "03 $T build   $j03  ($nb tasks)"
    fi
    j04=$($SUB 04-$T $(dep "$j03") --job-name=sb-04-$T --export=ALL,SB_TASK=$T --array=0-$((nq-1))%40 "$S/stage04_query.sh")
    echo "04 $T query   $j04  ($nq tasks)"
    j05=$($SUB 05-$T $(dep "$j04") --job-name=sb-05-$T --export=ALL,SB_TASK=$T "$S/stage05_metrics.sh")
    echo "05 $T metrics $j05"
    fig_deps="$fig_deps:$j05"
  done
  fig_deps="${fig_deps#:}"
  [ -n "$fig_deps" ] && { j06=$($SUB 06-figures $(dep "$fig_deps") --job-name=sb-06-figures "$S/stage06_figures.sh"); echo "06 figures    $j06"; }
fi
