#!/bin/bash
#SBATCH --job-name=sb-04-query
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
# Stage 04 — one array task = one (tool, dataset, params, threads, rep) query
# from data/$SB_TASK/manifest.tsv. Inputs and index are copied to $TMPDIR
# before the clock starts (rule 5); set SB_NO_TMPDIR_COPY=1 to read in place.
set -euo pipefail
source "${REPO_DIR:?}/scripts/hpc/config.sh"; source "$CONDA_BASE/etc/profile.d/conda.sh"; conda activate "$SB_ENV"
source "$REPO_DIR/scripts/tools/common.sh"
T="${SB_TASK:?}"
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SB_DATA/$T/manifest.tsv"); [ -n "$line" ] || exit 0
IFS=$'\t' read -r TOOL DATASET PARAMS THREADS REP <<< "$line"
OUT="$SB_RUNS/$T/$TOOL/$DATASET/${PARAMS//[,=]/_}/t${THREADS}_r${REP}"; mkdir -p "$OUT"
[ -f "$OUT/DONE" ] && { echo "done: $OUT"; exit 0; }
ds=$(dataset_row "$T" "$DATASET"); [ -n "$ds" ] || die "dataset $DATASET not in datasets.tsv"
IFS=$'\t' read -r _ KIND P1 P2 TRUTH META <<< "$ds"

W="${TMPDIR:-/tmp}/sb.$SLURM_JOB_ID.${SLURM_ARRAY_TASK_ID:-0}"; mkdir -p "$W"
trap 'rm -rf "$W"' EXIT
IDX_SRC=$(index_for "$T" "$TOOL" "$PARAMS")     # built in stage 03, or an existing/prebuilt index
if [ "${SB_NO_TMPDIR_COPY:-0}" = "1" ]; then
  IDX="$IDX_SRC"; IN1="$P1"; IN2="$P2"
else
  if [ -n "$IDX_SRC" ] && [ -e "$IDX_SRC" ]; then cp -r "$IDX_SRC" "$W/idx"; IDX="$W/idx"; else IDX="$IDX_SRC"; fi
  case "$KIND" in
    pairs)      mkdir -p "$W/data"
                if [ -d "$P1" ]; then cp "$P1"/*.fa.gz "$W/data/"; else while read -r f; do cp "$f" "$W/data/"; done < "$P1"; fi
                IN1="$W/data"; IN2="" ;;
    allvsall)   mkdir -p "$W/data"; while read -r f; do cp "$f" "$W/data/"; echo "$W/data/$(basename "$f")"; done < "$P1" > "$W/data.list"; IN1="$W/data.list"; IN2="" ;;
    *)          [ -n "$P1" ] && [ -f "$P1" ] && { cp "$P1" "$W/in_1.fq.gz"; IN1="$W/in_1.fq.gz"; } || IN1="$P1"
                if [ -n "$P2" ] && [ -f "$P2" ]; then cp "$P2" "$W/in_2.fq.gz"; IN2="$W/in_2.fq.gz"; else IN2=""; fi ;;
  esac
fi
export TASK="$T" TOOL DATASET PARAMS THREADS REP KIND IN1 IN2 TRUTH META IDX OUT W
echo "query $T $TOOL $DATASET $PARAMS t=$THREADS rep=$REP  idx=$IDX_SRC"
set +e
/usr/bin/time -v -o "$OUT/time.txt" \
  bash "$REPO_DIR/scripts/tools/$TOOL.sh" query "$T" "$DATASET" "$PARAMS" "$THREADS" "$W" "$IDX" "$OUT" > "$OUT/tool.log" 2>&1
rc=$?; set -e
[ $rc -eq 0 ] || { echo "tool failed rc=$rc; see $OUT/tool.log"; tail -20 "$OUT/tool.log"; exit $rc; }
[ "$T" = T1 ] && [ -f "$OUT/pairs.tsv" ] && python "$REPO_DIR/scripts/tools/py/names.py" "$OUT/pairs.tsv"   # one naming rule for all T1 tools
inputs=(); [ -n "$IN1" ] && [ -f "$IN1" ] && inputs+=("$IN1"); [ -n "$IN2" ] && [ -f "$IN2" ] && inputs+=("$IN2")
bp_sidecar=""; [ -f "$(dirname "$P1")/reads.bp" ] && bp_sidecar="$(dirname "$P1")/reads.bp"
python "$REPO_DIR/scripts/benchmark/timing_row.py" --task "$T" --tool "$TOOL" --stage query \
  --dataset "$DATASET" --params "$PARAMS" --threads "$THREADS" --rep "$REP" \
  --time-file "$OUT/time.txt" --index "$IDX" ${inputs[@]:+--input "${inputs[@]}"} \
  ${bp_sidecar:+--bp-sidecar "$bp_sidecar"} --sketch "$OUT/sketch" --out "$SB_DATA/$T/timing.csv"
rm -rf "$OUT/sketch"
touch "$OUT/DONE"
