#!/bin/bash
#SBATCH --job-name=sb-03-build
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
# Stage 03 — one array task = one (tool, refset, params) build from
# data/$SB_TASK/build_manifest.tsv. Timed; index bytes recorded.
# Kraken2 builds need --mem=200G: run_all.sh passes SB_TASK; override memory
# for that array with `--mem` on the sbatch line if the default is too small.
set -euo pipefail
source "${REPO_DIR:?}/scripts/hpc/config.sh"; source "$CONDA_BASE/etc/profile.d/conda.sh"; conda activate "$SB_ENV"
T="${SB_TASK:?SB_TASK (T1..T4) must be exported}"
MAN="$SB_DATA/$T/build_manifest.tsv"
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$MAN"); [ -n "$line" ] || exit 0
IFS=$'\t' read -r TOOL REF PARAMS THREADS <<< "$line"
IDX="$SB_IDX/$T/$TOOL/$REF/${PARAMS//[,=]/_}"; mkdir -p "$IDX"
[ -f "$IDX/DONE" ] && { echo "done: $IDX"; exit 0; }
export PARAMS THREADS TASK="$T" REF IDX
echo "build $T $TOOL $REF $PARAMS t=$THREADS -> $IDX"
/usr/bin/time -v -o "$IDX/time.txt" \
  bash "$REPO_DIR/scripts/tools/$TOOL.sh" build "$T" "$REF" "$PARAMS" "$THREADS" "$IDX" > "$IDX/build.log" 2>&1
python "$REPO_DIR/scripts/benchmark/timing_row.py" --task "$T" --tool "$TOOL" --stage index \
  --dataset "$REF" --params "$PARAMS" --threads "$THREADS" --rep 1 \
  --time-file "$IDX/time.txt" --index "$IDX" --out "$SB_DATA/$T/timing.csv"
touch "$IDX/DONE"
