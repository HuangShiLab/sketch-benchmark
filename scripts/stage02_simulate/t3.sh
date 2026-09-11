#!/bin/bash
#SBATCH --job-name=sb-02-t3
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-79
# T3 design: task 0 writes the design (abundance vectors + truth), then every
# task simulates its sample. Design is idempotent so racing tasks are harmless.
set -euo pipefail
source "${REPO_DIR:?}/scripts/hpc/config.sh"; source "$CONDA_BASE/etc/profile.d/conda.sh"; conda activate "$SB_ENV"
S="$REPO_DIR/scripts/stage02_simulate"; T="${SLURM_CPUS_PER_TASK:-4}"; D="$SB_SIM/T3"; mkdir -p "$D"
if [ "$SB_MINIMAL" = "1" ]; then
  [ -s "$D/samples.tsv" ] || python "$S/t3_design.py" --pool "$T3_POOL_LIST" --out-dir "$D" --per-group 2 --core 10 --shared 5 --depth-replicas 200000,500000 --seed "$SIM_SEED"
  SCALE=25   # 5M -> 200k pairs
else
  [ -s "$D/samples.tsv" ] || python "$S/t3_design.py" --pool "$T3_POOL_LIST" --out-dir "$D" --seed "$SIM_SEED"
  SCALE=1
fi
i="${SLURM_ARRAY_TASK_ID:-0}"
line=$(tail -n +2 "$D/samples.tsv" | sed -n "$((i+1))p"); [ -n "$line" ] || exit 0
read -r SID GRP NREADS_PAIRS BASE <<< "$line"
PAIRS=$((NREADS_PAIRS / SCALE)); NREADS=$((PAIRS * 2)); SEED=$((SIM_SEED + i))
OUT="$D/$SID"; mkdir -p "$OUT"; [ -f "$OUT/DONE" ] && exit 0
python "$S/community.py" t3 --abundance-json "$D/$SID.abund.json" --n-reads "$NREADS" --seed "$SEED" --out-prefix "$OUT/community"
cp "$OUT/community.truth.tsv" "$OUT/truth.tsv"; cp "$OUT/community.labels.tsv" "$OUT/labels.tsv"
bash "$S/iss_run.sh" "$OUT/community.fa.gz" "$OUT/community.abund.tsv" "$NREADS" "$SEED" "$ISS_MODEL" "$T" "$OUT"
printf 'dataset\tgroup\tpairs\tseed\tbase_sample\n%s\t%s\t%s\t%s\t%s\n' "$SID" "$GRP" "$PAIRS" "$SEED" "$BASE" > "$OUT/meta.tsv"
