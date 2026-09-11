#!/bin/bash
#SBATCH --job-name=sb-02-t2
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-8
# T2 communities. Task list (full): cov_10M_s{0,1,2} cov_50M_s{0,1,2} strain_10M_s{0,1,2}.
# Minimal: cov_500k_s0 strain_500k_s0 (tasks beyond the list exit 0).
set -euo pipefail
source "${REPO_DIR:?}/scripts/hpc/config.sh"; source "$CONDA_BASE/etc/profile.d/conda.sh"; conda activate "$SB_ENV"
S="$REPO_DIR/scripts/stage02_simulate"; T="${SLURM_CPUS_PER_TASK:-4}"
if [ "$SB_MINIMAL" = "1" ]; then
  LIST=("cov 500000 0" "strain 500000 0")
else
  LIST=("cov 10000000 0" "cov 10000000 1" "cov 10000000 2" "cov 50000000 0" "cov 50000000 1" "cov 50000000 2" "strain 10000000 0" "strain 10000000 1" "strain 10000000 2")
fi
i="${SLURM_ARRAY_TASK_ID:-0}"; [ "$i" -lt "${#LIST[@]}" ] || exit 0
read -r MODE PAIRS SEEDI <<< "${LIST[$i]}"
NREADS=$((PAIRS * 2)); SEED=$((SIM_SEED + SEEDI))
NAME="${MODE}_$((PAIRS/1000000))M_s${SEEDI}"; [ "$PAIRS" -lt 1000000 ] && NAME="${MODE}_$((PAIRS/1000))k_s${SEEDI}"
OUT="$SB_SIM/T2/$NAME"; mkdir -p "$OUT"; [ -f "$OUT/DONE" ] && exit 0
# background = 50 genomes from gtdb5k not among the targets (minimal: 10)
NBG=50; [ "$SB_MINIMAL" = "1" ] && NBG=10
grep -v -F -f "$T2_TARGETS_LIST" "$GTDB5K_LIST" | shuf --random-source=<(yes "$SEED") | head -n "$NBG" > "$OUT/background.list"
EXTRA=""
if [ "$MODE" = "strain" ]; then
  [ -f "$SB_SIM/T2/siblings/truth.tsv" ] || { echo "siblings missing; t1.sh must run first" >&2; exit 1; }
  EXTRA="--sibling-dir $SB_SIM/T2/siblings"
fi
python "$S/community.py" t2 --background "$OUT/background.list" --targets "$T2_TARGETS_LIST" \
    --n-reads "$NREADS" --seed "$SEED" --out-prefix "$OUT/community" $EXTRA
cp "$OUT/community.truth.tsv" "$OUT/truth.tsv"; cp "$OUT/community.labels.tsv" "$OUT/labels.tsv"
bash "$S/iss_run.sh" "$OUT/community.fa.gz" "$OUT/community.abund.tsv" "$NREADS" "$SEED" "$ISS_MODEL" "$T" "$OUT"
printf 'dataset\tmode\tpairs\tseed\n%s\t%s\t%s\t%s\n' "$NAME" "$MODE" "$PAIRS" "$SEED" > "$OUT/meta.tsv"
