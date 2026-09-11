#!/bin/bash
#SBATCH --job-name=sb-05-metrics
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
source "${REPO_DIR:?}/scripts/hpc/config.sh"; source "$CONDA_BASE/etc/profile.d/conda.sh"; conda activate "$SB_ENV"
T="${SB_TASK:?}"; t="${T,,}"
python "$REPO_DIR/scripts/benchmark/metrics_$T.py" --runs "$SB_RUNS/$T" --datasets "$SB_DATA/$T/datasets.tsv" --out-dir "$SB_DATA/$T"
python "$REPO_DIR/scripts/benchmark/matched.py" --task "$T" --data-dir "$SB_DATA/$T"
