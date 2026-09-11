#!/bin/bash
#SBATCH --job-name=sb-06-figures
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
source "${REPO_DIR:?}/scripts/hpc/config.sh"; source "$CONDA_BASE/etc/profile.d/conda.sh"; conda activate "$SB_ENV"
python "$REPO_DIR/scripts/benchmark/figures.py" --data-dir "$SB_DATA" --out-dir "$SB_FIG"
