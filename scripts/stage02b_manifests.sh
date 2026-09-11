#!/bin/bash
#SBATCH --job-name=sb-02b-manifests
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
# Writes the manifests from what stage 02 produced, then (unless SB_CONTINUE=0)
# submits phase run.
set -euo pipefail
source "${REPO_DIR:?}/scripts/hpc/config.sh"; source "$CONDA_BASE/etc/profile.d/conda.sh"; conda activate "$SB_ENV"
python "$REPO_DIR/scripts/benchmark/make_manifests.py" --tasks "${SB_TASKS:-T1,T2,T3,T4}"
if [ "${SB_CONTINUE:-1}" = "1" ]; then
  args=(--phase run --tasks "${SB_TASKS:-T1,T2,T3,T4}")
  [ "$SB_MINIMAL" = "1" ] && args+=(--minimal)
  nohup bash "$REPO_DIR/scripts/run_all.sh" "${args[@]}" > "$LOG_DIR/submit-run-$(date +%Y%m%d-%H%M%S).log" 2>&1 &
  echo "phase run submitted in background"
fi
