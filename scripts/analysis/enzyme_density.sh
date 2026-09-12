#!/bin/bash
#SBATCH --job-name=sb-enzyme-density
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
# Enzyme-density scan: type IIB tags per genome across the GTDB-5k pool (T1 genome list in
# minimal mode) for $DENSITY_ENZ. Writes data/enzyme_density/{per_genome,summary}.tsv — the
# data behind "BcgI ~ scaled 2000" and the GC-dependence claim. Stand-alone; not chained by run_all.sh.
#   sbatch scripts/analysis/enzyme_density.sh          (or: bash ... on a login node with SB_MINIMAL=1)
set -euo pipefail
source "${REPO_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}/scripts/hpc/config.sh"
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then source "$CONDA_BASE/etc/profile.d/conda.sh"; conda activate "$SB_ENV"; fi
LIST="$GTDB5K_LIST"; [ "$SB_MINIMAL" = "1" ] && LIST="$T1_GENOMES_LIST"
[ -s "$LIST" ] || { echo "genome list $LIST missing: run stage 00 first"; exit 1; }
OUT="$SB_DATA/enzyme_density"; mkdir -p "$OUT"
python "$REPO_DIR/scripts/tools/py/i2brad.py" density --enz "$DENSITY_ENZ" --threads "${SLURM_CPUS_PER_TASK:-4}" --genomes "$LIST" -o "$OUT"
column -t "$OUT/summary.tsv"
