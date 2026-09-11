#!/bin/bash
#SBATCH --job-name=sb-02-t1
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
# T1-sim: mutated genomes with exact truth; T2-strain siblings; T1-real pair truth (skani/FastANI/ANIm).
set -euo pipefail
source "${REPO_DIR:?}/scripts/hpc/config.sh"; source "$CONDA_BASE/etc/profile.d/conda.sh"; conda activate "$SB_ENV"
S="$REPO_DIR/scripts/stage02_simulate"; OUT="$SB_SIM/T1"; mkdir -p "$OUT"
P="0.001 0.005 0.01 0.02 0.05 0.10 0.15 0.20 0.25"; [ "$SB_MINIMAL" = "1" ] && P="0.01 0.05 0.20"
python "$S/t1_mutate.py" --genomes "$T1_GENOMES_LIST" --out-dir "$OUT/sim" --seed "$SIM_SEED" --p $P
# siblings for T2-strain (subs only) at the strain levels
python "$S/t1_mutate.py" --genomes "$T2_TARGETS_LIST" --out-dir "$SB_SIM/T2/siblings" --seed "$SIM_SEED" \
    --p 0.05 0.03 0.02 0.01 0.005 --series subs --truth "$SB_SIM/T2/siblings/truth.tsv"
# T1-real: consensus truth. FastANI and skani on all pairs; ANIm on the first 200 (or 20 minimal).
if [ -s "$T1_PAIRS" ]; then
  mkdir -p "$OUT/real"; cd "$OUT/real"
  tail -n +2 "$T1_PAIRS" | cut -f4 | sort -u > q.list; tail -n +2 "$T1_PAIRS" | cut -f5 | sort -u > r.list
  [ -s fastani.tsv ] || fastANI --ql q.list --rl r.list -t "$SLURM_CPUS_PER_TASK" -o fastani.tsv >/dev/null 2>&1 || true
  cat q.list r.list | sort -u > genomes.list
  [ -s skani.tsv ] || skani triangle -l genomes.list -t "$SLURM_CPUS_PER_TASK" --sparse -o skani.tsv >/dev/null 2>&1 || true
  NANIM=200; [ "$SB_MINIMAL" = "1" ] && NANIM=20
  mkdir -p anim
  tail -n +2 "$T1_PAIRS" | head -n "$NANIM" | while IFS=$'\t' read -r pid ga gb pa pb rank; do
    [ -s "anim/$pid.report" ] && continue
    zcat -f "$pa" > anim/a.fa; zcat -f "$pb" > anim/b.fa
    dnadiff -p "anim/$pid" anim/a.fa anim/b.fa >/dev/null 2>&1 || true
    rm -f anim/a.fa anim/b.fa anim/$pid.{1coords,1delta,mcoords,mdelta,qdiff,rdiff,snps,unqry,unref,delta}
  done
  python "$REPO_DIR/scripts/benchmark/t1_real_truth.py" --pairs "$T1_PAIRS" --fastani fastani.tsv --skani skani.tsv --anim-dir anim --out "$OUT/real/truth.tsv"
fi
touch "$OUT/DONE"
