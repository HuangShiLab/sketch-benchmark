#!/bin/bash
#SBATCH --job-name=sb-02-t4
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-14
# T4 sweeps around one 50%-host community (host = GRCh38, 30-genome background).
#   base_<model>_s<seed>   ISS models novaseq (3 seeds), hiseq, miseq
#   err<p>_s0              post-hoc substitutions on base_novaseq_s0
#   len<L>_s0              truncation of base_novaseq_s0
#   div_s<seed>            host = HG002 mat+pat (real haplotype) instead of GRCh38
set -euo pipefail
source "${REPO_DIR:?}/scripts/hpc/config.sh"; source "$CONDA_BASE/etc/profile.d/conda.sh"; conda activate "$SB_ENV"
S="$REPO_DIR/scripts/stage02_simulate"; T="${SLURM_CPUS_PER_TASK:-4}"; D="$SB_SIM/T4"; mkdir -p "$D"
PAIRS=10000000; [ "$SB_MINIMAL" = "1" ] && PAIRS=500000
if [ "$SB_MINIMAL" = "1" ]; then
  LIST=("base novaseq 0" "err 0.01 0" "len 75 0")
else
  LIST=("base novaseq 0" "base novaseq 1" "base novaseq 2" "base hiseq 0" "base miseq 0"
        "err 0.001 0" "err 0.005 0" "err 0.01 0" "err 0.02 0" "len 50 0" "len 75 0" "len 100 0"
        "div hg002 0" "div hg002 1" "div hg002 2")
fi
i="${SLURM_ARRAY_TASK_ID:-0}"; [ "$i" -lt "${#LIST[@]}" ] || exit 0
read -r KIND ARG SEEDI <<< "${LIST[$i]}"; SEED=$((SIM_SEED + SEEDI)); NREADS=$((PAIRS * 2))
simulate_base() { # simulate_base <name> <host.fa.gz|"a b"> <model> <seed>
  local name="$1" host="$2" model="$3" seed="$4" out="$D/$1"
  mkdir -p "$out"; [ -f "$out/DONE" ] && return 0
  python "$S/community.py" t4 --background "$T4_BACKGROUND_LIST" --host "$host" --host-frac 0.5 \
      --n-reads "$NREADS" --seed "$seed" --out-prefix "$out/community"
  cp "$out/community.truth.tsv" "$out/truth.tsv"; cp "$out/community.labels.tsv" "$out/labels.tsv"
  bash "$S/iss_run.sh" "$out/community.fa.gz" "$out/community.abund.tsv" "$NREADS" "$seed" "$model" "$T" "$out"
  printf 'dataset\tkind\targ\tseed\tread_len\terr_model\n%s\t%s\t%s\t%s\t%s\t%s\n' "$name" base "$model" "$seed" "$READ_LENGTH" "$model" > "$out/meta.tsv"
}
derive() { # derive <name> <base> <kind> <arg>  — err or len from an existing base dataset
  local name="$1" base="$D/$2" kind="$3" arg="$4" out="$D/$1"
  # array tasks run concurrently: wait (up to 10 h) for the base dataset
  local waited=0
  while [ ! -f "$base/DONE" ]; do sleep 300; waited=$((waited+300)); [ "$waited" -ge 36000 ] && { echo "base $2 never finished" >&2; exit 1; }; done
  mkdir -p "$out"; [ -f "$out/DONE" ] && return 0
  for m in 1 2; do
    if [ "$kind" = "err" ]; then python "$S/inject_errors.py" --p "$arg" --seed "$((SEED + m))" "$base/reads_$m.fq.gz" "$out/reads_$m.fq.gz"
    else python "$S/truncate_reads.py" --len "$arg" "$base/reads_$m.fq.gz" "$out/reads_$m.fq.gz"; fi
  done
  cp "$base/ground_truth_labels.txt" "$base/labels.tsv" "$base/truth.tsv" "$out/"
  seqkit stats -T -j "$T" "$out"/reads_?.fq.gz | awk 'NR>1{s+=$5} END{print s}' > "$out/reads.bp"
  local rl="$READ_LENGTH"; [ "$kind" = "len" ] && rl="$arg"
  printf 'dataset\tkind\targ\tseed\tread_len\terr_model\n%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$kind" "$arg" "$SEED" "$rl" "novaseq+$kind$arg" > "$out/meta.tsv"
  touch "$out/DONE"
}
case "$KIND" in
  base) simulate_base "base_${ARG}_s${SEEDI}" "$HUMAN_GENOME" "$ARG" "$SEED" ;;
  err)  derive "err${ARG}_s${SEEDI}" "base_novaseq_s0" err "$ARG" ;;
  len)  derive "len${ARG}_s${SEEDI}" "base_novaseq_s0" len "$ARG" ;;
  div)  if [ -s "$SB_REFS/HG002.mat.fa.gz" ] && [ -s "$SB_REFS/HG002.pat.fa.gz" ]; then
          HG="$D/hg002.fa.gz"; [ -s "$HG" ] || cat "$SB_REFS/HG002.mat.fa.gz" "$SB_REFS/HG002.pat.fa.gz" > "$HG"
          simulate_base "div_s${SEEDI}" "$HG" novaseq "$SEED"
        else echo "HG002 assemblies not fetched (config TODO); skipping div_s${SEEDI}"; fi ;;
esac
