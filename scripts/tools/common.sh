#!/bin/bash
# Shared helpers for tool wrappers. Sourced by stage04 and by every tools/*.sh.
#
# Wrapper interface (every tools/<tool>.sh implements both, or dies for tasks
# it does not take part in):
#   <tool>.sh build <task> <refset>  <params> <threads> <idx_dir>
#   <tool>.sh query <task> <dataset> <params> <threads> <work_dir> <idx_dir> <out_dir>
# Query environment (exported by stage04): KIND IN1 IN2 TRUTH META IDX OUT W
# Standard query outputs, one per task, written into <out_dir>:
#   T1 pairs.tsv  : genome_a  genome_b  j_est  ani_est  ci_lo  ci_hi     (NA allowed)
#   T2 detect.tsv : genome  score  called_present  ani_est  cov_est  abund_est
#   T3 dist.tsv   : sample_a  sample_b  d_est  estimand
#   T4 calls.tsv  : read_id  call            (call ∈ host|other; read_id normalised)
# A wrapper may leave a query-side sketch in <out_dir>/sketch/ so its bytes are
# recorded; stage04 deletes it after timing_row.py has measured it.

log() { echo "[$(date +%H:%M:%S)] $*" >&2; }
die() { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "$1 not on PATH (conda env $SB_ENV?)"; }

# param KEY [DEFAULT]  — from "$PARAMS" like "k=31,w=15,a=2"
param() {
  local key="$1" def="${2:-}" kv
  IFS=',' read -ra kvs <<< "${PARAMS:-}"
  for kv in "${kvs[@]}"; do [ "${kv%%=*}" = "$key" ] && { echo "${kv#*=}"; return; }; done
  echo "$def"
}

# verify_flags <tool> <flag>...  — fail early (W1) if the installed CLI drifted
verify_flags() {
  local tool="$1"; shift; local help
  help=$("$tool" --help 2>&1 || "$tool" -h 2>&1 || true)
  for f in "$@"; do echo "$help" | grep -q -- "$f" || die "$tool --help does not mention '$f'; CLI drift, fix tools/$(basename "$0")"; done
}

# dataset_row <task> <dataset>  — the line from datasets.tsv
dataset_row() { awk -F'\t' -v d="$2" '$1==d' "$SB_DATA/$1/datasets.tsv" | head -1; }

# index_for <task> <tool> <params>  — where stage 03 put it, or a prebuilt path
index_for() {
  local task="$1" tool="$2" p="$3" ref
  case "$task" in T1) ref="${DATASET_REFSET:-}";; T2) ref=gtdb5k;; T3) ref=samples;; T4) ref=t2t;; esac
  case "$tool" in
    kraken2)          [ "$task" = "T4" ] && { [ "$(PARAMS="$p" param db)" = "mixed" ] && echo "$KRAKEN2_DB_MIXED" || echo "$KRAKEN2_DB_T2T_ONLY"; } || echo "$SB_IDX/$task/kraken2/$ref/${p//[,=]/_}" ;;
    bowtie2)          echo "$(dirname "$BOWTIE2_INDEX")" ;;
    hostile)          echo "$HOSTILE_INDEX" ;;
    minimap2|minimap2_cov) [ "$task" = "T4" ] && echo "$MINIMAP2_INDEX" || echo "$SB_IDX/$task/minimap2_cov/$ref" ;;
    deacon_panhuman)  echo "$SB_IDX/deacon/panhuman-1.k31w15.idx" ;;
    skani|fastani)    echo "" ;;
    *)                if [ "$task" = "T1" ]; then echo "$SB_IDX/$task/$tool/$DATASET/${p//[,=]/_}"; else echo "$SB_IDX/$task/$tool/$ref/${p//[,=]/_}"; fi ;;
  esac
}

# --- read-id helpers (T4) -------------------------------------------------
# normalise: first token, strip /1 /2 and anything after '#'
norm_ids() { awk '{id=$1; sub(/^@/,"",id); sub(/\/[12]$/,"",id); sub(/#.*/,"",id); print id}'; }
fq_ids()   { zcat -f "$1" | awk 'NR%4==1' | norm_ids; }
# write_calls <all_ids> <host_ids> <out>  — every read once; host if in host_ids
write_calls() {
  awk 'NR==FNR{h[$1]=1; next} {print $1"\t"(($1 in h)?"host":"other")}' "$2" "$1" > "$3"
}
# all_ids <out_file>  — ids of IN1 (mate 1 defines the pair)
all_ids() { fq_ids "$IN1" | sort -u > "$1"; }

# --- T1 helpers -------------------------------------------------------------
# pairs_list <out>: pair_id path_a path_b genome_a genome_b, with paths remapped into $IN1
pairs_list() {
  tail -n +2 "$TRUTH" | awk -F'\t' -v d="$IN1" 'BEGIN{OFS="\t"} {
    na=$4; nb=$5; sub(/.*\//,"",na); sub(/.*\//,"",nb); print $1, d"/"na, d"/"nb, $2, $3 }' > "$1"
}
# genome id from a path, matching simlib.genome_id
gid() { python - "$1" <<'PY'
import sys, re
from pathlib import Path
n = Path(sys.argv[1]).name
for e in (".fna.gz",".fa.gz",".fasta.gz",".fna",".fa",".fasta"):
    if n.endswith(e): n = n[:-len(e)]; break
m = re.match(r"^(GC[AF]_\d+\.\d+)", n); print(m.group(1) if m else n.split("_genomic")[0])
PY
}

# refset_list <task> <refset> <out>  — the FASTA paths a stage-03 build sketches
refset_list() {
  case "$1/$2" in
    T1/sim)    ls "$SB_SIM/T1/sim"/*.fa.gz ;;
    T1/real)   cat "$SB_SIM/T1/real/genomes.list" ;;
    T1/scale)  cat "$GTDB5K_LIST" ;;
    T2/gtdb5k) { cat "$GTDB5K_LIST" "$T2_TARGETS_LIST"; [ -d "$SB_REFS/zymo_refs" ] && ls "$SB_REFS/zymo_refs"/*.f*a* ; } | sort -u ;;
    T4/t2t)    echo "$T2T_FASTA" ;;
    T3/samples) awk -F'\t' 'NR>1{print $1"\t"$3"\t"$4}' "$SB_DATA/T3/datasets.tsv" ;;   # sample  r1  r2
    *) die "unknown refset $1/$2" ;;
  esac > "$3"
}
