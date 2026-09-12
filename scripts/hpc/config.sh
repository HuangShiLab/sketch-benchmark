#!/bin/bash
# =============================================================================
# Sketch benchmark — HPC configuration (single source of truth)
# =============================================================================
# Source this in every job. Same conventions as rustyclean-paper/scripts/hpc/
# config.sh: paths on the shared filesystem, overridable by environment.
set -euo pipefail

# ---------------------------------------------------------------------------
# Project layout
# ---------------------------------------------------------------------------
export REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
export SB_DIR="${SB_DIR:-$REPO_DIR}"
export SB_REFS="${SB_REFS:-$SB_DIR/refs}"          # stage 00: reference collections
export SB_REAL="${SB_REAL:-$SB_DIR/data/real}"     # stage 00: real datasets
export SB_SIM="${SB_SIM:-$SB_DIR/data/sim}"        # stage 02: simulated datasets
export SB_IDX="${SB_IDX:-$SB_DIR/idx}"             # stage 03: sketches / indexes
export SB_RUNS="${SB_RUNS:-$SB_DIR/runs}"          # stage 04: raw tool output (gitignored)
export SB_DATA="${SB_DATA:-$SB_DIR/data}"          # stage 05: tidy tables (committed)
export SB_FIG="${SB_FIG:-$SB_DIR/figures}"
export LOG_DIR="${LOG_DIR:-$SB_DIR/logs}"
export SB_BIN="${SB_BIN:-$SB_DIR/env/bin}"         # source-built binaries
export SB_SRC="${SB_SRC:-$SB_DIR/env/src}"
export PATH="$SB_BIN:$PATH"

# Minimal mode (stage 01): tiny grids, 500k-read datasets, one rep, into data/minimal.
export SB_MINIMAL="${SB_MINIMAL:-0}"
if [ "$SB_MINIMAL" = "1" ]; then
  export SB_SIM="$SB_DIR/data/minimal/sim"
  export SB_IDX="$SB_DIR/idx/minimal"
  export SB_RUNS="$SB_DIR/runs/minimal"
  export SB_DATA="$SB_DIR/data/minimal"
  export SB_FIG="$SB_DIR/figures/minimal"
fi

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
export CONDA_BASE="${CONDA_BASE:-/group/aos_shihuang/conda}"
export SB_ENV="${SB_ENV:-sketch-benchmark}"
export SB_PARTITION="${SB_PARTITION:-amd}"
export SB_QOS="${SB_QOS:-normal}"

# ---------------------------------------------------------------------------
# Existing databases on this cluster (rustyclean-paper conventions)
# ---------------------------------------------------------------------------
export DB_ROOT="${DB_ROOT:-/lustre1/g/aos_shihuang/databases}"
export MICROBIAL_GENOME_DIR="${MICROBIAL_GENOME_DIR:-$DB_ROOT/GTDB/GTDBr202/GTDBr202_reference_genome}"
export T2T_FASTA="${T2T_FASTA:-$DB_ROOT/kraken2/kraken16/genomes/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz}"
export HUMAN_GENOME="${HUMAN_GENOME:-$DB_ROOT/kraken2/kraken16/genomes/GCF_000001405.40_GRCh38.p14_genomic.fna.gz}"
export KRAKEN2_DB_T2T_ONLY="${KRAKEN2_DB_T2T_ONLY:-$DB_ROOT/rustyclean_human_t2t_only/kraken2/t2t_only}"
export KRAKEN2_DB_MIXED="${KRAKEN2_DB_MIXED:-$DB_ROOT/kraken2/kraken16}"
export BOWTIE2_INDEX="${BOWTIE2_INDEX:-$DB_ROOT/rustyclean_human_t2t_only/bowtie2/t2t_only}"
export MINIMAP2_INDEX="${MINIMAP2_INDEX:-$DB_ROOT/rustyclean_human_t2t_only/minimap2/t2t_only.mmi}"
export HOSTILE_INDEX="${HOSTILE_INDEX:-$HOME/.local/share/hostile/human-t2t-hla}"
export SOURCE_TAXONOMY="${SOURCE_TAXONOMY:-/lustre1/g/aos_shihuang/tools/kraken2-standard-db/kraken_database/taxonomy}"
# Existing T4 panel from rustyclean-paper (18 datasets with per-read truth).
export RC_PANEL_DIR="${RC_PANEL_DIR:-/lustre1/g/aos_shihuang/rustyclean-paper/scratch/data/enhanced}"

# ---------------------------------------------------------------------------
# Downloads (stage 00). Fill the TODOs before running; stage 00 refuses to
# fetch a set whose URL/accession is still a TODO and records it as skipped.
# ---------------------------------------------------------------------------
export GTDB_METADATA_URL="${GTDB_METADATA_URL:-https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_metadata_r202.tar.gz}"
# deacon's own benchmark sets — Zenodo record ids from deacon/paper/benchmarking.md
export DEACON_ZENODO_CHM13="${DEACON_ZENODO_CHM13:-https://zenodo.org/records/15424966/files}"
export DEACON_ZENODO_ARGOS="${DEACON_ZENODO_ARGOS:-https://zenodo.org/records/15424142/files}"
export DEACON_ZENODO_RSV="${DEACON_ZENODO_RSV:-https://zenodo.org/records/15411280/files}"
export DEACON_PANHUMAN_URL="${DEACON_PANHUMAN_URL:-https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/3/panhuman-1.k31w15.idx}"
# TODO(W1): confirm current locations before stage 00.
export HG002_MAT_URL="${HG002_MAT_URL:-TODO}"          # HPRC HG002 maternal assembly (fa.gz)
export HG002_PAT_URL="${HG002_PAT_URL:-TODO}"          # HPRC HG002 paternal assembly (fa.gz)
export HG00438_R1_URL="${HG00438_R1_URL:-TODO}"        # HPRC HG00438 Illumina 2x150 (ENA fastq.gz)
export HG00438_R2_URL="${HG00438_R2_URL:-TODO}"
export ZYMO_EVEN_SRR="${ZYMO_EVEN_SRR:-TODO}"          # one 2x150 run of ZymoBIOMICS D6300
export ZYMO_LOG_SRR="${ZYMO_LOG_SRR:-TODO}"            # one 2x150 run of ZymoBIOMICS D6310
export ZYMO_REF_URL="${ZYMO_REF_URL:-TODO}"            # Zymo reference genomes bundle
export CAMI_MARINE_DIR="${CAMI_MARINE_DIR:-TODO}"      # staged by hand from the CAMI portal
export CAMI_STRAIN_DIR="${CAMI_STRAIN_DIR:-TODO}"
export HMP_SAMPLE_LIST="${HMP_SAMPLE_LIST:-$REPO_DIR/refs/hmp_samples.tsv}"   # sample_id  site  srr
export DASHING2_URL="${DASHING2_URL:-TODO}"            # github release asset, linux x86_64
export DEACON_UPSTREAM_URL="${DEACON_UPSTREAM_URL:-https://github.com/bede/deacon.git}"
# 2b-RAD tool rows (HuangShiLab); pin a commit before the full run. `TAG=main` builds the tip.
export SYN2BANI_URL="${SYN2BANI_URL:-https://github.com/HuangShiLab/Syn2bANI.git}";       export SYN2BANI_TAG="${SYN2BANI_TAG:-main}"
export FAST2BRAD_URL="${FAST2BRAD_URL:-https://github.com/HuangShiLab/Fast2bRAD-M.git}";  export FAST2BRAD_TAG="${FAST2BRAD_TAG:-main}"
export DEACON_UPSTREAM_TAG="${DEACON_UPSTREAM_TAG:-f2fa660d6dfc46dfa04d2b2dca1adeda8e2cfe27}"   # bede/deacon main @ 0.17.0, the commit env/deacon-syncmer.patch was made against

# ---------------------------------------------------------------------------
# Design constants (rule 3: one k per task)
# ---------------------------------------------------------------------------
export K_GENOME=21            # T1
export K_META=31              # T2 T3 T4
export SIM_SEED="${SIM_SEED:-42}"
export ISS_MODEL="${ISS_MODEL:-novaseq}"   # 151 bp on this cluster's ISS build
export READ_LENGTH=151

# Fixed reference lists — written once by stage 00, never regenerated.
export GTDB5K_LIST="$SB_REFS/gtdb5k.list"
export T1_GENOMES_LIST="$SB_REFS/t1_genomes.list"      # the 50 T1 base genomes
export T1_PAIRS="$SB_REFS/t1_pairs.tsv"                # 2,000 stratified real pairs
export T2_TARGETS_LIST="$SB_REFS/t2_targets.list"      # 20 spike-in targets (in GTDB-5k)
export T3_POOL_LIST="$SB_REFS/t3_pool.list"            # 200 genomes for T3 design
export T4_BACKGROUND_LIST="$SB_REFS/t4_background.list" # 30 microbial genomes for T4 sweeps

# ---------------------------------------------------------------------------
# Sketch-size sweep grids (rule 2). Bytes per 5 Mbp genome in comments.
# ---------------------------------------------------------------------------
export MASH_S="1000 4000 16000 64000"          # 8K 32K 128K 512K
export SOURMASH_SCALED="10000 2000 500 100"    # 4K 20K 80K 400K
export DASHING_LOG2M="10 12 14 16"             # 1K 4K 16K 64K
export BINDASH_S="1000 4000 16000 64000"
export SYLPH_C="1000 200 50"                   # c=200 default
export HULK_S="512 2048"
export I2BRAD_ENZ="BcgI BsaXI CspCI"            # motif-defined sampling: enzymes swept in T1 (density ~ scaled 2000-4000)
export DENSITY_ENZ="BcgI,BsaXI,CspCI,AlfI,CjePI,BplI"   # enzyme-density scan (scripts/analysis/enzyme_density.sh)
export DEACON_W="15"
export DEACON_A="1 2 3"
export DEACON_R="0 0.01 0.05"
export SYNCMER_S="16"                          # s-mer length; l = K_META
export KRAKEN2_CONF="0 0.1 0.3"
export BOWTIE2_MAPQ="0 10"

# Timing (rule 5)
export SB_THREADS="1 16"
export SB_REPS=3
if [ "$SB_MINIMAL" = "1" ]; then
  export MASH_S="1000 16000"; export SOURMASH_SCALED="2000 100"; export DASHING_LOG2M="12"
  export BINDASH_S="1000"; export SYLPH_C="200"; export HULK_S="512"; export I2BRAD_ENZ="BcgI"
  export DEACON_A="2"; export DEACON_R="0.01"; export KRAKEN2_CONF="0"; export BOWTIE2_MAPQ="0"
  export SB_THREADS="16"; export SB_REPS=1
fi

# T4 cautionary rows are expensive in Python; cap the reads they see.
export SB_T4_CAUTION_MAXREADS="${SB_T4_CAUTION_MAXREADS:-5000000}"

mkdir -p "$SB_REFS" "$SB_REAL" "$SB_SIM" "$SB_IDX" "$SB_RUNS" "$SB_DATA" "$SB_FIG" "$LOG_DIR" "$SB_BIN" 2>/dev/null || true
