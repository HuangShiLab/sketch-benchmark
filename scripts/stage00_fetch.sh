#!/bin/bash
#SBATCH --job-name=sb-00-fetch
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
# =============================================================================
# Stage 00 — environment manifest, reference lists, downloads.
# Idempotent: every product is guarded by a DONE marker or an existing file.
# A download whose URL/accession is still TODO is recorded as skipped in
# refs/manifest.tsv and the dependent datasets are dropped from the manifests.
# =============================================================================
set -euo pipefail
source "${REPO_DIR:?set REPO_DIR}/scripts/hpc/config.sh"
source "$CONDA_BASE/etc/profile.d/conda.sh"; conda activate "$SB_ENV"
cd "$SB_REFS"
MAN="$SB_REFS/manifest.tsv"
[ -f "$MAN" ] || printf 'item\tstatus\tpath\tbytes\tmd5\tnote\n' > "$MAN"
record() { printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "${4:-}" "${5:-}" "${6:-}" >> "$MAN"; }
fetch() { # fetch <item> <url> <dest>
  local item="$1" url="$2" dest="$3"
  if [ "$url" = "TODO" ]; then record "$item" skipped "$dest" "" "" "URL is TODO in config.sh"; return 0; fi
  if [ -s "$dest" ]; then record "$item" present "$dest" "$(stat -c %s "$dest")" ""; return 0; fi
  mkdir -p "$(dirname "$dest")"
  curl -L --retry 5 --retry-delay 30 -o "$dest.part" "$url" && mv "$dest.part" "$dest"
  record "$item" fetched "$dest" "$(stat -c %s "$dest")" "$(md5sum "$dest" | cut -d' ' -f1)"
}

# --- 1. tool versions --------------------------------------------------------
{
  printf 'tool\tversion\n'
  for t in mash sourmash sylph bindash dashing2 skani fastANI dnadiff kraken2 bowtie2 minimap2 samtools hostile deacon deacon-syncmer hulk simka iss seqkit taxonkit yacht gsearch; do
    v=$( { $t --version 2>&1 || $t version 2>&1 || $t -h 2>&1; } | head -1 | tr '\t' ' ' ) || v="MISSING"
    printf '%s\t%s\n' "$t" "${v:-MISSING}"
  done
} > "$SB_REFS/tool_versions.tsv"
python - <<'PY' >> "$SB_REFS/tool_versions.tsv"
import importlib
for m in ("sourmash","numpy","pandas","scipy"):
    try: print(f"py:{m}\t{importlib.import_module(m).__version__}")
    except Exception as e: print(f"py:{m}\tMISSING")
PY

# --- 2. GTDB metadata and the fixed genome lists ----------------------------
fetch gtdb_metadata "$GTDB_METADATA_URL" "$SB_REFS/bac120_metadata_r202.tar.gz"
if [ -s "$SB_REFS/bac120_metadata_r202.tar.gz" ] && [ ! -s "$SB_REFS/bac120_metadata_r202.tsv" ]; then
  tar -xzf "$SB_REFS/bac120_metadata_r202.tar.gz" -C "$SB_REFS"
  f=$(ls "$SB_REFS"/bac120_metadata_r202*.tsv | head -1); [ "$f" = "$SB_REFS/bac120_metadata_r202.tsv" ] || mv "$f" "$SB_REFS/bac120_metadata_r202.tsv"
fi
# Lists are written once. gtdb_lists.py refuses to overwrite an existing list.
python "$REPO_DIR/scripts/stage02_simulate/gtdb_lists.py" \
  --metadata "$SB_REFS/bac120_metadata_r202.tsv" --genome-dir "$MICROBIAL_GENOME_DIR" \
  --seed "$SIM_SEED" --out-dir "$SB_REFS" $([ "$SB_MINIMAL" = "1" ] && echo --minimal)
for l in gtdb5k.list t1_genomes.list t1_pairs.tsv t2_targets.list t3_pool.list t4_background.list; do
  record "list:$l" present "$SB_REFS/$l" "$(stat -c %s "$SB_REFS/$l")" ""
done

# --- 3. host references ------------------------------------------------------
[ -s "$T2T_FASTA" ] && record t2t present "$T2T_FASTA" || record t2t MISSING "$T2T_FASTA"
[ -s "$HUMAN_GENOME" ] && record grch38 present "$HUMAN_GENOME" || record grch38 MISSING "$HUMAN_GENOME"
fetch hg002_mat "$HG002_MAT_URL" "$SB_REFS/HG002.mat.fa.gz"
fetch hg002_pat "$HG002_PAT_URL" "$SB_REFS/HG002.pat.fa.gz"
fetch deacon_panhuman "$DEACON_PANHUMAN_URL" "$SB_IDX/deacon/panhuman-1.k31w15.idx"

# --- 4. real datasets --------------------------------------------------------
mkdir -p "$SB_REAL"/{T2,T3,T4}
fetch deacon_chm13_r1 "$DEACON_ZENODO_CHM13/chm13v2.r1.fastq.gz" "$SB_REAL/T4/chm13v2.r1.fastq.gz"
fetch deacon_chm13_r2 "$DEACON_ZENODO_CHM13/chm13v2.r2.fastq.gz" "$SB_REAL/T4/chm13v2.r2.fastq.gz"
fetch deacon_argos_r1 "$DEACON_ZENODO_ARGOS/argos988.r1.fastq.zst" "$SB_REAL/T4/argos988.r1.fastq.zst"
fetch deacon_argos_r2 "$DEACON_ZENODO_ARGOS/argos988.r2.fastq.zst" "$SB_REAL/T4/argos988.r2.fastq.zst"
fetch deacon_rsv_r1 "$DEACON_ZENODO_RSV/rsviruses17900.r1.fastq.gz" "$SB_REAL/T4/rsviruses17900.r1.fastq.gz"
fetch deacon_rsv_r2 "$DEACON_ZENODO_RSV/rsviruses17900.r2.fastq.gz" "$SB_REAL/T4/rsviruses17900.r2.fastq.gz"
for s in r1 r2; do f="$SB_REAL/T4/argos988.$s.fastq.zst"; [ -s "$f" ] && [ ! -s "${f%.zst}.gz" ] && { zstdcat "$f" | pigz -p 4 > "${f%.zst}.gz"; }; done
fetch hg00438_r1 "$HG00438_R1_URL" "$SB_REAL/T4/HG00438.r1.fastq.gz"
fetch hg00438_r2 "$HG00438_R2_URL" "$SB_REAL/T4/HG00438.r2.fastq.gz"
# dataset metadata: which real T4 sets are all-host / all-microbe
cat > "$SB_REAL/T4/datasets.tsv" <<TSV
dataset	r1	r2	label
chm13v2	$SB_REAL/T4/chm13v2.r1.fastq.gz	$SB_REAL/T4/chm13v2.r2.fastq.gz	host
argos988	$SB_REAL/T4/argos988.r1.fastq.gz	$SB_REAL/T4/argos988.r2.fastq.gz	other
rsviruses17900	$SB_REAL/T4/rsviruses17900.r1.fastq.gz	$SB_REAL/T4/rsviruses17900.r2.fastq.gz	other
HG00438	$SB_REAL/T4/HG00438.r1.fastq.gz	$SB_REAL/T4/HG00438.r2.fastq.gz	host
TSV

# Zymo mocks via SRA (needs sra-tools; accessions are TODO until chosen)
for pair in "zymo_even:$ZYMO_EVEN_SRR" "zymo_log:$ZYMO_LOG_SRR"; do
  name="${pair%%:*}"; srr="${pair##*:}"
  if [ "$srr" = "TODO" ]; then record "$name" skipped "" "" "" "SRR is TODO"; continue; fi
  if [ ! -s "$SB_REAL/T2/${name}_1.fastq.gz" ]; then
    fasterq-dump --split-files -e 4 -O "$SB_REAL/T2" "$srr" && pigz -p 4 "$SB_REAL/T2/${srr}"_*.fastq
    mv "$SB_REAL/T2/${srr}_1.fastq.gz" "$SB_REAL/T2/${name}_1.fastq.gz"; mv "$SB_REAL/T2/${srr}_2.fastq.gz" "$SB_REAL/T2/${name}_2.fastq.gz"
  fi
  record "$name" present "$SB_REAL/T2/${name}_1.fastq.gz" "$(stat -c %s "$SB_REAL/T2/${name}_1.fastq.gz")" ""
done
fetch zymo_refs "$ZYMO_REF_URL" "$SB_REFS/zymo_refs.tar.gz"

# CAMI II is staged by hand (portal login); record what is there.
for d in "marine:$CAMI_MARINE_DIR" "strain:$CAMI_STRAIN_DIR"; do
  name="cami_${d%%:*}"; dir="${d##*:}"
  if [ "$dir" = "TODO" ] || [ ! -d "$dir" ]; then record "$name" skipped "$dir" "" "" "stage by hand, set CAMI_*_DIR"; else record "$name" present "$dir" "$(du -sb "$dir" | cut -f1)" ""; fi
done

# HMP: sample list (sample_id  site  srr) is curated by hand; fetch + subsample to 5M pairs.
if [ -s "$HMP_SAMPLE_LIST" ] && [ "$(tail -n +2 "$HMP_SAMPLE_LIST" | grep -c '[^[:space:]]')" -gt 0 ]; then   # header-only = not curated yet
  mkdir -p "$SB_REAL/T3"
  tail -n +2 "$HMP_SAMPLE_LIST" | while IFS=$'\t' read -r sid site srr; do
    out="$SB_REAL/T3/${sid}_1.fastq.gz"
    [ -s "$out" ] && continue
    fasterq-dump --split-files -e 4 -O "$SB_REAL/T3/tmp" "$srr"
    seqkit sample -s "$SIM_SEED" -n 5000000 "$SB_REAL/T3/tmp/${srr}_1.fastq" | pigz -p 4 > "$out"
    seqkit sample -s "$SIM_SEED" -n 5000000 "$SB_REAL/T3/tmp/${srr}_2.fastq" | pigz -p 4 > "$SB_REAL/T3/${sid}_2.fastq.gz"
    rm -f "$SB_REAL/T3/tmp/${srr}"_*.fastq
  done
  record hmp present "$SB_REAL/T3" "$(du -sb "$SB_REAL/T3" | cut -f1)" ""
else
  record hmp skipped "$HMP_SAMPLE_LIST" "" "" "curate refs/hmp_samples.tsv"
fi

echo "stage 00 complete; see $MAN"
touch "$SB_REFS/DONE"
