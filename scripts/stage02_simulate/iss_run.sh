#!/bin/bash
# iss_run.sh <combined.fa.gz> <abund.tsv> <n_reads> <seed> <model> <threads> <out_dir>
# Runs InSilicoSeq, gzips, writes reads.bp and ground_truth_labels.txt (from labels.tsv).
set -euo pipefail
FA="$1"; AB="$2"; N="$3"; SEED="$4"; MODEL="$5"; T="$6"; OUT="$7"
mkdir -p "$OUT"; cd "$OUT"
[ -f DONE ] && exit 0
zcat -f "$FA" > genomes.fa
iss generate --genomes genomes.fa --abundance_file "$AB" --model "$MODEL" \
    --n_reads "$N" --output reads --cpus "$T" --seed "$SEED" > iss.log 2>&1
rm -f genomes.fa
pigz -p "$T" -f reads_R1.fastq; pigz -p "$T" -f reads_R2.fastq
mv reads_R1.fastq.gz reads_1.fq.gz; mv reads_R2.fastq.gz reads_2.fq.gz
seqkit stats -T -j "$T" reads_1.fq.gz reads_2.fq.gz | awk 'NR>1{s+=$5} END{print s}' > reads.bp
# per-read labels (host/other) from the genome-id prefix, same shape as rustyclean-paper's
if [ -f labels.tsv ]; then
  python - <<'PY'
import gzip
lab = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open("labels.tsv") if not l.startswith("genome\t")}
with gzip.open("reads_1.fq.gz", "rt") as fh, open("ground_truth_labels.txt", "w") as out:
    for i, line in enumerate(fh):
        if i % 4 == 0:
            rid = line[1:].split()[0].split("/")[0]
            gid = rid.rsplit("_", 1)[0]
            out.write(f"{rid}\t{lab.get(gid, 'other')}\n")
PY
fi
touch DONE
