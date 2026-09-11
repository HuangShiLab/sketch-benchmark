#!/usr/bin/env python
"""Per-read FracMinHash membership (cautionary rows for T4).
--mode sourmash: hash each read at (k, scaled) with sourmash's MinHash, count hashes present in the
  sorted uint64 array of reference hashes (np.searchsorted); host if count >= --abs. Mates pooled.
--mode sylph: write each read as a FASTA "genome", sketch with sylph, query the T2T .syldb;
  host if the read appears in the query output at all. Slow by design; capped."""
import argparse, gzip, os, subprocess, sys
import numpy as np
ap = argparse.ArgumentParser()
ap.add_argument("--mode", choices=["sourmash", "sylph"], required=True)
ap.add_argument("--reads", required=True); ap.add_argument("--reads2")
ap.add_argument("--db", required=True); ap.add_argument("--k", type=int, default=31); ap.add_argument("--scaled", type=int, default=8)
ap.add_argument("--abs", type=int, default=2); ap.add_argument("--max-reads", type=int, default=5_000_000)
ap.add_argument("--threads", type=int, default=8); ap.add_argument("--workdir", default="."); ap.add_argument("--out", required=True)
a = ap.parse_args()

def reads(path):
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h: return
            s = fh.readline().strip(); fh.readline(); fh.readline()
            yield h[1:].split()[0].split("/")[0].split("#")[0], s

if a.mode == "sourmash":
    import sourmash
    ref = np.load(a.db)
    r2 = reads(a.reads2) if a.reads2 else None
    with open(a.out, "w") as out:
        for n, (rid, s1) in enumerate(reads(a.reads)):
            if n >= a.max_reads: break
            mh = sourmash.MinHash(n=0, ksize=a.k, scaled=a.scaled)
            mh.add_sequence(s1, force=True)
            if r2 is not None:
                _, s2 = next(r2); mh.add_sequence(s2, force=True)
            hs = np.fromiter(mh.hashes, dtype=np.uint64, count=len(mh))
            hits = 0
            if len(hs):
                idx = np.searchsorted(ref, hs); idx[idx >= len(ref)] = 0
                hits = int((ref[idx] == hs).sum())
            out.write(f"{rid}\t{'host' if hits >= a.abs else 'other'}\n")
else:
    fa = os.path.join(a.workdir, "reads_as_genomes")
    os.makedirs(fa, exist_ok=True); ids = []
    for n, (rid, s) in enumerate(reads(a.reads)):
        if n >= a.max_reads: break
        with open(os.path.join(fa, f"r{n}.fa"), "w") as f: f.write(f">{rid}\n{s}\n")
        ids.append(rid)
    subprocess.run(f"ls {fa}/*.fa > {a.workdir}/g.list && sylph sketch -k {a.k} -c {a.scaled} -t {a.threads} -l {a.workdir}/g.list -o {a.workdir}/q >/dev/null 2>&1", shell=True, check=True)
    subprocess.run(f"sylph query -t {a.threads} {a.workdir}/q.syldb -d {a.db} -o {a.workdir}/q.tsv >/dev/null 2>&1 || true", shell=True)
    hit = set()
    if os.path.exists(f"{a.workdir}/q.tsv"):
        for l in open(f"{a.workdir}/q.tsv"):
            if l.startswith("Sample_file"): continue
            hit.add(os.path.basename(l.split("\t")[0]))
    with open(a.out, "w") as out:
        for n, rid in enumerate(ids):
            out.write(f"{rid}\t{'host' if f'r{n}.fa' in hit else 'other'}\n")
