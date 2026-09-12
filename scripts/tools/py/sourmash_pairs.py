#!/usr/bin/env python
"""Jaccard, containment and ANI (+CI) for listed pairs from a sourmash .zip (Python API).
pairs list: pair_id  path_a  path_b  genome_a  genome_b   (from common.sh pairs_list)
Writes pairs.tsv: genome_a genome_b j_est ani_est ci_lo ci_hi   (ani from max-containment, sourmash's estimator)
VERIFY(W1): MinHash.containment_ani(estimate_ci=True) exists in sourmash>=4.4; fields ani/ani_low/ani_high."""
import argparse, re, sys
from pathlib import Path
import sourmash
ap = argparse.ArgumentParser()
ap.add_argument("--sig", required=True); ap.add_argument("--k", type=int, required=True)
ap.add_argument("--pairs", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
import sys as _sys; from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parent))
from names import genome_name as gid
sigs = {}
for s in sourmash.load_file_as_signatures(a.sig, ksize=a.k):
    for key in {gid(s.filename or ""), gid(s.name or ""), str(s.name)}:
        if key: sigs.setdefault(key, s.minhash.flatten())
miss = 0
with open(a.out, "w") as out:
    out.write("genome_a\tgenome_b\tj_est\tani_est\tci_lo\tci_hi\n")
    for line in open(a.pairs):
        pid, pa, pb, ga, gb = line.rstrip("\n").split("\t")[:5]
        ka, kb = gid(pa), gid(pb)
        if ka not in sigs or kb not in sigs:
            miss += 1; out.write(f"{ga}\t{gb}\tNA\tNA\tNA\tNA\n"); continue
        A, B = sigs[ka], sigs[kb]
        j = A.jaccard(B)
        try:
            r1 = A.containment_ani(B, estimate_ci=True); r2 = B.containment_ani(A, estimate_ci=True)
            r = r1 if (r1.ani or 0) >= (r2.ani or 0) else r2
            ani, lo, hi = r.ani, r.ani_low, r.ani_high
        except Exception:
            ani = lo = hi = None
        f = lambda x: "NA" if x is None else f"{x:.6f}"
        out.write(f"{ga}\t{gb}\t{j:.6g}\t{f(ani)}\t{f(lo)}\t{f(hi)}\n")
if miss: print(f"WARN: {miss} pairs had no signature", file=sys.stderr)
