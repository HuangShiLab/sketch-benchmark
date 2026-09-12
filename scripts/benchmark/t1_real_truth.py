#!/usr/bin/env python
"""T1-real consensus truth: median of FastANI, skani and (where run) ANIm per pair.
Writes truth.tsv with the same columns as the simulated truth plus one column per reference method,
so the disagreement between references is a reportable noise floor rather than hidden."""
import argparse, csv, glob, os, re, statistics
from pathlib import Path
ap = argparse.ArgumentParser()
ap.add_argument("--pairs", required=True); ap.add_argument("--fastani"); ap.add_argument("--skani"); ap.add_argument("--anim-dir"); ap.add_argument("--out", required=True)
a = ap.parse_args()
import sys as _sys; from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[1] / 'tools' / 'py'))
from names import genome_name as gid
def key(x, y): return tuple(sorted((x, y)))
fast, sk, anim = {}, {}, {}
if a.fastani and os.path.exists(a.fastani):
    for l in open(a.fastani):
        q, r, ani, *_ = l.rstrip("\n").split("\t"); fast[key(gid(q), gid(r))] = float(ani) / 100
if a.skani and os.path.exists(a.skani):
    for r in csv.DictReader(open(a.skani), delimiter="\t"):
        sk[key(gid(r["Ref_file"]), gid(r["Query_file"]))] = float(r["ANI"]) / 100
if a.anim_dir:
    for rep in glob.glob(os.path.join(a.anim_dir, "*.report")):
        pid = Path(rep).stem; val = None
        for l in open(rep):
            if l.startswith("AvgIdentity"):
                val = float(l.split()[1]) / 100; break        # first AvgIdentity line = 1-to-1 alignments
        if val is not None: anim[pid] = val
with open(a.out, "w") as out:
    out.write("pair_id\tgenome_a\tgenome_b\tpath_a\tpath_b\tseries\tp_nominal\tani_true\tani_true_source\tani_fastani\tani_skani\tani_anim\tref_spread\n")
    for r in csv.DictReader(open(a.pairs), delimiter="\t"):
        k = key(r["genome_a"], r["genome_b"])
        vals = {"fastani": fast.get(k), "skani": sk.get(k), "anim": anim.get(r["pair_id"])}
        have = [v for v in vals.values() if v is not None]
        if not have: continue
        med = statistics.median(have); spread = (max(have) - min(have)) if len(have) > 1 else 0.0
        f = lambda v: "NA" if v is None else f"{v:.6f}"
        out.write(f"{r['pair_id']}\t{r['genome_a']}\t{r['genome_b']}\t{r['path_a']}\t{r['path_b']}\treal_{r['shared_rank']}\tNA\t{med:.6f}\tconsensus_{len(have)}\t{f(vals['fastani'])}\t{f(vals['skani'])}\t{f(vals['anim'])}\t{spread:.6f}\n")
