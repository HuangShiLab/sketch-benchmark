#!/usr/bin/env python
"""Similarity matrix (sourmash compare --csv, or dashing2 --cmpout) -> pairs.tsv (T1) or dist.tsv (T3).
T1: j_est = matrix value; ani_est = 1 - D where D = -(1/k) ln(2J/(1+J))  (Mash's formula).
T3: d_est = 1 - similarity; estimand from --estimand.
VERIFY(W1): dashing2 --cmpout layout; the parser accepts a CSV with a header of names
or a whitespace matrix (names then from --names, in list order)."""
import argparse, csv, math, re, sys
from pathlib import Path
ap = argparse.ArgumentParser()
ap.add_argument("--matrix", required=True); ap.add_argument("--csv", action="store_true")
ap.add_argument("--names"); ap.add_argument("--task", required=True); ap.add_argument("--k", type=int, default=21)
ap.add_argument("--estimand", default="jaccard"); ap.add_argument("--out", required=True)
a = ap.parse_args()

import sys as _sys; from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parent))
from names import genome_name as gid

names, M = None, []
txt = [l for l in open(a.matrix) if l.strip() and not l.startswith("#")]
if a.csv or "," in txt[0]:
    rows = list(csv.reader(txt)); names = rows[0]; M = [[float(x) for x in r[:len(names)]] for r in rows[1:]]
else:
    first = txt[0].split()
    if any(re.search(r"[A-Za-z]", t) for t in first[1:]):      # header row of names
        names = first; M = [[float(x) for x in l.split()[-len(names):]] for l in txt[1:]]
    else:
        M = [[float(x) for x in l.split() if re.match(r"^[-0-9.eE+]+$", x)] for l in txt]
        names = [l.split()[0] for l in txt] if len(M[0]) < len(txt[0].split()) else None
if names is None:
    names = [l.strip().split("\t")[0] for l in open(a.names) if l.strip()]
n = len(names)
if a.task == "T3" and a.names and names and "/" in names[0]:
    # dashing2 names are file paths; datasets.tsv sample ids come from --names (sample\tr1\tr2)
    r1_to_sid = {l.split("\t")[1].strip(): l.split("\t")[0] for l in open(a.names) if "\t" in l}
    names = [r1_to_sid.get(x, gid(x)) for x in names]
else:
    names = [gid(x) for x in names]
with open(a.out, "w") as out:
    if a.task == "T1":
        out.write("genome_a\tgenome_b\tj_est\tani_est\tci_lo\tci_hi\n")
        for i in range(n):
            for j in range(i + 1, n):
                J = M[i][j]
                ani = 1 - (-(1 / a.k) * math.log(2 * J / (1 + J))) if J > 0 else "NA"
                out.write(f"{names[i]}\t{names[j]}\t{J:.6g}\t{ani if ani=='NA' else f'{ani:.6f}'}\tNA\tNA\n")
    else:
        out.write("sample_a\tsample_b\td_est\testimand\n")
        for i in range(n):
            for j in range(i + 1, n):
                out.write(f"{names[i]}\t{names[j]}\t{1 - M[i][j]:.6f}\t{a.estimand}\n")
