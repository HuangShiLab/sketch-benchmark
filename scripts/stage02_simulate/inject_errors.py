#!/usr/bin/env python
"""Post-hoc substitution errors at rate p on FASTQ (gz in, gz out); qualities kept.
  inject_errors.py --p 0.01 --seed 42 in.fq.gz out.fq.gz"""
import argparse, gzip, random
ap = argparse.ArgumentParser()
ap.add_argument("--p", type=float, required=True); ap.add_argument("--seed", type=int, default=42)
ap.add_argument("inp"); ap.add_argument("out")
a = ap.parse_args(); rng = random.Random(a.seed)
alt = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}
with gzip.open(a.inp, "rt") as fi, gzip.open(a.out, "wt", compresslevel=4) as fo:
    while True:
        h = fi.readline()
        if not h: break
        s = fi.readline().rstrip("\n"); p = fi.readline(); q = fi.readline()
        s = "".join(rng.choice(alt[c]) if (c in alt and rng.random() < a.p) else c for c in s)
        fo.write(h); fo.write(s + "\n"); fo.write(p); fo.write(q)
