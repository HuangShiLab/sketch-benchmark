#!/usr/bin/env python
"""Truncate reads to L bases (gz in, gz out).  truncate_reads.py --len 75 in.fq.gz out.fq.gz"""
import argparse, gzip
ap = argparse.ArgumentParser(); ap.add_argument("--len", type=int, required=True)
ap.add_argument("inp"); ap.add_argument("out"); a = ap.parse_args()
with gzip.open(a.inp, "rt") as fi, gzip.open(a.out, "wt", compresslevel=4) as fo:
    while True:
        h = fi.readline()
        if not h: break
        s = fi.readline().rstrip("\n"); p = fi.readline(); q = fi.readline().rstrip("\n")
        fo.write(h); fo.write(s[:a.len] + "\n"); fo.write(p); fo.write(q[:a.len] + "\n")
