#!/usr/bin/env python
"""Kraken2 report -> detect.tsv per reference genome: score = reads assigned at or below the genome's taxid.
Several genomes may share a species taxid; they all receive that taxid's count (a known limitation of
taxid-level classification, reported as such)."""
import argparse, collections
ap = argparse.ArgumentParser(); ap.add_argument("--report", required=True); ap.add_argument("--map", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
clade = {}
for l in open(a.report):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 6: clade[f[4].strip()] = int(f[1])       # clade-level read count
g2t = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(a.map) if "\t" in l}
with open(a.out, "w") as out:
    out.write("genome\tscore\tcalled_present\tani_est\tcov_est\tabund_est\n")
    for g, t in g2t.items():
        n = clade.get(t, 0); out.write(f"{g}\t{n}\t{int(n >= 10)}\tNA\tNA\tNA\n")
