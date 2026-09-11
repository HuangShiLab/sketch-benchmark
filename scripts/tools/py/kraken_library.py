#!/usr/bin/env python
"""Concatenate reference genomes into one FASTA for kraken2-build / minimap2.
Headers become `<genome_id>|<contig_n> kraken:taxid|<ncbi_taxid>` (kraken2 needs the taxid
in the header); taxids come from GTDB metadata (ncbi_taxid, falling back to
ncbi_species_taxid). Writes <out>.map (contig -> genome) and, with taxids,
library.map (genome -> taxid) beside <out>. --no-taxid skips the metadata."""
import argparse, csv, gzip, re, sys
from pathlib import Path
ap = argparse.ArgumentParser()
ap.add_argument("--refs", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--metadata"); ap.add_argument("--no-taxid", action="store_true")
a = ap.parse_args()
def gid(p):
    n = Path(p).name
    for e in (".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta"):
        if n.endswith(e): n = n[:-len(e)]; break
    m = re.match(r"^(GC[AF]_\d+\.\d+)", n); return m.group(1) if m else n.split("_genomic")[0]
taxid = {}
if not a.no_taxid:
    if not a.metadata or not Path(a.metadata).exists(): sys.exit("metadata TSV required for taxids")
    for r in csv.DictReader(open(a.metadata), delimiter="\t"):
        acc = r["accession"].split("_", 1)[1] if r["accession"][:3] in ("RS_", "GB_") else r["accession"]
        t = r.get("ncbi_taxid") or r.get("ncbi_species_taxid") or ""
        if t and t not in ("none", "0"): taxid[acc] = t
out = open(a.out, "w"); cmap = open(a.out + ".map", "w"); lib = open(Path(a.out).parent / "library.map", "w")
missing = 0
for line in open(a.refs):
    p = line.strip();
    if not p: continue
    g = gid(p); t = taxid.get(g)
    if not a.no_taxid and not t: missing += 1; continue
    if not a.no_taxid: lib.write(f"{g}\t{t}\n")
    op = gzip.open if p.endswith(".gz") else open
    n = 0
    with op(p, "rt") as fh:
        for l in fh:
            if l.startswith(">"):
                n += 1; hdr = f"{g}|{n}"
                out.write(f">{hdr}" + (f" kraken:taxid|{t}" if t else "") + "\n"); cmap.write(f"{hdr}\t{g}\n")
            else: out.write(l)
if missing: print(f"WARN: {missing} genomes had no NCBI taxid and were skipped", file=sys.stderr)
