#!/usr/bin/env python
"""One rule for turning a sequence file name into the genome name used in truth tables.

  GCF_000123.1_ASM123v1_genomic.fna.gz  -> GCF_000123.1           (real assemblies)
  GCF_000123.1.orig.fa.gz               -> GCF_000123.1           (T1 collapsed original)
  GCF_000123.1.p0.01.subs.fa.gz         -> GCF_000123.1.p0.01.subs (T1 mutant; kept distinct)
  anything_else.fa                      -> anything_else

Used by t1_mutate.py (truth), t1_real_truth.py, the T1 parsers and, via
`python names.py pairs.tsv`, by stage 04 to normalise every T1 wrapper's output."""
import re, sys
from pathlib import Path

_EXT = (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta", ".fastq.gz", ".fq.gz", ".fastq", ".fq")
_MUT = re.compile(r"^(GC[AF]_\d+\.\d+)(\.p[0-9.]+\.[a-z]+)?$")
_GCF = re.compile(r"^(GC[AF]_\d+\.\d+)")

def genome_name(path):
    n = Path(str(path)).name
    for e in _EXT:
        if n.endswith(e): n = n[: -len(e)]; break
    if n.endswith(".orig"): n = n[:-5]
    if _MUT.match(n): return n
    m = _GCF.match(n)
    if m: return m.group(1)
    return n.split("_genomic")[0]

def normalise_pairs(path):
    """Rewrite genome_a/genome_b of a pairs.tsv in place through genome_name()."""
    lines = Path(path).read_text().splitlines()
    if not lines: return
    out = [lines[0]]
    for l in lines[1:]:
        f = l.split("\t")
        if len(f) >= 2: f[0] = genome_name(f[0]); f[1] = genome_name(f[1])
        out.append("\t".join(f))
    Path(path).write_text("\n".join(out) + "\n")

if __name__ == "__main__":
    for p in sys.argv[1:]: normalise_pairs(p)
