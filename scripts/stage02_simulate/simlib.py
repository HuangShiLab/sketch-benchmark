"""Shared helpers for the simulation stage. Pure Python + numpy; no Biopython.

Conventions
-----------
* A "genome" is one FASTA file. For simulation every genome is collapsed to a
  single record named by its genome id (filename stem up to the first '.'),
  contigs joined by SPACER_N Ns, so that InSilicoSeq — which treats each FASTA
  record as a separate genome — attributes reads to the genome, and read ids
  carry the genome id as their prefix (`<genome_id>_<n>/1`).
* Abundance files are ISS format: `<record_id>\t<fraction of reads>`.
* Everything random takes an explicit seed.
"""
from __future__ import annotations
import gzip, io, os, random, re
from pathlib import Path
import numpy as np

SPACER_N = 200
READ_LEN = int(os.environ.get("READ_LENGTH", "151"))


def opener(path):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")


def genome_id(path) -> str:
    """GCF_000005825.2_ASM582v2_genomic.fna.gz -> GCF_000005825.2"""
    name = Path(path).name
    for ext in (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    m = re.match(r"^(GC[AF]_\d+\.\d+)", name)
    return m.group(1) if m else name.split("_genomic")[0]


def read_fasta(path):
    """Yield (header, seq) with seq uppercased."""
    hdr, buf = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(buf).upper()
                hdr, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(buf).upper()


def collapse(path) -> tuple[str, str]:
    """One record per genome: contigs joined by N spacers."""
    gid = genome_id(path)
    seqs = [s for _, s in read_fasta(path)]
    return gid, ("N" * SPACER_N).join(seqs)


def write_fasta(records, out_path, width=80):
    op = gzip.open if str(out_path).endswith(".gz") else open
    with op(out_path, "wt") as fh:
        for hdr, seq in records:
            fh.write(f">{hdr}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


def genome_length(path) -> int:
    return sum(len(s) for _, s in read_fasta(path))


def lognormal_abundance(n, rng, sigma=1.0, mean=0.0):
    a = rng.lognormal(mean, sigma, size=n)
    return a / a.sum()


def write_abundance(pairs, out_path):
    with open(out_path, "w") as fh:
        for gid, frac in pairs:
            fh.write(f"{gid}\t{frac:.10g}\n")


def coverage_to_fraction(cov, length, n_reads, read_len=READ_LEN) -> float:
    """Fraction of reads a genome of `length` needs to reach coverage `cov`
    when the run has `n_reads` reads (mates counted separately) of `read_len`."""
    return cov * length / (n_reads * read_len)


def fraction_to_coverage(frac, length, n_reads, read_len=READ_LEN) -> float:
    return frac * n_reads * read_len / length


def mutate(seq: str, p_sub: float, p_indel: float, rng: random.Random, indel_mean=3):
    """Substitute at rate p_sub; insert/delete at rate p_indel each (per base),
    lengths geometric with the given mean. Returns (mutated, counts)."""
    out = []
    subs = ins_b = del_b = 0
    i, n = 0, len(seq)
    bases = "ACGT"
    while i < n:
        c = seq[i]
        r = rng.random()
        if c in bases and r < p_sub:
            out.append(rng.choice([b for b in bases if b != c])); subs += 1; i += 1
        elif p_indel > 0 and r < p_sub + p_indel:
            k = max(1, int(rng.expovariate(1 / indel_mean)))
            del_b += k; i += k                      # deletion: skip k bases
        elif p_indel > 0 and r < p_sub + 2 * p_indel:
            k = max(1, int(rng.expovariate(1 / indel_mean)))
            out.append("".join(rng.choice(bases) for _ in range(k))); ins_b += k
            out.append(c); i += 1                   # insertion before this base
        else:
            out.append(c); i += 1
    return "".join(out), dict(subs=subs, ins_bases=ins_b, del_bases=del_b, L_orig=n)


def ani_truth(counts: dict) -> dict:
    L, s, ins, dl = counts["L_orig"], counts["subs"], counts["ins_bases"], counts["del_bases"]
    return dict(
        ani_subs_only=1 - s / L,
        ani_gapped=1 - (s + ins + dl) / (L + ins),
    )


def read_list(path) -> list[str]:
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip() and not l.startswith("#")]
