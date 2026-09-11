#!/usr/bin/env python
"""Build one simulated community: combined FASTA, ISS abundance file, truth.

Modes
  t2     --background LIST --targets LIST --target-cov C1,C2,...  (spike-in panel)
         --sibling-dir DIR --sibling-ani A     (strain mode: targets replaced by
                                                mutated siblings; reference absent)
  t3     --pool LIST --abundance-json FILE     (abundance vector precomputed by t3_design.py)
  t4     --background LIST --host FASTA --host-frac F

Common: --n-reads N (mates counted) --seed S --out-prefix P
Writes P.fa.gz, P.abund.tsv, P.truth.tsv, P.labels.tsv (genome_id -> host/other)
"""
import argparse, json, random, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from simlib import (collapse, write_fasta, read_list, lognormal_abundance,
                    write_abundance, coverage_to_fraction, fraction_to_coverage, genome_id)

ap = argparse.ArgumentParser()
ap.add_argument("mode", choices=["t2", "t3", "t4"])
ap.add_argument("--background"); ap.add_argument("--targets"); ap.add_argument("--pool")
ap.add_argument("--target-cov", default="0.01,0.05,0.1,0.5,1,5,20")
ap.add_argument("--sibling-dir"); ap.add_argument("--sibling-p", default="0.05,0.03,0.02,0.01,0.005",
                help="strain mode: substitution rates of the siblings, rotated across targets")
ap.add_argument("--abundance-json")
ap.add_argument("--host"); ap.add_argument("--host-frac", type=float, default=0.5)
ap.add_argument("--n-reads", type=int, required=True)
ap.add_argument("--sigma", type=float, default=1.0)
ap.add_argument("--seed", type=int, required=True)
ap.add_argument("--out-prefix", required=True)
a = ap.parse_args()
rng = np.random.default_rng(a.seed); prng = random.Random(a.seed)
P = Path(a.out_prefix); P.parent.mkdir(parents=True, exist_ok=True)

records, abund, truth, labels = [], [], [], []

def add(path, frac, present=True, ani=1.0, label="other", gid_override=None, sim_from=None):
    gid, seq = collapse(sim_from or path)
    gid = gid_override or gid
    records.append((gid, seq)); abund.append((gid, frac))
    cov = fraction_to_coverage(frac, len(seq), a.n_reads)
    truth.append((gid, int(present), f"{cov:.6g}", f"{ani:.6f}", f"{frac:.6g}"))
    labels.append((gid, label))
    return len(seq)

if a.mode == "t2":
    bg = read_list(a.background); tg = read_list(a.targets)
    covs = [float(c) for c in a.target_cov.split(",")]
    # targets: rotate coverage levels; seed decides the rotation
    off = prng.randrange(len(covs))
    tfrac = []
    for i, t in enumerate(tg):
        cov = covs[(i + off) % len(covs)]
        gid = genome_id(t)
        if a.sibling_dir:
            ps = [float(x) for x in a.sibling_p.split(",")]
            psib = ps[(i + off) % len(ps)]
            sib = Path(a.sibling_dir) / f"{gid}.p{psib:g}.subs.fa.gz"
            if not sib.exists():
                sys.exit(f"no sibling {sib}; run t1_mutate.py on the targets with --p {a.sibling_p.replace(',', ' ')} --series subs")
            # reads come from the sibling; the *reference* genome id is recorded absent
            L = sum(len(s) for _, s in [collapse(sib)])
            frac = coverage_to_fraction(cov, L, a.n_reads)
            add(sib, frac, present=True, ani=1 - psib, label="other", gid_override=f"{gid}_sib")
            truth.append((gid, 0, "0", f"{1 - psib:.6f}", "0"))   # reference absent; sibling ANI recorded
        else:
            L = sum(len(s) for _, s in [collapse(t)])
            frac = coverage_to_fraction(cov, L, a.n_reads)
            add(t, frac, present=True, ani=1.0)
        tfrac.append(frac)
    rest = max(1e-6, 1 - sum(tfrac))
    bfrac = lognormal_abundance(len(bg), rng, a.sigma) * rest
    for b, f in zip(bg, bfrac):
        add(b, float(f))
elif a.mode == "t3":
    vec = json.load(open(a.abundance_json))     # {genome_path: fraction}
    for gpath, f in vec.items():
        if f > 0:
            add(gpath, float(f))
elif a.mode == "t4":
    from simlib import read_fasta
    bg = read_list(a.background)
    # Host stays one record per chromosome/contig (no 3 Gbp string), renamed
    # host_<i>; abundance proportional to length so host coverage is uniform.
    hrecs = [(f"host_{i}", s) for i, (_, s) in enumerate(read_fasta(a.host))]
    htot = sum(len(s) for _, s in hrecs)
    for hid, s in hrecs:
        frac = a.host_frac * len(s) / htot
        records.append((hid, s)); abund.append((hid, frac))
        truth.append((hid, 1, f"{fraction_to_coverage(frac, len(s), a.n_reads):.6g}", "1.000000", f"{frac:.6g}"))
        labels.append((hid, "host"))
    bfrac = lognormal_abundance(len(bg), rng, a.sigma) * (1 - a.host_frac)
    for b, f in zip(bg, bfrac):
        add(b, float(f))

# normalise (ISS requires fractions summing to 1)
tot = sum(f for _, f in abund)
abund = [(g, f / tot) for g, f in abund]
write_fasta(records, f"{P}.fa.gz")
write_abundance(abund, f"{P}.abund.tsv")
with open(f"{P}.truth.tsv", "w") as fh:
    fh.write("genome\tpresent\tcov_true\tani_true\tfrac_true\n")
    for row in truth: fh.write("\t".join(map(str, row)) + "\n")
with open(f"{P}.labels.tsv", "w") as fh:
    fh.write("genome\tlabel\n")
    for g, l in labels: fh.write(f"{g}\t{l}\n")
print(f"{P}: {len(records)} records, n_reads={a.n_reads}")
