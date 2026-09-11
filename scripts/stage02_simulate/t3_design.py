#!/usr/bin/env python
"""T3 design: 3 groups x n samples from a 200-genome pool.
Group g has a 60-genome core (group-specific lognormal means) + 20 genomes shared
by all groups; per-sample abundance = lognormal noise (sigma) around the group
vector. Writes one abundance JSON per sample, groups.tsv, and the true
Bray-Curtis / Jaccard matrices over relative read abundance.
"""
import argparse, itertools, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from simlib import read_list

ap = argparse.ArgumentParser()
ap.add_argument("--pool", required=True); ap.add_argument("--out-dir", required=True)
ap.add_argument("--groups", type=int, default=3); ap.add_argument("--per-group", type=int, default=20)
ap.add_argument("--core", type=int, default=60); ap.add_argument("--shared", type=int, default=20)
ap.add_argument("--sigma", type=float, default=0.6); ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--depth-replicas", default="1000000,5000000,20000000")
a = ap.parse_args()
rng = np.random.default_rng(a.seed)
pool = read_list(a.pool); out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
need = a.shared + a.groups * a.core
if len(pool) < need:
    print(f"WARN: pool {len(pool)} < {need}; cores will overlap", file=sys.stderr)
perm = rng.permutation(len(pool)).tolist()
shared = [pool[i] for i in perm[:a.shared]]
cores = [[pool[perm[(a.shared + g * a.core + i) % len(pool)]] for i in range(a.core)] for g in range(a.groups)]
group_mean = {}
for g in range(a.groups):
    genomes = shared + cores[g]
    mu = rng.normal(0, 1.0, size=len(genomes))
    group_mean[g] = dict(zip(genomes, mu))

samples, vectors = [], {}
for g in range(a.groups):
    for k in range(a.per_group):
        sid = f"g{g}s{k:02d}"
        genomes = list(group_mean[g])
        logv = np.array([group_mean[g][x] for x in genomes]) + rng.normal(0, a.sigma, size=len(genomes))
        v = np.exp(logv); v /= v.sum()
        vec = dict(zip(genomes, map(float, v)))
        json.dump(vec, open(out / f"{sid}.abund.json", "w"))
        samples.append((sid, g, 5000000)); vectors[sid] = vec
# depth replicas of the first sample of each group
for g in range(a.groups):
    base = f"g{g}s00"
    for d in a.depth_replicas.split(","):
        sid = f"{base}d{int(d)//1000000}M"
        json.dump(vectors[base], open(out / f"{sid}.abund.json", "w"))
        samples.append((sid, g, int(d))); vectors[sid] = vectors[base]

with open(out / "samples.tsv", "w") as fh:
    fh.write("sample\tgroup\tn_reads\tbase_sample\n")
    for sid, g, n in samples:
        fh.write(f"{sid}\t{g}\t{n}\t{sid.split('d')[0] if 'd' in sid[4:] else sid}\n")

# true distances on relative read abundance (what a k-mer sketch of reads sees)
allg = sorted({x for v in vectors.values() for x in v})
M = np.array([[vectors[s].get(x, 0.0) for x in allg] for s, _, _ in samples])
names = [s for s, _, _ in samples]
with open(out / "truth_dist.tsv", "w") as fh:
    fh.write("sample_a\tsample_b\tbc_true\tjac_true\n")
    for i, j in itertools.combinations(range(len(names)), 2):
        x, y = M[i], M[j]
        bc = np.abs(x - y).sum() / (x + y).sum()
        px, py = x > 0, y > 0
        jac = 1 - (px & py).sum() / max(1, (px | py).sum())
        fh.write(f"{names[i]}\t{names[j]}\t{bc:.6f}\t{jac:.6f}\n")
print(f"{len(samples)} samples ({a.groups} groups x {a.per_group} + depth replicas) -> {out}")
