#!/usr/bin/env python
"""Write the fixed genome lists from GTDB metadata + the local genome pool.
Refuses to overwrite: these lists are part of the benchmark's identity.

Outputs (in --out-dir):
  gtdb5k.list         5,000 species representatives (T2 reference DB, T3/T4 pools)
  t1_genomes.list     50 genomes stratified by GC and size (T1-sim bases)
  t1_pairs.tsv        2,000 real pairs: 500 each same-species/genus/family/order
  t2_targets.list     20 spike-in targets drawn from gtdb5k
  t3_pool.list        200 genomes for the T3 design
  t4_background.list  30 genomes for T4 sweeps
--minimal scales every count down 10x.
"""
import argparse, csv, glob, os, random, sys
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--metadata", required=True)
ap.add_argument("--genome-dir", required=True)
ap.add_argument("--out-dir", required=True)
ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--minimal", action="store_true")
a = ap.parse_args()
rng = random.Random(a.seed)
out = Path(a.out_dir)
N5K, NT1, NPAIR, NT2, NT3, NT4 = (5000, 50, 500, 20, 200, 30)
if a.minimal:
    N5K, NT1, NPAIR, NT2, NT3, NT4 = (500, 10, 25, 5, 40, 10)

# genome pool: accession (GCF_xxx.x) -> path
pool = {}
for p in glob.glob(os.path.join(a.genome_dir, "*")):
    n = os.path.basename(p)
    for pre in ("GCF_", "GCA_"):
        if n.startswith(pre):
            acc = n.split("_genomic")[0]
            acc = "_".join(acc.split("_")[:2])
            pool[acc] = p
if not pool:
    sys.exit(f"no genomes found under {a.genome_dir}")

# metadata: species reps with taxonomy, GC, size
reps = []
if os.path.exists(a.metadata):
    with open(a.metadata) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            if row.get("gtdb_representative", "t") not in ("t", "True", "TRUE"):
                continue
            acc = row["accession"].split("_", 1)[1] if row["accession"][:3] in ("RS_", "GB_") else row["accession"]
            if acc not in pool:
                continue
            tax = row["gtdb_taxonomy"].split(";")
            reps.append(dict(acc=acc, path=pool[acc], tax=tax,
                             gc=float(row.get("gc_percentage", 50) or 50),
                             size=int(float(row.get("genome_size", 0) or 0))))
else:
    print(f"WARN: metadata {a.metadata} missing; lists will be unstratified", file=sys.stderr)
    reps = [dict(acc=k, path=v, tax=[""] * 7, gc=50.0, size=0) for k, v in pool.items()]
reps.sort(key=lambda r: r["acc"])
rng.shuffle(reps)

def write_list(name, items):
    p = out / name
    if p.exists():
        print(f"keep existing {p}"); return
    with open(p, "w") as fh:
        for it in items:
            fh.write(it + "\n")
    print(f"wrote {p} ({len(items)})")

g5k = reps[:N5K]
write_list("gtdb5k.list", [r["path"] for r in g5k])

# T1 bases: stratify by GC tercile x size tercile within the 5k, round-robin
def terc(vals):
    s = sorted(vals); return s[len(s)//3], s[2*len(s)//3]
gq, sq = terc([r["gc"] for r in g5k]), terc([r["size"] for r in g5k])
cells = {}
for r in g5k:
    key = (sum(r["gc"] > t for t in gq), sum(r["size"] > t for t in sq))
    cells.setdefault(key, []).append(r)
t1 = []
while len(t1) < NT1 and any(cells.values()):
    for k in sorted(cells):
        if cells[k] and len(t1) < NT1:
            t1.append(cells[k].pop())
write_list("t1_genomes.list", [r["path"] for r in t1])

# T1 real pairs: stratified by lowest shared rank (species < genus < family < order)
by_rank = {}
for r in g5k:
    for lvl, name in ((6, "species"), (5, "genus"), (4, "family"), (3, "order")):
        key = (name, ";".join(r["tax"][:lvl + 1]))
        by_rank.setdefault(key, []).append(r)
pairs, seen = [], set()
for rank in ("species", "genus", "family", "order"):
    groups = [v for (n, _), v in by_rank.items() if n == rank and len(v) >= 2]
    rng.shuffle(groups)
    got = 0
    for g in groups:
        if got >= NPAIR: break
        x, y = rng.sample(g, 2)
        # lowest shared rank must be exactly `rank`: reject if they share a lower one
        deeper = {"order": 4, "family": 5, "genus": 6, "species": 7}[rank]
        if deeper <= 6 and x["tax"][:deeper] == y["tax"][:deeper]:
            continue
        key = tuple(sorted((x["acc"], y["acc"])))
        if key in seen: continue
        seen.add(key); got += 1
        pairs.append((f"{rank}_{got:04d}", x["acc"], y["acc"], x["path"], y["path"], rank))
    print(f"  {rank}: {got} pairs", file=sys.stderr)
p = out / "t1_pairs.tsv"
if p.exists():
    print(f"keep existing {p}")
else:
    with open(p, "w") as fh:
        fh.write("pair_id\tgenome_a\tgenome_b\tpath_a\tpath_b\tshared_rank\n")
        for row in pairs: fh.write("\t".join(row) + "\n")
    print(f"wrote {p} ({len(pairs)})")

rest = g5k[NT1:]
write_list("t2_targets.list", [r["path"] for r in rest[:NT2]])
write_list("t3_pool.list", [r["path"] for r in rest[NT2:NT2 + NT3]])
write_list("t4_background.list", [r["path"] for r in rest[NT2 + NT3:NT2 + NT3 + NT4]])
