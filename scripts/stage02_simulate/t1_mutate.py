#!/usr/bin/env python
"""T1 simulated pairs with exact truth.

For every genome in --genomes and every p in --p, write a mutated copy and a
truth row. Two series: subs-only (p_indel = 0) and indel (p_indel = p/10).
Also used by T2-strain to make siblings at a target ANI (--p 0.02 for 98%).

  t1_mutate.py --genomes t1_genomes.list --out-dir sim/T1 --seed 42 \
      --p 0.001 0.005 0.01 0.02 0.05 0.10 0.15 0.20 0.25
writes sim/T1/<gid>.p<p>.<series>.fa.gz and sim/T1/truth.tsv
"""
import argparse, random, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from simlib import collapse, write_fasta, mutate, ani_truth, read_list

ap = argparse.ArgumentParser()
ap.add_argument("--genomes", required=True)
ap.add_argument("--out-dir", required=True)
ap.add_argument("--p", type=float, nargs="+", required=True)
ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--series", nargs="+", default=["subs", "indel"])
ap.add_argument("--truth", default=None, help="truth tsv (default out-dir/truth.tsv)")
a = ap.parse_args()
out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
truth = Path(a.truth) if a.truth else out / "truth.tsv"
new = not truth.exists()
with open(truth, "a") as tf:
    if new:
        tf.write("pair_id\tgenome_a\tgenome_b\tpath_a\tpath_b\tseries\tp_nominal\tsubs\tins_bases\tdel_bases\tL_orig\tani_subs_only\tani_gapped\tani_true\tani_true_source\n")
    for gpath in read_list(a.genomes):
        gid, seq = collapse(gpath)
        # write the collapsed original once so tools compare like with like
        orig = out / f"{gid}.orig.fa.gz"
        if not orig.exists():
            write_fasta([(gid, seq)], orig)
        for p in a.p:
            for series in a.series:
                mp = out / f"{gid}.p{p:g}.{series}.fa.gz"
                if mp.exists():
                    continue
                rng = random.Random(f"{a.seed}:{gid}:{p}:{series}")
                mseq, c = mutate(seq, p, p / 10 if series == "indel" else 0.0, rng)
                write_fasta([(f"{gid}.p{p:g}.{series}", mseq)], mp)
                t = ani_truth(c)
                ani = t["ani_subs_only"] if series == "subs" else t["ani_gapped"]
                tf.write("\t".join(map(str, [
                    f"{gid}:p{p:g}:{series}", gid, f"{gid}.p{p:g}.{series}", orig, mp, series, p,
                    c["subs"], c["ins_bases"], c["del_bases"], c["L_orig"],
                    f"{t['ani_subs_only']:.6f}", f"{t['ani_gapped']:.6f}", f"{ani:.6f}", "simulation"])) + "\n")
print(f"truth: {truth}")
