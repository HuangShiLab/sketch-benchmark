#!/usr/bin/env python
"""Checks of the generic 2b-RAD tag sampler (scripts/tools/py/i2brad.py):
tag lengths match the Fast2bRAD-M table, tags are strand-invariant, containment
ANI recovers a known mutation rate, density/detect/bc subcommands produce sane numbers."""
import gzip, math, os, random, subprocess, sys, tempfile
from pathlib import Path
REPO = Path(__file__).resolve().parents[1]; PYD = REPO / "scripts" / "tools" / "py"
sys.path.insert(0, str(PYD))
import i2brad as ib

FAST2BRAD_TAG_LEN = {"CspCI": 33, "AloI": 27, "BsaXI": 27, "BaeI": 28, "BcgI": 32, "CjeI": 28, "PpiI": 27, "PsrI": 27,
                     "BplI": 27, "FalI": 27, "Bsp24I": 27, "HaeIV": 27, "CjePI": 27, "Hin4I": 27, "AlfI": 32}

def rand_seq(n, seed): rng = random.Random(seed); return "".join(rng.choices("ACGT", k=n))
def mutate(seq, p, seed):
    rng = random.Random(seed); s = list(seq)
    for i in range(len(s)):
        if rng.random() < p: s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)
def write_fa(path, recs):
    with gzip.open(path, "wt") as fh:
        for n, s in recs: fh.write(f">{n}\n{s}\n")
def run(*args):
    r = subprocess.run([sys.executable, str(PYD / "i2brad.py"), *map(str, args)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout

def test_tag_lengths():
    for e, k in FAST2BRAD_TAG_LEN.items(): assert ib.Enzyme(e).k == k, (e, ib.Enzyme(e).k, k)
    print("tag lengths ok (15 enzymes match Fast2bRAD-M)")

def test_planted_and_strand_invariance():
    g = rand_seq(200_000, 1)
    # plant a BcgI site on the top strand at 1000 and one on the bottom strand at 5000
    site = "CGA" + "ACGTAC" + "TGC"; g = g[:1000] + site + g[1012:]; g = g[:5000] + ib.rc(site) + g[5012:]
    e = ib.Enzyme("BcgI"); tags = e.tags(g)
    top = g[1000 - 10: 1012 + 10]; assert min(top, ib.rc(top)) in tags
    bot = g[5000 - 10: 5012 + 10]; assert min(bot, ib.rc(bot)) in tags
    # asymmetric enzyme: a site on the bottom strand must give the same tag as the top-strand site of the rc sequence
    for name in FAST2BRAD_TAG_LEN:
        en = ib.Enzyme(name); h = rand_seq(300_000, hash(name) % 1000)
        a, b = set(en.tags(h)), set(en.tags(ib.rc(h)))
        assert a == b and len(a) > 0, (name, len(a), len(b), len(a ^ b))
        assert all(len(t) == en.k for t in a)
    print("planted sites + strand invariance ok (15 enzymes)")

def test_ani_and_density(tmp):
    A = rand_seq(3_000_000, 7); B = mutate(A, 0.01, 8); C = rand_seq(3_000_000, 9)
    d = tmp / "g"; d.mkdir()
    write_fa(d / "GCF_000001.1.orig.fa.gz", [("A", A)]); write_fa(d / "GCF_000001.1.p0.01.subs.fa.gz", [("B", B)]); write_fa(d / "GCF_000002.1.orig.fa.gz", [("C", C)])
    run("pairs", "--enz", "BcgI", "--threads", 2, "--genomes", d, "-o", tmp / "pairs.tsv")
    rows = [l.split("\t") for l in (tmp / "pairs.tsv").read_text().splitlines()[1:]]
    ab = [r for r in rows if {r[0], r[1]} == {"GCF_000001.1", "GCF_000001.1.p0.01.subs"}]
    assert len(ab) == 1, rows
    ani, lo, hi = float(ab[0][3]), float(ab[0][4]), float(ab[0][5])
    # ~1,450 BcgI tags: sd(ANI) ~ 0.0008, so a single draw sits within +-0.003 of truth and the 95% CI is ~0.002 wide
    assert abs(ani - 0.99) < 0.004 and lo < ani < hi and (hi - lo) < 0.006, (ani, lo, hi)
    assert all({r[0], r[1]} != {"GCF_000001.1", "GCF_000002.1"} or int(r[6]) < 3 for r in rows)   # unrelated genomes share ~0 tags
    run("density", "--enz", "BcgI,BsaXI", "--threads", 2, "--genomes", d, "-o", tmp / "dens")
    summ = [l.split("\t") for l in (tmp / "dens" / "summary.tsv").read_text().splitlines()]
    bcg = dict(zip(summ[0], [r for r in summ[1:] if r[0] == "BcgI"][0]))
    se = float(bcg["median_scaled_equiv"]); assert 1400 < se < 3000, se          # random sequence: ~1 site / 2048 bp
    assert int(bcg["tag_len"]) == 32 and int(bcg["n_genomes"]) == 3
    print(f"ANI ok (est {ani:.4f} CI {lo:.4f}-{hi:.4f}); density ok (BcgI scaled-equivalent {se:.0f})")
    return A, C, d

def test_detect_and_bc(tmp, A, C, d):
    rng = random.Random(11); L = 151
    def reads(seq, n, err=0.0):
        out = []
        for i in range(n):
            s = rng.randrange(0, len(seq) - L); r = seq[s: s + L]
            if rng.random() < 0.5: r = ib.rc(r)
            out.append((f"r{i}", r))
        return out
    def write_fq(path, recs):
        with gzip.open(path, "wt") as fh:
            for n, s in recs: fh.write(f"@{n}\n{s}\n+\n{'I' * len(s)}\n")
    sA = reads(A, 20_000)                          # ~1x of 3 Mb
    write_fq(tmp / "sA.fq.gz", sA)
    run("build-db", "--enz", "BcgI", "--threads", 2, "--genomes", d / "GCF_000001.1.orig.fa.gz", d / "GCF_000002.1.orig.fa.gz", "-o", tmp / "db")
    run("detect", "--enz", "BcgI", "--db", tmp / "db", "--reads", tmp / "sA.fq.gz", "-o", tmp / "detect.tsv", "--read-len", L)
    det = {l.split("\t")[0]: l.split("\t") for l in (tmp / "detect.tsv").read_text().splitlines()[1:]}
    assert "GCF_000001.1" in det and det["GCF_000001.1"][2] == "1", det
    cov = float(det["GCF_000001.1"][4]); assert 0.6 < cov < 1.4, cov     # ~1x
    assert "GCF_000002.1" not in det or det["GCF_000002.1"][2] == "0", det
    # T3: two sketches, Bray-Curtis against a direct computation
    sB = reads(C, 20_000); write_fq(tmp / "sB.fq.gz", sB)
    run("sketch", "--enz", "BcgI", "--reads", tmp / "sA.fq.gz", "-o", tmp / "sA.tags.tsv.gz")
    run("sketch", "--enz", "BcgI", "--reads", tmp / "sB.fq.gz", "-o", tmp / "sB.tags.tsv.gz")
    mix = sA[:10_000] + sB[:10_000]; write_fq(tmp / "sM.fq.gz", mix); run("sketch", "--enz", "BcgI", "--reads", tmp / "sM.fq.gz", "-o", tmp / "sM.tags.tsv.gz")
    (tmp / "samples.tsv").write_text(f"sA\t{tmp/'sA.tags.tsv.gz'}\nsB\t{tmp/'sB.tags.tsv.gz'}\nsM\t{tmp/'sM.tags.tsv.gz'}\n")
    run("bc", "--samples", tmp / "samples.tsv", "--abund", 1, "-o", tmp / "dist.tsv")
    dist = {(l.split("\t")[0], l.split("\t")[1]): float(l.split("\t")[2]) for l in (tmp / "dist.tsv").read_text().splitlines()[1:]}
    a, b = ib.load_sketch(tmp / "sA.tags.tsv.gz"), ib.load_sketch(tmp / "sB.tags.tsv.gz")
    keys = set(a) | set(b); direct = sum(abs(a.get(t, 0) - b.get(t, 0)) for t in keys) / (sum(a.values()) + sum(b.values()))
    assert abs(dist[("sA", "sB")] - direct) < 1e-6 and dist[("sA", "sB")] > 0.99, (dist, direct)
    assert 0.3 < dist[("sA", "sM")] < 0.7, dist                                   # half the reads shared
    print(f"detect ok (cov_est {cov:.2f}); bc ok (A-B {dist[('sA','sB')]:.3f}, A-mix {dist[('sA','sM')]:.3f})")

if __name__ == "__main__":
    test_tag_lengths(); test_planted_and_strand_invariance()
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t); A, C, d = test_ani_and_density(tmp); test_detect_and_bc(tmp, A, C, d)
    print("ALL I2BRAD TESTS PASSED")
