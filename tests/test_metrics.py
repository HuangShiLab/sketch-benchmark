#!/usr/bin/env python
"""End-to-end checks of the metrics layer on tiny synthetic fixtures with known answers.
Run:  python tests/test_metrics.py      (or pytest tests/)
Covers metrics_T1..T4, matched.py, timing_row.py and make_manifests.py."""
import json, math, os, subprocess, sys, tempfile, textwrap
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
B = REPO / "scripts" / "benchmark"
PY = sys.executable

def run(script, *args, env=None):
    e = dict(os.environ); e.update(env or {})
    r = subprocess.run([PY, str(B / script), *map(str, args)], capture_output=True, text=True, env=e)
    assert r.returncode == 0, f"{script} failed:\n{r.stdout}\n{r.stderr}"
    return r.stdout

def w(path, text):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text).lstrip())

def approx(a, b, tol=1e-6): return abs(a - b) <= tol

# ------------------------------------------------------------------ T1
def test_T1(tmp):
    d = tmp / "T1"; sim = d / "sim"; runs = d / "runs"
    w(sim / "truth.tsv", """
    pair_id\tgenome_a\tgenome_b\tpath_a\tpath_b\tseries\tp_nominal\tsubs\tins_bases\tdel_bases\tL_orig\tani_subs_only\tani_gapped\tani_true\tani_true_source
    p1\tGA\tGA_p0.01_subs\t/x/GA.orig.fa.gz\t/x/GA.p0.01.subs.fa.gz\tsubs\t0.01\t0\t0\t0\t0\t0.99\t0.99\t0.990000\tsimulation
    p2\tGB\tGB_p0.02_subs\t/x/GB.orig.fa.gz\t/x/GB.p0.02.subs.fa.gz\tsubs\t0.02\t0\t0\t0\t0\t0.98\t0.98\t0.980000\tsimulation
    p3\tGC\tGC_p0.05_subs\t/x/GC.orig.fa.gz\t/x/GC.p0.05.subs.fa.gz\tsubs\t0.05\t0\t0\t0\t0\t0.95\t0.95\t0.950000\tsimulation
    p4\tGD\tGD_p0.20_subs\t/x/GD.orig.fa.gz\t/x/GD.p0.20.subs.fa.gz\tsubs\t0.20\t0\t0\t0\t0\t0.80\t0.80\t0.800000\tsimulation
    """)
    w(d / "datasets.tsv", f"dataset\tkind\tpath1\tpath2\ttruth\tmeta\nsim\tpairs\t{sim}\t\t{sim/'truth.tsv'}\t\n")
    r = runs / "toolA" / "sim" / "k_21_s_1000" / "t16_r1"
    # p1 exact with CI covering; p2 off by +0.01 with CI not covering; p3 exact no CI; p4 no estimate
    w(r / "pairs.tsv", """
    genome_a\tgenome_b\tj_est\tani_est\tci_lo\tci_hi
    GA_p0.01_subs\tGA\t0.5\t0.990000\t0.985\t0.995
    GB\tGB_p0.02_subs\t0.4\t0.990000\t0.988\t0.992
    GC\tGC_p0.05_subs\t0.2\t0.950000\tNA\tNA
    """)
    (r / "DONE").touch()
    run("metrics_T1.py", "--runs", runs, "--datasets", d / "datasets.tsv", "--out-dir", d)
    M = pd.read_csv(d / "metrics.csv"); S = pd.read_csv(d / "summary.csv")
    assert len(M) == 4 and M.no_estimate.sum() == 1, M
    allrow = S[S.ani_bin == "ALL"].iloc[0]
    assert approx(allrow.bias, (0 + 0.01 + 0) / 3), allrow.bias
    assert approx(allrow.rmse, math.sqrt((0.01 ** 2) / 3)), allrow.rmse
    assert approx(allrow.ci_coverage, 0.5), allrow.ci_coverage          # p1 covered, p2 not; p3 has no CI
    assert approx(allrow.no_estimate, 0.25)
    b95 = S[S.ani_bin == "0.950-0.990"].iloc[0]                      # p2 (0.98, err +0.01) and p3 (0.95, exact)
    assert b95.n == 2 and approx(b95.bias, 0.005) and approx(b95.rmse, math.sqrt(0.01 ** 2 / 2)) and approx(b95.ci_coverage, 0.0)
    b80 = S[S.ani_bin == "0.800-0.850"].iloc[0]; assert b80.n == 1 and b80.n_est == 0 and approx(b80.no_estimate, 1.0)
    b99 = S[S.ani_bin == "0.990-0.995"].iloc[0]; assert b99.n == 1 and approx(b99.bias, 0.0)
    print("T1 ok")

# ------------------------------------------------------------------ T2
def test_T2(tmp):
    d = tmp / "T2"; ds = d / "sim" / "cov_10M_s0"; runs = d / "runs"
    w(ds / "truth.tsv", """
    genome\tpresent\tcov_true\tani_true\tfrac_true
    G1\t1\t0.05\t1.000000\t0.001
    G2\t1\t1\t1.000000\t0.01
    G3\t1\t20\t1.000000\t0.2
    G4\t0\t0\t0.980000\t0
    """)
    dbl = d / "db.list"; w(dbl, "/db/G1_genomic.fna.gz\n/db/G2_genomic.fna.gz\n/db/G3_genomic.fna.gz\n/db/G4_genomic.fna.gz\n/db/G5_genomic.fna.gz\n/db/G6_genomic.fna.gz\n")
    w(d / "datasets.tsv", f"dataset\tkind\tpath1\tpath2\ttruth\tmeta\ncov_10M_s0\tsim\t{ds/'reads_1.fq.gz'}\t{ds/'reads_2.fq.gz'}\t{ds/'truth.tsv'}\t\n")
    r = runs / "toolB" / "cov_10M_s0" / "k_31_scaled_1000" / "t16_r1"
    # calls G2, G3 (TP), G5 (FP); misses G1 (FN); G4, G6 correctly absent
    w(r / "detect.tsv", """
    genome\tscore\tcalled_present\tani_est\tcov_est\tabund_est
    G2\t100\t1\t0.999\t1.2\t0.01
    G3\t2000\t1\t1.0\t18\t0.2
    G5\t5\t1\tNA\t0.1\tNA
    """)
    (r / "DONE").touch()
    run("metrics_T2.py", "--runs", runs, "--datasets", d / "datasets.tsv", "--out-dir", d, "--db-list", dbl)
    S = pd.read_csv(d / "summary.csv").iloc[0]; M = pd.read_csv(d / "metrics.csv")
    assert len(M) == 6, len(M)
    assert (S.tp, S.fp, S.fn, S.n_absent) == (2, 1, 1, 3), S
    assert approx(S.precision, 2 / 3) and approx(S.recall, 2 / 3) and approx(S.f1, 2 / 3)
    assert approx(S.fdr_absent, 1 / 3)
    assert approx(S["recall_cov_0.03-0.075"], 0.0) and approx(S["recall_cov_10-inf"], 1.0)
    assert approx(S.ani_mae, (0.001 + 0.0) / 2)              # G2 0.999 vs 1, G3 exact; G1 has no estimate
    # AUPRC by score: ranking G3(2000) G2(100) G5(5) then zeros -> P@1=1, P@2=1, third positive (G1) at score 0 tied with negatives
    assert 0.6 < S.auprc <= 1.0, S.auprc
    print("T2 ok")

# ------------------------------------------------------------------ T3
def test_T3(tmp):
    d = tmp / "T3"; sim = d / "sim"; runs = d / "runs"
    samples = ["g0s00", "g0s01", "g0s02", "g1s00", "g1s01", "g1s02"]; groups = [0, 0, 0, 1, 1, 1]
    # true BC: 0.1 within group, 0.8 between; jaccard: 0.05 within, 0.6 between
    lines = ["sample_a\tsample_b\tbc_true\tjac_true"]
    import itertools
    for (i, a), (j, b) in itertools.combinations(enumerate(samples), 2):
        same = groups[i] == groups[j]
        lines.append(f"{a}\t{b}\t{0.1 if same else 0.8}\t{0.05 if same else 0.6}")
    # depth replica of g0s00 at distance 0.02 from its base
    lines.append("g0s00\tg0s00d1M\t0.0\t0.0")
    for s in samples[1:]:
        same = groups[0] == groups[samples.index(s)]
        lines.append(f"g0s00d1M\t{s}\t{0.1 if same else 0.8}\t{0.05 if same else 0.6}")
    w(sim / "truth_dist.tsv", "\n".join(lines) + "\n")
    rows = ["dataset\tkind\tpath1\tpath2\ttruth\tmeta"]
    for s, g in zip(samples + ["g0s00d1M"], groups + [0]):
        base = "g0s00" if s == "g0s00d1M" else s
        w(sim / s / "meta.tsv", f"dataset\tgroup\tpairs\tseed\tbase_sample\n{s}\t{g}\t5000000\t42\t{base}\n")
        rows.append(f"{s}\tsim\t{sim/s/'reads_1.fq.gz'}\t{sim/s/'reads_2.fq.gz'}\t{sim/'truth_dist.tsv'}\t{sim/s/'meta.tsv'}")
    w(d / "datasets.tsv", "\n".join(rows) + "\n")
    r = runs / "toolC" / "all" / "k_31_scaled_1000_abund_1" / "t16_r1"
    # estimated distance = 0.5 * bc_true (perfect rank correlation), depth replica at 0.02
    est = ["sample_a\tsample_b\td_est\testimand"]
    for l in lines[1:]:
        a_, b_, bc, _ = l.split("\t"); v = 0.02 if {a_, b_} == {"g0s00", "g0s00d1M"} else 0.5 * float(bc)
        est.append(f"{a_}\t{b_}\t{v:.4f}\tangular")
    w(r / "dist.tsv", "\n".join(est) + "\n"); (r / "DONE").touch()
    run("metrics_T3.py", "--runs", runs, "--datasets", d / "datasets.tsv", "--out-dir", d, env={"HMP_SAMPLE_LIST": ""})
    S = pd.read_csv(d / "summary.csv").iloc[0]
    # p-value is a permutation floor at n=7 with heavy ties; only check it is a valid probability
    assert S.n_samples == 7 and approx(S.mantel_bc, 1.0) and approx(S.mantel_jac, 1.0) and 0 < S.mantel_bc_p <= 0.2, S
    assert approx(S.ari_ward, 1.0), S.ari_ward
    assert approx(S.depth_stability, 0.02), S.depth_stability
    assert S.procrustes_bc < 1e-6, S.procrustes_bc
    print("T3 ok")

# ------------------------------------------------------------------ T4
def test_T4(tmp):
    d = tmp / "T4"; ds = d / "sim" / "base_novaseq_s0"; runs = d / "runs"
    w(ds / "ground_truth_labels.txt", "\n".join([f"host_0_{i}/1\thost" for i in range(6)] + [f"GCF_1.1_{i}/1\tother" for i in range(4)]) + "\n")
    w(ds / "meta.tsv", "dataset\tkind\targ\tseed\tread_len\terr_model\nbase_novaseq_s0\tbase\tnovaseq\t42\t151\tnovaseq\n")
    w(d / "datasets.tsv", f"dataset\tkind\tpath1\tpath2\ttruth\tmeta\nbase_novaseq_s0\tsim\t{ds/'reads_1.fq.gz'}\t{ds/'reads_2.fq.gz'}\t{ds/'ground_truth_labels.txt'}\t{ds/'meta.tsv'}\n"
                          f"HG00438\treal\t/r/1.fq.gz\t/r/2.fq.gz\tlabel:host\t\n")
    r = runs / "deacon" / "base_novaseq_s0" / "k_31_w_15_a_2_r_0.01" / "t16_r1"
    # 5 of 6 host called host (1 FN), 1 of 4 other called host (1 FP); ids arrive with /1 and #suffix noise
    calls = [f"host_0_{i}\thost" for i in range(5)] + ["host_0_5\tother"] + ["GCF_1.1_0/1\thost"] + [f"GCF_1.1_{i}#x\tother" for i in (1, 2, 3)]
    w(r / "calls.tsv", "\n".join(calls) + "\n"); (r / "DONE").touch()
    r2 = runs / "deacon" / "HG00438" / "k_31_w_15_a_2_r_0.01" / "t16_r1"
    w(r2 / "calls.tsv", "\n".join([f"r{i}\thost" for i in range(98)] + ["r98\tother", "r99\tother"]) + "\n"); (r2 / "DONE").touch()
    run("metrics_T4.py", "--runs", runs, "--datasets", d / "datasets.tsv", "--out-dir", d)
    M = pd.read_csv(d / "metrics.csv").set_index("dataset")
    s = M.loc["base_novaseq_s0"]
    assert (s.tp, s.fp, s.tn, s.fn) == (5, 1, 3, 1), s
    assert approx(s.sens, 5 / 6) and approx(s.spec, 3 / 4) and approx(s.residual_host_per_M, 1e6 / 6) and approx(s.microbial_loss, 1 / 4)
    h = M.loc["HG00438"]; assert (h.tp, h.fn) == (98, 2) and approx(h.residual_host_per_M, 20000.0)
    print("T4 ok")

# ------------------------------------------------------------------ timing_row + matched
def test_timing_and_matched(tmp):
    d = tmp / "T1"
    tf = tmp / "time.txt"
    w(tf, """
    \tCommand being timed: "bash x"
    \tUser time (seconds): 12.50
    \tSystem time (seconds): 0.75
    \tElapsed (wall clock) time (h:mm:ss or m:ss): 1:05.25
    \tMaximum resident set size (kbytes): 123456
    """)
    idx = tmp / "idx"; w(idx / "ref.msh", "x" * 1000); bp = tmp / "reads.bp"; w(bp, "15000000000\n")
    for p in ("k_21_s_1000", "k_21_s_4000"):
        run("timing_row.py", "--task", "T1", "--tool", "toolA", "--stage", "query", "--dataset", "sim", "--params", p.replace("_s_", ",s=").replace("k_", "k="),
            "--threads", 16, "--rep", 1, "--time-file", tf, "--index", idx, "--bp-sidecar", bp, "--out", d / "timing.csv")
    T = pd.read_csv(d / "timing.csv"); row = T.iloc[0]
    assert approx(row.wall_s, 65.25) and approx(row.user_s, 12.5) and row.max_rss_kb == 123456 and row.input_bp == 15000000000 and row.index_bytes == 1000, row
    assert row.params == "k=21;s=1000"
    run("matched.py", "--task", "T1", "--data-dir", d, "--target-bytes", 900)
    Mt = pd.read_csv(d / "matched.csv"); assert len(Mt) == 1 and Mt.iloc[0].sketch_bytes == 1000, Mt
    print("timing/matched ok")

# ------------------------------------------------------------------ make_manifests
def test_manifests(tmp):
    sim, real, refs, data = tmp / "mm/sim", tmp / "mm/real", tmp / "mm/refs", tmp / "mm/data"
    (sim / "T2/cov_500k_s0").mkdir(parents=True); (sim / "T2/cov_500k_s0/DONE").touch()
    (sim / "T4/base_novaseq_s0").mkdir(parents=True); (sim / "T4/base_novaseq_s0/DONE").touch()
    (sim / "T3/g0s00").mkdir(parents=True); (sim / "T3/g0s00/DONE").touch(); (sim / "T3/truth_dist.tsv").touch()
    (sim / "T1/sim").mkdir(parents=True); (sim / "T1/sim/truth.tsv").touch()
    refs.mkdir(parents=True); (refs / "gtdb5k.list").touch()
    env = dict(SB_SIM=str(sim), SB_REAL=str(real), SB_REFS=str(refs), SB_DATA=str(data), SB_MINIMAL="1", K_GENOME="21", K_META="31",
               SB_THREADS="16", SB_REPS="1", MASH_S="1000 16000", SOURMASH_SCALED="2000 100", DASHING_LOG2M="12", BINDASH_S="1000",
               SYLPH_C="200", HULK_S="512", DEACON_W="15", DEACON_A="2", DEACON_R="0.01", SYNCMER_S="16", KRAKEN2_CONF="0", BOWTIE2_MAPQ="0",
               RC_PANEL_DIR="", CAMI_MARINE_DIR="TODO", CAMI_STRAIN_DIR="TODO", HMP_SAMPLE_LIST="")
    out = run("make_manifests.py", env=env)
    for t in ("T1", "T2", "T3", "T4"):
        m = pd.read_csv(data / t / "manifest.tsv", sep="\t", header=None); b = pd.read_csv(data / t / "build_manifest.tsv", sep="\t", header=None)
        assert m.shape[1] == 5 and len(m) > 0 and b.shape[1] == 4, (t, m.shape, b.shape)
    m4 = pd.read_csv(data / "T4/manifest.tsv", sep="\t", header=None, names=["tool", "dataset", "params", "threads", "rep"])
    assert set(m4.tool) >= {"deacon", "deacon_syncmer", "kraken2", "bowtie2", "sourmash_read"}, set(m4.tool)
    assert (m4.rep == 1).all() and (m4.threads == 16).all()
    m1 = pd.read_csv(data / "T1/manifest.tsv", sep="\t", header=None, names=["tool", "dataset", "params", "threads", "rep"])
    assert "scale" not in set(m1.dataset)            # minimal mode drops the scale set
    b1 = pd.read_csv(data / "T1/build_manifest.tsv", sep="\t", header=None, names=["tool", "ref", "params", "threads"])
    assert set(b1.ref) == {"sim"}, set(b1.ref)       # T1 builds only for refsets that actually exist
    print("manifests ok:", out.strip().replace("\n", " | "))

if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        test_T1(tmp); test_T2(tmp); test_T3(tmp); test_T4(tmp); test_timing_and_matched(tmp); test_manifests(tmp)
    print("ALL TESTS PASSED")
