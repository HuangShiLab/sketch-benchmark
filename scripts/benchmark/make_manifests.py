#!/usr/bin/env python
"""Write data/T*/{build_manifest,manifest,datasets}.tsv from the config grids
(read from the environment — run_all.sh sources config.sh first) and from the
datasets stage 02 actually produced.

Query rows: the full parameter grid once at 16 threads (accuracy is
deterministic), plus each tool's default parameters at threads x reps for
timing (rule 5). Minimal mode collapses to one rep at 16 threads.

Columns
  build_manifest.tsv : tool  dataset  params  threads
  manifest.tsv       : tool  dataset  params  threads  rep
  datasets.tsv       : dataset  kind  path1  path2  truth  meta
"""
import argparse, csv, glob, os, sys
from pathlib import Path

E = os.environ
def env(k, d=""): return E.get(k, d)
def grid(k): return env(k).split()
MIN = env("SB_MINIMAL") == "1"
K_G, K_M = env("K_GENOME", "21"), env("K_META", "31")
SIM, REAL, REFS = Path(env("SB_SIM")), Path(env("SB_REAL")), Path(env("SB_REFS"))
DATA = Path(env("SB_DATA"))
THREADS = [int(t) for t in grid("SB_THREADS")] or [16]
REPS = int(env("SB_REPS", "3"))
TMAX = max(THREADS)

# ---------------------------------------------------------------- tools/grids
def T1_tools():
    t = {
        "mash":     ([f"k={K_G},s={s}" for s in grid("MASH_S")] + ([] if MIN else [f"k=16,s=1000", f"k=31,s=1000"]), f"k={K_G},s=1000"),
        "bindash":  ([f"k={K_G},s={s}" for s in grid("BINDASH_S")], f"k={K_G},s=1000"),
        "dashing2": ([f"k={K_G},log2m={m},mode=setsketch" for m in grid("DASHING_LOG2M")] + ([] if MIN else [f"k={K_G},log2m=12,mode=bottomk", f"k={K_G},log2m=12,mode=pminhash"]), f"k={K_G},log2m=12,mode=setsketch"),
        "sourmash": ([f"k={K_G},scaled={s}" for s in grid("SOURMASH_SCALED")], f"k={K_G},scaled=1000"),
        "sylph":    ([f"k={K_G},c={c}" for c in grid("SYLPH_C")], f"k={K_G},c=200"),
        "skani":    (["default"], "default"),
        "fastani":  (["default"], "default"),
        "maxgeomhash": (["b=64"], "b=64"),
        "hypergen": (["default"], "default"),
        # motif-defined sampling (2b-RAD tags): mechanism row + tool row
        "i2brad":   ([f"enz={e}" for e in grid("I2BRAD_ENZ")] + ([] if MIN else ["enz=BcgI+AlfI"]), "enz=BcgI"),
        "syn2bani": (["panel=BcgI+AlfI+AloI+FalI"] + ([] if MIN else ["panel=BcgI"]), "panel=BcgI+AlfI+AloI+FalI"),
    }
    return t
def T2_tools():
    return {
        "sourmash_gather": ([f"k={K_M},scaled={s}" for s in ([1000] if MIN else [1000, 100])], f"k={K_M},scaled=1000"),
        "sylph_profile":   ([f"k={K_M},c={c}" for c in ([200] if MIN else [200, 100])], f"k={K_M},c=200"),
        "yacht":           ([f"k={K_M},scaled=1000,fpr=0.05"], f"k={K_M},scaled=1000,fpr=0.05"),
        "mash_screen":     ([f"k={K_M},s=10000,i=0.9"], f"k={K_M},s=10000,i=0.9"),
        "kraken2":         ([f"conf={c}" for c in ([0] if MIN else [0, 0.1])], "conf=0"),
        "minimap2_cov":    (["default"], "default"),
        "i2brad":          (["enz=BcgI,min_hits=5,min_frac=0.001"] + ([] if MIN else ["enz=BcgI,min_hits=10,min_frac=0.005"]), "enz=BcgI,min_hits=5,min_frac=0.001"),
        "fast2brad_m":     (["enz=BcgI,g=5"], "enz=BcgI,g=5"),
    }
def T3_tools():
    return {
        "sourmash_compare": ([f"k={K_M},scaled=1000,abund=1", f"k={K_M},scaled=1000,abund=0"], f"k={K_M},scaled=1000,abund=1"),
        "hulk":             ([f"k={K_M},s={s}" for s in grid("HULK_S")], f"k={K_M},s=512"),
        "dashing2":         ([f"k={K_M},log2m=12,mode=setsketch", f"k={K_M},log2m=12,mode=pminhash"], f"k={K_M},log2m=12,mode=setsketch"),
        "mash":             ([f"k={K_M},s=10000,m=1", f"k={K_M},s=10000,m=2"], f"k={K_M},s=10000,m=2"),
        "simka":            ([f"k={K_M}"], f"k={K_M}"),
        "sylph_profile_bc": ([f"k={K_M},c=200"], f"k={K_M},c=200"),
        "i2brad":           (["enz=BcgI,abund=1", "enz=BcgI,abund=0"], "enz=BcgI,abund=1"),
        "fast2brad_m":      (["enz=BcgI,g=5"], "enz=BcgI,g=5"),
    }
def T4_tools():
    dg = [f"k={K_M},w={w},a={a},r={r}" for w in grid("DEACON_W") for a in grid("DEACON_A") for r in grid("DEACON_R")]
    sg = [f"k={K_M},s={s},a={a},r={r}" for s in grid("SYNCMER_S") for a in grid("DEACON_A") for r in grid("DEACON_R")]
    return {
        "deacon":          (dg, f"k={K_M},w=15,a=2,r=0.01"),
        "deacon_panhuman": (dg, f"k={K_M},w=15,a=2,r=0.01"),
        "deacon_syncmer":  (sg, f"k={K_M},s=16,a=2,r=0.01"),
        "kraken2":         ([f"conf={c}" for c in grid("KRAKEN2_CONF")] + ([] if MIN else ["db=mixed,conf=0"]), "conf=0"),
        "bowtie2":         ([f"mapq={q}" for q in grid("BOWTIE2_MAPQ")], "mapq=0"),
        "hostile":         (["default"], "default"),
        "minimap2":        (["mapq=0"], "mapq=0"),
        "sylph_read":      ([f"k={K_M},c=200"], f"k={K_M},c=200"),          # cautionary row
        "sourmash_read":   ([f"k={K_M},scaled=8,a=2"], f"k={K_M},scaled=8,a=2"),  # cautionary row
    }
# tools whose index/sketch is built in stage 03 (others use existing indexes or need none)
BUILD = {
    "T1": {"mash", "bindash", "dashing2", "sourmash", "sylph", "maxgeomhash", "hypergen", "i2brad", "syn2bani"},
    "T2": {"sourmash_gather", "sylph_profile", "yacht", "mash_screen", "kraken2", "i2brad", "fast2brad_m"},
    "T3": {"sourmash_compare", "hulk", "dashing2", "mash", "sylph_profile_bc", "i2brad", "fast2brad_m"},
    "T4": {"deacon", "deacon_syncmer", "sylph_read", "sourmash_read"},
}
# build-time dataset (reference set) per task
BUILD_REF = {"T1": ["sim", "real", "scale"], "T2": ["gtdb5k"], "T3": ["samples"], "T4": ["t2t"]}

# ---------------------------------------------------------------- datasets
def datasets(task):
    rows = []  # (dataset, kind, path1, path2, truth, meta)
    if task == "T1":
        if (SIM / "T1/sim/truth.tsv").exists():
            rows.append(("sim", "pairs", str(SIM / "T1/sim"), "", str(SIM / "T1/sim/truth.tsv"), ""))
        if (SIM / "T1/real/truth.tsv").exists():
            rows.append(("real", "pairs", str(SIM / "T1/real/genomes.list"), "", str(SIM / "T1/real/truth.tsv"), str(REFS / "t1_pairs.tsv")))
        if (REFS / "gtdb5k.list").exists() and not MIN:
            rows.append(("scale", "allvsall", str(REFS / "gtdb5k.list"), "", "", ""))
    elif task == "T2":
        for d in sorted(glob.glob(str(SIM / "T2/*/DONE"))):
            dd = Path(d).parent
            if dd.name == "siblings": continue
            rows.append((dd.name, "sim", str(dd / "reads_1.fq.gz"), str(dd / "reads_2.fq.gz"), str(dd / "truth.tsv"), str(dd / "meta.tsv")))
        for name in ("zymo_even", "zymo_log"):
            r1 = REAL / f"T2/{name}_1.fastq.gz"
            if r1.exists():
                rows.append((name, "real", str(r1), str(REAL / f"T2/{name}_2.fastq.gz"), str(REFS / "zymo_truth.tsv"), ""))
        for cam, var in (("cami_marine", "CAMI_MARINE_DIR"), ("cami_strain", "CAMI_STRAIN_DIR")):
            d = env(var)
            if d and d != "TODO" and Path(d).exists():
                for s in sorted(glob.glob(f"{d}/*_1.fastq.gz"))[: (2 if MIN else 999)]:
                    sid = Path(s).name.replace("_1.fastq.gz", "")
                    rows.append((f"{cam}_{sid}", "real", s, s.replace("_1.fastq.gz", "_2.fastq.gz"), f"{d}/{sid}.truth.tsv", ""))
    elif task == "T3":
        for d in sorted(glob.glob(str(SIM / "T3/*/DONE"))):
            dd = Path(d).parent
            rows.append((dd.name, "sim", str(dd / "reads_1.fq.gz"), str(dd / "reads_2.fq.gz"), str(SIM / "T3/truth_dist.tsv"), str(dd / "meta.tsv")))
        for s in sorted(glob.glob(str(REAL / "T3/*_1.fastq.gz"))):
            sid = Path(s).name.replace("_1.fastq.gz", "")
            rows.append((sid, "real", s, s.replace("_1.fastq.gz", "_2.fastq.gz"), env("HMP_SAMPLE_LIST"), ""))
    elif task == "T4":
        for d in sorted(glob.glob(str(SIM / "T4/*/DONE"))):
            dd = Path(d).parent
            rows.append((dd.name, "sim", str(dd / "reads_1.fq.gz"), str(dd / "reads_2.fq.gz"), str(dd / "ground_truth_labels.txt"), str(dd / "meta.tsv")))
        rc = env("RC_PANEL_DIR")
        if rc and Path(rc).exists() and not MIN:
            for gt in sorted(glob.glob(f"{rc}/*/ground_truth_labels.txt")):
                dd = Path(gt).parent
                r1 = dd / "reads_R1.fastq.gz"; r2 = dd / "reads_R2.fastq.gz"; se = dd / "reads.fastq.gz"
                if r1.exists(): rows.append((f"rc_{dd.name}", "panel", str(r1), str(r2), gt, ""))
                elif se.exists(): rows.append((f"rc_{dd.name}", "panel", str(se), "", gt, ""))
        dt = REAL / "T4/datasets.tsv"
        if dt.exists() and not MIN:
            for row in csv.DictReader(open(dt), delimiter="\t"):
                if Path(row["r1"]).exists():
                    rows.append((row["dataset"], "real", row["r1"], row["r2"], f"label:{row['label']}", ""))
    return rows

# ---------------------------------------------------------------- write
ap = argparse.ArgumentParser(); ap.add_argument("--tasks", default="T1,T2,T3,T4"); a = ap.parse_args()
TOOLS = {"T1": T1_tools, "T2": T2_tools, "T3": T3_tools, "T4": T4_tools}
for task in a.tasks.split(","):
    out = DATA / task; out.mkdir(parents=True, exist_ok=True)
    ds = datasets(task)
    with open(out / "datasets.tsv", "w") as fh:
        fh.write("dataset\tkind\tpath1\tpath2\ttruth\tmeta\n")
        for r in ds: fh.write("\t".join(r) + "\n")
    tools = TOOLS[task]()
    nb = nq = 0
    with open(out / "build_manifest.tsv", "w") as fb, open(out / "manifest.tsv", "w") as fq:
        for tool, (params_list, default) in tools.items():
            if tool in BUILD[task]:
                # T1 sketches are built per dataset (the refset *is* the dataset): only emit refsets that exist
                refs = [d[0] for d in ds] if task == "T1" else BUILD_REF[task]
                for ref in refs:
                    for p in params_list:
                        fb.write(f"{tool}\t{ref}\t{p}\t{TMAX}\n"); nb += 1
            qds = (["all"] if ds else []) if task == "T3" else [d[0] for d in ds]   # T3 tools run all-vs-all over every sample
            if task == "T2" and tool == "minimap2_cov":
                qds = [d for d in qds if not d.startswith("cov_50M")]      # gold only on 10M sets
            if task == "T4" and tool in ("sylph_read", "sourmash_read"):
                qds = [d for d in qds if d.startswith(("base_novaseq_s0", "err", "len", "cov_500k"))]  # cautionary rows: sweeps only
            for d in qds:
                for p in params_list:
                    fq.write(f"{tool}\t{d}\t{p}\t{TMAX}\t1\n"); nq += 1
                for t in THREADS:
                    for rep in range(1, REPS + 1):
                        if t == TMAX and rep == 1: continue
                        fq.write(f"{tool}\t{d}\t{default}\t{t}\t{rep}\n"); nq += 1
    print(f"{task}: {len(ds)} datasets, {nb} build tasks, {nq} query tasks -> {out}")
