#!/usr/bin/env python
"""T2 metrics. Joins detect.tsv with per-dataset truth (genome, present, cov_true, ani_true).
Every reference genome in the database is a row; absent genomes not reported by a tool are
true negatives at the tool's default call. Writes
  metrics.csv  per (tool, params, dataset, genome)
  summary.csv  per (tool, params, dataset): precision, recall, f1, auprc, fdr_absent, n_absent,
               recall by coverage bin, ani_mae, cov_log_mae
Score ranking for AUPRC uses the tool's `score` column."""
import argparse, math, os
from pathlib import Path
import numpy as np, pandas as pd
from common import runs, read_tsv, fnum

COV_BINS = [(0, 0.03), (0.03, 0.075), (0.075, 0.3), (0.3, 0.75), (0.75, 3), (3, 10), (10, 1e9)]
def cbin(c):
    for lo, hi in COV_BINS:
        if lo <= c < hi: return f"{lo}-{hi if hi < 1e9 else 'inf'}"
    return "NA"
def auprc(y, s):
    y = np.asarray(y); s = np.asarray(s, dtype=float)
    if y.sum() == 0 or len(y) == 0: return float("nan")
    order = np.argsort(-s); y = y[order]; tp = np.cumsum(y); fp = np.cumsum(1 - y)
    prec = tp / (tp + fp); rec = tp / y.sum()
    ap = 0.0; prev = 0.0
    for p, r in zip(prec, rec):
        ap += p * (r - prev); prev = r
    return ap

ap = argparse.ArgumentParser(); ap.add_argument("--runs", required=True); ap.add_argument("--datasets", required=True); ap.add_argument("--out-dir", required=True)
ap.add_argument("--db-list", help="genome list of the reference DB (absent genomes); default $GTDB5K_LIST")
a = ap.parse_args()
ds = read_tsv(a.datasets)
dblist = a.db_list or os.environ.get("GTDB5K_LIST", "")
import re
def gid(p):
    n = Path(p).name
    for e in (".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta"):
        if n.endswith(e): n = n[:-len(e)]; break
    m = re.match(r"^(GC[AF]_\d+\.\d+)", n); return m.group(1) if m else n.split("_genomic")[0]
db = set(gid(l.strip()) for l in open(dblist) if l.strip()) if dblist and os.path.exists(dblist) else set()
truth = {}
for _, r in ds.iterrows():
    if r["truth"] and os.path.exists(r["truth"]):
        t = read_tsv(r["truth"]); truth[r["dataset"]] = t.set_index("genome")
rows, S = [], []
for tool, dataset, params, threads, rep, d in runs(a.runs, "detect.tsv"):
    if rep != 1 or dataset not in truth: continue
    T = truth[dataset]; E = read_tsv(d / "detect.tsv").drop_duplicates("genome").set_index("genome")
    genomes = set(T.index) | db | set(E.index)
    recs = []
    for g in genomes:
        t = T.loc[g] if g in T.index else None; e = E.loc[g] if g in E.index else None
        present = int(t["present"]) if t is not None else 0
        cov = fnum(t["cov_true"]) if t is not None else 0.0
        ani_t = fnum(t["ani_true"]) if t is not None else float("nan")
        called = int(fnum(e["called_present"]) == 1) if e is not None else 0
        score = fnum(e["score"]) if e is not None else 0.0
        recs.append(dict(tool=tool, params=params, dataset=dataset, genome=g, truth_present=present, truth_cov=cov, truth_ani=ani_t,
                         score=score if not math.isnan(score) else 0.0, called_present=called,
                         ani_est=fnum(e["ani_est"]) if e is not None else float("nan"),
                         cov_est=fnum(e["cov_est"]) if e is not None else float("nan"),
                         abund_est=fnum(e["abund_est"]) if e is not None else float("nan"), cov_bin=cbin(cov) if present else "absent"))
    R = pd.DataFrame(recs); rows.append(R)
    tp = int(((R.truth_present == 1) & (R.called_present == 1)).sum()); fp = int(((R.truth_present == 0) & (R.called_present == 1)).sum())
    fn = int(((R.truth_present == 1) & (R.called_present == 0)).sum()); n_abs = int((R.truth_present == 0).sum())
    prec = tp / (tp + fp) if tp + fp else float("nan"); rec = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * prec * rec / (prec + rec) if prec and rec and not math.isnan(prec) and not math.isnan(rec) and (prec + rec) else float("nan")
    pres = R[R.truth_present == 1]
    ani_mae = (pres.ani_est - pres.truth_ani).abs().mean() if pres.ani_est.notna().any() else float("nan")
    ok = pres[(pres.cov_est > 0) & (pres.truth_cov > 0)]
    cov_mae = (np.log10(ok.cov_est) - np.log10(ok.truth_cov)).abs().mean() if len(ok) else float("nan")
    row = dict(tool=tool, params=params, dataset=dataset, tp=tp, fp=fp, fn=fn, n_absent=n_abs, precision=prec, recall=rec, f1=f1,
               auprc=auprc(R.truth_present.values, R.score.values), fdr_absent=fp / n_abs if n_abs else float("nan"),
               ani_mae=ani_mae, cov_log10_mae=cov_mae)
    for b, g in pres.groupby("cov_bin"):
        row[f"recall_cov_{b}"] = g.called_present.mean()
    S.append(row)
os.makedirs(a.out_dir, exist_ok=True)
(pd.concat(rows) if rows else pd.DataFrame()).to_csv(Path(a.out_dir) / "metrics.csv", index=False)
pd.DataFrame(S).to_csv(Path(a.out_dir) / "summary.csv", index=False)
print(f"T2: {len(S)} runs -> {a.out_dir}")
