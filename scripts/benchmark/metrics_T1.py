#!/usr/bin/env python
"""T1 metrics. Joins each run's pairs.tsv with the dataset truth; writes
  metrics.csv  one row per (tool, params, dataset, pair)
  summary.csv  per (tool, params, dataset, ani_bin): n, bias, rmse, ci_coverage, no_estimate, and overall spearman
ANI bins: 75-80 ... 95-99, 99-99.5, 99.5-100."""
import argparse, math, os
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from common import runs, read_tsv, fnum

BINS = [(0.75, 0.80), (0.80, 0.85), (0.85, 0.90), (0.90, 0.95), (0.95, 0.99), (0.99, 0.995), (0.995, 1.0001)]
def abin(v):
    for lo, hi in BINS:
        if lo <= v < hi: return f"{lo:.3f}-{min(hi,1):.3f}"
    return "<0.75"

ap = argparse.ArgumentParser(); ap.add_argument("--runs", required=True); ap.add_argument("--datasets", required=True); ap.add_argument("--out-dir", required=True)
a = ap.parse_args()
ds = read_tsv(a.datasets); truth = {}
for _, r in ds.iterrows():
    if r["truth"] and os.path.exists(r["truth"]):
        t = read_tsv(r["truth"]); t["k"] = [tuple(sorted((x, y))) for x, y in zip(t["genome_a"], t["genome_b"])]
        truth[r["dataset"]] = t
rows = []
for tool, dataset, params, threads, rep, d in runs(a.runs, "pairs.tsv"):
    if rep != 1 or dataset not in truth: continue
    T = truth[dataset]; est = read_tsv(d / "pairs.tsv")
    est["k"] = [tuple(sorted((x, y))) for x, y in zip(est["genome_a"], est["genome_b"])]
    E = {}                                   # tuple key -> row dict (pandas .loc would misread the tuple as row/col)
    for _, e in est.iterrows(): E.setdefault(e["k"], e)
    for _, t in T.iterrows():
        e = E.get(t["k"])
        ani_est = fnum(e["ani_est"]) if e is not None else float("nan")
        rows.append(dict(tool=tool, params=params, dataset=dataset, pair_id=t["pair_id"], genome_a=t["genome_a"], genome_b=t["genome_b"],
                         ani_true=fnum(t["ani_true"]), ani_true_source=t.get("ani_true_source", ""), series=t.get("series", ""),
                         ani_est=ani_est, ci_lo=fnum(e["ci_lo"]) if e is not None else float("nan"), ci_hi=fnum(e["ci_hi"]) if e is not None else float("nan"),
                         j_est=fnum(e["j_est"]) if e is not None else float("nan"), no_estimate=int(e is None or math.isnan(ani_est))))
M = pd.DataFrame(rows)
os.makedirs(a.out_dir, exist_ok=True)
M.to_csv(Path(a.out_dir) / "metrics.csv", index=False)
if len(M):
    M["err"] = M["ani_est"] - M["ani_true"]; M["bin"] = M["ani_true"].map(abin)
    M["covered"] = ((M["ci_lo"] <= M["ani_true"]) & (M["ani_true"] <= M["ci_hi"])).where(M["ci_lo"].notna())
    S = []
    for (tool, params, dataset, b), g in M.groupby(["tool", "params", "dataset", "bin"]):
        ok = g.dropna(subset=["err"])
        S.append(dict(tool=tool, params=params, dataset=dataset, ani_bin=b, n=len(g), n_est=len(ok),
                      bias=ok["err"].mean() if len(ok) else np.nan, rmse=math.sqrt((ok["err"] ** 2).mean()) if len(ok) else np.nan,
                      mae=ok["err"].abs().mean() if len(ok) else np.nan,
                      ci_coverage=g["covered"].dropna().mean() if g["covered"].notna().any() else np.nan,
                      no_estimate=g["no_estimate"].mean()))
    for (tool, params, dataset), g in M.groupby(["tool", "params", "dataset"]):
        ok = g.dropna(subset=["err"])
        rho = spearmanr(ok["ani_est"], ok["ani_true"]).correlation if len(ok) > 2 else np.nan
        S.append(dict(tool=tool, params=params, dataset=dataset, ani_bin="ALL", n=len(g), n_est=len(ok),
                      bias=ok["err"].mean() if len(ok) else np.nan, rmse=math.sqrt((ok["err"] ** 2).mean()) if len(ok) else np.nan,
                      mae=ok["err"].abs().mean() if len(ok) else np.nan,
                      ci_coverage=g["covered"].dropna().mean() if g["covered"].notna().any() else np.nan,
                      no_estimate=g["no_estimate"].mean(), spearman=rho))
    pd.DataFrame(S).to_csv(Path(a.out_dir) / "summary.csv", index=False)
print(f"T1: {len(M)} pair rows -> {a.out_dir}")
