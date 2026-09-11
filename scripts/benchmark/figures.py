#!/usr/bin/env python
"""Draft figures from data/T*/summary.csv + timing.csv. One PNG per task plus the cross-task
bytes-vs-accuracy panel. Matplotlib only; publication styling is a later pass."""
import argparse, os
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument("--data-dir", required=True); ap.add_argument("--out-dir", required=True); a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True); D = Path(a.data_dir)
def load(t, f):
    p = D / t / f; return pd.read_csv(p) if p.exists() else None
def bytes_of(t):
    tm = load(t, "timing.csv")
    if tm is None: return {}
    tm["pdir"] = tm["params"].str.replace(";", "_").str.replace("=", "_"); q = tm[tm.stage == "query"]
    return {k: (v.sketch_bytes.median() if v.sketch_bytes.median() > 0 else v.index_bytes.median()) for k, v in q.groupby(["tool", "pdir"])}

fig, axes = plt.subplots(2, 2, figsize=(11, 8)); axes = axes.ravel()
spec = [("T1", "summary.csv", "rmse", "RMSE of ANI (95–99% bin)", lambda s: s[s.ani_bin == "0.950-0.990"] if "ani_bin" in s else s),
        ("T2", "summary.csv", "f1", "F1 (detection)", lambda s: s),
        ("T3", "summary.csv", "mantel_bc", "Mantel ρ vs true Bray–Curtis", lambda s: s),
        ("T4", "metrics.csv", "residual_host_per_M", "residual host reads per million", lambda s: s)]
for ax, (t, f, m, lab, sel) in zip(axes, spec):
    S = load(t, f); ax.set_title(f"{t}: {lab}"); ax.set_xlabel("bytes per sketch / index"); ax.set_xscale("log")
    if S is None or m not in S: ax.text(0.5, 0.5, "no data yet", ha="center", transform=ax.transAxes); continue
    S = sel(S); B = bytes_of(t)
    for tool, g in S.groupby("tool"):
        pts = sorted((B.get((tool, p), np.nan), g[g.params == p][m].mean()) for p in g.params.unique())
        pts = [(x, y) for x, y in pts if x == x and y == y]
        if pts: ax.plot(*zip(*pts), "o-", label=tool, ms=4)
    if t == "T4": ax.set_yscale("symlog")
    ax.legend(fontsize=7)
plt.tight_layout(); plt.savefig(Path(a.out_dir) / "fig_cross_task_bytes.png", dpi=150)
for t in ("T1", "T2", "T3", "T4"):
    tm = load(t, "timing.csv")
    if tm is None: continue
    q = tm[(tm.stage == "query") & (tm.input_bp > 0)].copy(); q["mbps"] = q.input_bp / q.wall_s / 1e6
    fig, ax = plt.subplots(figsize=(7, 4)); q.groupby(["tool", "threads"]).mbps.median().unstack().plot.bar(ax=ax)
    ax.set_ylabel("Mbp/s (median)"); ax.set_title(f"{t} query throughput"); plt.tight_layout(); plt.savefig(Path(a.out_dir) / f"fig_{t}_throughput.png", dpi=150); plt.close()
print(f"figures -> {a.out_dir}")
