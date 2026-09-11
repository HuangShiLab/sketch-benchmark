#!/usr/bin/env python
"""matched.csv: for every tool, the grid point whose query sketch (or index) bytes are nearest a
target budget — 32 KB for genome-scale tasks (T1), 25 MB for read-scale ones (T2-T4) — joined
to that point's headline accuracy. The "at equal size" table reviewers ask for (rule 2)."""
import argparse, math, os
from pathlib import Path
import pandas as pd
from common import sketch_bytes_lookup
ap = argparse.ArgumentParser(); ap.add_argument("--task", required=True); ap.add_argument("--data-dir", required=True)
ap.add_argument("--target-bytes", type=float, default=None); a = ap.parse_args()
target = a.target_bytes or (32e3 if a.task == "T1" else 25e6)
lookup = sketch_bytes_lookup(a.data_dir)
summ = Path(a.data_dir) / ("summary.csv" if a.task != "T4" else "metrics.csv")
if not summ.exists() or not lookup: print("matched: nothing to do"); raise SystemExit
S = pd.read_csv(summ)
head = {"T1": ("rmse", lambda df: df[df.get("ani_bin", "") == "0.950-0.990"] if "ani_bin" in df else df),
        "T2": ("f1", lambda df: df), "T3": ("mantel_bc", lambda df: df), "T4": ("residual_host_per_M", lambda df: df)}[a.task]
metric, sel = head; S = sel(S)
rows = []
for tool, g in S.groupby("tool"):
    best = None
    for params, gg in g.groupby("params"):
        b = lookup.get((tool, params), (0, 0))[0]
        if b <= 0: continue
        dist = abs(math.log10(b) - math.log10(target))
        if best is None or dist < best[0]: best = (dist, params, b, gg[metric].mean() if metric in gg else float("nan"))
    if best: rows.append(dict(tool=tool, params=best[1], sketch_bytes=best[2], target_bytes=target, metric=metric, value=best[3]))
pd.DataFrame(rows).to_csv(Path(a.data_dir) / "matched.csv", index=False)
print(f"matched: {len(rows)} tools -> {a.data_dir}/matched.csv")
