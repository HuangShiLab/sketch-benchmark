"""Shared pieces for metrics_T*.py: run discovery, param parsing, safe float."""
from __future__ import annotations
import glob, math, os
from pathlib import Path
import pandas as pd

def runs(runs_dir, task_file):
    """Yield (tool, dataset, params, threads, rep, path) for every finished query run."""
    for done in sorted(glob.glob(os.path.join(runs_dir, "*", "*", "*", "*", "DONE"))):
        d = Path(done).parent
        if not (d / task_file).exists(): continue
        tr = d.name                                  # t16_r1
        threads, rep = int(tr[1:].split("_r")[0]), int(tr.split("_r")[1])
        yield d.parents[2].name, d.parents[1].name, d.parent.name, threads, rep, d

def params_str(p):
    """t16 dir-safe params back to k=31;w=15 form (underscores in values are ambiguous; keep as is)."""
    return p

def fnum(x):
    try:
        v = float(x); return v if math.isfinite(v) else float("nan")
    except (TypeError, ValueError):
        return float("nan")

def read_tsv(path, **kw):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **kw)

def sketch_bytes_lookup(data_dir):
    """tool,params(dir-safe) -> median query sketch_bytes (or index_bytes) from timing.csv."""
    p = Path(data_dir) / "timing.csv"
    if not p.exists(): return {}
    t = pd.read_csv(p)
    t["pdir"] = t["params"].str.replace(";", "_").str.replace("=", "_")
    q = t[t["stage"] == "query"]
    out = {}
    for (tool, pdir), g in q.groupby(["tool", "pdir"]):
        sb = g["sketch_bytes"].median(); ib = g["index_bytes"].median()
        out[(tool, pdir)] = (sb if sb > 0 else ib, ib)
    return out
