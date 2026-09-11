#!/usr/bin/env python
"""T4 metrics. calls.tsv (read_id, call) vs per-read truth.
Truth sources: ground_truth_labels.txt (read_id<TAB>host|microbe|other), or `label:host` /
`label:other` for real all-host / all-microbial sets. Read ids are normalised identically on both sides.
  metrics.csv  per (tool, params, dataset, rep): tp fp tn fn sens spec f1 mcc residual_host_per_M microbial_loss
Only rep 1 is scored (depletion is deterministic); timing lives in timing.csv."""
import argparse, math, os, re
from pathlib import Path
import pandas as pd
from common import runs, read_tsv

def norm(rid): return rid.split()[0].lstrip("@").split("/")[0].split("#")[0]
def load_truth(spec):
    if spec.startswith("label:"): return None, spec.split(":")[1]
    host, other = set(), set()
    for l in open(spec):
        rid, lab = l.rstrip("\n").split("\t", 1)
        (host if lab.lower() == "host" else other).add(norm(rid))
    return (host, other), None

ap = argparse.ArgumentParser(); ap.add_argument("--runs", required=True); ap.add_argument("--datasets", required=True); ap.add_argument("--out-dir", required=True)
a = ap.parse_args()
ds = read_tsv(a.datasets).set_index("dataset")
meta = {}
for dname, r in ds.iterrows():
    if r["meta"] and os.path.exists(r["meta"]): meta[dname] = read_tsv(r["meta"]).iloc[0].to_dict()
cache, S = {}, []
for tool, dataset, params, threads, rep, d in runs(a.runs, "calls.tsv"):
    if rep != 1 or dataset not in ds.index: continue
    spec = ds.loc[dataset, "truth"]
    if spec not in cache: cache[spec] = load_truth(spec)
    sets, uniform = cache[spec]
    calls = {}
    for l in open(d / "calls.tsv"):
        rid, c = l.rstrip("\n").split("\t"); calls[norm(rid)] = c
    if uniform:
        n_host = sum(1 for c in calls.values() if c == "host"); n = len(calls)
        if uniform == "host": tp, fn, fp, tn = n_host, n - n_host, 0, 0
        else: tp, fn, fp, tn = 0, 0, n_host, n - n_host
    else:
        host, other = sets
        tp = sum(1 for r in host if calls.get(r) == "host"); fn = len(host) - tp
        fp = sum(1 for r in other if calls.get(r) == "host"); tn = len(other) - fp
        missing = len(host) + len(other) - len(calls)
    sens = tp / (tp + fn) if tp + fn else float("nan"); spec_ = tn / (tn + fp) if tn + fp else float("nan")
    prec = tp / (tp + fp) if tp + fp else float("nan")
    f1 = 2 * prec * sens / (prec + sens) if not (math.isnan(prec) or math.isnan(sens)) and prec + sens else float("nan")
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / den if den else float("nan")
    m = meta.get(dataset, {})
    S.append(dict(tool=tool, params=params, dataset=dataset, rep=rep, host_pct=m.get("host_frac", ""), read_len=m.get("read_len", ""),
                  err_model=m.get("err_model", ""), kind=ds.loc[dataset, "kind"], tp=tp, fp=fp, tn=tn, fn=fn, sens=sens, spec=spec_, f1=f1, mcc=mcc,
                  residual_host_per_M=1e6 * fn / (tp + fn) if tp + fn else float("nan"),
                  microbial_loss=fp / (tn + fp) if tn + fp else float("nan"), n_calls=len(calls)))
os.makedirs(a.out_dir, exist_ok=True)
pd.DataFrame(S).to_csv(Path(a.out_dir) / "metrics.csv", index=False)
print(f"T4: {len(S)} runs -> {a.out_dir}")
