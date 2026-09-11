#!/usr/bin/env python
"""T3 metrics. Each run's dist.tsv (all pairs) is scored against the true Bray-Curtis AND the
true Jaccard matrices (so the estimand mismatch is visible), plus group recovery and depth stability.
  metrics.csv  per (tool, params, pair): d_est, bc_true, jac_true, estimand
  summary.csv  per (tool, params): mantel_bc, mantel_jac (Spearman, 999 perms, p), ari_ward3,
               procrustes_bc, depth_stability (mean d between depth replicas of the same community)
Real (HMP) samples have no numeric truth; they get ari against site labels only."""
import argparse, itertools, math, os
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.spatial import procrustes
from common import runs, read_tsv, fnum

def mantel(D1, D2, perms=999, seed=0):
    n = D1.shape[0]; iu = np.triu_indices(n, 1)
    r = spearmanr(D1[iu], D2[iu]).correlation
    rng = np.random.default_rng(seed); cnt = 0
    for _ in range(perms):
        p = rng.permutation(n); rp = spearmanr(D1[iu], D2[np.ix_(p, p)][iu]).correlation
        if rp >= r: cnt += 1
    return r, (cnt + 1) / (perms + 1)
def ari(a, b):
    a = pd.Series(a).astype("category").cat.codes.values; b = pd.Series(b).astype("category").cat.codes.values
    n = len(a); ct = pd.crosstab(a, b).values
    def c2(x): return x * (x - 1) / 2
    sum_ij = c2(ct).sum(); sa = c2(ct.sum(1)).sum(); sb = c2(ct.sum(0)).sum(); tot = c2(n)
    exp = sa * sb / tot; mx = (sa + sb) / 2
    return (sum_ij - exp) / (mx - exp) if mx != exp else 1.0
def pcoa(D, k=2):
    n = D.shape[0]; J = np.eye(n) - np.ones((n, n)) / n; B = -0.5 * J @ (D ** 2) @ J
    w, v = np.linalg.eigh(B); idx = np.argsort(w)[::-1][:k]
    return v[:, idx] * np.sqrt(np.clip(w[idx], 0, None))

ap = argparse.ArgumentParser(); ap.add_argument("--runs", required=True); ap.add_argument("--datasets", required=True); ap.add_argument("--out-dir", required=True)
a = ap.parse_args()
ds = read_tsv(a.datasets)
sim = ds[ds.kind == "sim"]; truth_path = sim.truth.iloc[0] if len(sim) else None
meta = {}
for _, r in ds.iterrows():
    if r["meta"] and os.path.exists(r["meta"]):
        m = read_tsv(r["meta"]).iloc[0]; meta[r["dataset"]] = dict(group=m.get("group", ""), base=m.get("base_sample", r["dataset"]), pairs=fnum(m.get("pairs", "nan")))
truth = None
if truth_path and os.path.exists(truth_path):
    T = read_tsv(truth_path); truth = {tuple(sorted((x, y))): (fnum(b), fnum(j)) for x, y, b, j in zip(T.sample_a, T.sample_b, T.bc_true, T.jac_true)}
hmp = None
hl = os.environ.get("HMP_SAMPLE_LIST", "")
if hl and os.path.exists(hl):
    H = read_tsv(hl); hmp = dict(zip(H.iloc[:, 0], H.iloc[:, 1]))
rows, S = [], []
for tool, dataset, params, threads, rep, d in runs(a.runs, "dist.tsv"):
    if rep != 1: continue
    E = read_tsv(d / "dist.tsv"); E["d_est"] = E.d_est.map(fnum)
    est = E.estimand.iloc[0] if len(E) else "NA"
    names = sorted(set(E.sample_a) | set(E.sample_b)); ix = {s: i for i, s in enumerate(names)}
    D = np.zeros((len(names), len(names)))
    for x, y, v in zip(E.sample_a, E.sample_b, E.d_est):
        D[ix[x], ix[y]] = D[ix[y], ix[x]] = v
    D = np.nan_to_num(D, nan=1.0)
    simn = [s for s in names if truth and s in meta]
    row = dict(tool=tool, params=params, estimand=est, n_samples=len(names))
    if truth and len(simn) > 3:
        si = [ix[s] for s in simn]; Ds = D[np.ix_(si, si)]
        B = np.zeros_like(Ds); Jm = np.zeros_like(Ds)
        for i, j in itertools.combinations(range(len(simn)), 2):
            b, jc = truth.get(tuple(sorted((simn[i], simn[j]))), (np.nan, np.nan))
            B[i, j] = B[j, i] = b; Jm[i, j] = Jm[j, i] = jc
            rows.append(dict(tool=tool, params=params, sample_a=simn[i], sample_b=simn[j], d_est=Ds[i, j], bc_true=b, jac_true=jc, estimand=est))
        base = [s for s in simn if meta[s]["base"] == s]          # exclude depth replicas from structure metrics
        bi = [simn.index(s) for s in base]
        if len(bi) > 3:
            Db, Bb, Jb = Ds[np.ix_(bi, bi)], B[np.ix_(bi, bi)], Jm[np.ix_(bi, bi)]
            row["mantel_bc"], row["mantel_bc_p"] = mantel(Db, np.nan_to_num(Bb))
            row["mantel_jac"], row["mantel_jac_p"] = mantel(Db, np.nan_to_num(Jb))
            groups = [meta[s]["group"] for s in base]
            try:
                Z = linkage(squareform(Db, checks=False), "ward"); cl = fcluster(Z, len(set(groups)), "maxclust")
                row["ari_ward"] = ari(groups, cl)
            except Exception: row["ari_ward"] = np.nan
            try: row["procrustes_bc"] = procrustes(pcoa(np.nan_to_num(Bb)), pcoa(Db))[2]
            except Exception: row["procrustes_bc"] = np.nan
        # depth stability: distance between replicas of the same base community
        reps = [(s, meta[s]["base"]) for s in simn if meta[s]["base"] != s]
        dd = [Ds[simn.index(s), simn.index(b)] for s, b in reps if b in simn]
        row["depth_stability"] = float(np.mean(dd)) if dd else np.nan
        row["depth_stability_max"] = float(np.max(dd)) if dd else np.nan
    if hmp:
        hn = [s for s in names if s in hmp]
        if len(hn) > 3:
            hi = [ix[s] for s in hn]; Dh = D[np.ix_(hi, hi)]; sites = [hmp[s] for s in hn]
            try:
                Z = linkage(squareform(Dh, checks=False), "ward"); cl = fcluster(Z, len(set(sites)), "maxclust"); row["ari_hmp_site"] = ari(sites, cl)
            except Exception: row["ari_hmp_site"] = np.nan
    S.append(row)
os.makedirs(a.out_dir, exist_ok=True)
pd.DataFrame(rows).to_csv(Path(a.out_dir) / "metrics.csv", index=False)
pd.DataFrame(S).to_csv(Path(a.out_dir) / "summary.csv", index=False)
print(f"T3: {len(S)} runs -> {a.out_dir}")
