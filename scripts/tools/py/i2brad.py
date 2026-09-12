#!/usr/bin/env python
"""Generic motif-anchored k-mer sampler: in-silico type IIB restriction (2b-RAD) tags.

This is the *mechanism* row of the benchmark, deliberately tool-free: a k-mer is
selected iff it contains the enzyme's recognition motif at a fixed offset, the
same context-free predicate class as FracMinHash (hash below threshold) and
closed syncmers (minimal s-mer at a fixed position).  Fast2bRAD-M / Syn2bANI are
the tool rows and carry their own estimators.

Tag definition: the double-stranded core of the fragment, i.e. the recognition
site plus min(top, bottom) cut offset on each side (overhang bases excluded).
With the REBASE cut offsets below this reproduces the tag lengths of the
2bRAD-M/Fast2bRAD-M extractors for all 15 enzymes (BcgI 32, CspCI 33, BsaXI 27,
...).  Tags are strand-canonical (lexicographic min of tag / reverse complement).

Subcommands
  density   tags per genome per enzyme, GC, scaled-equivalent   (enzyme-density scan)
  pairs     T1: tag-set Jaccard / containment -> ANI per genome pair
  build-db  T2/T3: reference tag database (tag \\t genome), streamed at query time
  detect    T2: containment of each reference tag set in a read set
  sketch    T3 build: per-sample tag counter from reads
  bc        T3 query: Bray-Curtis (abund=1) or Jaccard distance (abund=0) between sketches
VERIFY(W1): recognition sites and cut offsets against REBASE for the enzymes you run."""
import argparse, gzip, itertools, math, os, re, subprocess, sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from names import genome_name

# name: (recognition site, (upstream top, upstream bottom), (downstream top, downstream bottom))  REBASE notation (a/b)SITE(c/d)
ENZYMES = {
    "CspCI":  ("CAANNNNNGTGG",  (11, 13), (12, 10)),
    "AloI":   ("GAACNNNNNNTCC", (7, 12),  (12, 7)),
    "BsaXI":  ("ACNNNNNCTCC",   (9, 12),  (10, 7)),
    "BaeI":   ("ACNNNNGTAYC",   (10, 15), (12, 7)),
    "BcgI":   ("CGANNNNNNTGC",  (10, 12), (12, 10)),
    "CjeI":   ("CCANNNNNNGT",   (8, 14),  (15, 9)),
    "PpiI":   ("GAACNNNNNCTC",  (7, 12),  (13, 8)),
    "PsrI":   ("GAACNNNNNNTAC", (7, 12),  (12, 7)),
    "BplI":   ("GAGNNNNNCTC",   (8, 13),  (13, 8)),
    "FalI":   ("AAGNNNNNCTT",   (8, 13),  (13, 8)),
    "Bsp24I": ("GACNNNNNNTGG",  (8, 13),  (12, 7)),
    "HaeIV":  ("GAYNNNNNRTC",   (7, 13),  (14, 9)),
    "CjePI":  ("CCANNNNNNNTC",  (7, 13),  (14, 8)),
    "Hin4I":  ("GAYNNNNNVTC",   (8, 13),  (13, 8)),
    "AlfI":   ("GCANNNNNNTGC",  (10, 12), (12, 10)),
}
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
         "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
DNA_COMP = str.maketrans("ACGT", "TGCA")

def rc(s): return s.translate(DNA_COMP)[::-1]
def rc_iupac(s): return s.translate(COMP)[::-1]

class Enzyme:
    def __init__(self, name):
        site, up, down = ENZYMES[name]
        self.name, self.site = name, site
        self.left, self.right = min(up), min(down)                 # double-stranded core
        self.k = self.left + len(site) + self.right
        fwd = "".join(IUPAC[c] for c in site); rev = "".join(IUPAC[c] for c in rc_iupac(site))
        self.pal = rc_iupac(site) == site
        # a palindromic site with symmetric offsets (AlfI) gives one window per site; with asymmetric offsets
        # (HaeIV) the two binding orientations give two different windows, both real tags -> scan both
        scan_rc = (not self.pal) or (self.left != self.right)
        self.re_f = re.compile(f"(?={fwd})"); self.re_r = re.compile(f"(?={rev})") if scan_rc else None
        self.L = len(site)
    def tags(self, seq):
        """Canonical tags of one sequence (uppercase, may contain N: such windows are skipped)."""
        out = []; n = len(seq); L = self.L
        for m in self.re_f.finditer(seq):
            s = m.start(); a, b = s - self.left, s + L + self.right
            if a >= 0 and b <= n:
                t = seq[a:b]
                if "N" not in t: r = rc(t); out.append(t if t <= r else r)
        if self.re_r is not None:
            for m in self.re_r.finditer(seq):                          # site on the bottom strand: offsets mirror
                s = m.start(); a, b = s - self.right, s + L + self.left
                if a >= 0 and b <= n:
                    t = seq[a:b]
                    if "N" not in t: r = rc(t); out.append(t if t <= r else r)
        return out
    def n_sites(self, seq):
        n = sum(1 for _ in self.re_f.finditer(seq))
        if not self.pal: n += sum(1 for _ in self.re_r.finditer(seq))       # physical sites: a palindrome is one site
        return n

def enzymes(spec):
    """'BcgI' | 'BcgI+AlfI' | 'BcgI,AlfI' -> [Enzyme]"""
    names = [x for x in re.split(r"[+,]", spec) if x]
    bad = [x for x in names if x not in ENZYMES]
    if bad: sys.exit(f"unknown enzyme(s) {bad}; known: {', '.join(ENZYMES)}")
    return [Enzyme(x) for x in names]

# ---------------------------------------------------------------- sequence IO
def _open(path):
    path = str(path)
    if path.endswith(".gz"):
        p = subprocess.Popen(["gzip", "-dc", path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        return p.stdout
    return open(path)

def fasta_records(path):
    name, buf = None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None: yield name, "".join(buf).upper()
                name, buf = line[1:].split()[0] if line[1:].split() else "", []
            else: buf.append(line.strip())
    if name is not None: yield name, "".join(buf).upper()

def read_seqs(path):
    """Sequences of a FASTA or FASTQ file (gz ok), sniffed from the first character."""
    with _open(path) as fh:
        first = fh.readline()
        if not first: return
        if first.startswith(">"):
            name, buf = first[1:].split()[0] if first[1:].split() else "", []
            for line in fh:
                if line.startswith(">"):
                    yield "".join(buf).upper(); buf = []
                else: buf.append(line.strip())
            yield "".join(buf).upper()
        else:                                                           # FASTQ: 4-line records
            yield fh.readline().strip().upper(); fh.readline(); fh.readline()
            while True:
                h = fh.readline()
                if not h: break
                yield fh.readline().strip().upper(); fh.readline(); fh.readline()

def genome_list(spec):
    """file of paths | directory | list of paths -> [paths]"""
    out = []
    for s in spec:
        p = Path(s)
        if p.is_dir(): out += sorted(str(x) for x in p.iterdir() if x.suffix in (".gz", ".fa", ".fna", ".fasta"))
        elif p.suffix in (".gz", ".fa", ".fna", ".fasta"): out.append(str(p))
        else: out += [l.split("\t")[-1].strip() if "\t" in l else l.strip() for l in open(p) if l.strip()]
    return out

# ---------------------------------------------------------------- workers
_ENZ = None
def _init(spec):
    global _ENZ; _ENZ = enzymes(spec)

def _genome_tags(path):
    tags = set()
    for _, seq in fasta_records(path):
        for e in _ENZ: tags.update(e.tags(seq))
    return genome_name(path), tags

def _genome_stats(path):
    L = gc = 0; sites = {e.name: 0 for e in _ENZ}; ntags = {e.name: set() for e in _ENZ}
    for _, seq in fasta_records(path):
        L += len(seq); gc += seq.count("G") + seq.count("C")
        for e in _ENZ: sites[e.name] += e.n_sites(seq); ntags[e.name].update(e.tags(seq))
    return genome_name(path), L, gc / L if L else float("nan"), sites, {k: len(v) for k, v in ntags.items()}

def _reads_counter(paths):
    c = Counter(); n_reads = 0
    for p in paths:
        for seq in read_seqs(p):
            n_reads += 1
            for e in _ENZ: c.update(e.tags(seq))
    return c, n_reads

def tag_len(spec):
    ks = {e.k for e in enzymes(spec)}
    return max(ks)                                                      # mixed panels: report the longest (used only for ANI transform)

# ---------------------------------------------------------------- subcommands
def cmd_density(a):
    paths = genome_list(a.genomes); enz = [e.name for e in enzymes(a.enz)]
    os.makedirs(a.out, exist_ok=True)
    rows = []
    with Pool(a.threads, initializer=_init, initargs=(a.enz,)) as pool:
        for gid, L, gc, sites, ntags in pool.imap_unordered(_genome_stats, paths, chunksize=4):
            rows.append((gid, L, gc, sites, ntags))
    with open(Path(a.out) / "per_genome.tsv", "w") as fh:
        fh.write("genome\tlength\tgc\t" + "\t".join(f"{e}_sites\t{e}_tags\t{e}_tags_per_mb\t{e}_scaled_equiv" for e in enz) + "\n")
        for gid, L, gc, sites, ntags in sorted(rows):
            cells = []
            for e in enz:
                k = Enzyme(e).k; n = ntags[e]
                cells += [str(sites[e]), str(n), f"{1e6 * n / L:.1f}" if L else "NA", f"{(L - k + 1) / n:.0f}" if n else "NA"]
            fh.write(f"{gid}\t{L}\t{gc:.4f}\t" + "\t".join(cells) + "\n")
    # summary per enzyme
    import numpy as np
    try: from scipy.stats import spearmanr
    except ImportError: spearmanr = None
    gcs = np.array([r[2] for r in rows]); Ls = np.array([r[1] for r in rows], float)
    with open(Path(a.out) / "summary.tsv", "w") as fh:
        fh.write("enzyme\ttag_len\tn_genomes\tmedian_tags\tmedian_tags_per_mb\tp10_tags_per_mb\tp90_tags_per_mb\tmedian_scaled_equiv\tfrac_lt_500_tags\tfrac_lt_1000_tags\tspearman_gc_density\n")
        for e in enz:
            n = np.array([r[4][e] for r in rows], float); d = 1e6 * n / Ls; k = Enzyme(e).k
            se = (Ls - k + 1) / np.maximum(n, 1)
            rho = spearmanr(gcs, d).correlation if spearmanr and len(rows) > 2 else float("nan")
            fh.write(f"{e}\t{k}\t{len(rows)}\t{np.median(n):.0f}\t{np.median(d):.1f}\t{np.percentile(d, 10):.1f}\t{np.percentile(d, 90):.1f}\t"
                     f"{np.median(se):.0f}\t{(n < 500).mean():.4f}\t{(n < 1000).mean():.4f}\t{rho:.3f}\n")
    print(f"density: {len(rows)} genomes x {len(enz)} enzymes -> {a.out}")

def wilson(x, n, z=1.959964):
    if n == 0: return float("nan"), float("nan")
    p = x / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)

def cmd_pairs(a):
    paths = genome_list(a.genomes); k = tag_len(a.enz)
    sets = {}
    with Pool(a.threads, initializer=_init, initargs=(a.enz,)) as pool:
        for gid, tags in pool.imap_unordered(_genome_tags, paths, chunksize=4): sets[gid] = tags
    names = sorted(sets)
    if a.pairs:                                                         # explicit pair list: pair_id path_a path_b genome_a genome_b
        todo = []
        for l in open(a.pairs):
            f = l.rstrip("\n").split("\t")
            if len(f) >= 3: todo.append((genome_name(f[1]), genome_name(f[2])))
    else:
        # all-vs-all through an inverted index: only pairs that share a tag get a row
        inv = {}
        for g in names:
            for t in sets[g]: inv.setdefault(t, []).append(g)
        shared = Counter()
        for gs in inv.values():
            if len(gs) > 1:
                for x, y in itertools.combinations(gs, 2): shared[(x, y) if x < y else (y, x)] += 1
        todo = list(shared)
    with open(a.out, "w") as fh:
        fh.write("genome_a\tgenome_b\tj_est\tani_est\tci_lo\tci_hi\tshared\tn_a\tn_b\tani_jaccard\n")
        for x, y in todo:
            A, B = sets.get(x, set()), sets.get(y, set()); s = len(A & B); u = len(A | B); n = min(len(A), len(B))
            j = s / u if u else float("nan"); c = s / n if n else 0.0
            if s == 0: ani, lo, hi, anij = "NA", "NA", "NA", "NA"
            else:
                ani = f"{c ** (1 / k):.6f}"; wl, wh = wilson(s, n); lo, hi = f"{wl ** (1 / k):.6f}", f"{wh ** (1 / k):.6f}"
                anij = f"{1 + math.log(2 * j / (1 + j)) / k:.6f}" if j > 0 else "NA"   # Mash transform, for the estimator comparison
            fh.write(f"{x}\t{y}\t{j:.6g}\t{ani}\t{lo}\t{hi}\t{s}\t{len(A)}\t{len(B)}\t{anij}\n")
    print(f"pairs: {len(names)} genomes, {len(todo)} pairs -> {a.out}")

def cmd_build_db(a):
    paths = genome_list(a.genomes); os.makedirs(a.out, exist_ok=True)
    n_tags = {}; n_lines = 0
    with gzip.open(Path(a.out) / "db.tsv.gz", "wt", compresslevel=4) as db, Pool(a.threads, initializer=_init, initargs=(a.enz,)) as pool:
        for gid, tags in pool.imap_unordered(_genome_tags, paths, chunksize=4):
            n_tags[gid] = len(tags)
            for t in tags: db.write(f"{t}\t{gid}\n"); n_lines += 1
    with open(Path(a.out) / "genomes.tsv", "w") as fh:
        fh.write("genome\tn_tags\n")
        for g in sorted(n_tags): fh.write(f"{g}\t{n_tags[g]}\n")
    Path(a.out, "enzymes.txt").write_text(a.enz + "\n")
    print(f"build-db: {len(n_tags)} genomes, {n_lines} tag records -> {a.out}")

def cmd_detect(a):
    _init(a.enz); k = tag_len(a.enz)
    counter, n_reads = _reads_counter(a.reads); total = sum(counter.values())
    hits = Counter(); depth = Counter()
    with gzip.open(Path(a.db) / "db.tsv.gz", "rt") as db:
        for line in db:
            t, g = line.rstrip("\n").split("\t")
            c = counter.get(t)
            if c: hits[g] += 1; depth[g] += c
    n_tags = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(Path(a.db) / "genomes.tsv").read().splitlines()[1:]}
    obs = max(a.read_len - k + 1, 1) / a.read_len                        # P(read fully contains a tag it overlaps)
    with open(a.out, "w") as fh:
        fh.write("genome\tscore\tcalled_present\tani_est\tcov_est\tabund_est\thits\tn_tags\n")
        for g, n in sorted(n_tags.items()):
            h = hits.get(g, 0)
            if h == 0: continue
            frac = h / n; called = int(h >= a.min_hits and frac >= a.min_frac)
            cov = depth[g] / n / obs; ab = depth[g] / total if total else 0.0
            fh.write(f"{g}\t{frac:.6g}\t{called}\t{'NA'}\t{cov:.4g}\t{ab:.4g}\t{h}\t{n}\n")
    print(f"detect: {n_reads} reads, {total} tag occurrences, {len(hits)} genomes with >=1 hit -> {a.out}")

def cmd_sketch(a):
    _init(a.enz); counter, n_reads = _reads_counter(a.reads)
    with gzip.open(a.out, "wt", compresslevel=4) as fh:
        fh.write(f"#reads={n_reads}\tenz={a.enz}\n")
        for t, c in counter.items(): fh.write(f"{t}\t{c}\n")
    print(f"sketch: {n_reads} reads -> {len(counter)} distinct tags -> {a.out}")

def load_sketch(path):
    c = {}
    with gzip.open(path, "rt") as fh:
        for l in fh:
            if l.startswith("#"): continue
            t, n = l.rstrip("\n").split("\t"); c[t] = int(n)
    return c

def cmd_bc(a):
    import numpy as np
    from scipy import sparse
    sids, files = [], []
    for l in open(a.samples):
        f = l.rstrip("\n").split("\t")
        if f and f[0]: sids.append(f[0]); files.append(f[1])
    idx = {}; rows, cols, vals = [], [], []
    for i, f in enumerate(files):
        for t, n in load_sketch(f).items():
            j = idx.setdefault(t, len(idx)); rows.append(i); cols.append(j); vals.append(1 if a.abund == 0 else n)
    M = sparse.csr_matrix((np.array(vals, float), (rows, cols)), shape=(len(files), len(idx)))
    tot = np.asarray(M.sum(1)).ravel()
    with open(a.out, "w") as fh:
        fh.write("sample_a\tsample_b\td_est\testimand\n")
        for i in range(len(sids)):
            for j in range(i + 1, len(sids)):
                if a.abund:
                    d = abs(M[i] - M[j]).sum() / (tot[i] + tot[j]) if tot[i] + tot[j] else float("nan"); est = "bray_curtis_tags"
                else:
                    inter = M[i].multiply(M[j]).nnz; uni = tot[i] + tot[j] - inter
                    d = 1 - inter / uni if uni else float("nan"); est = "jaccard_dist_tags"
                fh.write(f"{sids[i]}\t{sids[j]}\t{d:.6f}\t{est}\n")
    print(f"bc: {len(sids)} samples, {len(idx)} distinct tags -> {a.out}")

def cmd_tags(a):                                                        # debugging aid: print tags of a FASTA
    _init(a.enz)
    for name, seq in fasta_records(a.fasta):
        for e in _ENZ:
            for t in e.tags(seq): print(f"{name}\t{e.name}\t{t}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    def common(p, out=True):
        p.add_argument("--enz", default="BcgI", help="enzyme or panel, e.g. BcgI or BcgI+AlfI")
        p.add_argument("--threads", type=int, default=1)
        if out: p.add_argument("-o", "--out", required=True)
    p = sp.add_parser("density"); common(p); p.add_argument("--genomes", nargs="+", required=True); p.set_defaults(f=cmd_density)
    p = sp.add_parser("pairs"); common(p); p.add_argument("--genomes", nargs="+", required=True); p.add_argument("--pairs", help="pair list (pair_id path_a path_b ...); default all-vs-all"); p.set_defaults(f=cmd_pairs)
    p = sp.add_parser("build-db"); common(p); p.add_argument("--genomes", nargs="+", required=True); p.set_defaults(f=cmd_build_db)
    p = sp.add_parser("detect"); common(p); p.add_argument("--db", required=True); p.add_argument("--reads", nargs="+", required=True)
    p.add_argument("--min-hits", type=int, default=5); p.add_argument("--min-frac", type=float, default=0.001); p.add_argument("--read-len", type=int, default=151); p.set_defaults(f=cmd_detect)
    p = sp.add_parser("sketch"); common(p); p.add_argument("--reads", nargs="+", required=True); p.set_defaults(f=cmd_sketch)
    p = sp.add_parser("bc"); common(p); p.add_argument("--samples", required=True, help="sid \\t sketch.tsv.gz"); p.add_argument("--abund", type=int, default=1); p.set_defaults(f=cmd_bc)
    p = sp.add_parser("tags"); common(p, out=False); p.add_argument("fasta"); p.set_defaults(f=cmd_tags)
    a = ap.parse_args(); a.f(a)
