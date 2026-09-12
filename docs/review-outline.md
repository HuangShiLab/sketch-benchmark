# Review outline — k-mer sketching for biologists

> Superseded by the full draft in [manuscript.md](manuscript.md) (revision notes in [manuscript-changes.md](manuscript-changes.md)); kept as the planning record.

Working title: **k-mer sketching for biologists: mechanisms, size-matched benchmarks, and a build-your-own guide**
(the earlier "MinHash and its descendants" no longer covers the content once
motif-defined sampling is included).

One-sentence message: *ask what you are estimating before asking which tool is
fastest — sketch families differ in estimand and guarantee, and size-matched
benchmarks show the family matters more than the tool.*

Second message, for the agent era: *implementing a sketch is now cheap;
specifying its estimand, guarantee and test oracle is the scarce skill, and
this review supplies that vocabulary plus a reusable oracle
([sketch-benchmark](https://github.com/HuangShiLab/sketch-benchmark)).*

Target: benchmark-review with a tool catalogue — *Genome Biology* or
*Briefings in Bioinformatics*. If Section 6 is to lead, a *Nature Methods*
perspective is the shorter alternative (fewer results).

Length: ~7,000 words main text, 6 figures, 2 tables, 4 boxes, supplement.

---

## 1. Introduction — why "MinHash" is not one thing
- k-mer sets replaced alignment for scale; "MinHash" became an umbrella over
  ~30 tools estimating different quantities.
- Three failure modes seen in practice: Jaccard tools used for containment
  (metagenome vs genome); unweighted sketches used for abundance questions;
  set sketches applied to single reads (Box 3: a 150 bp read at `scaled=200`
  carries 0.6 expected sketch k-mers).
- What is new versus Rowe 2019, Marçais 2019, Zielezinski 2019, Ndiaye 2024:
  one framework across families, a size-matched benchmark, a build-your-own
  guide, and the first review to include an experimentally realizable sketch.

## 2. A unified framework (Fig 1, Box 1, Box 2)
- Pipeline: *estimand → k-mer (multi)set → selection rule → sketch structure →
  estimator → guarantee → cost*.
- The three quantities biologists actually want and their estimands: ANI /
  mutation rate (Jaccard or containment through the k-mer survival model);
  presence and coverage (containment); community dissimilarity (weighted
  Jaccard, Bray–Curtis proxies).
- **Box 1** "Six questions to ask of any sketch": what does it estimate; is the
  estimator unbiased; does it give a CI; does size scale with input; is it
  mergeable / streamable; is it comparable across tools.
- **Box 2** k-mer survival math: ANI ↔ Jaccard ↔ containment, the Mash
  transform, CIs (Blanca et al. 2022); why containment is the right estimand for
  asymmetric comparisons.

## 3. The mechanism families (Fig 2 gallery; Table 1; Table 2)
Same template per family: mechanism → estimand → guarantee → cost →
representative tools → where it breaks.

- **F1 fixed-size MinHash** — k-hash, bottom-k, one-permutation, b-bit; Mash,
  BinDash, Mash Screen.
- **F2 content-defined (context-free) sampling** — the k-mer decides for itself,
  so the same k-mer is selected in every genome and every read:
  - *hash threshold* — FracMinHash: sourmash, sylph, YACHT, CMash, MaxGeomHash;
  - *sub-hash position* — closed syncmers (also usable as an index sampler, §F6);
  - *sequence motif* — **2b-RAD tags**: a k-mer is selected iff it contains a
    type IIB recognition site at a fixed offset (BcgI: CGA·N6·TGC, 32 bp
    double-stranded core). Density ≈ 1 site / 2 kb, i.e. `scaled ≈ 2000`,
    set by the enzyme and tunable only discretely (enzyme choice, multi-enzyme
    panels). Non-uniform: motif frequency tracks GC and composition, so tag-poor
    genomes exist and uniform-sampling estimators need re-derivation (containment
    ANI holds, Jaccard ANI does not: a substitution inside the motif removes the
    tag from the denominator). The only selection rule an enzyme can execute
    (**Box 3**). Tools: Syn2bANI (ANI + SV), Fast2bRAD-M / MAP2B (profiling),
    Strain2bScan, sk2bGrow.
- **F3 weighted MinHash** — ICWS, BagMinHash, DartMinHash, ProbMinHash; HULK,
  GSearch. J_w ≠ J_P and why that matters for abundance.
- **F4 register sketches** — HLL, HyperMinHash, SetSketch; Dashing 1/2. Bytes
  per genome at fixed error.
- **F5 order-aware** — OrderMinHash, RabbitSketch.
- **F6 index samplers (not set sketches)** — minimizers, syncmers, strobemers,
  mod-minimizers; deacon, Kraken 2, minimap2. Window guarantees, context-freeness,
  why these are the right family for reads.
- **F7 vector / HD embeddings** — SimHash, HyperGen; LSH continuity, GPU fit.

**Table 1** family × property matrix: estimand; unbiased; CI; size ∝ input;
mergeable; streamable; context-free; density tunable; **experimentally
realizable**; implementation.
**Table 2** tool catalogue (~30 tools): year, language, family, defaults,
maintained.

## 4. Benchmark (Fig 3 accuracy vs bytes across tasks; Fig 4 per-task panels)
- Design: estimand-first eligibility; defaults **and** sketch-size-matched sweeps;
  two-tier truth (exact simulation; real-data consensus with spread); timing
  protocol ($TMPDIR copies, `/usr/bin/time -v`, 3 replicates, 1 and 16 threads).
- **T1** genome distance · **T2** detection in metagenomes · **T3** sample
  dissimilarity · **T4** read-level host depletion — one results subsection each.
- Cost axes: bytes per sketch, CPU time, and — new — *bases sequenced to obtain
  the sketch*: hash sketches need 100 % of the sample sequenced before they can
  be computed; a 2b-RAD sketch is sequenced directly at ~1 % of the bases.
- Cautionary rows: a family used outside its estimand (read-level sourmash /
  sylph in T4). Reported, never ranked as "losers".
- The syncmer-vs-minimizer result at matched index size and its dependence on
  sequencing error.
- 2b-RAD rows run under the same rules as everyone else (`i2brad` = mechanism,
  `syn2bani` / `fast2brad_m` = tools), size-matched against `scaled ≈ 2000`, with
  the negatives reported: GC dependence of density, no-estimate rate on tag-poor
  genomes (enzyme-density scan, `data/enzyme_density/`). T4 excluded — 2b-RAD
  handles host DNA by not sequencing it, not by filtering reads.

## 5. Decision guide for users (Fig 5 flowchart)
- Question → family → tool → parameters (k, s / scaled, w, thresholds).
- Pitfalls: k per task; hash-function compatibility across tools; sketch
  interchange; abundance weighting; CI availability; enzyme choice and tag-poor
  taxa for 2b-RAD.

## 6. Building your own sketch tool in the agent era (Fig 6; Box 4)
- The primitive catalogue: hashing (ntHash, xxhash), selection
  (`simd-minimizers`, syncmers, motif scanners), structures (sorted vectors, HLL
  registers, binary fuse filters), estimators, CIs — with the Rust / Python
  crates and APIs that compose them.
- Specification-first development: estimand, guarantee and test oracle before
  code; the benchmark harness as the oracle an agent runs.
- Worked example 1 — closed syncmers for deacon: what the human decided
  (estimand, threshold semantics, matched-size evaluation), what the agent
  produced (800-line patch, 107 upstream tests passing, smoke test), what the
  harness verified.
- Worked example 2 — one selection primitive, four tools: Syn2bANI (ANI + SV),
  Fast2bRAD-M (species), Strain2bScan (strain), sk2bGrow (growth dynamics).
  sk2bGrow is the clearest case in the review of a sketch used for an estimand
  beyond similarity.
- Agent-specific pitfalls: plausible-but-wrong CLI flags (`VERIFY(W1)` in this
  repo), silent estimand drift, tuning on the truth set, reproducing published
  numbers at unmatched sketch sizes.
- **Box 4** release checklist for a sketch tool.
- Reproducibility: pinned environments, manifests, one timing schema, idempotent
  stages.

## 7. Open problems and outlook
Long / noisy reads; pangenome-scale indexes and index algebra (union / subtract,
as in deacon's panhuman index); weighted containment with CIs; learned and GPU
sketches; a common sketch interchange format; machine-readable "mechanism
cards" agents can consume; sketches an instrument can measure directly.

## 8. Methods; data and code availability
HuangShiLab/sketch-benchmark, Zenodo snapshot, tool-version table.
Competing interests: the authors develop the 2b-RAD tools; the benchmark design
(this repository's first commit) predates any 2b-RAD result.

**Supplement**: full per-task tables, parameter sweeps, exact commands, tool
versions, minimal-run reproduction, enzyme-density scan.

---

## Figures, tables, boxes
| Item | Content |
|------|---------|
| Fig 1 | unified framework: estimand → selection → sketch → estimator → guarantee → cost |
| Fig 2 | mechanism gallery, one panel per family (F2 shows the three predicates side by side) |
| Fig 3 | accuracy vs bytes per sketch across T1–T3; second panel bases-sequenced axis |
| Fig 4 | per-task results T1–T4 |
| Fig 5 | decision flowchart |
| Fig 6 | spec → agent → oracle loop, with the two worked examples |
| Table 1 | family × property matrix |
| Table 2 | tool catalogue |
| Box 1 | six questions to ask of any sketch |
| Box 2 | k-mer survival math |
| Box 3 | one predicate, two executors: a hash or an enzyme (incl. the 150 bp read density argument) |
| Box 4 | release checklist |

## Placement rules for 2b-RAD (so the review does not read as self-promotion)
1. A subsection inside F2, one box, ≤ 10–15 % of the text; not in the title.
2. Benchmarked in silico under the same rules, size-matched; negatives reported.
3. Pre-registered design (repository history) and a competing-interests statement.
4. Tool series used as the *second* worked example in §6, not the first.

## Repository map
| Section | Where in the repo |
|---------|-------------------|
| §3 F2 motif sampling | `scripts/tools/py/i2brad.py` (enzyme table, tag definition, estimators) |
| §4 tasks | `scripts/stage02_simulate/`, `scripts/tools/*.sh`, `scripts/benchmark/metrics_T*.py` |
| §4 size matching | `scripts/benchmark/matched.py` |
| §4 enzyme-density scan | `scripts/analysis/enzyme_density.sh` → `data/enzyme_density/` |
| §6 example 1 | `env/deacon-syncmer.patch`, `tests/smoke_deacon_syncmer.sh` |
| §6 example 2 | `scripts/tools/syn2bani.sh`, `scripts/tools/fast2brad_m.sh` |
