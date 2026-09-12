# sketch-benchmark

A benchmark harness for "MinHash-family" sequence sketching tools, organised by
**sketching mechanism** rather than by tool. It backs a benchmark-style review and
the sketch-based host-depletion work in
[rustyclean](https://github.com/HuangShiLab/rustyclean).

Seven mechanism families are covered (F1 fixed-size MinHash, F2 content-defined
sampling — FracMinHash, syncmers and motif-defined 2b-RAD tags, F3 weighted MinHash,
F4 register sketches such as HyperLogLog/SetSketch, F5 order-aware sketches,
F6 index samplers such as minimizers, F7 vector/HD embeddings) across four tasks:

| Task | Question | Estimand | Truth | Standard output |
|------|----------|----------|-------|-----------------|
| **T1** | genome-vs-genome distance | ANI / Jaccard | mutated references (exact ANI); real pairs = median of FastANI, skani, ANIm | `pairs.tsv` |
| **T2** | genome detection in a metagenome | presence, coverage, ANI | spike-in simulations against the GTDB-5k pool | `detect.tsv` |
| **T3** | sample-vs-sample dissimilarity | Bray–Curtis / Jaccard rank structure | designed communities with known abundances; HMP body-site labels | `dist.tsv` |
| **T4** | read-level host depletion | per-read host / other | simulated human + microbial reads with labels; real human runs | `calls.tsv` |

Every tool is run at its **defaults** and along a **sketch-size sweep**, so accuracy can
be reported against bytes per sketch (`matched.csv`) and not only at the vendors'
defaults. Timing uses `/usr/bin/time -v` on inputs copied to `$TMPDIR`, three
replicates, at 1 and 16 threads.

The design document (datasets, ground truth, metrics, eligibility rules) is the
companion work plan; this repository is its executable form. The outline of the
review it feeds is in [docs/review-outline.md](docs/review-outline.md).

## Layout

```
env/                 conda environment, cargo-built tools, deacon syncmer patch
scripts/hpc/         config.sh (single source of truth) and the sbatch wrapper
scripts/run_all.sh   chains every SLURM stage with dependencies
scripts/stage00_fetch.sh          references, public datasets, tool versions
scripts/stage02_simulate/         T1 mutations, T2 spike-ins, T3 communities, T4 sweeps (ISS)
scripts/stage02b_manifests.sh     discovers datasets, writes build/query manifests
scripts/stage03_build.sh          index/sketch building   (array over build_manifest.tsv)
scripts/stage04_query.sh          querying + timing       (array over manifest.tsv)
scripts/stage05_metrics.sh        metrics_T1..T4.py, matched.py
scripts/stage06_figures.sh        cross-task figures
scripts/tools/<tool>.sh           one wrapper per tool: build + query, standard outputs
scripts/tools/py/                 parsers shared by the wrappers; i2brad.py (in-silico 2b-RAD tag sampler); names.py (one genome-naming rule)
scripts/analysis/                 stand-alone analyses (enzyme-density scan)
scripts/benchmark/                manifests, timing rows, metrics, figures
docs/                             review outline
tests/                            metrics tests on synthetic fixtures; deacon syncmer smoke test
```

Large trees (`refs/`, `idx/`, `runs/`, `data/sim`, `data/real`) never enter git;
tidy tables under `data/` and `figures/` do.

## Quick start (HKU HPC2021 / SLURM)

```bash
git clone git@github.com:HuangShiLab/sketch-benchmark.git
cd sketch-benchmark
conda env create -f env/environment.yml        # -> sketch-benchmark
bash env/cargo-tools.sh                        # deacon (+ syncmer fork), gsearch, dashing2
```

Edit `scripts/hpc/config.sh`: the `TODO` download URLs / SRA accessions, the existing
database paths under `DB_ROOT`, partition and QoS. Then:

```bash
bash scripts/run_all.sh --minimal --dry-run    # print every sbatch line, submit nothing
bash scripts/run_all.sh --minimal              # 10x-smaller end-to-end run (data/minimal, idx/minimal, runs/minimal)
bash scripts/run_all.sh                        # full run
```

`run_all.sh` options: `--phase {all,sim,run}`, `--tasks T1,T2` (subset), `--after JOBID`
(chain behind an existing job), `--no-continue` (do not let stage 02b submit the run
phase). Job ids are appended to `logs/jobs.tsv`.

Stages are idempotent: each dataset, index and run directory carries a `DONE` marker and
is skipped when present, so a failed array task is re-run by resubmitting the stage.

## Conventions

* `scripts/hpc/config.sh` is the only place paths, grids and seeds live. `SB_MINIMAL=1`
  redirects to `data/minimal`, `idx/minimal`, `runs/minimal` and shrinks every grid.
* Wrapper interface: `build <task> <refset> <params> <threads> <idx_dir>` and
  `query <task> <dataset> <params> <threads> <work_dir> <idx_dir> <out_dir>`.
  `params` is a `k=31,w=15` string, parsed with `param KEY DEFAULT` from
  `scripts/tools/common.sh`.
* Standard outputs (one file per run directory, read by the metrics scripts):
  `pairs.tsv` (`genome_a genome_b j_est ani_est ci_lo ci_hi`),
  `detect.tsv` (`genome score called_present ani_est cov_est abund_est`),
  `dist.tsv` (`sample_a sample_b d_est estimand`),
  `calls.tsv` (`read_id call`, call in `host|other`).
* Timing rows (`data/<task>/timing.csv`) share one schema:
  `task,tool,tool_version,family,stage,dataset,params,threads,rep,wall_s,user_s,sys_s,max_rss_kb,input_bp,sketch_bytes,index_bytes,node,cpu_model,timestamp`.
* k = 21 for T1 and k = 31 for T2–T4; simulation seed 42; ISS `novaseq` model, 151 bp.
* Tools that estimate a different quantity from the task's estimand are still run but
  flagged as cautionary rows (e.g. read-level sketching in T4), never as "losers".

## Motif-defined sampling (2b-RAD tags)

`scripts/tools/py/i2brad.py` is a tool-free sampler for type IIB restriction tags:
a k-mer is selected iff it contains the enzyme's recognition motif at a fixed offset
(15 enzymes; the double-stranded-core tag definition reproduces the Fast2bRAD-M tag
lengths for all of them, and tags are strand-invariant — `tests/test_i2brad.py`).
It is the *mechanism* row in T1–T3 (`i2brad`); the *tool* rows are `syn2bani`
(T1) and `fast2brad_m` (T2/T3), built by `env/cargo-tools.sh`. The enzyme-density
scan (`sbatch scripts/analysis/enzyme_density.sh`) writes tags per genome, GC
dependence and the FracMinHash `scaled`-equivalent for `$DENSITY_ENZ` over the
GTDB-5k pool — on random sequence BcgI comes out at `scaled ≈ 2100`.

## The deacon syncmer fork

`env/deacon-syncmer.patch` adds a `--scheme {minimizer,syncmer}` option to
`deacon index build` (closed syncmers via `simd-minimizers`, index format version 4,
`deacon index info` reports the scheme). It applies cleanly to upstream commit
`f2fa660d6dfc46dfa04d2b2dca1adeda8e2cfe27` (`DEACON_UPSTREAM_TAG` in `config.sh`);
all 107 upstream tests pass with the patch, and `tests/smoke_deacon_syncmer.sh`
checks the two schemes against each other on mutated reads. `env/cargo-tools.sh`
builds it as `$SB_BIN/deacon-syncmer`.

## Tests

```bash
python tests/test_metrics.py            # metrics_T1..T4, matched.py, timing_row.py, make_manifests.py on synthetic fixtures
python tests/test_i2brad.py             # 2b-RAD sampler: tag lengths, strand invariance, containment ANI, detect, Bray-Curtis
bash tests/smoke_deacon_syncmer.sh [path/to/deacon-syncmer]   # needs the patched binary
```

## Before the first full run (W1 checklist)

Nothing in this repository has been executed on the cluster yet. Items marked
`VERIFY(W1)` in the wrappers are CLI details that must be checked against the
installed versions, and `TODO` entries in `config.sh` are download locations:

* `bindash`, `dashing2`, `hulk`, `yacht`, `sylph` flag names and output columns
  (`scripts/tools/*.sh`, `scripts/tools/py/matrix_to_pairs.py`).
* `sourmash` Python API: `MinHash.containment_ani(estimate_ci=True)` field names
  (`scripts/tools/py/sourmash_pairs.py`).
* `maxgeomhash.sh` and `hypergen.sh` are stubs that exit 3 until the tools are installed.
* `syn2bani.sh` (ANI units and `std_err` scale, `sketch --enzymes`) and `fast2brad_m.sh`
  (columns of `*.GCF_detected.xls`, `quantify -l` list format); pin `SYN2BANI_TAG` /
  `FAST2BRAD_TAG` in `config.sh`. Recognition sites and cut offsets in `i2brad.py` against REBASE.
* `config.sh` TODOs: HPRC HG002 assemblies and HG00438 reads, ZymoBIOMICS runs and
  reference bundle, CAMI marine/strain-madness staging directories, dashing2 binary URL;
  `refs/hmp_samples.tsv` (columns `sample_id site srr`) must be curated by hand.

Run `bash scripts/run_all.sh --minimal` first; the minimal run exercises every wrapper
on small inputs and surfaces CLI drift before the full grids are submitted.

## License

MIT, see `LICENSE`.
