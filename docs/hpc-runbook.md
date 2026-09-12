# HPC runbook — running the benchmark on HKU HPC2021

Operational plan for one person on the cluster. Every command is meant to be
pasted; every number comes from `scripts/hpc/config.sh` or from a dry run of
the manifests. Sizes and times are estimates to plan queue and quota with,
not measurements.

## 0. Goals, and what changed with the manuscript revision

The four tasks and the seven rules of `docs/manuscript.md` §5 are unchanged.
The revision did not re-aim the benchmark; it re-ordered it and made three
things explicit:

1. **A pre-step.** The enzyme-density scan must run before the T1/T2 motif
   rows are interpreted: it fixes the FracMinHash `scaled` that matches BcgI
   per genome (manuscript §2.3.4, Q8) and the no-estimate rate (§3.7).
2. **A new axis, already recorded.** Bases inside tags (Rule 7) are written
   by the `i2brad` rows; nothing extra to run.
3. **Deferred by decision, not by accident** (manuscript §5.8): long reads in
   T4, repeat-fraction strata in T1, the k ablation beyond Mash, GSearch /
   BinDash 2.0 / kssd / Kmer-db / Vclust / ska2, and the wet-lab concordance
   test. Do not add these to the first full run.

Run order, by what the paper leans on most per CPU-hour:

| Priority | What | Answers | Why first |
|---|---|---|---|
| A | minimal run (all tasks) | nothing scientific | surfaces every CLI-drift failure in ~2 h instead of mid-grid |
| A | enzyme-density scan | Q8 setup, §3.7 numbers | 1 job, ~1 h, feeds T1/T2 interpretation |
| A | T4 | §7.3 scheme contrast, T2T vs panhuman, host-depletion table; also rustyclean | cheapest task per result; every run is minutes |
| B | T1 sim + real | Q1 (CI coverage), Q8 (motif vs hash, size-matched), raw vs converted ANI | small data, many tools, Figure 3 |
| C | T2 | Q3 (safe scale factor), strain specificity, FDR | heaviest data (50 M-pair sets) |
| D | T3 | Q2 (J_w vs J_P), Q4 (depth) | Simka and per-sample sketches are the slow part |

## 1. Sizes to plan with

Full mode (`SB_MINIMAL` unset), from a dry manifest build:

| Task | Datasets | Build tasks | Query tasks | Simulated data | Notes |
|---|---|---|---|---|---|
| T1 | sim, real, scale | 93 | ~240 | ~2 GB (950 mutant FASTAs) | scale = 5,000 genomes all-vs-all, sketch tools + i2brad + syn2bani only |
| T2 | 9 sim + 2 Zymo (+ CAMI if staged) | 11 | 554 | ~45 GB (six 10 M-pair and three 50 M-pair sets) | InSilicoSeq is the slow step |
| T3 | 60 + 9 depth replicas (+ HMP) | 12 | 53 | ~75 GB (5 M pairs/sample) | Simka needs ~48 GB RAM |
| T4 | 15 sim + 4 real | 20 | 1,426 | ~30 GB sim + real sets | deacon variants are 800 of the 1,426; each is minutes |

Working space: reserve **~400 GB** on `/lustre1/g/aos_shihuang` for
`data/sim`, `data/real`, `idx` and `runs` (the Kraken 2 database built from
the 5,000-genome pool alone is ~20 GB). Cap the real T4 human set to 10 M
pairs (stage 00 does this for HMP; do the same for HG00438 by pointing the
URL at one lane or subsampling after download).

Compute, order of magnitude: ~5,000–8,000 core-hours in total; wall time
1–2 days per task with the array concurrency limits below, dominated by
InSilicoSeq (simulation), Bowtie2/hostile (T4), and Kraken 2 on 50 M pairs (T2).

## 2. One-time setup (login node)

```bash
cd /lustre1/g/aos_shihuang            # or your project space
git clone git@github.com:HuangShiLab/sketch-benchmark.git
cd sketch-benchmark
source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda env create -f env/environment.yml           # env name: sketch-benchmark (30-60 min)
conda activate sketch-benchmark
module load rust 2>/dev/null || curl https://sh.rustup.rs -sSf | sh -s -- -y   # cargo for deacon-syncmer, syn2bani, fast2bRAD-M, gsearch
bash env/cargo-tools.sh                           # builds into env/bin; warns (does not stop) on optional failures
ls env/bin                                        # expect: deacon-syncmer syn2bani fast2bRAD-M [gsearch] [dashing2]
```

Then edit `scripts/hpc/config.sh`:

| Variable | Set to |
|---|---|
| `SB_PARTITION`, `SB_QOS` | your partition/QoS (defaults `amd`, `normal`) |
| `DB_ROOT`, `MICROBIAL_GENOME_DIR`, `T2T_FASTA`, `HUMAN_GENOME`, `KRAKEN2_DB_*`, `BOWTIE2_INDEX`, `MINIMAP2_INDEX`, `HOSTILE_INDEX` | the existing database paths (rustyclean-paper conventions) |
| `HG002_MAT_URL`, `HG002_PAT_URL` | HPRC HG002 maternal/paternal assemblies (`.fa.gz`) — T4 divergence panel |
| `HG00438_R1_URL`, `HG00438_R2_URL` | one Illumina lane of an HPRC individual — T4 real host set |
| `ZYMO_EVEN_SRR`, `ZYMO_LOG_SRR`, `ZYMO_REF_URL` | one 2×150 run each of D6300 and D6310, and the reference bundle — T2 real tier |
| `CAMI_MARINE_DIR`, `CAMI_STRAIN_DIR` | directories you stage by hand from the CAMI portal (`*_1.fastq.gz`, `*_2.fastq.gz`, `<sid>.truth.tsv`), or leave `TODO` |
| `DASHING2_URL` | a linux x86_64 release asset (or leave `TODO`: Dashing 2 rows are skipped) |
| `SYN2BANI_TAG`, `FAST2BRAD_TAG`, `DEACON_UPSTREAM_TAG` | pin commits before the full run (`main` is fine for the minimal run) |
| `refs/hmp_samples.tsv` | `sample_id  site  srr` rows for the HMP tier of T3, or leave header-only to skip |

Anything left as `TODO` is recorded as skipped in `refs/manifest.tsv`; the
pipeline runs without it.

## 3. W1: interface verification (30–60 min, login node)

Every wrapper marked `VERIFY(W1)` guessed a CLI detail. Check them against
the installed versions before submitting anything:

```bash
grep -rn "VERIFY(W1)" scripts/tools | cut -d: -f1 | sort -u
bindash sketch --help | head -40          # scripts/tools/bindash.sh: flag names, sketchsize64
dashing2 sketch --help | head -40         # scripts/tools/dashing2.sh: -k, -S, --cmpout layout
hulk smash --help                         # scripts/tools/hulk.sh: -m weightedjaccard, output matrix name
yacht --help; yacht run --help            # scripts/tools/yacht.sh
sylph query --help                        # scripts/tools/sylph.sh: genome-as-query mode, column names
python -c "import sourmash; m=sourmash.MinHash(0,21,scaled=1000); print(m.containment_ani)"   # sourmash_pairs.py
syn2bani dist --help; syn2bani sketch --help          # scripts/tools/syn2bani.sh: --enzymes, units of ani/std_err
fast2bRAD-M quantify --help; fast2bRAD-M extract --help   # scripts/tools/fast2brad_m.sh: -l list format, output columns
```

Fix the wrapper where the installed CLI differs, delete the `VERIFY(W1)`
comment, and run the two local test suites (they do not need the tools):

```bash
python tests/test_metrics.py && python tests/test_i2brad.py
bash tests/smoke_deacon_syncmer.sh env/bin/deacon-syncmer       # ~1 min; must print PASS
```

## 4. Minimal run (Priority A; ~2–4 h wall)

```bash
bash scripts/run_all.sh --minimal --dry-run       # prints every sbatch line, submits nothing
bash scripts/run_all.sh --minimal                 # data/minimal, idx/minimal, runs/minimal
tail -f logs/jobs.tsv                             # one line per submitted job: stage, job id
squeue -u $USER
```

It runs stage 00 (downloads), 02 (all simulations, 10× smaller), 02b
(manifests; this job submits the run phase itself), 03/04/05/06 per task.
When it finishes:

```bash
ls data/minimal/T*/metrics.csv data/minimal/T*/summary.csv data/minimal/T*/timing.csv
grep -l "tool failed" logs/sb-04-*.out | head          # any wrapper that died
awk -F'\t' '$1=="T4"' data/minimal/T4/timing.csv | head
```

Every tool must produce a metrics row here. A wrapper that fails in minimal
mode will fail 100 times in the full grid; fix it now. Failed array tasks are
re-run by resubmitting the stage — completed runs carry `DONE` and are skipped:

```bash
bash scripts/run_all.sh --minimal --phase run --tasks T2      # only T2, only build→query→metrics
```

## 5. Enzyme-density scan (Priority A; 1 job)

```bash
sbatch scripts/analysis/enzyme_density.sh          # 16 CPUs, ~1 h on 5,000 genomes
column -t data/enzyme_density/summary.tsv           # tag_len, median tags, tags/Mb, scaled-equivalent, GC Spearman, tag-poor fraction
```

Needs `refs/gtdb5k.list`, written by the full-mode stage 00 (`--phase sim`); a
minimal run writes its own lists under `refs/minimal/` and `SB_MINIMAL=1 bash
scripts/analysis/enzyme_density.sh` scans the 10 T1 genomes as a smoke test. The
summary feeds manuscript §2.3.4/§3.7 directly and is the first result worth
reading.

## 6. Full run

Phase `sim` first (downloads + simulations), then `run` per task in priority
order so that results arrive while later simulations are still running.

```bash
bash scripts/run_all.sh --dry-run                            # check the full plan
bash scripts/run_all.sh --phase sim --no-continue            # 00 fetch → 02 t1/t2/t3/t4 → 02b manifests (does not start queries)
```

Simulation wall times (16 CPUs): T1 minutes; T4 ~40 min per dataset (15
datasets, array); T3 ~20 min per sample (69, array of 80, %-limited); T2 10 M
pairs ~30 min, 50 M pairs ~3 h. Expect the sim phase to take most of a day.

Then, as each task's simulations finish (`ls data/sim/T4/*/DONE | wc -l`):

```bash
bash scripts/run_all.sh --phase run --tasks T4               # 20 builds → 1,426 queries (%40) → metrics
bash scripts/run_all.sh --phase run --tasks T1
bash scripts/run_all.sh --phase run --tasks T2
bash scripts/run_all.sh --phase run --tasks T3
bash scripts/run_all.sh --phase run --tasks T1,T2,T3,T4      # re-run anything that failed; DONE runs are skipped
```

`--after JOBID` chains a run phase behind a still-running simulation job.
Stage 03 arrays are limited to 20 concurrent tasks, stage 04 to 40; raise
`%20`/`%40` in `scripts/run_all.sh` if the partition is idle.

Monitoring:

```bash
squeue -u $USER -o "%.10i %.20j %.8T %.10M %.6D %R"
tail -3 logs/jobs.tsv
find runs/T4 -name DONE | wc -l; wc -l < data/T4/manifest.tsv        # progress = DONE / manifest
grep -L DONE $(find runs/T4 -name tool.log | sed 's/tool.log//' | sed 's#$#DONE#') 2>/dev/null | head   # runs without DONE
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS,ExitCode | tail
```

## 7. Where results land and how they map to the manuscript

| File | Content | Manuscript |
|---|---|---|
| `data/<T>/metrics.csv` | one row per (tool, params, dataset, unit) | supplement |
| `data/<T>/summary.csv` | per (tool, params, dataset): bias/RMSE/MAE per bin, CI coverage, no-estimate, Spearman (T1); P/R/F1, AUPRC, FDR, recall per coverage bin (T2); Mantel vs BC/Jaccard/J_P, ARI, Procrustes, depth stability (T3); sens/spec/F1/MCC, residual host per M (T4) | Table 5, Figure 3 |
| `data/<T>/timing.csv` | one schema for all tasks: wall/user/sys, max RSS, input bp, sketch and index bytes, node, CPU | Table 5 timing columns, Figure 2 x-axis |
| `data/<T>/matched.csv` | accuracy at the size-matched budget (32 KB genome, 25 MB sample) | Rule 1 comparisons |
| `data/enzyme_density/summary.tsv` | tags per genome, GC dependence, scaled-equivalent | §2.3.4, §3.7, Q8 |
| `runs/T2/i2brad/*/detect.tsv.stat`, T3 sketch headers | bases inside tags / bases scanned | Rule 7 panel of Figure 2 |
| `figures/` | stage 06 cross-task panels | Figure 2 draft |

`bash scripts/stage05_metrics.sh` and `stage06_figures.sh` can be re-run on
a login node at any time (`SB_TASK=T1 bash scripts/stage05_metrics.sh`).

## 8. Failure playbook

| Symptom | Cause | Fix |
|---|---|---|
| `tool --help does not mention '<flag>'` in tool.log | CLI drift (`verify_flags`) | edit the wrapper; resubmit the stage |
| many `sb-04-*` tasks end in seconds with rc≠0 | one wrapper broken | `grep -h "tool failed" logs/sb-04-T*.out | sort | uniq -c` to find it |
| stage 02 T2 tasks time out | 50 M-pair InSilicoSeq | raise `--time` in `t2.sh` (default 12 h) or run seeds 1–2 later |
| `DONE` missing but no error | job killed by walltime | resubmit; partial outputs are overwritten |
| `sbatch: error: ... QOSMaxSubmitJobPerUserLimit` | queue limit | `submit.sh` retries every 120 s automatically; nothing to do |
| `siblings missing; t1.sh must run first` (T2 strain) | T2 started before T1 | rerun `--phase sim` for T2 after T1's sim job finishes |
| Kraken 2 T2 build OOM | 5,000-genome DB | raise `--mem` in `stage03_build.sh` to 128G for `kraken2` |
| Simka OOM | 69 samples | `-max-memory` in `simka.sh`, node with 192 GB |

## 9. Hand-off checklist

Before numbers go into the manuscript: every `VERIFY(W1)` comment removed;
`refs/tool_versions.tsv` committed; `data/*/summary.csv`, `timing.csv`,
`matched.csv`, `enzyme_density/` committed; `runs/` and `idx/` archived to
Zenodo (they are gitignored); the repository tagged (`git tag v0.1-results`).
