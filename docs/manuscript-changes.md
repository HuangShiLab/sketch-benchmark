# Manuscript revision notes

`docs/manuscript.md` (and its export `docs/manuscript.docx`, made by
`docs/tools/md2docx.py`) is a revision of the first draft
(`minhash-review-draft.docx`, 16.1k words, 115 references). This file records
what changed, what must be verified before submission, and what still has to
be decided. The revision is 21.4k words; a cut list to journal length is at the
end.

## What was kept

The draft's spine is sound and was kept verbatim wherever it was not touched:
the selection-rule/estimator split (§2.1), the seven families, the property
matrix (Table 2), the failure-mode catalogue (§3.1–3.6), the eligibility
matrix with a reason per cell (§4.2), the benchmark rules (§5.1), the
LLM-deduplication cross-pollination (§6), the open questions and the
conclusions. Citation numbering [1]–[115] is unchanged; new references are
[116]–[126].

## What changed, by section

| Where | Change |
|---|---|
| Title, authors | Added a working title and an author/correspondence placeholder. |
| Abstract | Rewritten: the three context-free predicates (hash threshold, sub-hash position, sequence motif), the sequencing-cost axis, and the build-your-own primer. |
| §1.3–1.5 | Scope now covers any selection rule statable in the same terms, names restriction-tag sketching as unreviewed; five contributions (added the primer, §7); organization updated. |
| §2.1 | Context-freedom defined as "a predicate on the k-mer alone"; the three predicates named; Figure 1 gains panel (e) (motif rule); Table 1 F2 row extended; intersection-commuting statement now covers threshold, syncmer and motif rules. |
| §2.3.4 (new) | *Motif-defined sampling: restriction tags as a context-free predicate* — mechanism, the three differences from hashing (discrete density, non-uniform sampling, enzyme-executable), containment-versus-Jaccard estimator argument, the tool cluster, and the competing-interest note. Former §2.3.4 (kssd) is now §2.3.5. |
| Table 2 | Three rows added: in-silico 2b-RAD tags (mechanism), Syn2bANI, Fast2bRAD-M/MAP2B. Reading paragraph adjusted. |
| §2.7.3, Table 4 | Cross-reference to motif tags; timeline rows for 2012 (2b-RAD), 2022 (2bRAD-M) and 2026 (i2bRAD-M, Syn2bANI, Fast2bRAD-M). |
| §3.7 (new) | *Non-uniform sampling: when the predicate is a motif* — tag-poor genomes as a no-estimate failure, why ANI must go through containment, and the computed-vs-sequenced concordance question. Cites our simulation: 0.990 ± 0.001 at 1% divergence from ~1,450 BcgI tags (`tests/test_i2brad.py`). |
| §4.1.1–4.1.4 | One paragraph each on how the motif family meets or fails the task's estimand. |
| Table 3 | F2 split into F2h (hash) and F2m (motif); reasons s and t added; the "read by column" paragraph updated. |
| §4.3 | 150 bp read: ~120 k-mers at k = 31, 0.12 expected hashes (was 150 / 0.15). |
| §5 (rewritten) | Now matches the implemented pipeline: seven rules (Rule 7 = sequencing-cost axis), the actual grids (divergence rates, coverages, community design, T4 sweeps), the real-data tiers that are configured, the entrants that have wrappers, the metrics the scripts compute, size-matched budgets (32 KB / 25 MB), the deacon scheme contrast, Figure 2's fifth panel, Table 5's bases-sequenced column, a concrete reproducibility section, and §5.8 *Planned extensions* listing what the draft promised but the pipeline does not yet do (long reads, repeat stratification, k ablation beyond Mash, extra comparators, wet-lab concordance). |
| §7 (new) | *Building your own sketch tool in the era of coding agents*: primitive catalogue, specification template, worked example 1 (deacon closed syncmers: decisions, patch, 107 tests, smoke numbers), worked example 2 (one predicate → four 2b-RAD tools), five failure modes met while building the harness, Box 4 release checklist. |
| §8 | Questions numbered Q1–Q9; Q8 (motif vs hash at matched density) and Q9 (computed vs sequenced sketch concordance) added; Q5/Q6 reconciled with what the pipeline measures. |
| §9 | Conclusions and limitations updated (third limitation = competing interest); *Competing interests* and *Data and code availability* sections added. |
| References | [116] 2b-RAD, [117] 2bRAD-M, [118] i2bRAD-M protocol, [119] Syn2bANI, [120] Fast2bRAD-M, [121] MAP2B placeholder, [122] Strain2bScan, [123] sk2bGrow, [124] Blanca et al. 2022, [125] REBASE, [126] this repository. |

## Must verify before submission

Items the reviser could not confirm from the sources at hand. None was
changed; each is flagged here rather than in the text.

Draft claims and citations

- [12] BinDash 2.0, [28] MaxGeomHash "accepted at RECOMB 2026", [29] dN/dS preprint, [30] Pilea (Microbiome 2026), [78] HyperSketch and [80] BioData Mining 2026, [84] ReSkmer (Genome Biol 2026), [85] Wu & Medvedev 2026: dates, venues, DOIs (several carry a 10.64898 prefix; confirm it is the current bioRxiv prefix).
- [25] sylph volume/pages (2025;43(8)); [48] ExaLogLog (EDBT 2025 pages); [52] SubseqSketch; [51] RabbitSketch; [71] GreedyMini; [98] EvANI; [97] Ponsero et al.
- §1.2 / §3.3: the Belbasi worked example "true Jaccard 0.90 → estimate 0.44" and the "up to 14% relative error in MashMap divergence" figure — check against [60].
- §2.3.3: "CMash deprecated in favour of sourmash and YACHT" [113] — check the repository README.
- §2.4.4, §5.1 Rule 5, §6.4: datasketch v2.0.0 (July 2026) affine32 change [101,114] — check the release notes.
- §6.2: RefinedWeb "9000 hashes as 450 buckets of 20" — the paper may state 20 buckets of 450; the LSH-curve argument in the same paragraph only holds for 450 × 20, so the two must agree.
- §6.4: SEDD "1.2-trillion-token corpus in ~2.9 h, up to 158× faster" [100].
- §2.5.3: UltraLogLog 17% and ExaLogLog 43% space figures [47,48].
- §3.1: the numerical consequences of equation (8) quoted from [19] (f = 1 at 10⁴ elements, ε = 0.07; errors beyond ±5% at scaled = 1000 and 10⁷ elements).

New material

- [117] 2bRAD-M author list and pages; [118] i2bRAD-M status and author list; [119], [122], [123] status (preprint DOIs if available); [121] whether MAP2B has a standalone citation; [125] REBASE year/volume.
- §2.3.4 and `i2brad.py`: BcgI cut offsets (10/12)…(12/10) and the other 14 enzymes against REBASE. The tag lengths reproduce Fast2bRAD-M's table for all 15, which is strong but indirect evidence.
- §7.3 numbers (800-line patch, 107 tests, 3,979/3,969 and 3,867/3,804 of 4,000) are from this repository's `env/deacon-syncmer.patch` and `tests/smoke_deacon_syncmer.sh`; rerun before quoting.
- §3.7 / §5.2: "0.990 ± 0.001 at 1% divergence from ~1,450 BcgI tags" is from six seeds in `tests/test_i2brad.py`-style simulation on random 3 Mb sequence; the benchmark's real genomes will replace it.

## Decisions for the authors

1. Title — the working title is long; alternatives: "Sketching by estimand: …" or "Choose the estimator, not the name: …".
2. §7 names no coding agent. Decide whether to name the tool(s) used and how to describe the human/agent division of labour; the section is written so either choice works.
3. Author list, affiliations, and the exact competing-interests wording (§9).
4. Journal. As written (21k words, 5 tables, 3 figures, 4 boxes) the piece is too long for any venue; see the cut list.
5. Whether Figure 1 is drawn with five panels (as the caption now says) — no figure files exist yet.

## Cut list to ~10k words

In order of least loss: (i) Table 4 and §2.10 to the supplement (−600); (ii) §2.2.2, §2.5.3, §2.6.2 to two sentences each (−700); (iii) §6 to one page, keeping §6.2's LSH-curve point and the datasketch reproducibility point (−900); (iv) §5.6 figure/table specifications to the supplement, leaving one paragraph (−900); (v) §4.3 merged into §3 as a closing paragraph (−500); (vi) §5.2–5.5 entrant lists into Table 5's Notes column (−600); (vii) §7.5 to a bulleted half page (−400); (viii) §1.4/1.5 halved (−300). That reaches ~16k; the last step to 10k is moving §5's task designs to a Methods/Supplement and leaving §5 as rules + placeholders (−4,000), which fits the "protocol is the deliverable" stance.

## Repository changes made for the manuscript

- `scripts/stage02_simulate/t3_design.py`: probability-Jaccard distance (`jp_true`) added to the T3 truth, so Q2's cross-scoring has a ground truth.
- `scripts/benchmark/metrics_T3.py`: Mantel correlation against `jp_true`.
- `scripts/tools/py/i2brad.py`: bases inside tags recorded for the sequencing-cost axis (Rule 7) — `.stat` sidecar in T2, header line in T3 sketches.
- `docs/tools/md2docx.py`: Markdown → .docx exporter used for `docs/manuscript.docx`.
