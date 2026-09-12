# Which quantity does your sketch estimate? A mechanism-first guide, benchmark protocol and build-your-own primer for k-mer sketching

Authors: [to be completed] — HuangShiLab, Faculty of Dentistry, The University of Hong Kong, and co-authors.

Correspondence: Shi Huang (huanglab07@gmail.com)

# Abstract

K-mer sketching has become the standard substrate for approximate genomic comparison, and the label "MinHash" now covers tools that estimate different quantities from the same data: Jaccard similarity, containment, two distinct weighted overlaps, an edit-distance proxy, or nothing at all. Choosing a tool by name is therefore unsafe; the documented failure cases are silent, returning a plausible number that answers a question other than the one posed.

This review organizes k-mer sketching by what each method estimates. Every sketch decomposes into a selection rule, which decides which k-mers are retained, and an estimator, which converts retained k-mers into a number. The two axes vary independently, and the estimator axis determines fitness for a task. We present a seven-family taxonomy along these axes and show that the context-free family — whose rule is a predicate evaluated on the k-mer itself — admits three predicates: a hash threshold (FracMinHash), a sub-hash position (syncmers), and a sequence motif (type IIB restriction tags, the basis of 2bRAD-M). The motif predicate is the only selection rule an enzyme can execute, so the sketch can be sequenced directly instead of computed from a full shotgun library; this adds a cost axis, bases sequenced to obtain the sketch, that no hash-based sketch can occupy.

Three artifacts carry the argument. A property matrix states each method's estimands and algebraic properties in falsifiable form. A task-eligibility matrix marks every family eligible, partial, or excluded for four common tasks — genome-to-genome average nucleotide identity, genome detection in a metagenome, metagenome-to-metagenome comparison, and read-level filtering — with a stated reason for every exclusion. A pre-registered, executable benchmark protocol specifies the data, two-tier ground truth, size-matched parameter sweeps and statistical plan needed to convert the matrix's untested cells into measurements. A closing primer addresses tool builders in the era of coding agents, where implementing a sketch has become cheap and specifying its estimand, guarantee and test oracle is the scarce skill; two worked examples — adding a syncmer scheme to a host-depletion tool, and a family of tools built on one restriction-tag predicate — show what the human must decide and what the agent can be left to do.

The practical shift is from choosing tools by name to choosing, and building, estimators by question.

# 1. Introduction

## 1.1 The sketch explosion in genomics

Sequencing archives grow faster than the compute budgets of the groups that analyze them, and whole-database questions — which genomes resemble this one, which species sit in this metagenome, which reads came from this reference — are increasingly answered without alignment. Sketching is the enabling compromise: replace each dataset with a small summary, built by a fixed rule, that supports approximate answers to a specific similarity query. By 2019 the approach had matured enough to support three surveys from different angles. Rowe's practical guide organized the field by algorithm family, pairing MinHash, HyperLogLog, Bloom filters, and Count-Min sketches with the genomic software built on each [92]. Marçais and colleagues took the theory view, covering compressed indexes, membership structures, locality-sensitive hashing, and minimizers as four routes to sublinear space [93]. Zielezinski and colleagues demonstrated the demand side directly, benchmarking 74 alignment-free comparison methods — most of them k-mer-based — across five application areas [94]. The adoption trigger in genomics was Mash, which showed that a few thousand hash values per genome suffice to estimate mutation distances across tens of thousands of genomes on commodity hardware [7]; the k-mer sketch became a default pre-filter, index, and sometimes the final answer.

The field has not stood still since those surveys. FracMinHash (preprint) [17] revived threshold-based sampling and made containment a first-class estimand, anchoring the sourmash ecosystem for metagenome search [21]. sylph paired scaled sampling with a zero-inflated Poisson model that debiases containment-based average nucleotide identity (ANI) at low coverage [25]. MaxGeomHash (preprint, accepted at RECOMB 2026) changed the size class itself, producing sketches that grow sublinearly with the set while remaining mergeable [28]. In parallel, a theory-correction line has hardened what was folklore: confidence intervals for FracMinHash-derived mutation rates [18], a safe scale-factor formula showing that widely used default thresholds can exceed their nominal error [19], and repeat-robust estimators that remove the bias induced by duplicated k-mers (preprint) [85]. The result is a larger, more diverse toolbox than any existing review describes.

## 1.2 The naming problem: one label, many estimators

Growth has outpaced terminology. The label "k-mer sketching" now covers tools whose estimators answer fundamentally different questions: Jaccard similarity and mutation distance (Mash [7]), containment of a small genome in a large mixture [15,16], two distinct weighted Jaccard quantities [37,38], edit-distance proxies [50], and — in the case of minimizer schemes — no similarity estimator at all, only an indexing structure [60]. Two tools that both advertise "MinHash" may therefore estimate different quantities from the same data, and the label alone says nothing about which question a tool can answer. Section 2 makes this precise by separating the two design axes that the naming convention collapses: the selection rule, which decides which k-mers a sketch retains, and the estimator, which converts retained k-mers into a number. The introduction needs only the observation that the two axes vary independently, and that the estimator axis determines fitness for a task.

The failures this produces are silent: the software runs, returns a number, and the number answers a different question than the user asked. Three documented cases frame the problem. First, estimator mismatch by size class. A FracMinHash sketch retains an expected fraction  $1/s$  of distinct k-mers [17]; a single short read contains only a few hundred distinct k-mers, so at the widely used default  $s=1000$  the read-level sketch is usually empty and a containment query returns zero hits. This is not a bug (the sketch correctly implements a threshold rule whose sampling density cannot support read-scale sets), but it defeats the intuition of users who expect "a MinHash" to work at any input size. Second, a biased estimator hiding behind a familiar construction. Belbasi and colleagues proved that the Jaccard estimator built on minimizer sketches is biased and inconsistent, because minimizer selection depends on k-mer layout along the sequence; in their worked example, a true Jaccard of 0.90 yields an estimate of only 0.44 [60]. Third, two estimands under one name. HULK advertises weighted-Jaccard ( $J_{w}$ ) estimation for microbiome profiles [39,40], but Ertl proved that the collision probability of its underlying HistoSketch algorithm equals the probability Jaccard  $J_{P}$  instead [38], a scale-invariant quantity that is provably at least as large as  $J_{w}$  on L1-normalized inputs. In each case the failure is a mismatch between the question and the estimator, not an implementation error, and in each case it went unflagged until someone checked the math.

## 1.3 Gap and scope

Existing reviews leave this mismatch axis unmapped. The 2019 surveys predate FracMinHash, sylph, and the entire theory-correction line [92,93,94]. Zheng, Marçais, and Kingsford wrote a guide restricted to that single family: minimizer sketches [95]. The nearest neighbor to the present work is Ndiaye and colleagues' 2024 review of minimizer, syncmer, and strobemer sketching across five application areas [96]; it touches on context-free selection properties, but organizes the field by method family rather than by estimand, includes no benchmark, and treats MinHash and FracMinHash only peripherally. On the benchmark side, EvANI compares ANI and evolutionary-distance estimators across tools [98], and Ponsero and colleagues compared twelve k-mer-based tools on metagenome sample comparison [97]; each is a genuine multi-tool evaluation, and each is confined to a single task and organized by tool rather than by estimator. A check of the 2023–2026 literature found no review organized by estimand, no estimator-by-task decision framework, and no cross-family benchmark spanning genome comparison, metagenome detection, sample comparison, and read filtering. No review has treated restriction-site tags — the reduced-representation sketches of 2b-RAD and 2bRAD-M [116,117] — as a member of the same design space, although, as §2.3 shows, they are a context-free sampling rule in exactly the sense that makes FracMinHash work, with the added property that the rule can be executed by an enzyme.

That matrix is the scope of this review. We cover k-mer sketching on nucleotide sequences, treating the estimator — what quantity a sketch provably approximates — as the primary axis and the biological task as the second; "MinHash-family" is used for the hash-based lineage, and the taxonomy admits any selection rule whose properties can be stated in the same terms. Protein-level and hyperdimensional sketches enter only where they change the estimator argument, and alignment-based methods appear only as ground-truth anchors. The claim is deliberately narrow and falsifiable: for each of four common tasks, exactly which sketch families produce a defensible estimate, and which fail silently.

## 1.4 Contributions

This review makes five contributions, each producing a concrete artifact.

First, a seven-family taxonomy of k-mer sketching organized by what each method estimates and how it selects k-mers, with properties stated precisely enough to be falsified: mergeability, intersection-commutativity, size class, context-freedom (§2; Tables 1 and 2). The taxonomy separates selection rules from estimators and shows that most "MinHash variants" differ on at least one of the two axes. Within the context-free family it identifies three interchangeable predicates — hash threshold, sub-hash position, and sequence motif — and places restriction-tag sketching in the design space for the first time.

Second, an estimator-by-task eligibility matrix (Table 3) covering four tasks — genome–genome ANI, genome detection in metagenomes, metagenome–metagenome comparison, and read-level filtering. Every exclusion or partial entry carries a stated reason, so the matrix functions as a decision procedure rather than a league table.

Third, an executable benchmarking protocol (§5) with a pre-specified statistical plan — cluster bootstrap over genomes and mixed models, because pairwise distances are not independent samples — size-matched parameter sweeps, a sequencing-cost axis, and version-locked software and data archives, implemented as a public pipeline [126] whose design was committed before any result existed.

Fourth, a cross-domain analysis (§6) of genomic sketching and large language model (LLM) training-data deduplication, two fields that descend from the same 1997 minwise-hashing construction [1], have converged on similar parameter choices independently, and barely cite each other.

Fifth, a build-your-own primer for the era of coding agents (§7): a specification template that turns the taxonomy's questions into a checkable request, two worked examples of agent-assisted tool development with the human decisions made explicit, the failure modes we encountered, and a release checklist (Box 4).

## 1.5 Organization

Section 2 develops the taxonomy and its property matrix. Section 3 collects the known failure modes of sketch estimators, pairing each with its theoretical boundary. Section 4 defines the four tasks and derives the eligibility matrix. Section 5 specifies the benchmarking protocol; Section 6 presents the cross-domain analysis; Section 7 is the primer for tool builders; Section 8 lists the open questions, most of which the benchmark is designed to settle, and Section 9 concludes.

# 2. A taxonomy of sketch selection rules

## 2.1 Formalism: a sketch is a selection rule plus an estimator

Every tool in this review runs the same pipeline. A genome, read set, or metagenome becomes a set  $A$  of  $n$  k-mers; a hash function  $h$  maps each k-mer into  $[0,H_{max}]$ ; a selection rule decides which hashed values to keep; an estimator converts one or two retained sets into a number. The classical targets are Jaccard similarity and the containment of  $A$  in  $B$ :

$J(A,B)=(|A∩B|)/(|A∪B|),  C(A,B)=(|A∩B|)/(|A|).$ 	(1)

Everything that follows varies these two design choices independently — what to keep, and what to compute from what was kept. Tools sharing the label "MinHash" differ on both axes, and the estimator axis, not the label, determines which question a tool can answer. Three properties locate any sketch in this taxonomy: the estimand (Jaccard, containment, mutation distance, weighted or probability Jaccard, an edit-distance proxy, cardinality, or nothing at all); the size class (fixed budget independent of  $n$ ; linear growth  $Θ(n/s)$ ; or sublinear growth); and whether selection is context-free — decided by a predicate on the k-mer alone — or context-dependent on neighboring k-mers, in Edgar's terminology [61]. Three context-free predicates appear in this review: a hash threshold (§2.3.2), the position of a minimal s-mer (syncmers, §2.7.3), and a sequence motif (restriction tags, §2.3.4); they share every algebraic property and differ in density control, uniformity, and in whether an enzyme can execute them. The third property decides which set-algebra operations survive sketching and separates estimation from indexing. Throughout,  $k$  denotes k-mer length; sketch sizes are denoted  $s$  or  $m$  where noted.

Figure 1 makes the selection axis concrete: one 14-element set under bottom- $k$ , a threshold rule, a register scheme, a window-minimum rule, and a motif rule. The five sketches disagree about which elements represent the set, and that disagreement, not implementation detail, propagates into downstream estimates.

Figure 1: One shared 14-element k-mer set under five selection rules; selected elements in dark blue. (a) Bottom- $k$  keeps a fixed count of smallest hashes. (b) A threshold (scaled) rule keeps all hashes below  $H_{max}/s$ ; sketch size tracks set size. (c) A register scheme keeps one extreme value per hash-space bin. (d) A window-minimum (minimizer) rule keeps the smallest hash per window of  $w=4$  consecutive k-mers, making selection depend on sequence context. (e) A motif rule keeps every k-mer containing a fixed recognition motif at a fixed offset (a type IIB restriction site); like (b) it is context-free, but its density is set by the motif rather than by a threshold.

Read across the panels, the point is divergence. Panels (a) and (b) select identical elements for small sets, but their selections diverge as  $n$  grows: the fixed budget of (a) thins its coverage of a large set, while (b) grows with the set and tracks it. Panel (c) retains extrema, not elements, so the set-algebra interpretation of (a) and (b) does not apply; its registers feed cardinality and likelihood-based estimators instead. Panel (d) is the outlier: its selection changes when the underlying sequence is edited, which is why window-minimum sketches serve indexing rather than estimation. Panel (e) selects like (b) — each element decides for itself — but the elements it selects are fixed by sequence, not by a random hash, so two genomes share a tag only where they share the motif and its flanks. Estimator design therefore begins with the selection rule.

Table 1 organizes the field into seven families along these axes. Two algebraic columns need defining. A sketch is mergeable if the sketch of a union can be computed from the sketches alone,  $S(A∪B)=f(S(A),S(B))$ , which underlies distributed sketching and database construction. It is intersection-commuting if  $S(A)∩S(B)=S(A∩B)$ , a strictly stronger property that makes intersections (and hence subtraction of a reference from a metagenome sketch) computable directly on sketches. Intersection-commuting holds only for context-free schemes — threshold sketches, syncmers, and motif tags — a fact that explains the gather and index-subtraction operations of §2.3.

Table 1. Seven sketch families defined by selection rule and estimator.

| Family | Selection rule | Estimand | Size class | Mergeable | ∩-commuting | Weights | Asymmetric sets | Representative methods |
|---|---|---|---|---|---|---|---|---|
| F1 Fixed-size MinHash (bottom- $k$ ) | $k$  smallest hashes | Jaccard; mutation distance | Fixed  $k$ | Yes (union-homomorphic) | No; Jaccard via union sketch | No | No | Broder (SEQUENCES 1997 [1]); KMV [2,3]; Mash (Genome Biol 2016 [7]); BinDash (Bioinformatics 2019 [11]) |
| F2 Threshold/scaled (context-free predicate) | All  $h(x)≤H_{max}/s$ ; or any predicate on  $x$  alone: minimal s-mer position (syncmer), recognition motif (restriction tag) | Containment; Jaccard; ANI | $Θ(n/s)$ ; motif: set by enzyme density | Yes | Yes | No | Yes (design target) | Modulo sketch [1]; containment MinHash (Appl Math Comput 2019 [15]); Mash Screen (Genome Biol 2019 [16]); FracMinHash (preprint [17]); sourmash [20,21]; 2b-RAD tags (Nat Methods 2012 [116]; 2bRAD-M, Genome Biol 2022 [117]) |
| F3 Weighted | Consistent weighted sampling | $J_{w}$  or  $J_{P}$  (distinct quantities) | Fixed  $m$ | Not standardized | No | Yes | Partial | CWS [32]; ICWS (ICDM 2010 [33]); BagMinHash (KDD 2018 [35]); ProbMinHash (IEEE TKDE 2022 [38]) |
| F4 Register/extreme-value | One extremum per hash-space bin | Cardinality; Jaccard via MLE | Fixed register array | Yes (element-wise max) | No | Via multiplicity-aware variants | Partial | HyperLogLog (AofA 2007 [5]); HyperMinHash (IEEE TKDE 2022 [42]); SetSketch (PVLDB 2021 [43]); Dashing (Genome Biol 2019 [44]) |
| F5 Order-sensitive | Minima of ordered k-mer subsequences | Edit-distance proxy (weighted Jaccard of positional multisets) | Fixed,  $∝mℓ$ | No (per sequence) | No | Implicit (multiplicities) | No | OrderMinHash (Bioinformatics 2019 [50]); tensor sketch (preprint [53]); SubseqSketch (WABI 2025 [52]) |
| F6 Sampling schemes (indexing) | Minimum hash per window of  $w$  k-mers | None — indexing seeds; Jaccard estimator biased and inconsistent [60] | $Θ(2n/(w+1))$ | Partial (context-dependent) | No for minimizers; yes for syncmers | No | n/a | Minimizers [58,59]; syncmers (PeerJ 2021 [61]); strobemers (Genome Res 2021 [63]); mod-minimizer (WABI 2024 [67]) |
| F7 Vector/hyperdimensional | No element selection; random projection and superposition | Inner product → Jaccard → ANI | Fixed dimension  $D$ | Additive (caveat duplicates) | No | No | Partial | SimHash (STOC 2002 [4]); DotHash (KDD 2023 [76]); HyperGen (Bioinformatics 2024 [77]) |

The table is a claim about priorities, not a league table. The estimand column is the primary sort key: tools in different rows answer different questions even when their k-mer plumbing is identical. The size-class column dictates feasibility — a fixed- $k$  sketch behaves identically on a bacterial genome and a 100 Gbp metagenome, while a threshold sketch does not. The mergeability and intersection columns dictate which distributed and subtractive workflows are expressible at all. The F6 row is deliberately deflationary: minimizer schemes dominate indexing practice, but theirs is the only row whose estimand reads "none," anticipating the bias results of §2.7.

## 2.2 F1: Fixed-size MinHash (bottom-k)

### 2.2.1 Lineage: from document resemblance to genomic distance

The family originates with Broder's resemblance sketches for near-duplicate document detection, which showed that a fixed-size sample of minimal hash values suffices to estimate  $J$  [1]. The bottom- $k$  formulation (retain exactly the  $k$  smallest hash values; the notation is historical, with  $k$  here a sketch size distinct from k-mer length) was analyzed as k-minimum-values (KMV) sampling in the streaming literature [2] and given its modern estimation theory by Cohen and Kaplan [3]. Mash carried the construction into genomics, pairing it with a mutation model: under a Poisson model of random, independent point mutation acting on unique k-mers, a fraction  $e^{−kd}$  of shared k-mers survives divergence  $d$ , and solving for  $d$  in terms of the Jaccard estimate  $j$  yields the Mash distance [7]

$D=−(1)/(k)ln(2j)/(1+j),$ 	(2)

with  $D=0$  when  $j=1$  and  $D$  diverging as  $j→0$ . The Poisson assumption is a modeling choice, not a theorem (repeats, indels, and non-uniform mutation violate it), but nearly every sketch-based average nucleotide identity (ANI) number passes through equation (2) or a close relative.

### 2.2.2 Bit-width compression

A second line of work shrinks each retained hash rather than their number. b-bit MinHash keeps only the lowest  $b$  bits per value, accepting a small estimator correction for large space savings [8,9]. One permutation hashing (OPH) replaces  $k$  independent permutations with one permutation plus binning, cutting sketching cost toward  $O(n)$  at the price of empty bins requiring densification [10]. BinDash ported the b-bit-plus-OPH stack to rolling k-mer hashes for laptop-scale genome comparison [11]; BinDash 2.0 (preprint) [12] adds SIMD densification for trillion-scale search. These methods inherit the bottom- $k$  estimand and differ only in encoding.

### 2.2.3 Exact property statement

Bottom- $k$  sketches are homomorphic under union: with a shared hash function,  $S(A∪B)=bottom-k(S(A)∪S(B))$ , so sketches merge losslessly [1]. They are not homomorphic under intersection: the intersection of two sketches is not the sketch of the intersection, and  $|A∩B|$  cannot be read off directly. Jaccard must instead be estimated through the union sketch, as Mash writes explicitly [7]:

$J(A,B)≈(|S(A∪B)∩S(A)∩S(B)|)/(|S(A∪B)|).$ 	(3)

With the merged sketch of size  $k$  as denominator, equation (3) degrades exactly when  $|A∪B|≫k$ : fixed-size sketches price their error in units of union size, which is why asymmetric comparisons defeat them — a theme revisited in §3.

### 2.2.4 Application side

Two applications show the family's range. finch wraps Mash-style sketches in dynamic abundance filtering, excluding low-coverage k-mers before sketching and extending bottom- $k$  to raw read sets [13]. PopPUNK fits Jaccard estimates at multiple k values and decomposes them into core- and accessory-genome components, turning a scalar estimand into a two-dimensional strain-relationship score [14]. Both leave the estimator untouched and innovate at the selection and modeling layers.

## 2.3 F2: Threshold/scaled sketches

### 2.3.1 The modulo legacy and the containment index

Broder's 1997 paper already contained a second scheme beside bottom- $k$ : the modulo sketch, which retains every hash congruent to  $0 mod M$  and grows linearly with the set. Broder himself cautioned that estimating the containment of a small set in a much larger one is "rather error prone due to the paucity of samples," and the idea lay dormant in genomics for two decades [1]. Koslicki and Zabeti revived the asymmetric problem from the theory side, showing that MinHash's relative error for containment explodes when  $|A|≪|B|$  and deriving exact error formulas for a containment-index estimator; tellingly, the analysis appeared in an applied-mathematics journal [15], a venue gap that recurs in §2.10. Mash Screen engineered the asymmetric design directly: only the reference  $A$  is sketched, and the query  $B$  is streamed through it, giving [16]

$c_{k}(A,B)≈(|S(A)∩π(B)|)/(|S(A)|),$ 	(4)

where  $π(B)$  denotes the hashes of all k-mers of  $B$  (no sketch of  $B$  is needed). This is an unbiased estimate of the containment index with expected error  $O(1/sqrt(m))$  for sketch size  $m$ . Two Mash Screen outputs must not be conflated: the containment index of equation (4) is a set-overlap fraction, whereas the containment score  $M_{c}(A,B)≈c_{k}^{ 1/k}$  is a derived per-nucleotide identity estimate under a k-mer survival model [16]; only the second carries modeling assumptions.

### 2.3.2 FracMinHash and the intersection-commuting property

FracMinHash (preprint) [17] reinstates the modulo idea with a MinHash-flavored threshold. For a scale factor  $s$ ,

$FRAC_{s}(A)={ x∈A:h(x)≤H_{max}/s },$ 	(5)

so the sketch retains an expected fraction  $1/s$  of distinct k-mers, with size floating as  $Θ(n/s)$  — "a mix of MinHash and ModHash" in the authors' words [17]. Its decisive property is intersection-commuting: because retention depends only on each element's own hash,  $FRAC_{s}(A)∩FRAC_{s}(B)=FRAC_{s}(A∩B)$  exactly. This makes the sourmash ecosystem's signature operations (gather, which greedily subtracts the best-matching reference from a metagenome sketch, and index subtraction generally) well-defined set algebra rather than heuristic [20,21]. The infrastructure spans the Rust-core sourmash v4 [21] and the branchwater plugin ecosystem for petabyte-scale search (preprint) [22]. On the theory side, Rahman Hera and Koslicki have formalized which measures admit sound FracMinHash estimation and under what scale factors [19]; their safe-scale-factor critique of fixed default thresholds is taken up in §3.

### 2.3.3 The tool cluster

The F2 application cluster is large and growing. CMash estimates Jaccard and containment at multiple k from a single truncated sketch [23]; its authors now deprecate it in favor of sourmash and YACHT.[113] Metalign uses containment MinHash only as a pre-filter before alignment-based profiling [24]. YACHT recasts detection as a hypothesis test, returning a presence/absence decision with confidence rather than an abundance point estimate [26]. fmh-funprofiler applies FracMinHash protein sketches to functional profiling [27]; a preprint extends the same machinery to alignment-free dN/dS estimation [29]; Pilea uses FracMinHash containment to profile bacterial growth dynamics [30]. Two members sit at the family's edge. sylph is F2-adjacent rather than a FracMinHash derivative, pairing FracMinHash-style scaled sampling with its own abundance-corrected sketch and a zero-inflated Poisson model that debiases containment ANI at low coverage [25]. MaxGeomHash (preprint, accepted at RECOMB 2026) generalizes the size class itself, producing variable-size sketches of  $blg(n/b)+O(b)$  elements (sublinear in  $n$ , between MinHash's constant and FracMinHash's linear footprint) while remaining mergeable and asymptotically unbiased [28].

### 2.3.4 Motif-defined sampling: restriction tags as a context-free predicate

Threshold sketches owe their intersection-commuting property (§2.3.2) to a single fact: whether a k-mer is retained depends on the k-mer alone. Syncmers share the property through a different predicate — the position of the minimal s-mer inside the k-mer (§2.7.3) — and neither predicate has anything to do with hashing as such. A third predicate has been in the wet lab for over a decade. Type IIB restriction enzymes cut on both sides of their recognition site and release iso-length fragments; 2b-RAD used these fragments as genome-wide genotyping tags [116], and 2bRAD-M sequenced only them to profile low-biomass and degraded microbiomes from roughly 1% of the bases of a shotgun library [117,118]. Computationally, a tag is a k-mer selected if and only if it contains the recognition motif at a fixed offset. For BcgI (CGA·N6·TGC, REBASE cut offsets (10/12)…(12/10) [125]) the double-stranded core is a 32-mer; applying the same rule in silico to a reference — an "i2bRAD" digestion — yields the reference's tag set, and applying it to reads yields the sample's. The rule is context-free, so tag sets commute with intersection exactly as FracMinHash sketches do, and the same k-mer is selected in every genome and every read that contains it. Figure 1(e) adds the motif rule to the selection gallery.

Three properties separate the motif predicate from the hash predicate. (i) Density is fixed by the enzyme and tunable only discretely — by enzyme choice or by a multi-enzyme panel — where the scale factor of equation (5) is continuous. For a site with six specified bases the expected density on random sequence is  $2/4^{6}$  per position, about one tag per 2 kb: the BcgI tag set of a 5 Mb genome has roughly 2,400 members, the size of a FracMinHash sketch at  $s≈2000$  (on random sequence our sampler measures a scaled-equivalent of 2,100; §5.2). (ii) Sampling is not uniform over k-mer space: motif frequency tracks GC content and dinucleotide usage, so tags per megabase vary across taxa and some genomes are tag-poor. The consequence is a failure mode different from hashing — a no-estimate rather than a wrong estimate — and estimators derived under uniform sampling must be re-derived (§3.7). (iii) The predicate can be executed by an enzyme. Every other selection rule in Table 1 is computed after sequencing, so a sketch's byte cost is paid on top of a full library; a restriction digest performs the selection before sequencing, so the library is the sketch. This adds a cost axis that no hash-based sketch can occupy — bases sequenced to obtain the sketch — and it is why the family's tools were built for specimens shotgun sequencing handles poorly: formalin-fixed sections, picogram inputs, and samples with more than 99% host DNA [117,118].

The estimator changes with the predicate. Under the substitution model of §2.2, a tag of length  $k$  survives divergence  $d$  with probability  $(1−d)^{k}$  whether it was selected by hash or by motif, so the containment estimator  $C(A,B)=|A∩B|/|A|$  and its transform  $d=1−C^{1/k}$  remain unbiased; the Jaccard route of equation (2) does not, because a substitution inside the recognition site deletes the tag from the mutated genome instead of replacing it with a freshly sampled k-mer (§3.7). The tool cluster follows the tasks of §4. MAP2B and its Rust port Fast2bRAD-M profile metagenomes from tags, deciding presence on species-specific tag subsets with a breadth-times-depth score (the G-score) and a learned false-positive filter [117,118,120]; Syn2bANI estimates ANI between genomes from tag anchors — because tags occupy fixed positions, matches chain without a seeding step, and aligned fraction and structural-variant calls come out of the same pass [119]; Strain2bScan and sk2bGrow carry the same tag sets to strain-level and growth-dynamics estimands [122,123]. In this review the motif predicate enters the benchmark twice: as a tool-free mechanism row (in-silico tags scored with containment and Bray–Curtis estimators, so that the sampling rule is measured apart from any tool) and as the tool rows. We state the competing interest here as well as in §9: several authors develop these tools, and the benchmark design was fixed before any tag-based result existed.

### 2.3.5 An adjacent route

kssd represents a neighboring design: it samples a subspace of k-mer substring space rather than thresholding hash values, and supports set operations on its sketches [31]. It is not a MinHash-family method, but it competes for the same containment-and-distance questions and is a natural benchmark comparator.

## 2.4 F3: Weighted sketches —  $J_{w}$  is not  $J_{P}$

### 2.4.1 Lineage

When k-mer multiplicities matter (metagenome abundance profiles being the canonical case), the estimand shifts from set overlap to distributional overlap, and a separate lineage applies. Consistent weighted sampling (CWS) extended minwise hashing to real-valued weight vectors [32]; ICWS made sampling worst-case constant time per non-zero weight [33]; 0-bit CWS compressed signatures to single bits [34]; BagMinHash removed the remaining speed bottleneck [35]; DartMinHash (preprint) cut costs further for sparse vectors [36]. In parallel, Moulton and Jiang defined a second generalization, the probability Jaccard  $J_{P}$ , with a matching sampler [37], and Ertl's ProbMinHash provided a full locality-sensitive hash family for it [38].

### 2.4.2 Two quantities called "weighted Jaccard"

The two estimands are different functions. For non-negative weight vectors  $w_{A},w_{B}$  over a domain  $D$ ,

$J_{w}(A,B)=(∑_{d∈D}^{​} min(w_{A}(d),w_{B}(d)))/(∑_{d∈D}^{​} max(w_{A}(d),w_{B}(d))),$ 	(6)

$J_{P}(A,B)=∑_{d∈D}^{​} (∑_{d′∈D}^{​} max​((w_{A}(d′))/(w_{A}(d)),(w_{B}(d′))/(w_{B}(d))))^{−1}.$ 	(7)

The practical difference is scale invariance:  $J_{P}$  is unchanged when either vector is rescaled, treating inputs as probability distributions, while  $J_{w}$  treats them as weighted sets whose absolute mass matters [37,38]. The two are bracketed tightly: under L1 normalization,  $J_{P}∈[J_{w}, 2J_{w}/(1+J_{w})]$  [37], so they agree at high similarity and diverge most in the mid-range — where metagenomic profile comparisons live. A benchmark reporting "weighted Jaccard" without specifying which quantity is reporting an ambiguous number.

### 2.4.3 The HULK boundary case

HistoSketch was proposed for streaming histograms under a weighted-Jaccard interpretation [39], and HULK applied it to microbiome k-mer spectra [40]. Ertl then proved that the final HistoSketch algorithm's collision probability equals  $J_{P}$ , not the  $J_{w}$  its design intended [38]. A widely used microbiome tool therefore estimates a different quantity than its documentation states. GSearch makes the distinction load-bearing by running both: its ProbMinHash module estimates  $J_{P}$ , which better discriminates closely related genomes, while its SetMinHash module (SuperMinHash, densified one-permutation hashing, or SetSketch) estimates unweighted  $J$  [41].

### 2.4.4 Merge operations are not standardized

One property remains genuinely open. Unweighted MinHash signatures merge by component-wise minimum, but the weighted literature standardizes no merge operation: the CWS, BagMinHash, and ProbMinHash papers define none, and datasketch's WeightedMinHash exposes none [101]. Whether a weighted signature can be merged depends on the intended union semantics (max-union of weights versus additive bag-union), and the field has not settled it. Version drift compounds this: datasketch v2.0.0 (July 2026) changed the default MinHash permutation scheme to affine32, repairing an overestimation bias on large sets but breaking hash compatibility with earlier versions [101].[114] Weighted-sketch results should be reported with implementation and version pinned.

## 2.5 F4: Register/extreme-value sketches

### 2.5.1 From cardinality to joint estimation

HyperLogLog (HLL) replaced retained samples with a fixed register array, each register recording the most extreme hash mapped to its bin, estimating cardinality with standard error about  $1.04/sqrt(m)$  in 1.5 KB [5]. HyperMinHash layered Jaccard estimation onto an HLL scaffold, compressing MinHash buckets from  $O(logn)$  to  $O(loglogn)$  bits while retaining streaming updates and unions [42]; SetSketch unified the two ends, interpolating between MinHash-like and HLL-like behavior via an adjustable base, with clean maximum-likelihood machinery [43]. The price of compactness is paid at estimation time: there is no closed-form collision probability as for MinHash, inclusion–exclusion approximations perform poorly, and the best joint estimates come from maximum-likelihood estimation (Ertl's HLL-MLE, preprint) [6,43].

### 2.5.2 Genomics adoption

Dashing brought the register family into genomics, pairing HLL with Ertl's MLE and joint-MLE estimators to sketch and distance over 87,000 genomes in minutes [44]. Dashing 2 rebuilt the system on SetSketch, reclaiming the bits that HLL's leading-zero-count wastes, added multiplicity-aware sketching via ProbMinHash, and integrated locality-sensitive hashing for all-pairs comparison at million-genome scale [45]. ntCard marks the family's boundary: a streaming sampler that reconstructs the k-mer frequency histogram, and hence cardinality, but offers no set-similarity estimand [46]. It belongs to the family by construction while answering a different question — a compact case for sorting by estimand rather than mechanism.

### 2.5.3 Theory successors not yet in genomics

The register line has continued in the theory literature: UltraLogLog reduces space by 17% at equal error [47], ExaLogLog extends distinct counting to the exa-scale at 43% less space [48], and HyperLogLogLog compresses the register distribution itself [49]. None has yet been adopted by a genomics tool. The lag fits the family's history (twelve years separated HLL from Dashing) and represents headroom for the next Dashing rather than a theory deficiency.

## 2.6 F5: Order-sensitive sketches

### 2.6.1 OrderMinHash

OrderMinHash (OMH) is a gapped locality-sensitive hash family for the edit distance: each sketch entry concatenates  $ℓ$  minimum-hash k-mers in textual order, making the sketch sensitive to k-mer content and relative order alike, with the weighted Jaccard of positional k-mer multisets as the proxy estimand rather than the edit distance itself [50]. Encoding order imposes three constraints: the sketch is defined per sequence and cannot be merged across a multi-sequence file; canonical k-mers are forbidden (reverse-complement collapsing scrambles positional information); and each entry costs roughly  $ℓ$  times a MinHash entry. The authors acknowledge a large gap in the theoretical guarantee, a criticism sharpened by McCauley (preprint): OMH bounds neither collision threshold in the worst case and so yields no similarity-search guarantee [57].

### 2.6.2 Implementations and new members

RabbitSketch is engineering substrate rather than a new method: a multithreaded library implementing MinHash, OMH, kssd, and HLL, whose parallel OMH improves on the reference implementation by about two orders of magnitude [51]. Two recent methods extend the family. SubseqSketch sketches random subsequences, and its sketches' cosine similarity tracks edit similarity more closely than OMH (Pearson 0.918 versus 0.726 in their evaluation) [52]. Tensor sketching (preprint) embeds edit distance into an  $ℓ_{2}$  metric via tensor products [53], and its long-seed-sketch application improves alignment of distant sequences to graphs [54]; the theoretical underpinnings trace to streaming edit-distance embeddings [55] and locality-sensitive bucketing functions [56]. F5 remains the smallest family in genomics use, in part because its mergeability and canonicalization restrictions clash with standard k-mer pipelines.

## 2.7 F6: Sampling schemes for indexing — not estimators

### 2.7.1 Dual origin and defining properties

The minimizer was proposed independently twice: as winnowing for document fingerprinting, with a proved asymptotic density of  $2/(w+1)$  and a  $1.5/(w+1)$  lower bound for local algorithms [59], and as minimizers for biological sequence comparison, with the same density and the window guarantee that any shared substring of length  $w+k−1$  shares a minimizer [58]. Selection is context-dependent in Edgar's terminology: retention depends on a k-mer's neighbors within the window, not on the k-mer alone [61]. That is what makes minimizers well suited to seeding (the window guarantee) and unusable as similarity sketches.

### 2.7.2 The Belbasi theorem

Belbasi et al. proved that the minimizer-based Jaccard estimator is biased and inconsistent: in their worked example, a true Jaccard of 0.90 yields an estimate of about 0.44 [60]. MashMap3's minmers repair the estimator, generalizing minimizers to a scheme with unbiased local Jaccard estimation [66]. The mechanism and magnitude of the bias are treated in §3.3.

### 2.7.3 The scheme family

The scheme family has been redesigned repeatedly. Syncmers select k-mers whose smallest contained s-mer sits at a fixed position, making selection context-free and commutative with set intersection, the property otherwise unique to threshold and motif sketches (§2.3), and parameterized syncmers (preprint) [62] improve long-read mapping [61]. Strobemers link several short k-mers spaced across a window, tolerating substitutions and small indels [63]. On density, the mod-minimizer approaches the asymptotic optimum  $1/w$  for long k-mers [67], its open-closed extension covers small  $k$  [68], and a near-tight lower bound delimits what any forward scheme can achieve [69]. Engineering keeps pace: SIMD implementations reach about 500 Mbp/s [70], and greedy constructions approach the lower bounds at the cost of storing explicit orders [71].

### 2.7.4 Tools

The family's tools are indexers and aligners, not estimators. minimap2 seeds alignment with a minimizer index [72]; Kraken 2 stores minimizers with spaced-seed masking in a compact hash table [73]; strobealign aligns reads with flexible-size randstrobe seeds [74]; and deacon (preprint) queries a pangenome minimizer index for host depletion at over 250 Mbp/s on a laptop [75]. Each exploits the window guarantee; none reports a Jaccard estimate, and by the Belbasi result none should.

## 2.8 F7: Vector/hyperdimensional sketches

The final family abandons element selection altogether. SimHash signs random hyperplane projections so that bit-signature collisions encode cosine similarity [4], a construction engineered at web scale by Manku et al. [79]. DotHash replaced signs with superposition of random high-dimensional vectors, whose inner product unbiasedly estimates  $|A∩B|$  [76]. HyperGen carries this into genomics: k-mers are FracMinHash-sampled, encoded as quasi-orthogonal hypervectors, and summed into a fixed-width sketch; estimation then runs inner product to intersection size to Jaccard to ANI via equation (2), with accuracy claims confined to ANI above 85% [77]. The chain matters: the sketch's native currency is an inner product, and the ANI number inherits every conversion assumption. HyperSketch (preprint) encodes de Bruijn graph topology rather than flat k-mer sets [78], and a vector-symbolic architecture has been applied to viral pangenome classification [80].

## 2.9 A property matrix of representative methods

Table 2 compresses the taxonomy into a method-level matrix over nine properties: four estimand questions (Jaccard, containment, distance/ANI, abundance), cardinality, and four operational properties (mergeability, intersection-commuting, streaming construction, asymmetric-set support). Cells marked "partial," "claimed," or "untested" mark thin published evidence — the cells the benchmark of §5 is designed to fill with measurements.

Table 2. Property matrix for representative methods. yes/no from cited sources; partial = restricted scope; claimed = asserted by authors without independent verification; untested = no published evidence located.

| Method (family) | Jaccard | Containment | Distance/ANI | Abundance | Cardinality | Mergeable | ∩-commuting | Streaming | Asymmetric sets |
|---|---|---|---|---|---|---|---|---|---|
| Mash (F1) [7] | yes | no | yes | no | no | yes | no | yes | no |
| BinDash (F1) [11] | yes | no | yes | no | no | untested | no | yes | no |
| Mash Screen (F2) [16] | no | yes | yes (score) | no | no | yes | no (reference sketch is fixed-size MinHash, not scaled) | yes (query side) | yes |
| sourmash / FracMinHash (F2) [17,21] | yes | yes | yes [19] | partial | partial | yes | yes | yes | yes |
| sylph (F2-adjacent) [25] | no | yes | yes | yes (debiased) | no | untested | untested | yes | yes |
| in-silico 2b-RAD tags (F2, motif predicate) [116,117] | yes (biased conversion, §3.7) | yes | yes (containment route) | yes (tag counts) | no | yes | yes | yes | yes |
| Syn2bANI (F2 motif + chaining) [119] | no | partial (screen) | yes (MLE on chained anchors, with AF) | no | no | yes (tag database) | yes | no | partial |
| Fast2bRAD-M / MAP2B (F2 motif) [117,118,120] | no | yes (species-specific tags) | no | yes (G-score, relative abundance) | no | yes | yes | yes | yes |
| YACHT (F2) [26] | no | partial (test only) | yes (ANI-based test) | no | no | untested | yes | yes | yes |
| Dashing 2 (F4) [45] | yes (MLE) | partial | yes | yes (multiplicity) | yes | yes | no | yes | partial |
| GSearch (F3/F4) [41] | yes | partial | yes | yes ( $J_{P}$  path) | no | untested | no | yes | partial |
| HULK (F3, boundary) [40] | no | no | no | yes (estimates  $J_{P}$  [38]) | no | no (no merge defined) | no | yes | no |
| OrderMinHash (F5) [50] | no | no | partial (edit proxy) | implicit | no | no | no | no | no |
| HyperGen (F7) [77] | yes | untested | yes (ANI > 85%) | no | partial | partial (duplicate caveat) | no | yes | untested |
| minimizer schemes (F6) [58,61] | biased [60] | biased [60] | no | no | no | partial | no (yes for syncmers) | yes | n/a |

Three readings matter for what follows. First, the estimand columns separate cleanly by family: no tool earns a "yes" in both the Jaccard and abundance columns without running two distinct estimators, as GSearch does explicitly and as the tag row does (presence Jaccard and count Bray–Curtis are two estimators over one sketch). Second, the operational columns cluster as Table 1 predicts (the F2 rows, hash and motif alike, and the syncmer corner of F6 are the intersection-commuting entries) because these are properties of the selection rule, not of the implementation. Third, the "untested" cells are not rhetorical gaps: for several widely used tools, mergeability and asymmetric-set behavior have never been measured under controlled conditions, and sylph's hybrid design inherits properties from scaled sampling only by analogy until tested. The benchmark section treats every "claimed" and "untested" cell as a hypothesis.

## 2.10 Timeline: theory precedes application, with a shrinking lag

Table 4 lays out the two publication streams — sketch theory and genomics application — from 1997 to 2026.

Table 4. Theory versus genomics application, 1997–2026.

| Year | Theory | Genomics application |
|---|---|---|
| 1997 | MinHash and modulo sketch (SEQUENCES) [1] | — |
| 2002 | KMV bottom- $k$  (RANDOM) [2]; SimHash (STOC) [4] | — |
| 2003 | Winnowing (SIGMOD) [59] | — |
| 2004 | — | Minimizers (Bioinformatics) [58] |
| 2007 | HyperLogLog (AofA) [5]; SimHash at web scale (WWW) [79] | — |
| 2010 | CWS [32]; ICWS (ICDM) [33]; b-bit MinHash (WWW) [8] | — |
| 2012 | One permutation hashing (NeurIPS) [10] | 2b-RAD restriction-tag genotyping (Nat Methods) [116] |
| 2015 | 0-bit CWS (KDD) [34] | — |
| 2016 | — | Mash (Genome Biol) [7]; sourmash (JOSS) [20] |
| 2017 | HistoSketch (ICDM) [39]; HLL maximum-likelihood estimators (preprint) [6] | ntCard (Bioinformatics) [46]; MashMap (RECOMB) [64] |
| 2018 | BagMinHash (KDD) [35]; probability Jaccard  $J_{P}$  (ICDM) [37] | minimap2 [72]; FastANI [82]; finch (JOSS) [13] |
| 2019 | Containment index (Appl Math Comput) [15]; OMH (Bioinformatics/ISMB) [50] | Mash Screen [16]; Dashing [44]; BinDash [11]; PopPUNK [14]; HULK [40]; Kraken 2 [73] |
| 2020 | ProbMinHash early access [38]; DartMinHash (preprint) [36] | Metalign (Genome Biol) [24] |
| 2021 | SetSketch (PVLDB) [43] | Syncmers (PeerJ) [61]; strobemers (Genome Res) [63]; kssd (Genome Biol) [31] |
| 2022 | ProbMinHash and HyperMinHash journal versions (IEEE TKDE) [38,42]; HyperLogLogLog (KDD) [49]; minimizer Jaccard bias proof (Bioinformatics/ISMB) [60] | FracMinHash and gather (preprint) [17]; CMash [23]; strobealign [74]; branchwater (preprint) [22]; 2bRAD-M tag-based profiling (Genome Biol) [117] |
| 2023 | DotHash (KDD) [76]; FracMinHash confidence-interval theory (Genome Res) [18] | Dashing 2 [45]; MashMap3 minmers [66]; long seed sketches [54]; skani (Nat Methods) [81] |
| 2024 | Mod-minimizer (WABI) [67]; UltraLogLog (PVLDB) [47]; forward-scheme density lower bound [69] | sylph (Nat Biotechnol) [25]; YACHT [26]; GSearch (NAR) [41]; HyperGen [77]; sourmash v4 [21]; fmh-funprofiler [27]; BinDash 2.0 (preprint) [12] |
| 2025 | Safe scale factors for FracMinHash (Algorithms Mol Biol) [19]; open-closed mod-minimizer [68]; ExaLogLog (EDBT) [48]; SubseqSketch (WABI) [52]; MaxGeomHash (preprint) [28] | RabbitSketch [51]; deacon (preprint) [75]; SimdMinimizers (SEA) [70]; sketch-based dN/dS (preprint) [29] |
| 2026 | MaxGeomHash accepted at RECOMB 2026 [28] | Pilea (Microbiome) [30]; HyperSketch (preprint) [78]; vector-symbolic pangenome classification [80]; i2bRAD-M protocol, Syn2bANI, Fast2bRAD-M [118–120] |

Nineteen years separate MinHash from Mash; twelve separate HyperLogLog from Dashing; nine separate consistent weighted sampling from HULK. Even the minimizer — born twice within a year on both sides — waited until 2022 for its misuse as a Jaccard estimator to be ruled out [60] and 2023 for the repair [66]. Two forces now compress the lag. Venue mixing has increased: containment theory reached genomics through an applied-mathematics journal [15], and the FracMinHash line moved from a 2022 preprint to confidence-interval theory [18], a scale-factor framework [19], and a RECOMB-accepted generalization [28] in four years. The application column has become theory-literate: Dashing 2 adopted SetSketch within two years of publication [43,45]. The 2024–2026 rows are the densest on both sides; the organizing question — which quantity a sketch estimates — now decides correctness.

# 3. When sketch estimators silently fail

Section 2 organized the field by what each family estimates and what its selection rule preserves. None of those properties guarantees that a particular estimate is trustworthy. The failure modes collected here share a signature: the software runs to completion and returns a plausible number, and the error is visible only from outside — from theory, from simulation, or from comparison against alignment. Each subsection isolates one mechanism and its consequence; where the empirical magnitude of a failure has not been measured, we say so and defer the measurement to the protocol of §5.

## 3.1 Scale factor and sampling sufficiency

Threshold sketches (family F2) size themselves to the data, which moves the failure point from a fixed sketch budget to the scale factor. Rahman Hera and Koslicki derived the conditions under which FracMinHash estimation is sound: for similarity and distance measures expressible through intersection and set sizes — cosine, Bray–Curtis, Kulczynski-1, and Sørensen among them — a retained fraction  $f=1/s∈(0,1]$  achieves relative error  $ε$  with confidence  $1−α$  only if [19]

$f ≥ (3 (2+ε)^{2} ln​(6/α))/(ε^{2} min{|A|,|B|}).$ 	(8)

The bound prices accuracy in units of the smaller set, complementing the fixed-size case where equation (3) prices error in units of union size. Its consequences are uncomfortable for current defaults. With sets of  $10^{4}$  distinct elements, demanding  $ε=0.07$  at  $α=0.05$  makes the right-hand side of (8) exceed 1, so the theory requires  $f=1$  — that is, no subsampling at all [19]. The same paper shows that a preset scale factor can produce errors larger than users expect: at the sourmash default scaled=1000 ( $f=10^{−3}$ ) and  $10^{7}$  distinct elements, many estimates fall outside  $±5%$  of the true value [19,21]. The scale factor is therefore an accuracy parameter set by the smaller input and the tolerated error, not a compression convenience to be fixed once and forgotten; both the chosen value and the set sizes should travel with any reported estimate.

## 3.2 Confidence intervals in practice

A point estimate without an interval invites overinterpretation, yet interval machinery is recent and thinly validated. Hera, Pierce-Ward, and Koslicki built an early such theory for FracMinHash: a debiased containment estimator  $C$  that is asymptotically normal, converted to a mutation rate by the delta method [18],

$d=1−C^{ 1/k},  se(d)=(1)/(k) C^{ 1/k−1} se(C),$ 	(9)

with  $k$  the k-mer length; the construction holds across a wider range of evolutionary distances than earlier fixed-sketch analyses [18]. sylph attaches confidence intervals to its containment ANI from the sampling variance of its zero-inflated Poisson model, which debiases the estimate at low sequencing coverage [25]. Both intervals are nominal: they inherit asymptotic normality and the correctness of the underlying k-mer survival model. Whether a nominal 95% interval from either pipeline actually covers the truth 95% of the time across realistic genomes — with repeats, indels, compositional heterogeneity, and shallow coverage — has not been measured by any benchmark we are aware of. Realized coverage, not only point accuracy, is a scoring axis that §5 builds into the benchmarking protocol.

## 3.3 Sequence-layout dependence

Sampling schemes (family F6) fail through a different route: the retained set depends on where k-mers sit, not only on which k-mers exist. Belbasi and colleagues proved that the minimizer Jaccard estimator is both biased and inconsistent — the expected gap between estimate and truth does not vanish even as sequence length grows — because the estimator depends on the layout of shared k-mers along the two sequences, while the Jaccard index depends only on how many k-mers are shared [60]. The distortion is not subtle: in their worked example, a true Jaccard of 0.90 yields an estimate near 0.44, and the bias propagates into MashMap's reported divergence with errors of up to 14% relative (about half a percentage point of divergence in absolute terms) [60]. The repair came from within the same lineage. MashMap3 replaces minimizers with minmers — a generalization combining rolling MinHash with multiple retained samples per window — and indexing minmers rather than minimizers removes the Jaccard estimator bias present in earlier MashMap versions [66]. Two practical lessons follow. For indexing, the bias is harmless, since seed selection needs no unbiased estimand; it bites only when an indexing sketch is repurposed as a similarity estimator. And the failure is silent precisely because the output looks like an ordinary similarity value — nothing in the number announces that it measures layout as much as overlap.

## 3.4 Repeats and k-mer redundancy

Set-based estimators count each distinct k-mer once, so repetitive sequence breaks the link between k-mer overlap and nucleotide divergence: a substitution in a high-copy element removes one distinct k-mer no matter how many copies it appears in, and the repeat content of a genome, not just its divergence, shapes the estimate. Three lines of work attack this at different layers. Skmer corrects Jaccard-derived distances for the low coverage and sequencing error of genome skims, but models repeat structure only implicitly [83]. ReSkmer derives the expected k-mer intersection between two repetitive genomes directly, jointly modeling repeat spectra, coverage, and error, and recovers accurate population-genomic distances even for highly repetitive inputs [84]. Most recently, Wu and Medvedev (preprint) define repeat-robust mutation-rate estimators built on "novel" k-mers — those arising from new mutations — and prove that the estimators combine with FracMinHash sampling without introducing systematic bias [85]. The corrections target different repeat regimes, and their comparative accuracy on real repeat-rich genomes has not been measured; that comparison is a second benchmark target carried into §5.

## 3.5 Divergence versus missing sequence

Containment answers "how much of set  $A$  appears in set  $B$ ," and nothing more. The Mash Screen authors state the caveat themselves: the identity formula "does not attempt to separate divergence from missing sequence" [16]. The ambiguity is real in both directions. A genome present in a sample but substantially diverged, and a genome absent except for one shared mobile element, can yield the same containment value; so can a genuinely diverged strain and a complete match recovered at 40% because the assembly is incomplete or the sequencing shallow. Section 2 separated the containment index — a set-overlap fraction — from the containment score, a derived per-nucleotide identity; the boundary here is semantic rather than algebraic. Low containment is not evidence of divergence, and high containment is not evidence of presence, until a second signal rules the alternative out. Screening results therefore need coverage breadth or aligned fraction alongside containment before any biological interpretation, a point skani makes constructively by estimating identity only over the regions it can chain [81].

## 3.6 ANI conversion model error

Even a perfect set-similarity estimate can yield a wrong ANI, because the conversion itself carries assumptions. The Mash distance of equation (2) rests on a Poisson model — random, independent point mutations acting on unique k-mers, with a fraction  $e^{−kd}$  of shared k-mers surviving divergence  $d$  [7]. Each clause of that model has a genomic counterexample. Indels violate independence in bulk: a single deletion destroys up to  $k$  overlapping k-mers at once, so indel-rich divergence reads as far more substitution than occurred. Compositional heterogeneity violates uniformity, since mutation is not equally likely across a genome whose k-mer distribution is skewed. Repeats violate uniqueness, compounding §3.4. Converted ANI is thus a model extrapolation layered on top of sampling error — accurate near the assumptions and systematically biased away from them, in directions that depend on the genome pair rather than the sketch. The discipline this demands is cheap: raw sketch estimates (Jaccard or containment) and converted ANI values are different objects with different error sources, and both should be reported — the raw estimate with  $k$  and scale factor, the converted value labeled with the model that produced it [7,18]. Reporting only the converted number hides the one quantity another researcher could recompute correctly.

## 3.7 Non-uniform sampling: when the predicate is a motif

Threshold and motif predicates share context-freedom but differ in what they assume about k-mer space. A hash threshold samples uniformly at random: every k-mer in every genome has the same retention probability, and that uniformity underwrites the estimators of §2.3 and the safe-scale-factor bound of equation (8). A motif samples deterministically and non-uniformly: retention is a property of sequence composition, density varies with GC content and k-mer usage, and the retained fraction cannot be set by the user below the enzyme's own. Three consequences follow. First, tag-poor genomes: a genome with few recognition sites yields a small sketch whatever its size, so the sample-size term of equation (8) is governed by the tag count rather than by a chosen scale factor. The protocol therefore records tags per genome for each enzyme and reports a no-estimate rate alongside accuracy (§5.2); a motif sketch that cannot estimate says so, which is the benign end of the failure spectrum of this section. Second, the ANI conversion must use containment, not Jaccard. Under the substitution model of §2.2 a tag survives divergence  $d$  with probability  $(1−d)^{k}$  regardless of how it was selected, so  $C=|A∩B|/|A|$  and  $d=1−C^{1/k}$  are unbiased; the Jaccard route of equation (2) is not, because a substitution inside the recognition site removes the tag from the mutated genome's set instead of replacing it with a freshly sampled k-mer, so the union shrinks and Jaccard is inflated relative to the hash case. On simulated genome pairs at 1% divergence the containment estimator recovers 0.990 ± 0.001 across replicate seeds from about 1,450 BcgI tags (§5.2); the Jaccard route is scored beside it in the benchmark so that the size of its bias is measured rather than asserted. Third, the physical and the computational sketch can disagree. Incomplete digestion, PCR duplicates, and sequencing error act on a wet-lab tag set but not on the in-silico one, so the agreement between the two is itself a measurable quantity — the only case in this review in which a sketch can be validated against a laboratory measurement rather than against another algorithm (§8, Q9).

# 4. Four tasks, four estimators

Section 2 sorted sketches by what they estimate; this section turns that ordering into a decision. The common genomic sketch tasks reduce to four estimation problems, and each requires a specific quantity. A family is eligible for a task only when its native estimand — not its label — matches that quantity. Where the match is exact we mark the family eligible; where the mechanism can serve only inside a larger pipeline we mark it partial; where the estimand itself is wrong we exclude it and give the reason. The question here is prior to accuracy: §3 asks how eligible estimators break, while this section asks which quantity the task demands at all.

## 4.1 Task definitions and the estimand each requires

### 4.1.1 T1: genome–genome ANI

The task compares two assemblies of comparable size, symmetrically, and reports average nucleotide identity (ANI) together with the aligned fraction (AF) — the share of each genome that participates in the comparison. The pure-sketch reference is Mash, which converts a bottom- $k$  Jaccard estimate into a distance through equation (2) [7]. That route returns no AF: a shared-k-mer fraction over the union conflates point divergence with regions that simply do not align, and for genomes shaped by horizontal transfer or rearrangement the two components must be separated. The accurate tools are therefore hybrids rather than pure sketches. The MashMap line used winnowed MinHash to estimate local Jaccard along mapped intervals [64,65], and MashMap3 replaced minimizers with minmers to make those local estimates unbiased [66]; FastANI segments the query, maps the segments MashMap-style, and averages identity over orthologous regions alone [82]; skani sparsely chains FracMinHash-sampled markers and estimates ANI strictly on shared regions, reporting AF alongside, which is what makes it robust on fragmented metagenome-assembled genomes [81]. The task's estimand is thus a symmetric distance computed over the alignable fraction. Every pure sketch supplies the first term and none the second, so pure sketches enter Table 3 as partial — the chaining or mapping layer supplies the rest.

Motif-defined tags supply the same two pieces by a different route. Because a tag's position is fixed by its recognition site, tag matches chain without a seeding step, and Syn2bANI reports ANI on the chained fraction together with AF and structural-variant calls [119]; the tool-free mechanism row — tags with a containment estimator — enters as partial, like the other pure sketches, since containment supplies the identity term and nothing supplies AF.

### 4.1.2 T2: detecting a genome in a metagenome

The task asks whether genome  $A$  is present in read set  $B$  with  $|A|≪|B|$ . The required quantity is the asymmetric containment  $C(A,B)$  of equation (1), followed by a statistical decision. Koslicki and Zabeti derived exact error formulas for a containment MinHash estimator aimed at differently sized sets [15]; Mash Screen engineered the asymmetry directly, sketching only the reference and streaming the mixture past it [16]. Detection, though, is a decision under uncertainty rather than a point estimate. sylph debiases containment ANI at low coverage with a zero-inflated Poisson model [25], and YACHT abandons the point estimate altogether for a hypothesis test that returns presence or absence with calibrated confidence [26]. Fixed-size bottom- $k$  sketches are structurally disfavoured here: their error is priced in units of union size (§2.2), so a 5 Mbp genome inside a 50 Gbp metagenome contributes almost no samples to the merged sketch. Eligibility for T2 therefore tracks the size-class column of Table 1 as much as the estimand column.

Tag-based profilers make the same decision on a species-specific subset of the reference tag set — tags unique to one species in the database — scored by breadth times depth (the G-score) and passed through a learned false-positive filter [117,118]. The sample-size regime is the one that governs scaled sketches (§3.1): at 0.01× coverage a 2,400-tag genome contributes about 24 tag observations, so the T2 coverage sweep tests the hash and motif predicates on the same curve.

### 4.1.3 T3: metagenome–metagenome comparison

The comparison is symmetric again, but the objects are distributions: relative abundances are the signal, and the required estimands are weighted — the weighted Jaccard  $J_{w}$ , the probability Jaccard  $J_{P}$ , and Bray–Curtis dissimilarity, whose exact multiset definitions genomics inherited from k-mer counting [86]. The three are not interchangeable.  $J_{w}$  and Bray–Curtis are monotone reparameterizations of each other on count vectors —  $J_{w}=S/(2−S)$  for Bray–Curtis similarity  $S$  — while  $J_{P}$  is a genuinely different quantity, scale-invariant where  $J_{w}$  is not (§2.4). Sketch implementations run on two mechanisms. HistoSketch and its microbiome application HULK brought consistent weighted sampling to k-mer spectra [39,40], with the caveat that the collision probability was later proved to equal  $J_{P}$  rather than the  $J_{w}$  the framing claims [38]. Dashing 2 instead pairs register sketches with BagMinHash and ProbMinHash to sketch multiplicities directly [45]; SimkaMin sketches the multiset distances themselves [87]; and kWIP's weighted inner product attacks the same distributional question from a count-sketch base [88]. Any benchmark that reports "weighted Jaccard" without saying which one is reporting an ambiguous number.

Tag counts are a weighted sketch by construction — each tag occurrence is a read — so tag-count vectors support Bray–Curtis directly, and species profiles built from them support the profile-level Bray–Curtis that most microbiome studies compute; both enter T3.

### 4.1.4 T4: read-level classification and filtering

The task makes a per-read membership decision: does this read derive from an indexed reference set? No set-similarity quantity is involved. What is required is a per-read hit guarantee — any shared substring of length  $w+k−1$  must yield a shared seed — which is exactly the window guarantee of minimizers [58,59]; syncmers keep a positional guarantee while making selection context-free [61]. deacon (preprint) shows the task in its purest current form, querying a human-pangenome minimizer index to discard host reads at over 250 Mbp/s on a laptop [75]. Similarity sketches are an estimator mismatch for T4. The task asks a membership question, and a bottom- $k$ , threshold, or register sketch answers a different one — how similar are these two sets — that nobody posed; at read length the answer is usually not even computable, because the read's sketch is empty (§4.3). The mismatch runs in both directions: the window schemes that dominate T4 are themselves disqualified as similarity estimators, their Jaccard estimator being biased and inconsistent [60]. Both directions are enforced by theorem, not by benchmark taste.

Motif tags are excluded for the same reason as threshold sketches: most reads contain no recognition site — about 6% of 151 bp reads contain a BcgI core — so there is no per-read guarantee. A 2b-RAD library handles host DNA by not sequencing most of it, which is a library-design answer to the host problem rather than a filtering one.

## 4.2 Table 3: the task × family eligibility matrix

Table 3 assembles the four task definitions against the seven families of §2, with one line of justification for every exclusion or partial match.

Table 3. Task × family eligibility matrix. ✓ = native estimand matches the task; partial = usable only with caveats or as a pipeline component; ✗ = estimand mismatch. Letters refer to the reason lines below. F2 is split into its hash-threshold (F2h) and motif (F2m) predicates because their cost models differ.

| Task (required estimand) | F1 bottom- $k$ | F2h threshold | F2m motif | F3 weighted | F4 register | F5 order | F6 sampling | F7 vector |
|---|---|---|---|---|---|---|---|---|
| T1 genome–genome ANI (symmetric distance + AF) | partialᵃ | partialᵃ | partialˢ | ✗ᵇ | partialᵃ | ✗ᶜ | partialᵈ | partialᵉ |
| T2 genome detection in a metagenome (asymmetric containment + test) | ✗ᶠ | ✓ | ✓ | ✗ᵍ | partialʰ | ✗ᶜ | ✗ⁱ | ✗ʲ |
| T3 metagenome–metagenome comparison (weighted distance) | ✗ᵏ | partialˡ | ✓ | ✓ | ✓ | ✗ᶜ | ✗ⁱ | ✗ᵐ |
| T4 read classification/filtering (membership + window guarantee) | ✗ⁿ | ✗ᵒ | ✗ᵗ | ✗ᵖ | ✗ᵠ | ✗ʳ | ✓ | ✗ᵠ |

a. Symmetric Jaccard/distance estimand matches, but no aligned fraction; AF requires a chaining or mapping layer (MashMap line, RECOMB 2017 [64]; Bioinformatics 2018 [65]; Bioinformatics 2023 [66]; FastANI, Nat Commun 2018 [82]; skani, Nat Methods 2023 [81]).

b. Multiplicity estimands answer distributional questions T1 does not pose; weights add no signal to genome–genome comparison.

c. Edit-distance proxy defined per sequence; sketches cannot merge across genomes and no ANI estimand exists (Bioinformatics 2019 [50]).

d. Window-minimum selection is a biased, inconsistent Jaccard estimator [60]; eligible only as the seed layer of sketch-chaining hybrids [64,81].

e. Inner-product→ANI chain validated only above ~85% ANI and reports no AF (Bioinformatics 2024 [77]).

f. Error priced in units of union size; the small genome's hashes vanish in the metagenome union [15,16].

g. Abundance weighting conflates detection with coverage; no presence/absence test formulation exists.

h. Containment via maximum-likelihood estimation exists, but asymmetric-set error is uncharacterized (Genome Res 2023 [45]).

i. Minimizer-based Jaccard and containment estimation is biased and inconsistent [60]; no multiplicity channel.

j. Containment untested; the inner-product estimand targets symmetric comparison [77].

k. Unweighted set overlap discards the abundance signal that defines the task.

l. Abundance-tracking variants carry counts, but the native estimands remain unweighted; weighted distances are derived, not estimated (JOSS 2024 [21]).

m. No validated weighted estimand [77].

n. Context-free global sampling gives no per-read hit guarantee; a fixed- $k$  sketch of a ~150-k-mer read is nearly empty.

o. Sketch size  $Θ(n/s)$  leaves most reads with zero retained hashes at standard scale factors (FracMinHash preprint [17]).

p. Single-read weight vectors are too small for consistent sampling; no membership semantics.

q. The representation answers set-level queries only; there is no per-k-mer membership query.

r. Order sketches cost ~ $ℓ×$  per entry and answer edit similarity, not membership [50].

s. Containment ANI is unbiased for motif tags, but AF requires chaining on tag anchors (Syn2bANI [119]); density is fixed by the enzyme and tag-poor genomes return no estimate (§3.7).

t. Most reads carry no recognition site, so there is no per-read membership query (§4.1.4); the family's answer to host DNA is at the library, not the read.

Read by column, no family is eligible everywhere, and the ✓ cells are nearly diagonal: T2 admits the two context-free predicates, T3 the weighted mechanisms (F3, F4, and the motif tags, whose counts are weights by construction), T4 only the index samplers, with the single exception of T1, where the best answers are hybrids built on F6 seeds or on tag anchors. That exception is the table's main lesson: ANI on the alignable fraction is not a set-overlap quantity, so no pure selection rule estimates it. Read by row, the partial cells cluster in T1, and F1 exemplifies the pattern: its estimand is right for T1 but its mechanism carries no positional information, and its fixed budget disqualifies it wherever set sizes differ. The F6 column compresses the whole argument into one diagonal pair — ✗ as an estimator, ✓ as an index. Section 5 benchmarks only the ✓ and partial entries of this matrix; the ✗ cells are excluded by construction, not by measurement.

## 4.3 The cost of mismatch: three recurring errors

Three mismatch patterns recur in the literature and in practice. The first is applying a T1 instrument to T2: a Jaccard-optimized sketch used to detect a genome in a metagenome. The Mash Screen authors state the failure directly — Mash "cannot reliably estimate the containment of a genome within a metagenome" [16] — and the theory explains it: the relative error of MinHash-based estimation grows sharply when  $|A|≪|B|$  [15]. The failure is silent, because the tool returns a small number indistinguishable from true absence.

The second is reporting a converted ANI without the raw estimate beneath it. Mash Screen itself separates the assumption-light containment index from the model-laden containment score  $c_{k}^{ 1/k}$  [16], yet pipelines that publish only the converted number discard exactly the information needed to check it, and every conversion assumption — k-mer survival, Poisson mutation (§3.6) — is carried invisibly. sylph's low-coverage corrections show how large the modeling adjustment can be before an ANI figure is trustworthy [25]. Raw estimate and converted quantity should travel together.

The third is read-level similarity search with threshold sketches. At a scale factor of  $s=1000$ , a 150 bp read contributes about 120 k-mers at  $k=31$  and an expected 0.12 retained hashes, so most reads have empty sketches and the search returns zero hits — an empty result indistinguishable from "no match" (FracMinHash preprint [17]). The safe-scale-factor analysis quantifies how many elements a threshold sketch must retain before its estimates carry any guarantee at all [19]. The correct instrument for the question is a membership index with a window guarantee [58,75]. All three errors share one structure: the tool computed a serviceable estimate of the wrong quantity, and no improvement to the estimator repairs that.

# 5. A benchmarking protocol

Section 4 closed with an eligibility matrix (Table 3): only families whose native estimand matches a task — the ✓ and partial cells — are meaningful contestants for it. This section turns that screen into an executable benchmark for four tasks (T1–T4), each with fixed data, two-tier ground truth, a parameterization rule, and a statistical plan. The protocol is implemented as a public SLURM pipeline [126]; every number below is a value fixed in its configuration file, the design was committed before any result existed, and results await execution on HPC resources and appear as marked placeholders (§5.6). Each design answers an open measurement question from §3 — confidence interval (CI) coverage (§3.2), scale-factor safety (§3.1), the motif-versus-hash comparison at matched density (§3.7), and sequencing-depth stability, posed as its own protocol question via the §5.4 depth sweep (§8, Q4).

## 5.1 Design principles

Seven rules govern all four tasks.

Rule 1 — eligibility by estimand; defaults and a size sweep; a bytes axis. A tool enters a task only from a ✓ or partial cell of Table 3; ✗ cells are excluded by construction, not by measurement (§4.2). Every entrant is run under its documented default parameters and along a sketch-size sweep (Mash  $s∈{1000, 4000, 16000, 64000}$ ; FracMinHash  $s∈{10000, 2000, 500, 100}$ ; Dashing  $log_{2}m∈{10,12,14,16}$ ; sylph  $c∈{1000, 200, 50}$ ; HULK  $s∈{512, 2048}$ ; deacon thresholds  $a∈{1,2,3}$ ,  $r∈{0, 0.01, 0.05}$ ; the motif rows over enzymes BcgI, BsaXI, CspCI and a BcgI + AlfI panel). All accuracy-versus-cost results are plotted against bytes per sketch rather than sketch count or sampling parameter, because the families scale differently: bottom- $k$  sketches store a fixed number of hashes [7], threshold sketches grow linearly with distinct-element count [19], register sketches trade register count against bits per register [44], and motif sketches grow with the number of recognition sites. Bytes is the only currency in which these are commensurable. Size-matched comparisons are additionally reported at a fixed budget per task — the grid point nearest 32 KB per genome in T1 and 25 MB per sample in T3 — so that "more accurate" is never read across a difference in size.

Rule 2 — one  $k$  per task. T1 uses  $k=21$ , the value at which the Mash distance model was calibrated [7]; T2–T4 use  $k=31$ , the value shared by the containment profilers [25] and the filtering indexes [75] that define those tasks. Motif tags have their own  $k$  — the double-stranded core length of the enzyme (32 for BcgI, 27 for BsaXI, 33 for CspCI) — which is reported with every estimate. A  $k$  ablation is run for Mash in T1 only ( $k∈{16, 21, 31}$ ), where the Jaccard-to-ANI conversion is most  $k$ -sensitive (§3.6).

Rule 3 — two-tier ground truth. Each task pairs a simulated tier, where the true estimand is known exactly by construction, with a real-data tier anchored to an alignment-based or exact-counting consensus. Where reference methods disagree, as alignment-based ANI estimates do near species boundaries [82], the disagreement is reported as a spread, not resolved into one number. A benchmark that validates sketch methods against other sketch methods cannot detect a bias shared by the whole family, so the anchor tier must be estimand-faithful even where expensive: the median of FastANI [82], skani [81], and ANIm [110] with its spread for T1; designed abundance vectors plus exact  $k$ -mer multiset counts [86] for T3; mock communities of known composition [108] and the CAMI II truth sets [105] for T2; per-read labels from simulation and single-organism real runs for T4.

Rule 4 — timing protocol. Wall-clock time, user and system time, peak resident memory, and on-disk size are recorded separately for index construction and query, with input bases recorded for throughput. Inputs and indexes are copied to node-local storage ($TMPDIR) before the clock starts so that network filesystem latency does not contaminate timings. Every default-parameter measurement is repeated three times at 1 and 16 threads and the median reported; the sweep grids are timed once at 16 threads. The node and CPU model are recorded with every row; earlier benchmarks pooled stages or used a single thread count [94].

Rule 5 — version locking and data archiving. Every tool is pinned to an exact version — a conda environment for packaged tools, a commit hash for source-built ones — and a version table is written before any run. Sketches, configurations, logs, and result tables are archived. The need is not hypothetical: datasketch v2.0.0 changed its default MinHash permutation scheme to affine32 to fix a similarity over-estimation bias on large sets, rendering hash values incompatible with earlier releases [101,114]. A comparison that pools results across such a change measures the library version, not the algorithm.

Rule 6 — the statistical unit is the genome, not the pair. Fifty genomes yield 1,225 pairwise distances, but those pairs share genomes and are not independent samples. All tool-versus-tool differences are therefore tested with a bootstrap clustered at the genome level, resampling genomes with replacement and recomputing every induced pair under both tools, and confirmed with a mixed model carrying crossed random effects (random intercepts for genome  $i$  and genome  $j$ ). Resampling a genome twice induces a self-pair; self-pairs are excluded from each bootstrap replicate. For a task estimate  $\hat{θ}_{i}$  with reference value  $θ_{i}$  over  $n$  evaluation units (pairs, samples, or reads), the protocol reports

$bias=(1)/(n)∑_{i=1}^{n} (\hat{θ}_{i}−θ_{i}),  RMSE=sqrt((1)/(n)∑_{i=1}^{n} (\hat{θ}_{i}−θ_{i})^{2}),$ 	(10)

and, where a tool emits a nominal 95% CI  $I_{i}$ ,

$coverage=(1)/(n)∑_{i=1}^{n} 1{θ_{i}∈I_{i}}.$ 	(11)

A gap between tools is called a difference only when the clustered-bootstrap interval for the paired RMSE contrast excludes zero; untested numerical gaps are described as numerical gaps, never as trends. Tools that fail to compile or run under the pinned environment are excluded and reported as excluded, not silently dropped — the deprecated CMash repository, marked by its authors as superseded by sourmash and YACHT [23,113], is the standing case.

Rule 7 — a sequencing-cost axis. Bytes per sketch measure what a sketch costs to store; they do not measure what it costs to obtain. For every entrant the protocol records the input bases from which the sketch was computed, and for the motif rows additionally the fraction of those bases that lie inside tags — the bases a restriction library would have sequenced to produce the same sketch. Figure 2 carries a second cost panel on this axis. Hash-based sketches sit at 100% of input bases by construction; the panel therefore measures the one property of §2.3.4 that byte accounting cannot see.

## 5.2 T1: genome–genome ANI

The simulated tier mutates 50 GTDB R202 species representatives [104], stratified by GC content and genome size, at nine divergence rates (0.1, 0.5, 1, 2, 5, 10, 15, 20, and 25%) in two parallel series — substitution-only, and substitution plus indels at one indel per ten substitutions with geometric lengths — recording the realized substitution and indel counts so that the true identity is exact rather than nominal, in both its substitution-only and gapped forms. The two series separate the failure modes of §3.5 and §3.6: indels remove sequence, substitutions only diverge it, and a sketch that conflates the two will drift between series where an AF-aware method will not. The real tier consists of 2,000 GTDB R202 pairs, 500 at each lowest shared rank (species, genus, family, order), so that the 75–100% ANI range is populated without selecting on any sketch tool's output; the anchor is the median of FastANI, skani, and ANIm per Rule 3, with the spread between them reported. A 5,000-genome all-versus-all run measures throughput at database scale. The motif rows are preceded by an enzyme-density scan over the same 5,000 genomes — tags per genome and per megabase for six enzymes, the Spearman correlation of density with GC content, the fraction of tag-poor genomes, and the FracMinHash scaled-equivalent  $(L−k+1)/n_{tags}$  — which is the data behind the claim that BcgI sits near  $s≈2000$  and the basis for the size-matched comparison of Rule 1.

Entrants are the ✓ and partial cells of Table 3: Mash [7], BinDash [11], Dashing 2 [45], sourmash [21], sylph in genome mode [25], HyperGen [77] and MaxGeomHash (preprint, RECOMB 2026) [28] where installable, the sketch-chaining hybrid skani [81], and the motif rows — in-silico tags with both the containment and the Jaccard estimator reported so that the bias of §3.7 is measured, and Syn2bANI with its default four-enzyme panel and with BcgI alone [119]. FastANI [82] and ANIm serve as reference baselines. Metrics are bias, RMSE and mean absolute error per ANI bin (75–80, 80–85, 85–90, 90–95, 95–99, 99–99.5, 99.5–100%; equation 10), CI coverage where a tool emits intervals (equation 11, targeting §3.2), the no-estimate rate, the Spearman correlation with truth, the detection floor — the lowest consensus ANI at which at least 95% of true pairs are called — and RMSE versus bytes under defaults and along the sweep. Raw-estimate scoring — Jaccard or containment RMSE against exact simulated values, kept separate from converted ANI — isolates conversion error from sketch error (§3.6).

## 5.3 T2: genome detection in a metagenome

The simulated tier spikes 20 target genomes into a 30-genome log-normal background at seven coverages — 0.01, 0.05, 0.1, 0.5, 1, 5, and 20× — rotated across targets so that every coverage is represented in every community, with reads generated by InSilicoSeq [107] under its NovaSeq model at 2 × 151 bp, at 10 million and 50 million read pairs and three seeds each. This spans the regime where scale-factor safety fails (§3.1): a threshold sketch of a 0.01× genome retains only a handful of hashes, and a motif sketch a few dozen tags. A strain panel spikes siblings of the targets at 0.5, 1, 2, 3 and 5% divergence while the indexed reference itself is absent, measuring strain-level specificity; the 5,000-genome GTDB pool serves as the reference database and yields an empirical false discovery rate (FDR) on the genomes known to be absent. The real tier uses the ZymoBIOMICS even and log-distributed community standards (D6300, D6310) [108] and the CAMI II marine and strain-madness datasets [105]. Entrants are the containment family: Mash Screen [16], sourmash gather (preprint) [17], sylph [25], YACHT [26], the in-silico tag row (containment of each reference tag set with hit-count and fraction thresholds), and Fast2bRAD-M with its G-score decision [120]; MaxGeomHash [28] and Dashing 2 [45] enter where installable to test the hypotheses behind their partial cells. Kraken 2 [73] classification is reported as a non-sketch context comparator, and minimap2 read mapping supplies the coverage anchor on the 10-million-pair sets. A scale-factor sweep on Rule 1's grid records observed containment error at each scale factor, for plotting against the predicted safe-scale-factor boundary of equation (8) (§3.1).

Two estimands are scored separately and never averaged: species-level sensitivity (calling a genome present when it is) and strain-level specificity (not calling the indexed reference when only its sibling is present). Collapsing them rewards tools that are liberal at both levels. Metrics are precision, recall, F1 and the area under the precision–recall curve over the tool's score, FDR on absent genomes, recall per coverage bin, and — for tools that report them — the mean absolute error of ANI and of log10 coverage. For the decision tools, the scored quantity is empirical FDR at each nominal confidence level — the fraction of calls at that level that are absent genomes — answering whether sketch-based uncertainty statements mean what they claim (§3.2).

## 5.4 T3: metagenome–metagenome comparison

The simulated tier is a  $3×20$  design: three groups of 20 samples drawn from a 200-genome pool, each group built on a 60-genome core with group-specific log-normal mean abundances plus 20 genomes shared by all groups, with per-sample log-normal noise ( $σ=0.6$ ) around the group vector, at 5 million read pairs per sample. Ground truth is two-tiered within the task: Bray–Curtis, Jaccard and probability-Jaccard distances computed exactly from the designed abundance vectors (what a profiler tries to recover) and exact  $k$ -mer multiset Bray–Curtis computed by Simka from the reads [86] (what a  $k$ -mer sketch literally estimates). Depth stability is measured by re-simulating a subset of samples at 1, 5, and 20 million pairs, probing whether estimates move with sequencing effort (§3.1). The real tier uses Human Microbiome Project samples across body sites [106], with body-site labels as the external anchor.

Every entrant's output is cross-scored against both weighted ground truths,  $J_{w}$  and  $J_{P}$  (equations 6–7). The reason is concrete: the collision probability of HistoSketch, which HULK implements, equals  $J_{P}$ , not the weighted Jaccard its framing claims [38,40]. Cross-scoring turns that correction into a measurement: a tool estimating  $J_{w}$  scores well against the  $J_{w}$  column and poorly against  $J_{P}$ , and vice versa. Entrants are sourmash compare with and without abundance tracking [21], HULK [40], Dashing 2 in SetSketch and ProbMinHash modes [45], Mash on reads with and without a minimum-count filter [7], sylph profiles converted to Bray–Curtis [25], the in-silico tag row (tag-count Bray–Curtis and tag-presence Jaccard), and Fast2bRAD-M species profiles converted to Bray–Curtis [120], with exact Simka [86] as the reference baseline. Metrics are the Mantel (Spearman) correlation against each truth, the adjusted Rand index of Ward clustering against the design groups, a Procrustes distance between ordinations, the depth-stability displacement between depth replicas, and, on the real tier, the adjusted Rand index against body site.

## 5.5 T4: read-level classification and filtering

The simulated tier mixes reads from a human genome (GRCh38 or T2T-CHM13) at 50% with a 30-genome microbial background under InSilicoSeq's NovaSeq model (three seeds) and its HiSeq and MiSeq models, then derives sweeps from the NovaSeq base: post-hoc substitution error at 0.1, 0.5, 1 and 2% per base, and truncation to 50, 75 and 100 bp, so that the window guarantee qualifying this family (§4.1.4) is stressed by error and by read length separately. An inter-individual divergence panel replaces the reference with the HG002 maternal and paternal assemblies [111] as the source of host reads while the index stays reference-built, measuring how sensitivity erodes as query and index diverge. The real tier uses read sets whose provenance fixes their label: CHM13 reads (all host), the FDA-ARGOS and RSV read sets released with deacon (all non-host) [75,112], and an HPRC individual's short reads (host). The central manipulation holds tool, scheme, and parameters constant while only index content changes — deacon under a T2T-only index and under panhuman-1 (Zenodo record 17288185 [112]) — and, separately, holds tool, index content, and parameters constant while only the selection scheme changes: deacon with  $w=15$  minimizers against a fork of deacon selecting closed syncmers with  $s=16$ , which yields the same index size (§7.3). Because one factor varies at a time, the index effect and the scheme effect are separated — a distinction that is confounded whenever tool, scheme, and index content vary together. Kraken 2 [73] with T2T-only and mixed databases, Bowtie2, hostile, and minimap2 [72] serve as non-sketch reference points; read-level sylph and sourmash are run as cautionary rows, their ✗ cell of Table 3 measured rather than assumed. Metrics are per-read sensitivity, specificity, F1 and Matthews correlation, residual host reads per million input reads, microbial reads lost per million, throughput, index bytes, and build time under Rule 4.

## 5.6 Results placeholders

The benchmark has not yet run on HPC resources; the three result objects are specified below so outputs drop in without restructuring. Figure 2 will carry the accuracy–cost trade-off for all four tasks.

[PLACEHOLDER: pending HPC results — Figure 2: four-task RMSE-versus-bytes panels plus a bases-sequenced panel; axes, legend and panel layout specified below]

Figure 2 is a four-panel row, one panel per task, with a fifth panel for the sequencing-cost axis of Rule 7. The x-axis of every task panel is bytes per sketch (log scale, 10²–10⁷); the y-axis is RMSE against the tier-appropriate ground truth (T1: ANI RMSE in percent; T2: sensitivity at fixed empirical FDR; T3: Bray–Curtis RMSE; T4: false-negative rate). Default-parameter runs are open markers, sweep runs filled, one color per family, with the two F2 predicates distinguished by marker shape. The T2 panel overlays the predicted safe-scale-factor boundary of equation (8), against which the §5.3 sweep is read. The fifth panel plots T1 and T2 accuracy against bases sequenced to obtain the sketch, on which every hash sketch sits at 100% of input and the motif rows at the tag fraction. The figure answers four questions: which family sits on the Pareto frontier for each task; whether the frontier shifts between default and matched parameters, indicating miscalibrated published defaults; whether the byte-for-byte ordering of families is stable across tasks or reverses, quantifying how far "which sketch is best" is task-dependent; and what accuracy costs when the sketch is sequenced rather than computed.

Figure 2: RMSE versus bytes per sketch for tasks T1–T4, and versus bases sequenced for T1–T2. Open markers, default parameters; filled markers, sweep parameters. (Placeholder — HPC results pending.)

Figure 3 addresses uncertainty calibration rather than point accuracy.

[PLACEHOLDER: pending HPC results — Figure 3: CI coverage panels; axes, legend and panel layout specified below]

Figure 3 plots empirical CI coverage (equation 11) for every tool that emits uncertainty statements in T1 and T2, with nominal 95% coverage as a horizontal reference line. The x-axis is the truth-stratifying variable per task (consensus ANI bin for T1, spike-in coverage for T2), so that undercoverage can be localized: the open question from §3.2 is not whether intervals are wrong on average but where they fail, and theory predicts failure exactly at low coverage and high divergence [18]. A reader can read off, for each tool, where its reported intervals can be trusted and where they should be ignored. Tools that emit no intervals appear in a side panel of point-estimate standard errors, so the absence of uncertainty quantification is itself visible as a gap.

Figure 3: Empirical coverage of nominal 95% confidence intervals against the T1/T2 ground truths. (Placeholder — HPC results pending.)

Table 5 will collect the headline number per task and tool; its structure is fixed below.

[PLACEHOLDER: pending HPC results — Table 5: task × tool main results; column headers specified below, rows grouped by task]

Table 5. Main benchmark results by task and tool. One row per task × tool × parameter set. RMSE/bias from equation (10); coverage from equation (11); times under Rule 4 (median of three repeats, 16 threads). For T1 rows the raw-estimate column carries Jaccard/containment RMSE against exact simulated values beside converted-ANI RMSE, separating sketch error from conversion error (§3.6). T3 rows report RMSE against both  $J_{w}$  and  $J_{P}$  ground truths (§5.4). The bases-sequenced column is 100% for every hash sketch and the tag fraction for motif rows (Rule 7). All cells pending HPC execution.

| Task | Tool (pinned version) | Family | Parameters (default / sweep) | Bytes per sketch | Bases sequenced (% of input) | Primary accuracy (RMSE or sensitivity/specificity) | Raw-estimate RMSE (T1 only) | CI coverage (nominal 95%) | Build time (s) | Query time (s) | Notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| T1 | [pending] |  |  |  |  |  |  |  |  |  |  |
| T2 | [pending] |  |  |  |  |  |  |  |  |  |  |
| T3 | [pending] |  |  |  |  |  |  |  |  |  |  |
| T4 | [pending] |  |  |  |  |  |  |  |  |  |  |

Table 5 is the lookup table a practitioner reads first: it compresses each entrant to the byte cost it charges, the error it returns at that cost, whether its uncertainty statements are calibrated, and what its build-versus-query profile implies for deployment. Because every row carries both parameter sets, the table shows how much of each tool's published performance is tuning rather than the sketch itself. Rows for excluded tools are retained with the exclusion reason in Notes, so the table documents the full eligibility set, not only survivors.

## 5.7 Reproducibility

The pipeline [126] is a set of SLURM stages chained by job dependencies — fetch, simulate, manifest, build, query, metrics, figures — driven by one configuration file that holds every path, grid, and seed. A uniform wrapper per tool implements two calls, build and query, with a fixed argument list and one of four standard output tables (pair estimates, detections, sample distances, per-read calls), so that numeric differences cannot be traced to wrapper idiosyncrasies; a single genome-naming rule is applied to every output before scoring. Every dataset, index, and run directory carries a completion marker and is skipped when present, so a failed array task is re-run by resubmitting the stage rather than by restarting the benchmark. Timing rows from all tasks share one schema. A tenfold-smaller "minimal" mode exercises every wrapper end-to-end on small inputs and is the first thing run on a new cluster, because the most common failure of a benchmark harness is an interface that drifted between the version the wrapper was written against and the version installed. Sketches, configurations, logs, and result tables are archived on Zenodo with a DOI per release. Executing the archived pipeline populates the §5.6 placeholders; re-running it on new tools extends the benchmark without redesign. The protocol, not any result table, is the deliverable.

## 5.8 Planned extensions

Five measurements are specified but not yet in the pipeline, and are listed so that their absence is not mistaken for a design decision: long reads (10–15 kb) in T4, where the window guarantee behaves differently; a repeat-fraction stratification of the T1 genomes, so that the repeat-robust estimators of §3.4 — Skmer [83], ReSkmer [84], and the novelty-based estimators [85] — can be scored in the regime they target; a  $k$  ablation beyond Mash; the adjacent-family and non-sketch comparators GSearch [41], BinDash 2.0 [12], kssd [31], Kmer-db [89], Vclust [90], and ska2 [91] where installable; and the wet-lab concordance test of §8, Q9, which needs a mock community sequenced both by shotgun and by 2bRAD-M.

# 6. Cross-pollination: genomic sketching and LLM data deduplication

## 6.1 One ancestor, two lineages

Genomic sketching and the deduplication pipelines of large language model (LLM) training corpora descend from a single publication: Broder's 1997 definition of document resemblance and containment, written to find near-duplicate web pages for the AltaVista search engine [1]. The genomic lineage begins with Mash, which in 2016 recast minwise hashing as a k-mer sketch for genome distance estimation [7]. The LLM lineage starts earlier than is often assumed: GPT-3 already fuzzily deduplicated its training documents with Spark's MinHashLSH implementation using 10 hashes [102], and the practice became load-bearing infrastructure in the RefinedWeb (preprint) [103] and FineWeb (preprint) [99] corpora, where MinHash deduplication over Common Crawl snapshots is a first-class stage of the data pipeline. Each lineage rebuilt the vocabulary. Shingles became k-mers; resemblance became the Jaccard index and then average nucleotide identity; near-duplicate documents became near-identical genomes. The underlying estimator, in both worlds, is still Broder's.

## 6.2 The parallel parameter problem

The two fields face the same design choice — how many hashes, and how to split them into bands and rows for locality-sensitive hashing (LSH) — and resolve it the same way: by convention. On the LLM side, FineWeb deduplicates with word-level 5-gram MinHash, 112 hash functions split into 14 buckets of 8 rows, targeting documents at least 75% similar [99]; the SlimPajama configuration inherited by SEDD uses 128 hashes as  $b=16$  bands of  $r=8$  rows [100]; RefinedWeb (preprint) spent 9000 hashes as 450 buckets of 20 [103]. On the genomic side, the defaults are  $k=21$  with sketch size  $s=1000$  in Mash [7] and scaled sketching at  $s=1000$  in sourmash [20]. The banding choice controls the match probability through the standard LSH curve

$P(match∣J)=1−(1−J^{r})^{b},$

where  $J$  is the true Jaccard similarity (not to be confused with the sketch-size  $s$  of §2.1),  $r$  the rows per band, and  $b$  the number of bands. All three LLM configurations place the steepest point near  $J≈0.71$ – $0.74$ , but RefinedWeb's 80-fold larger sketch buys a far sharper transition, effectively trading memory for a crisp acceptance boundary, while FineWeb accepts a soft shoulder to keep signatures small [99,103]. Neither literature benchmarks its chosen operating point against measured estimator error for its own workload; the numbers persist because earlier pipelines used them. This shared reliance on folklore (the same parameter-selection rule, inherited rather than derived, on both sides; §3.1) is the strongest evidence that the two fields are solving one problem twice.

## 6.3 Mutual non-citation, with one engineering exception

Despite the common ancestor, the literatures do not talk to each other. In the corpora we examined, FineWeb's 80 references contain no genomic sketching work, SEDD's reference list likewise stops at Broder 1997, and the Mash-lineage papers — Mash [7], sourmash [20], GSearch [41] — cite no LLM deduplication literature; a direct search for "sourmash" against "FineWeb" or "datatrove" returns nothing. We cannot exclude isolated cross-references outside these corpora, but the pattern is consistent enough to treat as structural. The one verified bridge is code, not citation: GSearch [41] is built on its authors' Rust HNSW crate hnsw_rs, which its maintainers report is reused in LLM vector-database and semantic-search applications [115]. Engineering artifacts cross the boundary freely; ideas, apparently, do not.

## 6.4 What each side could take from the other

The transfer is asymmetric and mostly untried. Genomics has already imported from the broader sketching literature: GSearch combines ProbMinHash, Ertl's locality-sensitive scheme for the probability Jaccard similarity of weighted sets [38], and the register-based SetSketch [43] with graph indexing. The reverse direction is open. Weighted sketches such as ProbMinHash address multiplicity directly, which is exactly the structure of web text, where boilerplate and templated content repeat at scale; LLM deduplication currently handles this with set-based MinHash and post-hoc heuristics rather than weighted estimators. In the opposite direction, the LLM side has solved an engineering problem genomics still treats as exotic: SEDD (preprint) executes the full MinHash-LSH deduplication of a 1.2-trillion-token corpus on GPUs in about 2.9 hours, up to  $158×$  faster than a CPU datasketch baseline, using dense-bucket strategies and tiled all-pairs comparison kernels [100] — machinery directly applicable to billion-genome-scale sketch search. Finally, both worlds share a reproducibility failure mode: the widely used datasketch library changed the default MinHash permutation scheme to affine32 in a 2026 major release, silently invalidating sketches computed under earlier defaults [101,114]. Genomic sketch formats, which embed version and parameter metadata in the signature (§2.4.4, §5 Rule 5), offer a ready-made answer to that problem.

# 7. Building your own sketch tool in the era of coding agents

## 7.1 Why the primitive catalogue matters now

Table 2 lists some thirty tools, and most are recombinations of about ten primitives: a rolling hash (ntHash, xxhash, or a fixed random permutation), a selection predicate (bottom- $k$ , threshold, syncmer, minimizer, motif), a container (sorted vector, register array, hash set, binary fuse or Bloom filter), an estimator (Jaccard, containment, maximum likelihood on registers), a conversion model (equation 2 or its containment form, equation 9), and a decision rule (threshold, hypothesis test, learned filter). The primitives are well engineered and openly available — SIMD minimizer and syncmer selection at hundreds of megabases per second [70], compact membership filters, and the sourmash and Dashing libraries expose them through stable interfaces [21,45] — and coding agents have made their recombination cheap. In our own experience (§7.3, §7.4), an agent can produce a working, tested sketch tool or extension in hours when the specification is complete. What an agent does not supply is the specification: the estimand, the guarantee, and the test that would reveal a silent failure of the kind catalogued in §3. Those are exactly the three things this review's taxonomy makes explicit, and they are where a biologist's effort now belongs. The shift is from writing the sketch to specifying it.

## 7.2 A specification template

The six questions of Box 1 turn into a specification when each is answered before any code exists.

1. Estimand. Which quantity of equation (1), (6)–(7), or none of them (membership)? Name the column of Table 2 the tool will occupy. "Similarity" is not an answer; an agent asked for similarity will implement Jaccard.
2. Size class and asymmetry. Will the two inputs be comparable in size (T1, T3) or differ by orders of magnitude (T2)? This decides between fixed-size and threshold sampling before any accuracy question arises (§2.2.3, §3.1).
3. Guarantee. Unbiased? With a confidence interval? With a window guarantee? Under which mutation model? The guarantee names the theory the tool inherits and the assumptions §3.6 says must travel with every estimate.
4. Algebra. Will sketches be merged, subtracted, or intersected? If so, the selection rule must be context-free (Table 1), which excludes minimizers and fixes the choice within F2.
5. Cost budget. Bytes per sketch, seconds per query, and — if the sketch is to be measured rather than computed — bases sequenced (Rule 7).
6. Test oracle. A simulation in which the estimand is known exactly, and a real-data anchor computed by an independent method. The oracle is written before the tool, and the tool is not allowed to read it.

A request phrased in these terms is checkable: the agent's output either occupies the specified cell of Table 2 with the specified guarantee, verified by the specified oracle, or it does not. The benchmark harness of §5 is one such oracle, built to be reusable: a new tool enters it by implementing the two-call wrapper interface of §5.7 and is scored on the same truth as every existing entrant.

## 7.3 Worked example 1: adding a selection scheme to an existing tool

deacon depletes host reads by querying a minimizer index of human assemblies [75]. Minimizers are context-dependent (§2.7.1): a sequencing error inside a window can change which  $k$ -mer is selected, so a read can lose a match without losing the underlying  $k$ -mer. Closed syncmers are context-free and were reported to conserve more  $k$ -mers under mutation [61], which suggested a testable hypothesis: at matched index size, a syncmer index should retain more host reads at realistic error rates. The human decisions were four. The estimand is per-read membership (T4), so the evaluation metric is host reads retained and non-host reads falsely retained, not any similarity. The comparison is size-matched: at  $k=31$ , a minimizer window of  $w=15$  and closed syncmers with  $s=16$  select the same expected fraction of  $k$ -mers, so the two indexes have the same size and the contrast isolates the scheme. The decision rule — deacon's absolute and relative hit thresholds — is left unchanged, so the contrast is not confounded by tuning. And the expected effect size is small and error-dependent, so the test must sweep sequencing error.

Given that specification, an agent produced a patch adding a --scheme option to deacon's index builder, with a versioned index header so that an index records which scheme built it, in roughly 800 lines against a pinned upstream commit; all 107 upstream tests pass with the patch applied, and a smoke test with reads mutated at 1% and 2% per base retains 3,979 versus 3,969 and 3,867 versus 3,804 of 4,000 host reads for syncmers versus minimizers, with zero of 4,000 random reads retained under either scheme [126]. The harness then decides what the smoke test cannot: whether the gain survives across error models, read lengths, and a pangenome index, and whether it exists at all at 0.1% error, where a pre-analysis simulation predicts none. Three features of the exercise generalize. The hypothesis, the matched-size rule, and the metric were fixed before the code; the agent's work was verified against the tool's own test suite and an independent smoke test, not by inspection; and the benchmark, not the smoke test, is the arbiter of the claim.

## 7.4 Worked example 2: one predicate, a family of tools

The motif predicate of §2.3.4 has been specified once — an enzyme, a tag definition, and a canonicalization — and reused across four estimands: Fast2bRAD-M profiles species from tags [120], Syn2bANI estimates ANI and structural variation from tag anchors [119], Strain2bScan resolves strains [122], and sk2bGrow infers growth dynamics from tag coverage along the chromosome [123]. Two lessons follow. First, a validated primitive is a compounding asset: the tag definition was checked against the wet lab once — the tags a digest produces are the tags the software predicts — and every tool built on it inherits that check. Second, the primitive's constraints compound as well: every tool in the family inherits the fixed density and composition dependence of §3.7, and no downstream estimator repairs a tag-poor genome. When we built the mechanism row of §5 as a tool-free sampler from published REBASE cut offsets, it reproduced the tag lengths of the family's extractor for all fifteen supported enzymes and passed a strand-invariance test on every one — a test oracle independent of the tool being tested, which is what §7.2 asks for.

## 7.5 Failure modes specific to agent-built tools

Building the harness of §5 with agent assistance surfaced five failure modes, each of which we now guard against explicitly.

Plausible but wrong interfaces. Agents produce command-line flags and output column names that look right for the tool and version they were trained on. Every wrapper in the harness carries a verification marker until its interface has been checked against the installed version, and the minimal-mode run of §5.7 exists to surface this drift before the full grids are submitted.

Silent estimand drift. Asked for a distance, an agent returns the Mash distance; asked for detection, it returns Jaccard. The specification of §7.2 names the estimand, and the standard output tables of §5.7 force it into a named column.

Tests that cannot fail. An agent will write tests that pass on the code it wrote. The valuable tests are those with answers known independently: synthetic fixtures with exact truth, invariances (strand, order, merge), and comparison against a second implementation. One such fixture found that our T1 scoring had collapsed every mutant of a genome onto a single identifier — a naming bug invisible to inspection and fatal to every result.

Tuning on the truth. A tool under development must not read the benchmark's truth files, and defaults reported as "recommended" must have been set before the tool met the benchmark. Rule 1's size sweep exists so that a default tuned to one benchmark is visible as a point on a curve rather than as a headline.

Unmatched comparisons. "More accurate" at default parameters compares sketch sizes as much as methods (§5.1, Rule 1); "faster" compares thread counts and file systems (Rule 4). Both are cheap to control and, in our experience, the first thing an agent-generated benchmark omits.

Box 4. Release checklist for a sketch tool.

| Item | What to state or provide |
|---|---|
| Estimand | One sentence naming the quantity (Jaccard, containment,  $J_{w}$ ,  $J_{P}$ , membership); the column of Table 2 the tool occupies |
| Guarantee | Unbiased or not; CI or not; the mutation model and the assumptions behind any converted quantity (§3.6) |
| Selection rule | Context-free or not; size class; how density is set and what it should be for the smallest intended input (§3.1) |
| Hash and  $k$  | Hash function, seed, canonicalization; compatibility with other tools' sketches or explicit incompatibility |
| Algebra | Merge and subtraction semantics, or a statement that none are defined (§2.4.4) |
| Test oracle | A simulation with exact truth and a real-data anchor, shipped with the tool and runnable by others |
| Reporting | Raw estimate and converted quantity both emitted, with  $k$  and sampling parameters (§3.6) |
| Format | Version-stamped sketch files that record scheme, parameters and hash function (§5, Rule 5) |
| Benchmark row | A wrapper for a shared harness, so the tool is scored on the same truth as its comparators |
| Interests | For tools benchmarked by their authors, the design date of the benchmark and a competing-interests statement |

# 8. Open questions

The benchmark protocol of §5 was designed backward from a list of questions that neither theory nor published benchmarks currently answer. Nine such questions recur across the taxonomy. Each is stated below in falsifiable form, together with the specific figure, table, or design block of §5 intended to resolve it.

Q1. Confidence interval (CI) coverage. Do nominal 95% intervals achieve their stated coverage in practice? FracMinHash-based estimators ship analytic intervals for mutation rate [18], and sylph derives containment intervals from a per-k-mer model [25], but each was validated only on its authors' own simulated mixtures; no published measurement exists of realized coverage across estimator families on shared reference data. Figure 3 reports realized versus nominal coverage for every estimator that emits an interval.

Q2. $J_{w}$  versus  $J_{P}$  disagreement. How large is the gap between the two quantities marketed as "weighted Jaccard"? ProbMinHash collides with probability equal to the probability Jaccard similarity  $J_{P}$ [38], whereas streaming-histogram sketches such as HULK target the weighted-set Jaccard  $J_{w}$ [40], and the two coincide only under restrictive conditions on the abundance distributions. The T3 cross-scoring block of Table 5 scores every tool claiming "weighted Jaccard" against both ground truths, turning a labeling ambiguity into a measured quantity.

Q3. Empirical boundary of the safe scale factor. Does the error-controlled scale-factor formula hold on real metagenomes? Hera and Koslicki derived a lower bound on the FracMinHash scale factor needed to keep distance error within a target tolerance and showed that the widely used default  $s=1000$  can violate it at deep evolutionary distances (Algorithms for Molecular Biology 20:8) [19]; the derivation, however, assumes a Poisson mutation model, and its boundary has never been measured on communities with heterogeneous coverage. The T2 design sweeps scale factors on CAMI-style mixtures, and the corresponding panel of Figure 2 plots observed error against the predicted boundary.

Q4. Depth sensitivity. At what sequencing depth do weighted sketches stabilize? Multiplicity-aware estimators are motivated by abundance-weighted comparison [40], yet no published measurement exists of the read depth at which their estimates stop drifting for a fixed community composition. The T3 depth sweep (1M, 5M, and 20M reads per sample) estimates each tool's stabilization point directly.

Q5. Repeat robustness. Which estimators degrade least on repeat-rich genomes? ReSkmer showed that explicitly modeling repeats repairs population-genomic distances that standard k-mer estimators get wrong [84], and a recent preprint derives repeat-robust mutation-rate estimators (preprint) [85], but neither study benchmarks mainstream sketching tools, so the cross-family ranking is unknown. The repeat-fraction stratification of the T1 genomes (§5.8, planned) partitions the benchmark by genomic repeat content and reports per-bin accuracy; until it runs, T1 reports accuracy per ANI bin only.

Q6. End-to-end cost of sublinear sketches. Does the asymptotic win of sublinear-size sketches survive full-system accounting? MaxGeomHash (preprint; accepted at RECOMB 2026) proves that sketch size can grow sublinearly with database size at fixed containment error [28], but its reported advantages rest on sampling-time arguments, and build, merge, and I/O costs have not been benchmarked against fixed-size baselines. The T1/T2 timing protocol measures wall-clock build and query time and on-disk size alongside accuracy (Rule 4), so the crossover point — if one exists — becomes an empirical number rather than an asymptotic claim.

Q7. Protein-level sketching. Does the seven-family taxonomy transfer to amino-acid alphabets? Early evidence is suggestive: FracMinHash containment has been used to estimate dN/dS (preprint) [29], indicating that scaled sketches remain statistically tractable away from nucleotide data, but no systematic treatment of protein-alphabet estimators exists. This question and Q9 are the two that §5 does not answer; it is recorded here as a scope extension, since the alphabet change alters k-mer space size, mutation dynamics, and the validity of every family-specific estimator assumption at once.

Q8. Motif versus hash at matched density. Does non-uniform sampling cost accuracy? The motif predicate and the hash threshold can be compared at equal tags per genome — BcgI sits near  $s≈2000$  on random sequence (§2.3.4), and the enzyme-density scan of §5.2 fixes the matched scale factor per genome. T1 and T2 report both predicates on the same curves; any gap is the price of enzymatic executability, and its sign is not known in advance, since motif tags are deterministic and duplicate-free where hash samples are not.

Q9. Concordance of computed and sequenced sketches. Does the in-silico tag set of a sample agree with the tag set a 2bRAD-M library of the same sample yields? The motif family is the only one in which the sketch can be measured, and mock communities sequenced both ways exist [117,118]. The concordance bounds the wet-lab noise — incomplete digestion, amplification, sequencing error — that every downstream estimator inherits, and it is the one question in this list that no amount of computation settles. It is recorded as a scope extension (§5.8).

# 9. Conclusions

## 9.1 From a tool directory to an estimator decision framework

K-mer sketching is no longer one method with many implementations but a design space with two independent axes — the selection rule and the estimator — and only the second axis determines which biological question a sketch can answer. This review organized the field along that axis and delivered the corresponding artifacts: Figure 1 and Tables 1 and 2 define the seven families and state their properties in falsifiable form, and place restriction-tag sketching in that space as a third context-free predicate beside hash thresholds and syncmers; Table 4 dates the transfer between sketch theory and genomics application; Table 3 converts the taxonomy into an eligibility matrix over four common tasks, with a stated reason for every exclusion; §5 turns the thin-evidence cells of that matrix into an executable, version-locked benchmarking protocol with a cost axis for sketches that are sequenced rather than computed; and §7 turns the same taxonomy into a specification a tool builder can hand to a coding agent. The takeaway reduces to one sentence: choose — or build — a sketch by the estimator it implements and the task at hand, not by the tool name. The silent failures of §3, from empty read-level sketches at default scale factors to a biased minimizer Jaccard estimator [60], all trace back to the reverse choice. What remains is measurement, and the open questions of §8 specify exactly which measurements, from realized confidence-interval coverage to the concordance of a computed sketch with a sequenced one.

## 9.2 Limitations

Three limitations bound the conclusions. First, the scope is nucleotide k-mer sketches; protein-level sketching enters only as an open question (§8, Q7), because the alphabet change alters k-mer space size, mutation dynamics, and every family-specific estimator assumption at once. Second, the benchmark results are pending: §5 is a protocol with placeholders, and the empirical claims it is designed to settle will land when the compute runs complete; until then the eligibility matrix rests on theory and on the tools' own validations. Third, the motif-defined family is developed by several of the authors. We have addressed this by fixing the benchmark design before any tag-based result existed, by scoring the mechanism apart from the tools, by reporting no-estimate rates and the Jaccard-route bias rather than only the favourable containment route, and by stating the interest below; readers should nonetheless weigh the family's rows with that in mind.

## Competing interests

[Authors] develop 2bRAD-M, Fast2bRAD-M, Syn2bANI, Strain2bScan and sk2bGrow. The benchmark protocol (§5) and its implementation [126] were committed before any result for these tools existed, and the design and the results will be released together with their version history. The other tools benchmarked are not developed by the authors.

## Data and code availability

The benchmark pipeline, tool wrappers, simulation scripts, the tool-free tag sampler, the deacon syncmer patch, and the tests described in §7 are available at https://github.com/HuangShiLab/sketch-benchmark [126] under the MIT license. Sketches, configurations, logs, and result tables will be archived on Zenodo with a DOI per release when the benchmark has run.

# References

[1] Broder AZ. On the resemblance and containment of documents. In: Compression and Complexity of SEQUENCES 1997. IEEE; 1997. p. 21–29. doi:10.1109/SEQUEN.1997.666900

[2] Bar-Yossef Z, Jayram T, Kumar R, Sivakumar D, Trevisan L. Counting distinct elements in a data stream. In: Proceedings of RANDOM 2002 (Lecture Notes in Computer Science, vol 2483). Springer; 2002. p. 1–10. doi:10.1007/3-540-45726-7_1

[3] Cohen E, Kaplan H. Summarizing data using bottom-k sketches. In: Proceedings of PODC 2007. ACM; 2007. p. 225–234. doi:10.1145/1281100.1281133

[4] Charikar MS. Similarity estimation techniques from rounding algorithms. In: Proceedings of STOC 2002. ACM; 2002. p. 380–388. doi:10.1145/509907.509965

[5] Flajolet P, Fusy É, Gandouet O, Meunier F. HyperLogLog: the analysis of a near-optimal cardinality estimation algorithm. In: Proceedings of the 2007 Conference on Analysis of Algorithms (DMTCS Proceedings AH). 2007. p. 137–156. doi:10.46298/dmtcs.3545

[6] Ertl O. New cardinality estimation algorithms for HyperLogLog sketches. arXiv. 2017. arXiv:1702.01284 (preprint).

[7] Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM. Mash: fast genome and metagenome distance estimation using MinHash. Genome Biol. 2016;17:132. doi:10.1186/s13059-016-0997-x

[8] Li P, König AC. b-Bit minwise hashing. In: Proceedings of WWW 2010. ACM; 2010. p. 671–680. doi:10.1145/1772690.1772759

[9] Li P, König AC. Theory and applications of b-bit minwise hashing. Commun ACM. 2011;54(8):101–109. doi:10.1145/1978542.1978566

[10] Li P, Owen AB, Zhang CH. One permutation hashing. In: Advances in Neural Information Processing Systems 25 (NeurIPS 2012). 2012. p. 3122–3130.

[11] Zhao X. BinDash, software for fast genome distance estimation on a typical personal laptop. Bioinformatics. 2019;35(4):671–673. doi:10.1093/bioinformatics/bty651

[12] Zhao J, Zhao X, Pierre-Both J, Konstantinidis KT. BinDash 2.0: new MinHash scheme allows ultra-fast and accurate genome search and comparisons. bioRxiv. 2024. doi:10.1101/2024.03.13.584875 (preprint).

[13] Bovee R, Greenfield N. Finch: a tool adding dynamic abundance filtering to genomic MinHashing. J Open Source Softw. 2018;3(22):505. doi:10.21105/joss.00505

[14] Lees JA, Harris SR, Tonkin-Hill G, et al. Fast and flexible bacterial genomic epidemiology with PopPUNK. Genome Res. 2019;29(2):304–316. doi:10.1101/gr.241455.118

[15] Koslicki D, Zabeti H. Improving MinHash via the containment index with applications to metagenomic analysis. Appl Math Comput. 2019;354:206–215. doi:10.1016/j.amc.2019.02.018

[16] Ondov BD, Starrett GJ, Sappington A, Kostic A, Koren S, Buck CB, Phillippy AM. Mash Screen: high-throughput sequence containment estimation for genome discovery. Genome Biol. 2019;20:232. doi:10.1186/s13059-019-1841-x

[17] Irber LC, Brooks PT, Reiter TE, Pierce-Ward NT, Hera MR, Koslicki D, Brown CT. Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers. bioRxiv. 2022. doi:10.1101/2022.01.11.475838 (preprint).

[18] Rahman Hera M, Pierce-Ward NT, Koslicki D. Deriving confidence intervals for mutation rates across a wide range of evolutionary distances using FracMinHash. Genome Res. 2023;33(7):1061–1068. doi:10.1101/gr.277651.123

[19] Rahman Hera M, Koslicki D. Estimating similarity and distance using FracMinHash. Algorithms Mol Biol. 2025;20:8. doi:10.1186/s13015-025-00276-8

[20] Brown CT, Irber L. sourmash: a library for MinHash sketching of DNA. J Open Source Softw. 2016;1(5):27. doi:10.21105/joss.00027

[21] Irber L, Pierce-Ward NT, Abuelanin M, et al. sourmash v4: a multitool to quickly search, compare, and analyze genomic and metagenomic data sets. J Open Source Softw. 2024;9(98):6830. doi:10.21105/joss.06830

[22] Irber L, Pierce-Ward NT, Brown CT. Sourmash branchwater enables lightweight petabyte-scale sequence search. bioRxiv. 2022. doi:10.1101/2022.11.02.514947 (preprint).

[23] Liu S, Koslicki D. CMash: fast, multi-resolution estimation of k-mer-based Jaccard and containment indices. Bioinformatics. 2022;38(Suppl 1):i28–i35. doi:10.1093/bioinformatics/btac237

[24] LaPierre N, Alser M, Eskin E, Koslicki D, Mangul S. Metalign: efficient alignment-based metagenomic profiling via containment min hash. Genome Biol. 2020;21:242. doi:10.1186/s13059-020-02159-0

[25] Shaw J, Yu YW. Rapid species-level metagenome profiling and containment estimation with sylph. Nat Biotechnol. 2025;43(8):1348–1359. doi:10.1038/s41587-024-02412-y

[26] Koslicki D, White S, Ma C, Novikov A. YACHT: an ANI-based statistical test to detect microbial presence/absence in a metagenomic sample. Bioinformatics. 2024;40(2):btae047. doi:10.1093/bioinformatics/btae047

[27] Rahman Hera M, Liu S, Wei W, Rodriguez JS, Ma C, Koslicki D. Metagenomic functional profiling: to sketch or not to sketch? Bioinformatics. 2024;40(Suppl 2):ii165–ii173. doi:10.1093/bioinformatics/btae397

[28] Rahman Hera M, Koslicki D, Martínez C. MaxGeomHash: an algorithm for variable-size random sampling of distinct elements. bioRxiv. 2025. doi:10.1101/2025.11.11.687920 (preprint).

[29] Rodriguez JS, Rahman Hera M, Koslicki D. Leveraging FracMinHash containment for genomic dN/dS. bioRxiv. 2025. doi:10.1101/2025.11.12.688019 (preprint).

[30] Chen X, et al. Pilea: profiling bacterial growth dynamics from metagenomes with sketching. Microbiome. 2026. doi:10.1186/s40168-026-02374-0

[31] Yi H, Lin Y, Lin C, Jin W. Kssd: sequence dimensionality reduction by k-mer substring space sampling enables real-time large-scale datasets analysis. Genome Biol. 2021;22:84. doi:10.1186/s13059-021-02303-4

[32] Manasse M, McSherry F, Talwar K. Consistent weighted sampling. Microsoft Research Technical Report MSR-TR-2010-73; 2010.

[33] Ioffe S. Improved consistent sampling, weighted Minhash and L1 sketching. In: Proceedings of ICDM 2010. IEEE; 2010. p. 246–255. doi:10.1109/ICDM.2010.80

[34] Li P. 0-Bit consistent weighted sampling. In: Proceedings of KDD 2015. ACM; 2015. p. 665–674. doi:10.1145/2783258.2783406

[35] Ertl O. BagMinHash – minwise hashing algorithm for weighted sets. In: Proceedings of KDD 2018. ACM; 2018. p. 1368–1377. doi:10.1145/3219819.3220089

[36] Christiani T. DartMinHash: fast sketching for weighted sets. arXiv. 2020. arXiv:2005.11547 (preprint).

[37] Moulton R, Jiang Y. Maximally consistent sampling and the Jaccard index of probability distributions. In: Proceedings of ICDM 2018. IEEE; 2018. p. 347–356. doi:10.1109/ICDM.2018.00050

[38] Ertl O. ProbMinHash – a class of locality-sensitive hash algorithms for the (probability) Jaccard similarity. IEEE Trans Knowl Data Eng. 2022;34(7):3491–3506. doi:10.1109/TKDE.2020.3021176

[39] Yang D, Li B, Rettig L, Cudre-Mauroux P. HistoSketch: fast similarity-preserving sketching of streaming histograms with concept drift. In: Proceedings of ICDM 2017. IEEE; 2017. p. 545–554. doi:10.1109/ICDM.2017.64

[40] Rowe WPM, Carrieri AP, Alcon-Giner C, et al. Streaming histogram sketching for rapid microbiome analytics. Microbiome. 2019;7:40. doi:10.1186/s40168-019-0653-2

[41] Zhao J, Both JP, Rodriguez-R LM, Konstantinidis KT. GSearch: ultra-fast and scalable genome search by combining K-mer hashing with hierarchical navigable small world graphs. Nucleic Acids Res. 2024;52(16):e74. doi:10.1093/nar/gkae609

[42] Yu YW, Weber GM. HyperMinHash: MinHash in LogLog space. IEEE Trans Knowl Data Eng. 2022;34(1):328–339. doi:10.1109/TKDE.2020.2981311

[43] Ertl O. SetSketch: filling the gap between MinHash and HyperLogLog. Proc VLDB Endow. 2021;14(11):2244–2257. doi:10.14778/3476249.3476276

[44] Baker DN, Langmead B. Dashing: fast and accurate genomic distances with HyperLogLog. Genome Biol. 2019;20:265. doi:10.1186/s13059-019-1875-0

[45] Baker DN, Langmead B. Genomic sketching with multiplicities and locality-sensitive hashing using Dashing 2. Genome Res. 2023;33(7):1218–1227. doi:10.1101/gr.277655.123

[46] Mohamadi H, Khan H, Birol I. ntCard: a streaming algorithm for cardinality estimation in genomics data. Bioinformatics. 2017;33(9):1324–1330. doi:10.1093/bioinformatics/btw832

[47] Ertl O. UltraLogLog: a practical and more space-efficient alternative to HyperLogLog for approximate distinct counting. Proc VLDB Endow. 2024;17(7):1655–1668. doi:10.14778/3654621.3654632

[48] Ertl O. ExaLogLog: space-efficient and practical approximate distinct counting up to the exa-scale. In: Proceedings of the 28th International Conference on Extending Database Technology (EDBT 2025). OpenProceedings.org; 2025. p. 829–841.

[49] Karppa M, Pagh R. HyperLogLogLog: cardinality estimation with one log more. In: Proceedings of KDD 2022. ACM; 2022. p. 753–761. doi:10.1145/3534678.3539246

[50] Marçais G, DeBlasio D, Pandey P, Kingsford C. Locality-sensitive hashing for the edit distance. Bioinformatics. 2019;35(14):i127–i135. doi:10.1093/bioinformatics/btz354

[51] Zhang T, Yin Z, Xu X, et al. RabbitSketch: a high-performance sketching library for genome analysis. Bioinformatics. 2025;41(5):btaf249. doi:10.1093/bioinformatics/btaf249

[52] Chen K, Pattar V, Shao M. Sequence similarity estimation by random subsequence sketching. In: 25th International Workshop on Algorithms in Bioinformatics (WABI 2025), LIPIcs vol 344. 2025. p. 7:1–7:16. doi:10.4230/LIPIcs.WABI.2025.7

[53] Joudaki A, Rätsch G, Kahles A. Fast alignment-free similarity estimation by tensor sketching. bioRxiv. 2020. doi:10.1101/2020.11.13.381814 (preprint).

[54] Joudaki A, Meterez A, Mustafa H, Groot Koerkamp R, Kahles A, Rätsch G. Aligning distant sequences to graphs using long seed sketches. Genome Res. 2023;33(7):1208–1217. doi:10.1101/gr.277659.123

[55] Chakraborty D, Goldenberg E, Koucký M. Streaming algorithms for embedding and computing edit distance in the low distance regime. In: Proceedings of STOC 2016. ACM; 2016. doi:10.1145/2897518.2897577

[56] Chen K, Shao M. Locality-sensitive bucketing functions for the edit distance. Algorithms Mol Biol. 2023;18:7. doi:10.1186/s13015-023-00234-2

[57] McCauley S. Approximate similarity search under edit distance using locality-sensitive hashing. arXiv. 2019. arXiv:1907.01600 (preprint).

[58] Roberts M, Hayes W, Hunt BR, Mount SM, Yorke JA. Reducing storage requirements for biological sequence comparison. Bioinformatics. 2004;20(18):3363–3369. doi:10.1093/bioinformatics/bth408

[59] Schleimer S, Wilkerson DS, Aiken A. Winnowing: local algorithms for document fingerprinting. In: Proceedings of SIGMOD 2003. ACM; 2003. p. 76–85. doi:10.1145/872757.872770

[60] Belbasi M, Blanca A, Harris RS, Koslicki D, Medvedev P. The minimizer Jaccard estimator is biased and inconsistent. Bioinformatics. 2022;38(Suppl 1):i169–i176. doi:10.1093/bioinformatics/btac244

[61] Edgar R. Syncmers are more sensitive than minimizers for selecting conserved k-mers in biological sequences. PeerJ. 2021;9:e10805. doi:10.7717/peerj.10805

[62] Dutta A, Pellow D, Shamir R. Parameterized syncmer schemes improve long-read mapping. bioRxiv. 2022. doi:10.1101/2022.01.10.475696 (preprint).

[63] Sahlin K. Effective sequence similarity detection with strobemers. Genome Res. 2021;31(11):2080–2094. doi:10.1101/gr.275648.121

[64] Jain C, Dilthey A, Koren S, Aluru S, Phillippy AM. A fast approximate algorithm for mapping long reads to large reference databases. In: RECOMB 2017 (Lecture Notes in Computer Science, vol 10229). Springer; 2017. p. 66–81. doi:10.1007/978-3-319-56970-3_5

[65] Jain C, Koren S, Dilthey A, Phillippy AM, Aluru S. A fast adaptive algorithm for computing whole-genome homology maps. Bioinformatics. 2018;34(17):i748–i756. doi:10.1093/bioinformatics/bty597

[66] Kille B, Garrison E, Treangen TJ, Phillippy AM. Minmers are a generalization of minimizers that enable unbiased local Jaccard estimation. Bioinformatics. 2023;39(9):btad512. doi:10.1093/bioinformatics/btad512

[67] Groot Koerkamp R, Pibiri GE. The mod-minimizer: a simple and efficient sampling algorithm for long k-mers. In: 24th International Workshop on Algorithms in Bioinformatics (WABI 2024), LIPIcs vol 312. 2024. p. 11:1–11:23. doi:10.4230/LIPIcs.WABI.2024.11

[68] Groot Koerkamp R, Liu D, Pibiri GE. The open-closed mod-minimizer algorithm. Algorithms Mol Biol. 2025;20:4. doi:10.1186/s13015-025-00270-0

[69] Kille B, Groot Koerkamp R, McAdams D, Liu A, Treangen TJ. A near-tight lower bound on the density of forward sampling schemes. Bioinformatics. 2025;41(1):btae736. doi:10.1093/bioinformatics/btae736

[70] Groot Koerkamp R, Martayan I. SimdMinimizers: computing random minimizers, fast. In: 23rd International Symposium on Experimental Algorithms (SEA 2025), LIPIcs vol 338. 2025. p. 20:1–20:14. doi:10.4230/LIPIcs.SEA.2025.20

[71] Golan S, Tziony I, Kraus M, Orenstein Y, Shur A. GreedyMini: generating low-density DNA minimizers. Bioinformatics. 2025;41(Suppl 1):i275–i284. doi:10.1093/bioinformatics/btaf251

[72] Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094–3100. doi:10.1093/bioinformatics/bty191

[73] Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biol. 2019;20:257. doi:10.1186/s13059-019-1891-0

[74] Sahlin K. Strobealign: flexible seed size enables ultra-fast and accurate read alignment. Genome Biol. 2022;23:260. doi:10.1186/s13059-022-02831-7

[75] Constantinides B, Lees J, Crook DW. Deacon: fast sequence filtering and contaminant depletion. bioRxiv. 2025. doi:10.1101/2025.06.09.658732 (preprint).

[76] Nunes I, Heddes M, Vergés P, et al. DotHash: estimating set similarity metrics for link prediction and document deduplication. In: Proceedings of KDD 2023. ACM; 2023. p. 1758–1769. doi:10.1145/3580305.3599314

[77] Xu W, Hsu PK, Moshiri N, Yu S, Rosing T. HyperGen: compact and efficient genome sketching using hyperdimensional vectors. Bioinformatics. 2024;40(7):btae452. doi:10.1093/bioinformatics/btae452

[78] Cumbo F, Dhillon K, Najafi MH, Aygun S, Blankenberg D. HyperSketch: de Bruijn graph sketching for genomic similarity estimation with hyperdimensional computing. bioRxiv. 2026. doi:10.64898/2026.09.01.748726 (preprint).

[79] Manku GS, Jain A, Das Sarma A. Detecting near-duplicates for web crawling. In: Proceedings of WWW 2007. ACM; 2007. p. 141–150. doi:10.1145/1242572.1242592

[80] Cumbo F, Dhillon K, Joshi J, Chicco D, Aygun S, Blankenberg D. A novel vector-symbolic architecture for graph encoding and its application to viral pangenome-based species classification. BioData Min. 2026;19:54. doi:10.1186/s13040-026-00561-1

[81] Shaw J, Yu YW. Fast and robust metagenomic sequence comparison through sparse chaining with skani. Nat Methods. 2023;20(11):1661–1665. doi:10.1038/s41592-023-02018-3

[82] Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. Nat Commun. 2018;9:5114. doi:10.1038/s41467-018-07641-9

[83] Sarmashghi S, Bohmann K, Gilbert MTP, Bafna V, Mirarab S. Skmer: assembly-free and alignment-free sample identification using genome skims. Genome Biol. 2019;20:34. doi:10.1186/s13059-019-1632-4

[84] Charvel E, Thomas I, Alves Monteiro HJ, Sarmashghi S, Dunshea G, Bafna V, Mirarab S. ReSkmer: modeling repeats allows k-mer-based alignment-free methods to calculate population genomic distances. Genome Biol. 2026;27:233. doi:10.1186/s13059-026-04108-9

[85] Wu H, Medvedev P. The gift of novelty: repeat-robust k-mer-based estimators of mutation rates. bioRxiv. 2026. doi:10.64898/2026.04.01.715966 (preprint).

[86] Benoit G, Peterlongo P, Mariadassou M, et al. Multiple comparative metagenomics using multiset k-mer counting. PeerJ Comput Sci. 2016;2:e94. doi:10.7717/peerj-cs.94

[87] Benoit G, Mariadassou M, Robin S, Schbath S, Peterlongo P, Lemaitre C. SimkaMin: fast and resource frugal de novo comparative metagenomics. Bioinformatics. 2020;36(4):1275–1276. doi:10.1093/bioinformatics/btz685

[88] Murray KD, Webers C, Ong CS, Borevitz J, Warthmann N. kWIP: the k-mer weighted inner product, a de novo estimator of genetic similarity. PLoS Comput Biol. 2017;13(10):e1005727. doi:10.1371/journal.pcbi.1005727

[89] Deorowicz S, Gudyś A, Długosz M, Kokot M, Danek A. Kmer-db: instant evolutionary distance estimation. Bioinformatics. 2019;35(1):133–136. doi:10.1093/bioinformatics/bty610

[90] Zielezinski A, et al. Ultrafast and accurate sequence alignment and clustering of viral genomes. Nat Methods. 2025;22:1191–1194. doi:10.1038/s41592-025-02701-7

[91] Derelle R, von Wachsmann J, Maklin T, et al. Seamless, rapid, and accurate analyses of outbreak genomic data using split k-mer analysis (ska2). Genome Res. 2024;34(10):1661–1673. doi:10.1101/gr.279449.124

[92] Rowe WPM. When the levee breaks: a practical guide to sketching algorithms for processing the flood of genomic data. Genome Biol. 2019;20:199. doi:10.1186/s13059-019-1809-x

[93] Marçais G, Solomon B, Patro R, Kingsford C. Sketching and sublinear data structures in genomics. Annu Rev Biomed Data Sci. 2019;2:93–118. doi:10.1146/annurev-biodatasci-072018-021156

[94] Zielezinski A, Girgis HZ, Bernard G, et al. Benchmarking of alignment-free sequence comparison methods. Genome Biol. 2019;20:144. doi:10.1186/s13059-019-1755-7

[95] Zheng H, Marçais G, Kingsford C. Creating and using minimizer sketches in computational genomics. J Comput Biol. 2023;30(12):1251–1276. doi:10.1089/cmb.2023.0094

[96] Ndiaye M, Prieto-Baños S, Fitzgerald LM, et al. When less is more: sketching with minimizers in genomics. Genome Biol. 2024;25:270. doi:10.1186/s13059-024-03414-4

[97] Ponsero AJ, Miller M, Hurwitz BL. Comparison of k-mer-based de novo comparative metagenomic tools and approaches. Microbiome Res Rep. 2023;2(4):27. doi:10.20517/mrr.2023.26

[98] Majidian S, Hwang S, Zakeri M, Langmead B. EvANI benchmarking workflow for evolutionary distance estimation. Brief Bioinform. 2025;26(3):bbaf267. doi:10.1093/bib/bbaf267

[99] Penedo G, Kydlíček H, Ben Allal L, Lozhkov A, Mitchell M, Raffel C, Von Werra L, Wolf T. The FineWeb datasets: decanting the web for the finest text data at scale. arXiv. 2024. arXiv:2406.17557 (preprint).

[100] Son Y, Kim C, Lee J. SEDD: scalable and efficient dataset deduplication with GPUs. arXiv. 2025. arXiv:2501.01046 (preprint).

[101] Zhu E. datasketch: big data looks small. Zenodo; 2017. doi:10.5281/zenodo.598238

[102] Brown TB, et al. Language models are few-shot learners. In: Advances in Neural Information Processing Systems 33 (NeurIPS 2020). 2020. p. 1877–1901.

[103] Penedo G, Malartic Q, Hesslow D, Cojocaru R, Cappelli A, Alobeidli H, Pannier B, Almazrouei E, Launay J. The RefinedWeb dataset for Falcon LLM: outperforming curated corpora with web data, and web data only. arXiv. 2023. arXiv:2306.01116 (preprint).

[104] Parks DH, Chuvochina M, Rinke C, Mussig AJ, Chaumeil PA, Hugenholtz P. GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. Nucleic Acids Res. 2022;50(D1):D785–D794. doi:10.1093/nar/gkab776

[105] Meyer F, Fritz A, Deng ZL, et al. Critical Assessment of Metagenome Interpretation: the second round of challenges. Nat Methods. 2022;19(4):429–440. doi:10.1038/s41592-022-01431-4

[106] Human Microbiome Project Consortium. Structure, function and diversity of the healthy human microbiome. Nature. 2012;486(7402):207–214. doi:10.1038/nature11234

[107] Gourlé H, Karlsson-Lindsjö O, Hayer J, Bongcam-Rudloff E. Simulating Illumina metagenomic data with InSilicoSeq. Bioinformatics. 2019;35(3):521–522. doi:10.1093/bioinformatics/bty630

[108] Nicholls SM, Quick JC, Tang S, Loman NJ. Ultra-deep, long-read nanopore sequencing of mock microbial community standards. GigaScience. 2019;8(5):giz043. doi:10.1093/gigascience/giz043

[109] Richter M, Rosselló-Móra R. Shifting the genomic gold standard for the prokaryotic species definition. Proc Natl Acad Sci USA. 2009;106(45):19126–19131. doi:10.1073/pnas.0906412106

[110] Pritchard L, Glover RH, Humphris S, Elphinstone JG, Toth IK. Genomics and taxonomy in diagnostics for food security: soft-rotting enterobacterial plant pathogens. Anal Methods. 2016;8(1):12–24. doi:10.1039/C5AY02550H

[111] Zook JM, Catoe D, McDaniel J, et al. Extensive sequencing of seven human genomes to characterize benchmark reference materials. Sci Data. 2016;3:160025. doi:10.1038/sdata.2016.25

[112] Constantinides B. Panhuman DNA minimizer index for Deacon version 1 (format version 3). Zenodo; 2025. doi:10.5281/zenodo.17288185

[113] Koslicki D, Liu S. CMash: fast and accurate set similarity estimation via containment min hash (repository deprecated; README directs users to sourmash and YACHT). GitHub. https://github.com/dkoslicki/CMash

[114] Zhu E. datasketch v2.0.0 release notes: default MinHash permutation scheme changed to affine32. GitHub; 2026. https://github.com/ekzhu/datasketch/releases/tag/v2.0.0

[115] Both JP. hnsw_rs: Rust implementation of HNSW approximate nearest-neighbor search. GitHub. https://github.com/jean-pierreBoth/hnswlib-rs

[116] Wang S, Meyer E, McKay JK, Matz MV. 2b-RAD: a simple and flexible method for genome-wide genotyping. Nat Methods. 2012;9(8):808–810. doi:10.1038/nmeth.2023

[117] Sun Z, Huang S, Zhang M, Zhu Q, Haiminen N, Carrieri AP, Vázquez-Baeza Y, Parida L, Kim HC, Knight R, Liu YY. Species-resolved sequencing of low-biomass or degraded microbiomes using 2bRAD-M. Genome Biol. 2022;23:36. doi:10.1186/s13059-021-02576-9

[118] Geng Q, Lv J, Sun Z, Huang Y, et al. Integrated 2bRAD-M approach for comprehensive and cost-efficient metagenomic profiling of challenging environmental and biomedical specimens. Nat Protoc. 2026; in revision.

[119] HuangShiLab. Syn2bANI: strain-level ANI estimation and structural comparison via fixed restriction-site anchors. Manuscript in preparation; code: https://github.com/HuangShiLab/Syn2bANI

[120] HuangShiLab. Fast2bRAD-M: a Rust reimplementation of the 2bRAD-M profiling pipeline. https://github.com/HuangShiLab/Fast2bRAD-M

[121] MAP2B: the computational engine of 2bRAD-M profiling with machine-learning false-positive removal; see [118] for the current description. [verify whether a standalone MAP2B citation exists]

[122] HuangShiLab. Strain2bScan: strain-level scanning from 2b-RAD tags. Manuscript in preparation; code: https://github.com/HuangShiLab/Strain2bScan

[123] HuangShiLab. sk2bGrow: growth-dynamics prediction from 2b-RAD tag coverage. Manuscript in preparation; https://github.com/HuangShiLab/sk2bGrow-paper

[124] Blanca A, Harris RS, Koslicki D, Medvedev P. The statistics of k-mers from a sequence undergoing a simple mutation process without spurious matches. J Comput Biol. 2022;29(2):155–168. doi:10.1089/cmb.2021.0431

[125] Roberts RJ, Vincze T, Posfai J, Macelis D. REBASE: a database for DNA restriction and modification: enzymes, genes and genomes. Nucleic Acids Res. 2023;51(D1):D629–D630. doi:10.1093/nar/gkac975

[126] HuangShiLab. sketch-benchmark: a mechanism-first benchmark harness for k-mer sketching tools. 2026. https://github.com/HuangShiLab/sketch-benchmark
