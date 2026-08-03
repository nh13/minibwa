# minibwa `--meth` confident wrong-chromosome tail — investigation RESULTS (negative)

**Date:** 2026-07-20
**Companion to:** `2026-06-25-minibwa-core-pe-noalt-investigation-results.md` (non-ALT WGS) and
`2026-06-25-minibwa-alt-pe-mapq-investigation-results.md` (ALT arm).
**Question:** why does minibwa place more reads on the **wrong chromosome at confident mapQ** than the
reference aligner on directional EM-seq, and is it fixable (placement and/or mapQ)?
**Verdict:** **Not a defect, and not fixable by mapQ.** minibwa's PE mapQ on this population is
well calibrated (empirical Q≈30.5 where it reports 20–49). On **total** error at mapQ>=20 minibwa
already beats the reference aligner (+82,093 correct, −602 errors). The wrong-chromosome excess is a
base-rate artifact of a 99.9%-correct population, not over-confidence. A candidate patch was built,
tuned and measured, then **rejected**: it is strictly dominated by simply raising the mapQ filter.

**Worktree:** `pe-mapq-repeat/` (branch `fix/pe-mapq-low-copy-repeats`, commit `cd88c2e`) — **do not merge.**
**Artifacts:** `/Volumes/scratch-00001/tmp/minibwa/` (SAMs, grades, feature dumps, `PR_DRAFT.md`).

## Harness / method

- 10.7M directional EM-seq pairs (holodeck sim, GRCh38 truth in read name), `minibwa map --meth -t12`,
  GATK hg38 index (decoys + ALTs). Grader: `analysis/tail_eval.py` / `tail_eval2.py` vs `golden.bam`.
- Real-data arm: **SRR4235788** (NA12878 directional WGBS, the run cited in `tex/minibwa.tex`), first
  5M pairs, 125 bp. Streamed ~0.8% of the 79 GB run via `curl | gzip -dc | head`.
- Non-`--meth` control: HG002 WGS-1M (Zenodo 19703025).
- Baseline is r409, which already contains r404's PE-mapQ rework and its `score_se2` cap.

## Baseline

| mapQ>= | wrongCHR | wrongPOS | correct |
|---|---:|---:|---:|
| 60 | 0   | 939   | 9,057,031 |
| 40 | 19  | 1,091 | 9,660,270 |
| 30 | 40  | 1,296 | 9,696,028 |
| 20 | 136 | 1,473 | 9,775,818 |

## Finding 1 — the errors are near-ties, and `frac_high` cannot see them

For all 136 confident wrong-chromosome ends, the score margin to the next-best locus:

| margin (best − 2nd, pts; 1 mismatch ≈ 10) | reads |
|---|---:|
| 0 (exact tie) | 53 |
| 10 | 69 |
| 20 | 9 |
| 30 | 5 |

122/136 (90%) are within one mismatch. Composition: ~23 X/Y, ~29 ALT/decoy/random, ~84 primary-chromosome
paralogs. The repeat down-weight that should catch these is driven by `frac_high`, and
`mb_cal_high_cov()` (map-algo.c:175) only counts a seed when `sai[i].size > max_occ` (**250**). Low-copy
repeats seed 2–20 times, so `frac_high == 0` for essentially all of them — instrumented `mb_pair` dumps
confirm `fh=0` and `n_sub=1` on these pairs. **The blind spot is real.**

## Finding 2 — but the population `frac_high` misses is 99.9% correct

Instrumented every pair where a low-copy-ambiguity term fires (630,590 pairs), labelled against truth,
and split the ones a demotion would move across mapQ 20:

- **A** = demoted **and** wrong-chromosome: **58 pairs / 124 ends**
- **B** = demoted **and** correct: **65,055 pairs / 149,187 ends**
- ratio **1 : 1,122**

`P(wrong | ambiguity fires) = 58/65,113 = 0.089%` → **Q ≈ 30.5**. minibwa reports 20–49 for these.
**It is already telling the truth.** Any rule that demotes them below 20 is a calibration *regression*,
which is exactly what the measured calibration showed (bin 10–19 observed Q moved 17.5 → 27.9, i.e. from
honest to badly under-confident).

## Finding 3 — nothing in the alignment separates A from B

AUC of A-vs-B for every internal feature (0.5 = coin flip):

| feature | AUC | feature | AUC |
|---|---:|---|---:|
| max_nhit | 0.612 | se_margin | 0.535 |
| min_nhit | 0.611 | mapq_raw | 0.528 |
| n_pp | 0.591 | min_end_margin | 0.523 |
| pair_margin | 0.567 | f_high | 0.504 |
| max_end_margin | 0.545 | min_mapq_se | 0.504 |

Cost curve for the best feature is **flat** — the hallmark of no signal:

| catch | sacrifices | cost |
|---|---|---|
| 50% of A | 46,416 B | 1,601 correct per wrong |
| 80% of A | 55,497 B | 1,196 per wrong |
| 95% of A | 61,164 B | 1,110 per wrong |

The "is there a competing *concordant pair*" hypothesis is refuted outright: 96.6% of wrong pairs have
one — and so do 97.7% of correct pairs.

## Finding 4 — the three external-information fixes also fail

Compared reference at the chosen locus vs the best XA alternative, in read orientation:

| hypothesis | A (wrong) | B (correct) | verdict |
|---|---:|---:|---|
| **H3** mean base qual @ discriminating pos | 25.8 | 26.9 | refuted |
| **H3** min base qual @ discriminating pos | 23.0 | 25.2 | refuted |
| **H2** margin *entirely* bisulfite-illusory | 21.8% | 19.9% | refuted |
| **H1** ALT/random/Un contig involved | 16.4% | 20.6% | refuted (*lower* in wrong) |

**H2** deserves a note because it was the most promising `--meth`-specific idea: under directional BS a
{C,T} difference (R1) or {G,A} (R2) between two loci cannot be resolved by the read, so such a margin is
illusory. Wrong and correct pairs carry the *same* fraction of illusory positions — it does not separate.
**H1** would cover only ~16% of the errors even if a `.alt` track were wired in.

The one real stratifier is same-chrom vs cross-chrom competitor (A 3.6% vs B 39.1% same-chrom), but it is
partly definitional (a wrong-*chromosome* read's true locus is cross-chromosome by construction) and even
its best stratum sits at 684 correct-per-wrong.

## Finding 5 — the premise was cherry-picked; minibwa beats the reference on total error

The original framing counted only wrong-chromosome. At mapQ>=20:

| aligner | correct | wrongPOS | wrongCHR | total wrong | err rate | obsQ |
|---|---:|---:|---:|---:|---:|---:|
| minibwa r409 | 9,775,818 | 1,473 | 136 | **1,609** | 1.65e-04 | **37.8** |
| reference aligner | 9,693,725 | 2,189 | 22 | **2,211** | 2.28e-04 | 36.4 |

minibwa has **+82,093 more correct and 602 fewer total errors**. The reference wins on wrong-chromosome
(22 vs 136) and loses on wrong-position (2,189 vs 1,473) — a difference in error *composition*, not
calibration quality. Global mapQ AUC is 0.9848 for both before and after any patch (Δ = −3e-5).

## The rejected patch

A minimal change was built and fully measured: per-end best-vs-second-best *locus* gap
(`dp_max_se − dp_max_se2`) ramped over `pe_rep_span` mismatches, combined across ends with `min()` so
mate rescue survives, folded into the `frac_high` factor. It works as designed — placements byte-identical,
mapQ monotonically lowered, demoted reads 99.7% XA-bearing and enriched 8.0x on `chrUn*` / 7.4x on
`*_random` on real SRR4235788 data. It is still **not worth taking**:

| | correct | total wrong |
|---|---:|---:|
| **patched**, filter mapQ>=20 | 9,626,759 | **1,336** |
| **unpatched**, filter mapQ>=30 | 9,696,028 | **1,336** |

Identical error count; the unpatched aligner keeps **69,269 more correct reads**. Anyone wanting a lower
error rate simply raises their threshold and comes out ahead — no aligner change required. `pe_rep_span`
only slides along this curve (span=2: 37 wrongCHR / 9,687,582 correct; span=3: 12 / 9,626,759).

## Conclusions

1. `frac_high`'s `max_occ` blind spot for low-copy repeats is **real** and worth recording, but closing it
   does not improve mapQ — the population behind the blind spot is 99.9% correctly placed.
2. The confident wrong-chromosome tail is a **base-rate artifact**: ~65k reads that are each ~99.9% likely
   right will produce ~58 wrong pairs by arithmetic alone. mapQ ≈ 30 is the honest label for them, and it
   is what minibwa already assigns.
3. **No mapQ-side fix exists.** Best separating feature AUC 0.61; cost pinned at ~1,100 correct sacrificed
   per wrong removed at every threshold.
4. This corroborates the 2026-06-25 non-ALT conclusion ("structural, not a tunable defect") from the other
   direction: that study found the truth copy *absent* from the candidate set 192/192 times; here the truth
   locus is present as a near-tie in 47% of cases and still unrecoverable. Both point at **seeding /
   candidate-set completeness**, not the mapQ formula.
5. Do not compare aligners on a single error class. On total error at a fixed threshold minibwa is the
   better-calibrated of the two here.

## If this is revisited

The only lever with real headroom is the **candidate set** (finding the true locus), not the scoring of
the candidates you already have — consistent with the June non-ALT study. A `.alt`/segdup annotation is the
standard remedy but caps out at ~16% of these errors. Anything score-based is bounded by Finding 3.

## Reproduction

```sh
D=/Volumes/scratch-00001/tmp/minibwa
minibwa map --meth -t12 --xa=200 $D/index/Homo_sapiens_assembly38.fasta \
  $D/data/place_r1.fq.gz $D/data/place_r2.fq.gz > out.sam
python3 $D/analysis/tail_eval2.py out.sam $D/data/golden.bam   # correct/wrongPOS/wrongCHR by mapQ
python3 $D/analysis/mapq_auc.py  $D/data/golden.bam before:a.sam after:b.sam   # AUC + calibration
python3 $D/analysis/sep2.py      $D/feat.tsv a.sam b.sam $D/data/golden.bam    # A-vs-B separability
python3 $D/analysis/why.py       a.sam b.sam $D/data/golden.bam $D/index/*.fasta  # H1/H2/H3
```

Real WGBS arm: `SRR4235788` (ENA `ftp.sra.ebi.ac.uk/vol1/fastq/SRR423/008/SRR4235788/`), stream the head
of both mates rather than fetching all 79 GB.
