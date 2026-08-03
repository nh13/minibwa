# minibwa mate-rescue window-narrowing gate — ablation RESULTS (negative)

**Date:** 2026-08-01
**Question:** `pe.c:352` narrows the mate-rescue search window to `2*len` around the k-mer anchor
only when `max_ug >= 10 && max_ug >= len>>1 && n_good == 1`. bwa-mem3's ROC-validated k-mer rescue
narrows 94.8% of the time; measured on 14,077 real rescue decisions this gate narrows only 42.3%.
Loosening it is a large, well-defined speedup — is it free?
**Verdict:** **No. Rejected.** Every ablation level that narrows the window more aggressively is
faster and also moves confident placements on real data. The Hi-C speedup (−44.9% at the loosest
level) is real and large, but the correct fix is not a looser global window — it is skipping
rescue outright when the insert-size distribution is uninformative, which is a different, untried
idea.

**Worktree:** `perf/rescue-window-gate/` (branch `perf/rescue-window-gate`) — **do not merge.**
**Commits:** `9d0f121` (SMEM depth tunables, prerequisite plumbing) and `5985142` (the ablation
ladder itself).

## Harness / method

`pe.c:352`'s gate has three conditions, each independently blocking narrowing on real rescue
decisions (14,077 sampled):

| condition | blocks |
|---|---:|
| `max_ug >= 10` | 0.2% |
| `max_ug >= len>>1` | 70.2% |
| `n_good == 1` | 74.8% |

Only **42.3%** of rescues narrow under the stock gate.

An `MB_RESCUE_GATE` environment variable selects one of five ablation levels from a **single
build**, so no rung differs from any other by codegen or code layout — only by which branch of the
gate is taken at runtime:

| level | change |
|---|---|
| 0 | stock (all three conditions active) |
| 1 | drop `n_good == 1` |
| 2 | relax `len>>1` → `len>>2` |
| 3 | both (1) and (2) |
| 4 | keep only `max_ug >= 10` |

Level 0 was verified **md5-identical** to r416 (`ebc59ea`) stock output before any other level was
trusted.

Two arms were run: a simulated PE ROC (hg38, 150 bp, insert 500±100, 1,072,259 pairs) and real data
(WGS and Hi-C, 2,000,000 primary records each) compared against stock, order-balanced with
direction alternated each round.

## Finding 1 — mapping-phase speed

| level | simulated | real WGS | real Hi-C |
|---|---:|---:|---:|
| L1 | +0.5% | −1.47% | −6.80% |
| L2 | +1.0% | −0.42% | −1.94% |
| L3 | −0.2% | −3.99% | −16.24% |
| L4 | −1.3% | −5.54% | −44.92% |

Hi-C absolute: L0 164.9 s → L4 90.8 s. The simulated arm barely moves at any level; the real-data
Hi-C column is the one large effect.

## Finding 2 — the simulated ROC passed, but it was the wrong experiment

All five levels passed the simulated ROC cleanly: every MAPQ bin ≥10 carried identical mismap
counts across all five levels, on 1,072,259 pairs. That result does not generalize, because 150 bp
reads at insert 500±100 drawn directly from the reference are the regime where mate rescue barely
matters in the first place — there is little room for a looser rescue window to change anything,
so a pass here is close to a null test.

## Finding 3 — on real data, every level moves confident placements

Compared against stock over 2,000,000 primary records each:

| | speed | any move | MAPQ ≥1 | MAPQ ≥10 | MAPQ 60+ |
|---|---:|---:|---:|---:|---:|
| WGS L1 | −1.47% | 0.0083% | 1 in 19,600 | 1 in 333,000 | none |
| WGS L3 | −3.99% | 0.0174% | 1 in 9,200 | 1 in 36,400 | none |
| WGS L4 | −5.54% | 0.0946% | 1 in 1,350 | 1 in 4,300 | 1 in 667,000 |
| HiC L1 | −6.80% | 0.0150% | 1 in 9,900 | 1 in 105,000 | 1 in 2,000,000 |
| HiC L3 | −16.24% | 0.0256% | 1 in 5,900 | 1 in 51,300 | 1 in 2,000,000 |
| HiC L4 | −44.92% | 0.1069% | 1 in 1,460 | 1 in 5,390 | 1 in 222,000 |

There is no truth set on real data, so "moved" is not the same as "wrong." But a MAPQ ≥10
placement is a claim the aligner makes to the user, and every level tested makes some of those
claims differently from stock, including at MAPQ 60 for the loosest levels on both libraries.

WGS L3 is the one interesting cell: −3.99% speed, zero MAPQ 60+ movement, and the smallest
MAPQ≥10 movement rate of the non-trivial levels (1 in 36,400). It is noted here rather than
folded into the general rejection, since it is a materially different risk profile from L4 on
either library.

## Why it was abandoned

The gate as written blocks narrowing on 74.8% of real rescues via `n_good == 1` alone and another
70.2% via `len>>1`, so the 42.3% narrowing rate is much lower than bwa-mem3's 94.8%. That gap looks
like free performance sitting on the table. It is not: on real data, every ablation level that
closes the gap trades some MAPQ ≥10 placements against stock, and the sharpest levels move MAPQ 60
placements too — small in absolute rate, but there is no threshold tested here that reproduces
bwa-mem3's rate at zero movement. The simulated ROC that motivated the investigation could not
have caught this, because insert 500±100 read-from-reference simulation is not a regime where
rescue does much work either way.

## What survives

The Hi-C column is the finding worth keeping: −44.92% at L4 says the current gate spends enormous
time attempting rescues that were never going to succeed, because for Hi-C the mate is effectively
uniformly distributed across the genome — there is no informative insert-size window to narrow
around in the first place. The right fix is not to loosen the window gate globally for every
library type; it is to detect that the insert-size distribution is uninformative and skip rescue
in that regime. That costs WGS nothing (its insert-size distribution is informative, so the
detector would not fire) and is defensible on first principles rather than tuned to a benchmark.
This is tracked as a separate issue, not implemented here.

## Reproduction

```sh
cd perf/rescue-window-gate
make   # single build; level selected at runtime

# level 0 must be byte-identical to stock before trusting any other level:
MB_RESCUE_GATE=0 ./minibwa map -t12 hg38 wgs_1.fq wgs_2.fq > l0.sam
md5sum l0.sam   # compare against r416 (ebc59ea) stock output

# select an ablation level (1-4) and re-run the same command:
MB_RESCUE_GATE=4 ./minibwa map -t12 hg38 wgs_1.fq wgs_2.fq > l4.sam
```

Real-data arms: WGS and Hi-C, 2,000,000 primary records each, order-balanced against stock with
direction alternated each round to control for machine/thermal drift.
