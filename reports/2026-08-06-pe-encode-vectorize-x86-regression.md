# minibwa `pe-encode-vectorize` — x86 regression RESULTS (negative)

**Date:** 2026-08-06
**Question:** `pe-encode-vectorize` splits the paired-end query encode in `mb_matesw` into two
loops so the reverse-complement auto-vectorizes. The 2026-08-07 ladders show it `~flat` on arm64
and a cost on x86 (+0.36pp c6a16/gcc, +0.46pp c6a16/clang). Is that cost real, is it standalone,
and is there any regime where the feature pays?
**Verdict:** **No regime. Retired.** The cost is real and reproducible at 5–11 standard errors,
it is standalone, and the feature's own premise is sound but worthless: the region it optimizes
does not appear in a profile at all. Flat on the architecture it was written for, a measurable
cost on the other.

**Worktree:** `simd-revcomp/` (branch `perf/pe-encode-vectorize`, commit `8e1e79b`) — **do not
merge.**

## Method

No new measurement was needed. Three sources, all free:

1. The per-rep timings already published for the four 2026-08-07 ladder cells, at
   `s3://nh13-minibwa-bench/runs/ebc59eaf…/c96816dd075f9274/704544e0…/{gcc-15,clang-21}/{c8g16,c6a16}/map-t16/wgs-1M/bench.db`.
2. Both `pe.c` variants compiled for x86-64 **and** arm64 in the pinned image
   `704544e0…` (gcc 15.2 / clang 21.1) under the real ladder flags
   (`-std=c99 -O3 -Wall`, plus `-msse4.2 -mpopcnt` on x86), comparing vectorizer reports and
   disassembly.
3. A local `sample(1)` profile of the mapping phase (`wgs-1M`, `-t 8`, SAM to `/dev/null`).

## The cost is real

Pooling all 10 reps per rung. The `11 vs 9` column is a control: if rung 11 inherits the step,
the change at rung 10 is a durable level shift rather than a transient.

| cell | rung 9 | rung 10 | 10 vs 9 | Welch *t* | Cohen *d* | 11 vs 9 |
|---|---:|---:|---:|---:|---:|---:|
| gcc-15 / c8g16 | 8.7360 s | 8.7317 s | −0.049 pp | −1.35 | −0.60 | −0.031 pp |
| clang-21 / c8g16 | 8.5107 s | 8.5137 s | +0.035 pp | +1.36 | +0.61 | +0.061 pp |
| gcc-15 / c6a16 | 12.5739 s | 12.6160 s | **+0.335 pp** | **+5.25** | **+2.35** | +0.328 pp |
| clang-21 / c6a16 | 11.8025 s | 11.8658 s | **+0.536 pp** | **+11.26** | **+5.03** | +0.532 pp |

Per-rep CV on these rungs is 0.05–0.18%, so the standard error on a 10-rep mean is 0.02–0.06%
and the x86 shift is 6–11 SEs out. The earlier "comparable to rung 10's own 0.30% scatter"
caveat compared a shift in a *mean* against *per-rep* spread — the wrong denominator. arm64 is
flat with opposite signs across the two compilers, which is noise, not a small win.

## It is standalone, not an adjacency artifact

Rungs 8, 9 and 10 touch disjoint files — `format.c`, `map-main.c`, `pe.c` — and minibwa builds
without LTO, so rung 9 and rung 10 are separate translation units that cannot affect each
other's codegen. The feature's hunk also applies cleanly to `pe.c` at the ladder base
`ebc59eaf`, confirming rung 9's `pe.c` is stock. No manifest reordering was required.

## The premise is correct and the mechanism is the opposite of the one expected

The x86 vectorizer does **not** bail. It vectorizes exactly as the commit message claims, and
the arm64 reports are identical line for line:

```
gcc   rung 9    pe.c:421 basic block part vectorized using 8 byte vectors
gcc   rung 10   pe.c:429 loop vectorized using 16 byte vectors
                pe.c:429 loop versioned for vectorization because of possible aliasing
clang rung 9    pe.c:423 loop not vectorized
clang rung 10   pe.c:429 vectorized loop (vectorization width: 16, interleaved count: 2)
```

`mb_matesw` instruction counts grow comparably on both architectures — x86 gcc 453→618, clang
402→649; arm64 gcc 412→540, clang 355→497. The codegen change is architecture-neutral. What
differs is how the hardware handles it.

The split introduces a round trip through memory the fused loop did not have. Loop 1 is a
double-indirect table gather, never vectorized on either architecture, so it fills `qs[r][0]`
with **scalar 1-byte stores**. Loop 2 then immediately reads that same buffer back with
**16-byte unaligned vector loads**:

```
1f80: movdqu (%rax,%r8,1),%xmm5
1f86: movdqu 0x10(%rax,%r8,1),%xmm6
1fce: pshufb %xmm4,%xmm6
1fd3: movdqu %xmm6,-0xf(%rdx,%rdi,1)     # unaligned, cache-line-splitting
```

On x86 a wide load overlapping narrower stores still in the store buffer cannot be
store-forwarded; the load blocks until they drain to L1 (~12–19 cycles on Zen3), and it recurs
every vector iteration because loop 1 ran immediately before over the same ~150 bytes. Rung 9's
fused loop kept the value in a register and never round-tripped.

This accounts for all four cells and for the deleted `map-t1` axis: **c6a16** is Zen3 with 8
physical cores × 2 SMT threads, so siblings contend for the store buffer; **c8g16** is
Graviton4 with 16 physical cores, no SMT, and more tolerant forwarding; **`map-t1`** was flat
everywhere because one thread per core leaves the store buffer uncontended.

*Caveat:* this dataset cannot separate "x86" from "SMT" — c6a16 is both and c8g16 is neither.
Both readings point the same way and both describe the shipped configuration.

## There is no regime where it pays

A local profile of the mapping phase settles it. Self time, as a share of compute:

| self% | symbol |
|---:|---|
| 27.8 | `ksw_extd2_sse` |
| 17.8 | `ksw_ll_u8_core` |
| 11.1 | `mb_bwt_sa_batch` |
| 6.2 | `mb_bwt_rank1a` |
| 5.4 | `mb_bwt_smem_batch` |
| 4.4 | `mb_bwt_rank2a` |

SW kernels are 49.1% of compute and BWT/SMEM 27.0%. **The nt4 encode does not appear at all** —
no symbol matches `nt4` or `encode`, and the lookups inline into callers that are themselves
near-zero (`mb_map_sai` 0.21%, `process_batch` 0.15%, `kseq2bseq` 0.14%). Arithmetic agrees:
`wgs-1M` is ~300 MB of sequence translated roughly twice, which at 1–3 bytes/cycle/core across
16 threads is 4–13 ms of a 12.6 s run.

That the feature is *flat* rather than positive on arm64, where the vectorization succeeds and
costs nothing, is the tell: the region is too small a share of runtime for a local speedup to
surface. A regime with more mate rescue (Hi-C) amplifies the gain and the stall together, and
on x86 the stall dominates. More threads per core makes it strictly worse, so `map-t96` would
increase SMT pressure rather than relieve it.

Note the irony for anyone revisiting this: `ksw_ll_u8_core` at 17.8% *is* mate rescue — the
Farrar kernel called from `mb_matesw`, the same function this feature modified. The feature was
in the right function and optimized the one part of it that costs nothing.

## If someone retries this

The fix is not a loop split. It is a `pshufb`/`tbl` nt4 translation that vectorizes *loop 1*,
eliminating both the scalar store pass and the round trip — `s2n-lite.h` already provides the
SSE→NEON shim to write it once. But the profile above says the whole encode path is under 0.1%
of runtime, so that work is unmeasurable on this ladder and should not be undertaken for
performance reasons.
