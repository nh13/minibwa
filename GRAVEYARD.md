# Graveyard

Changes investigated and deliberately **not** shipped. They are recorded so the
same idea is not re-litigated, and their branches are kept so each disproof stays
reproducible: the branch is the artifact you can check out and re-measure.
Nothing here is a candidate for upstreaming.

The report path on each entry names the maintainer's own analysis notes. Those
are not part of this repository — the branch and this summary are what travel
with the distribution.

## `pe-mapq-low-copy-repeats`

- **Branch:** `fix/pe-mapq-low-copy-repeats`
- **Report:** `reports/2026-07-20-minibwa-meth-pe-mapq-wrongchrom-negative-result.md`

AUC flat at -2.6e-5 and strictly dominated by raising the mapQ threshold: unpatched at mapQ>=30 gives the same 1,336 errors as patched at mapQ>=20, with 69,269 more correct reads. The wrong-chromosome tail is a base-rate artifact, not a calibration defect.

## `rescue-window-gate`

- **Branch:** `perf/rescue-window-gate`
- **Report:** `reports/2026-08-01-mate-rescue-window-gate.md`

Large speedups (-5.5% WGS, -44.9% Hi-C) but every ablation level moves confident placements on real data (WGS L4: 1 in 4,300 at MAPQ>=10). The simulated ROC that passed was the wrong experiment. The Hi-C half is void too: that arm ran without --hic, which already skips all mate rescue, so it timed a configuration nobody should use for the assay (5.9-7.9x on 100k HG002 pairs; proper-paired 88,411 -> 0). The idea it was said to leave behind -- skip rescue when the insert-size distribution is uninformative -- is already implemented twice: --hic for the whole assay, and the per-read futility gate at pe.c:358, whose constants (10, 0.33) match bwa-mem3's --rescue-skip. See reports/2026-08-08-hic-mate-rescue-futility.md.

## `pe-encode-vectorize`

- **Branch:** `perf/pe-encode-vectorize`
- **Report:** `reports/2026-08-06-pe-encode-vectorize-x86-regression.md`

Flat on arm64 and a reproducible cost on x86 (+0.34pp gcc t=5.3, +0.54pp clang t=11.3), standalone rather than an adjacency artifact. The premise holds -- the revcomp does auto-vectorize, identically on both architectures -- but the split makes loop 2 reload with 16-byte vector loads what scalar loop 1 just wrote byte-by-byte, which defeats x86 store-to-load forwarding. Moot regardless: the nt4 encode is absent from a profile, under 0.1% of runtime, so nothing here was worth optimizing.

