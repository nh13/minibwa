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

Large speedups (-5.5% WGS, -44.9% Hi-C) but every ablation level moves confident placements on real data (WGS L4: 1 in 4,300 at MAPQ>=10). The simulated ROC that passed was the wrong experiment. The Hi-C result survives as a different idea: skip rescue when the insert-size distribution is uninformative.

