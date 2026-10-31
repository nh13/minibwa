#!/bin/sh
# Static guard for the index-threads-v2 feature.
#
# index-threads-v2 parallelizes the fused SSA-sampling / BWT-inversion loop in
# mb_bwt_libsais() with OpenMP. A downstream merge (e.g. against a feature that
# refactors the same region) can silently DROP that parallelization: the SA and
# BWT it produces are byte-identical either way, so the distribution's
# byte-identity gate passes while `dist` quietly ships a correct-but-slower index
# build. The byte-identity gate structurally cannot catch a same-output slowdown;
# this suite is the thing that does.
#
# It is a STATIC check (grep the assembled source), not a timing measurement, so
# it is deterministic and never perf-flaky. Run against the assembled checkout,
# whose path is passed as $1 (the convention the in-tree suites use).
set -eu

REPO="${1:-.}"
SRC="$REPO/index.c"
[ -f "$SRC" ] || { echo "FAIL: '$SRC' not found (assembled checkout path wrong?)"; exit 1; }

# Extract just mb_bwt_libsais so the pragmas we require are the ones inside THIS
# function, not some unrelated OpenMP loop elsewhere in the file.
fn="$(awk '/mb_bwt_t \*mb_bwt_libsais\(/{f=1} f{print} f&&/^}/{exit}' "$SRC")"
[ -n "$fn" ] || { echo "FAIL: mb_bwt_libsais() not found in index.c"; exit 1; }

# The parallel-for drives the fused loop; the atomic write is unique to that
# loop's single-primary store and is the tell that the parallel inversion (not
# just the seq fill) survived. Require both.
printf '%s\n' "$fn" | grep -q '#pragma omp parallel for' || {
	echo "FAIL: mb_bwt_libsais lost '#pragma omp parallel for' -- the index-threads-v2 SA-to-BWT parallelization was dropped in a merge"
	exit 1
}
printf '%s\n' "$fn" | grep -q '#pragma omp atomic write' || {
	echo "FAIL: mb_bwt_libsais lost the fused-inversion '#pragma omp atomic write' -- the parallel BWT inversion was dropped in a merge"
	exit 1
}

echo "PASS: index-threads-v2 OpenMP parallelization intact in mb_bwt_libsais"
