#!/bin/sh
# Fixture for Task 4 (reconciliation pass), case (b): the PARALOG GUARD.
#
# This is the BLOCKER fixture: it must stress the MB_LIFT_TOL boundary, not be
# trivially far apart.  Two near-identical primary loci whose ALT twins lift to
# lifted_st values JUST OUTSIDE MB_LIFT_TOL (=10) must NOT be merged into one
# group — the read must stay a low-MAPQ multi-mapper.
#
# Reference layout on chrP (primary):
#   [pad1 200bp][COPY1 (150bp)][gap 220bp][COPY2 (150bp, == COPY1)][pad2 200bp]
#   COPY1 1-based start = 201            (0-based 200)
#   COPY2 1-based start = 201+150+220=571 (0-based 570)
#   => the two paralog loci on chrP are 370bp apart.
#
# ALT twins (separate ALT contigs, each a single COPY-sized contig):
#   chrP_altA  = COPY (== COPY1 == COPY2).  .alt aligns it to chrP at POS=201
#               (1-based) => lifts to lifted_st ~= 200 (group of COPY1).
#   chrP_altB  = COPY.  .alt aligns it to chrP at POS=221 (1-based) =>
#               lifts to lifted_st ~= 220.  This is altA's lift (200) + 20,
#               i.e. exactly TWICE MB_LIFT_TOL away — just outside tolerance.
#
# A read from COPY maps to:
#   chrP@201   (COPY1, lifted_st ~200)
#   chrP@571   (COPY2, lifted_st ~570)
#   chrP_altA  (lifts to ~200)
#   chrP_altB  (lifts to ~220)
#
# Correct grouping with MB_LIFT_TOL=10:
#   {chrP@201 (~200), chrP_altA (~200)}  and  {chrP@571 (~570)} and
#   {chrP_altB (~220)} — altB is 20bp from the COPY1 group's 200, > tol, so it
#   does NOT merge.  The point: even the CLOSEST competing placements (~200 vs
#   ~220) stay distinct at the 2*tol boundary, so multiple genuine loci survive
#   => the read is a multi-mapper => MAPQ 0 (or very low).  If grouping were
#   sloppy and merged altB into the COPY1 group, we would still have COPY2 as a
#   distinct competitor, so the multi-mapper verdict holds regardless; the
#   assertion that proves the boundary is checked by the placement dump.
set -eu
. "$(dirname "$0")/fixlib.sh"
d="$1"; mkdir -p "$d"

# Unique-ish building blocks (LCG PRNG) so each region is internally unique but
# COPY1==COPY2 exactly (the paralog).

PAD1=$(gen 1 200)
COPY=$(gen 2 150)
GAP=$(gen  3 220)
PAD2=$(gen 4 200)

CHRP="${PAD1}${COPY}${GAP}${COPY}${PAD2}"
LEN=${#CHRP}

# ALT twins: each contig is just COPY (150bp).
printf '>chrP\n%s\n'      "$CHRP" > "$d/ref.fa"
printf '>chrP_altA\n%s\n' "$COPY">> "$d/ref.fa"
printf '>chrP_altB\n%s\n' "$COPY">> "$d/ref.fa"

# .alt: altA aligns to chrP at POS 201 (the COPY1 locus); altB at POS 221
# (20bp downstream — 2*MB_LIFT_TOL).  Both full-length 150M.
LC=${#COPY}
printf 'chrP_altA\t0\tchrP\t201\t60\t%dM\t*\t0\t0\t*\t*\n' "$LC"  > "$d/ref.fa.alt"
printf 'chrP_altB\t0\tchrP\t221\t60\t%dM\t*\t0\t0\t*\t*\n' "$LC" >> "$d/ref.fa.alt"

# read: the COPY sequence (matches both COPY1 and COPY2, and both ALT twins).
QUAL=$(printf "%${LC}s" '' | tr ' ' 'I')
printf '@r-para\n%s\n+\n%s\n' "$COPY" "$QUAL" > "$d/reads.fq"

echo "fixture(paralog): chrP=$LEN COPY1@201 COPY2@571 altA->201 altB->221 (Δ20 > 2*tol)" >&2
