#!/bin/sh
# Fixture for Task 7 fix round 1, Fix 3: REJECT length-changing lifted spans.
#
# A seed that spans an indel (or two adjacent .alt lift blocks) lifts to a primary
# span whose length differs from the seed length.  Injecting it as a contiguous
# len-bp anchor would corrupt the chain coordinates, so mb_anchor_project_alt must
# REJECT it (pri_en-pri_st+1 != q->len -> continue).
#
# Construction:
#   chrP_altI is an ALT contig of LEFT(75) + RIGHT(75) = 150bp.  Its .alt CIGAR is
#   75M5D75M: the two 75bp blocks lift to primary with a 5bp DELETION from the ALT
#   between them, so the primary span of a seed covering BOTH blocks is 5bp LONGER
#   than the ALT seed.
#   The read == chrP_altI (150bp), unique to the ALT contig, so its full-length
#   SMEM (which SPANS the 5bp deletion, covering both blocks) survives and reaches
#   the projection.  Because that span lifts to length seed_len+5 != seed_len, the
#   projection must drop it -> NO MB_PROJ line for the spanning seed.
#
# Observable assertion (test-indel.sh, via the MB_PROJ_TRACE probe seam):
#   With Fix 3:    the full-length seed (len ~150, spanning the deletion) produces
#                  NO projected anchor (rejected).  Any projected anchor that DOES
#                  appear has len == its lifted span length (single-block seeds
#                  only).
#   Without Fix 3: the spanning seed would be injected as a len-bp anchor at a
#                  primary span that is actually len+5 bp -- a coordinate bug.
#
# Usage: mkfixture-indel.sh <out-dir>
set -eu
d="$1"; mkdir -p "$d"

gen() { python3 -c "
import random
random.seed($1)
print(''.join(random.choice('ACGT') for _ in range($2)))
"; }

LEFT=$(gen 11 75)          # 75bp block A (shared ALT<->primary)
RIGHT=$(gen 22 75)         # 75bp block B (shared ALT<->primary)
GAP=$(gen 33 5)            # 5bp present on primary only (D in the ALT->primary CIGAR)
ALT="${LEFT}${RIGHT}"      # 150bp ALT contig == read
LA=${#ALT}

# chrP (primary): PAD0 + LEFT + GAP + RIGHT + PADZ.  The .alt deletes GAP from the
# ALT, so the two ALT blocks lift to primary positions 5bp farther apart than they
# are on the ALT -> a seed covering both blocks changes length under the lift.
PAD0=$(gen 1 200)
PADZ=$(gen 5 200)
CHRP="${PAD0}${LEFT}${GAP}${RIGHT}${PADZ}"
LIFT_OFF=${#PAD0}
LIFT_POS1=$(( LIFT_OFF + 1 ))    # 1-based POS where block A (alt 0) lifts

printf '>chrP\n%s\n'      "$CHRP" > "$d/ref.fa"
printf '>chrP_altI\n%s\n' "$ALT"  >> "$d/ref.fa"

# .alt: chrP_altI -> chrP at LIFT_POS1, forward, 75M5D75M (5bp deletion from ALT).
printf 'chrP_altI\t0\tchrP\t%d\t60\t75M5D75M\t*\t0\t0\t*\t*\n' "$LIFT_POS1" \
    > "$d/ref.fa.alt"

QUAL=$(printf "%${LA}s" '' | tr ' ' 'I')
printf '@r-indel\n%s\n+\n%s\n' "$ALT" "$QUAL" > "$d/reads.fq"

echo "LIFT_POS1=$LIFT_POS1" > "$d/meta.txt"
echo "fixture(indel): chrP len=${#CHRP} alt=chrP_altI(${LA}) CIGAR=75M5D75M LIFT@${LIFT_POS1}" >&2
