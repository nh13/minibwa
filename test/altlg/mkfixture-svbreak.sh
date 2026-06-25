#!/bin/sh
# Fixture for SV-breakpoint-aware ALT grouping (the W2.2 hi_vs_z gap).
#
# An ALT contig whose .alt CIGAR carries an SV-scale indel describes a structural
# alignment in which two M blocks lift to primary positions thousands of bp apart.
# A read whose ALT footprint spans BOTH blocks (across the breakpoint) lifts to a
# placement whose min/max-COLLAPSED [lifted_st, lifted_en] is dragged thousands of
# bp from where the read's PRIMARY twin actually sits.  Under the old single-
# interval grouping that ALT hit fails to co-locate with its primary twin -> they
# compete -> a FALSE MAPQ 0.  The multi-interval placement exposes a sub-interval
# AT the twin's position, so they group and the read keeps its confident MAPQ.
#
# Construction (forward and reverse variants):
#   ALT contig chrP_altV = LEFT(75) + RIGHT(75) = the 150bp read.  The read maps
#   150M to chrP_altV (unique to that contig).
#
#   .alt CIGAR = 75M 2000D 75M (a 2000bp DELETION from the ALT, i.e. an SV-scale
#   indel): the LEFT block lifts to primary [200,275); the RIGHT block lifts to
#   primary [2275,2350) -- 2000bp away.  The read's full footprint spans both
#   blocks, so the OLD collapse yields lifted_st = 200 (the min).
#
#   The PRIMARY TWIN (LEFT+RIGHT contiguous, == the read) is placed at primary
#   offset 2275 -- exactly the RIGHT block's lifted position.  (.alt lifting is
#   purely coordinate-based, so the twin's bases need not equal the RIGHT block's
#   declared bases.)  The read maps 150M to the twin at 1-based POS 2276.
#
#   OLD behavior (RED):  ALT collapsed lifted_st = 200; twin lifted_st = 2275;
#                        |Δ| = 2075 >> MB_LIFT_TOL=10 -> no co-location -> compete
#                        -> primary MAPQ 0.
#   NEW behavior (GREEN): ALT sub-placements = {200 (LEFT), 2275 (RIGHT)}; the 2275
#                        sub-interval co-locates with the twin -> they group ->
#                        the primary twin is sam_pri with MAPQ>0, the ALT secondary.
#
#   The forward variant uses .alt FLAG 0 (read maps forward to the twin); the
#   reverse variant uses .alt FLAG 0x10 and a reverse-complemented twin (read maps
#   reverse to the twin), exercising the reverse-fold lift paths.
#
# Usage: mkfixture-svbreak.sh <out-dir>
set -eu
d="$1"; mkdir -p "$d"

gen() { python3 -c "
import random
random.seed($1)
print(''.join(random.choice('ACGT') for _ in range($2)))
"; }
rc() { python3 -c "
import sys
s=sys.argv[1]
print(s.translate(str.maketrans('ACGT','TGCA'))[::-1])
" "$1"; }

A=75            # LEFT block length
B=75            # RIGHT block length
N=2000          # SV-scale deletion from the ALT (the breakpoint gap on primary)

LEFT=$(gen 11 "$A")
RIGHT=$(gen 22 "$B")
READ="${LEFT}${RIGHT}"      # 150bp read == ALT contig (forward)
LR=${#READ}

# .alt LEFT block lifts to primary 0-based 200 (1-based POS 201); RIGHT block to
# 0-based 200+A+N.  Place the read twin (READ, contiguous) at that RIGHT position.
PAD0=$(gen 1 200)
P0=${#PAD0}                       # 200 (0-based start of the LEFT block)
TWIN0=$(( P0 + A + N ))           # 0-based start of the RIGHT block == twin start
ALT_POS1=$(( P0 + 1 ))            # 1-based POS the .alt block aligns LEFT at
TWIN_POS1=$(( TWIN0 + 1 ))        # 1-based POS the read twin maps at
JUNK_LEN=$(( TWIN0 - P0 ))        # filler between PAD0 and the twin
JUNK=$(gen 99 "$JUNK_LEN")
PADZ=$(gen 5 200)

QUAL=$(printf "%${LR}s" '' | tr ' ' 'I')

# ---- forward variant ----
#   chrP_fwd = PAD0 + JUNK + READ(twin, fwd) + PADZ ; .alt FLAG 0.
CHRP_FWD="${PAD0}${JUNK}${READ}${PADZ}"
printf '>chrP\n%s\n'      "$CHRP_FWD" >  "$d/ref.fa"
printf '>chrP_altV\n%s\n' "$READ"     >> "$d/ref.fa"
printf 'chrP_altV\t0\tchrP\t%d\t60\t%dM%dD%dM\t*\t0\t0\t*\t*\n' \
    "$ALT_POS1" "$A" "$N" "$B" > "$d/ref.fa.alt"
printf '@r-sv-fwd\n%s\n+\n%s\n' "$READ" "$QUAL" > "$d/reads.fq"

# ---- reverse variant (separate index dir) ----
#   chrP_rev = PAD0 + JUNK + rc(READ) (twin maps reverse) + PADZ ; .alt FLAG 0x10.
mkdir -p "$d/rev"
RCREAD=$(rc "$READ")
CHRP_REV="${PAD0}${JUNK}${RCREAD}${PADZ}"
printf '>chrP\n%s\n'      "$CHRP_REV" >  "$d/rev/ref.fa"
printf '>chrP_altV\n%s\n' "$READ"     >> "$d/rev/ref.fa"
printf 'chrP_altV\t16\tchrP\t%d\t60\t%dM%dD%dM\t*\t0\t0\t*\t*\n' \
    "$ALT_POS1" "$A" "$N" "$B" > "$d/rev/ref.fa.alt"
printf '@r-sv-rev\n%s\n+\n%s\n' "$READ" "$QUAL" > "$d/rev/reads.fq"

{
    echo "ALT_POS1=$ALT_POS1"
    echo "TWIN_POS1=$TWIN_POS1"
} > "$d/meta.txt"
echo "fixture(svbreak): .alt=${A}M${N}D${B}M  ALT block@${ALT_POS1}  twin@${TWIN_POS1} (Δ$(( N + A )) >> MB_LIFT_TOL)" >&2
