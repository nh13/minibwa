#!/bin/sh
# Fixture for test-projrev.sh: REVERSE-strand ALT-seed -> primary projection,
# SAM-observable via MAPQ differential (no probe required).
#
# MECHANISM
# ---------
# The read's full 150bp SMEM has SA=2: one occurrence on the ALT contig
# (forward) and one on the primary LIFT locus (reverse complement).  Running
# with -c 1 (max_occ=1) forces the strided sampler to take exactly ONE of the
# two hits from the SA interval.  The SA-array lexicographic order (determined
# by the sequence seeds below) places the ALT contig hit first, so:
#
#   WITHOUT projection (MB_NO_ALT_PROJECT=1):
#     Sampler takes the ALT-contig hit -> ALT anchor -> alignment to chrP_altS
#     only; no primary anchor at LIFT -> chrP_altS is SAM primary (MAPQ 60),
#     NO chrP record at LIFT_POS1.
#
#   WITH projection (default):
#     Sampler takes the ALT-contig hit -> ALT anchor -> mb_anchor_project_alt
#     injects a reverse-strand native primary anchor at LIFT_POS1 (folded_rev =
#     blk_rev ^ alt_rev = 1 ^ 0 = 1).  DP at LIFT: read maps perfectly (150M,
#     NM=0, score=300) because the primary stores RC(read) at LIFT_POS1.
#     Reconciliation groups {chrP_altS, chrP LIFT}: scores tie at 300; non-ALT
#     preference promotes chrP -> chrP is SAM primary (FLAG=16, MAPQ=60).
#     chrP_altS becomes secondary (FLAG=256, MAPQ=0).
#
# WHY THE EXISTING segdup-rev FIXTURE CAN'T DO THIS
# --------------------------------------------------
# In segdup-rev the LIFT locus stores RC(CORE)+TAIL1 where TAIL1 != RC(ALTTAIL),
# so the primary alignment has ~50 bp of mismatches -> DP score ~202 < ALT
# score 300 -> ALT wins as group rep -> chrP LIFT is secondary (MAPQ=0).
# Making TAIL1 = RC(ALTTAIL) (perfect match) would create junction k-mers
# with SA=2 that are always fully sampled (SA <= max_occ=4), bypassing
# projection entirely.  The trick here is using -c 1: the SA=2 full-read SMEM
# IS subsampled (SA > max_occ=1), and the SA-array ordering deterministically
# puts the ALT contig hit first for the chosen sequence seeds.
#
# REVERSE-STRAND PROJECTION CODE PATH
# ------------------------------------
# The ALT anchor is forward on chrP_altS (alt_rev=0).  The .alt block is RC
# (FLAG 0x10, blk_rev=1).  mb_anchor_project_alt computes:
#   folded_rev = blk_rev ^ alt_rev = 1 ^ 0 = 1   (reverse primary strand)
# This exercises the same Fix 1 + Fix 2 coordinate path tested in
# test-segdup-rev.sh (via MB_PROJ_TRACE), but here the result is directly
# observable in the SAM output.
#
# REFERENCE LAYOUT (chrP, primary)
# ---------------------------------
#   [PAD0(200bp)]  [LIFT: RC(ALTTAIL)(50bp) + RC(CORE)(100bp)]  [decoy copies x N]  [PADZ(200bp)]
#    <- LIFT_POS1=201 (1-based) ^
#
#   The LIFT region is RC(read), so the read aligns perfectly in REVERSE.
#   Decoy copies store RC(CORE)+TAILi (unique random tails), providing the
#   sub-SMEM anchors that are also subsampled by -c 1.
#
# Usage: mkfixture-projrev.sh <out-dir> [N-copies]
set -eu
. "$(dirname "$0")/fixlib.sh"
d="$1"; mkdir -p "$d"
N="${2:-20}"


# Sequence seeds chosen so the SA-array lexicographic order places the
# ALT-contig hit BEFORE the primary LIFT hit in the 150bp SMEM's SA interval,
# making max_occ=1 (-c 1) deterministically pick the ALT contig -> projection
# is the ONLY path to the primary LIFT locus.
CORE=$(gen 500 100)        # 100bp shared seed core (repetitive via RC on primary)
ALTTAIL=$(gen 5777 50)     # 50bp tail UNIQUE to the ALT contig (and RC is unique to LIFT locus)
READ="${CORE}${ALTTAIL}"   # 150bp read == ALT contig sequence (forward)
LR=${#READ}
RCCORE=$(rc "$CORE")       # RC(CORE): stored at each decoy and at the 2nd half of LIFT
RCALTTAIL=$(rc "$ALTTAIL") # RC(ALTTAIL): stored at the 1st half of LIFT only

# LIFT region = RC(ALTTAIL)(50bp) + RC(CORE)(100bp) = RC(read) -- perfect
# reverse match for the read.  The .alt FLAG=0x10 block maps ALT -> primary in
# RC, so primary[LIFT_POS1-1 .. LIFT_POS1+LR-2] = RC(ALT) = RC(read).
LIFT_REGION="${RCALTTAIL}${RCCORE}"

# chrP: PAD0 + LIFT + (decoy copies: PADi + RC(CORE) + TAILi) + PADZ
PAD0=$(gen 1 200)
CHRP="${PAD0}${LIFT_REGION}"
LIFT_OFF=${#PAD0}
LIFT_POS1=$(( LIFT_OFF + 1 ))    # 1-based POS of the LIFT locus

i=2
while [ "$i" -le "$N" ]; do
    PADI=$(gen $(( 2000 + i )) 70)
    TAILI=$(gen $(( 4000 + i )) 50)
    CHRP="${CHRP}${PADI}${RCCORE}${TAILI}"
    i=$(( i + 1 ))
done
PADZ=$(gen 5 200)
CHRP="${CHRP}${PADZ}"

printf '>chrP\n%s\n'      "$CHRP"  > "$d/ref.fa"
printf '>chrP_altS\n%s\n' "$READ" >> "$d/ref.fa"

# .alt: chrP_altS -> chrP at the LIFT locus, RC (FLAG 0x10), full length.
# The block strand (blk_rev=1) XOR the forward ALT-seed strand (alt_rev=0)
# gives folded_rev=1: the projected primary anchor is on the REVERSE strand.
printf 'chrP_altS\t16\tchrP\t%d\t60\t%dM\t*\t0\t0\t*\t*\n' "$LIFT_POS1" "$LR" \
    > "$d/ref.fa.alt"

QUAL=$(printf "%${LR}s" '' | tr ' ' 'I')
printf '@r-projrev\n%s\n+\n%s\n' "$READ" "$QUAL" > "$d/reads.fq"

echo "LIFT_POS1=$LIFT_POS1" > "$d/meta.txt"
echo "fixture(projrev): chrP len=${#CHRP} N=$N core=100 alttail=50 LIFT@${LIFT_POS1} RC ALT=chrP_altS(${LR}) SA_full=2 max_occ=1" >&2
