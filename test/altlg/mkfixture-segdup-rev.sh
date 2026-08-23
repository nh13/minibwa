#!/bin/sh
# Fixture for Task 7 fix round 1: REVERSE-strand ALT-seed -> primary projection.
#
# Same segdup-recovery mechanism as mkfixture-segdup.sh, but the ALT contig is
# REVERSE-COMPLEMENT aligned to the primary (.alt FLAG 0x10).  This exercises the
# reverse code paths in mb_anchor_project_alt that the forward-only segdup fixture
# never touched:
#   * Fix 1: recover the ALT-contig forward span for a reverse-aligned block.
#   * Fix 2: emit the projected primary tpos in the FORWARD contig frame (do NOT
#            re-fold for the projected strand) so the lifted POS is correct.
#
# Reference layout on chrP (primary) -- mirrors the forward segdup fixture but the
# shared core is stored REVERSE-COMPLEMENTED on the primary:
#   [PAD0(200)] [RC(CORE)(100) + TAIL1(50)]          <- LIFT locus (1-based POS 201)
#   [PADi(70)]  [RC(CORE)(100) + TAILi(50)]  x (N-1) <- segdup decoy copies
#   [PADZ(200)]
#
#   The primary stores RC(CORE).  The read (= ALT contig) is CORE+ALTTAIL forward;
#   its CORE matches the primary copies only on the REVERSE strand.  Each copy has
#   a UNIQUE 50bp tail (TAIL1, TAILi, ...) so the read's full 150bp is NOT present
#   on the primary -- only the 100bp CORE is shared/repetitive.  With N primary
#   copies plus the ALT, that CORE seed's SA interval is size N+1; `-c <max_occ>`
#   with N+1 > max_occ subsamples it and drops the LIFT copy's own seed.
#
# ALT twin:
#   chrP_altS = CORE + ALTTAIL (150bp), forward.  Because ALTTAIL is unique to the
#   ALT contig, the read's full 150bp SMEM is UNIQUE to chrP_altS (SA size 1) -> the
#   ALT seed is NEVER subsampled and always produces an ALT anchor.
#   .alt maps chrP_altS onto chrP at the LIFT locus with FLAG 0x10 (RC), full
#   length (150M).  (Bases past CORE mismatch the primary's TAIL1; that is fine --
#   .alt alignments carry mismatches and the projected anchor is only a SEED.)
#
# Reverse-lift coordinate check (what Fix 1+2 must reproduce):
#   The .alt block is alt[0,150) -> primary[LIFT0, LIFT0+150), rev=1, where
#   LIFT0 = LIFT_POS1-1 (0-based).  Under l2b_lift's reverse rule
#   (pri = pri_en-1-(alt-alt_st)) the surviving full-read seed's forward ALT span
#   alt[0..149] lifts to primary span [LIFT0, LIFT0+149]; its FORWARD last base is
#   LIFT0+149, so mb_hit_set_coor's ts = (LIFT0+149)+1-150 = LIFT0 -> 1-based
#   POS = LIFT_POS1, on the REVERSE strand (FLAG&0x10).
#
# Expected projected primary record (asserted in test-segdup-rev.sh):
#   a REVERSE-strand chrP record at POS == LIFT_POS1.
#   BEFORE Fix 1+2 the projected POS is wrong (mirror-image / dropped);
#   AFTER, the reverse-strand chrP record appears at POS LIFT_POS1.
#
# Usage: mkfixture-segdup-rev.sh <out-dir> [N-copies]
set -eu
. "$(dirname "$0")/fixlib.sh"
d="$1"; mkdir -p "$d"
N="${2:-20}"


CORE=$(gen 100 100)        # 100bp shared seed core (segdup): on every copy + ALT
ALTTAIL=$(gen 777 50)      # 50bp tail UNIQUE to the ALT contig (makes ALT SMEM unique)
READ="${CORE}${ALTTAIL}"   # 150bp read == ALT contig sequence (forward)
LR=${#READ}
LC=${#CORE}
RCCORE=$(rc "$CORE")       # reverse complement of CORE: what the primary stores

# chrP: PAD0 + (LIFT copy: RC(CORE)+TAIL1) + decoys + PADZ
PAD0=$(gen 1 200)
TAIL1=$(gen 4001 50)       # the LIFT copy's own tail (differs from ALTTAIL)
CHRP="${PAD0}${RCCORE}${TAIL1}"
LIFT_OFF=${#PAD0}
LIFT_POS1=$(( LIFT_OFF + 1 ))    # 1-based POS where the .alt block (and lifted seed) lands

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

# .alt: chrP_altS -> chrP at the LIFT locus, RC (FLAG 0x10), full length (150M).
printf 'chrP_altS\t16\tchrP\t%d\t60\t%dM\t*\t0\t0\t*\t*\n' "$LIFT_POS1" "$LR" \
    > "$d/ref.fa.alt"

QUAL=$(printf "%${LR}s" '' | tr ' ' 'I')
printf '@r-segdup-rev\n%s\n+\n%s\n' "$READ" "$QUAL" > "$d/reads.fq"

echo "LIFT_POS1=$LIFT_POS1" > "$d/meta.txt"
echo "fixture(segdup-rev): chrP len=${#CHRP} N=$N core=$LC LIFT@${LIFT_POS1} RC ALT=chrP_altS(150)" >&2
