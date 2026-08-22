#!/bin/sh
# Fixture for Task 7 (ALT-seed -> primary anchor projection; segdup recovery).
#
# Demonstrates the core failure mode: a SEGDUPLICATED primary locus whose own
# seed is dropped by max_occ subsampling, so the ONLY way a primary candidate
# exists at that locus is by PROJECTING the surviving ALT-contig seed back onto
# the primary assembly (DRAGEN's mechanism).
#
# Reference layout on chrP (primary):
#   [PAD0(200)] [CORE(100) + TAIL1(50)]            <- LIFT locus (1-based POS 201)
#   [PADi(70)]  [CORE(100) + TAILi(50)]  x (N-1)   <- segdup decoy copies
#   [PADZ(200)]
#
#   * CORE (100bp) is byte-identical across the LIFT copy and every decoy copy
#     AND the ALT contig -> they all share the same 100bp seed.  With N copies on
#     the primary plus the ALT, that seed's SA interval has size N+1; run minibwa
#     with `-c <max_occ>` small enough (N+1 > max_occ) that the interval is
#     SUBSAMPLED.  The strided sampler drops the LIFT copy's own seed position.
#   * Each copy has a UNIQUE 50bp tail (TAIL1, TAILi, ...) so the only shared,
#     repetitively-seeded region is CORE.
#
# ALT twin:
#   chrP_altS = CORE + ALTTAIL  (150bp).  Because ALTTAIL is unique to the ALT
#   contig, the read's full 150bp SMEM is UNIQUE to chrP_altS (SA size 1) -> the
#   ALT seed is NEVER subsampled and always produces an ALT anchor.
#   .alt maps chrP_altS onto chrP at the LIFT locus, full length (150M).  (The
#   tail bases mismatch the primary's own TAIL1; that is fine -- .alt alignments
#   carry mismatches, and the projected anchor is only a SEED: DP re-aligns it.)
#
# Read = chrP_altS sequence (CORE + ALTTAIL, 150bp).
#
# Behaviour (asserted in test-segdup.sh, with `-c 4 -N 20`):
#   WITHOUT projection: the LIFT-locus primary candidate is absent (its CORE seed
#       was subsampled out) -> NO chrP record at POS 201.
#   WITH projection: the surviving ALT seed lifts to the LIFT locus and injects a
#       primary anchor there -> a chrP record at POS 201 appears.
#
# Usage: mkfixture-segdup.sh <out-dir> [N-copies]
set -eu
d="$1"; mkdir -p "$d"
N="${2:-20}"   # number of primary copies of CORE (>= the max_occ used in tests)

gen() { python3 -c "
import random
random.seed($1)
print(''.join(random.choice('ACGT') for _ in range($2)))
"; }

CORE=$(gen 100 100)        # 100bp shared seed core (segdup): on every copy + ALT
ALTTAIL=$(gen 777 50)      # 50bp tail UNIQUE to the ALT contig (makes ALT SMEM unique)
READ="${CORE}${ALTTAIL}"   # 150bp read == ALT contig sequence
LR=${#READ}
LC=${#CORE}

# chrP: PAD0 + (LIFT copy: CORE+TAIL1) + (decoy copies: PADi + CORE + TAILi) + PADZ
PAD0=$(gen 1 200)
TAIL1=$(gen 4001 50)       # the LIFT copy's own tail (differs from ALTTAIL)
CHRP="${PAD0}${CORE}${TAIL1}"
LIFT_OFF=${#PAD0}
LIFT_POS1=$(( LIFT_OFF + 1 ))    # 1-based POS of the LIFT locus (where .alt maps)

i=2
while [ "$i" -le "$N" ]; do
    PADI=$(gen $(( 2000 + i )) 70)
    TAILI=$(gen $(( 4000 + i )) 50)
    CHRP="${CHRP}${PADI}${CORE}${TAILI}"
    i=$(( i + 1 ))
done
PADZ=$(gen 5 200)
CHRP="${CHRP}${PADZ}"

printf '>chrP\n%s\n'      "$CHRP"  > "$d/ref.fa"
printf '>chrP_altS\n%s\n' "$READ" >> "$d/ref.fa"

# .alt: chrP_altS -> chrP at the LIFT locus, full length (150M).
printf 'chrP_altS\t0\tchrP\t%d\t60\t%dM\t*\t0\t0\t*\t*\n' "$LIFT_POS1" "$LR" \
    > "$d/ref.fa.alt"

QUAL=$(printf "%${LR}s" '' | tr ' ' 'I')
printf '@r-segdup\n%s\n+\n%s\n' "$READ" "$QUAL" > "$d/reads.fq"

echo "LIFT_POS1=$LIFT_POS1" > "$d/meta.txt"
echo "fixture(segdup): chrP len=${#CHRP} N=$N core=$LC LIFT@${LIFT_POS1} ALT=chrP_altS(150,unique tail)" >&2
