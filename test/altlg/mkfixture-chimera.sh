#!/bin/sh
# Fixture for Task 4 (reconciliation pass), the CHIMERA case — locks Fix 1.
#
# Fix 1 tightened mb_reconcile_alt's group-merge predicate: two hits are the
# SAME liftover group only if they are the same read alignment on primary vs
# ALT — same pri_tid, same strand, |Δlifted_st| <= MB_LIFT_TOL, AND their query
# spans OVERLAP.  Before Fix 1 the merge ignored query overlap, so a chimeric
# read whose two DIFFERENT query segments happen to lift near the SAME primary
# locus would be (wrongly) collapsed into one group: the second segment's hit
# would be demoted to a same-group subordinate (secondary, 0x100), removed as a
# competitor, and the first segment's MAPQ inflated.
#
# This fixture constructs exactly that adversarial geometry:
#
#   chrP (primary) is one 900bp unique sequence (LCG PRNG, seed 11).
#
#   SEG1 = chrP[120,240)  (120bp)  -- the read's FIRST query segment.
#   SEG2 = chrP[560,680)  (120bp)  -- a DIFFERENT, distant region; the read's
#          SECOND query segment.  SEG1 and SEG2 share no homology (distinct
#          windows of a unique sequence), so they are genuine chimeric segments
#          of one read, not competitors for the same locus.
#
#   The read = SEG1 immediately followed by SEG2 (240bp).  Query span of the
#   SEG1 alignment is ~[0,120); of the SEG2 alignment ~[120,240) — DISJOINT.
#
#   ALT contig chrP_altE = a copy of chrP[560,760) (200bp) that BEGINS with SEG2
#   (so SEG2's hit starts at ALT offset ~0).  Its .alt line aligns chrP_altE to
#   chrP at POS=121 (1-based) — i.e. it LIES that this ALT contig corresponds to
#   the chrP region [120,320), the SAME primary window SEG1 lifts into.  Because
#   SEG2 sits at the very start of the ALT contig, SEG2's best hit (on chrP_altE)
#   lifts to lifted_st ~= 120 — within MB_LIFT_TOL of SEG1's own chrP hit
#   (lifted_st ~120) — but at a DISJOINT query span.
#
# Result:
#   SEG1 hit:  chrP,       query ~[0,120),   lifts to lifted_st ~120
#   SEG2 hit:  chrP_altE,  query ~[120,240), lifts to lifted_st ~120  (same!)
#
#   They are co-located (same pri_tid/strand, |Δlst| <= tol) but query-DISJOINT.
#   BEFORE Fix 1: merged -> SEG2's chrP_altE hit is secondary (0x100), competitor
#                 suppressed, SEG1 MAPQ inflated.
#   AFTER  Fix 1: NOT merged (query spans don't overlap) -> the chrP_altE hit
#                 stays an independent representative (NOT 0x100 secondary); the
#                 genuine second locus is preserved.
#
# The observable, fixture-stable signal asserted by the test is the SAM flag of
# the chrP_altE record: secondary (0x100 set) == wrongly merged (RED, pre-Fix1);
# NOT secondary == correctly kept separate (GREEN, post-Fix1).
set -eu
d="$1"; mkdir -p "$d"

# 900bp unique sequence (LCG PRNG, seed 11).
S=$(python3 -c "
import random
random.seed(11)
print(''.join(random.choice('ACGT') for _ in range(900)))
")

# SEG1 = chrP[120,240) 1-based 121..240; SEG2 = chrP[560,680) 1-based 561..680.
SEG1=$(printf '%s' "$S" | cut -c121-240)
SEG2=$(printf '%s' "$S" | cut -c561-680)

# ALT contig chrP_altE = chrP[560,760) 1-based 561..760 (200bp), contains SEG2.
ALTE=$(printf '%s' "$S" | cut -c561-760)
LE=${#ALTE}

# --- reference ---
printf '>chrP\n%s\n'      "$S"    > "$d/ref.fa"
printf '>chrP_altE\n%s\n' "$ALTE">> "$d/ref.fa"

# --- .alt: chrP_altE deliberately aligned to chrP at POS=121 (the SEG1 window)
#     so SEG2's ALT hit lifts onto SEG1's primary locus (lifted_st ~120). ---
printf 'chrP_altE\t0\tchrP\t121\t60\t%dM\t*\t0\t0\t*\t*\n' "$LE" > "$d/ref.fa.alt"

# --- read: SEG1, a short non-aligning SPACER, then SEG2.  The spacer (a window
#     of chrP far from both SEG1 and SEG2, reverse-complemented so it does not
#     extend either flank) guarantees the SEG1 and SEG2 alignments have a CLEAN
#     gap between their query spans — strictly DISJOINT, not merely abutting —
#     so mb_qspan_overlap is unambiguously false. ---
SPACER=$(printf '%s' "$S" | cut -c801-830 | rev | tr 'ACGTacgt' 'TGCAtgca')
RCHIM="${SEG1}${SPACER}${SEG2}"
LR=${#RCHIM}
QUAL=$(printf "%${LR}s" '' | tr ' ' 'I')
printf '@r-chimera\n%s\n+\n%s\n' "$RCHIM" "$QUAL" > "$d/reads.fq"

echo "fixture(chimera): chrP=900 SEG1=chrP[120,240) SEG2=chrP[600,720) altE=chrP[560,760)->POS121" >&2
