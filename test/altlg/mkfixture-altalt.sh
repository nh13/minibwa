#!/bin/sh
# Fixture for the ALT-ALTERNATE fold-in (mb_reconcile_alt step 2b).
#
# Root cause it locks: an ALT contig's .alt CIGAR carries STRUCTURAL INDELS
# (the duplicated / divergent segments the alt-aligner could not co-linearly
# align).  Those indels OFFSET an ALT twin's lifted placement from the read's
# true primary copy by up to the indel size -- beyond MB_LIFT_TOL -- OR drop it
# into an insertion HOLE (unliftable).  The precise per-base lift then fails to
# co-locate the ALT twin with its primary hit, so the twin survives as a
# co-equal group representative: it both steals the SAM-primary slot (read
# placed on the ALT contig) AND dilutes the primary's MAPQ to 0.
#
# The fix recognizes -- via the .alt CORRESPONDENCE (this ALT contig is, by
# construction, an alternate of THIS primary region) -- that an ALT hit sharing
# read bases with a primary hit inside the ALT's own primary region is an
# ALTERNATE PLACEMENT of that locus, and folds it into the primary's group.
# Placement-based, not "discard all ALT": only ALT-vs-non-ALT folds, so two
# genuine primary loci never merge (paralog safety; see mkfixture-paralog.sh).
#
# Geometry (chrP = 2000bp unique sequence, LCG-free python random seed 23):
#
#   R-disp  = chrP[600,750)            -- true origin chrP:601 (1-based).
#   chrP_altF = R-disp ++ chrP[1000,1200)
#     .alt:  chrP_altF -> chrP POS 451, CIGAR 150M 400D 200M
#       150M : ALT[0,150)   -> chrP[450,600)   (R-disp's ALT hit lifts to chrP:451,
#                                               DISPLACED 150bp from its true 601)
#       400D : chrP[600,1000) skipped
#       200M : ALT[150,350) -> chrP[1000,1200)
#     contig primary span chrP[450,1200) COVERS the true origin chrP:601
#     (which sits in the 400D gap) -> fold-in applies.
#
#   R-hole  = chrP[1500,1650)          -- true origin chrP:1501 (1-based).
#   chrP_altG = chrP[1300,1450) ++ R-hole ++ chrP[1600,1700)
#     .alt:  chrP_altG -> chrP POS 1301, CIGAR 150M 150I 150D 100M
#       150M : ALT[0,150)   -> chrP[1300,1450)
#       150I : ALT[150,300) = R-hole -> INSERTION HOLE (unliftable)
#       150D : chrP[1450,1600) skipped
#       100M : ALT[300,400) -> chrP[1600,1700)
#     contig primary span chrP[1300,1700) COVERS the true origin chrP:1501
#     (in the 150D gap) -> fold-in applies even though the ALT hit is unliftable.
#
# Each read therefore has exactly two equal-scoring hits: the true chrP hit and
# an ALT hit that the per-base lift CANNOT co-locate with it (displaced / hole).
#
# Usage: mkfixture-altalt.sh <outdir>
set -eu
d="${1:?usage: mkfixture-altalt.sh <outdir>}"
mkdir -p "$d"

S=$(python3 -c "
import random
random.seed(23)
print(''.join(random.choice('ACGT') for _ in range(2000)))
")

sub() { printf '%s' "$S" | cut -c"$1"-"$2"; }

RDISP=$(sub 601 750)            # chrP[600,750)
ALTF=$(sub 601 750)$(sub 1001 1200)   # R-disp ++ chrP[1000,1200)
RHOLE=$(sub 1501 1650)         # chrP[1500,1650)
ALTG=$(sub 1301 1450)$(sub 1501 1650)$(sub 1601 1700)  # ++ R-hole ++

# --- reference ---
printf '>chrP\n%s\n'       "$S"    > "$d/ref.fa"
printf '>chrP_altF\n%s\n'  "$ALTF">> "$d/ref.fa"
printf '>chrP_altG\n%s\n'  "$ALTG">> "$d/ref.fa"

# --- .alt: structural-indel CIGARs that displace / hole the ALT twins ---
printf 'chrP_altF\t0\tchrP\t451\t60\t150M400D200M\t*\t0\t0\t*\t*\n'      > "$d/ref.fa.alt"
printf 'chrP_altG\t0\tchrP\t1301\t60\t150M150I150D100M\t*\t0\t0\t*\t*\n'>> "$d/ref.fa.alt"

# --- reads ---
QUAL150=$(printf '%150s' '' | tr ' ' 'I')
printf '@r-disp\n%s\n+\n%s\n' "$RDISP" "$QUAL150"  > "$d/reads.fq"
printf '@r-hole\n%s\n+\n%s\n' "$RHOLE" "$QUAL150" >> "$d/reads.fq"

echo "fixture(altalt): chrP=2000  altF=disp(150M400D200M)  altG=hole(150M150I150D100M)" >&2
