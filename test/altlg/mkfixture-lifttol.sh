#!/bin/sh
# Fixture for the runtime-tunable ALT lift tolerance (--alt-lift-tol).
#
# Models a near-twin whose lifted placement sits 15bp away from where the read
# aligns directly on the primary -- i.e. lift/placement JITTER of 15bp (more than
# the default MB_LIFT_TOL=10, less than a real paralog separation).  Such jitter
# is what a too-tight tolerance fails to absorb: the twin and the primary fall
# into separate groups, compete, and the read is wrongly driven to MAPQ 0.
#
#   chrP            800bp unique sequence (LCG PRNG, seed 11).
#   read r-near     chrP[300,420) (1-based 301..420), 120bp.  Aligns directly to
#                   chrP at POS 301  => primary hit, lifted_st = 300 (0-based).
#   chrP_altT       EXACT copy of chrP[300,420) (120bp), so the read also aligns
#                   to it perfectly => an equal-scoring ALT hit.
#   .alt            chrP_altT aligned to chrP at POS 316 (120M).  The ALT hit
#                   therefore LIFTS to chrP 0-based 315 -- 15bp from the primary's
#                   300.  |Δlifted_st| = 15.
#
# Expected:
#   --alt-lift-tol 10 (default): 15 > 10 => twin is a SEPARATE group => equal-score
#       competitor => r-near MAPQ 0 (conservative default; under-merges jitter).
#   --alt-lift-tol 20          : 15 <= 20 => twin GROUPS with the primary => chrP is
#       sam_pri with MAPQ>0, the twin secondary (recovered).
set -eu
d="$1"; mkdir -p "$d"

# 800bp unique sequence (LCG PRNG, seed 11).
S=$(python3 -c "
import random
random.seed(11)
print(''.join(random.choice('ACGT') for _ in range(800)))
")

# ALT twin: EXACT copy of chrP[300,420) (1-based 301..420), 120bp.
ALTT=$(printf '%s' "$S" | cut -c301-420)
LT=${#ALTT}

# --- reference ---
printf '>chrP\n%s\n'      "$S"    > "$d/ref.fa"
printf '>chrP_altT\n%s\n' "$ALTT">> "$d/ref.fa"

# --- .alt: altT placed 15bp downstream of its true primary window (POS 316, not
#     301) so the lift lands at chrP 0-based 315 -- 15bp of jitter vs the read's
#     direct primary alignment at 300. ---
printf 'chrP_altT\t0\tchrP\t316\t60\t%dM\t*\t0\t0\t*\t*\n' "$LT" > "$d/ref.fa.alt"

# --- read ---
QUAL120=$(printf '%120s' '' | tr ' ' 'I')
RNEAR=$(printf '%s' "$S" | cut -c301-420)
printf '@r-near\n%s\n+\n%s\n' "$RNEAR" "$QUAL120" > "$d/reads.fq"

echo "fixture(lifttol): chrP=800  altT=chrP[300,420) placed +15bp (lifted_st jitter=15)" >&2
