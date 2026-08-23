#!/bin/sh
# Fixture for the reverse-multi-block and multi-record/hard-clip .alt lift paths
# (regression for the l2b_set_alt / l2b_lift coordinate bugs).  Every OTHER
# fixture emits a single full-length M block per ALT contig, which is exactly the
# coincidentally-correct regime; these two contigs are the cases that regime hid:
#
#   chrPrev_alt (150bp) : reverse record (FLAG 16) with an internal deletion,
#                         CIGAR 75M2000D75M @ chrP POS 201.  A reverse record's
#                         CIGAR walks the reverse-complement of the ALT contig, so
#                         the block coordinates need an RC->forward remap; before
#                         the fix the two blocks were swapped.
#   chrPhc_alt  (100bp) : two records for one ALT contig — a primary 50M @ POS
#                         1001 (FLAG 0) and a supplementary 50H50M @ POS 2001
#                         (FLAG 2048).  The leading hard clip is the offset into
#                         the ALT contig; before the fix H consumed nothing and
#                         the array was left unsorted, so the supplementary span
#                         lifted to the wrong place (or read as a hole).
#
# Emits ref.fa + ref.fa.alt only; the test drives ex-lift-check in check mode.
set -eu
. "$(dirname "$0")/fixlib.sh"
d="$1"; mkdir -p "$d"

# Sequence content is irrelevant to liftover (pure coordinate arithmetic from the
# .alt CIGARs); gen <seed> <length> comes from fixlib.sh.
{
    printf '>chrP\n%s\n'        "$(gen 71 3000)"
    printf '>chrPrev_alt\n%s\n' "$(gen 72 150)"
    printf '>chrPhc_alt\n%s\n'  "$(gen 73 100)"
} > "$d/ref.fa"

{
    printf 'chrPrev_alt\t16\tchrP\t201\t60\t75M2000D75M\t*\t0\t0\t*\t*\n'
    printf 'chrPhc_alt\t0\tchrP\t1001\t60\t50M\t*\t0\t0\t*\t*\n'
    printf 'chrPhc_alt\t2048\tchrP\t2001\t60\t50H50M\t*\t0\t0\t*\t*\n'
} > "$d/ref.fa.alt"
