#!/bin/sh
# Regression test for reverse-multi-block and multi-record/hard-clip .alt lifting.
# Drives ex-lift-check in check mode: each "<ctg> <pos> <expected>" triple asserts
# l2b_lift(ctg, pos) maps to primary <expected>, and that the contig's lift[] is
# sorted by alt_st (the invariant l2b_lift's binary search relies on).
#
# Expected primary coordinates (0-based), derived by hand from the .alt CIGARs:
#   chrPrev_alt (reverse, 75M2000D75M @ chrP POS 201 -> 0-based 200):
#     forward-ALT base 0   = RC base 149 -> primary 2349
#     forward-ALT base 74  = RC base 75  -> primary 2275
#     forward-ALT base 75  = RC base 74  -> primary 274
#     forward-ALT base 149 = RC base 0   -> primary 200
#   chrPhc_alt (50M @1001 + 50H50M @2001, both forward; POS 0-based 1000/2000):
#     base 0  -> 1000 ; base 49 -> 1049 ; base 50 -> 2000 ; base 99 -> 2049
# Usage: test/altlg/test-revmulti.sh [<minibwa-dir>]
set -eu

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MKFIXTURE="$MDIR/test/altlg/mkfixture-revmulti.sh"
MINIBWA="$MDIR/minibwa"
EX_LIFT="$MDIR/api-test/ex-lift-check"

TMPD=$(mktemp -d /tmp/altlg-test.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

echo "[test-revmulti] building fixture in $TMPD ..."
/bin/sh "$MKFIXTURE" "$TMPD"

echo "[test-revmulti] indexing with minibwa ..."
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null

echo "[test-revmulti] running ex-lift-check (check mode) ..."
# ex-lift-check exits non-zero on any MISMATCH or unsorted lift[].
OUT=$("$EX_LIFT" "$TMPD/ref.fa" \
        chrPrev_alt 0 2349 \
        chrPrev_alt 74 2275 \
        chrPrev_alt 75 274 \
        chrPrev_alt 149 200 \
        chrPhc_alt 0 1000 \
        chrPhc_alt 49 1049 \
        chrPhc_alt 50 2000 \
        chrPhc_alt 99 2049)
echo "$OUT" | sed 's/^/    /'

if echo "$OUT" | grep -q MISMATCH; then
    echo "FAIL: at least one lift MISMATCH (see above)"
    exit 1
fi
if echo "$OUT" | grep -q 'sorted=0'; then
    echo "FAIL: a contig's lift[] is not sorted by alt_st"
    exit 1
fi

echo "[test-revmulti] PASS: reverse-multi-block and hard-clip multi-record .alt lift correctly"
exit 0
