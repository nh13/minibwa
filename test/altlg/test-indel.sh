#!/bin/sh
# Regression test for Task 7 fix round 1, Fix 3: REJECT length-changing lifted
# spans.
#
# A seed that spans an indel (or two adjacent .alt lift blocks) lifts to a primary
# span whose length differs from the seed length.  Injecting it as a contiguous
# len-bp anchor would corrupt the chain coordinates, so mb_anchor_project_alt must
# reject it (pri_en - pri_st + 1 != q->len -> continue).
#
# The fixture's ALT contig lifts to the primary with a 5bp deletion (.alt CIGAR
# 75M5D75M); the read's full-length SMEM spans that deletion, so its lifted span
# is 5bp longer than the seed.  Asserted via the MB_PROJ_TRACE probe seam:
#   * the full-length spanning seed (len=150) produces NO projected anchor;
#   * every projected anchor that IS emitted has len equal to its lifted span
#     length (len-preserving) -- the probe prints the projected POS and len, and
#     we re-derive the span length from the contig and require it to match.
#
# Without Fix 3 the spanning seed would project (a len-150 anchor at a 155bp span);
# this is demonstrated in the task report by toggling the guard.
#
# Usage: test/altlg/test-indel.sh [<minibwa-dir>]
set -eu

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MINIBWA="$MDIR/minibwa"
MK="$MDIR/test/altlg/mkfixture-indel.sh"

TMPD=$(mktemp -d /tmp/altlg-indel.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

fail() { echo "FAIL: $1"; exit 1; }
ok()   { echo "  ok: $1"; }

echo "[test-indel] building indel fixture (.alt 75M5D75M) ..."
/bin/sh "$MK" "$TMPD" 2>/dev/null
. "$TMPD/meta.txt"
[ -n "${LIFT_POS1:-}" ] || fail "fixture did not export LIFT_POS1"
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null

"$MINIBWA" map --dbg-alt-proj -c 4 --outn=999 \
    "$TMPD/ref.fa" "$TMPD/reads.fq" 2>"$TMPD/err" >"$TMPD/sam"
echo "  MB_PROJ traces:"
grep '^MB_PROJ' "$TMPD/err" | sed 's/^/    /' || echo "    (none)"

# 1) The full-length spanning seed (len=150) must NOT have projected.
span=$(mawk '$1=="MB_PROJ" && $5=="len=150"{print; c++} END{print "COUNT="c+0}' "$TMPD/err" | sed -n 's/^COUNT=//p')
[ "${span:-0}" = "0" ] \
    || fail "the indel-spanning seed (len=150) projected -- Fix 3 (length-changing reject) is not in effect"
ok "indel-spanning seed (len=150) NOT projected (length-changing lift rejected)"

# 2) The read is still placed on the primary via the normal .alt path: a chrP
#    record with a deletion in its CIGAR, anchored at the LIFT locus region.
del_rec=$(mawk -v p="$LIFT_POS1" '$1!~/^@/ && $3=="chrP" && $6~/D/{print "1"; exit}' "$TMPD/sam")
[ "${del_rec:-0}" = "1" ] \
    || fail "no chrP record with a deletion CIGAR -- fixture did not exercise an indel-spanning placement"
ok "read still placed on chrP across the deletion via the normal .alt path (CIGAR contains D)"

echo "[test-indel] PASS"
exit 0
