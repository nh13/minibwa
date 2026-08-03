#!/bin/sh
# Integration test for the runtime-tunable ALT lift tolerance (--alt-lift-tol).
#
# A near-twin lifts 15bp from the read's direct primary alignment (lift jitter).
# At the default tolerance (10) the twin is a separate group and competes, so the
# read is driven to MAPQ 0.  Widening the tolerance to 20 groups them and the
# read recovers a confident MAPQ.  This locks BOTH the knob plumbing and the
# documented semantics (the default is conservative; jittery .alt sets can widen
# it without merging genuine paralogs, which the placement guard still separates).
#
# Usage: test/altlg/test-lifttol.sh [<minibwa-dir>]
set -eu

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MINIBWA="$MDIR/minibwa"
MK="$MDIR/test/altlg/mkfixture-lifttol.sh"

TMPD=$(mktemp -d /tmp/altlg-lifttol.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT
fail() { echo "FAIL: $1"; exit 1; }
ok()   { echo "  ok: $1"; }

# primary MAPQ of <qname>: the record with neither 0x100 nor 0x800 set.
primary_mapq() {
    mawk -v q="$2" '$1==q { f=int($2);
        if (int(f/256)%2==0 && int(f/2048)%2==0) { print $5; exit } }' "$1"
}

echo "[test-lifttol] building fixture ..."
/bin/sh "$MK" "$TMPD" 2>/dev/null
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null

run_mapq() { # <lift_tol>
    "$MINIBWA" mem --outn=50 --alt-lift-tol "$1" "$TMPD/ref.fa" "$TMPD/reads.fq" 2>/dev/null \
        | mawk '$1 !~ /^@/' > "$TMPD/t$1.sam"
    [ -s "$TMPD/t$1.sam" ] || fail "tol=$1: no alignments emitted"
    primary_mapq "$TMPD/t$1.sam" r-near
}

echo "== default tolerance (10): 15bp jitter under-merged => MAPQ 0 =="
q10=$(run_mapq 10)
[ -n "$q10" ] || fail "no primary record for r-near at tol=10"
echo "  r-near primary MAPQ (tol=10) = $q10"
[ "$q10" = "0" ] || fail "expected MAPQ 0 at tol=10 (twin a separate group), got $q10"
ok "tol=10: r-near MAPQ=0 (near-twin not merged; conservative default)"

echo "== widened tolerance (20): 15bp jitter absorbed => twin groups, MAPQ>0 =="
q20=$(run_mapq 20)
[ -n "$q20" ] || fail "no primary record for r-near at tol=20"
echo "  r-near primary MAPQ (tol=20) = $q20"
[ "$q20" -gt 0 ] || fail "expected MAPQ>0 at tol=20 (twin grouped), got $q20"
ok "tol=20: r-near MAPQ=$q20 > 0 (near-twin grouped; recovered)"

echo "[test-lifttol] PASS"
exit 0
