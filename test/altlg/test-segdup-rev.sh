#!/bin/sh
# Regression test for Task 7 fix round 1, Fix 1 + Fix 2: REVERSE-strand ALT-seed
# -> primary projection.
#
# The original forward-only segdup fixture never exercised a reverse .alt block,
# so two reverse coordinate bugs slipped through:
#   Fix 1: the ALT contig-local forward span recovery was forward-only.
#   Fix 2: the projected primary tpos was re-folded for the reverse strand, but
#          process_batch keeps the position component in the FORWARD contig frame
#          for both strands -> the lifted POS was the mirror image (wrong end of
#          the contig) instead of the true locus.
#
# This test uses a REVERSE .alt block (FLAG 0x10) and asserts, via the
# MB_PROJ_TRACE probe seam, that the surviving full-read ALT seed projects to the
# CORRECT primary coordinate on the REVERSE strand:  chrP POS=LIFT_POS1, strand -.
#
# Why a probe instead of a pure SAM assertion: under a reverse (RC) alignment the
# shared CORE collapses to one of many identical-scoring paralog loci, so the
# projected reverse anchor is masked at SAM level by paralog selection (the same
# happens with projection OFF -- the LIFT locus simply isn't one of the few
# secondaries minibwa surfaces).  The probe asserts the load-bearing coordinate
# the production projection actually computes, which is exactly what Fix 1+2 set.
# RED (projection off) emits NO probe line; GREEN (default) emits the correct one.
#
# Usage: test/altlg/test-segdup-rev.sh [<minibwa-dir>]
set -eu

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MINIBWA="$MDIR/minibwa"
MK="$MDIR/test/altlg/mkfixture-segdup-rev.sh"

C=4
N=20

TMPD=$(mktemp -d /tmp/altlg-segdup-rev.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

fail() { echo "FAIL: $1"; exit 1; }
ok()   { echo "  ok: $1"; }

# proj_line <stderr-trace> <strand>  -> the MB_PROJ POS for chrP on that strand,
# for the FULL-read seed (len=150); empty if none.
proj_pos() { mawk -v s="$2" '$1=="MB_PROJ" && $2=="chrP" && $4==s && $5=="len=150"{print $3}' "$1"; }

echo "[test-segdup-rev] building reverse segdup fixture (N=$N) ..."
/bin/sh "$MK" "$TMPD" "$N" 2>/dev/null
. "$TMPD/meta.txt"
[ -n "${LIFT_POS1:-}" ] || fail "fixture did not export LIFT_POS1"
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null
echo "[test-segdup-rev] expected projected locus: chrP POS=$LIFT_POS1 strand=- (full-read seed)"

# --- RED: projection disabled -> no projected anchor at all ---
echo "== RED: MB_NO_ALT_PROJECT=1 (no projection) =="
MB_NO_ALT_PROJECT=1 MB_PROJ_TRACE=1 "$MINIBWA" mem -c "$C" --outn=999 \
    "$TMPD/ref.fa" "$TMPD/reads.fq" 2>"$TMPD/red.err" >/dev/null
red_pos=$(proj_pos "$TMPD/red.err" -)
[ -z "$red_pos" ] || fail "RED: projection emitted a trace with seam off (got POS=$red_pos)"
ok "RED: no projected anchor with projection disabled"

# --- GREEN: projection enabled (default) ---
echo "== GREEN: projection enabled (default) =="
MB_PROJ_TRACE=1 "$MINIBWA" mem -c "$C" --outn=999 \
    "$TMPD/ref.fa" "$TMPD/reads.fq" 2>"$TMPD/green.err" >/dev/null
echo "  MB_PROJ traces (GREEN):"
grep '^MB_PROJ' "$TMPD/green.err" | sed 's/^/    /' || true
green_pos=$(proj_pos "$TMPD/green.err" -)
[ -n "$green_pos" ] || fail "GREEN: no reverse-strand projected anchor for the full-read seed"
echo "  expected POS=$LIFT_POS1 strand=-, actual POS=$green_pos strand=-"
[ "$green_pos" = "$LIFT_POS1" ] \
    || fail "GREEN: reverse projection landed at POS=$green_pos, expected $LIFT_POS1 (Fix 1/2 coordinate bug)"
ok "GREEN: reverse ALT seed projected to the CORRECT primary locus (POS=$LIFT_POS1, strand -)"

# Determinism.
MB_PROJ_TRACE=1 "$MINIBWA" mem -c "$C" --outn=999 \
    "$TMPD/ref.fa" "$TMPD/reads.fq" 2>"$TMPD/green2.err" >/dev/null
[ "$(proj_pos "$TMPD/green2.err" -)" = "$green_pos" ] \
    || fail "GREEN: reverse projected POS not deterministic across runs"
ok "GREEN: reverse projected POS deterministic across runs"

echo "[test-segdup-rev] PASS"
exit 0
