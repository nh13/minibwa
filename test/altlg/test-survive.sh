#!/bin/sh
# Integration test for Task 3 (survival guard in mb_select_sub).
#
# Fixture: chrP (400bp unique sequence) + chrP_alt (ALT = chrP with 1bp SNP at
# 1-based position 126, putting the SNP at read position 55 of the 100bp read
# chrP[70,170)).
#
# With -p 0.9 (pri_ratio=0.9) the chrP_alt twin chains LOWER (sub-SMEM score
# ~55 vs primary ~100) and fails both the score-ratio and min-diff checks in
# mb_select_sub (pre-DP at :617), so it is DROPPED without the guard.
# With the survival guard the twin is force-kept because its lifted placement
# co-locates with the kept chrP primary (same pri_tid=chrP, same rev, |Δlst|=0).
#
# GREEN assertion (guard present): chrP_alt IS in SAM output.
# RED assertion (guard ablated): with --dbg-no-alt-survive the guard is skipped
# and chrP_alt must be ABSENT.  Asserting both states proves the fixture actually
# exercises the guard's drop path (not that chrP_alt survives for some unrelated
# reason), so the test cannot silently stop stressing the guard.
#
# Baseline (chrM, no .alt): guard never fires; output byte-identical without
# guard because l2b has no ALT contigs (guard gated on l2b ALT presence).
#
# Usage: test/altlg/test-survive.sh [<minibwa-dir>]
set -eu
. "$(dirname "$0")/lib.sh"

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MKFIXTURE="$MDIR/test/altlg/mkfixture-survive.sh"
MINIBWA="$MDIR/minibwa"

# Shared no-.alt baseline check (see lib-baseline.sh).
. "$(dirname "$0")/lib-baseline.sh"

TMPD=$(mktemp -d /tmp/altlg-survive.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT


# ---------------------------------------------------------------------------
# Part 1: survival fixture (unique chrP + chrP_alt with 1bp SNP at read pos 55)
# ---------------------------------------------------------------------------
echo "[test-survive] building fixture ..."
/bin/sh "$MKFIXTURE" "$TMPD" 2>/dev/null

echo "[test-survive] indexing ..."
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null

# Map with -p 0.9 (high pri_ratio to drop the ~55bp chain) --outn=50
echo "[test-survive] mapping with -p 0.9 --outn=50 ..."
"$MINIBWA" map -p 0.9 --outn=50 "$TMPD/ref.fa" "$TMPD/reads.fq" 2>/dev/null \
    > "$TMPD/out.sam"

echo "----- SAM output (-p 0.9 --outn=50) -----"
grep -v "^@" "$TMPD/out.sam"
echo "----- end -----"

# Primary alignment on chrP must always be present
pri_present=$(mawk '$1=="r-survive" && $3=="chrP" && ($2+0 == 0) {found=1} END{print found+0}' \
    "$TMPD/out.sam")
if [ "$pri_present" != "1" ]; then
    fail "primary chrP alignment should always be present (got: $pri_present)"
fi
ok "primary chrP alignment present (RNAME=chrP, FLAG=0)"

# chrP_alt alignment must be PRESENT: the guard force-kept it.
# (Without the guard, this fails: RED state.)
alt_present=$(mawk '$1=="r-survive" && $3=="chrP_alt" {found=1} END{print found+0}' \
    "$TMPD/out.sam")
if [ "$alt_present" != "1" ]; then
    fail "chrP_alt alignment should survive (guard keeps it); got: absent [RED: guard not working]"
fi
ok "chrP_alt alignment survived through guard (GREEN)"

# Also verify the lifted placement of chrP_alt co-locates with the primary:
# both map to chrP at approximately the same position.
pri_pos=$(mawk '$1=="r-survive" && $3=="chrP" && ($2+0==0) {print $4}' "$TMPD/out.sam")
alt_pos=$(mawk '$1=="r-survive" && $3=="chrP_alt" {print $4}' "$TMPD/out.sam" | head -1)
if [ -n "$pri_pos" ] && [ -n "$alt_pos" ]; then
    d=$(( alt_pos - pri_pos ))
    [ "$d" -lt 0 ] && d=$(( -d ))
    if [ "$d" -gt 10 ]; then
        fail "chrP_alt POS=$alt_pos not co-located with chrP POS=$pri_pos (|diff|=$d > 10)"
    fi
    ok "chrP_alt POS=$alt_pos co-locates with chrP POS=$pri_pos (|diff|=$d <= 10)"
fi

# --- RED state: ablate ONLY the survival guard and confirm chrP_alt is dropped.
# This proves the fixture geometry genuinely drives the twin below the
# score-ratio/min-diff thresholds, so the GREEN assertion above is meaningful. ---
echo "[test-survive] mapping with --dbg-no-alt-survive (guard ablated) ..."
"$MINIBWA" map -p 0.9 --outn=50 --dbg-no-alt-survive "$TMPD/ref.fa" "$TMPD/reads.fq" 2>/dev/null \
    > "$TMPD/out.red.sam"
echo "----- SAM output (guard ablated) -----"
grep -v "^@" "$TMPD/out.red.sam"
echo "----- end -----"
red_pri=$(mawk '$1=="r-survive" && $3=="chrP" && ($2+0==0) {found=1} END{print found+0}' "$TMPD/out.red.sam")
[ "$red_pri" = "1" ] || fail "(RED) primary chrP alignment should still be present with the guard ablated"
red_alt=$(mawk '$1=="r-survive" && $3=="chrP_alt" {found=1} END{print found+0}' "$TMPD/out.red.sam")
if [ "$red_alt" != "0" ]; then
    fail "(RED) chrP_alt should be DROPPED with the guard ablated, but it survived; the fixture is not exercising the guard's drop path"
fi
ok "chrP_alt dropped when guard ablated (RED): fixture provably stresses the survival guard"

# ---------------------------------------------------------------------------
# Part 2: chrM baseline — guard never fires (no .alt -> l2b has no ALT ctgs)
# ---------------------------------------------------------------------------
echo "[test-survive] chrM baseline (no .alt): byte-identical to stock ..."
chrm_baseline "(BASELINE) chrM" se --outn=5

echo "[test-survive] PASS"
exit 0
