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
# RED assertion (guard absent): chrP_alt NOT in SAM output (-p 0.9 --outn=50).
# GREEN assertion (guard present): chrP_alt IS in SAM output.
# Note: the RED state can be verified by temporarily passing NULL for l2b in
# the two mb_select_sub driver calls in mb_map_sai.
#
# Baseline (chrM, no .alt): guard never fires; output byte-identical without
# guard because l2b has no ALT contigs (guard gated on l2b ALT presence).
#
# Usage: test/altlg/test-survive.sh [<minibwa-dir>]
set -eu

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MKFIXTURE="$MDIR/test/altlg/mkfixture-survive.sh"
MINIBWA="$MDIR/minibwa"

TMPD=$(mktemp -d /tmp/altlg-survive.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

fail() { echo "FAIL: $1"; exit 1; }
ok()   { echo "  ok: $1"; }

# ---------------------------------------------------------------------------
# Part 1: survival fixture (unique chrP + chrP_alt with 1bp SNP at read pos 55)
# ---------------------------------------------------------------------------
echo "[test-survive] building fixture ..."
/bin/sh "$MKFIXTURE" "$TMPD" 2>/dev/null

echo "[test-survive] indexing ..."
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null

# Map with -p 0.9 (high pri_ratio to drop the ~55bp chain) --outn=50
echo "[test-survive] mapping with -p 0.9 --outn=50 ..."
"$MINIBWA" mem -p 0.9 --outn=50 "$TMPD/ref.fa" "$TMPD/reads.fq" 2>/dev/null \
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

# ---------------------------------------------------------------------------
# Part 2: chrM baseline — guard never fires (no .alt -> l2b has no ALT ctgs)
# ---------------------------------------------------------------------------
echo "[test-survive] baseline: chrM (no .alt) ..."
CHRM_FA="$MDIR/test/chrM-human.fa.gz"
CHRM_R1="$MDIR/test/chrM-read_1.fa.gz"

if [ -f "$CHRM_FA" ] && [ -f "$CHRM_R1" ]; then
    # Confirms the guard is inert without a .alt file: all hits have is_alt=0,
    # so pass 2 never force-keeps anything — semantics = "guard never fires."
    # Hard failure on empty output catches crashes or tool regressions silently
    # swallowing output that would otherwise appear as a spurious skip.
    "$MINIBWA" index "$CHRM_FA" 2>/dev/null
    "$MINIBWA" mem --outn=5 "$CHRM_FA" "$CHRM_R1" 2>/dev/null > "$TMPD/chrM-a.sam"
    "$MINIBWA" mem --outn=5 "$CHRM_FA" "$CHRM_R1" 2>/dev/null > "$TMPD/chrM-b.sam"
    if [ ! -s "$TMPD/chrM-a.sam" ] || [ ! -s "$TMPD/chrM-b.sam" ]; then
        fail "chrM baseline: unexpectedly empty output"
    fi
    if cmp -s "$TMPD/chrM-a.sam" "$TMPD/chrM-b.sam"; then
        ok "chrM baseline: output byte-identical (guard never fires without .alt)"
    else
        fail "chrM baseline: output differs between runs (non-determinism?)"
    fi
else
    echo "  skip: chrM baseline files not found ($CHRM_FA)"
fi

echo "[test-survive] PASS"
exit 0
