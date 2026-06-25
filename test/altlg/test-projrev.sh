#!/bin/sh
# Integration test: REVERSE-strand ALT-seed -> primary anchor projection,
# observable end-to-end via MAPQ in the SAM output (no probe required).
#
# WHAT IS TESTED
# --------------
# Proves that mb_anchor_project_alt's reverse-block code path (blk_rev ^ alt_rev,
# Fix 1 + Fix 2 in seed.c) recovers a primary alignment that is SAM-observable
# as MAPQ > 0 on the primary assembly, not just through the MB_PROJ_TRACE
# diagnostic probe used by test-segdup-rev.sh.
#
# MECHANISM (see mkfixture-projrev.sh for the full derivation)
# ------------------------------------------------------------
# The read's 150bp SMEM has SA=2: the ALT contig (forward) and the primary LIFT
# locus (reverse complement).  -c 1 (max_occ=1) subsamples this SA=2 interval
# to 1 hit; the SA-array lexicographic ordering for the chosen sequences places
# the ALT-contig hit first, so the strided sampler takes only that hit.
#
# RED  (MB_NO_ALT_PROJECT=1, projection disabled):
#   Only the ALT-contig anchor exists.  The read aligns to chrP_altS (SAM
#   primary, FLAG=0, MAPQ=60).  No chrP record at LIFT_POS1 is emitted.
#
# GREEN (projection enabled, default):
#   The ALT anchor is lifted through the RC .alt block:
#     blk_rev=1, alt_rev=0 -> folded_rev = blk_rev ^ alt_rev = 1 (reverse).
#   A native reverse-strand primary anchor is injected at LIFT_POS1.  DP at
#   LIFT: read maps perfectly (150M, NM=0, score=300) because the primary
#   stores RC(read) at LIFT_POS1.  Reconciliation groups {chrP_altS, chrP LIFT};
#   scores tie at 300; non-ALT preference promotes chrP -> chrP is SAM primary
#   (FLAG=16, MAPQ=60).  chrP_altS becomes secondary (FLAG=256, MAPQ=0).
#
# Assertions:
#   RED:   chrP POS=LIFT_POS1 absent (primary seed subsampled, no projection).
#   GREEN: chrP POS=LIFT_POS1 present with MAPQ > 0 (reverse projection works).
#   GREEN: chrP POS=LIFT_POS1 is the SAM primary record (FLAG & 0x900 == 0).
#   Also: chrM byte-identical baseline (no .alt => projection is a strict no-op).
#
# Usage: test/altlg/test-projrev.sh [<minibwa-dir>]
set -eu

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MINIBWA="$MDIR/minibwa"
MK="$MDIR/test/altlg/mkfixture-projrev.sh"

# max_occ=1: subsamples the 150bp SMEM (SA=2) to 1 hit; the SA-array ordering
# for the chosen sequence seeds deterministically picks the ALT-contig hit,
# leaving the primary LIFT locus unreachable without projection.
C=1
N=20

TMPD=$(mktemp -d /tmp/altlg-projrev.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

fail() { echo "FAIL: $1"; exit 1; }
ok()   { echo "  ok: $1"; }

# chrP POS lookup helpers (same style as test-segdup.sh).
has_pos()     { mawk -v p="$2" 'BEGIN{f=0} $1!~/^@/ && $3=="chrP" && $4==p{f=1} END{print f}' "$1"; }
mapq_at_pos() { mawk -v p="$2" '$1!~/^@/ && $3=="chrP" && $4==p{print $5; exit}' "$1"; }
flag_at_pos() { mawk -v p="$2" '$1!~/^@/ && $3=="chrP" && $4==p{print $2; exit}' "$1"; }

echo "[test-projrev] building reverse-projection fixture (N=$N copies) ..."
/bin/sh "$MK" "$TMPD" "$N" 2>/dev/null
. "$TMPD/meta.txt"   # sets LIFT_POS1
[ -n "${LIFT_POS1:-}" ] || fail "fixture did not export LIFT_POS1"
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null

echo "[test-projrev] LIFT locus is chrP POS=$LIFT_POS1; max_occ (-c) = $C"
echo "[test-projrev] 150bp SMEM SA=2 (ALT-contig + LIFT-reverse); -c 1 subsamples to ALT hit first"

# =========================================================================
# --- RED: projection disabled ---
echo "== RED: MB_NO_ALT_PROJECT=1 (no projection) =="
MB_NO_ALT_PROJECT=1 "$MINIBWA" mem -c "$C" --outn=999 "$TMPD/ref.fa" "$TMPD/reads.fq" \
    2>/dev/null > "$TMPD/red.sam"

red_hit=$(has_pos "$TMPD/red.sam" "$LIFT_POS1")
[ "$red_hit" = "0" ] || fail "RED: chrP record exists at LIFT POS=$LIFT_POS1 without projection; the primary LIFT locus is reachable by direct seeds (SA ordering or max_occ changed - check sequence seeds and -c value)"
ok "RED: no chrP record at LIFT POS=$LIFT_POS1 (primary seed subsampled out, ALT hit sampled instead)"

# Sanity: the ALT contig IS present (the ALT seed was taken by the sampler).
alt_present=$(mawk '$1!~/^@/ && $3=="chrP_altS"{print "1"; exit}' "$TMPD/red.sam")
[ "${alt_present:-0}" = "1" ] || fail "RED: chrP_altS alignment missing (ALT seed not taken by sampler)"
ok "RED: chrP_altS alignment present (ALT seed was the sampled hit, as expected)"

# =========================================================================
# --- GREEN: projection enabled (default) ---
echo "== GREEN: projection enabled (default) =="
"$MINIBWA" mem -c "$C" --outn=999 "$TMPD/ref.fa" "$TMPD/reads.fq" \
    2>/dev/null > "$TMPD/green.sam"

green_hit=$(has_pos "$TMPD/green.sam" "$LIFT_POS1")
[ "$green_hit" = "1" ] || fail "GREEN: projection did not produce a chrP record at LIFT POS=$LIFT_POS1"
ok "GREEN: chrP record present at LIFT POS=$LIFT_POS1 (reverse ALT seed projected to primary)"

# MAPQ > 0: the projected reverse alignment is the SAM primary and has a
# confident score because the primary stores RC(read) at LIFT_POS1 (perfect
# 150M match; score=300, same as the ALT alignment -> non-ALT preference
# promotes chrP to group representative).
green_mapq=$(mapq_at_pos "$TMPD/green.sam" "$LIFT_POS1")
[ -n "$green_mapq" ] || fail "GREEN: could not read MAPQ from chrP POS=$LIFT_POS1 record"
[ "$green_mapq" -gt 0 ] || fail "GREEN: chrP POS=$LIFT_POS1 has MAPQ=$green_mapq (not > 0); ALT contig may have won group rep (check scores)"
ok "GREEN: chrP POS=$LIFT_POS1 has MAPQ=$green_mapq > 0 (confident primary alignment via reverse projection)"

# The chrP LIFT record must be the SAM PRIMARY (FLAG & 0x900 == 0), not a
# secondary.  non-ALT preference makes it the group representative.
green_flag=$(flag_at_pos "$TMPD/green.sam" "$LIFT_POS1")
primary_mask=$(( green_flag & 0x900 ))
[ "$primary_mask" -eq 0 ] || fail "GREEN: chrP POS=$LIFT_POS1 has FLAG=$green_flag (0x$(printf '%x' $green_flag)); expected SAM primary (FLAG & 0x900 == 0)"
ok "GREEN: chrP POS=$LIFT_POS1 is the SAM primary record (FLAG=$green_flag, FLAG & 0x900 = 0)"

# The reverse strand bit must be set (FLAG & 0x10): the projected anchor is
# reverse (folded_rev = blk_rev ^ alt_rev = 1 ^ 0 = 1).
rev_bit=$(( green_flag & 0x10 ))
[ "$rev_bit" -ne 0 ] || fail "GREEN: chrP POS=$LIFT_POS1 is NOT on the reverse strand (FLAG=$green_flag); expected FLAG&0x10 set (folded_rev=1)"
ok "GREEN: chrP POS=$LIFT_POS1 is on the reverse strand (FLAG&0x10 set, folded_rev=blk_rev^alt_rev=1)"

# The ALT twin chrP_altS must be present as a secondary (subordinate to chrP).
alt_in_green=$(mawk '$1!~/^@/ && $3=="chrP_altS"{print "1"; exit}' "$TMPD/green.sam")
[ "${alt_in_green:-0}" = "1" ] || fail "GREEN: chrP_altS secondary missing (ALT seed should still produce an ALT alignment)"
ok "GREEN: chrP_altS present as secondary subordinate in the liftover group"

# Determinism: chrP POS=LIFT_POS1 MAPQ is stable across runs.
"$MINIBWA" mem -c "$C" --outn=999 "$TMPD/ref.fa" "$TMPD/reads.fq" 2>/dev/null > "$TMPD/green2.sam"
green2_mapq=$(mapq_at_pos "$TMPD/green2.sam" "$LIFT_POS1")
[ "$green2_mapq" = "$green_mapq" ] \
    || fail "GREEN: chrP POS=$LIFT_POS1 MAPQ not deterministic (run1=$green_mapq run2=$green2_mapq)"
ok "GREEN: MAPQ=$green_mapq at LIFT POS=$LIFT_POS1 deterministic across runs"

# =========================================================================
# --- baseline: chrM (no .alt) -- projection is a strict no-op; SE + PE ---
echo "[test-projrev] chrM baseline (no .alt; projection inert) ..."
CHRM_FA="$MDIR/test/chrM-human.fa.gz"
CHRM_R1="$MDIR/test/chrM-read_1.fa.gz"
CHRM_R2="$MDIR/test/chrM-read_2.fa.gz"
if [ -f "$CHRM_FA" ] && [ -f "$CHRM_R1" ]; then
    "$MINIBWA" index "$CHRM_FA" 2>/dev/null
    "$MINIBWA" mem --outn=5 "$CHRM_FA" "$CHRM_R1" 2>/dev/null > "$TMPD/chrM-se-on.sam"
    MB_NO_ALT_PROJECT=1 "$MINIBWA" mem --outn=5 "$CHRM_FA" "$CHRM_R1" 2>/dev/null > "$TMPD/chrM-se-off.sam"
    [ -s "$TMPD/chrM-se-on.sam" ] || fail "chrM SE baseline: empty output"
    cmp -s "$TMPD/chrM-se-on.sam" "$TMPD/chrM-se-off.sam" \
        || fail "chrM SE baseline: projection changed output without .alt"
    ok "chrM SE baseline: byte-identical with/without projection (no .alt => inert)"
    if [ -f "$CHRM_R2" ]; then
        "$MINIBWA" mem --outn=5 "$CHRM_FA" "$CHRM_R1" "$CHRM_R2" 2>/dev/null > "$TMPD/chrM-pe-on.sam"
        MB_NO_ALT_PROJECT=1 "$MINIBWA" mem --outn=5 "$CHRM_FA" "$CHRM_R1" "$CHRM_R2" 2>/dev/null > "$TMPD/chrM-pe-off.sam"
        cmp -s "$TMPD/chrM-pe-on.sam" "$TMPD/chrM-pe-off.sam" \
            || fail "chrM PE baseline: projection changed output without .alt"
        ok "chrM PE baseline: byte-identical with/without projection (no .alt => inert)"
    else
        echo "  skip: chrM R2 not found ($CHRM_R2)"
    fi
else
    echo "  skip: chrM baseline files not found ($CHRM_FA)"
fi

echo "[test-projrev] PASS"
exit 0
