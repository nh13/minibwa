#!/bin/sh
# Integration test for Hook C UNLIFTABLE-demotion path (pe.c mb_distinct_lifted_group).
#
# MECHANISM:
#   mb_distinct_lifted_group() returns 0 when EITHER hit is unliftable (!pa.liftable
#   || !pb.liftable), so it falls through to NORMAL chimeric demotion.  Only two
#   LIFTABLE hits that are genuinely distinct primary groups get the "continue" skip.
#   An unliftable is_alt hit is NOT a protected distinct group -- shielding it would
#   preserve a spurious competitor and wrongly affect MAPQ.
#
# FIXTURE (mkfixture-hookc.sh):
#   chrHC        4000bp primary
#   chrHC_altU   200bp ALT contig:
#     bases [0,100)   = chrHC[500,600) -- covered by .alt block (liftable)
#     bases [100,200) = chrHC[700,800) -- NO .alt block -> LIFT HOLE (unliftable)
#   .alt: chrHC_altU -> chrHC POS 501, CIGAR 100M  (covers only alt[0,100))
#
#   R1 = chrHC[300,420) fwd (unique anchor)
#   R2 = revcomp(chrHC[700,800)):
#     - maps to chrHC@700        (primary, liftable=1)
#     - maps to chrHC_altU@100   (ALT, footprint [100,200) in HOLE -> liftable=0)
#
# WHAT IS PROVEN:
#   (1) ex-place-check confirms the chrHC_altU hit is unliftable (liftable=0) --
#       fails loudly if liftable=1 (wrong fixture; test proves nothing).
#   (2) The unliftable ALT hit is demoted to secondary (0x100) -- NOT shielded.
#       Would FAIL if mb_distinct_lifted_group returned 1 for unliftable hits,
#       because the demotion "continue" would fire and the hit would stay as a
#       supplementary (0x800) separate-group representative.
#   (3) The primary chrHC R2 hit is a proper-pair primary with MAPQ > 0.
#   (4) chrM PE baseline: no .alt -> byte-identical across runs.
#
# Usage: test/altlg/test-hookc.sh [<minibwa-dir>]
set -eu

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MINIBWA="$MDIR/minibwa"

# Shared no-.alt baseline check (see lib-baseline.sh).
. "$(dirname "$0")/lib-baseline.sh"
MKFIXTURE="$MDIR/test/altlg/mkfixture-hookc.sh"
EX_PLACE="$MDIR/api-test/ex-place-check"

TMPD=$(mktemp -d /tmp/altlg-hookc.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

fail() { echo "FAIL: $1"; exit 1; }
ok()   { echo "  ok: $1"; }

has_bit() { mawk -v f="$1" -v b="$2" 'BEGIN{ f=int(f); b=int(b);
    printf "%d\n", (int(f/b)%2==1)?1:0; }'; }
flag_of() { mawk -v q="$2" -v c="$3" -v p="${4:-}" \
    '$1==q && $3==c && (p=="" || $4==p){print $2; exit}' "$1"; }
mapq_of() { mawk -v q="$2" -v c="$3" -v p="${4:-}" \
    '$1==q && $3==c && (p=="" || $4==p){print $5; exit}' "$1"; }

# =========================================================================
echo "[test-hookc] building fixture ..."
/bin/sh "$MKFIXTURE" "$TMPD/fix" 2>/dev/null

echo "[test-hookc] indexing ..."
"$MINIBWA" index "$TMPD/fix/ref.fa" 2>/dev/null

# =========================================================================
# (1) PROVE the ALT hit is unliftable via ex-place-check.
#     r2_probe.fq contains only R2 so ex-place-check can call mb_map on it.
echo "[test-hookc] running ex-place-check to prove unliftability ..."
"$EX_PLACE" "$TMPD/fix/ref.fa" "$TMPD/fix/r2_probe.fq" all 2>/dev/null \
    > "$TMPD/place.txt"
echo "----- ex-place-check output (R2 probe) -----"
cat "$TMPD/place.txt"
echo "---------------------------------------------"

# Extract liftable flag for the chrHC_altU hit (is_alt=1 line).
alt_liftable=$(mawk '$2=="chrHC_altU"{
    for(i=1;i<=NF;i++){
        if(index($i,"liftable=")==1){print substr($i,10); exit}
    }}' "$TMPD/place.txt")
[ -n "$alt_liftable" ] || fail "ex-place-check: no chrHC_altU hit found for R2 probe (fixture broken)"
[ "$alt_liftable" = "0" ] || \
    fail "ex-place-check: chrHC_altU liftable=$alt_liftable (expected 0 -- fix is in HOLE); fixture does not exercise the intended path"
ok "chrHC_altU ALT hit: liftable=0 confirmed (footprint in lift HOLE)"

# Also confirm the primary chrHC hit IS liftable (sanity check).
pri_liftable=$(mawk '$2=="chrHC"{
    for(i=1;i<=NF;i++){
        if(index($i,"liftable=")==1){print substr($i,10); exit}
    }}' "$TMPD/place.txt")
[ "$pri_liftable" = "1" ] || \
    fail "ex-place-check: chrHC primary hit liftable=$pri_liftable (expected 1)"
ok "chrHC primary hit: liftable=1 confirmed"

# =========================================================================
# (2)+(3) Run PE mapping and assert on flags/MAPQ.
echo "[test-hookc] running PE minibwa mem ..."
"$MINIBWA" mem --outn=50 "$TMPD/fix/ref.fa" \
    "$TMPD/fix/reads_1.fq" "$TMPD/fix/reads_2.fq" 2>/dev/null \
    | mawk '$1 !~ /^@/' > "$TMPD/pe.sam"
[ -s "$TMPD/pe.sam" ] || fail "PE mapping produced no alignments"

echo "----- r-hookc SAM records -----"
mawk '$1=="r-hookc"{printf "%s flag=%d %s pos=%s mapq=%s\n",$1,$2,$3,$4,$5}' \
    "$TMPD/pe.sam"
echo "-------------------------------"

# --- (2a) The unliftable ALT hit (chrHC_altU, last mate) must be secondary (0x100). ---
# It is demoted because mb_distinct_lifted_group returns 0 for the unliftable hit,
# so the chimeric demotion fires normally.
alt_flag=$(flag_of "$TMPD/pe.sam" r-hookc chrHC_altU 101)
[ -n "$alt_flag" ] || fail "no r-hookc record on chrHC_altU pos=101 (ALT hit absent)"
[ "$(has_bit "$alt_flag" 256)" = "1" ] || \
    fail "chrHC_altU pos=101 flag=$alt_flag NOT secondary (0x100): unliftable ALT hit was SHIELDED by Hook C (mb_distinct_lifted_group wrongly returned 1 for unliftable hit)"
ok "chrHC_altU pos=101 (unliftable ALT hit) is secondary 0x100 (flag=$alt_flag): demoted normally, not shielded"

# --- (2b) The unliftable ALT hit must NOT be supplementary (0x800). ---
# If it were 0x800 it would mean it was kept as a separate-group representative --
# that is the "shielded" failure mode (Hook C wrongly fired the continue).
[ "$(has_bit "$alt_flag" 2048)" = "0" ] || \
    fail "chrHC_altU pos=101 flag=$alt_flag is supplementary (0x800): unliftable ALT hit was treated as a distinct group (wrongly shielded)"
ok "chrHC_altU pos=101 (unliftable ALT hit) is NOT supplementary: correctly NOT treated as a distinct group"

# --- (3a) R2 primary (chrHC pos=701) must be a proper-pair primary with MAPQ > 0. ---
r2_flag=$(flag_of "$TMPD/pe.sam" r-hookc chrHC 701)
[ -n "$r2_flag" ] || fail "no r-hookc record on chrHC pos=701 (R2 primary absent)"
[ "$(has_bit "$r2_flag" 256)" = "0" ] || \
    fail "R2 primary chrHC pos=701 flag=$r2_flag is secondary -- pairing failed"
[ "$(has_bit "$r2_flag" 2)" = "1" ] || \
    fail "R2 primary chrHC pos=701 flag=$r2_flag is NOT a proper pair"
r2_mapq=$(mapq_of "$TMPD/pe.sam" r-hookc chrHC 701)
[ "$r2_mapq" -gt 0 ] || \
    fail "R2 primary chrHC pos=701 MAPQ=$r2_mapq (expected >0): unliftable competitor inflated or collapsed MAPQ"
ok "R2 primary: chrHC pos=701 proper-pair primary, MAPQ=$r2_mapq > 0 (flag=$r2_flag)"

# --- (3b) R1 primary (chrHC pos=301) must be a proper-pair primary with MAPQ > 0. ---
r1_flag=$(flag_of "$TMPD/pe.sam" r-hookc chrHC 301)
[ -n "$r1_flag" ] || fail "no r-hookc record on chrHC pos=301 (R1 primary absent)"
[ "$(has_bit "$r1_flag" 256)" = "0" ] || \
    fail "R1 primary chrHC pos=301 flag=$r1_flag is secondary"
[ "$(has_bit "$r1_flag" 2)" = "1" ] || \
    fail "R1 primary chrHC pos=301 flag=$r1_flag is NOT a proper pair"
r1_mapq=$(mapq_of "$TMPD/pe.sam" r-hookc chrHC 301)
[ "$r1_mapq" -gt 0 ] || \
    fail "R1 primary chrHC pos=301 MAPQ=$r1_mapq (expected >0)"
ok "R1 primary: chrHC pos=301 proper-pair primary, MAPQ=$r1_mapq > 0 (flag=$r1_flag)"

# =========================================================================
# (4) Baseline: chrM PE with no .alt -> byte-identical across runs.
echo "[test-hookc] chrM baseline (no .alt): byte-identical to stock ..."
chrm_baseline "(BASELINE) chrM" pe --outn=5

echo "[test-hookc] PASS"
exit 0
