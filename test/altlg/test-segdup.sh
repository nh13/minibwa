#!/bin/sh
# Integration test for Task 7: ALT-seed -> primary anchor projection.
#
# Proves the segdup-recovery mechanism end-to-end on a single binary by toggling
# the projection with the MB_NO_ALT_PROJECT test seam:
#
#   RED  (MB_NO_ALT_PROJECT=1): a segduplicated primary locus whose CORE seed is
#        subsampled out by `-c` produces NO primary candidate -> there is NO chrP
#        record at the ALT-lifted locus (POS 201).
#   GREEN (default): the surviving, unique ALT seed lifts to that locus and
#        injects a primary anchor -> a chrP record at POS 201 appears.
#
# Also asserts the chrM baseline is byte-identical with/without the seam (no .alt
# => no ALT anchors => projection is a strict no-op), for SE and PE.
#
# Usage: test/altlg/test-segdup.sh [<minibwa-dir>]
set -eu

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MINIBWA="$MDIR/minibwa"
MK_SEG="$MDIR/test/altlg/mkfixture-segdup.sh"

C=4         # max_occ: < N+1 so the CORE seed interval is subsampled
N=20        # primary copies of CORE

TMPD=$(mktemp -d /tmp/altlg-segdup.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

fail() { echo "FAIL: $1"; exit 1; }
ok()   { echo "  ok: $1"; }

# chrP_pos_list <sam>  -> sorted space-separated list of chrP POS values
chrP_pos_list() { mawk '$1!~/^@/ && $3=="chrP"{print $4}' "$1" | sort -n | tr '\n' ' '; }
# has_pos <sam> <pos>  -> "1" if any chrP record has that POS, else "0"
has_pos() { mawk -v p="$2" 'BEGIN{f=0} $1!~/^@/ && $3=="chrP" && $4==p{f=1} END{print f}' "$1"; }

echo "[test-segdup] building segdup fixture (N=$N copies) ..."
/bin/sh "$MK_SEG" "$TMPD" "$N" 2>/dev/null
. "$TMPD/meta.txt"   # sets LIFT_POS1
[ -n "${LIFT_POS1:-}" ] || fail "fixture did not export LIFT_POS1"
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null

echo "[test-segdup] LIFT locus is chrP POS=$LIFT_POS1; max_occ (-c) = $C"

# --- RED: projection disabled ---
echo "== RED: MB_NO_ALT_PROJECT=1 (no projection) =="
"$MINIBWA" map --dbg-no-alt-proj -c "$C" --outn=999 "$TMPD/ref.fa" "$TMPD/reads.fq" \
    2>/dev/null > "$TMPD/red.sam"
mawk '$1 !~ /^@/' "$TMPD/red.sam" > "$TMPD/red.body.sam"
[ -s "$TMPD/red.body.sam" ] || fail "RED: no alignments emitted"
echo "  chrP POS (RED):   $(chrP_pos_list "$TMPD/red.sam")"
red_hit=$(has_pos "$TMPD/red.sam" "$LIFT_POS1")
[ "$red_hit" = "0" ] || fail "RED: a chrP record exists at LIFT POS=$LIFT_POS1 without projection; fixture is not forcing the segdup-subsample drop (try a smaller -c or larger N)"
ok "RED: no chrP candidate at LIFT POS=$LIFT_POS1 (segdup seed subsampled out)"

# --- GREEN: projection enabled (default) ---
echo "== GREEN: projection enabled (default) =="
"$MINIBWA" map -c "$C" --outn=999 "$TMPD/ref.fa" "$TMPD/reads.fq" \
    2>/dev/null > "$TMPD/green.sam"
echo "  chrP POS (GREEN): $(chrP_pos_list "$TMPD/green.sam")"
green_hit=$(has_pos "$TMPD/green.sam" "$LIFT_POS1")
[ "$green_hit" = "1" ] || fail "GREEN: projection did not produce a chrP candidate at LIFT POS=$LIFT_POS1"
ok "GREEN: chrP candidate recovered at LIFT POS=$LIFT_POS1 (ALT seed projected to primary)"

# The ALT twin itself must still be present (sanity: the ALT seed survived).
alt_present=$(mawk '$1!~/^@/ && $3=="chrP_altS"{print "1"; exit}' "$TMPD/green.sam")
[ "${alt_present:-0}" = "1" ] || fail "GREEN: chrP_altS ALT hit missing (ALT seed did not survive)"
ok "GREEN: chrP_altS ALT hit present (the projected anchor's source seed)"

# Determinism: GREEN chrP POS list is stable across runs.
"$MINIBWA" map -c "$C" --outn=999 "$TMPD/ref.fa" "$TMPD/reads.fq" 2>/dev/null > "$TMPD/green2.sam"
[ "$(chrP_pos_list "$TMPD/green.sam")" = "$(chrP_pos_list "$TMPD/green2.sam")" ] \
    || fail "GREEN: chrP POS list not deterministic across runs"
ok "GREEN: chrP POS list deterministic across runs"

# =========================================================================
# --- baseline: chrM (no .alt) -- projection is a strict no-op; SE + PE ---
echo "[test-segdup] chrM baseline (no .alt; projection inert) ..."
CHRM_FA="$MDIR/test/chrM-human.fa.gz"
CHRM_R1="$MDIR/test/chrM-read_1.fa.gz"
CHRM_R2="$MDIR/test/chrM-read_2.fa.gz"
if [ -f "$CHRM_FA" ] && [ -f "$CHRM_R1" ]; then
    # Copy and index inside TMPD rather than beside the shared fixture (same reason as
    # lib-baseline.sh): keeps the check hermetic, avoids racing another script indexing
    # the same chrM, and leaves no untracked .mbw/.l2b in test/.
    cp "$CHRM_FA" "$TMPD/chrM-ref.fa.gz"; CHRM_FA="$TMPD/chrM-ref.fa.gz"
    "$MINIBWA" index "$CHRM_FA" 2>/dev/null
    # SE: default vs seam-off must be byte-identical (no ALT anchors either way).
    "$MINIBWA" map --outn=5 "$CHRM_FA" "$CHRM_R1" 2>/dev/null > "$TMPD/chrM-se-on.sam"
    "$MINIBWA" map --dbg-no-alt-proj --outn=5 "$CHRM_FA" "$CHRM_R1" 2>/dev/null > "$TMPD/chrM-se-off.sam"
    [ -s "$TMPD/chrM-se-on.sam" ] || fail "chrM SE baseline: empty output"
    grep -v '^@' "$TMPD/chrM-se-on.sam" > "$TMPD/chrM-se-on.body" || true
    grep -v '^@' "$TMPD/chrM-se-off.sam" > "$TMPD/chrM-se-off.body" || true
    cmp -s "$TMPD/chrM-se-on.body" "$TMPD/chrM-se-off.body" \
        || fail "chrM SE baseline: projection changed alignments without .alt"
    ok "chrM SE baseline: identical alignments with/without projection (no .alt => inert)"
    if [ -f "$CHRM_R2" ]; then
        "$MINIBWA" map --outn=5 "$CHRM_FA" "$CHRM_R1" "$CHRM_R2" 2>/dev/null > "$TMPD/chrM-pe-on.sam"
        "$MINIBWA" map --dbg-no-alt-proj --outn=5 "$CHRM_FA" "$CHRM_R1" "$CHRM_R2" 2>/dev/null > "$TMPD/chrM-pe-off.sam"
        grep -v '^@' "$TMPD/chrM-pe-on.sam" > "$TMPD/chrM-pe-on.body" || true
        grep -v '^@' "$TMPD/chrM-pe-off.sam" > "$TMPD/chrM-pe-off.body" || true
        cmp -s "$TMPD/chrM-pe-on.body" "$TMPD/chrM-pe-off.body" \
            || fail "chrM PE baseline: projection changed alignments without .alt"
        ok "chrM PE baseline: identical alignments with/without projection (no .alt => inert)"
    else
        echo "  skip: chrM R2 not found ($CHRM_R2)"
    fi
else
    echo "  skip: chrM baseline files not found ($CHRM_FA)"
fi

echo "[test-segdup] PASS"
exit 0
