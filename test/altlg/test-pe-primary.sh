#!/bin/sh
# Integration test for the opt-in PE-pair-primary selection (--pe-pair-primary).
#
# r-para's R1 maps equally to two paralog copies (copyA=chrP[500,650),
# copyB=chrP[2500,2650)); its mate R2 is unique near copyA, so the PE pairing
# chooses copyA.  Assertions:
#   - WITH --pe-pair-primary: R1 SAM primary is at copyA (chrP pos in [490,660)),
#     i.e. the mate-consistent copy -- NOT copyB.
#   - default (no flag): chrM PE output is byte-identical across runs (the opt-in
#     selection is off, so plain `minibwa mem` is unchanged).
#
# Usage: test/altlg/test-pe-primary.sh [<minibwa-dir>]
set -eu
MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MINIBWA="$MDIR/minibwa"
MK="$MDIR/test/altlg/mkfixture-pe-primary.sh"
TMPD=$(mktemp -d /tmp/altlg-pe-primary.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT
fail() { echo "FAIL: $1"; exit 1; }
ok()   { echo "  ok: $1"; }

# R1 primary record (no 0x100/0x800) position for qname/1
r1_primary_pos() {
    mawk -v q="$2" '$1==q { f=int($2);
        if (int(f/64)%2==1 && int(f/256)%2==0 && int(f/2048)%2==0) { print $4; exit } }' "$1"
}
r1_primary_rname() {
    mawk -v q="$2" '$1==q { f=int($2);
        if (int(f/64)%2==1 && int(f/256)%2==0 && int(f/2048)%2==0) { print $3; exit } }' "$1"
}

echo "[test-pe-primary] building fixture ..."
/bin/sh "$MK" "$TMPD" 2>/dev/null
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null

in_copyA() { [ "$1" -ge 490 ] && [ "$1" -lt 660 ]; }   # copyA window (~chrP:501)
in_copyB() { [ "$1" -ge 2490 ] && [ "$1" -lt 2660 ]; } # copyB window (~chrP:2501)

echo "[test-pe-primary] mapping DEFAULT (no flag) ..."
"$MINIBWA" mem --outn=50 "$TMPD/ref.fa" "$TMPD/reads_1.fq" "$TMPD/reads_2.fq" \
    2>/dev/null | mawk '$1 !~ /^@/' > "$TMPD/off.sam"
[ -s "$TMPD/off.sam" ] || fail "no alignments emitted (default)"
echo "----- r-para records (default) -----"
mawk '$1=="r-para"{printf "  flag=%d %s:%s mapq=%s\n",$2,$3,$4,$5}' "$TMPD/off.sam"
off_pos=$(r1_primary_pos "$TMPD/off.sam" r-para)
[ -n "$off_pos" ] || fail "no R1 primary record for r-para (default)"
# The bug: default per-read order emits the higher-scoring (exact) copyB, which
# is NOT mate-consistent.  Locks that the flag is actually doing something.
if in_copyB "$off_pos"; then
    ok "default emits copyB (chrP:$off_pos) -- the mate-INconsistent copy (bug present without flag)"
else
    fail "default R1 primary at chrP:$off_pos; fixture expected copyB (~2501) to win by score without the flag"
fi

echo "[test-pe-primary] mapping WITH --pe-pair-primary ..."
"$MINIBWA" mem --outn=50 --pe-pair-primary "$TMPD/ref.fa" "$TMPD/reads_1.fq" "$TMPD/reads_2.fq" \
    2>/dev/null | mawk '$1 !~ /^@/' > "$TMPD/on.sam"
[ -s "$TMPD/on.sam" ] || fail "no alignments emitted"
echo "----- r-para records (--pe-pair-primary) -----"
mawk '$1=="r-para"{printf "  flag=%d %s:%s mapq=%s\n",$2,$3,$4,$5}' "$TMPD/on.sam"
rn=$(r1_primary_rname "$TMPD/on.sam" r-para)
pos=$(r1_primary_pos "$TMPD/on.sam" r-para)
[ -n "$pos" ] || fail "no R1 primary record for r-para"
[ "$rn" = "chrP" ] || fail "r-para/1 primary on '$rn'; expected chrP"
if in_copyA "$pos"; then
    ok "r-para/1 primary at copyA (chrP:$pos) -- mate-consistent copy chosen by --pe-pair-primary"
else
    fail "r-para/1 primary at chrP:$pos; expected copyA in [490,660) (mate-consistent), not copyB"
fi

echo "[test-pe-primary] default (no flag) chrM PE byte-identical ..."
CHRM_FA="$MDIR/test/chrM-human.fa.gz"; CHRM_R1="$MDIR/test/chrM-read_1.fa.gz"; CHRM_R2="$MDIR/test/chrM-read_2.fa.gz"
if [ -f "$CHRM_FA" ] && [ -f "$CHRM_R1" ] && [ -f "$CHRM_R2" ]; then
    "$MINIBWA" index "$CHRM_FA" 2>/dev/null
    "$MINIBWA" mem --outn=5 "$CHRM_FA" "$CHRM_R1" "$CHRM_R2" 2>/dev/null > "$TMPD/m-a.sam"
    "$MINIBWA" mem --outn=5 "$CHRM_FA" "$CHRM_R1" "$CHRM_R2" 2>/dev/null > "$TMPD/m-b.sam"
    cmp -s "$TMPD/m-a.sam" "$TMPD/m-b.sam" \
        && ok "chrM PE byte-identical across default runs (opt-in off)" \
        || fail "chrM PE differs between default runs"
else
    echo "  skip: chrM PE baseline files not found"
fi

echo "[test-pe-primary] PASS"
exit 0
