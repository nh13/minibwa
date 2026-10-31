#!/bin/sh
# Integration test for the ALT-ALTERNATE fold-in (mb_reconcile_alt step 2b).
#
# An ALT twin that the per-base lift cannot co-locate with its primary hit --
# because an .alt-internal structural indel DISPLACES it beyond MB_LIFT_TOL
# (r-disp) or drops it into an insertion HOLE (r-hole) -- must still be folded
# into the primary's group, recognized as an alternate placement of that locus
# via the .alt correspondence.  After the fold:
#   - the SAM primary is the non-ALT chrP hit (read NOT placed on the ALT contig);
#   - the primary MAPQ is > 0 (the ALT twin no longer dilutes it);
#   - the ALT twin is secondary (0x100) or dropped, and at the grouping level is
#     a subordinate (parent != its own id).
#
# Baseline: chrM (no .alt) => mb_any_alt gate off => reconcile inert => output
# byte-identical across runs.
#
# Usage: test/altlg/test-altalt.sh [<minibwa-dir>]
set -eu
. "$(dirname "$0")/lib.sh"

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MINIBWA="$MDIR/minibwa"

# Shared no-.alt baseline check (see lib-baseline.sh).
. "$(dirname "$0")/lib-baseline.sh"
EX_GROUP="$MDIR/api-test/ex-group-check"
MK="$MDIR/test/altlg/mkfixture-altalt.sh"

TMPD=$(mktemp -d /tmp/altlg-altalt.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT


# primary record (no 0x100, no 0x800) fields for a qname
primary_rname() { mawk -v q="$2" '$1==q{f=int($2); if(int(f/256)%2==0 && int(f/2048)%2==0){print $3; exit}}' "$1"; }

# group-check field for the ALT contig hit of a qname (is_alt=1 row on <ctg>)
grp_field() {
    mawk -v q="$2" -v c="$3" -v k="$4" '
        $1==q { ctg=""; want=""; for(i=1;i<=NF;i++){ n=index($i,"=");
            if(n>0){ key=substr($i,1,n-1); val=substr($i,n+1);
                     if(key=="ctg") ctg=val; if(key==k) want=val; } }
          if(ctg==c){ print want; exit } }' "$1"
}

echo "[test-altalt] building fixture ..."
/bin/sh "$MK" "$TMPD" 2>/dev/null
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null

echo "[test-altalt] mapping (--outn=50) ..."
"$MINIBWA" map --outn=50 "$TMPD/ref.fa" "$TMPD/reads.fq" 2>/dev/null > "$TMPD/full.sam"
mawk '$1 !~ /^@/' "$TMPD/full.sam" > "$TMPD/aln.sam"
[ -s "$TMPD/aln.sam" ] || fail "no alignments emitted"
echo "----- altalt SAM -----"
mawk '{printf "%s flag=%d %s pos=%s mapq=%s\n",$1,$2,$3,$4,$5}' "$TMPD/aln.sam"
echo "----------------------"

[ -x "$EX_GROUP" ] || fail "ex-group-check probe not built ($EX_GROUP); run 'make -C api-test'"
"$EX_GROUP" "$TMPD/ref.fa" "$TMPD/reads.fq" 2>/dev/null > "$TMPD/grp.txt"
echo "----- altalt grouping -----"
cat "$TMPD/grp.txt"
echo "---------------------------"

check_read() {   # <qname> <alt-ctg>
    q="$1"; altc="$2"
    pr=$(primary_rname "$TMPD/aln.sam" "$q")
    pf=$(primary_flag  "$TMPD/aln.sam" "$q")
    pq=$(primary_mapq  "$TMPD/aln.sam" "$q")
    [ -n "$pr" ] || fail "($q) no primary record"
    [ "$pr" = "chrP" ] || fail "($q) SAM primary on '$pr'; expected chrP (ALT twin stole the primary slot)"
    ok "($q) SAM primary is chrP (flag=$pf)"
    [ "$pq" -gt 0 ] || fail "($q) primary MAPQ=$pq; expected >0 (ALT twin must not dilute it)"
    ok "($q) primary MAPQ=$pq > 0"
    af=$(flag_of "$TMPD/aln.sam" "$q" "$altc")
    if [ -n "$af" ]; then
        s=$(has_bit "$af" 256)
        [ "$s" = "1" ] || fail "($q) $altc present but not secondary (flag=$af)"
        ok "($q) $altc twin is secondary (flag=$af)"
    else
        ok "($q) $altc twin absent (collapsed into primary group)"
    fi
    # grouping level: the ALT hit must be a subordinate (parent != its own id)
    gid=$(grp_field "$TMPD/grp.txt" "$q" "$altc" id)
    gpa=$(grp_field "$TMPD/grp.txt" "$q" "$altc" parent)
    [ -n "$gid" ] && [ -n "$gpa" ] || fail "($q) no $altc hit in grouping dump"
    [ "$gpa" != "$gid" ] || fail "($q) $altc is its own group rep (parent=$gpa == id=$gid); ALT alternate must be folded into the primary group"
    ok "($q) $altc folded as subordinate (parent=$gpa != id=$gid)"
}

echo "== r-disp: deletion-displaced ALT twin folds into chrP =="
check_read r-disp chrP_altF
echo "== r-hole: insertion-hole (unliftable) ALT twin folds into chrP =="
check_read r-hole chrP_altG

echo "[test-altalt] chrM baseline (no .alt): byte-identical to stock ..."
chrm_baseline "(BASELINE) chrM" se --outn=5

echo "[test-altalt] PASS"
exit 0
