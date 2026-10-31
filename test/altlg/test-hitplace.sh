#!/bin/sh
# Integration test for mb_hit_place (Task 2: per-hit lifted PLACEMENT).
#
# Builds the placement fixture (chrP primary; chrP_alt forward ALT with a 5bp
# insertion @300; chrR_alt reverse-strand ALT = revcomp of chrP[250,450)),
# indexes it, maps the reads, and asserts on the placement printed by
# api-test/ex-place-check.  Each ex-place-check line looks like:
#
#   <qname> <ctg> ts=.. te=.. rev=.. is_alt=.. cg=.. || place: pri=<ctg> \
#       lst=<lifted_st> len=<lifted_en> prev=<rev> liftable=<0|1>
#
# Cases (c) and (d) are the BLOCKER branches: a reverse .alt block (strand fold +
# min/max over lifted outputs) and a footprint end sitting in a hole (walk to the
# first liftable base).  A naive lifted_st = lift(alt_st) FAILS both.
#
# Usage: test/altlg/test-hitplace.sh [<minibwa-dir>]
set -eu
. "$(dirname "$0")/lib.sh"

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MKFIXTURE="$MDIR/test/altlg/mkfixture-place.sh"
MINIBWA="$MDIR/minibwa"
EX_PLACE="$MDIR/api-test/ex-place-check"

TMPD=$(mktemp -d /tmp/altlg-hitplace.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

echo "[test-hitplace] building fixture in $TMPD ..."
/bin/sh "$MKFIXTURE" "$TMPD"

echo "[test-hitplace] indexing with minibwa ..."
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null

echo "[test-hitplace] running ex-place-check (exact, post-DP) ..."
"$EX_PLACE" "$TMPD/ref.fa" "$TMPD/reads.fq" all 2>/dev/null > "$TMPD/place.txt"
echo "[test-hitplace] running ex-place-check (coarse, h->p forced NULL) ..."
"$EX_PLACE" "$TMPD/ref.fa" "$TMPD/reads.fq" coarse 2>/dev/null > "$TMPD/place-coarse.txt"

echo "----- exact placement -----"
cat "$TMPD/place.txt"
echo "----- coarse placement -----"
cat "$TMPD/place-coarse.txt"


# field extractor: print value of key=<val> on the line matching <qname> <ctg> [rev=<r>]
# usage: field <file> <qname> <ctg> <key>            (first matching hit)
#        field_rev <file> <qname> <ctg> <hrev> <key> (matching hit on read-strand hrev)
field() {
	mawk -v q="$2" -v c="$3" -v k="$4" '
		$1==q && $2==c {
			for(i=1;i<=NF;i++){n=index($i,k"="); if(n==1){print substr($i,length(k)+2); exit}}
		}' "$1"
}
field_rev() {
	mawk -v q="$2" -v c="$3" -v hr="rev="$4 -v k="$5" '
		$1==q && $2==c {
			ok=0; for(i=1;i<=NF;i++) if($i==hr) ok=1;
			if(ok){ for(i=1;i<=NF;i++){n=index($i,k"="); if(n==1){print substr($i,length(k)+2); exit}} }
		}' "$1"
}

eq() { # eq <desc> <got> <want>
	if [ "$2" != "$3" ]; then fail "$1: expected '$3', got '$2'"; fi
	echo "  ok: $1 = $2"
}
within() { # within <desc> <got> <want> <tol>
	d=$(( $2 - $3 )); [ "$d" -lt 0 ] && d=$(( -d ))
	if [ "$d" -gt "$4" ]; then fail "$1: |$2 - $3| = $d > tol $4"; fi
	echo "  ok: $1 = $2 (within $4 of $3)"
}

P="$TMPD/place.txt"
C="$TMPD/place-coarse.txt"

echo "== case (b): primary chrP hit is identity passthrough =="
eq "pa-fwd/chrP pri"      "$(field "$P" pa-fwd chrP pri)"      chrP
eq "pa-fwd/chrP lst==ts"  "$(field "$P" pa-fwd chrP lst)"      120
eq "pa-fwd/chrP len==te"  "$(field "$P" pa-fwd chrP len)"      270
eq "pa-fwd/chrP prev"     "$(field "$P" pa-fwd chrP prev)"     0
eq "pa-fwd/chrP liftable" "$(field "$P" pa-fwd chrP liftable)" 1

echo "== case (a): forward ALT hit lifts to the SAME primary footprint =="
eq "pa-fwd/chrP_alt pri"      "$(field "$P" pa-fwd chrP_alt pri)"      chrP
eq "pa-fwd/chrP_alt prev"     "$(field "$P" pa-fwd chrP_alt prev)"     0
eq "pa-fwd/chrP_alt liftable" "$(field "$P" pa-fwd chrP_alt liftable)" 1
within "pa-fwd/chrP_alt lst groups with chrP" "$(field "$P" pa-fwd chrP_alt lst)" "$(field "$P" pa-fwd chrP lst)" 10

echo "== case (c): reverse-strand ALT folds; groups with the forward primary =="
# c-fwd: chrP hit (primary, read forward) vs chrR_alt hit (read on minus strand).
# block_rev(chrR_alt)=1 XOR h->rev(1) = 0 -> same folded strand as the chrP hit.
eq "c-fwd/chrP prev"           "$(field_rev "$P" c-fwd chrP 0 prev)"     0
eq "c-fwd/chrP lst"            "$(field_rev "$P" c-fwd chrP 0 lst)"      300
eq "c-fwd/chrR_alt pri"        "$(field_rev "$P" c-fwd chrR_alt 1 pri)"  chrP
eq "c-fwd/chrR_alt liftable"   "$(field_rev "$P" c-fwd chrR_alt 1 liftable)" 1
eq "c-fwd/chrR_alt folded rev" "$(field_rev "$P" c-fwd chrR_alt 1 prev)" 0
within "c-fwd/chrR_alt lst groups with chrP (strand-folded)" \
	"$(field_rev "$P" c-fwd chrR_alt 1 lst)" "$(field_rev "$P" c-fwd chrP 0 lst)" 10

echo "== case (d): one footprint end in a hole; walk to first liftable base =="
# hole-end on chrP_alt starts at ALT 300 (inside the 5I insertion -> a hole).
# Naive lift(alt_st=300) returns a hole -> liftable=0; correct impl walks to the
# first liftable base (ALT 305 -> chrP 300) so liftable=1 and lst==300.
eq "hole-end/chrP_alt pri"      "$(field "$P" hole-end chrP_alt pri)"      chrP
eq "hole-end/chrP_alt liftable" "$(field "$P" hole-end chrP_alt liftable)" 1
eq "hole-end/chrP_alt prev"     "$(field "$P" hole-end chrP_alt prev)"     0
eq "hole-end/chrP_alt lst (first liftable base)" "$(field "$P" hole-end chrP_alt lst)" 300

echo "== twin with small (<=tol) read-vs-ALT indel still groups =="
within "twin-indel chrP_alt groups with chrP primary" \
	"$(field "$P" twin-indel chrP_alt lst)" "$(field "$P" twin-indel chrP lst)" 10

echo "== coarse (pre-DP, h->p==NULL) agrees with exact on lifted_st (within tol) =="
within "coarse hole-end/chrP_alt lst == exact" \
	"$(field "$C" hole-end chrP_alt lst)" "$(field "$P" hole-end chrP_alt lst)" 10
eq     "coarse hole-end/chrP_alt liftable" "$(field "$C" hole-end chrP_alt liftable)" 1
within "coarse c-fwd/chrR_alt lst == exact (reverse fold)" \
	"$(field_rev "$C" c-fwd chrR_alt 1 lst)" "$(field_rev "$P" c-fwd chrR_alt 1 lst)" 10

echo "[test-hitplace] PASS"
exit 0
