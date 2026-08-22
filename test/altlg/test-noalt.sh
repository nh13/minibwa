#!/bin/sh
# Integration test for --no-alt, and for .alt resolution being independent of
# how the index was loaded.
#
# The fixture (mkfixture.sh) is a primary chrP plus an ALT twin chrP_alt that
# differs only by a 5bp insertion, and a read landing in the shared flank.  The
# read therefore hits both contigs equally well, which is precisely the case
# ALT-awareness exists for: with the .alt loaded the two hits are recognised as
# one locus and the primary keeps MAPQ 60; without it they look like a genuine
# repeat and MAPQ collapses to 0.  That 60-vs-0 split is the signal every check
# below keys on -- it is unambiguous and it is the user-visible consequence.
#
# Verifies:
#   1. Liveness   -- ALT-on and --no-alt really do differ (MAPQ 60 vs 0).
#   2. Equivalence-- --no-alt is byte-identical to deleting the .alt outright.
#   3. Precedence -- --no-alt wins over an explicit --alt FILE, either order.
#   4. mmap parity-- --mmap agrees with the normal loader, ALT on AND off.
#      Before .alt resolution moved into the caller, mb_idx_load auto-detected
#      <prefix>.alt but mb_idx_load_mmap did not, so --mmap silently disabled
#      ALT-awareness on the very same index.  This check is that regression.
#
# Usage: test/altlg/test-noalt.sh [<minibwa-dir>]
set -eu

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MINIBWA="$MDIR/minibwa"

TMPD=$(mktemp -d /tmp/altlg-noalt.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

fail() { echo "FAIL: $1"; exit 1; }
ok()   { echo "  ok: $1"; }

# Alignment records only: @PG carries the command line, which legitimately
# differs between two runs we expect to align identically.
sam_body() { mawk '$1!~/^@/{print}' "$1"; }

# MAPQ of the primary (neither secondary 0x100 nor supplementary 0x800).
# mawk has no and(), so test the two bits by arithmetic.
primary_mapq() {
	mawk '$1!~/^@/ && int($2/256)%2==0 && int($2/2048)%2==0 {print $5; exit}' "$1"
}

# Run `map` with the fixture and the given extra flags, into $TMPD/$1.sam.
run() {
	out="$1"; shift
	"$MINIBWA" map --outn=5 "$@" "$TMPD/ref.fa" "$TMPD/r1.fq" \
		2>/dev/null > "$TMPD/$out.sam" \
		|| fail "minibwa map failed for '$out' (flags: $*)"
	[ -s "$TMPD/$out.sam" ] || fail "empty output for '$out' (flags: $*)"
	sam_body "$TMPD/$out.sam" > "$TMPD/$out.body"
}

echo "[test-noalt] building fixture ..."
sh "$MDIR/test/altlg/mkfixture.sh" "$TMPD"
[ -f "$TMPD/ref.fa.alt" ] || fail "fixture did not produce ref.fa.alt"
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null || fail "indexing the fixture failed"

# ============================================================
# 1. Liveness: the flag actually does something
# ============================================================
echo "== liveness =="
run alt-on
run alt-off --no-alt

mq_on=$(primary_mapq "$TMPD/alt-on.sam")
mq_off=$(primary_mapq "$TMPD/alt-off.sam")
[ "${mq_on:-}" = "60" ] \
	|| fail "ALT on: primary MAPQ is '${mq_on:-<none>}', expected 60 (is the .alt being loaded at all?)"
ok "ALT on: primary MAPQ 60 (ALT twin recognised, MAPQ preserved)"
[ "${mq_off:-}" = "0" ] \
	|| fail "--no-alt: primary MAPQ is '${mq_off:-<none>}', expected 0 (ALT-awareness was not suppressed)"
ok "--no-alt: primary MAPQ 0 (ALT twin looks like a repeat, as it should)"

cmp -s "$TMPD/alt-on.body" "$TMPD/alt-off.body" \
	&& fail "--no-alt produced byte-identical output to ALT-on; the flag is a no-op"
ok "--no-alt changes the alignment records"

# ============================================================
# 2. Equivalence: --no-alt == no .alt file on disk
# ============================================================
echo "== equivalence with a deleted .alt =="
mv "$TMPD/ref.fa.alt" "$TMPD/hidden.alt"
run alt-absent
mv "$TMPD/hidden.alt" "$TMPD/ref.fa.alt"

cmp -s "$TMPD/alt-off.body" "$TMPD/alt-absent.body" \
	|| fail "--no-alt differs from running against an index with no .alt beside it"
ok "--no-alt is byte-identical to having no .alt at all"

# ============================================================
# 3. Precedence: --no-alt beats an explicit --alt, in either order
# ============================================================
echo "== precedence over --alt FILE =="
cp "$TMPD/ref.fa.alt" "$TMPD/elsewhere.alt"

run alt-explicit --alt "$TMPD/elsewhere.alt"
cmp -s "$TMPD/alt-explicit.body" "$TMPD/alt-on.body" \
	|| fail "--alt FILE with the same content differs from the auto-detected .alt"
ok "--alt FILE reproduces the auto-detected result (control)"

run alt-then-no --alt "$TMPD/elsewhere.alt" --no-alt
run no-then-alt --no-alt --alt "$TMPD/elsewhere.alt"
cmp -s "$TMPD/alt-then-no.body" "$TMPD/alt-off.body" \
	|| fail "'--alt FILE --no-alt' did not suppress ALT-awareness"
ok "'--alt FILE --no-alt' suppresses ALT-awareness"
cmp -s "$TMPD/no-then-alt.body" "$TMPD/alt-off.body" \
	|| fail "'--no-alt --alt FILE' did not suppress ALT-awareness (order-dependent!)"
ok "'--no-alt --alt FILE' suppresses ALT-awareness (order-independent)"

# ============================================================
# 4. mmap parity, ALT on and ALT off
# ============================================================
echo "== --mmap parity =="
run mmap-on  --mmap
run mmap-off --mmap --no-alt

mq_mmap=$(primary_mapq "$TMPD/mmap-on.sam")
[ "${mq_mmap:-}" = "60" ] \
	|| fail "--mmap with an adjacent .alt: primary MAPQ is '${mq_mmap:-<none>}', expected 60 -- the mmap loader is not seeing the .alt"
cmp -s "$TMPD/mmap-on.body" "$TMPD/alt-on.body" \
	|| fail "--mmap changed the alignment records with ALT on (loader-dependent .alt resolution)"
ok "--mmap matches the normal loader with ALT on"

cmp -s "$TMPD/mmap-off.body" "$TMPD/alt-off.body" \
	|| fail "--mmap changed the alignment records under --no-alt"
ok "--mmap matches the normal loader under --no-alt"

echo "[test-noalt] PASS"
exit 0
