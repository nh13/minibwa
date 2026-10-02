#!/bin/sh
# Malformed .alt records must be rejected, not half-loaded.
#
# l2b_set_alt() parses .alt lines out of a user-supplied SAM. Every coordinate it
# derives is attacker- or accident-controlled, and one of them is a SUBTRACTION:
# a reverse record's ALT span is remapped to forward as `ctg->len - alt_en`. If
# the record's CIGAR query span runs past the ALT contig it names -- a .alt built
# against a different reference revision is enough -- that wraps on uint64 to a
# near-UINT64_MAX coordinate, the block sorts to the end of lift[], and
# l2b_lift()'s binary search reads it as valid. The failure is silent and it
# inflates confidence: a garbage lift makes an ALT twin look co-located with its
# primary, and MAPQ goes UP.
#
# Verifies:
#   1. Over-long reverse record -> rejected (MAPQ not inflated by a garbage lift).
#   2. POS < 1                  -> rejected.
#   3. Unknown CIGAR operator   -> rejected.
#   4. A WELL-FORMED reverse record still loads -- the guards must reject only
#      what is actually malformed, or they would silently disable ALT support.
#
# Usage: test/altlg/test-badalt.sh [<minibwa-dir>]
set -eu
. "$(dirname "$0")/lib.sh"
. "$(dirname "$0")/fixlib.sh"

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MINIBWA="$MDIR/minibwa"

TMPD=$(mktemp -d /tmp/altlg-badalt.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

# chrP (400bp primary) + chrP_alt, the reverse complement of chrP[100,250).
# A reverse .alt record is the one that exercises the subtraction.
PRI=$(gen 7 400)
ALT=$(rc "$(printf '%s' "$PRI" | cut -c101-250)")
printf '>chrP\n%s\n>chrP_alt\n%s\n' "$PRI" "$ALT" > "$TMPD/ref.fa"
printf '@r1\n%s\n+\n%s\n' \
    "$(printf '%s' "$PRI" | cut -c121-220)" "$(printf '%100s' '' | tr ' ' 'I')" > "$TMPD/reads.fq"
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null || fail "indexing the fixture failed"

# Run with the given .alt content and report the primary's MAPQ.
mapq_with() {
    printf '%s' "$1" > "$TMPD/ref.fa.alt"
    "$MINIBWA" map --outn=5 "$TMPD/ref.fa" "$TMPD/reads.fq" 2>/dev/null > "$TMPD/out.sam" \
        || fail "minibwa map failed on .alt <$1>"
    [ -s "$TMPD/out.sam" ] || fail "empty output for .alt <$1>"
    primary_mapq "$TMPD/out.sam"
}

# ============================================================
# 4 first: establish that a well-formed reverse record DOES load.
# Without this the rejection checks below could all pass for the
# trivial reason that ALT support is broken outright.
# ============================================================
echo "== well-formed reverse record (control) =="
good=$(mapq_with 'chrP_alt	16	chrP	101	60	150M	*	0	0	*	*
')
[ "$good" = "60" ] \
    || fail "well-formed reverse .alt: primary MAPQ is '${good:-<none>}', expected 60 -- ALT support is not working, so the rejection checks below would be vacuous"
ok "well-formed reverse record loads: primary MAPQ $good"

# The same read with NO .alt at all: the baseline a rejected record must match.
rm -f "$TMPD/ref.fa.alt"
"$MINIBWA" map --outn=5 "$TMPD/ref.fa" "$TMPD/reads.fq" 2>/dev/null > "$TMPD/none.sam"
none=$(primary_mapq "$TMPD/none.sam")
[ "$none" != "$good" ] \
    || fail "no-.alt MAPQ ($none) equals the loaded one ($good); this fixture cannot tell the two apart"
ok "no .alt at all: primary MAPQ $none (the baseline a rejected record must match)"

# ============================================================
# 1. Over-long reverse record: 300M on a 150bp ALT contig.
# ============================================================
echo "== malformed records are rejected =="
span=$(mapq_with 'chrP_alt	16	chrP	101	60	300M	*	0	0	*	*
')
[ "$span" = "$none" ] \
    || fail "over-long reverse record: primary MAPQ is '$span', expected '$none' -- the span was accepted and a wrapped coordinate lifted it"
ok "reverse record longer than its ALT contig is rejected (MAPQ $span, not inflated)"

# ============================================================
# 2. POS < 1 (a mapped record must be 1-based >= 1).
# ============================================================
pos=$(mapq_with 'chrP_alt	16	chrP	0	60	150M	*	0	0	*	*
')
[ "$pos" = "$none" ] \
    || fail "POS=0 record: primary MAPQ is '$pos', expected '$none'"
ok "POS < 1 is rejected (MAPQ $pos)"

# ============================================================
# 3. Unknown CIGAR operator.
# ============================================================
cig=$(mapq_with 'chrP_alt	16	chrP	101	60	100M50Z	*	0	0	*	*
')
[ "$cig" = "$none" ] \
    || fail "unknown CIGAR op: primary MAPQ is '$cig', expected '$none'"
ok "unknown CIGAR operator is rejected (MAPQ $cig)"

echo "[test-badalt] PASS"
exit 0
