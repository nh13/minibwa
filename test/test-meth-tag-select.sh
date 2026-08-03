#!/bin/sh
# Tests for --meth-tags, the selector for which of the Bismark XR/XG/XM tags
# --meth emits.
#
# The properties under test:
#  - "all" is the default, and a bare --meth is unchanged by this option.
#  - Every spec form (all, none, inclusion list, ^exclusion list) selects the
#    tag set it names.
#  - A spec is a *set*: order and repetition are immaterial, names are
#    case-insensitive, and ^XM is exactly XR,XG.
#  - Inclusion and exclusion forms cannot be mixed; a mixed or otherwise
#    malformed spec exits non-zero instead of quietly falling back.
#  - Selection is an emission filter and not an alignment knob: SAM fields
#    1-11 and the NM/MD/AS tags are byte-identical under every spec.
#
# Uses the bundled chrM reference and read pairs, so the invariance check runs
# over a thousand real alignments rather than a hand-built fixture.
#
# Run from the repo root after `make`.

set -e

MINIBWA="${MINIBWA:-./minibwa}"
TEST_DIR="$(dirname "$0")"
TMP="${TMPDIR:-/tmp}/minibwa-meth-tag-select.$$"
mkdir -p "$TMP"
trap "rm -rf $TMP" EXIT

gzip -dc "$TEST_DIR/chrM-human.fa.gz" > "$TMP/chrM.fa"
gzip -dc "$TEST_DIR/chrM-read_1.fa.gz" > "$TMP/r1.fa"
gzip -dc "$TEST_DIR/chrM-read_2.fa.gz" > "$TMP/r2.fa"

"$MINIBWA" index --meth "$TMP/chrM.fa" 2>/dev/null

fail() { echo "FAIL: $*" >&2; exit 1; }

# Map with the given --meth-tags arguments and write the SAM records (no
# header) to $2. -b MD is on so the invariance check can see an MD tag.
map_tags() {
	out="$1"; shift
	"$MINIBWA" map --meth -b MD "$@" -a "$TMP/chrM.fa" "$TMP/r1.fa" "$TMP/r2.fa" \
		> "$TMP/raw.sam" 2>/dev/null || fail "'minibwa map --meth $*' exited non-zero"
	# grep exits 1 on an all-header SAM, which is not a test failure here.
	grep -v '^@' "$TMP/raw.sam" > "$out" || :
}

# The sorted, de-duplicated set of Bismark tag names present in a SAM file,
# as a space-separated string; empty if there are none.
tags_in() {
	grep -oE 'X[RGM]:Z:' "$1" 2>/dev/null | sed 's/:Z://' | sort -u | tr '\n' ' ' \
		| sed 's/ $//'
}

# SAM fields 1-11 plus the NM/MD/AS tags: everything that must not move when
# only the methylation tag selection changes.
project_invariant() {
	awk -F'\t' '{
		printf "%s", $1
		for (i = 2; i <= 11; ++i) printf "\t%s", $i
		for (i = 12; i <= NF; ++i)
			if ($i ~ /^NM:/ || $i ~ /^MD:/ || $i ~ /^AS:/) printf "\t%s", $i
		print ""
	}' "$1"
}

# Assert that a spec succeeds and yields exactly the given tag set.
check_tags() {
	spec="$1"; want="$2"
	map_tags "$TMP/sel.sam" --meth-tags "$spec"
	got=$(tags_in "$TMP/sel.sam")
	[ "$got" = "$want" ] || fail "--meth-tags '$spec': expected tags [$want], got [$got]"
	# ... and that it changed nothing but the tags.
	project_invariant "$TMP/sel.sam" > "$TMP/sel.inv"
	cmp -s "$TMP/base.inv" "$TMP/sel.inv" || \
		fail "--meth-tags '$spec' altered SAM fields 1-11 or NM/MD/AS"
}

# Assert that a spec is rejected: non-zero exit, and no SAM on stdout.
check_rejected() {
	spec="$1"
	rc=0
	"$MINIBWA" map --meth --meth-tags "$spec" -a "$TMP/chrM.fa" "$TMP/r1.fa" \
		> "$TMP/rej.sam" 2>/dev/null || rc=$?
	[ "$rc" -ne 0 ] || fail "--meth-tags '$spec' should have exited non-zero"
	[ ! -s "$TMP/rej.sam" ] || fail "--meth-tags '$spec' exited non-zero but still wrote SAM"
}

# ---------------------------------------------------------------------------
# The default is "all", and it is unchanged by the existence of this option.
map_tags "$TMP/base.sam"
project_invariant "$TMP/base.sam" > "$TMP/base.inv"
[ -s "$TMP/base.sam" ] || fail "baseline --meth run produced no alignments"

got=$(tags_in "$TMP/base.sam")
[ "$got" = "XG XM XR" ] || fail "bare --meth should emit all three tags, got [$got]"

map_tags "$TMP/all.sam" --meth-tags all
cmp -s "$TMP/base.sam" "$TMP/all.sam" || fail "'--meth-tags all' differs from a bare --meth"

# ---------------------------------------------------------------------------
# Every spec form selects the set it names. (tags_in sorts, so the expected
# set is always in XG XM XR order regardless of emission order.)
check_tags "all"      "XG XM XR"
check_tags "none"     ""
check_tags "XR"       "XR"
check_tags "XG"       "XG"
check_tags "XM"       "XM"
check_tags "XR,XG"    "XG XR"
check_tags "XR,XG,XM" "XG XM XR"
check_tags "^XM"      "XG XR"
check_tags "^XR"      "XG XM"
check_tags "^XR,XG"   "XM"
check_tags "^XR,XG,XM" ""

# ---------------------------------------------------------------------------
# Set semantics: order and repetition are immaterial, and the exclusion form
# is just another spelling of the same set.
map_tags "$TMP/incl.sam" --meth-tags XR,XG
for equivalent in "XG,XR" "XR,XG,XR" "XR,XR,XG" "^XM"; do
	map_tags "$TMP/eq.sam" --meth-tags "$equivalent"
	cmp -s "$TMP/incl.sam" "$TMP/eq.sam" || \
		fail "--meth-tags '$equivalent' should be identical to 'XR,XG'"
done

# Tag names and keywords are case-insensitive.
for spelling in "xr,xg" "Xr,xG" "XR,Xg"; do
	map_tags "$TMP/case.sam" --meth-tags "$spelling"
	cmp -s "$TMP/incl.sam" "$TMP/case.sam" || \
		fail "--meth-tags '$spelling' should be identical to 'XR,XG'"
done
for spelling in "ALL" "All"; do
	map_tags "$TMP/case.sam" --meth-tags "$spelling"
	cmp -s "$TMP/base.sam" "$TMP/case.sam" || \
		fail "--meth-tags '$spelling' should be identical to 'all'"
done
check_tags "NONE" ""
check_tags "^xm"  "XG XR"

# ---------------------------------------------------------------------------
# Mixing inclusion and exclusion is rejected rather than resolved: "XR,^XM"
# could mean "only XR" or "everything but XM", and guessing either way would
# silently produce the wrong tag set.
check_rejected "XR,^XM"
check_rejected "^XR,^XM"
check_rejected "XR,XG,^XM"

# Malformed specs exit non-zero; none of them may fall back to the default.
check_rejected ""          # empty spec
check_rejected "^"         # bare exclusion marker
check_rejected "XQ"        # not a methylation tag
check_rejected "NM"        # a real SAM tag, but not one this option selects
check_rejected "X"         # too short
check_rejected "XRR"       # too long
check_rejected "XR,"       # trailing comma
check_rejected ",XR"       # leading comma
check_rejected "XR,,XG"    # empty list element
check_rejected "XR,all"    # keywords are whole-spec, not list members
check_rejected "^all"
check_rejected "^none"
check_rejected "XR XG"     # the separator is a comma, never a space

echo "PASS: --meth-tags selects XR/XG/XM and defaults to all"
echo "PASS: --meth-tags specs are sets: order, repetition and case immaterial"
echo "PASS: --meth-tags rejects mixed and malformed specs with a non-zero exit"
echo "PASS: --meth-tags leaves SAM fields 1-11 and NM/MD/AS byte-identical"
