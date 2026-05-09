#!/bin/sh
# End-to-end smoke test for --meth SAM tag emission.
#
# Builds a small synthetic reference with two CpG-rich blocks and a PE
# read pair simulating an OT-source fragment with methylated CpGs:
# - R1 reads the top strand directly (BS preserves Cs at methylated CpGs).
# - R2 reads the bottom strand (revcomp displayed in SAM).
#
# Verifies XR/XG (Bismark OT encoding) and XM (Z at every CpG cytosine).
#
# Run from the repo root after `make`.

set -e

MINIBWA="${MINIBWA:-./minibwa}"
TMP="${TMPDIR:-/tmp}/minibwa-meth-tags.$$"
mkdir -p "$TMP"
trap "rm -rf $TMP" EXIT

cat > "$TMP/ref.fa" <<'EOF'
>chr1
AAAAAAAAAACGCGCGCGAAAAAAAAAACGCGCGCGAAAAAAAAAACGCGCGCG
EOF

cat > "$TMP/r1.fq" <<'EOF'
@read_meth/1
AAAAAAAAAACGCGCGCGAAAAAAAAAACGCGCGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

cat > "$TMP/r2.fq" <<'EOF'
@read_meth/2
CGCGCGCGTTTTTTTTTTCGCGCGCGTTTTTTTTTT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

"$MINIBWA" index --meth "$TMP/ref.fa" 2>/dev/null
"$MINIBWA" map --meth -a "$TMP/ref.fa" "$TMP/r1.fq" "$TMP/r2.fq" 2>/dev/null > "$TMP/out.sam"

# minibwa emits SAM records in segment order (R1 then R2 per query), so
# the first non-header record is R1 and the second is R2.
R1=$(grep -v '^@' "$TMP/out.sam" | sed -n '1p')
R2=$(grep -v '^@' "$TMP/out.sam" | sed -n '2p')

# OT fragment, R1: XR=CT, XG=CT.
echo "$R1" | grep -q 'XR:Z:CT' || { echo "FAIL: R1 XR != CT" >&2; echo "$R1"; exit 1; }
echo "$R1" | grep -q 'XG:Z:CT' || { echo "FAIL: R1 XG != CT" >&2; echo "$R1"; exit 1; }

# OT fragment, R2: XR=GA, XG=CT (Bismark CTOT-style encoding for R2 of OT).
echo "$R2" | grep -q 'XR:Z:GA' || { echo "FAIL: R2 XR != GA" >&2; echo "$R2"; exit 1; }
echo "$R2" | grep -q 'XG:Z:CT' || { echo "FAIL: R2 XG != CT" >&2; echo "$R2"; exit 1; }

# XM: Z at top-strand C of every CpG (positions 10,12,14,16, 27,29,31,33).
EXPECTED_XM='..........Z.Z.Z.Z...........Z.Z.Z.Z'
for rec in "$R1" "$R2"; do
	XM=$(echo "$rec" | grep -oE 'XM:Z:[^[:space:]]+' | sed 's/XM:Z://')
	# trailing-position behavior is reference-end dependent; check the prefix.
	prefix=$(echo "$XM" | cut -c1-${#EXPECTED_XM})
	if [ "$prefix" != "$EXPECTED_XM" ]; then
		echo "FAIL: XM prefix mismatch" >&2
		echo "  expected: $EXPECTED_XM" >&2
		echo "  got:      $prefix" >&2
		exit 1
	fi
done

echo "PASS: --meth XR/XG/XM tags emit correctly for a methylated OT fragment."
