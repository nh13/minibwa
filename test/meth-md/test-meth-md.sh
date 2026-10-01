#!/bin/sh
# MD must list exactly the differences NM counts, with and without --meth.
#
# With --meth the scoring matrix treats a read's bisulfite/enzymatic conversion
# (T at a reference C, or A at a reference G) as a match, and NM counts only
# differences the matrix penalizes, plus ambiguous bases. MD must list exactly
# those bases, so CIGAR + SEQ + MD reconstructs the converted reference that NM
# was computed against.
#
# This converts the bundled chrM reads as a fully unmethylated directional
# library (R1 C->T, R2 G->A) and adds single-end reads over the chrM reference N
# (a read N on it, a base on it, and a read N elsewhere). Every mapped primary
# record is checked against the reference:
#   - NM = MD mismatches + MD deleted bases + CIGAR inserted bases;
#   - every base MD lists is the reference base there, and every ambiguous base
#     (read or reference) is listed;
#   - every difference MD leaves out is the read's conversion for its strand
#     under --meth (R1 forward / R2 reverse: C->T; R1 reverse / R2 forward:
#     G->A, in reference orientation), and there are none without --meth.
# Mapping with --eqx must not change MD. The unconverted reads mapped without
# --meth are checked the same way as a control.
#
# Usage: sh test/meth-md/test-meth-md.sh [<minibwa-dir>]

set -e

DIR="${1:-.}"
MINIBWA="${MINIBWA:-$DIR/minibwa}"
TMP="${TMPDIR:-/tmp}/minibwa-meth-md.$$"
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

gzip -dc "$DIR/test/chrM-human.fa.gz" > "$TMP/ref.fa"
gzip -dc "$DIR/test/chrM-read_1.fa.gz" > "$TMP/r_1.fa"
gzip -dc "$DIR/test/chrM-read_2.fa.gz" > "$TMP/r_2.fa"
# Convert sequence lines only (FASTA = alternating header/sequence lines).
awk 'NR%2==0{gsub(/[Cc]/,"T")}1' "$TMP/r_1.fa" > "$TMP/bs_1.fa"
awk 'NR%2==0{gsub(/[Gg]/,"A")}1' "$TMP/r_2.fa" > "$TMP/bs_2.fa"
# The reference on one line, uppercased, for the checks below.
grep -v '^>' "$TMP/ref.fa" | tr -d '\n' | tr 'acgtn' 'ACGTN' > "$TMP/ref.txt"

# Single-end reads over the reference N at 1-based 3107: one with N on it, one
# with A on it, and one with an extra read N 20 bases later. The --meth copies
# are converted as R1 (C->T), leaving the N.
awk -v p=3107 '{
	s = substr($0, p - 75, 151)
	print ">nn"; print s
	print ">na"; print substr(s, 1, 75) "A" substr(s, 77)
	print ">rn"; print substr(s, 1, 95) "N" substr(s, 97)
}' "$TMP/ref.txt" > "$TMP/ns.fa"
awk 'NR%2==0{gsub(/C/,"T")}1' "$TMP/ns.fa" > "$TMP/ns_bs.fa"

"$MINIBWA" index "$TMP/ref.fa" 2>/dev/null
"$MINIBWA" index --meth "$TMP/ref.fa" 2>/dev/null
for eqx in "" --eqx; do
	"$MINIBWA" map --meth -b MD $eqx "$TMP/ref.fa" "$TMP/bs_1.fa" "$TMP/bs_2.fa" 2>/dev/null > "$TMP/meth$eqx.sam"
done
cp "$TMP/meth.sam" "$TMP/meth-pe.sam"
"$MINIBWA" map --meth -b MD "$TMP/ref.fa" "$TMP/ns_bs.fa" 2>/dev/null | grep -v '^@' >> "$TMP/meth.sam"
"$MINIBWA" map -b MD "$TMP/ref.fa" "$TMP/r_1.fa" "$TMP/r_2.fa" 2>/dev/null > "$TMP/plain.sam"
"$MINIBWA" map -b MD "$TMP/ref.fa" "$TMP/ns.fa" 2>/dev/null | grep -v '^@' >> "$TMP/plain.sam"

# Prints "<checked> <failing> <mean MD edits> <N/N records> <first failure>" for a SAM file;
# meth=1 allows the read's strand conversion to be left out of MD.
md_check () {
	grep -v '^@' "$1" | awk -v meth="$2" -v reffile="$TMP/ref.txt" '
	BEGIN { getline ref < reffile }
	function fail(why) { ++bad; if (!ex) ex = $1 " " why " NM=" nm " MD=" md " CIGAR=" $6 }
	{
		flag = $2
		if (int(flag/4)%2 || int(flag/256)%2 || int(flag/2048)%2) next
		nm = -1; md = ""
		for (i = 12; i <= NF; ++i) {
			if ($i ~ /^NM:i:/) nm = substr($i, 6) + 0
			else if ($i ~ /^MD:Z:/) md = substr($i, 6)
		}
		++n
		if (nm < 0 || md == "") { fail("missing NM/MD"); next }
		# Expand MD: "." per matching base, the reference letter per mismatch, and
		# the deleted reference letters queued separately.
		mdx = ""; dels = ""; m = md
		while (m != "") {
			if (match(m, /^[0-9]+/)) {
				for (k = substr(m, 1, RLENGTH) + 0; k > 0; --k) mdx = mdx "."
			} else if (match(m, /^\^[A-Z]+/)) {
				dels = dels substr(m, 2, RLENGTH - 1)
			} else if (match(m, /^[A-Z]/)) {
				mdx = mdx substr(m, 1, 1)
			} else { fail("bad MD"); next }
			m = substr(m, RLENGTH + 1)
		}
		edits = length(dels); ins = 0
		letters = mdx; gsub(/\./, "", letters); edits += length(letters)
		top = (int(flag/16)%2 == int(flag/128)%2)
		rp = $4; qp = 1; mi = 1; di = 1; cig = $6; nn = 0
		while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
			len = substr(cig, 1, RLENGTH - 1) + 0; op = substr(cig, RLENGTH, 1)
			cig = substr(cig, RLENGTH + 1)
			if (op == "M" || op == "=" || op == "X") {
				for (k = 0; k < len; ++k) {
					r = substr(ref, rp + k, 1); q = substr($10, qp + k, 1); e = substr(mdx, mi + k, 1)
					if (r == "N" && q == "N" && !nn++) ++nnrec
					if (e != ".") {
						if (e != r) { fail("MD letter " e " is not reference " r " at " rp + k); next }
					} else if (r == "N" || q == "N") {
						fail("ambiguous base at " rp + k " not in MD"); next
					} else if (q != r && !(meth && (top ? r == "C" && q == "T" : r == "G" && q == "A"))) {
						fail("difference " r ">" q " at " rp + k " not in MD"); next
					}
				}
				rp += len; qp += len; mi += len
			} else if (op == "D") {
				if (substr(dels, di, len) != substr(ref, rp, len)) { fail("MD deletion is not the reference"); next }
				rp += len; di += len
			} else if (op == "N") {
				rp += len
			} else if (op == "I") {
				qp += len; ins += len
			} else if (op == "S") {
				qp += len
			}
		}
		if (mi - 1 != length(mdx) || di - 1 != length(dels)) { fail("MD does not cover the CIGAR"); next }
		if (nm != edits + ins) { fail("NM != MD edits + insertions"); next }
		total += edits
	} END { printf "%d %d %.3f %d %s", n, bad + 0, (n ? total/n : 999), nnrec + 0, ex }'
}

for mode in meth plain; do
	meth=0; [ "$mode" = meth ] && meth=1
	set -- $(md_check "$TMP/$mode.sam" $meth)
	n=$1; bad=$2; mean=$3; nnrec=$4; shift 4
	echo "  $mode: $n primary alignments ($nnrec with N on reference N), $bad failing, mean MD edits $mean"
	if [ "$n" -lt 1900 ]; then
		echo "FAIL: $mode mapped only $n reads" >&2; exit 1
	fi
	if [ "$nnrec" -lt 1 ]; then
		echo "FAIL: $mode did not align the read with N on the reference N" >&2; exit 1
	fi
	if [ "$bad" -ne 0 ]; then
		echo "FAIL: $mode MD wrong, e.g. $*" >&2; exit 1
	fi
	# A converted 151 bp read carries ~30 conversions; listing them would put
	# the mean far above this bound.
	if ! awk -v v="$mean" 'BEGIN{exit !(v < 2.0)}'; then
		echo "FAIL: $mode mean MD edits ($mean) too high" >&2; exit 1
	fi
done

# --eqx rewrites the CIGAR (X still marks every raw difference) but not MD.
md_by_read () { grep -v '^@' "$1" | awk '{ for (i = 12; i <= NF; ++i) if ($i ~ /^MD:Z:/) print $1, $2, $i }' | sort; }
md_by_read "$TMP/meth-pe.sam" > "$TMP/md.txt"
md_by_read "$TMP/meth--eqx.sam" > "$TMP/md-eqx.txt"
if ! cmp -s "$TMP/md.txt" "$TMP/md-eqx.txt"; then
	echo "FAIL: --eqx changed MD under --meth" >&2; exit 1
fi

echo "PASS: MD lists exactly the differences NM counts, with and without --meth"
