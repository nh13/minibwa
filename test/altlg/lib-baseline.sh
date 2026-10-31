# Shared no-.alt baseline check for the ALT liftover-group suite.
#
# What the ALT feature must guarantee on a reference with no .alt file is that its
# output is *identical to stock minibwa's* -- the hooks are all gated, so an unpatched
# and a patched binary should agree byte for byte.
#
# Running this build twice and comparing cannot show that: it only shows the aligner is
# deterministic, which is true whether or not the hooks leak. Seven scripts here used to
# do exactly that and reported it as "inert without .alt". The comparison needs a second,
# unpatched binary, so it is opt-in: set MB_STOCK to a stock minibwa and the check runs;
# leave it unset and the check reports itself skipped rather than passing vacuously.
#
# Callers must define: MINIBWA, MDIR, TMPD, and the ok()/fail() helpers.
#
#   chrm_baseline "(BASELINE) chrM PE" pe --outn=5
#   chrm_baseline "(BASELINE) chrM SE" se

chrm_baseline() {
	_label="$1"; _mode="$2"; shift 2
	_fa="$MDIR/test/chrM-human.fa.gz"
	_r1="$MDIR/test/chrM-read_1.fa.gz"
	_r2="$MDIR/test/chrM-read_2.fa.gz"

	if [ ! -f "$_fa" ] || [ ! -f "$_r1" ]; then
		echo "  skip: $_label -- chrM fixture not found"
		return 0
	fi
	if [ "$_mode" = pe ] && [ ! -f "$_r2" ]; then
		echo "  skip: $_label -- chrM mate file not found"
		return 0
	fi
	if [ -z "${MB_STOCK:-}" ]; then
		echo "  skip: $_label -- set MB_STOCK=<stock minibwa> to compare against unpatched output"
		return 0
	fi
	if [ ! -x "$MB_STOCK" ]; then
		fail "$_label -- MB_STOCK='$MB_STOCK' is not executable"
	fi

	# Copy and index inside TMPD rather than beside the shared fixture: keeps the check
	# hermetic, avoids racing another script indexing the same chrM, and leaves no
	# untracked .mbw/.l2b in test/.
	cp "$_fa" "$TMPD/base-ref.fa.gz"
	if [ "$_mode" = pe ]; then set -- "$@" "$TMPD/base-ref.fa.gz" "$_r1" "$_r2"
	else set -- "$@" "$TMPD/base-ref.fa.gz" "$_r1"; fi

	"$MINIBWA" index "$TMPD/base-ref.fa.gz" 2>/dev/null
	"$MINIBWA" map "$@" 2>/dev/null | grep -v '^@PG' > "$TMPD/base-patched.sam"
	"$MB_STOCK" map "$@" 2>/dev/null | grep -v '^@PG' > "$TMPD/base-stock.sam"

	[ -s "$TMPD/base-patched.sam" ] || fail "$_label -- empty output from this build"
	[ -s "$TMPD/base-stock.sam" ]   || fail "$_label -- empty output from MB_STOCK"

	# @PG carries the command line and version, so it differs by construction and is
	# stripped above; every other line must match.
	if cmp -s "$TMPD/base-patched.sam" "$TMPD/base-stock.sam"; then
		ok "$_label byte-identical to stock (hooks provably inert without .alt)"
	else
		fail "$_label differs from stock -- an ALT hook is firing on a no-.alt reference"
	fi
}
