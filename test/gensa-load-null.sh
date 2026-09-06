#!/bin/sh
# Regression test: `minibwa gensa` (non-raw) must not crash when mb_bwt_load()
# returns NULL. mb_bwt_load() returns NULL for an unreadable/missing file or a
# corrupt BWT (bad magic, short header, ...). Before the guard, main_gensa passed
# that NULL straight into mb_bwt_gen_sa(), which dereferences it -> SIGSEGV (exit
# 139). With the guard it reports the offending path and returns 1. The raw
# loader (`-r` -> mb_bwt_load_raw) has a different, non-NULL contract and is
# intentionally out of scope here. Generates all inputs inline (no committed data).
set -eu

MINIBWA="${MINIBWA:-./minibwa}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail() { echo "FAIL: $*" >&2; exit 1; }

# Run gensa and capture the exit status without letting `set -e` abort us.
run_gensa() { # <in.bwt> <out.bwt>
	set +e
	"$MINIBWA" gensa "$1" "$2" >/dev/null 2>&1
	rc=$?
	set -e
}

echo "# (a) a missing input BWT must exit 1, not crash"
run_gensa "$WORK/does-not-exist.bwt" "$WORK/out-a.bwt"
[ "$rc" -eq 139 ] && fail "(a) gensa crashed (SIGSEGV) on a missing BWT file"
[ "$rc" -eq 1 ] || fail "(a) gensa on a missing BWT should exit 1, got $rc"

echo "# (b) a corrupt input BWT (bad magic) must exit 1, not crash"
printf 'NOTAMAGIC and some junk bytes that are not a valid minibwa BWT header' > "$WORK/corrupt.bwt"
run_gensa "$WORK/corrupt.bwt" "$WORK/out-b.bwt"
[ "$rc" -eq 139 ] && fail "(b) gensa crashed (SIGSEGV) on a corrupt BWT file"
[ "$rc" -eq 1 ] || fail "(b) gensa on a corrupt BWT should exit 1, got $rc"

echo "# (c) a valid BWT still gensa's cleanly (exit 0) -- guard doesn't break the happy path"
ref="$WORK/ref.fa"
awk 'BEGIN{
	b="ACGT"; x=1; printf ">chrA\n";
	for(i=0;i<2000;i++){ x=(1103515245*x+12345)%2147483648; printf "%s", substr(b,(x%4)+1,1);
		if((i+1)%60==0) printf "\n"; }
	printf "\n";
}' > "$ref"
"$MINIBWA" index "$ref" "$WORK/idx" >/dev/null 2>&1 || fail "(c) index build failed"
run_gensa "$WORK/idx.mbw" "$WORK/out-c.mbw"
[ "$rc" -eq 0 ] || fail "(c) gensa on a valid BWT should exit 0, got $rc"

echo "PASS: gensa NULL-load guard"
