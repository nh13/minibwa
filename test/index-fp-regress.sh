#!/bin/sh
# Regression test for the reference-fingerprint guard on incremental .mbw reuse
# (index.c). Proves that `minibwa index` refuses to reuse a .mbw built from a
# different reference (the silent wrong-coordinate bug), while still reusing it
# for a legitimate same-reference add-density run.
#
# Generates all FASTA content inline (no committed test data). Run from a build
# tree: `sh test/index-fp-regress.sh` (uses ./minibwa) or MINIBWA=/path/to/minibwa.
set -eu

MINIBWA="${MINIBWA:-./minibwa}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail() { echo "FAIL: $*" >&2; exit 1; }

# Deterministic pseudo-random DNA of a given length and integer seed.
gen_fa() { # <len> <seed> <name> <out.fa>
	awk -v len="$1" -v seed="$2" -v name="$3" 'BEGIN{
		b="ACGT"; x=seed+1; printf ">%s\n", name;
		for(i=0;i<len;i++){ x=(1103515245*x+12345)%2147483648; printf "%s", substr(b,(x%4)+1,1);
			if((i+1)%60==0) printf "\n"; }
		if(len%60!=0) printf "\n";
	}' > "$4"
}

mbw_md5() { if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'; else md5 -q "$1"; fi; }

refA="$WORK/refA.fa"; refB="$WORK/refB.fa"; refAp="$WORK/refAp.fa"
gen_fa 2000 1 chrA "$refA"      # reference A
gen_fa 3000 2 chrB "$refB"      # different length AND content
# Same length as A, different content (flip to a different seed): tests that a
# length-only check (rejected design A) would miss but the content hash catches.
gen_fa 2000 9 chrA "$refAp"

pfx="$WORK/idx"

echo "# (a) differing-length reference swap must be rejected"
"$MINIBWA" index "$refA" "$pfx" >/dev/null 2>&1
[ -f "$pfx.mbw.fp" ] || fail "(a) no .mbw.fp written on fresh build"
md5_A="$(mbw_md5 "$pfx.mbw")"
err="$("$MINIBWA" index "$refB" "$pfx" 2>&1 >/dev/null)"
echo "$err" | grep -q "doesn't match the current reference" || fail "(a) no mismatch warning on ref swap"
[ "$(mbw_md5 "$pfx.mbw")" != "$md5_A" ] || fail "(a) .mbw not rebuilt after ref swap"

echo "# (b) same-length, different-content swap must be rejected (content hash, not length)"
rm -f "$pfx".* ; "$MINIBWA" index "$refA" "$pfx" >/dev/null 2>&1
md5_A="$(mbw_md5 "$pfx.mbw")"
err="$("$MINIBWA" index "$refAp" "$pfx" 2>&1 >/dev/null)"
echo "$err" | grep -q "doesn't match the current reference" || fail "(b) no mismatch warning on same-length swap"
[ "$(mbw_md5 "$pfx.mbw")" != "$md5_A" ] || fail "(b) .mbw not rebuilt after same-length swap"

echo "# (c) legitimate same-reference add-density reuse must NOT rebuild"
rm -f "$pfx".* ; "$MINIBWA" index -u 4 "$refA" "$pfx" >/dev/null 2>&1
md5_A="$(mbw_md5 "$pfx.mbw")"
err="$("$MINIBWA" index -u 3,4 "$refA" "$pfx" 2>&1 >/dev/null)"
echo "$err" | grep -q "doesn't match" && fail "(c) spurious mismatch warning on same-ref reuse"
[ "$(mbw_md5 "$pfx.mbw")" = "$md5_A" ] || fail "(c) .mbw rebuilt despite unchanged reference"
[ -f "$pfx.sa.u3" ] || fail "(c) new sidecar sa.u3 not created on reuse"

echo "# (d) a .mbw with no fingerprint must rebuild conservatively"
rm -f "$pfx.mbw.fp"
err="$("$MINIBWA" index "$refA" "$pfx" 2>&1 >/dev/null)"
echo "$err" | grep -q "fingerprint is missing" || fail "(d) missing .fp did not force a rebuild"
[ -f "$pfx.mbw.fp" ] || fail "(d) .fp not recreated"

echo "# (e) a corrupt/truncated fingerprint must rebuild conservatively"
printf 'XX' > "$pfx.mbw.fp"
err="$("$MINIBWA" index "$refA" "$pfx" 2>&1 >/dev/null)"
echo "$err" | grep -q "rebuilding the BWT" || fail "(e) corrupt .fp did not force a rebuild"

# (f) The low-memory (-l) path rebuilds the .mbw from scratch, so it must refresh
# the fingerprint and drop stale sidecars too. Otherwise a -l rebuild from refB
# leaves the refA .fp in place, and a later (non -l) refA run sees a matching
# fingerprint and reuses refB's BWT -- the silent wrong-coordinate bug. -l forbids
# multi-density, so the sa.u4 sidecar is seeded by a non -l build first.
if "$MINIBWA" index -l "$refA" "$WORK/lm" >/dev/null 2>&1; then
	echo "# (f) -l (low-memory) rebuild must refresh the fingerprint and drop stale sidecars"
	rm -f "$pfx".* ; "$MINIBWA" index -u 3,4 "$refA" "$pfx" >/dev/null 2>&1
	[ -f "$pfx.sa.u4" ] || fail "(f) setup: expected sa.u4 sidecar from -u 3,4"
	"$MINIBWA" index -l "$refB" "$pfx" >/dev/null 2>&1
	[ -f "$pfx.mbw.fp" ] || fail "(f) -l rebuild wrote no .mbw.fp"
	[ ! -f "$pfx.sa.u4" ] || fail "(f) -l rebuild left a stale sa.u4 sidecar"
	md5_B="$(mbw_md5 "$pfx.mbw")"
	err="$("$MINIBWA" index "$refA" "$pfx" 2>&1 >/dev/null)"
	echo "$err" | grep -q "doesn't match the current reference" || fail "(f) stale fingerprint after -l let refB's .mbw be reused for refA"
	[ "$(mbw_md5 "$pfx.mbw")" != "$md5_B" ] || fail "(f) .mbw not rebuilt for refA after -l refB"
	rm -f "$WORK/lm".*
else
	echo "# (f) skipped: this minibwa was built without -l (USE_GPL) support"
fi

# (g) A fresh build must drop stale sidecars even when no .mbw is present. A deleted
# .mbw (or a partially-copied index dir) can leave a .sa.u<N> from a previous
# reference; without removal it would be attached at map time with wrong coordinates.
echo "# (g) fresh build with a deleted .mbw must still drop stale sidecars"
rm -f "$pfx".* ; "$MINIBWA" index -u 3,4 "$refA" "$pfx" >/dev/null 2>&1
[ -f "$pfx.sa.u4" ] || fail "(g) setup: expected sa.u4 sidecar from -u 3,4"
rm -f "$pfx.mbw" "$pfx.mbw.fp"                       # .mbw gone, stale sa.u4 (refA) remains
"$MINIBWA" index "$refB" "$pfx" >/dev/null 2>&1      # fresh build (no .mbw) from a DIFFERENT reference
[ ! -f "$pfx.sa.u4" ] || fail "(g) stale sa.u4 (from refA) survived a fresh build with no .mbw"

# (h) A fingerprint match with an UNREADABLE .mbw falls through to a rebuild
# (mb_bwt_load_nosa returns NULL). That fallback is a fresh build too, so it must
# drop stale sidecars first -- otherwise a leftover .sa.u<N> is later attached to
# the rebuilt BWT with wrong coordinates. The pre-load removal is gated on
# !reused, so the reused-but-load-failed path needs its own removal.
echo "# (h) reused-but-unreadable .mbw must drop stale sidecars on the rebuild fallback"
rm -f "$pfx".* ; "$MINIBWA" index -u 3,6 "$refA" "$pfx" >/dev/null 2>&1
[ -f "$pfx.sa.u6" ] || fail "(h) setup: expected sa.u6 sidecar from -u 3,6"
printf 'GARBAGE-not-a-valid-BWT' > "$pfx.mbw"        # corrupt .mbw; keep the matching .mbw.fp
"$MINIBWA" index -u 3 "$refA" "$pfx" >/dev/null 2>&1  # same ref: fp matches -> reused, load fails -> rebuild
[ ! -f "$pfx.sa.u6" ] || fail "(h) stale sa.u6 survived a rebuild triggered by an unreadable reused .mbw"

echo "PASS: index fingerprint guard"
