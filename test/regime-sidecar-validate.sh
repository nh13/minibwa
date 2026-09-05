#!/bin/sh
# Regression test for sidecar validation in mb_regime_discover (regime.c). The
# discovery step used to register a <prefix>.sa.u* sidecar after reading only its
# sa_bit field, so a malformed sidecar (wrong magic, truncated SA payload) became
# a selectable regime. mb_regime_pick() can prefer a sparser sidecar (higher
# speed_rank), and mb_bwt_load_sa() then fails on the corrupt file with no
# fallback -- map aborts. discover now mirrors the loader's checks (MB_SA_MAGIC,
# the SA count expected from seq_len, and a complete on-disk payload), so a
# corrupt sidecar is skipped at discovery and never listed or selected.
# Generates all FASTA inline (no committed test data).
set -eu

MINIBWA="${MINIBWA:-./minibwa}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail() { echo "FAIL: $*" >&2; exit 1; }

gen_fa() { # <len> <seed> <name> <out.fa>
	awk -v len="$1" -v seed="$2" -v name="$3" 'BEGIN{
		b="ACGT"; x=seed+1; printf ">%s\n", name;
		for(i=0;i<len;i++){ x=(1103515245*x+12345)%2147483648; printf "%s", substr(b,(x%4)+1,1);
			if((i+1)%60==0) printf "\n"; }
		if(len%60!=0) printf "\n";
	}' > "$4"
}

ref="$WORK/ref.fa"; gen_fa 2000 1 chrA "$ref"

# A two-density index: sa_bits[0]=4 is bundled into .mbw, sa_bits[1]=6 becomes the
# sidecar .sa.u6 (regime sa64).
build_index() { # <prefix>
	"$MINIBWA" index -u 4,6 "$ref" "$1" >/dev/null 2>&1
	[ -f "$1.sa.u6" ] || fail "setup: expected sidecar $1.sa.u6 to exist"
}

echo "# (a) an intact sidecar is discovered: both sa16 (bundled) and sa64 (sidecar) list"
build_index "$WORK/ok"
out="$("$MINIBWA" map --list-regimes "$WORK/ok" 2>/dev/null)"
echo "$out" | grep -q 'sa16 ' || fail "(a) expected bundled regime sa16"
echo "$out" | grep -q 'sa64 ' || fail "(a) expected sidecar regime sa64"

echo "# (b) a sidecar with a corrupted magic is skipped, not registered"
build_index "$WORK/badmagic"
printf 'XXXX' | dd of="$WORK/badmagic.sa.u6" bs=1 count=4 conv=notrunc >/dev/null 2>&1
out="$("$MINIBWA" map --list-regimes "$WORK/badmagic" 2>/dev/null)"
echo "$out" | grep -q 'sa16 ' || fail "(b) bundled regime sa16 must still list"
echo "$out" | grep -q 'sa64 ' && fail "(b) corrupt-magic sidecar must NOT be registered as sa64"

echo "# (c) a sidecar with a truncated SA payload is skipped, not registered"
build_index "$WORK/short"
head -c 20 "$WORK/short.sa.u6" > "$WORK/short.sa.u6.tmp" && mv "$WORK/short.sa.u6.tmp" "$WORK/short.sa.u6"
out="$("$MINIBWA" map --list-regimes "$WORK/short" 2>/dev/null)"
echo "$out" | grep -q 'sa16 ' || fail "(c) bundled regime sa16 must still list"
echo "$out" | grep -q 'sa64 ' && fail "(c) truncated sidecar must NOT be registered as sa64"

echo "PASS: regime sidecar validation"
