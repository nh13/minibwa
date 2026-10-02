#!/bin/sh
# Regression test for the signed `1 << sa_bit` overflow in regime naming and
# --list-regimes (regime.c: fill_bwt_regime, mb_regime_list_print). A valid but
# extreme density (-u 31, accepted by the index side's [0,32) contract) makes the
# 32-bit signed shift 1<<31 overflow into the sign bit -> a negative regime name
# ("sa-2147483648") and rate ("1/-2147483648"). The unsigned shift prints the
# correct positive value. Generates all FASTA inline (no committed test data).
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

echo "# (a) an extreme but valid density (-u 31) must name the regime with an unsigned shift"
"$MINIBWA" index -u 31 "$ref" "$WORK/x31" >/dev/null 2>&1
out="$("$MINIBWA" map --list-regimes "$WORK/x31" 2>/dev/null)"
echo "$out" | grep -q -- '-2147483648' && fail "(a) signed 1<<31 overflow: negative density in --list-regimes"
echo "$out" | grep -q 'sa2147483648' || fail "(a) expected regime name sa2147483648"
echo "$out" | grep -q '1/2147483648' || fail "(a) expected sampling rate 1/2147483648"

echo "# (b) a normal density (-u 4) still prints sa16 / 1/16"
"$MINIBWA" index -u 4 "$ref" "$WORK/x4" >/dev/null 2>&1
out="$("$MINIBWA" map --list-regimes "$WORK/x4" 2>/dev/null)"
echo "$out" | grep -q 'sa16 ' || fail "(b) expected regime name sa16"
echo "$out" | grep -q '1/16 ' || fail "(b) expected sampling rate 1/16"

echo "PASS: regime sa_bit shift guard"
