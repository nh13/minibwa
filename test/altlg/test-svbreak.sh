#!/bin/sh
# Integration test: SV-breakpoint-aware ALT grouping (multi-interval mb_place_t).
#
# THE BUG (verified on real WGS): an ALT contig whose .alt CIGAR carries an
# SV-scale indel describes two M blocks lifting to primary positions thousands of
# bp apart.  A read whose ALT footprint spans the breakpoint collapsed (under the
# old single-interval placement) to a lifted_st thousands of bp from the read's
# PRIMARY twin -> the ALT hit failed to co-locate with the twin -> they competed
# -> a FALSE MAPQ 0 (bwa-postalt keeps these confident).
#
# THE FIX (Option 1, multi-interval): mb_hit_place records one sub-placement per
# overlapping .alt lift block; two hits co-locate iff ANY sub-interval pair shares
# (pri_tid, rev, |Δst| <= MB_LIFT_TOL).  The breakpoint-spanning ALT hit then keeps
# a sub-interval AT the twin's position -> they group -> MAPQ recovered.
#
# Cases:
#  (fwd) forward .alt block: read maps forward to its primary twin.  RED before
#        the fix (primary MAPQ 0); GREEN after (primary MAPQ>0, ALT secondary).
#  (rev) reverse .alt block (FLAG 0x10): read maps reverse to its primary twin.
#        Exercises the reverse-fold lift paths.  Same RED->GREEN assertion.
#  (par) PARALOG SAFETY: two distinct primary loci 7kb apart (same contig/strand,
#        overlapping query) + an ALT twin co-located with ONE of them.  The OTHER
#        primary stays a distinct competing group -> MAPQ MUST remain 0.  Proves
#        the multi-interval change does NOT merge two far-apart primary loci.
#  Baseline: chrM (no .alt) — mb_any_alt gate => pass never runs => byte-identical.
#
# Usage: test/altlg/test-svbreak.sh [<minibwa-dir>]
set -eu

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MINIBWA="$MDIR/minibwa"

# Shared no-.alt baseline check (see lib-baseline.sh).
. "$(dirname "$0")/lib-baseline.sh"
MK_SV="$MDIR/test/altlg/mkfixture-svbreak.sh"

TMPD=$(mktemp -d /tmp/altlg-svbreak.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

fail() { echo "FAIL: $1"; exit 1; }
ok()   { echo "  ok: $1"; }

# primary_flag/_mapq <sam> <qname> -> the SAM-primary record (no 0x100, no 0x800)
primary_flag() {
    mawk -v q="$2" '$1==q { f=int($2);
        if (int(f/256)%2==0 && int(f/2048)%2==0) { print $2; exit } }' "$1"
}
primary_mapq() {
    mawk -v q="$2" '$1==q { f=int($2);
        if (int(f/256)%2==0 && int(f/2048)%2==0) { print $5; exit } }' "$1"
}
# secondary present on a given contig (0x100 set)?  prints the FLAG of the first.
alt_secondary_flag() {
    mawk -v q="$2" -v c="$3" '$1==q && $3==c { f=int($2);
        if (int(f/256)%2==1) { print $2; exit } }' "$1"
}

# =========================================================================
echo "[test-svbreak] building svbreak fixture ..."
/bin/sh "$MK_SV" "$TMPD" 2>/dev/null

# ---------------- forward variant ----------------
echo "== case (fwd): breakpoint-spanning ALT groups with forward primary twin =="
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null
"$MINIBWA" map --outn=50 "$TMPD/ref.fa" "$TMPD/reads.fq" 2>/dev/null > "$TMPD/fwd.full.sam"
mawk '$1 !~ /^@/' "$TMPD/fwd.full.sam" > "$TMPD/fwd.sam"
[ -s "$TMPD/fwd.sam" ] || fail "(fwd) no alignments emitted"
echo "----- fwd SAM -----"
mawk '{printf "%s flag=%d %s pos=%s mapq=%s cig=%s\n",$1,$2,$3,$4,$5,$6}' "$TMPD/fwd.sam"
echo "-------------------"
fwd_flag=$(primary_flag "$TMPD/fwd.sam" r-sv-fwd)
fwd_mapq=$(primary_mapq "$TMPD/fwd.sam" r-sv-fwd)
[ -n "$fwd_mapq" ] || fail "(fwd) no SAM-primary record for r-sv-fwd"
# The SAM-primary must be the primary-contig twin (chrP), forward strand.
fwd_rname=$(mawk -v q=r-sv-fwd '$1==q { f=int($2);
    if (int(f/256)%2==0 && int(f/2048)%2==0) { print $3; exit } }' "$TMPD/fwd.sam")
[ "$fwd_rname" = "chrP" ] || fail "(fwd) SAM-primary on $fwd_rname, expected chrP twin"
ok "(fwd) SAM-primary is chrP twin (flag=$fwd_flag)"
[ "$fwd_mapq" -gt 0 ] || fail "(fwd) chrP twin MAPQ=$fwd_mapq; expected >0 (breakpoint ALT must group, not compete)"
ok "(fwd) chrP twin MAPQ=$fwd_mapq > 0 (grouped with breakpoint-spanning ALT)"
fwd_alt=$(alt_secondary_flag "$TMPD/fwd.sam" r-sv-fwd chrP_altV)
[ -n "$fwd_alt" ] || fail "(fwd) chrP_altV not present as secondary; expected grouped/demoted"
ok "(fwd) chrP_altV ALT twin demoted to secondary (flag=$fwd_alt)"

# ---------------- reverse variant ----------------
echo "== case (rev): breakpoint-spanning ALT groups with reverse primary twin =="
"$MINIBWA" index "$TMPD/rev/ref.fa" 2>/dev/null
"$MINIBWA" map --outn=50 "$TMPD/rev/ref.fa" "$TMPD/rev/reads.fq" 2>/dev/null > "$TMPD/rev.full.sam"
mawk '$1 !~ /^@/' "$TMPD/rev.full.sam" > "$TMPD/rev.sam"
[ -s "$TMPD/rev.sam" ] || fail "(rev) no alignments emitted"
echo "----- rev SAM -----"
mawk '{printf "%s flag=%d %s pos=%s mapq=%s cig=%s\n",$1,$2,$3,$4,$5,$6}' "$TMPD/rev.sam"
echo "-------------------"
rev_flag=$(primary_flag "$TMPD/rev.sam" r-sv-rev)
rev_mapq=$(primary_mapq "$TMPD/rev.sam" r-sv-rev)
[ -n "$rev_mapq" ] || fail "(rev) no SAM-primary record for r-sv-rev"
rev_rname=$(mawk -v q=r-sv-rev '$1==q { f=int($2);
    if (int(f/256)%2==0 && int(f/2048)%2==0) { print $3; exit } }' "$TMPD/rev.sam")
[ "$rev_rname" = "chrP" ] || fail "(rev) SAM-primary on $rev_rname, expected chrP twin"
ok "(rev) SAM-primary is chrP twin (flag=$rev_flag)"
# confirm the twin maps on the reverse strand (0x10 set)
rev_is_rev=$(mawk -v f="$rev_flag" 'BEGIN{ print (int(int(f)/16)%2==1)?1:0 }')
[ "$rev_is_rev" = "1" ] || fail "(rev) chrP twin not reverse-strand (flag=$rev_flag); fixture/lift mismatch"
ok "(rev) chrP twin is reverse-strand (flag=$rev_flag, 0x10 set)"
[ "$rev_mapq" -gt 0 ] || fail "(rev) chrP twin MAPQ=$rev_mapq; expected >0 (reverse breakpoint ALT must group)"
ok "(rev) chrP twin MAPQ=$rev_mapq > 0 (grouped via reverse-folded sub-interval)"
rev_alt=$(alt_secondary_flag "$TMPD/rev.sam" r-sv-rev chrP_altV)
[ -n "$rev_alt" ] || fail "(rev) chrP_altV not present as secondary; expected grouped/demoted"
ok "(rev) chrP_altV ALT twin demoted to secondary (flag=$rev_alt)"

# =========================================================================
# --- case (par): paralog safety — two primary loci 7kb apart + ALT ---
echo "== case (par): two primary loci 7kb apart + ALT must STAY MAPQ 0 =="
PARD="$TMPD/par"; mkdir -p "$PARD"
python3 - "$PARD" <<'PY'
import random, sys
d = sys.argv[1]
def gen(seed, n):
    random.seed(seed); return ''.join(random.choice('ACGT') for _ in range(n))
COPY = gen(2, 150)               # identical paralog COPY
PAD0 = gen(1, 200)
GAP  = gen(3, 7000)              # the two primary loci are 7kb apart
PADZ = gen(5, 200)
CHRP = PAD0 + COPY + GAP + COPY + PADZ
open(d + "/ref.fa", "w").write(">chrP\n%s\n>chrP_altA\n%s\n" % (CHRP, COPY))
# ALT co-locates with COPY1 (1-based POS 201).  COPY2 (7kb away) stays distinct.
open(d + "/ref.fa.alt", "w").write("chrP_altA\t0\tchrP\t201\t60\t150M\t*\t0\t0\t*\t*\n")
open(d + "/reads.fq", "w").write("@r-par7k\n%s\n+\n%s\n" % (COPY, "I" * 150))
PY
"$MINIBWA" index "$PARD/ref.fa" 2>/dev/null
"$MINIBWA" map --outn=50 "$PARD/ref.fa" "$PARD/reads.fq" 2>/dev/null > "$PARD/par.full.sam"
mawk '$1 !~ /^@/' "$PARD/par.full.sam" > "$PARD/par.sam"
[ -s "$PARD/par.sam" ] || fail "(par) no alignments emitted"
echo "----- par SAM -----"
mawk '{printf "%s flag=%d %s pos=%s mapq=%s cig=%s\n",$1,$2,$3,$4,$5,$6}' "$PARD/par.sam"
echo "-------------------"
# Sanity: the read must produce >= 2 primary-contig hits (the two paralog loci) so
# the multi-mapper verdict is genuinely about distinct primaries, not a missed hit.
n_chrP=$(mawk -v q=r-par7k '$1==q && $3=="chrP"{c++} END{print c+0}' "$PARD/par.sam")
[ "$n_chrP" -ge 2 ] || fail "(par) only $n_chrP chrP hits; fixture must expose two distinct primary loci"
ok "(par) two distinct chrP primary loci present ($n_chrP hits)"
par_mapq=$(primary_mapq "$PARD/par.sam" r-par7k)
[ -n "$par_mapq" ] || fail "(par) no SAM-primary record for r-par7k"
[ "$par_mapq" = "0" ] || fail "(par) r-par7k MAPQ=$par_mapq; expected 0 (distinct primary loci must NOT merge)"
ok "(par) r-par7k MAPQ=$par_mapq (paralogs not merged; multi-interval change is locus-safe)"

# =========================================================================
# --- baseline: chrM (no .alt) — mb_any_alt gate => byte-identical ---
echo "[test-svbreak] chrM baseline (no .alt): byte-identical to stock ..."
chrm_baseline "(BASELINE) chrM" se --outn=5

echo "[test-svbreak] PASS"
exit 0
