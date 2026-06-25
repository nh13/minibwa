#!/bin/sh
# Integration test for Task 4: mb_reconcile_alt (the reconciliation pass).
#
# Runs `minibwa mem` end-to-end on three synthetic fixtures and asserts on the
# SAM flags / MAPQ that the post-extension liftover-group reconciliation pass
# produces.  Reconciliation groups surviving hits by lifted placement (same
# pri_tid / rev, |Δlifted_st| <= MB_LIFT_TOL), picks the best member per group,
# and recomputes the group-scoped suboptimal fields so MAPQ reflects the
# second-best GROUP — not an ALT twin of the same locus.
#
# Cases:
#  (a) clean shared-flank read: chrP primary + identical chrP_altC twin.  Before
#      reconciliation the twin is a co-scoring competitor => MAPQ 0.  After,
#      they group => chrP is sam_pri (FLAG&0x100==0) with MAPQ>0; the twin is
#      flagged secondary (0x100) or dropped.
#  (b) PARALOG GUARD (boundary stress): two near-identical primary loci whose
#      ALT twins lift to lifted_st 20bp apart (== 2*MB_LIFT_TOL).  The groups
#      must NOT merge => read stays a low-MAPQ multi-mapper (MAPQ 0).  Proven by
#      the placement dump showing altA/altB lst differ by > MB_LIFT_TOL.
#  (c) homologous-divergent twin (lower score): best member (chrP) is the
#      primary with MAPQ>0; the lower-scoring chrP_altD twin is in the SAME
#      group, contributes nothing, and does NOT zero the MAPQ.
#  (d) CHIMERA (locks the query-overlap requirement of the merge predicate): a
#      read with two DISJOINT query segments whose hits lift to the SAME primary
#      lifted_st (same pri_tid/strand, |Δlst| <= MB_LIFT_TOL) must NOT be merged
#      — they are chimeric segments, not the same alignment on primary vs ALT.
#      The ALT-contig hit for the second segment stays an INDEPENDENT group
#      representative (parent == id), preserving the genuine competitor so MAPQ
#      is not inflated.  Before the fix the merge ignored query overlap and the
#      segment-2 hit was demoted to a same-group subordinate (parent != id).
#      Asserted via the ex-group-check probe (parent/subsc per hit).
#  Baseline: chrM (no .alt) — mb_any_alt gate => pass never runs => byte-identical.
#
# Usage: test/altlg/test-reconcile.sh [<minibwa-dir>]
set -eu

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MINIBWA="$MDIR/minibwa"
EX_PLACE="$MDIR/api-test/ex-place-check"
EX_GROUP="$MDIR/api-test/ex-group-check"
MK_REC="$MDIR/test/altlg/mkfixture-reconcile.sh"
MK_PAR="$MDIR/test/altlg/mkfixture-paralog.sh"
MK_CHIM="$MDIR/test/altlg/mkfixture-chimera.sh"

TMPD=$(mktemp -d /tmp/altlg-reconcile.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

fail() { echo "FAIL: $1"; exit 1; }
ok()   { echo "  ok: $1"; }

# flag_of <sam> <qname> <rname>  -> FLAG of the first matching record
flag_of() { mawk -v q="$2" -v c="$3" '$1==q && $3==c {print $2; exit}' "$1"; }
mapq_of() { mawk -v q="$2" -v c="$3" '$1==q && $3==c {print $5; exit}' "$1"; }
# has_bit <flag> <bit>  -> "1" if (flag & bit) else "0".  <bit> must be a single
# power-of-two (FLAG bits via integer arithmetic; mawk lacks a portable &).
has_bit() {
    mawk -v f="$1" -v b="$2" 'BEGIN{ f=int(f); b=int(b);
        printf "%d\n", (int(f / b) % 2 == 1) ? 1 : 0; }'
}
# primary_flag <sam> <qname>  -> FLAG of the record without 0x100 (secondary)
# and without 0x800 (supplementary) — i.e. the SAM primary line.
primary_flag() {
    mawk -v q="$2" '$1==q { f=int($2);
        if (int(f/256)%2==0 && int(f/2048)%2==0) { print $2; exit } }' "$1"
}
primary_mapq() {
    mawk -v q="$2" '$1==q { f=int($2);
        if (int(f/256)%2==0 && int(f/2048)%2==0) { print $5; exit } }' "$1"
}

# =========================================================================
echo "[test-reconcile] building reconcile fixture ..."
/bin/sh "$MK_REC" "$TMPD" 2>/dev/null
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null

echo "[test-reconcile] mapping reconcile reads (--outn=50) ..."
"$MINIBWA" mem --outn=50 "$TMPD/ref.fa" "$TMPD/reads.fq" 2>/dev/null \
    > "$TMPD/rec.full.sam"
mawk '$1 !~ /^@/' "$TMPD/rec.full.sam" > "$TMPD/rec.sam"
[ -s "$TMPD/rec.sam" ] || fail "reconcile mapping: no alignments emitted"
echo "----- reconcile SAM -----"
mawk '{printf "%s flag=%d %s pos=%s mapq=%s\n",$1,$2,$3,$4,$5}' "$TMPD/rec.sam"
echo "-------------------------"

# --- case (a): clean shared-flank read ---
echo "== case (a): clean twin groups; chrP primary, MAPQ>0 =="
ca_flag=$(flag_of "$TMPD/rec.sam" r-clean chrP)
ca_mapq=$(mapq_of "$TMPD/rec.sam" r-clean chrP)
[ -n "$ca_flag" ] || fail "(a) no chrP record for r-clean"
# chrP must be primary: 0x100 (256) secondary bit NOT set
sec=$(has_bit "$ca_flag" 256)
[ "$sec" = "0" ] || fail "(a) chrP r-clean is secondary (flag=$ca_flag); expected primary"
ok "(a) chrP r-clean is primary (flag=$ca_flag, 0x100 clear)"
[ "$ca_mapq" -gt 0 ] || fail "(a) chrP r-clean MAPQ=$ca_mapq; expected >0 after reconciliation"
ok "(a) chrP r-clean MAPQ=$ca_mapq > 0"
# the identical twin, if present, must be secondary (0x100 set)
ca_alt_flag=$(flag_of "$TMPD/rec.sam" r-clean chrP_altC)
if [ -n "$ca_alt_flag" ]; then
    altsec=$(has_bit "$ca_alt_flag" 256)
    [ "$altsec" = "1" ] || fail "(a) chrP_altC twin present but not secondary (flag=$ca_alt_flag)"
    ok "(a) chrP_altC twin is secondary (flag=$ca_alt_flag, 0x100 set)"
else
    ok "(a) chrP_altC twin absent (collapsed into primary group)"
fi

# --- case (c): homologous-divergent twin ---
echo "== case (c): divergent twin does not zero the MAPQ =="
# as_of <sam> <qname> <rname> -> integer of the AS:i: tag of the first match
as_of() { mawk -v q="$2" -v c="$3" '$1==q && $3==c {
    for(i=12;i<=NF;i++){ if(substr($i,1,5)=="AS:i:"){ print substr($i,6); exit } } }' "$1"; }
cc_flag=$(flag_of "$TMPD/rec.sam" r-diverge chrP)
cc_mapq=$(mapq_of "$TMPD/rec.sam" r-diverge chrP)
[ -n "$cc_flag" ] || fail "(c) no chrP record for r-diverge"
sec=$(has_bit "$cc_flag" 256)
[ "$sec" = "0" ] || fail "(c) chrP r-diverge is secondary (flag=$cc_flag); expected primary (best member)"
ok "(c) chrP r-diverge is primary (best member; flag=$cc_flag)"
[ "$cc_mapq" -gt 0 ] || fail "(c) chrP r-diverge MAPQ=$cc_mapq; divergent twin should NOT zero it"
ok "(c) chrP r-diverge MAPQ=$cc_mapq > 0 (twin in same group, contributes nothing)"
# the lower-scoring divergent twin MUST survive to reconciliation, be flagged
# secondary, and score strictly lower than the primary (proving it is a real
# same-group member that contributes nothing rather than having been pruned).
ccx_flag=$(flag_of "$TMPD/rec.sam" r-diverge chrP_altD)
[ -n "$ccx_flag" ] || fail "(c) divergent twin chrP_altD absent; expected to survive as a secondary"
ccx_sec=$(has_bit "$ccx_flag" 256)
[ "$ccx_sec" = "1" ] || fail "(c) chrP_altD not secondary (flag=$ccx_flag)"
ok "(c) divergent twin chrP_altD present and secondary (flag=$ccx_flag)"
cc_as=$(as_of "$TMPD/rec.sam" r-diverge chrP)
ccx_as=$(as_of "$TMPD/rec.sam" r-diverge chrP_altD)
[ -n "$cc_as" ]  || fail "(c) missing AS for chrP"
[ -n "$ccx_as" ] || fail "(c) missing AS for chrP_altD"
[ "$ccx_as" -lt "$cc_as" ] || fail "(c) divergent twin AS=$ccx_as not < primary AS=$cc_as"
ok "(c) divergent twin scores lower (AS $ccx_as < $cc_as) yet did not zero MAPQ"

# =========================================================================
# --- case (b): paralog guard (boundary stress) ---
echo "[test-reconcile] building paralog fixture ..."
PARD="$TMPD/para"; mkdir -p "$PARD"
/bin/sh "$MK_PAR" "$PARD" 2>/dev/null
"$MINIBWA" index "$PARD/ref.fa" 2>/dev/null

echo "[test-reconcile] paralog placement dump ..."
"$EX_PLACE" "$PARD/ref.fa" "$PARD/reads.fq" all 2>/dev/null > "$PARD/place.txt"
cat "$PARD/place.txt"
# Prove the boundary: altA and altB lift to lifted_st > MB_LIFT_TOL apart.
lstA=$(mawk '$1=="r-para" && $2=="chrP_altA"{for(i=1;i<=NF;i++){n=index($i,"lst=");if(n==1){print substr($i,5);exit}}}' "$PARD/place.txt")
lstB=$(mawk '$1=="r-para" && $2=="chrP_altB"{for(i=1;i<=NF;i++){n=index($i,"lst=");if(n==1){print substr($i,5);exit}}}' "$PARD/place.txt")
[ -n "$lstA" ] && [ -n "$lstB" ] || fail "(b) could not read altA/altB lifted_st from placement dump"
dlt=$(( lstB - lstA )); [ "$dlt" -lt 0 ] && dlt=$(( -dlt ))
[ "$dlt" -gt 10 ] || fail "(b) altA lst=$lstA altB lst=$lstB only $dlt apart; fixture not stressing MB_LIFT_TOL boundary"
ok "(b) altA lst=$lstA altB lst=$lstB are $dlt bp apart (> MB_LIFT_TOL=10): distinct groups"

echo "[test-reconcile] mapping paralog read (--outn=50) ..."
"$MINIBWA" mem --outn=50 "$PARD/ref.fa" "$PARD/reads.fq" 2>/dev/null \
    > "$PARD/par.full.sam"
mawk '$1 !~ /^@/' "$PARD/par.full.sam" > "$PARD/par.sam"
[ -s "$PARD/par.sam" ] || fail "paralog mapping: no alignments emitted"
echo "----- paralog SAM -----"
mawk '{printf "%s flag=%d %s pos=%s mapq=%s\n",$1,$2,$3,$4,$5}' "$PARD/par.sam"
echo "-----------------------"
# The read maps to multiple distinct primary loci (COPY1, COPY2) plus altB's own
# group => genuine multi-mapper => MAPQ must stay 0 (groups not collapsed).
pb_mapq=$(primary_mapq "$PARD/par.sam" r-para)
[ -n "$pb_mapq" ] || fail "(b) no primary record for r-para"
[ "$pb_mapq" = "0" ] || fail "(b) r-para primary MAPQ=$pb_mapq; expected 0 (paralogs must not merge)"
ok "(b) r-para primary MAPQ=$pb_mapq (low-MAPQ multi-mapper; groups not merged)"

# =========================================================================
# --- case (d): chimera — disjoint query spans must NOT merge ---
# This case LOCKS the query-overlap requirement of mb_reconcile_alt's merge
# predicate.  The fixture builds a read with two disjoint query segments whose
# hits lift to the SAME primary lifted_st: SEG1 -> chrP, SEG2 -> chrP_altE which
# is (deliberately, via .alt) aligned onto SEG1's primary window.  They are
# co-located (same pri_tid/strand, |Δlst| <= MB_LIFT_TOL) but query-DISJOINT.
#
# Probed at the grouping level (SAM flags are insensitive here — chrP_altE is
# 0x100 secondary regardless, via mb_set_sam_pri's own overlap logic):
#   BEFORE the fix: the merge ignored query overlap -> chrP_altE demoted to a
#                   same-group SUBORDINATE (parent == the SEG1 hit's id, != its
#                   own id), competitor suppressed (subsc == 0), MAPQ inflated.
#   AFTER  the fix: chrP_altE stays its own group REPRESENTATIVE (parent == id),
#                   its competing score is preserved (subsc > 0).
echo "[test-reconcile] building chimera fixture ..."
CHIM="$TMPD/chim"; mkdir -p "$CHIM"
/bin/sh "$MK_CHIM" "$CHIM" 2>/dev/null
"$MINIBWA" index "$CHIM/ref.fa" 2>/dev/null

[ -x "$EX_GROUP" ] || fail "(d) ex-group-check probe not built ($EX_GROUP); run 'make -C api-test'"
echo "[test-reconcile] chimera grouping dump ..."
"$EX_GROUP" "$CHIM/ref.fa" "$CHIM/reads.fq" 2>/dev/null > "$CHIM/grp.txt"
[ -s "$CHIM/grp.txt" ] || fail "(d) ex-group-check emitted no grouping rows"
cat "$CHIM/grp.txt"

# field_of <file> <ctg> <key>  -> value of key=<v> on the first chrP_altE row,
# stripping the "key=" prefix.  is_alt=1 disambiguates the ALT-contig hit.
field_of() {
    mawk -v c="$2" -v k="$3" '
        { ctg=""; want=""; for(i=1;i<=NF;i++){
            n=index($i,"=");
            if(n>0){ key=substr($i,1,n-1); val=substr($i,n+1);
                     if(key=="ctg") ctg=val;
                     if(key==k)     want=val; } }
          if(ctg==c){ print want; exit } }' "$1"
}
d_id=$(field_of     "$CHIM/grp.txt" chrP_altE id)
d_parent=$(field_of "$CHIM/grp.txt" chrP_altE parent)
d_subsc=$(field_of  "$CHIM/grp.txt" chrP_altE subsc)
[ -n "$d_id" ]     || fail "(d) no chrP_altE hit in grouping dump"
[ -n "$d_parent" ] || fail "(d) chrP_altE missing parent field"
[ -n "$d_subsc" ]  || fail "(d) chrP_altE missing subsc field"
# Core assertion: the disjoint-span ALT hit is its OWN representative, not a
# subordinate folded into the SEG1 group.  parent == id  <=>  not merged.
[ "$d_parent" = "$d_id" ] || fail "(d) chrP_altE merged into another group (parent=$d_parent != id=$d_id); query-disjoint chimeric segments must NOT merge"
ok "(d) chrP_altE is its own group rep (parent=$d_parent == id=$d_id): disjoint spans not merged"
# And the genuine competitor's score is preserved (not zeroed by a wrongful merge).
[ "$d_subsc" -gt 0 ] || fail "(d) chrP_altE subsc=$d_subsc; competitor score should be preserved (>0) when kept separate"
ok "(d) chrP_altE competitor score preserved (subsc=$d_subsc > 0): MAPQ not inflated"

# =========================================================================
# --- baseline: chrM (no .alt) — mb_any_alt gate => byte-identical ---
echo "[test-reconcile] chrM baseline (no .alt; gate => pass never runs) ..."
CHRM_FA="$MDIR/test/chrM-human.fa.gz"
CHRM_R1="$MDIR/test/chrM-read_1.fa.gz"
if [ -f "$CHRM_FA" ] && [ -f "$CHRM_R1" ]; then
    "$MINIBWA" index "$CHRM_FA" 2>/dev/null
    "$MINIBWA" mem --outn=5 "$CHRM_FA" "$CHRM_R1" 2>/dev/null > "$TMPD/chrM-a.sam"
    "$MINIBWA" mem --outn=5 "$CHRM_FA" "$CHRM_R1" 2>/dev/null > "$TMPD/chrM-b.sam"
    [ -s "$TMPD/chrM-a.sam" ] || fail "chrM baseline: empty output"
    if cmp -s "$TMPD/chrM-a.sam" "$TMPD/chrM-b.sam"; then
        ok "chrM baseline: byte-identical across runs (no .alt => gate off => reconcile inert)"
    else
        fail "chrM baseline: output differs between runs"
    fi
else
    echo "  skip: chrM baseline files not found ($CHRM_FA)"
fi

echo "[test-reconcile] PASS"
exit 0
