#!/bin/sh
# Integration test for Task 8: PE liftover-group reconciliation integration.
#
# Runs `minibwa mem` in PAIRED-END mode on synthetic fixtures and asserts on the
# SAM flags / MAPQ that the three PE hooks produce.  WGS is always PE, so this is
# what makes the ALT liftover-group feature actually work on the real target.
#
# Cases:
#  (HAPPY) Hook B: a pair whose BOTH mates lie inside a primary region with an
#      identical ALT twin contig.  On primary the pair is best; on the ALT twin an
#      equally-good second PAIR forms (pairs form only within one contig).  Without
#      Hook B that twin pair depresses paux.sub_sc => mapq_pe ~= 0 and the ALT twin
#      can even win sam_pri.  With Hook B (exclude non-rep group members from pair
#      enumeration) the twin pair never forms: chrP is the proper-pair primary for
#      BOTH mates with MAPQ>0, the chrP_altF twin records are secondary (0x100).
#      [RED proven: pre-Hook-B mapq_pe collapses and chrP_altF wins the pair.]
#  (PARALOG) Hook C: a pair where one mate is an ambiguous paralog (COPY1/COPY2 on
#      chrP whose ALT twins lift just outside MB_LIFT_TOL).  The chimeric demotion
#      loop (pe.c:524-535) would, by query overlap alone, fold the COPY2 hit under
#      the chosen COPY1 pair-hit -- merging two DISTINCT lifted groups and undoing
#      the Task-4 paralog guard.  Hook C skips the demotion across distinct lifted
#      groups, so COPY2 stays its own group representative: it appears as a
#      supplementary (0x800) record, NOT a demoted secondary (0x100).
#      [RED proven: pre-Hook-C COPY2 is 0x100 (merged).]
#  (RESCUE) is_alt stamp + Hook A: one mate is peppered so it cannot seed and is
#      recovered only by mate rescue, landing inside an ALT-twinned window.  The
#      rescued ALT hit must carry is_alt (mb_matesw_align memset clears it) so the
#      reconciliation groups it.  Asserts: the unseedable mate is rescue-only
#      (unmapped single-end), the rescued chrP_altG hit is grouped (secondary,
#      0x100) rather than a second proper-pair primary, chrP is the proper-pair
#      primary with MAPQ>0.
#  (BASELINE) chrM PE with NO .alt: every hook is gated (mb_any_alt false / no ALT
#      hits), so output is byte-identical across runs.
#
# Usage: test/altlg/test-pe.sh [<minibwa-dir>]
set -eu

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MINIBWA="$MDIR/minibwa"

# Shared no-.alt baseline check (see lib-baseline.sh).
. "$(dirname "$0")/lib-baseline.sh"
MK_HAPPY="$MDIR/test/altlg/mkfixture-pe-happy.sh"
MK_PARA="$MDIR/test/altlg/mkfixture-pe-paralog.sh"
MK_RESC="$MDIR/test/altlg/mkfixture-pe-rescue.sh"

TMPD=$(mktemp -d /tmp/altlg-pe.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

fail() { echo "FAIL: $1"; exit 1; }
ok()   { echo "  ok: $1"; }

# has_bit <flag> <bit> -> "1" if (flag & bit) else "0" (single power-of-two bit).
has_bit() { mawk -v f="$1" -v b="$2" 'BEGIN{ f=int(f); b=int(b);
    printf "%d\n", (int(f/b)%2==1)?1:0; }'; }
# flag_of <sam> <qname> <rname> <pos> -> FLAG of the matching record (pos optional).
flag_of() { mawk -v q="$2" -v c="$3" -v p="${4:-}" '$1==q && $3==c && (p=="" || $4==p){print $2; exit}' "$1"; }
mapq_of() { mawk -v q="$2" -v c="$3" -v p="${4:-}" '$1==q && $3==c && (p=="" || $4==p){print $5; exit}' "$1"; }
# proper-pair primary FLAG for a given mate-bit (0x40 first / 0x80 last).
prim_flag_mate() { mawk -v q="$2" -v mb="$3" '$1==q { f=int($2);
    if (int(f/256)%2==0 && int(f/2048)%2==0 && int(f/mb)%2==1) { print $2; exit } }' "$1"; }
prim_mapq_mate() { mawk -v q="$2" -v mb="$3" '$1==q { f=int($2);
    if (int(f/256)%2==0 && int(f/2048)%2==0 && int(f/mb)%2==1) { print $5; exit } }' "$1"; }

run_pe() { # <fixture-dir> <out.sam>
    "$MINIBWA" index "$1/ref.fa" 2>/dev/null
    "$MINIBWA" mem --outn=50 "$1/ref.fa" "$1/reads_1.fq" "$1/reads_2.fq" 2>/dev/null \
        | mawk '$1 !~ /^@/' > "$2"
    [ -s "$2" ] || fail "mapping produced no alignments ($1)"
}

# =========================================================================
echo "[test-pe] (HAPPY / Hook B) building fixture ..."
HD="$TMPD/happy"; /bin/sh "$MK_HAPPY" "$HD" 2>/dev/null
run_pe "$HD" "$TMPD/happy.sam"
echo "----- happy r-twin SAM -----"
mawk '$1=="r-twin"{printf "%s flag=%d %s pos=%s mapq=%s\n",$1,$2,$3,$4,$5}' "$TMPD/happy.sam"
echo "----------------------------"

# chrP must be the proper-pair primary for BOTH mates with MAPQ>0.
for mb in 64 128; do  # 0x40 first mate, 0x80 last mate
    pf=$(prim_flag_mate "$TMPD/happy.sam" r-twin "$mb")
    pm=$(prim_mapq_mate "$TMPD/happy.sam" r-twin "$mb")
    [ -n "$pf" ] || fail "(HAPPY) no proper-pair primary for mate bit $mb"
    # primary must be on chrP (not the ALT twin) and a proper pair (0x2)
    rn=$(mawk -v f="$pf" '$1=="r-twin" && $2==f{print $3; exit}' "$TMPD/happy.sam")
    [ "$rn" = "chrP" ] || fail "(HAPPY) mate $mb primary on $rn, expected chrP (ALT twin won the pair)"
    [ "$(has_bit "$pf" 2)" = "1" ] || fail "(HAPPY) mate $mb primary not a proper pair (flag=$pf)"
    [ "$pm" -gt 0 ] || fail "(HAPPY) mate $mb MAPQ=$pm; twin depressed mapq_pe (expected >0)"
    ok "(HAPPY) mate $mb: chrP proper-pair primary, MAPQ=$pm > 0"
done
# every chrP_altF twin record must be secondary (0x100).  Run the assertion in a
# single mawk that exits nonzero on the first non-secondary record; a `while read`
# pipeline would `exit 1` only from its subshell under /bin/sh, falsely PASSing.
if ! mawk '$1=="r-twin" && $3=="chrP_altF" && int($2/256)%2==0 { bad=1; exit } END{ exit (bad?1:0) }' "$TMPD/happy.sam"; then
    fail "(HAPPY) found chrP_altF record that is not secondary (0x100)"
fi
ok "(HAPPY) all chrP_altF twin records are secondary (0x100)"

# =========================================================================
echo "[test-pe] (PARALOG / Hook C) building fixture ..."
PD="$TMPD/para"; /bin/sh "$MK_PARA" "$PD" 2>/dev/null
run_pe "$PD" "$TMPD/para.sam"
echo "----- paralog r-para SAM -----"
mawk '$1=="r-para"{printf "%s flag=%d %s pos=%s mapq=%s\n",$1,$2,$3,$4,$5}' "$TMPD/para.sam"
echo "------------------------------"
# COPY1 is the chosen pair-hit (first mate, the unique R2 anchors it near COPY1).
# COPY2 lives at chrP pos 2371 (0-based 2370).  Hook C must keep it a DISTINCT
# lifted group: it must NOT be demoted to 0x100 secondary (which is what the
# cross-group merge would do).  A separate group rep that is not sam_pri shows as
# 0x800 supplementary.
c2_flag=$(flag_of "$TMPD/para.sam" r-para chrP 2371)
[ -n "$c2_flag" ] || fail "(PARALOG) no chrP record at pos 2371 (COPY2)"
c2_sec=$(has_bit "$c2_flag" 256)
[ "$c2_sec" = "0" ] || fail "(PARALOG) COPY2 (pos 2371) is 0x100 secondary (flag=$c2_flag): merged into COPY1's group (Hook C failed)"
# Pin the CURRENT state: a distinct group rep that is not sam_pri shows as 0x800
# supplementary.  NOTE: whether distinct paralog loci should be 0x800 vs 0x100
# (per bwa-postalt) is a DEFERRED open question for Task 6; this assertion locks
# in today's behavior so any future change here is caught and reviewed.
[ "$(has_bit "$c2_flag" 2048)" = "1" ] || fail "(PARALOG) COPY2 expected supplementary (0x800), got flag=$c2_flag"
ok "(PARALOG) COPY2 (pos 2371, flag=$c2_flag) NOT demoted to secondary: distinct lifted group preserved (0x800 supplementary)"

# =========================================================================
echo "[test-pe] (RESCUE / is_alt stamp + Hook A) building fixture ..."
RD="$TMPD/resc"; /bin/sh "$MK_RESC" "$RD" 2>/dev/null
"$MINIBWA" index "$RD/ref.fa" 2>/dev/null
# (1) the peppered mate is rescue-only: unmapped when mapped single-end.
mawk '/^@r-rescue/{p=4} p>0{print;p--}' "$RD/reads_2.fq" > "$RD/r2.fq"
r2se=$("$MINIBWA" mem --outn=50 "$RD/ref.fa" "$RD/r2.fq" 2>/dev/null | mawk '$1=="r-rescue/2"{print $2; exit}')
[ -n "$r2se" ] || fail "(RESCUE) no single-end record for R2"
[ "$(has_bit "$r2se" 4)" = "1" ] || fail "(RESCUE) R2 mapped single-end (flag=$r2se); fixture must be rescue-only"
ok "(RESCUE) R2 is rescue-only (unmapped single-end, flag=$r2se)"
# (2) full PE run.
run_pe "$RD" "$TMPD/resc.sam"
echo "----- rescue r-rescue SAM -----"
mawk '$1=="r-rescue"{printf "%s flag=%d %s pos=%s mapq=%s\n",$1,$2,$3,$4,$5}' "$TMPD/resc.sam"
echo "-------------------------------"
# chrP must be the proper-pair primary for the rescued (last) mate with MAPQ>0.
rp_flag=$(prim_flag_mate "$TMPD/resc.sam" r-rescue 128)
rp_mapq=$(prim_mapq_mate "$TMPD/resc.sam" r-rescue 128)
[ -n "$rp_flag" ] || fail "(RESCUE) no proper-pair primary for the rescued mate"
rp_rn=$(mawk -v f="$rp_flag" '$1=="r-rescue" && $2==f{print $3; exit}' "$TMPD/resc.sam")
[ "$rp_rn" = "chrP" ] || fail "(RESCUE) rescued-mate primary on $rp_rn, expected chrP"
[ "$rp_mapq" -gt 0 ] || fail "(RESCUE) rescued-mate MAPQ=$rp_mapq; expected >0"
ok "(RESCUE) rescued mate: chrP proper-pair primary, MAPQ=$rp_mapq > 0"
# the rescued ALT hit must be GROUPED (secondary, 0x100) -- not a 2nd proper-pair
# primary on the ALT contig.
rg_flag=$(mawk '$1=="r-rescue" && $3=="chrP_altG" && int($2/128)%2==1{print $2; exit}' "$TMPD/resc.sam")
[ -n "$rg_flag" ] || fail "(RESCUE) rescued ALT hit chrP_altG absent for the rescued mate"
[ "$(has_bit "$rg_flag" 256)" = "1" ] || fail "(RESCUE) rescued chrP_altG flag=$rg_flag not secondary; ALT mate not grouped"
[ "$(has_bit "$rg_flag" 2)" = "0" ] || fail "(RESCUE) rescued chrP_altG flag=$rg_flag is a proper pair; ALT mate competed"
ok "(RESCUE) rescued chrP_altG hit is grouped (secondary 0x100, not a competing pair)"

# =========================================================================
echo "[test-pe] (BASELINE) chrM PE, no .alt: byte-identical to stock ..."
chrm_baseline "(BASELINE) chrM PE" pe --outn=5

echo "[test-pe] PASS"
exit 0
