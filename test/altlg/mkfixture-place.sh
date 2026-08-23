#!/bin/sh
# Fixture for Task 2 (mb_hit_place — per-hit lifted placement).
#
# Reference (ref.fa):
#   chrP       primary contig (length L)
#   chrP_alt   forward ALT  = chrP with a 5bp insertion ("GGGGG") at offset 300,
#              so the .alt CIGAR is 300M5I<rest>M, FLAG 0  (POS 1).
#   chrR_alt   reverse-strand ALT = reverse complement of chrP[250,450);
#              .alt line is FLAG 16, RNAME chrP, POS 251, CIGAR 200M.
#
# Reads (one FASTA each so QNAME survives):
#   pa-fwd     chrP[120,270)  forward, inside the shared flank -> case (a)/(b):
#              primary chrP hit (identity) and forward chrP_alt hit lift to the
#              same primary lifted_st (~120), rev 0, liftable 1.
#   c-fwd      chrP[300,400)  forward      -> case (c) primary side: chrP hit, rev 0
#   c-rev-twin reverse complement of chrP[300,400): maps to chrR_alt on the minus
#              strand; block_rev(1) XOR h->rev(1) = 0, footprint folds back to
#              chrP[300,400) so lifted_st ~300, rev 0 -> groups with c-fwd. case (c)
#   hole-end   chrP_alt read that STARTS inside the 5bp "GGGGG" insertion (a hole)
#              and continues into liftable bases -> case (d): liftable 1, lifted_st
#              taken from the first liftable base.
#   twin-indel chrP read with a tiny (<=tol) deletion vs chrP_alt: still groups.
set -eu
. "$(dirname "$0")/fixlib.sh"
d="$1"; mkdir -p "$d"

S=$(printf '%s' 'GATCCTAGCATGCTAGGCTAACGTTAGCCGATCGTAGCTAGGCATCGATCGTAGCTAGCTAGGCATCGATTACGATCGGCTAATCGATCGTAGCTGATCGA'\
'TCGTAGCTAGCATCGATCGTAGCATCGGCTAGCATCGATCGATTACGCATCGATCGTAGGCTAGCATCGATCGTAGCTAGCATCGGCTAGCATCGATTACG'\
'ATCGGCTAATCGATCGTAGCTGATCGATCGTAGCTAGCATCGATCGTAGCATCGGCTAGCATCGATCGATTACGCATCGATCGTAGGCTAGCATCGATCGT'\
'AGCTAGCATCGGCTAGCATCGATTACGATCGGCTAATCGATCGTAGCTGATCGAGATCCTAGCATGCTAGGCTAACGTTAGCCGATCGTAGCTAGGCATCG'\
'ATCGTAGCTAGCTAGGCATCGATTACGATCGGCTAATCGATCGTAGCTGATCGATCGTAGCTAGCATCGATCGTAGCATCGGCTAGCATCGATCGATTACG')
L=${#S}

# revcomp helper (reads stdin, writes revcomp to stdout)

# --- contigs ---
printf '>chrP\n%s\n' "$S" > "$d/ref.fa"
ALT=$(printf '%s' "$S" | mawk '{print substr($0,1,300) "GGGGG" substr($0,301)}')   # 5bp insertion at 300
printf '>chrP_alt\n%s\n' "$ALT" >> "$d/ref.fa"
# chrR_alt = revcomp of chrP[250,450) (1-based cut 251..450)
WIN=$(printf '%s' "$S" | cut -c251-450)
RWIN=$(printf '%s' "$WIN" | rc)
printf '>chrR_alt\n%s\n' "$RWIN" >> "$d/ref.fa"

# --- .alt ---
rest=$(( L - 300 ))
{
  printf 'chrP_alt\t0\tchrP\t1\t60\t300M5I%dM\t*\t0\t0\t*\t*\n' "$rest"
  printf 'chrR_alt\t16\tchrP\t251\t60\t200M\t*\t0\t0\t*\t*\n'
} > "$d/ref.fa.alt"

Q=$(printf '%150s' '' | tr ' ' 'I')   # 150 I quals (long enough for any read)
emit() {  # emit <name> <seq>
  n="$1"; s="$2"; ql=$(printf '%s' "$Q" | cut -c1-${#s})
  printf '@%s\n%s\n+\n%s\n' "$n" "$s" "$ql" >> "$d/reads.fq"
}
: > "$d/reads.fq"

# (a)/(b) forward read in the shared flank: chrP[120,270) (1-based 121..270)
emit pa-fwd "$(printf '%s' "$S" | cut -c121-270)"

# (c) primary-side forward read chrP[300,400) (1-based 301..400)
CF=$(printf '%s' "$S" | cut -c301-400)
emit c-fwd "$CF"
# (c) reverse-complement twin of the SAME read -> hits chrR_alt on minus strand
emit c-rev-twin "$(printf '%s' "$CF" | rc)"

# (d) one-end-in-a-hole: read whose first 5 bases are the inserted "GGGGG"
# (positions 301..305 on chrP_alt, a hole) then 145 liftable bases from chrP_alt
# 306.. (== chrP 301..).  Built directly from chrP_alt so it aligns there.
HOLE=$(printf '%s' "$ALT" | cut -c301-450)   # GGGGG + 145 liftable bases
emit hole-end "$HOLE"

# small (<=tol) read-vs-ALT indel twin: chrP read with a 3bp deletion at offset 60
# relative to chrP_alt's content; placement absorbs it and it still groups (~chrP 120)
TW=$(printf '%s' "$S" | mawk -v s=121 -v e=270 '{x=substr($0,s,e-s+1); print substr(x,1,60) substr(x,64)}')
emit twin-indel "$TW"

echo "fixture: L=$L  chrP_alt(+5I@300)  chrR_alt(rev of chrP[250,450))" >&2
