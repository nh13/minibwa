#!/bin/sh
# Fixture for Task 4 (reconciliation pass mb_reconcile_alt), cases (a) and (c).
#
# chrP is a single 800bp unique (non-repetitive) sequence (LCG PRNG, seed 7) so
# the ONLY competitor for any read is its ALT twin, never a paralogous primary
# locus.  Two well-separated regions, each with exactly ONE ALT twin.  Crucially
# each ALT is a SHORT contig covering only its own region (not a full-length
# copy of chrP), so the region-A twin cannot also compete in region D.
#
#   Region A  chrP[200,350)  -- read r-clean.  Twin chrP_altC = chrP[150,400)
#             IDENTICAL (250bp), aligned to chrP at POS=151 (250M).  The read
#             maps equally well to chrP and chrP_altC: WITHOUT reconciliation the
#             twin is a co-scoring competitor => MAPQ 0; WITH it they group =>
#             chrP is sam_pri, MAPQ>0, chrP_altC secondary (case a).
#
#   Region D  chrP[500,620)  -- read r-diverge.  Twin chrP_altD = chrP[450,670)
#             (220bp) with ONE SNP inside [500,620), aligned to chrP at POS=451
#             (220M).  The read scores STRICTLY LOWER on chrP_altD (AS ~230 vs
#             the primary's 240) yet still lifts to the SAME locus and survives
#             the score-ratio prune (Task 3 guard) -- one SNP costs score without
#             fragmenting every seed (2+ closely-spaced SNPs would kill seeding
#             on the short ALT contig and the twin would never reach this pass).
#             The lower-scoring twin is a real same-group member and must NOT
#             zero the MAPQ (case c).
set -eu
d="$1"; mkdir -p "$d"

# 800bp unique sequence (LCG PRNG, seed 7).
S=$(python3 -c "
import random
random.seed(7)
print(''.join(random.choice('ACGT') for _ in range(800)))
")

# Region-A twin: identical 250bp window chrP[150,400) (1-based 151..400).
ALTC=$(printf '%s' "$S" | cut -c151-400)
LC=${#ALTC}

# Region-D twin: 220bp window chrP[450,670) (1-based 451..670) with ONE SNP inside
# the r-diverge window chrP[500,620).  Window-local 0-based offset 90 maps to
# chrP 540 (read position 40 of 120) -- one mismatch, AS ~230 < the primary's 240.
ALTD=$(printf '%s' "$S" | cut -c451-670 | python3 -c "
import sys
s=list(sys.stdin.read().strip())
flip={'A':'T','C':'G','G':'C','T':'A','a':'t','c':'g','g':'c','t':'a'}
s[90] = flip[s[90]]   # window-local 0-based -> chrP 540 (inside [500,620))
print(''.join(s))
")
LD=${#ALTD}

# --- reference ---
printf '>chrP\n%s\n'      "$S"    > "$d/ref.fa"
printf '>chrP_altC\n%s\n' "$ALTC">> "$d/ref.fa"
printf '>chrP_altD\n%s\n' "$ALTD">> "$d/ref.fa"

# --- .alt: each ALT aligns full-length forward to chrP at its window offset ---
printf 'chrP_altC\t0\tchrP\t151\t60\t%dM\t*\t0\t0\t*\t*\n' "$LC"  > "$d/ref.fa.alt"
printf 'chrP_altD\t0\tchrP\t451\t60\t%dM\t*\t0\t0\t*\t*\n' "$LD" >> "$d/ref.fa.alt"

# --- reads ---
QUAL150=$(printf '%150s' '' | tr ' ' 'I')
QUAL120=$(printf '%120s' '' | tr ' ' 'I')
# r-clean: chrP[200,350) (1-based 201..350) -- region A, identical twin.
RCLEAN=$(printf '%s' "$S" | cut -c201-350)
printf '@r-clean\n%s\n+\n%s\n' "$RCLEAN" "$QUAL150" > "$d/reads.fq"
# r-diverge: chrP[500,620) (1-based 501..620) -- region D, divergent twin.
RDIV=$(printf '%s' "$S" | cut -c501-620)
printf '@r-diverge\n%s\n+\n%s\n' "$RDIV" "$QUAL120" >> "$d/reads.fq"

echo "fixture(reconcile): chrP=800  altC=identical chrP[150,400)  altD=2SNP chrP[450,670)" >&2
