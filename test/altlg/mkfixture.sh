#!/bin/sh
# Synthetic ref: chrP (primary) + chrP_alt (ALT, = chrP with a 5bp insertion at
# offset 300, so the .alt CIGAR is 300M5I<rest>M).
# Emits ref.fa, ref.fa.alt (SAM: chrP_alt aligned to chrP), and reads.
set -eu
d="$1"; mkdir -p "$d"
S=$(printf '%s' 'GATCCTAGCATGCTAGGCTAACGTTAGCCGATCGTAGCTAGGCATCGATCGTAGCTAGCTAGGCATCGATTACGATCGGCTAATCGATCGTAGCTGATCGA'\
'TCGTAGCTAGCATCGATCGTAGCATCGGCTAGCATCGATCGATTACGCATCGATCGTAGGCTAGCATCGATCGTAGCTAGCATCGGCTAGCATCGATTACG'\
'ATCGGCTAATCGATCGTAGCTGATCGATCGTAGCTAGCATCGATCGTAGCATCGGCTAGCATCGATCGATTACGCATCGATCGTAGGCTAGCATCGATCGT'\
'AGCTAGCATCGGCTAGCATCGATTACGATCGGCTAATCGATCGTAGCTGATCGAGATCCTAGCATGCTAGGCTAACGTTAGCCGATCGTAGCTAGGCATCG'\
'ATCGTAGCTAGCTAGGCATCGATTACGATCGGCTAATCGATCGTAGCTGATCGATCGTAGCTAGCATCGATCGTAGCATCGGCTAGCATCGATCGATTACG')
printf '>chrP\n%s\n' "$S" > "$d/ref.fa"
ALT=$(printf '%s' "$S" | mawk '{print substr($0,1,300) "GGGGG" substr($0,301)}')   # 5bp insertion at 300
printf '>chrP_alt\n%s\n' "$ALT" >> "$d/ref.fa"
# .alt: chrP_alt (query) aligned to chrP (ref): 300M5I<rest>M ; POS=1, MAPQ 60
rest=$(( ${#S} - 300 ))
printf 'chrP_alt\t0\tchrP\t1\t60\t300M5I%dM\t*\t0\t0\t*\t*\n' "$rest" > "$d/ref.fa.alt"
# read fully inside the shared flank (offset 120, 150bp) -> lifts cleanly
printf '@r1\n%s\n+\n%s\n' "$(printf '%s' "$S" | cut -c121-270)" "$(printf '%150s' '' | tr ' ' 'I')" > "$d/r1.fq"
