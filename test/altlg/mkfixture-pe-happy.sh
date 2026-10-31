#!/bin/sh
# Fixture for Task 8 (PE reconciliation integration), the HAPPY-PATH case.
#
# Goal: a paired read whose BOTH mates land inside a primary region that also has
# an identical ALT twin contig.  On the primary the pair (R1@chrP, R2@chrP) is the
# best PAIR; on the ALT twin the pair (R1@altF, R2@altF) is an EQUALLY-GOOD second
# PAIR (same insert, same score, identical sequence).  In mb_pair_hits, pairs only
# form within ONE contig (pe.c:174 requires hi->tid == hk->tid), so the twin pair
# is a real competitor:
#
#   BEFORE Hook B:  paux.sub_sc ~= paux.score  =>  mapq_pe ~= 0 (the pair looks
#                   like a 2-way multi-mapper even though both placements are the
#                   SAME locus on primary vs ALT).
#   AFTER  Hook B:  the ALT (non-rep group member) hits are excluded from pair
#                   enumeration  =>  the twin pair never forms  =>  paux.sub_sc is
#                   the true second-best (none) => mapq_pe > 0; chrP is sam_pri for
#                   BOTH mates, the chrP_altF twin records are secondary.
#
# pestat needs >= MIN_DIR_CNT (20) unique FR pairs to not fail, so we emit a block
# of FILLER pairs from unique chrP windows (no ALT twin) at a consistent ~480bp FR
# insert.  Those establish the insert distribution; the single test pair (r-twin)
# exercises the feature.  Filler reads are 120bp, well-separated, each unique.
#
# Reference layout (chrP, 6000bp unique LCG sequence, seed 21):
#   - chrP[2400,2940) (540bp) is duplicated as ALT contig chrP_altF; .alt aligns
#     chrP_altF full-length forward to chrP at POS=2401.
#   - r-twin: R1 = chrP[2480,2600) fwd, R2 = revcomp(chrP[2820,2940)).  Insert
#     (5' to 5') ~= 460bp, FR orientation, well within the filler distribution.
#     Both mates lie inside [2400,2940) so both also map to chrP_altF.
set -eu
d="$1"; mkdir -p "$d"

python3 - "$d" <<'PY'
import sys, random
d = sys.argv[1]
random.seed(21)
N = 6000
S = ''.join(random.choice('ACGT') for _ in range(N))

def revcomp(s):
    c = {'A':'T','C':'G','G':'C','T':'A'}
    return ''.join(c[b] for b in reversed(s))

# ALT-twinned window on chrP: [2400, 2940) (0-based), 540bp.
ALT_ST, ALT_EN = 2400, 2940
ALTF = S[ALT_ST:ALT_EN]
LF = len(ALTF)

# --- reference ---
with open(f"{d}/ref.fa", "w") as f:
    f.write(f">chrP\n{S}\n")
    f.write(f">chrP_altF\n{ALTF}\n")

# --- .alt: chrP_altF aligns full-length forward to chrP at 1-based POS=2401 ---
with open(f"{d}/ref.fa.alt", "w") as f:
    f.write(f"chrP_altF\t0\tchrP\t{ALT_ST+1}\t60\t{LF}M\t*\t0\t0\t*\t*\n")

RL = 120
QUAL = 'I' * RL

def fq(name, r1, r2):
    return (f"@{name}/1\n{r1}\n+\n{QUAL}\n", f"@{name}/2\n{r2}\n+\n{QUAL}\n")

reads1, reads2 = [], []

# --- FILLER pairs: 30 unique FR pairs spread across chrP, insert ~480bp, all
#     OUTSIDE the ALT window so each is a clean unique pair for pestat. ---
filler_starts = []
pos = 200
n_filler = 30
for _ in range(n_filler):
    # keep filler fragments clear of the ALT window [2400,2940)
    while ALT_ST - 120 < pos < ALT_EN + 120:
        pos += 600
    if pos + 480 + RL > N - 200:
        pos = 200
    filler_starts.append(pos)
    pos += 150
jit = random.Random(99)
for i, st in enumerate(filler_starts):
    # Jitter the 5'-5' insert so the proper-pair window has real width (a fixed
    # insert collapses pestat's std.dev to 0 and the [lo,hi] window to a single
    # value, which would reject the test pair).  Inserts span ~400..520.
    gap = jit.randint(280, 400)            # 5'-5' insert = gap + RL = 400..520
    if st + gap + RL > N - 200:
        gap = 280
    r1 = S[st:st+RL]
    r2 = revcomp(S[st+gap:st+gap+RL])
    a, b = fq(f"f{i:02d}", r1, r2)
    reads1.append(a); reads2.append(b)

# --- TEST pair r-twin: both mates inside the ALT-twinned window. ---
# R1 = chrP[2480,2600) fwd ; R2 = revcomp(chrP[2820,2940)).  5'-5' insert = 460,
# squarely inside the filler distribution (400..520).  Both mates lie in the ALT
# window [2400,2940) so each also maps to chrP_altF.
r1 = S[2480:2600]
r2 = revcomp(S[2820:2940])
a, b = fq("r-twin", r1, r2)
reads1.append(a); reads2.append(b)

with open(f"{d}/reads_1.fq", "w") as f:
    f.writelines(reads1)
with open(f"{d}/reads_2.fq", "w") as f:
    f.writelines(reads2)

sys.stderr.write(
    f"fixture(pe-happy): chrP={N} altF=chrP[{ALT_ST},{ALT_EN})->POS{ALT_ST+1} "
    f"filler={n_filler}FR(ins480) r-twin(ins460,both mates in ALT window)\n")
PY

echo "fixture(pe-happy) written to $d" >&2
