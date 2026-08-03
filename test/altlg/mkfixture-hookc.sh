#!/bin/sh
# Fixture for Hook C UNLIFTABLE-demotion test (test-hookc.sh).
#
# Goal: a PE read whose mate (R2) has an is_alt hit on chrHC_altU whose ENTIRE
# footprint falls in a LIFT HOLE (no .alt block covers those ALT bases).  Hook C
# in pe.c skips chimeric demotion only for DISTINCT LIFTABLE groups; an unliftable
# hit returns 0 from mb_distinct_lifted_group and falls through to NORMAL demotion.
# This fixture proves that path: the unliftable ALT hit must become secondary
# (parent != its id) and must NOT remain a shielded independent competitor.
#
# Geometry:
#   chrHC          4000 bp primary (random seed 41)
#   chrHC_altU     200 bp ALT contig:
#     bases  [0,100)  = chrHC[500,600)  -- covered by .alt block -> liftable
#     bases [100,200) = chrHC[700,800)  -- NO .alt block covers this -> HOLE
#   .alt record: chrHC_altU -> chrHC POS 501, CIGAR 100M  (covers only alt[0,100))
#
#   R1 = chrHC[300,420) fwd  (unique anchor, 120bp)
#   R2 = revcomp(chrHC[700,800)) = revcomp(chrHC_altU[100,200))
#      -> maps to chrHC@700 (primary, liftable) AND to chrHC_altU@100 (ALT, HOLE -> unliftable)
#   Insert (5'-5'): R1@300 fwd, R2@700 rev -> ~520 bp, within filler distribution.
#
#   Filler pairs: 30 unique FR pairs from chrHC windows outside [400,900),
#   jittered insert 400..520 bp, establish pestat.
#
# Unliftability proof: the chrHC_altU footprint [100,200) overlaps NO .alt lift
# block (the single block covers alt[0,100), i.e. alt_en=100 <= fp_st=100), so
# mb_hit_place returns liftable=0.  ex-place-check will confirm this.
set -eu
d="$1"; mkdir -p "$d"

python3 - "$d" <<'PY'
import sys, random
d = sys.argv[1]
random.seed(41)
N = 4000
S = ''.join(random.choice('ACGT') for _ in range(N))

def revcomp(s):
    c = {'A':'T','C':'G','G':'C','T':'A'}
    return ''.join(c[b] for b in reversed(s))

# ALT contig: first 100bp = chrHC[500,600); second 100bp = chrHC[700,800)
ALTL = S[500:600]   # liftable part (covered by .alt block)
ALTH = S[700:800]   # hole part (NOT covered by .alt block)
ALTU = ALTL + ALTH  # 200bp total

# .alt block covers only the first 100bp of chrHC_altU (alt[0,100) -> chrHC POS 501)
# so alt[100,200) is a HOLE -> any hit with footprint inside [100,200) -> liftable=0

with open(f"{d}/ref.fa", "w") as f:
    f.write(f">chrHC\n{S}\n")
    f.write(f">chrHC_altU\n{ALTU}\n")
with open(f"{d}/ref.fa.alt", "w") as f:
    # QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
    # covers chrHC_altU[0,100) -> chrHC POS 501 (1-based)
    f.write("chrHC_altU\t0\tchrHC\t501\t60\t100M\t*\t0\t0\t*\t*\n")

RL1, RL2 = 120, 100  # read lengths

def fq(name, r1, r2):
    q1 = 'I' * len(r1); q2 = 'I' * len(r2)
    return (f"@{name}/1\n{r1}\n+\n{q1}\n", f"@{name}/2\n{r2}\n+\n{q2}\n")

reads1, reads2 = [], []

# Filler pairs: unique FR windows from chrHC, clear of the ALT-relevant region [400,900)
# jittered insert so pestat has real width (inserts ~400..520)
jit = random.Random(77)
pos = 50
n_filler = 30
for i in range(n_filler):
    while 400 - RL1 < pos < 900 + RL2:
        pos += 400
    if pos + 520 + RL2 > N - 100:
        pos = 50
    gap = jit.randint(280, 400)  # 5'-5' insert = gap + RL1 = 400..520
    r1 = S[pos:pos+RL1]
    r2 = revcomp(S[pos+gap:pos+gap+RL2])
    a, b = fq(f"f{i:02d}", r1, r2)
    reads1.append(a); reads2.append(b)
    pos += 100

# Test pair r-hookc:
# R1 = chrHC[300,420) fwd (unique anchor)
# R2 = revcomp(chrHC[700,800)) -- this is also chrHC_altU[100,200) revcomp'd
#   -> primary hit: chrHC@700 (liftable=1)
#   -> ALT hit: chrHC_altU@100 (footprint [100,200) in HOLE -> liftable=0)
# 5'-5' insert: (700+100) - 300 = 500bp, within the filler distribution (400..520)
r1 = S[300:420]
r2 = revcomp(S[700:800])
a, b = fq("r-hookc", r1, r2)
reads1.append(a); reads2.append(b)

with open(f"{d}/reads_1.fq", "w") as f: f.writelines(reads1)
with open(f"{d}/reads_2.fq", "w") as f: f.writelines(reads2)

# Also write R2 alone for the ex-place-check probe (single-read FASTQ)
with open(f"{d}/r2_probe.fq", "w") as f:
    f.write(f"@r-hookc-r2\n{r2}\n+\n{'I'*len(r2)}\n")

sys.stderr.write(
    f"fixture(hookc): chrHC={N} altU=200bp(liftable[0,100)=chrHC[500,600),"
    f"hole[100,200)=chrHC[700,800)) "
    f"R1=chrHC[300,420) R2=revcomp(chrHC[700,800)) insert~500 filler={n_filler}\n")
PY

echo "fixture(hookc) written to $d" >&2
