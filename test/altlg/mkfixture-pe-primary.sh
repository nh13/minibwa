#!/bin/sh
# Fixture for the opt-in PE-pair-primary selection (--pe-pair-primary).
#
# Bug it locks: minibwa's PE pairing (mb_pair_hits) correctly chooses the
# mate-consistent endpoint among near-equal paralog copies, but mb_set_sam_pri
# then emits the SAM primary by per-read DP-score/hash order -- which can be the
# OTHER (mate-inconsistent) copy.  --pe-pair-primary makes mb_set_sam_pri honor
# the pair-chosen endpoint.  (No ALT contigs here: this is a primary-vs-primary
# paralog scenario, exercising the general PE path.)
#
# The bug ONLY manifests when the read also carries an ALT hit (r_any_alt true):
# pe.c's Hook C skips demoting a distinct-lifted-group competitor (the paralog-
# safety guard), so the mate-INconsistent paralog survives as a representative,
# and the default mb_set_sam_pri picks it by score order.  Without an ALT hit,
# Hook C demotes it and the bug is masked -- so the fixture includes an ALT
# contig to reproduce the real (subtelomeric, ALT-region) scenario.
#
# Geometry (chrP = 4000bp unique LCG sequence, seed 29; 150bp motif M):
#   copyA = chrP[500,650)   <- M with ONE mismatch (DP score slightly lower)
#   copyB = chrP[2500,2650) <- EXACT M (highest score -> sorts first)
#   chrP_altX               <- EXACT M; .alt maps it onto copyB (POS 2501), so
#                              R1's altX hit folds under copyB and makes
#                              r_any_alt true for R1.
#   R1 (r-para/1) = M       -> hits copyA, copyB (+ folded altX).
#   R2 (r-para/2) = revcomp(chrP[900,1050)) -> UNIQUE, ~550bp 5'-5' from copyA;
#                              forms a proper FR pair with copyA only (copyB is
#                              ~1600bp away, outside the proper-pair window).
#   => pairing chooses copyA; copyA and copyB are distinct lifted groups so Hook C
#      keeps copyB a representative.  Default emits copyB (higher score) as R1
#      primary -- WRONG.  --pe-pair-primary emits copyA (mate-consistent).
#
# 30 filler unique FR pairs (insert ~400-520, jittered) establish pestat.
set -eu
d="${1:?usage: mkfixture-pe-primary.sh <outdir>}"; mkdir -p "$d"

python3 - "$d" <<'PY'
import sys, random
d = sys.argv[1]
random.seed(29)
N = 4000
S = list(''.join(random.choice('ACGT') for _ in range(N)))

def revcomp(s):
    c = {'A':'T','C':'G','G':'C','T':'A'}
    return ''.join(c[b] for b in reversed(s))

# Place a 150bp motif M at copyA=[500,650) and copyB=[2500,2650).  copyB is an
# EXACT copy of M; copyA carries ONE mismatch so its DP score is strictly lower
# (but the gap stays well under pen_unpair, so the proper pair is still applied).
# R1 = M then sorts copyB FIRST by score, so the default (per-read order) emits
# copyB as primary -- the WRONG, mate-inconsistent copy.  --pe-pair-primary
# instead emits copyA (the mate-consistent, pair-chosen endpoint).
M = ''.join(random.choice('ACGT') for _ in range(150))
S[2500:2650] = list(M)                 # copyB: exact
flip = {'A':'C','C':'G','G':'T','T':'A'}
A = list(M)
A[75] = flip[A[75]]                    # 1 mismatch, mid-motif (seeds still hit)
S[500:650] = A                          # copyA: 1 mismatch vs M
S = ''.join(S)

with open(f"{d}/ref.fa", "w") as f:
    f.write(f">chrP\n{S}\n")
    f.write(f">chrP_altX\n{M}\n")       # ALT contig: exact copy of M
# .alt: chrP_altX aligns full-length forward onto copyB (1-based POS=2501), so
# R1's altX hit lifts onto copyB and folds under it (=> r_any_alt true for R1).
with open(f"{d}/ref.fa.alt", "w") as f:
    f.write(f"chrP_altX\t0\tchrP\t2501\t60\t150M\t*\t0\t0\t*\t*\n")

RL = 150
QUAL = 'I' * RL
def fq(name, r1, r2):
    return (f"@{name}/1\n{r1}\n+\n{QUAL}\n", f"@{name}/2\n{r2}\n+\n{QUAL}\n")

reads1, reads2 = [], []
# FILLER: 30 unique FR pairs outside [500,650) and [2500,2650], insert ~400-520.
jit = random.Random(7)
starts = []
pos = 1000
for _ in range(30):
    while (480 < pos < 700) or (2480 < pos < 2700):
        pos += 300
    if pos + 520 + RL > N - 100:
        pos = 1000
    starts.append(pos); pos += 70
for i, st in enumerate(starts):
    gap = jit.randint(250, 370)          # 5'-5' insert = gap + RL = 400..520
    r1 = S[st:st+RL]
    r2 = revcomp(S[st+gap:st+gap+RL])
    a, b = fq(f"f{i:02d}", r1, r2)
    reads1.append(a); reads2.append(b)

# TEST pair: R1 = motif M (maps to copyA and copyB); R2 unique near copyA.
r1 = M
r2 = revcomp(S[900:1050])
a, b = fq("r-para", r1, r2)
reads1.append(a); reads2.append(b)

with open(f"{d}/reads_1.fq", "w") as f: f.writelines(reads1)
with open(f"{d}/reads_2.fq", "w") as f: f.writelines(reads2)
sys.stderr.write("fixture(pe-primary): chrP=4000 M@[500,650)&[2500,2650) "
                 "R2 unique near copyA; 30 filler FR pairs\n")
PY
echo "fixture(pe-primary) written to $d" >&2
