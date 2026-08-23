#!/bin/sh
# Fixture for Task 8 (PE reconciliation integration), the RESCUED-ALT-MATE case.
#
# Goal: one mate (R2) cannot seed on its own (peppered with mismatches so no
# k-mer seed survives), so it is recovered ONLY by mate rescue (mb_matesw_core).
# The rescued mate lands inside a primary region that has an identical ALT twin
# contig, so mate rescue produces hits on BOTH chrP and the ALT twin.  For the
# liftover-group reconciliation (Hook A) to fold the rescued ALT hit under its
# primary, the rescued hit must carry is_alt -- but mb_matesw_align does
# memset(h,0), clearing is_alt.  The is_alt stamp in mb_matesw_core fixes that.
#
#   WITHOUT the stamp + Hook A: the rescued ALT hit has is_alt=0, so reconcile is
#       blind to it; it competes as a distinct locus (or mis-groups), depressing
#       the pair MAPQ / mis-flagging the twin.
#   WITH both: the rescued chrP_altG hit is grouped under chrP (secondary, 0x100),
#       chrP is the proper-pair primary, MAPQ stays high.
#
# R1 is clean and unique (anchors the pair).  A block of filler pairs (jittered
# insert) establishes pestat.  R2's true locus is inside the ALT-twinned window.
set -eu
d="$1"; mkdir -p "$d"

python3 - "$d" <<'PY'
import sys, random
d = sys.argv[1]
random.seed(31)
N = 6000
S = ''.join(random.choice('ACGT') for _ in range(N))

def revcomp(s):
    c = {'A':'T','C':'G','G':'C','T':'A'}
    return ''.join(c[b] for b in reversed(s))

def pepper(s, period, seed):
    # flip every `period`-th base so no k-mer (default 11..19) survives a window
    r = random.Random(seed)
    flip = {'A':'C','C':'G','G':'T','T':'A'}
    out = list(s)
    for i in range(period//2, len(out), period):
        out[i] = flip[out[i]]
    return ''.join(out)

# ALT-twinned window on chrP: [3000, 3600) (600bp), duplicated as chrP_altG.
ALT_ST, ALT_EN = 3000, 3600
ALTG = S[ALT_ST:ALT_EN]
LG = len(ALTG)

with open(f"{d}/ref.fa", "w") as f:
    f.write(f">chrP\n{S}\n")
    f.write(f">chrP_altG\n{ALTG}\n")
with open(f"{d}/ref.fa.alt", "w") as f:
    f.write(f"chrP_altG\t0\tchrP\t{ALT_ST+1}\t60\t{LG}M\t*\t0\t0\t*\t*\n")

RL = 120
reads1, reads2 = [], []
def fq(name, r1, r2):
    return (f"@{name}/1\n{r1}\n+\n{'I'*len(r1)}\n", f"@{name}/2\n{r2}\n+\n{'I'*len(r2)}\n")

# Filler pairs (jittered FR insert 400..520) from unique windows clear of the ALT
# window, to establish pestat.
jit = random.Random(77)
pos = 200
n_filler = 30
for i in range(n_filler):
    while ALT_ST - 200 < pos < ALT_EN + 200:
        pos += 600
    if pos + 520 + RL > N - 200:
        pos = 200
    gap = jit.randint(280, 400)
    r1 = S[pos:pos+RL]
    r2 = revcomp(S[pos+gap:pos+gap+RL])
    a, b = fq(f"f{i:02d}", r1, r2)
    reads1.append(a); reads2.append(b)
    pos += 150

# Test pair r-rescue: BOTH mates lie inside the ALT-twinned window [3000,3600), so
# each also has a chrP_altG hit.  R1 = chrP[3080,3200) seeds cleanly (anchors and
# also maps to chrP_altG).  R2's true locus chrP[3440,3560) is peppered so it
# cannot seed and is recovered ONLY by mate rescue -- landing on both chrP and
# chrP_altG.  Because R1 also has a chrP_altG hit, the rescued R2 chrP_altG hit can
# form a competing (altG,altG) PAIR.  Only if the rescued mate carries is_alt does
# Hook A group it (parent != id) so Hook B excludes it from the second pairing;
# without the is_alt stamp the (altG,altG) pair competes and depresses mapq_pe.
r1 = S[3080:3200]
r2_true = revcomp(S[3440:3560])
r2 = pepper(r2_true, 11, 5)
a, b = fq("r-rescue", r1, r2)
reads1.append(a); reads2.append(b)

with open(f"{d}/reads_1.fq", "w") as f: f.writelines(reads1)
with open(f"{d}/reads_2.fq", "w") as f: f.writelines(reads2)

sys.stderr.write(
    f"fixture(pe-rescue): chrP={N} altG=chrP[{ALT_ST},{ALT_EN})->POS{ALT_ST+1} "
    f"R1=chrP[3080,3200) clean, R2=chrP[3440,3560) peppered(period11)->rescue-only\n")
PY

echo "fixture(pe-rescue) written to $d" >&2
