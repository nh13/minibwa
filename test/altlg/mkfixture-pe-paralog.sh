#!/bin/sh
# Fixture for Task 8 (PE reconciliation integration), the PARALOG GUARD case
# (Hook C): the chimeric demotion loop at pe.c:524-535 demotes a hit `p` under the
# chosen pair-hit `h[r]` when their QUERY spans overlap a lot (ol > mask_level).
# Without Hook C that demotion fires on two hits that are DIFFERENT lifted groups
# (different primary loci / paralogs), merging them and inflating MAPQ — undoing
# the Task-4(b) paralog guard for PE.  Hook C skips the demotion when `p` and `h[r]`
# are different lifted groups (different pri_tid/rev or |Δlifted_st| > MB_LIFT_TOL).
#
# Geometry: chrP carries two near-identical paralog copies of COPY (COPY1, COPY2)
# whose ALT twins lift just OUTSIDE MB_LIFT_TOL apart (the same 2*tol boundary the
# SE paralog fixture stresses).  One mate (R1) of the pair IS the ambiguous COPY
# sequence, so it has competing hits at COPY1, COPY2, altA, altB — overlapping
# query spans but DIFFERENT lifted groups.  The other mate (R2) is unique and
# anchors the proper pair near COPY1.  A block of filler pairs establishes pestat.
#
#   BEFORE Hook C: :526 demotes the COPY2/altB-group hit under the COPY1 pair-hit
#                  (query spans fully overlap) => groups merge => the read looks
#                  uniquely placed => MAPQ inflated.
#   AFTER  Hook C: the cross-group demotion is skipped => COPY2 stays a distinct
#                  competing locus => R1 remains a low-MAPQ multi-mapper.
set -eu
d="$1"; mkdir -p "$d"

python3 - "$d" <<'PY'
import sys, random
d = sys.argv[1]

def gen(seed, n):
    r = random.Random(seed)
    return ''.join(r.choice('ACGT') for _ in range(n))

def revcomp(s):
    c = {'A':'T','C':'G','G':'C','T':'A'}
    return ''.join(c[b] for b in reversed(s))

PAD1 = gen(1, 2000)
COPY = gen(2, 150)
GAP  = gen(3, 220)
PAD2 = gen(4, 2000)
CHRP = PAD1 + COPY + GAP + COPY + PAD2
N = len(CHRP)
# COPY1 0-based start = 2000 ; COPY2 0-based start = 2000+150+220 = 2370.
C1 = len(PAD1)                    # 2000
C2 = len(PAD1) + len(COPY) + len(GAP)  # 2370
LC = len(COPY)

# ALT twins: each contig is just COPY (150bp); altA->COPY1 locus, altB 20bp
# downstream (2*MB_LIFT_TOL), just outside tolerance -> distinct groups.
with open(f"{d}/ref.fa", "w") as f:
    f.write(f">chrP\n{CHRP}\n")
    f.write(f">chrP_altA\n{COPY}\n")
    f.write(f">chrP_altB\n{COPY}\n")
with open(f"{d}/ref.fa.alt", "w") as f:
    f.write(f"chrP_altA\t0\tchrP\t{C1+1}\t60\t{LC}M\t*\t0\t0\t*\t*\n")
    f.write(f"chrP_altB\t0\tchrP\t{C1+1+20}\t60\t{LC}M\t*\t0\t0\t*\t*\n")

QUAL = 'I' * LC
reads1, reads2 = [], []

def fq(name, r1, r2):
    q1 = 'I' * len(r1); q2 = 'I' * len(r2)
    return (f"@{name}/1\n{r1}\n+\n{q1}\n", f"@{name}/2\n{r2}\n+\n{q2}\n")

# Filler pairs from unique windows of PAD1/PAD2 to establish pestat.  Jitter the
# insert so pestat's proper-pair window has real width (a fixed insert collapses
# std.dev to 0 and the [lo,hi] window to a single value).  Inserts span ~400..520.
n_filler = 30
jit = random.Random(99)
starts = []
pos = 100
for _ in range(n_filler):
    starts.append(pos); pos += 55
for i, st in enumerate(starts):
    gap = jit.randint(280, 400)            # 5'-5' insert = gap + 120 = 400..520
    r1 = CHRP[st:st+120]
    r2 = revcomp(CHRP[st+gap:st+gap+120])
    a, b = fq(f"f{i:02d}", r1, r2)
    reads1.append(a); reads2.append(b)

# Test pair r-para: R1 = COPY (ambiguous: COPY1/COPY2/altA/altB).  R2 = unique
# window downstream of COPY1, revcomp, giving an FR pair (insert 460) anchored
# near COPY1, squarely inside the filler distribution.
r1 = COPY
r2 = revcomp(CHRP[C1+340:C1+340+120])
a, b = fq("r-para", r1, r2)
reads1.append(a); reads2.append(b)

with open(f"{d}/reads_1.fq", "w") as f: f.writelines(reads1)
with open(f"{d}/reads_2.fq", "w") as f: f.writelines(reads2)

sys.stderr.write(
    f"fixture(pe-paralog): chrP={N} COPY1@{C1} COPY2@{C2} altA->{C1+1} altB->{C1+21} "
    f"(Δ20>2*tol) filler={n_filler} r-para(R1=COPY ambiguous)\n")
PY

echo "fixture(pe-paralog) written to $d" >&2
