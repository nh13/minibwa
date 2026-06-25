#!/bin/sh
# Integration test for Task 5: --alt-records flag and --alt FILE CLI.
#
# Verifies three things:
#   1. WITHOUT --alt-records, no ALT-contig record appears in SAM output
#      (default out_n=0 drops secondaries/supplements on ALT contigs).
#   2. WITH --alt-records, the ALT-contig record IS emitted with FULL-LENGTH SEQ
#      (length == read length, not '*').
#   3. --alt FILE loads a .alt from a non-adjacent path and enables ALT mapping.
#   4. chrM baseline (no .alt): output is byte-identical with/without --alt-records.
#
# Usage: test/altlg/test-output.sh [<minibwa-dir>]
set -eu

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MINIBWA="$MDIR/minibwa"

TMPD=$(mktemp -d /tmp/altlg-output.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

fail() { echo "FAIL: $1"; exit 1; }
ok()   { echo "  ok: $1"; }

# ============================================================
# Build a minimal fixture:
#   chrP (300bp primary) + chrP_altA (150bp ALT contig)
#   .alt maps chrP_altA -> chrP at POS 101 (1-based), 150M
#   read = chrP_altA sequence (unique; maps cleanly to the ALT contig)
# ============================================================
echo "[test-output] building fixture ..."

python3 -c "
import random, sys
random.seed(42)
def seq(n): return ''.join(random.choice('ACGT') for _ in range(n))

pad1  = seq(100)
core  = seq(150)   # shared core: appears in chrP at offset 100, and IS chrP_altA
pad2  = seq(50)

chrp  = pad1 + core + pad2   # 300bp primary
alt   = core                  # 150bp ALT contig == core

print(f'>chrP\n{chrp}')
print(f'>chrP_altA\n{alt}')
" > "$TMPD/ref.fa"

# .alt: chrP_altA (0-based [0,150)) -> chrP POS 101 (1-based), 150M
printf 'chrP_altA\t0\tchrP\t101\t60\t150M\t*\t0\t0\t*\t*\n' > "$TMPD/ref.fa.alt"

# Read: the ALT contig sequence
python3 -c "
import random
random.seed(42)
def seq(n): return ''.join(random.choice('ACGT') for _ in range(n))
seq(100)   # consume pad1 seed state
core = seq(150)
qual = 'I' * 150
print(f'@r-alt')
print(core)
print('+')
print(qual)
" > "$TMPD/reads.fq"

"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null

READ_LEN=$(python3 -c "
with open('$TMPD/reads.fq') as f:
    next(f); seq=next(f).strip(); print(len(seq))
")
echo "[test-output] read length = $READ_LEN; ALT contig = chrP_altA; lifted to chrP POS 101"

# ============================================================
# 1. WITHOUT --alt-records: no chrP_altA record in SAM
# ============================================================
echo "== WITHOUT --alt-records =="
"$MINIBWA" mem "$TMPD/ref.fa" "$TMPD/reads.fq" 2>/dev/null > "$TMPD/no-flag.sam"
alt_lines=$(mawk '$1!~/^@/ && $3=="chrP_altA"' "$TMPD/no-flag.sam" | wc -l | tr -d ' ')
[ "$alt_lines" -eq 0 ] \
    || fail "without --alt-records: $alt_lines chrP_altA record(s) emitted (expected 0)"
ok "without --alt-records: no chrP_altA record (ALT hits suppressed by default)"

# ============================================================
# 2. WITH --alt-records: chrP_altA IS emitted with full-length SEQ
# ============================================================
echo "== WITH --alt-records =="
"$MINIBWA" mem --alt-records "$TMPD/ref.fa" "$TMPD/reads.fq" 2>/dev/null > "$TMPD/with-flag.sam"

# Check the ALT record exists
alt_lines=$(mawk '$1!~/^@/ && $3=="chrP_altA"' "$TMPD/with-flag.sam" | wc -l | tr -d ' ')
[ "$alt_lines" -gt 0 ] \
    || fail "with --alt-records: no chrP_altA record emitted"
ok "with --alt-records: chrP_altA record is present ($alt_lines line(s))"

# Check the SEQ field is full-length (not '*')
alt_seq=$(mawk '$1!~/^@/ && $3=="chrP_altA"{print $10; exit}' "$TMPD/with-flag.sam")
[ "$alt_seq" != "*" ] \
    || fail "with --alt-records: SEQ field for chrP_altA is '*' (expected full SEQ)"
alt_seq_len=${#alt_seq}
[ "$alt_seq_len" -eq "$READ_LEN" ] \
    || fail "with --alt-records: SEQ length $alt_seq_len != read length $READ_LEN"
ok "with --alt-records: SEQ is full-length ($alt_seq_len == $READ_LEN)"

# Report the FLAG bit for Task 6 reference
alt_flag=$(mawk '$1!~/^@/ && $3=="chrP_altA"{print $2; exit}' "$TMPD/with-flag.sam")
echo "  (info) ALT record FLAG = $alt_flag (0x$(printf '%x' $alt_flag)) -- noted for Task 6"

# ============================================================
# 3. --alt FILE: load .alt from a non-adjacent path
# ============================================================
echo "== --alt FILE from a non-adjacent path =="
mkdir -p "$TMPD/altdir"
cp "$TMPD/ref.fa.alt" "$TMPD/altdir/custom.alt"

# Build index without adjacent .alt, then specify it via --alt
python3 -c "
import random
random.seed(99)
def seq(n): return ''.join(random.choice('ACGT') for _ in range(n))
pad = seq(50); core = seq(100); pad2 = seq(50)
chrp = pad + core + pad2
alt  = core
print(f'>chrQ\n{chrp}')
print(f'>chrQ_altX\n{alt}')
" > "$TMPD/ref2.fa"
printf 'chrQ_altX\t0\tchrQ\t51\t60\t100M\t*\t0\t0\t*\t*\n' > "$TMPD/altdir/ref2.alt"
python3 -c "
import random
random.seed(99)
def seq(n): return ''.join(random.choice('ACGT') for _ in range(n))
seq(50); core = seq(100); qual = 'I'*100
print('@r-alt2'); print(core); print('+'); print(qual)
" > "$TMPD/reads2.fq"
"$MINIBWA" index "$TMPD/ref2.fa" 2>/dev/null

# With --alt pointing to the non-adjacent file, ALT record should appear
"$MINIBWA" mem --alt-records --alt "$TMPD/altdir/ref2.alt" \
    "$TMPD/ref2.fa" "$TMPD/reads2.fq" 2>/dev/null > "$TMPD/alt-file.sam"
alt2_lines=$(mawk '$1!~/^@/ && $3=="chrQ_altX"' "$TMPD/alt-file.sam" | wc -l | tr -d ' ')
[ "$alt2_lines" -gt 0 ] \
    || fail "--alt FILE: no chrQ_altX record (--alt did not load the non-adjacent .alt)"
ok "--alt FILE: chrQ_altX ALT record present ($alt2_lines line(s)) from non-adjacent .alt"

# WITHOUT --alt FILE and no adjacent .alt: chrQ_altX is not is_alt, so it
# appears as a plain primary (FLAG has no 0x100/0x800 bits set), not an ALT hit.
# With --alt FILE: chrQ_altX IS is_alt -> FLAG has 0x100 or 0x800.
"$MINIBWA" mem --alt-records "$TMPD/ref2.fa" "$TMPD/reads2.fq" 2>/dev/null > "$TMPD/no-alt-file.sam"
alt2_lines_without=$(mawk '$1!~/^@/ && $3=="chrQ_altX"' "$TMPD/no-alt-file.sam" | wc -l | tr -d ' ')
[ "$alt2_lines_without" -gt 0 ] \
    || fail "--alt FILE sanity: without --alt, chrQ_altX record missing; cannot compare FLAG classification"
flag_without=$(mawk '$1!~/^@/ && $3=="chrQ_altX"{print $2; exit}' "$TMPD/no-alt-file.sam")
flag_with=$(mawk '$1!~/^@/ && $3=="chrQ_altX"{print $2; exit}' "$TMPD/alt-file.sam")
# Without --alt: chrQ_altX is primary (not secondary/supplementary)
[ $(( ${flag_without:-0} & 0x900 )) -eq 0 ] \
    || fail "--alt FILE sanity: without --alt, chrQ_altX FLAG=$flag_without has 0x100/0x800 set (unexpected)"
ok "--alt FILE sanity: without --alt, chrQ_altX FLAG=$flag_without (primary, not ALT-classified)"
# With --alt: chrQ_altX is secondary or supplementary (0x100 or 0x800)
[ $(( ${flag_with:-0} & 0x900 )) -ne 0 ] \
    || fail "--alt FILE: with --alt, chrQ_altX FLAG=$flag_with lacks 0x100/0x800 (not treated as ALT)"
ok "--alt FILE: with --alt, chrQ_altX FLAG=$flag_with (secondary/supplementary = ALT-classified)"

# ============================================================
# 4. chrM baseline: byte-identical with/without --alt-records (no .alt => inert)
# ============================================================
echo "== chrM baseline (no .alt; --alt-records inert) =="
CHRM_FA="$MDIR/test/chrM-human.fa.gz"
CHRM_R1="$MDIR/test/chrM-read_1.fa.gz"
CHRM_R2="$MDIR/test/chrM-read_2.fa.gz"
# Filter helper: strip header lines for body-only comparison
sam_body() { mawk '$1!~/^@/{print}' "$1"; }

if [ -f "$CHRM_FA" ] && [ -f "$CHRM_R1" ]; then
    "$MINIBWA" mem --outn=5 "$CHRM_FA" "$CHRM_R1" 2>/dev/null > "$TMPD/chrM-se-no-flag.sam"
    "$MINIBWA" mem --outn=5 --alt-records "$CHRM_FA" "$CHRM_R1" 2>/dev/null > "$TMPD/chrM-se-with-flag.sam"
    [ -s "$TMPD/chrM-se-no-flag.sam" ] || fail "chrM SE baseline: empty output"
    sam_body "$TMPD/chrM-se-no-flag.sam"   > "$TMPD/chrM-se-no-flag.body"
    sam_body "$TMPD/chrM-se-with-flag.sam" > "$TMPD/chrM-se-with-flag.body"
    cmp -s "$TMPD/chrM-se-no-flag.body" "$TMPD/chrM-se-with-flag.body" \
        || fail "chrM SE baseline: --alt-records changed alignment records without .alt"
    ok "chrM SE baseline: alignment records byte-identical with/without --alt-records (no .alt => inert)"
    if [ -f "$CHRM_R2" ]; then
        "$MINIBWA" mem --outn=5 "$CHRM_FA" "$CHRM_R1" "$CHRM_R2" 2>/dev/null > "$TMPD/chrM-pe-no-flag.sam"
        "$MINIBWA" mem --outn=5 --alt-records "$CHRM_FA" "$CHRM_R1" "$CHRM_R2" 2>/dev/null > "$TMPD/chrM-pe-with-flag.sam"
        sam_body "$TMPD/chrM-pe-no-flag.sam"   > "$TMPD/chrM-pe-no-flag.body"
        sam_body "$TMPD/chrM-pe-with-flag.sam" > "$TMPD/chrM-pe-with-flag.body"
        cmp -s "$TMPD/chrM-pe-no-flag.body" "$TMPD/chrM-pe-with-flag.body" \
            || fail "chrM PE baseline: --alt-records changed alignment records without .alt"
        ok "chrM PE baseline: alignment records byte-identical with/without --alt-records (no .alt => inert)"
    else
        echo "  skip: chrM R2 not found ($CHRM_R2)"
    fi
else
    echo "  skip: chrM baseline files not found ($CHRM_FA)"
fi

echo "[test-output] PASS"
exit 0
