#!/bin/sh
# Integration test for l2b_set_alt / l2b_lift (Task 1: span-lift index).
# Builds a synthetic fixture, indexes it with minibwa, and verifies:
#   1. chrP_alt is flagged ALT (ex-lift-check exits 0).
#   2. lift(100) succeeded (ok=1) and maps to pri=100 (identity, in block 0).
#   3. lift(400) succeeded (ok=1) and maps to pri=395 (5bp insertion shifts by 5).
# Usage: test/altlg/test-lift.sh [<minibwa-dir>]
# <minibwa-dir> defaults to two levels above this script's directory.
set -eu

MDIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
MKFIXTURE="$MDIR/test/altlg/mkfixture.sh"
MINIBWA="$MDIR/minibwa"
EX_LIFT="$MDIR/api-test/ex-lift-check"

TMPD=$(mktemp -d /tmp/altlg-test.XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

echo "[test-lift] building fixture in $TMPD ..."
/bin/sh "$MKFIXTURE" "$TMPD"

echo "[test-lift] indexing with minibwa ..."
"$MINIBWA" index "$TMPD/ref.fa" 2>/dev/null

echo "[test-lift] running ex-lift-check ..."
OUT=$("$EX_LIFT" "$TMPD/ref.fa" 2>&1)
echo "[test-lift] output: $OUT"

# --- assertion 1: chrP_alt is flagged ALT (program exits non-zero and prints error if not) ---
# exit-code already checked by set -e above.

# Output format: lift(100)=<ok1> pri=<pp1>  lift(400)=<ok2> pri=<pp2>

# --- assertion 2a: lift(100) succeeded (ok1=1) ---
OK1=$(printf '%s\n' "$OUT" | mawk '{if(match($0,/lift\(100\)=[0-9]+/)) print substr($0,RSTART+10,RLENGTH-10)}')
if [ "$OK1" != "1" ]; then
    echo "FAIL: lift(100) expected ok=1, got '$OK1'"
    exit 1
fi

# --- assertion 2b: lift(100) maps to pri=100 (identity: pos 100 is in block 0, no shift) ---
PRI1=$(printf '%s\n' "$OUT" | mawk '{if(match($0,/lift\(100\)=[0-9]+ pri=[0-9]+/)) {s=substr($0,RSTART,RLENGTH); if(match(s,/pri=[0-9]+/)) print substr(s,RSTART+4,RLENGTH-4)}}')
if [ "$PRI1" != "100" ]; then
    echo "FAIL: lift(100) expected pri=100, got '$PRI1'"
    exit 1
fi

# --- assertion 3a: lift(400) succeeded (ok2=1) ---
OK2=$(printf '%s\n' "$OUT" | mawk '{if(match($0,/lift\(400\)=[0-9]+/)) print substr($0,RSTART+10,RLENGTH-10)}')
if [ "$OK2" != "1" ]; then
    echo "FAIL: lift(400) expected ok=1, got '$OK2'"
    exit 1
fi

# --- assertion 3b: lift(400) maps to pri=395 (5bp insertion shifts coordinates by 5) ---
PRI2=$(printf '%s\n' "$OUT" | mawk '{if(match($0,/lift\(400\)=[0-9]+ pri=[0-9]+/)) {s=substr($0,RSTART,RLENGTH); if(match(s,/pri=[0-9]+/)) print substr(s,RSTART+4,RLENGTH-4)}}')
if [ "$PRI2" != "395" ]; then
    echo "FAIL: lift(400) expected pri=395, got '$PRI2'"
    exit 1
fi

echo "[test-lift] PASS: chrP_alt flagged ALT, lift(100) ok=1 pri=100, lift(400) ok=1 pri=395"
exit 0
