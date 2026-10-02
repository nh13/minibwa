# Shared helpers for the ALT liftover-group test suite (sourced by test-*.sh).
#
# Every helper here was previously copy-pasted across many test-*.sh scripts,
# sometimes with small and occasionally INCOMPATIBLE per-file variations.  The
# signatures below are supersets that are backward-compatible with every prior
# call site: trailing arguments are optional and default to "match anything", so
# both e.g. `primary_mapq f` (any primary) and `primary_mapq f q` (that read's
# primary) work, and both `flag_of f q c` and `flag_of f q c pos` work.

ok()   { echo "  ok: $1"; }
fail() { echo "FAIL: $1"; exit 1; }

# SAM body: drop @ header lines.
sam_body() { mawk '$1!~/^@/{print}' "$1"; }

# Is bit <b> set in flag <f>?  (mawk has no and(); test by arithmetic.)
has_bit() { mawk -v f="$1" -v b="$2" 'BEGIN{ f=int(f); b=int(b); printf "%d\n", (int(f/b)%2==1)?1:0 }'; }

# FLAG / MAPQ of the <qname> record on contig <ctg>, optionally at 1-based POS <pos>.
flag_of() { mawk -v q="$2" -v c="$3" -v p="${4:-}" '$1==q && $3==c && (p=="" || $4==p){print $2; exit}' "$1"; }
mapq_of() { mawk -v q="$2" -v c="$3" -v p="${4:-}" '$1==q && $3==c && (p=="" || $4==p){print $5; exit}' "$1"; }

# AS:i: tag of the <qname> record on contig <ctg>.
as_of() { mawk -v q="$2" -v c="$3" '$1==q && $3==c {
    for(i=12;i<=NF;i++){ if(substr($i,1,5)=="AS:i:"){ print substr($i,6); exit } } }' "$1"; }

# FLAG / MAPQ of the PRIMARY record (neither secondary 0x100 nor supplementary
# 0x800).  With <qname>, restrict to that read; without it, the first primary.
primary_flag() { mawk -v q="${2:-}" '$1!~/^@/ && (q=="" || $1==q){f=int($2); if(int(f/256)%2==0 && int(f/2048)%2==0){print $2; exit}}' "$1"; }
primary_mapq() { mawk -v q="${2:-}" '$1!~/^@/ && (q=="" || $1==q){f=int($2); if(int(f/256)%2==0 && int(f/2048)%2==0){print $5; exit}}' "$1"; }

# Is there a record on chrP at 1-based POS <pos>?  (chrP is the primary contig in
# every fixture that uses this.)
has_pos() { mawk -v p="$2" 'BEGIN{f=0} $1!~/^@/ && $3=="chrP" && $4==p{f=1} END{print f}' "$1"; }
