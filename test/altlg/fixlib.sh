# Shared fixture-builder helpers for the ALT liftover-group suite (sourced by
# mkfixture-*.sh).  Previously copy-pasted (gen into 7 scripts, rc into 4).

# Deterministic random DNA of a given length: gen <seed> <length>.
gen() { python3 -c "
import random
random.seed($1)
print(''.join(random.choice('ACGT') for _ in range($2)))
"; }

# Reverse-complement.  With an argument it complements that string; with none it
# acts as a stdin filter.  Handles upper and lower case.
rc() {
    if [ "$#" -ge 1 ]; then printf '%s' "$1" | rev | tr 'ACGTacgt' 'TGCAtgca'
    else rev | tr 'ACGTacgt' 'TGCAtgca'; fi
}
