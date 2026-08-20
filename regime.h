#ifndef MB_REGIME_H
#define MB_REGIME_H

#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Which SA-lookup backend a regime uses. cp_occ is a later PR (M3); M2 only
 * ever discovers MB_BACKEND_BWT regimes. */
typedef enum { MB_BACKEND_BWT, MB_BACKEND_CP_OCC } mb_backend_t;

typedef struct {
	mb_backend_t backend;
	int sa_bit;
	uint64_t est_ram;
	int speed_rank;
	uint32_t mode_mask;
	char name[16];      /* "sa16","sa8","sa8-b2" */
	char sa_path[1024]; /* "" = bundled in .mbw */
} mb_regime_t;

#define MB_MODE_SRPE 1u
#define MB_MODE_METH 2u
#define MB_MODE_HIC  4u
#define MB_MODE_LR   8u

/* Discover which SA-density regimes exist on disk for `prefix` (i.e.
 * `<prefix>.l2b` and `<prefix>.mbw`, plus any `<prefix>.sa.u*` sidecars).
 * Writes up to `max` regimes into `out` and returns the count (0 if the
 * base index files are missing). `b2_available` is accepted for forward
 * compatibility with M3's cp_occ regimes; M2 always discovers zero of them. */
int mb_regime_discover(const char *prefix, int b2_available, mb_regime_t *out, int max);

/* Detect the usable memory budget: min(host available, cgroup limit,
 * user_cap_bytes). A zero cap/limit/host value means "no such constraint"
 * and is ignored when computing the min. Pass 0 for user_cap_bytes to
 * ignore the user cap entirely. */
uint64_t mb_mem_budget(uint64_t user_cap_bytes);

/* Pick the fastest regime (by speed_rank) that is eligible for `mode` and
 * fits within `budget` bytes. If `forced` is non-NULL, return the index of
 * the regime whose name matches it exactly (ignoring budget), or -1 if no
 * such eligible regime exists. Returns -1 if nothing fits/matches. */
int mb_regime_pick(const mb_regime_t *r, int n, uint64_t budget, uint32_t mode, const char *forced);

/* Print a simple aligned table of the discovered regimes to `fp`. */
void mb_regime_list_print(FILE *fp, const mb_regime_t *r, int n);

#ifdef __cplusplus
}
#endif

#endif /* MB_REGIME_H */
