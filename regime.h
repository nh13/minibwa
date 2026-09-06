#ifndef MB_REGIME_H
#define MB_REGIME_H

#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque forward declaration: avoids pulling the public minibwa.h (and its
 * whole API surface) into this internal header just to spell the return
 * type of mb_idx_load_regime() below. Callers get the real definition of
 * mb_idx_t from minibwa.h, which they already include. */
struct mb_idx_s;

/* Which SA-lookup backend a regime uses: minibwa's native classic-BWT (bwt.c)
 * or bwa-mem2/bwa-mem3's checkpointed cp_occ FM-index (b2idx.cpp; compile-
 * optional, requires MB_HAVE_B2). */
typedef enum { MB_BACKEND_BWT, MB_BACKEND_CP_OCC } mb_backend_t;

typedef struct {
	mb_backend_t backend;
	int sa_bit;
	uint64_t est_ram;
	int speed_rank;
	uint32_t mode_mask;
	int is_meth;        /* 1 = methylated index (<prefix>.meth.mbw); 0 = normal */
	char name[16];      /* "sa16","sa8","sa8-b2" */
	char sa_path[1024]; /* "" = bundled in .mbw */
} mb_regime_t;

#define MB_MODE_SRPE 1u
#define MB_MODE_METH 2u
#define MB_MODE_HIC  4u
#define MB_MODE_LR   8u

/* Discover which SA-density regimes exist on disk for `prefix`. For a normal
 * index (is_meth==0) that is `<prefix>.l2b` and `<prefix>.mbw` plus any
 * `<prefix>.sa.u*` sidecars; for a methylated index (is_meth!=0) it is
 * `<prefix>.l2b` and `<prefix>.meth.mbw` (a single bundled regime, no sidecars).
 * Writes up to `max` regimes into `out` and returns the count (0 if the base
 * index files are missing). When `b2_available` and a co-located bwa-mem3
 * cp_occ index is present, also discovers one MB_BACKEND_CP_OCC regime for it
 * (compile-optional: requires MB_HAVE_B2). */
int mb_regime_discover(const char *prefix, int is_meth, int b2_available, mb_regime_t *out, int max);

/* Detect the usable memory budget: min(host available, cgroup limit,
 * user_cap_bytes). A zero cap/limit/host value means "no such constraint"
 * and is ignored when computing the min. Pass 0 for user_cap_bytes to
 * ignore the user cap entirely. */
uint64_t mb_mem_budget(uint64_t user_cap_bytes);

/* Pick the fastest regime (by speed_rank) that is eligible for `mode` and
 * fits within `budget` bytes. When `use_mmap` is set, a bundled-SA regime is
 * demand-paged (near-zero resident RAM) and so is exempt from `budget`; a
 * sidecar-SA regime still heap-loads its SA and is gated by `budget` even under
 * mmap. If `forced` is non-NULL, return the index of the regime whose name
 * matches it exactly (ignoring budget), or -1 if no such eligible regime
 * exists. Returns -1 if nothing fits/matches. */
int mb_regime_pick(const mb_regime_t *r, int n, uint64_t budget, uint32_t mode, int use_mmap, const char *forced);

/* Print a simple aligned table of the discovered regimes to `fp`. */
void mb_regime_list_print(FILE *fp, const mb_regime_t *r, int n);

/* Load an index for a specific SA regime `rg` (see mb_regime_discover above).
 * `prefix` is the index prefix (as passed to `minibwa index`); `use_mmap`/
 * `preload` behave as in mb_idx_load_mmap(). Returns NULL on any failure.
 * Defined in map-algo.c; declared here (not in the public minibwa.h) since
 * it takes an mb_regime_t, an internal type callers outside the aligner
 * have no business constructing. */
struct mb_idx_s *mb_idx_load_regime(const char *prefix, const mb_regime_t *rg, int use_mmap, int preload);

#ifdef __cplusplus
}
#endif

#endif /* MB_REGIME_H */
