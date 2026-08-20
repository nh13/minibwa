/* b2idx — minibwa seeding backed by the bwa-mem2/bwa-mem3 FM-index.
 *
 * THE PROTOTYPE. Swaps minibwa's native .mbw FM-index for bwa-mem2's
 * .bwt.2bit.64 checkpoint-occ index for the SEEDING step only (SMEM search +
 * suffix-array lookup). Everything downstream — chaining, extension, pairing,
 * SAM — is unchanged and keeps using minibwa's .l2b reference layer, because
 * the two indexes share the same 2N (forward + reverse-complement) SA
 * coordinate space (validated: resolved positions matched 99.25% incl. rbeg).
 *
 * The four entry points mirror the minibwa FM-index seeding API (bwt.h) that
 * seed.c calls, so the swap is a branch on idx->b2:
 *   mb_bwt_cache      <-> (built inside mb_b2_load)
 *   mb_bwt_smem_batch <-> mb_b2_smem_batch   (per-entry min_len/min_occ/st/en)
 *   mb_bwt_sa_batch   <-> mb_b2_sa_batch     (SA index -> 2N position, in place)
 *
 * Implementation is C++ (bwa-mem3's FMI_search is C++); this header is the
 * C-callable seam included by seed.c / map-algo.c. */
#ifndef MB_B2IDX_H
#define MB_B2IDX_H

#include "bwt.h"   /* mb_sai_t, mb_sai_v, mb_smem_entry_t (minibwa types) */

#ifdef __cplusplus
extern "C" {
#endif

/* Load the bwa-mem2 FM-index at `prefix` (.bwt.2bit.64 + .amb/.ann/.pac) and
 * build the 10-mer SMEM cache over its cp_occ. Returns an opaque handle, or
 * NULL on failure. */
void *mb_b2_load(const char *prefix);

/* Free the handle (deletes the FMI_search, which frees cp_occ / sa_*). */
void  mb_b2_destroy(void *b2);

/* Batched SMEM collection: exact structural port of mb_bwt_smem_batch, honoring
 * each entry's min_len / min_occ / st / en, appending mb_sai_t intervals to
 * a[i].v via the caller's kalloc pool `km`. Intervals are in the bwa-mem2 SA
 * space (x[0]/x[1] are bwa-mem2 SA-interval bounds). */
void  mb_b2_smem_batch(void *km, void *b2, int32_t n, mb_smem_entry_t *a);

/* Resolve n suffix-array indices in place: a[i] (bwa-mem2 SA index) -> a[i]
 * (position in the 2N concatenated reference), matching mb_bwt_sa_batch. Uses
 * bwa-mem2's prefetch-pooled get_sa_entries_prefetch. */
void  mb_b2_sa_batch(void *km, void *b2, int64_t n, uint64_t *a);

/* Cheaply read the SA sampling rate of the cp_occ index at `prefix`
 * (<prefix>.bwt.2bit.64) WITHOUT doing a full mb_b2_load. Returns the SA
 * interval (1<<sa_compx, e.g. 8 or 16), or -1 if the file is absent. Used by
 * regime discovery to label/rank a cp_occ regime without paying the cost of
 * loading it. */
int mb_b2_peek_sa_intv(const char *prefix);

#ifdef __cplusplus
}
#endif
#endif
