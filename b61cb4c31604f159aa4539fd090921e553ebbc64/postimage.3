#ifndef MINIBWA_H
#define MINIBWA_H

#include <stdint.h>

#include "l2bit.h"

#define MB_VERSION "0.7-r421"

#define MB_F_PAF              (0x1LL)       // output in the PAF format
#define MB_F_NO_UNMAP         (0x2LL)       // output unmapped query sequences
#define MB_F_COPY_COMMENT     (0x4LL)       // copy FASTX comments to output
#define MB_F_PE               (0x8LL)       // paired-end mode
#define MB_F_LONG             (0x10LL)      // long-sequence mode
#define MB_F_EQX              (0x20LL)      // = in CIGAR
#define MB_F_NO_KALLOC        (0x40LL)      // disable kalloc
#define MB_F_NO_ALN           (0x80LL)      // skip base alignment
#define MB_F_PE_PREDEF        (0x100LL)     // use predefined PE
#define MB_F_WRITE_DS         (0x200LL)     // write ds:Z
#define MB_F_WRITE_CS         (0x400LL)     // write cs:Z
#define MB_F_WRITE_MD         (0x800LL)     // write MD:Z
#define MB_F_2ND_SEQ          (0x1000LL)    // in SAM, write SEQ for secondary alignments
#define MB_F_SUPP_SOFT        (0x2000LL)    // in SAM, use soft-clips for supplementary alignments
#define MB_F_ADAP             (0x4000LL)    // adaptive mode
#define MB_F_PRIMARY5         (0x8000LL)    // for Hi-C
#define MB_F_NO_PAIRING       (0x10000LL)   // don't pair reads
#define MB_F_METH             (0x20000LL)   // methylation mode
#define MB_F_ALT_RECORDS      (0x40000LL)   // emit ALT-contig hits with full SEQ

#define MB_CIGAR_MATCH      0
#define MB_CIGAR_INS        1
#define MB_CIGAR_DEL        2
#define MB_CIGAR_N_SKIP     3
#define MB_CIGAR_SOFTCLIP   4
#define MB_CIGAR_HARDCLIP   5
#define MB_CIGAR_PADDING    6
#define MB_CIGAR_EQ_MATCH   7
#define MB_CIGAR_X_MISMATCH 8

#define MB_CIGAR_STR  "MIDNSHP=XB"

typedef struct {
	uint64_t flag;
	// seeding options
	int32_t min_len; // min seed length
	int32_t max_sub_occ; // look for shorter seed if smem occ below this value
	int32_t max_occ; // max interval occurrence
	// general algorithm options
	int32_t bw, bw_long; // bandwidth
	int32_t max_gap; // break a chain if there are no seeds in a max_gap window
	int32_t max_sr_len; // in the adaptive sr mode, treat reads longer than this as long reads
	// chaining options
	int32_t max_chain_skip;
	int32_t max_chain_iter;
	int32_t min_chain_score; // min chaining score
	float chain_gap_scale;
	// hit processing options
	float mask_level;
	int32_t mask_len;
	float pri_ratio;
	int32_t best_n;
	// alignment options
	int32_t a, b;     // match, mismatch
	int32_t b_ts;     // transition mismatch
	int32_t b_ambi;   // ambiguous mismatch
	int32_t q, q2;    // gap open, long gap open
	int32_t e, e2;    // gap extension, long gap extension
	int32_t end_bonus;
	int32_t min_dp_max; // min_dp_max*a is the min score
	int32_t zdrop;
	int32_t zdrop_inv;
	int32_t min_ksw_len;
	// pairing options
	int32_t max_pe_ins;
	int32_t max_rescue;
	int32_t pen_unpair;
	int32_t pe_avg, pe_std, pe_lo, pe_hi;
	// input/output options
	int32_t sb_len;   // number of bases for batch smem
	int32_t sb_seq;   // number of sequences for batch smem
	int32_t n_thread; // number of worker threads, excluding I/O threads
	int32_t out_n;    // max number of secondary alignments to output
	float out_s;
	int32_t seed;
	int32_t xa_max;
	int64_t mb_size;  // mini-batch size
	int64_t max_mb_size;
	int64_t max_sw_mat;
	int64_t cap_kalloc;
	int32_t lift_tol;  // ALT liftover-group co-location tolerance in bp (default MB_LIFT_TOL)
} mb_opt_t;

struct mb_idx_s;
typedef struct mb_idx_s mb_idx_t;

typedef struct {
	uint32_t cap;               // the capacity of cigar[]
	int32_t dp_score, dp_max0;  // DP score; score of the max-scoring segment
	int32_t dp_max, dp_max2;    // adjusted score and second best score for mapQ
	uint32_t n_ambi:31, cs:1;   // number of ambiguous bases;
	int32_t n_cigar;            // number of cigar operations in cigar[]
	uint32_t cigar[];           // cs/MD is appended at the end
} mb_extra_t;

#define MB_PARENT_UNSET   (-1)
#define MB_PARENT_TMP_PRI (-2)

typedef struct {
	int64_t tid;            // target ID (the original tid, NOT stranded)
	int64_t ts, te;         // target start and end
	int32_t id;             // ID for internal uses
	int32_t cnt;            // number of anchors
	int32_t score, score0;  // chaining score; score0 is the original chaining score
	int32_t as;             // offset in the a[] array (for internal uses only)
	int32_t qs, qe;         // query start and end
	int32_t parent, n_sub, subsc;
	int32_t mlen, blen;
	int32_t mapq;
	uint32_t hash;
	uint32_t rev:1, proper_pair:1, sam_pri:1, flt:1, inv:1, split:2, split_inv:1, rescued:1, frac_high:8, seed_ratio:8, is_alt:1, dummy:6;
	mb_extra_t *p;
} mb_hit_t;

/* Co-location tolerance (bp) on the lifted footprint start used to group hits
 * by primary locus.  Two hits are the same locus iff they share pri_tid, rev,
 * and their lifted_st differ by no more than the tolerance.  This is the
 * COMPILE-TIME DEFAULT; the effective value is opt->lift_tol (runtime-tunable
 * via --alt-lift-tol), threaded into the survival guard, reconciliation, and the
 * PE demotion guard.  Raise it for .alt files whose ALT-to-primary CIGARs carry
 * larger indels (lifting a read start across an indel can drift lifted_st by up
 * to the indel size); the query-span-overlap requirement guards against merging
 * genuine paralogs even at a looser tolerance. */
#define MB_LIFT_TOL 10

/* Maximum number of per-block lifted sub-placements retained in mb_place_t.
 * A read footprint (~150 bp pre/post DP) overlaps at most this many .alt lift
 * blocks in practice; if MORE blocks overlap (a pathologically fragmented .alt
 * CIGAR over the footprint) the FIRST MB_MAX_SUBPL are kept and the rest spill
 * (documented in mb_hit_place).  Spilling only ever DROPS candidate co-location
 * intervals -- it can never invent a spurious match -- so it is conservatively
 * safe for paralog isolation. */
#define MB_MAX_SUBPL 8

/* One lifted sub-placement: where a single overlapping .alt lift block maps the
 * footprint onto primary coordinates.
 *   st       representative (min) primary coordinate of this block's lifted span
 *   pri_tid  primary contig this block lands on
 *   rev      strand of this block's footprint on primary (.alt block strand XOR h->rev) */
typedef struct {
	int64_t st;
	int64_t pri_tid;
	uint8_t rev;
} mb_subpl_t;

/* The lifted PLACEMENT of one hit: the primary footprint it occupies over its
 * liftable portion.  Computed by mb_hit_place().
 *
 * MULTI-INTERVAL placement (SV-breakpoint-aware grouping): instead of collapsing
 * every overlapping .alt lift block into a single [lifted_st, lifted_en], the
 * placement records ONE sub-placement per overlapping block in subpl[].  A
 * breakpoint-spanning ALT hit whose footprint straddles an SV-scale indel then
 * exposes BOTH the near-breakpoint primary position AND the far one as separate
 * sub-placements, so it can still co-locate with its primary twin via the
 * matching sub-interval (mb_places_colocate) instead of being dragged thousands
 * of bp away by a min/max collapse.  Co-location requires a SHARED sub-interval,
 * so two distinct primary loci that happen to land in one inflated span are NOT
 * merged (paralog safety).
 *
 *   pri_tid    REPRESENTATIVE primary contig (== subpl[0].pri_tid; == h->tid for
 *              non-ALT hits).  Kept for back-compat readers.
 *   lifted_st  REPRESENTATIVE primary coordinate (== subpl[0].st).  Back-compat
 *              grouping key for any reader not yet on the multi-interval API.
 *   lifted_en  max primary coordinate over all sub-placements (cosmetic: nothing
 *              reads it for grouping decisions).
 *   rev        REPRESENTATIVE strand (== subpl[0].rev; .alt block strand XOR h->rev).
 *   liftable   1 iff n_subpl >= 1 (at least one block of the footprint lifts);
 *              0 iff the ENTIRE footprint falls in holes (ALT-specific -> own group).
 *   n_subpl    number of valid sub-placements (1..MB_MAX_SUBPL; 0 when !liftable).
 *   subpl      the per-block lifted sub-placements (first n_subpl entries valid). */
typedef struct {
	int64_t pri_tid;
	int64_t lifted_st, lifted_en;
	uint8_t rev, liftable;
	int n_subpl;
	mb_subpl_t subpl[MB_MAX_SUBPL];
} mb_place_t;

struct mb_tbuf_s;
typedef struct mb_tbuf_s mb_tbuf_t;

#ifdef __cplusplus
extern "C" {
#endif

mb_idx_t *mb_idx_load(const char *prefix, int32_t is_meth);
mb_idx_t *mb_idx_load_mmap(const char *prefix, int32_t is_meth, int preload);
void mb_idx_destroy(mb_idx_t *idx);
void mb_idx_set_alt(mb_idx_t *idx, const char *fn);
const char *mb_idx_ctg_name(const mb_idx_t *idx, int32_t tid);
int64_t mb_idx_ctg_len(const mb_idx_t *idx, int32_t tid);

/**
 * Compute the lifted PLACEMENT (primary footprint) of one hit.
 *
 * For a non-ALT hit this is the identity placement on its own contig.  For an
 * ALT hit it lifts the aligned footprint (the chain interval [ts,te) pre-DP, or
 * the exact CIGAR span post-DP) through the .alt span-lift to primary
 * coordinates, taking the min/max over the LIFTED (primary) outputs so a reverse
 * .alt block folds correctly and a footprint end sitting in a hole walks inward
 * to the first/last liftable base.  See mb_place_t.
 *
 * @param l2b  the span-lift index (must have .alt loaded for ALT hits)
 * @param h    the hit; h->p may be NULL (pre-DP, coarse) or non-NULL (post-DP, exact)
 * @return     the placement; liftable==0 if the whole footprint is in holes
 */
mb_place_t mb_hit_place(const l2b_t *l2b, const mb_hit_t *h);

void mb_opt_init(mb_opt_t *opt);
int mb_opt_preset(mb_opt_t *opt, const char *preset);

mb_tbuf_t *mb_tbuf_init(int no_kalloc);
void mb_tbuf_destroy(mb_tbuf_t *b);
int32_t mb_tbuf_reset(mb_tbuf_t *b, int64_t max_block_size);

/**
 * Align one sequence
 *
 * @param opt        options, typically initialized by mb_opt_init()
 * @param idx        index
 * @param qlen       query length
 * @param seq        query sequence, ASCII or 01/2/3 encoded
 * @param mt         methylation type: 0 for unmethylated, 1 for read1 (C-to-T) and 2 for read2 (G-to-A)
 * @param n_hit      (out) number of hits
 * @param b          thread buffer; can be NULL
 * @param qname      query name
 *
 * @return hit array
 */
mb_hit_t *mb_map(const mb_opt_t *opt, const mb_idx_t *idx, int32_t qlen, const char *seq, int32_t mt, int32_t *n_hit, mb_tbuf_t *b, const char *qname);

/**
 * Align a set of sequences in batch
 *
 * @param opt        options, typically initialized by mb_opt_init()
 * @param idx        index
 * @param n_seq      number of sequences
 * @param qlen       query lengths, of size n_seq
 * @param seq        query sequences, ASCII or 01/2/3 encoded, of size n_seq
 * @param n_hit      (out) number of hits, of size n_seq
 * @param b          thread buffer; can be NULL
 * @param qname      query name, of size n_seq
 *
 * @return hits, of size n_seq
 */
mb_hit_t **mb_map_batch(const mb_opt_t *opt, const mb_idx_t *idx, int32_t n_seq, const int32_t *qlen, const char **seq, int32_t *n_hit, mb_tbuf_t *b, const char **qname);

#ifdef __cplusplus
}
#endif

#endif
