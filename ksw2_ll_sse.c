#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "ksw2.h"

#if defined(__ARM_NEON)
#include "s2n-lite.h"
#elif defined(__SSE2__)
#include <xmmintrin.h>
#else
#error "Missing SSE2 or NEON intrinsics"
#endif

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

#ifndef Kgrow
#define Kgrow(km, type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = (type*)krealloc(km, ptr, (__m) * sizeof(type)); \
		} \
	} while (0)
#endif

typedef struct {
	int qlen, slen;
	uint8_t shift, mdiff, max, size;
	__m128i *qp, *H0, *H1, *E, *Hmax;
	void *km;
} kswq_t;

/**
 * Initialize the query data structure
 *
 * @param size   Number of bytes used to store a score; valid valures are 1 or 2
 * @param qlen   Length of the query sequence
 * @param query  Query sequence
 * @param m      Size of the alphabet
 * @param mat    Scoring matrix in a one-dimension array
 *
 * @return       Query data structure
 */
void *ksw_ll_qinit(void *km, int size, int qlen, const uint8_t *query, int m, const int8_t *mat)
{
	kswq_t *q;
	int slen, a, tmp, p;

	size = size > 1? 2 : 1;
	p = 8 * (3 - size); // # values per __m128i
	slen = (qlen + p - 1) / p; // segmented length
#if defined(__ARM_NEON)
	q = (kswq_t*)kmalloc(km, sizeof(kswq_t) + 256 + 16 * slen * (m + 6)); // a single block of memory
#else
	q = (kswq_t*)kmalloc(km, sizeof(kswq_t) + 256 + 16 * slen * (m + 4)); // a single block of memory
#endif
	q->qp = (__m128i*)(((size_t)q + sizeof(kswq_t) + 15) >> 4 << 4); // align memory
#if defined(__ARM_NEON)
	q->H0 = q->qp + slen * m; // on arm64 the u8 kernel interleaves E' with H' in H0/H1: {H,E} per stripe
	q->H1 = q->H0 + 2 * slen;
	q->E  = q->H1 + 2 * slen;
#else
	q->H0 = q->qp + slen * m;
	q->H1 = q->H0 + slen;
	q->E  = q->H1 + slen;
#endif
	q->Hmax = q->E + slen;
	q->slen = slen; q->qlen = qlen; q->size = size;
	q->km = km;
	// compute shift
	tmp = m * m;
	for (a = 0, q->shift = 127, q->mdiff = 0; a < tmp; ++a) { // find the minimum and maximum score
		if (mat[a] < (int8_t)q->shift) q->shift = mat[a];
		if (mat[a] > (int8_t)q->mdiff) q->mdiff = mat[a];
	}
	q->max = q->mdiff;
	q->shift = 256 - q->shift; // NB: q->shift is uint8_t
	q->mdiff += q->shift; // this is the difference between the min and max scores
	// An example: p=8, qlen=19, slen=3 and segmentation:
	//  {{0,3,6,9,12,15,18,-1},{1,4,7,10,13,16,-1,-1},{2,5,8,11,14,17,-1,-1}}
	if (size == 1) {
		int8_t *t = (int8_t*)q->qp;
#if defined(__ARM_NEON)
		const int bias = 0; // the arm64 u8 kernel adds the signed profile with USQADD, so no bias
#else
		const int bias = q->shift;
#endif
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? -1 : ma[query[k]]) + bias;
		}
	} else {
		int16_t *t = (int16_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? -1 : ma[query[k]]);
		}
	}
	return q;
}

static inline int ksw_le_u8(__m128i a, __m128i b)
{
#if defined(__ARM_NEON)
	// narrow the 16 compare lanes to a 64-bit mask: cheaper and shorter than a umaxv reduction
	uint8x16_t gt = vcgtq_u8(a, b);
	return vget_lane_u64(vreinterpret_u64_u8(vshrn_n_u16(vreinterpretq_u16_u8(gt), 4)), 0) == 0;
#elif defined(__SSE2__)
	return _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_subs_epu8(a, b), _mm_setzero_si128())) == 0xffff;
#endif
}

#if defined(__ARM_NEON)
// H'(i-1,j-1)+S(i,j) saturated to [0,255] in one step from the SIGNED profile: identical to
// adding the biased profile and subtracting the bias except in cells past 255-shift, and a
// row holding such a cell trips the gmax+shift>=255 exit in both forms with the same result
static inline __m128i ksw_usqadd(__m128i h, __m128i s) { return vsqaddq_u8(h, vreinterpretq_s8_u8(s)); }
// {H',E'} of one stripe with a single two-register store/load
static inline void ksw_st2(__m128i *p, __m128i h, __m128i e) { uint8x16x2_t v; v.val[0] = h; v.val[1] = e; vst1q_u8_x2((uint8_t*)p, v); }
#endif

static inline int ksw_max_u8(__m128i x)
{
#if defined(__ARM_NEON)
	return vmaxvq_u8(x);
#elif defined(__SSE2__)
	x = _mm_max_epu8(x, _mm_srli_si128(x, 8));
	x = _mm_max_epu8(x, _mm_srli_si128(x, 4));
	x = _mm_max_epu8(x, _mm_srli_si128(x, 2));
	x = _mm_max_epu8(x, _mm_srli_si128(x, 1));
	return _mm_extract_epi16(x, 0) & 0x00ff;
#endif
}

ksw_llrst_t ksw_ll_u8_core(void *q_, int tlen, const uint8_t *target, int _gapo, int _gape, int xtra)
{
	kswq_t *q = (kswq_t*)q_;
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m128i gapoe, gape, *H0, *H1, *Hmax;
#if defined(__ARM_NEON)
	__m128i gape2x, hl, *qp = q->qp;
	int gapo_pos = (uint8_t)(_gapo + _gape) > (uint8_t)_gape; // o>0 as the 8-bit vectors see it
	int hmax_pending = 0, qshift = q->shift;
	hl = _mm_setzero_si128(); // H'(i-1,slen-1) carried in a register (H0 starts all zero)
#else
	__m128i shift, *E;
#endif
	ksw_llrst_t r = { 0, -1, -1, -1, -1 };

	// initialization
	minsc = (xtra&KSW_LL_SUBO)? xtra&0xffff : 0x10000;
	endsc = (xtra&KSW_LL_STOP)? xtra&0xffff : 0x10000;
	m_b = n_b = 0; b = 0;
	gapoe = _mm_set1_epi8(_gapo + _gape);
	gape = _mm_set1_epi8(_gape);
#if !defined(__ARM_NEON)
	shift = _mm_set1_epi8(q->shift);
#endif
	H0 = q->H0; H1 = q->H1; Hmax = q->Hmax;
	slen = q->slen;
#if defined(__ARM_NEON)
	memset(H0,   0, 2 * slen * sizeof(__m128i)); // the {H,E} pairs read by the first row
#else
	E = q->E;
	memset(E,    0, slen * sizeof(__m128i));
	memset(H0,   0, slen * sizeof(__m128i));
#endif
	memset(Hmax, 0, slen * sizeof(__m128i));
#if defined(__ARM_NEON)
	{ // two gap extensions, saturated as the two single steps would be
		int gape2 = 2 * (int)(uint8_t)_gape;
		gape2x = _mm_set1_epi8(gape2 < 255? gape2 : 255);
	}
#endif
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, imax;
#if defined(__ARM_NEON)
		__m128i e, h, t, f, max, mf, me, *S = qp + target[i] * slen; // s is the 1st score vector
#else
		__m128i e, h, t, f, max, mf, me, *S = q->qp + target[i] * slen; // s is the 1st score vector
#endif
		f = max = _mm_setzero_si128();
#if defined(__ARM_NEON)
		if (LIKELY(gapo_pos)) {
			/* arm64 fast path, exact for o>0 (the recurrence is spelled out in the x86 loop below):
			 *  - H'(i-1,slen-1), the diagonal input of the first stripe, comes from the register
			 *    hl rather than a store->load round trip (reloaded only if lazy-F touched it)
			 *  - M is one USQADD of the signed profile (see ksw_usqadd)
			 *  - as o>=0, F' reads max(M,E) in place of the reduced H'(i,j): the E'-o-r term
			 *    it drops is dominated by E'-r. this takes H' off the loop-carried f chain
			 *  - two stripes per iteration, with F'(i,j+2) taken directly from F'(i,j):
			 *      F'(i,j+2) = max{F'(i,j)-2r, max(M,E)(i,j)-o-2r, max(M,E)(i,j+1)-o-r}
			 *    bit-identical to the two single steps as saturating subtraction distributes
			 *    over max. F'(i,j+1) is still computed, but off the chain, which this halves
			 *  - H' and E' of a stripe are adjacent ({H,E} at H1+2j) and stored with one
			 *    two-register store, which costs a single vector-pipe slot instead of two,
			 *    except the last stripe of the last pair (below), which stores H' alone so
			 *    hl can be set before E' is computed, then stores E' separately */
			h = _mm_slli_si128(hl, 1); // h=H(i-1,-1); << instead of >> because x64 is little-endian
			j = 0;
			if (slen & 1) { // the first stripe alone when slen is odd
				h = ksw_usqadd(h, _mm_load_si128(S));
				e = _mm_load_si128(H0 + 1);
				me = _mm_max_epu8(h, e);
				h = _mm_max_epu8(me, f);
				max = _mm_max_epu8(max, me);
				e = _mm_subs_epu8(e, gape);
				t = _mm_subs_epu8(h, gapoe);
				e = _mm_max_epu8(e, t);
				_mm_store_si128(H1 + 0, h);
				_mm_store_si128(H1 + 1, e);
				hl = h;
				t = _mm_subs_epu8(me, gapoe);
				f = _mm_subs_epu8(f, gape);
				f = _mm_max_epu8(f, t); // f=F'(i,1)
				h = _mm_load_si128(H0 + 0);
				j = 1;
			}
			for (; LIKELY(j + 2 < slen); j += 2) { // every pair but the last
				__m128i f1, t0;
				// stripe j
				h = ksw_usqadd(h, _mm_load_si128(S + j)); // h=M=H'(i-1,j-1)+S(i,j), saturated to [0,255]
				e = _mm_load_si128(H0 + 2 * j + 1); // e=E'(i,j)
				me = _mm_max_epu8(h, e); // me=max(M,E)
				h = _mm_max_epu8(me, f); // h=H'(i,j)
				max = _mm_max_epu8(max, me); // == max over h: every F' is below an earlier max(M,E)
				e = _mm_subs_epu8(e, gape);
				t = _mm_subs_epu8(h, gapoe); // H'(i,j) - o - r
				e = _mm_max_epu8(e, t); // e=E'(i+1,j)
				ksw_st2(H1 + 2 * j, h, e);
				t0 = _mm_subs_epu8(me, gapoe); // max(M,E)(i,j) - o - r
				f1 = _mm_subs_epu8(f, gape);
				f1 = _mm_max_epu8(f1, t0); // f1=F'(i,j+1)
				h = _mm_load_si128(H0 + 2 * j);
				// stripe j+1
				h = ksw_usqadd(h, _mm_load_si128(S + j + 1));
				e = _mm_load_si128(H0 + 2 * j + 3);
				me = _mm_max_epu8(h, e);
				h = _mm_max_epu8(me, f1); // h=H'(i,j+1)
				max = _mm_max_epu8(max, me);
				e = _mm_subs_epu8(e, gape);
				t = _mm_subs_epu8(h, gapoe);
				e = _mm_max_epu8(e, t); // e=E'(i+1,j+1)
				ksw_st2(H1 + 2 * j + 2, h, e);
				t = _mm_subs_epu8(me, gapoe); // max(M,E)(i,j+1) - o - r
				t0 = _mm_subs_epu8(t0, gape); // max(M,E)(i,j) - o - 2r
				t = _mm_max_epu8(t, t0);
				f = _mm_subs_epu8(f, gape2x);
				f = _mm_max_epu8(f, t); // f=F'(i,j+2)
				h = _mm_load_si128(H0 + 2 * j + 2);
			}
			if (j < slen) { // the last pair: H'(i,slen-1) is stored alone so the next row's first stripe does not wait on E'
				__m128i f1, t0;
				// stripe j
				h = ksw_usqadd(h, _mm_load_si128(S + j)); // h=M=H'(i-1,j-1)+S(i,j), saturated to [0,255]
				e = _mm_load_si128(H0 + 2 * j + 1); // e=E'(i,j)
				me = _mm_max_epu8(h, e); // me=max(M,E)
				h = _mm_max_epu8(me, f); // h=H'(i,j)
				max = _mm_max_epu8(max, me); // == max over h: every F' is below an earlier max(M,E)
				e = _mm_subs_epu8(e, gape);
				t = _mm_subs_epu8(h, gapoe); // H'(i,j) - o - r
				e = _mm_max_epu8(e, t); // e=E'(i+1,j)
				ksw_st2(H1 + 2 * j, h, e);
				t0 = _mm_subs_epu8(me, gapoe); // max(M,E)(i,j) - o - r
				f1 = _mm_subs_epu8(f, gape);
				f1 = _mm_max_epu8(f1, t0); // f1=F'(i,j+1)
				h = _mm_load_si128(H0 + 2 * j);
				// stripe j+1
				h = ksw_usqadd(h, _mm_load_si128(S + j + 1));
				e = _mm_load_si128(H0 + 2 * j + 3);
				me = _mm_max_epu8(h, e);
				h = _mm_max_epu8(me, f1); // h=H'(i,j+1)
				max = _mm_max_epu8(max, me);
				_mm_store_si128(H1 + 2 * j + 2, h);
				hl = h;
				e = _mm_subs_epu8(e, gape);
				t = _mm_subs_epu8(h, gapoe);
				e = _mm_max_epu8(e, t); // e=E'(i+1,j+1)
				_mm_store_si128(H1 + 2 * j + 3, e);
				t = _mm_subs_epu8(me, gapoe); // max(M,E)(i,j+1) - o - r
				t0 = _mm_subs_epu8(t0, gape); // max(M,E)(i,j) - o - 2r
				t = _mm_max_epu8(t, t0);
				f = _mm_subs_epu8(f, gape2x);
				f = _mm_max_epu8(f, t); // f=F'(i,slen)
			}
			// lazy-F loop. the exit test reads H'(i,j)-o-r of the cell BEFORE the update: with
			// o>0, f-r<=max(H',f)-o-r holds iff f-r<=H'-o-r, so it no longer waits on the max
			if (ksw_le_u8(f, _mm_setzero_si128())) goto end_loop_u8; // F'(i,slen)==0: the first lazy step is a no-op and its test holds
			for (k = 0; LIKELY(k < 16); ++k) {
				f = _mm_slli_si128(f, 1);
				// the first stripe alone (the usual exit), then two stripes per test:
				// once the test holds at a stripe, the rest of the pass cannot change H', so
				// testing after a block is exact as long as the pass boundary is respected
				h = _mm_load_si128(H1 + 0);
				t = _mm_subs_epu8(h, gapoe);
				h = _mm_max_epu8(h, f);
				_mm_store_si128(H1 + 0, h);
				f = _mm_subs_epu8(f, gape);
				if (UNLIKELY(ksw_le_u8(f, t))) goto lazy_done_u8;
				for (j = 1; LIKELY(j + 1 < slen); j += 2) {
					uint8x16_t g0, g1;
					h = _mm_load_si128(H1 + 2 * j);
					t = _mm_subs_epu8(h, gapoe);
					h = _mm_max_epu8(h, f);
					_mm_store_si128(H1 + 2 * j, h);
					f = _mm_subs_epu8(f, gape);
					g0 = vcgtq_u8(f, t);
					h = _mm_load_si128(H1 + 2 * j + 2);
					t = _mm_subs_epu8(h, gapoe);
					h = _mm_max_epu8(h, f);
					_mm_store_si128(H1 + 2 * j + 2, h);
					f = _mm_subs_epu8(f, gape);
					g1 = vcgtq_u8(f, t);
					{ // exit if either stripe's test holds: g0==0 or g1==0
						uint64_t m = vget_lane_u64(vreinterpret_u64_u8(vshrn_n_u16(vreinterpretq_u16_u8(vpmaxq_u8(g0, g1)), 4)), 0);
						if (UNLIKELY((uint32_t)m == 0 || (m >> 32) == 0)) goto lazy_done_u8;
					}
				}
				if (j < slen) {
					h = _mm_load_si128(H1 + 2 * j);
					t = _mm_subs_epu8(h, gapoe);
					h = _mm_max_epu8(h, f);
					_mm_store_si128(H1 + 2 * j, h);
					f = _mm_subs_epu8(f, gape);
					if (UNLIKELY(ksw_le_u8(f, t))) goto lazy_done_u8;
				}
			}
lazy_done_u8:
			// the register copy of H'(i,slen-1) is stale once a lazy pass reached the last stripe
			if (k > 0 || j + 3 >= slen - 1) hl = _mm_load_si128(H1 + 2 * (slen - 1)); // j is the start of the block that exited
		} else {
			// o<=0 as the vectors see it: the reassociated single-stripe loop, whose lazy-F test
			// f-r<=max(H',f)-o-r always holds, so that loop stops at the first stripe
			h = _mm_slli_si128(hl, 1);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = ksw_usqadd(h, _mm_load_si128(S + j));
				e = _mm_load_si128(H0 + 2 * j + 1);
				me = _mm_max_epu8(h, e); // me=max(M,E)
				mf = _mm_max_epu8(h, f); // mf=max(M,F)
				h = _mm_max_epu8(me, f); // h=H'(i,j)
				max = _mm_max_epu8(max, me); // == max over h: every F' is below an earlier max(M,E)
				_mm_store_si128(H1 + 2 * j, h);
				e = _mm_subs_epu8(e, gape);
				t = _mm_subs_epu8(mf, gapoe);
				e = _mm_max_epu8(e, t);
				_mm_store_si128(H1 + 2 * j + 1, e);
				f = _mm_subs_epu8(f, gape);
				t = _mm_subs_epu8(me, gapoe);
				f = _mm_max_epu8(f, t);
				h = _mm_load_si128(H0 + 2 * j);
			}
			f = _mm_slli_si128(f, 1);
			_mm_store_si128(H1 + 0, _mm_max_epu8(_mm_load_si128(H1 + 0), f));
			hl = _mm_load_si128(H1 + 2 * (slen - 1));
		}
#else
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 1); // h=H(i-1,-1); << instead of >> because x64 is little-endian
		for (j = 0; LIKELY(j < slen); ++j) {
			/* SW cells are computed in the following order:
			 *   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			 *   E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
			 *   F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
			 */
			// compute H'(i,j); note that at the beginning, h=H'(i-1,j-1)
			h = _mm_adds_epu8(h, _mm_load_si128(S + j));
			h = _mm_subs_epu8(h, shift); // h=M=H'(i-1,j-1)+S(i,j)
			e = _mm_load_si128(E + j); // e=E'(i,j)
			me = mf = h = _mm_max_epu8(_mm_max_epu8(h, e), f); // h=H'(i,j)
			max = _mm_max_epu8(max, h); // set max
			_mm_store_si128(H1 + j, h); // save to H'(i,j)
			// now compute E'(i+1,j)
			e = _mm_subs_epu8(e, gape); // e=E'(i,j) - e_del
			t = _mm_subs_epu8(mf, gapoe); // max(M,F) - o_del - e_del
			e = _mm_max_epu8(e, t); // e=E'(i+1,j)
			_mm_store_si128(E + j, e); // save to E'(i+1,j)
			// now compute F'(i,j+1)
			f = _mm_subs_epu8(f, gape);
			t = _mm_subs_epu8(me, gapoe); // max(M,E) - o_ins - e_ins
			f = _mm_max_epu8(f, t);
			// get H'(i-1,j) and prepare for the next j
			h = _mm_load_si128(H0 + j); // h=H'(i-1,j)
		}
		// NB: we do not need to set E(i,j) as we disallow adjecent insertion and then deletion
		for (k = 0; LIKELY(k < 16); ++k) { // this block mimics SWPS3; NB: H(i,j) updated in the lazy-F loop cannot exceed max
			f = _mm_slli_si128(f, 1);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epu8(h, f); // h=H'(i,j)
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu8(h, gapoe);
				f = _mm_subs_epu8(f, gape);
				if (UNLIKELY(ksw_le_u8(f, h))) goto end_loop_u8;
			}
		}
#endif
end_loop_u8:
		imax = ksw_max_u8(max); // imax is the maximum number in max
		if (imax >= minsc) { // write the b array; this condition adds branching unfornately
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) { // then append
				if (n_b == m_b) Kgrow(q->km, uint64_t, b, n_b, m_b);
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i; // te is the end position on the target
#if defined(__ARM_NEON)
			// keep the H1 vector lazily: a run of rows each beating the last is copied once,
			// when a row fails to beat it (its H1 is then H0) or when the loop ends
			hmax_pending = 1;
			if (gmax + qshift >= 255 || gmax >= endsc) {
				for (j = 0; LIKELY(j < slen); ++j)
					_mm_store_si128(Hmax + j, _mm_load_si128(H1 + 2 * j));
				hmax_pending = 0;
				break;
			}
		} else if (hmax_pending) {
			for (j = 0; LIKELY(j < slen); ++j)
				_mm_store_si128(Hmax + j, _mm_load_si128(H0 + 2 * j));
			hmax_pending = 0;
		}
#else
			for (j = 0; LIKELY(j < slen); ++j) // keep the H1 vector
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax + q->shift >= 255 || gmax >= endsc) break;
		}
#endif
		S = H1; H1 = H0; H0 = S; // swap H0 and H1
	}
#if defined(__ARM_NEON)
	if (hmax_pending) // the last row set the max; after the swap its H1 is H0
		for (i = 0; i < slen; ++i) _mm_store_si128(Hmax + i, _mm_load_si128(H0 + 2 * i));
#endif
	r.score = gmax + q->shift < 255? gmax : 255;
	r.te = te;
	if (r.score != 255) { // get a->qe, the end of query match; find the 2nd best score
		int max = -1, tmp, low, high, qlen = slen * 16;
		uint8_t *t = (uint8_t*)Hmax;
		for (i = 0; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / 16 + i % 16 * slen;
			else if ((int)*t == max && (tmp = i / 16 + i % 16 * slen) < r.qe) r.qe = tmp; 
		if (b) {
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}
	kfree(q->km, b);
	return r;
}

static inline int ksw_le_epi16(__m128i a, __m128i b)
{
#if defined(__ARM_NEON)
	uint16x8_t gt = vcgtq_s16(vreinterpretq_s16_u8(a), vreinterpretq_s16_u8(b));
	return vget_lane_u64(vreinterpret_u64_u8(vshrn_n_u16(gt, 4)), 0) == 0;
#elif defined(__SSE2__)
	return _mm_movemask_epi8(_mm_cmpgt_epi16(a, b)) == 0;
#endif
}

static inline int ksw_max_i16(__m128i x)
{
#if defined(__ARM_NEON)
	return vmaxvq_s16(vreinterpretq_s16_u8(x));
#elif defined(__SSE2__)
	x = _mm_max_epi16(x, _mm_srli_si128(x, 8));
	x = _mm_max_epi16(x, _mm_srli_si128(x, 4));
	x = _mm_max_epi16(x, _mm_srli_si128(x, 2));
	return _mm_extract_epi16(x, 0);
#endif
}

ksw_llrst_t ksw_ll_i16_core(void *q_, int tlen, const uint8_t *target, int _gapo, int _gape, int xtra)
{
	kswq_t *q = (kswq_t*)q_;
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m128i gapoe, gape, *H0, *H1, *E, *Hmax;
#if defined(__ARM_NEON)
	__m128i gape2x, *qp = q->qp;
	int gapo_pos = (uint16_t)(_gapo + _gape) > (uint16_t)_gape; // o>0 as the 16-bit vectors see it
	int hmax_pending = 0;
#endif
	ksw_llrst_t r = { 0, -1, -1, -1, -1 };

	// initialization
	minsc = (xtra&KSW_LL_SUBO)? xtra&0xffff : 0x10000;
	endsc = (xtra&KSW_LL_STOP)? xtra&0xffff : 0x10000;
	m_b = n_b = 0; b = 0;
	gapoe = _mm_set1_epi16(_gapo + _gape);
	gape = _mm_set1_epi16(_gape);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	memset(E,    0, slen * sizeof(__m128i));
	memset(H0,   0, slen * sizeof(__m128i));
	memset(Hmax, 0, slen * sizeof(__m128i));
#if defined(__ARM_NEON)
	{ // two gap extensions, saturated as the two single steps would be
		int gape2 = 2 * (int)(uint16_t)_gape;
		gape2x = _mm_set1_epi16(gape2 < 65535? gape2 : 65535);
	}
#endif
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, imax;
#if defined(__ARM_NEON)
		__m128i e, t, h, f, max, mf, me, *S = qp + target[i] * slen; // s is the 1st score vector
#else
		__m128i e, t, h, f, max, mf, me, *S = q->qp + target[i] * slen; // s is the 1st score vector
#endif
		f = max = _mm_setzero_si128();
#if defined(__ARM_NEON)
		if (LIKELY(gapo_pos)) {
			// the arm64 fast path of ksw_ll_u8_core. every e/f/h/me stays in [0,32767], so the
			// signed max distributes over the unsigned saturating subtraction just as in u8
			h = _mm_load_si128(H0 + slen - 1);
			h = _mm_slli_si128(h, 2);
			for (j = 0; LIKELY(j + 1 < slen); j += 2) {
				__m128i f1, t0;
				h = _mm_adds_epi16(h, _mm_load_si128(S + j)); // h=M (match score)
				e = _mm_load_si128(E + j);
				me = _mm_max_epi16(h, e); // me=max(M,E)
				h = _mm_max_epi16(me, f); // h=H'(i,j)
				max = _mm_max_epi16(max, me); // == max over h: every F' is below an earlier max(M,E)
				_mm_store_si128(H1 + j, h);
				e = _mm_subs_epu16(e, gape);
				t = _mm_subs_epu16(h, gapoe); // H'(i,j) - gapoe
				e = _mm_max_epi16(e, t);
				_mm_store_si128(E + j, e);
				t0 = _mm_subs_epu16(me, gapoe); // max(M,E)(i,j) - gapoe
				f1 = _mm_subs_epu16(f, gape);
				f1 = _mm_max_epi16(f1, t0); // f1=F'(i,j+1)
				h = _mm_load_si128(H0 + j);
				h = _mm_adds_epi16(h, _mm_load_si128(S + j + 1));
				e = _mm_load_si128(E + j + 1);
				me = _mm_max_epi16(h, e);
				h = _mm_max_epi16(me, f1); // h=H'(i,j+1)
				max = _mm_max_epi16(max, me);
				_mm_store_si128(H1 + j + 1, h);
				e = _mm_subs_epu16(e, gape);
				t = _mm_subs_epu16(h, gapoe);
				e = _mm_max_epi16(e, t);
				_mm_store_si128(E + j + 1, e);
				t = _mm_subs_epu16(me, gapoe); // max(M,E)(i,j+1) - gapoe
				t0 = _mm_subs_epu16(t0, gape); // max(M,E)(i,j) - gapoe - gape
				t = _mm_max_epi16(t, t0);
				f = _mm_subs_epu16(f, gape2x);
				f = _mm_max_epi16(f, t); // f=F'(i,j+2)
				h = _mm_load_si128(H0 + j + 1);
			}
			if (j < slen) { // the last stripe when slen is odd
				h = _mm_adds_epi16(h, _mm_load_si128(S + j));
				e = _mm_load_si128(E + j);
				me = _mm_max_epi16(h, e);
				h = _mm_max_epi16(me, f);
				max = _mm_max_epi16(max, me);
				_mm_store_si128(H1 + j, h);
				e = _mm_subs_epu16(e, gape);
				t = _mm_subs_epu16(h, gapoe);
				e = _mm_max_epi16(e, t);
				_mm_store_si128(E + j, e);
				t = _mm_subs_epu16(me, gapoe);
				f = _mm_subs_epu16(f, gape);
				f = _mm_max_epi16(f, t);
			}
			if (ksw_le_epi16(f, _mm_setzero_si128())) goto end_loop_i16; // F'(i,slen)==0: the first lazy step is a no-op and its test holds
			for (k = 0; LIKELY(k < 16); ++k) { // exit test from the cell before the update; see ksw_ll_u8_core
				f = _mm_slli_si128(f, 2);
				for (j = 0; LIKELY(j < slen); ++j) {
					h = _mm_load_si128(H1 + j);
					t = _mm_subs_epu16(h, gapoe);
					h = _mm_max_epi16(h, f);
					_mm_store_si128(H1 + j, h);
					f = _mm_subs_epu16(f, gape);
					if (UNLIKELY(ksw_le_epi16(f, t))) goto end_loop_i16;
				}
			}
		} else { // o<=0 as the vectors see it; see ksw_ll_u8_core
			h = _mm_load_si128(H0 + slen - 1);
			h = _mm_slli_si128(h, 2);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_adds_epi16(h, _mm_load_si128(S + j)); // h=M (match score)
				e = _mm_load_si128(E + j);
				me = _mm_max_epi16(h, e); // me=max(M,E)
				mf = _mm_max_epi16(h, f); // mf=max(M,F)
				h = _mm_max_epi16(me, f); // h=H'(i,j)
				max = _mm_max_epi16(max, h);
				_mm_store_si128(H1 + j, h);
				e = _mm_subs_epu16(e, gape);
				t = _mm_subs_epu16(mf, gapoe);
				e = _mm_max_epi16(e, t);
				_mm_store_si128(E + j, e);
				f = _mm_subs_epu16(f, gape);
				t = _mm_subs_epu16(me, gapoe);
				f = _mm_max_epi16(f, t);
				h = _mm_load_si128(H0 + j);
			}
			f = _mm_slli_si128(f, 2);
			_mm_store_si128(H1, _mm_max_epi16(_mm_load_si128(H1), f));
		}
#else
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 2);
		for (j = 0; LIKELY(j < slen); ++j) {
			h = _mm_adds_epi16(h, _mm_load_si128(S++)); // h=M (match score)
			e = _mm_load_si128(E + j);
			me = mf = h = _mm_max_epi16(_mm_max_epi16(h, e), f); // h=H'(i,j)
			max = _mm_max_epi16(max, h);
			_mm_store_si128(H1 + j, h);
			e = _mm_subs_epu16(e, gape);
			t = _mm_subs_epu16(mf, gapoe); // max(M,F) - gapoe
			e = _mm_max_epi16(e, t);
			_mm_store_si128(E + j, e);
			f = _mm_subs_epu16(f, gape);
			t = _mm_subs_epu16(me, gapoe); // max(M,E) - gapoe
			f = _mm_max_epi16(f, t);
			h = _mm_load_si128(H0 + j);
		}
		for (k = 0; LIKELY(k < 16); ++k) {
			f = _mm_slli_si128(f, 2);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epi16(h, f);
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu16(h, gapoe);
				f = _mm_subs_epu16(f, gape);
				if (UNLIKELY(ksw_le_epi16(f, h))) goto end_loop_i16;
			}
		}
#endif
end_loop_i16:
		imax = ksw_max_i16(max);
		if (imax >= minsc) {
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) {
				if (n_b == m_b) Kgrow(q->km, uint64_t, b, n_b, m_b);
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i;
#if defined(__ARM_NEON)
			hmax_pending = 1; // see ksw_ll_u8_core
			if (gmax >= endsc) {
				for (j = 0; LIKELY(j < slen); ++j)
					_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
				hmax_pending = 0;
				break;
			}
		} else if (hmax_pending) {
			for (j = 0; LIKELY(j < slen); ++j)
				_mm_store_si128(Hmax + j, _mm_load_si128(H0 + j));
			hmax_pending = 0;
		}
#else
			for (j = 0; LIKELY(j < slen); ++j)
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax >= endsc) break;
		}
#endif
		S = H1; H1 = H0; H0 = S;
	}
#if defined(__ARM_NEON)
	if (hmax_pending)
		for (i = 0; i < slen; ++i) _mm_store_si128(Hmax + i, _mm_load_si128(H0 + i));
#endif
	r.score = gmax; r.te = te;
	{
		int max = -1, tmp, low, high, qlen = slen * 8;
		uint16_t *t = (uint16_t*)Hmax;
		for (i = 0, r.qe = -1; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / 8 + i % 8 * slen;
			else if ((int)*t == max && (tmp = i / 8 + i % 8 * slen) < r.qe) r.qe = tmp;
		if (b) {
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}
	kfree(q->km, b);
	return r;
}

int ksw_ll_i16(void *q_, int tlen, const uint8_t *target, int _gapo, int _gape, int *qe, int *te)
{
	ksw_llrst_t r;
	r = ksw_ll_i16_core(q_, tlen, target, _gapo, _gape, 0);
	*qe = r.qe, *te = r.te;
	return r.score;
}
