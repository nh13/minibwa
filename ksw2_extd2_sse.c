#include <string.h>
#include <assert.h>
#include "ksw2.h"

#if defined(__ARM_NEON)
#define __SSE4_1__
#include "s2n-lite.h"
#elif defined(__SSE4_1__)
#include <smmintrin.h>
#elif defined(__SSE2__)
#include <xmmintrin.h>
#else
#error "Missing SSE2 or NEON intrinsics"
#endif

/* KSW_EXTD2_NEON_KERNEL selects one of two DP rail layouts. Both produce the
 * same scores and CIGARs; only the speed differs, and it differs by ISA.
 *
 * 1 (arm64): x, v and x2 are double-buffered by row parity, so the t-1
 *   shifted rails are plain unaligned loads from the previous row's buffer
 *   (the boundary byte is inserted with ksw_insert0) instead of one vext per
 *   rail per block. The rails themselves store the SAME unbiased values as
 *   the original kernel below -- no bias is added or subtracted anywhere in
 *   this path -- so the two kernels are byte-identical by construction over
 *   the full scoring-parameter domain, not just realistic ones. This
 *   measures +4-9%.
 *
 * 0 (x86 and everything else): the original single-buffer kernel, where the
 *   shift is one _mm_alignr_epi8 per rail. The double buffers' extra memory
 *   traffic does not pay there (0.97-0.98x at short reads), so x86 keeps the
 *   original code as it was and inherits its output and speed by construction. */
#ifndef KSW_EXTD2_NEON_KERNEL
#if defined(__ARM_NEON)
#define KSW_EXTD2_NEON_KERNEL 1
#else
#define KSW_EXTD2_NEON_KERNEL 0
#endif
#endif

/* d with the k bits taken from g. On arm64 this is one BIT; clang folds the
 * portable spelling back to and/orr, so the instruction is written out. The flag
 * form needs d's k bits already clear; the sel form is a full blend. */
#if defined(__ARM_NEON)
static inline __m128i ksw_bit(__m128i d, __m128i g, __m128i k)
{
	__asm__("bit %0.16b, %1.16b, %2.16b" : "+w"(d) : "w"(g), "w"(k));
	return d;
}
static inline __m128i ksw_bitins_flag(__m128i d, __m128i g, __m128i k) { return ksw_bit(d, g, k); }
static inline __m128i ksw_bitsel(__m128i d, __m128i g, __m128i k) { return ksw_bit(d, g, k); }
#else
static inline __m128i ksw_bitins_flag(__m128i d, __m128i g, __m128i k) { return _mm_or_si128(d, _mm_and_si128(g, k)); }
static inline __m128i ksw_bitsel(__m128i d, __m128i g, __m128i k) { return _mm_blendv_epi8(d, g, k); }
#endif

/* not a general _mm_shuffle_epi8: vqtbl1q_u8 zeroes any index >= 16 where SSSE3
 * only zeroes on bit 7. Every index built below is <= 12, so the two agree. */
#if defined(__ARM_NEON)
static inline __m128i ksw_shuffle_epi8(__m128i a, __m128i b) { return vqtbl1q_u8(a, b); }
static inline __m128i ksw_xor_si128(__m128i a, __m128i b) { return veorq_u8(a, b); }
#else
static inline __m128i ksw_shuffle_epi8(__m128i a, __m128i b) { return _mm_shuffle_epi8(a, b); }
static inline __m128i ksw_xor_si128(__m128i a, __m128i b) { return _mm_xor_si128(a, b); }
#endif

/* a compare mask (0x00/0xff) to 0/1; on arm64 without a constant register: |-1| = 1 */
#if defined(__ARM_NEON)
static inline __m128i ksw_mask01(__m128i m) { return vreinterpretq_u8_s8(vabsq_s8(vreinterpretq_s8_u8(m))); }
#else
static inline __m128i ksw_mask01(__m128i m) { return _mm_and_si128(m, _mm_set1_epi8(1)); }
#endif

#if KSW_EXTD2_NEON_KERNEL

/* v with lane 0 replaced by b: the row boundary byte goes into the first
 * shifted rail vector. Only compiled under KSW_EXTD2_NEON_KERNEL, which
 * implies __ARM_NEON (see the default above), so this is NEON-only. */
static inline __m128i ksw_insert0(__m128i v, int8_t b) { return vreinterpretq_u8_s8(vsetq_lane_s8(b, vreinterpretq_s8_u8(v), 0)); }

/* sign-extend 16 int8 to four int32 vectors; the i16->i32 step folds into the
 * accumulate */
static inline void ksw_widen_i8x16_pair(const int8_t *p, __m128i *w)
{
#if defined(__ARM_NEON)
	int8x16_t b = vld1q_s8(p);
	int16x8_t lo = vmovl_s8(vget_low_s8(b)), hi = vmovl_s8(vget_high_s8(b));
	w[0] = vreinterpretq_u8_s32(vmovl_s16(vget_low_s16(lo)));
	w[1] = vreinterpretq_u8_s32(vmovl_s16(vget_high_s16(lo)));
	w[2] = vreinterpretq_u8_s32(vmovl_s16(vget_low_s16(hi)));
	w[3] = vreinterpretq_u8_s32(vmovl_s16(vget_high_s16(hi)));
#elif defined(__SSE4_1__)
	__m128i b = _mm_loadu_si128((const __m128i*)p);
	w[0] = _mm_cvtepi8_epi32(b);
	w[1] = _mm_cvtepi8_epi32(_mm_srli_si128(b,  4));
	w[2] = _mm_cvtepi8_epi32(_mm_srli_si128(b,  8));
	w[3] = _mm_cvtepi8_epi32(_mm_srli_si128(b, 12));
#else
	int k;
	for (k = 0; k < 4; ++k) w[k] = ksw_i8x4_to_i32x4(p + k * 4);
#endif
}

#else /* the original helpers */

static inline __m128i ksw_alignr15(__m128i cur, __m128i prev) /* {prev[15], cur[0..14]} */
{
#if defined(__ARM_NEON)
	return vextq_u8(prev, cur, 15);
#elif defined(__SSSE3__) || defined(__SSE4_1__)
	return _mm_alignr_epi8(cur, prev, 15);
#else
	return _mm_or_si128(_mm_slli_si128(cur, 1), _mm_srli_si128(prev, 15)); /* only prev[15] survives */
#endif
}

/* sign-extend 16 int8 to four int32 vectors; the i16->i32 step folds into the accumulate */
static inline void ksw_widen_i8x16_pair(const int8_t *p, __m128i *w)
{
#if defined(__ARM_NEON)
	int8x16_t b = vld1q_s8(p);
	int16x8_t lo = vmovl_s8(vget_low_s8(b)), hi = vmovl_s8(vget_high_s8(b));
	w[0] = vreinterpretq_u8_s32(vmovl_s16(vget_low_s16(lo)));
	w[1] = vreinterpretq_u8_s32(vmovl_s16(vget_high_s16(lo)));
	w[2] = vreinterpretq_u8_s32(vmovl_s16(vget_low_s16(hi)));
	w[3] = vreinterpretq_u8_s32(vmovl_s16(vget_high_s16(hi)));
#elif defined(__SSE4_1__)
	__m128i b = _mm_loadu_si128((const __m128i*)p);
	w[0] = _mm_cvtepi8_epi32(b);
	w[1] = _mm_cvtepi8_epi32(_mm_srli_si128(b,  4));
	w[2] = _mm_cvtepi8_epi32(_mm_srli_si128(b,  8));
	w[3] = _mm_cvtepi8_epi32(_mm_srli_si128(b, 12));
#else
	int k;
	for (k = 0; k < 4; ++k) w[k] = ksw_i8x4_to_i32x4(p + k * 4);
#endif
}

#endif /* KSW_EXTD2_NEON_KERNEL */

static inline __m128i ksw_i8x4_to_i32x4(const int8_t *x)
{
#if defined(__ARM_NEON)
    return vreinterpretq_u8_s32(vmovl_s16(vget_low_s16(vmovl_s8(vcreate_s8((uint64_t)(uint32_t)*(int32_t*)x)))));
#elif defined(__SSE4_1__)
	return _mm_cvtepi8_epi32(_mm_cvtsi32_si128(*(int32_t*)x));
#else
	return _mm_setr_epi32(x[0], x[1], x[2], x[3]);
#endif
}

void ksw_extd2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
// The two kernels share the loops below and differ in these macros: how the
// t-1 shifted rails xt1/vt1/x2t1 are produced (block1/block3), where u and v are
// stored (block2), what the four rail stores subtract (__dp_store_*), and how
// the scalar readers of u8[]/v8[] undo any bias (__dp_u8/__dp_v8 and the
// widening of v8[] into H[]).
#if KSW_EXTD2_NEON_KERNEL
// x, v and x2 are double-buffered by row parity: xo/vo/x2o hold row r-1 and are
// only read, xn/vn/x2n receive row r. The t-1 shifted rails xt1, vt1 and x2t1
// are then plain unaligned loads from the row r-1 buffers: block3 loads the
// next block's at the end of an iteration (the loop is entered with the first
// block's, boundary byte inserted), so nothing is shifted or copied in the loop.
#define __dp_code_block1 \
	z = _mm_load_si128(&s[t]); \
	a = _mm_add_epi8(xt1, vt1);                      /* a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14] */ \
	ut = _mm_load_si128(&u[t]);                      /* ut <- u[t..t+15] */ \
	b = _mm_add_epi8(_mm_load_si128(&y[t]), ut);     /* b <- y[r-1][t..t+15] + u[r-1][t..t+15] */ \
	a2= _mm_add_epi8(x2t1, vt1); \
	b2= _mm_add_epi8(_mm_load_si128(&y2[t]), ut);

#define __dp_code_block3 \
	xt1 = _mm_loadu_si128((const __m128i*)(xo8  + t * 16 + 15)); /* xt1 <- x[r-1][t+15..t+30]; one block past en_ lands in the pad */ \
	vt1 = _mm_loadu_si128((const __m128i*)(vo8  + t * 16 + 15)); \
	x2t1= _mm_loadu_si128((const __m128i*)(x2o8 + t * 16 + 15));

// Rails store the same unbiased values as the original kernel below -- no
// bias is added or subtracted anywhere here, so the store macros are
// identical to the x86 ones. The overflow constraint is exactly the
// original's: x[t-1]+v[t-1] reaches -2(q+e), so it must fit int8, i.e.
// 2(q+e) <= 128.
#define __dp_code_block2 \
	_mm_store_si128(&u[t], _mm_sub_epi8(z, vt1));    /* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
	_mm_store_si128(&vn[t], _mm_sub_epi8(z, ut));    /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
	tmp = _mm_sub_epi8(z, q_); \
	a = _mm_sub_epi8(a, tmp); \
	b = _mm_sub_epi8(b, tmp); \
	tmp = _mm_sub_epi8(z, q2_); \
	a2= _mm_sub_epi8(a2, tmp); \
	b2= _mm_sub_epi8(b2, tmp);

#define __dp_store_x(V)  _mm_store_si128(&xn[t],  _mm_sub_epi8((V), qe_))
#define __dp_store_y(V)  _mm_store_si128(&y[t],   _mm_sub_epi8((V), qe_))
#define __dp_store_x2(V) _mm_store_si128(&x2n[t], _mm_sub_epi8((V), qe2_))
#define __dp_store_y2(V) _mm_store_si128(&y2[t],  _mm_sub_epi8((V), qe2_))

#define __dp_u8(I) u8[I]
#define __dp_v8(I) v8[I]
#define __dp_v8_i32x16(P, W) ksw_widen_i8x16_pair((P), (W))
#define __dp_v8_i32x4(P) ksw_i8x4_to_i32x4(P)
#else
// x1_, v1_ and x21_ carry the previous rail vector, not its last byte, so the
// shift-shift-or collapses to one alignr per rail
#define __dp_code_block1 \
	z = _mm_load_si128(&s[t]); \
	tmp = _mm_load_si128(&x[t]);                     /* tmp <- x[r-1][t..t+15] */ \
	xt1 = ksw_alignr15(tmp, x1_);                   /* xt1 <- x[r-1][t-1..t+14] */ \
	x1_ = tmp; \
	tmp = _mm_load_si128(&v[t]);                     /* tmp <- v[r-1][t..t+15] */ \
	vt1 = ksw_alignr15(tmp, v1_);                   /* vt1 <- v[r-1][t-1..t+14] */ \
	v1_ = tmp; \
	a = _mm_add_epi8(xt1, vt1);                      /* a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14] */ \
	ut = _mm_load_si128(&u[t]);                      /* ut <- u[t..t+15] */ \
	b = _mm_add_epi8(_mm_load_si128(&y[t]), ut);     /* b <- y[r-1][t..t+15] + u[r-1][t..t+15] */ \
	tmp = _mm_load_si128(&x2[t]); \
	x2t1= ksw_alignr15(tmp, x21_); \
	x21_= tmp; \
	a2= _mm_add_epi8(x2t1, vt1); \
	b2= _mm_add_epi8(_mm_load_si128(&y2[t]), ut);

#define __dp_code_block3 ((void)0)                   /* block1 shifts the rails itself */

#define __dp_code_block2 \
	_mm_store_si128(&u[t], _mm_sub_epi8(z, vt1));    /* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
	_mm_store_si128(&v[t], _mm_sub_epi8(z, ut));     /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
	tmp = _mm_sub_epi8(z, q_); \
	a = _mm_sub_epi8(a, tmp); \
	b = _mm_sub_epi8(b, tmp); \
	tmp = _mm_sub_epi8(z, q2_); \
	a2= _mm_sub_epi8(a2, tmp); \
	b2= _mm_sub_epi8(b2, tmp);

#define __dp_store_x(V)  _mm_store_si128(&x[t],  _mm_sub_epi8((V), qe_))
#define __dp_store_y(V)  _mm_store_si128(&y[t],  _mm_sub_epi8((V), qe_))
#define __dp_store_x2(V) _mm_store_si128(&x2[t], _mm_sub_epi8((V), qe2_))
#define __dp_store_y2(V) _mm_store_si128(&y2[t], _mm_sub_epi8((V), qe2_))

#define __dp_u8(I) u8[I]
#define __dp_v8(I) v8[I]
#define __dp_v8_i32x16(P, W) ksw_widen_i8x16_pair((P), (W))
#define __dp_v8_i32x4(P) ksw_i8x4_to_i32x4(P)
#endif

	int r, t, qe0 = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, last_max_H = 0, wl, wr, max_sc, min_sc, long_thres, long_diff;
	int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
	int32_t *H = 0, H0 = 0, last_H0_t = 0;
	uint8_t *qr, *sf, *mem, *mem2 = 0;
	__m128i q_, q2_, zero_, sc_mch_, pmat_;
	int use_lut;
	__m128i *u, *y, *y2, *s, *p = 0;
	__m128i qe_, qe2_;
#if KSW_EXTD2_NEON_KERNEL
	__m128i m1_;
	__m128i *vb[2], *xb[2], *x2b[2];
#else
	__m128i *v, *x, *x2;
#endif

	ksw_reset_extz(ez);
	if (m <= 1 || qlen <= 0 || tlen <= 0) return;

	if (q2 + e2 < q + e) t = q, q = q2, q2 = t, t = e, e = e2, e2 = t; // make sure q+e no larger than q2+e2

	zero_   = _mm_set1_epi8(0);
	q_      = _mm_set1_epi8(q);
	q2_     = _mm_set1_epi8(q2);
	qe_     = _mm_set1_epi8(q + e);
	qe2_    = _mm_set1_epi8(q2 + e2);
#if KSW_EXTD2_NEON_KERNEL
	m1_     = _mm_set1_epi8(-1); // a > -1 is a >= 0, the right-alignment flag test
#endif
	sc_mch_ = _mm_set1_epi8(mat[0]);
		// XOR-indexed substitution LUT; KSW_EZ_GENERIC_SC keeps the scalar mat[] lookup
	use_lut = !(flag & KSW_EZ_GENERIC_SC);
	if (use_lut) {
		int8_t pmat[16];
		int8_t w_N = mat[m*m-1] == 0? (int8_t)-e2 : mat[m*m-1];
		int lt;
		pmat[0] = mat[0];                                  /* match            */
		pmat[1] = pmat[2] = pmat[3] = mat[1];              /* ACGT mismatch    */
		for (lt = 4; lt < 16; ++lt) pmat[lt] = w_N;        /* anything with N  */
		pmat_ = _mm_loadu_si128((const __m128i*)pmat);
	} else pmat_ = zero_;

	if (w < 0) w = tlen > qlen? tlen : qlen;
	wl = wr = w;
	tlen_ = (tlen + 15) / 16;
	n_col_ = qlen < tlen? qlen : tlen;
	n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 15) / 16 + 1;
	qlen_ = (qlen + 15) / 16;
	for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
		max_sc = max_sc > mat[t]? max_sc : mat[t];
		min_sc = min_sc < mat[t]? min_sc : mat[t];
	}
	if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

	long_thres = e != e2? (q2 - q) / (e - e2) - 1 : 0;
	if (q2 + e2 + long_thres * e2 > q + e + long_thres * e)
		++long_thres;
	long_diff = long_thres * (e - e2) - (q2 - q) - e2;

#if KSW_EXTD2_NEON_KERNEL
	// nine padded rails (u, y, y2 and two buffers each of v, x, x2), then s, sf
	// and qr; the 16-byte pad in front of each rail keeps the t-1 loads in bounds
	mem = (uint8_t*)kcalloc(km, (size_t)(tlen_ + 1) * 9 + tlen_ * 2 + qlen_ + 2, 16);
	u = (__m128i*)(((size_t)mem + 15) >> 4 << 4) + 1; // 16-byte aligned, after the pad
	vb[0]  = u + tlen_ + 1,     vb[1]  = vb[0] + tlen_ + 1;
	xb[0]  = vb[1] + tlen_ + 1, xb[1]  = xb[0] + tlen_ + 1;
	y      = xb[1] + tlen_ + 1;
	x2b[0] = y + tlen_ + 1,     x2b[1] = x2b[0] + tlen_ + 1;
	y2     = x2b[1] + tlen_ + 1;
	s = y2 + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 16;
	memset(u,  -q  - e,  tlen_ * 16);
	memset(vb[0],  -q  - e,  tlen_ * 16), memset(vb[1],  -q  - e,  tlen_ * 16);
	memset(xb[0],  -q  - e,  tlen_ * 16), memset(xb[1],  -q  - e,  tlen_ * 16);
	memset(y,  -q  - e,  tlen_ * 16);
	memset(x2b[0], -q2 - e2, tlen_ * 16), memset(x2b[1], -q2 - e2, tlen_ * 16);
	memset(y2, -q2 - e2, tlen_ * 16);
#else
	mem = (uint8_t*)kcalloc(km, tlen_ * 8 + qlen_ + 1, 16);
	u = (__m128i*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
	v = u + tlen_, x = v + tlen_, y = x + tlen_, x2 = y + tlen_, y2 = x2 + tlen_;
	s = y2 + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 16;
	memset(u,  -q  - e,  tlen_ * 16);
	memset(v,  -q  - e,  tlen_ * 16);
	memset(x,  -q  - e,  tlen_ * 16);
	memset(y,  -q  - e,  tlen_ * 16);
	memset(x2, -q2 - e2, tlen_ * 16);
	memset(y2, -q2 - e2, tlen_ * 16);
#endif
	if (!approx_max) {
		H = (int32_t*)kmalloc(km, tlen_ * 16 * 4);
		for (t = 0; t < tlen_ * 16; ++t) H[t] = KSW_NEG_INF;
	}
	if (with_cigar) {
		mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 16);
		p = (__m128i*)(((size_t)mem2 + 15) >> 4 << 4);
		off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
		off_end = off + qlen + tlen - 1;
	}

	if (use_lut) {
		/* query-N -> 8, which keeps every sf ^ qrr index <= 12. Confined to the
		 * prepass: qr[] has no other reader. */
		for (t = 0; t < qlen; ++t) {
			uint8_t c = query[qlen - 1 - t];
			qr[t] = c == m - 1? 8 : c;
		}
	} else {
		for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
	}
	memcpy(sf, target, tlen);

	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1, st0, en0, st_, en_;
		int8_t x1, x21, v1;
		uint8_t *qrr = qr + (qlen - 1 - r);
		__m128i xt1, x2t1, vt1;
#if KSW_EXTD2_NEON_KERNEL
		__m128i *xn = xb[r&1], *vn = vb[r&1], *x2n = x2b[r&1]; // row r
		int8_t *u8 = (int8_t*)u, *v8 = (int8_t*)vn;
		int8_t *xo8 = (int8_t*)xb[(r+1)&1], *vo8 = (int8_t*)vb[(r+1)&1], *x2o8 = (int8_t*)x2b[(r+1)&1]; // row r-1
#else
		int8_t *u8 = (int8_t*)u, *v8 = (int8_t*)v, *x8 = (int8_t*)x, *x28 = (int8_t*)x2;
		__m128i x1_, x21_, v1_;
#endif
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-wr+1)>>1) st = (r-wr+1)>>1; // take the ceil
		if (en > (r+wl)>>1) en = (r+wl)>>1; // take the floor
		if (st > en) {
			ez->zdropped = 1;
			break;
		}
		st0 = st, en0 = en;
		st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
#if KSW_EXTD2_NEON_KERNEL
		// set boundary conditions: same unbiased values as the original kernel
		// below; only the (r-1,s-1) carry-in source differs (the previous
		// row's double-buffer slot rather than the single shared rail).
		if (st > 0) {
			if (st - 1 >= last_st && st - 1 <= last_en) {
				x1 = xo8[st - 1], x21 = x2o8[st - 1], v1 = vo8[st - 1]; // (r-1,s-1) calculated in the last round
			} else {
				x1 = -q - e, x21 = -q2 - e2;
				v1 = -q - e;
			}
		} else {
			x1 = -q - e, x21 = -q2 - e2;
			v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
		}
		if (en >= r) {
			((int8_t*)y)[r] = -q - e, ((int8_t*)y2)[r] = -q2 - e2;
			u8[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
		}
#else
		// set boundary conditions
		if (st > 0) {
			if (st - 1 >= last_st && st - 1 <= last_en) {
				x1 = x8[st - 1], x21 = x28[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
			} else {
				x1 = -q - e, x21 = -q2 - e2;
				v1 = -q - e;
			}
		} else {
			x1 = -q - e, x21 = -q2 - e2;
			v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
		}
		if (en >= r) {
			((int8_t*)y)[r] = -q - e, ((int8_t*)y2)[r] = -q2 - e2;
			u8[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
		}
#endif
		// loop fission: set scores first
		if (use_lut) {
			// 5 ops -> 2: one XOR, one byte shuffle
			for (t = st0; t <= en0; t += 16) {
				__m128i sq, st;
				sq = _mm_loadu_si128((__m128i*)&sf[t]);
				st = _mm_loadu_si128((__m128i*)&qrr[t]);
				_mm_storeu_si128((__m128i*)((int8_t*)s + t),
								 ksw_shuffle_epi8(pmat_, ksw_xor_si128(sq, st)));
			}
		} else {
			for (t = st0; t <= en0; ++t)
				((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
		}
		// core loop
		st_ = st / 16, en_ = en / 16;
#if KSW_EXTD2_NEON_KERNEL
		xt1  = ksw_insert0(_mm_loadu_si128((const __m128i*)(xo8  + st - 1)), x1);  /* x[r-1][st-1..st+14] */
		vt1  = ksw_insert0(_mm_loadu_si128((const __m128i*)(vo8  + st - 1)), v1);
		x2t1 = ksw_insert0(_mm_loadu_si128((const __m128i*)(x2o8 + st - 1)), x21);
#else
		// lane 15: ksw_alignr15() reads the carry from the top byte of the previous vector
		x1_  = _mm_slli_si128(_mm_cvtsi32_si128((uint8_t)x1),  15);
		x21_ = _mm_slli_si128(_mm_cvtsi32_si128((uint8_t)x21), 15);
		v1_  = _mm_slli_si128(_mm_cvtsi32_si128((uint8_t)v1),  15);
#endif
		assert(en_ - st_ + 1 <= n_col_);
		if (!with_cigar) { // score only
			for (t = st_; t <= en_; ++t) {
				__m128i z, a, b, a2, b2, ut, tmp;
				__dp_code_block1;
#ifdef __SSE4_1__
				z = _mm_max_epi8(z, a);
				z = _mm_max_epi8(z, b);
				z = _mm_max_epi8(z, a2);
				z = _mm_max_epi8(z, b2);
				z = _mm_min_epi8(z, sc_mch_);
				__dp_code_block2; // save u[] and v[]; update a, b, a2 and b2
				__dp_store_x(_mm_max_epi8(a,  zero_));
				__dp_store_y(_mm_max_epi8(b,  zero_));
				__dp_store_x2(_mm_max_epi8(a2, zero_));
				__dp_store_y2(_mm_max_epi8(b2, zero_));
				__dp_code_block3;
#else
				tmp = _mm_cmpgt_epi8(a,  z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(b,  z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(a2, z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a2));
				tmp = _mm_cmpgt_epi8(b2, z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b2));
				tmp = _mm_cmplt_epi8(sc_mch_, z);
				z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
				__dp_code_block2;
				tmp = _mm_cmpgt_epi8(a, zero_);
				__dp_store_x(_mm_and_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(b, zero_);
				__dp_store_y(_mm_and_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(a2, zero_);
				__dp_store_x2(_mm_and_si128(tmp, a2));
				tmp = _mm_cmpgt_epi8(b2, zero_);
				__dp_store_y2(_mm_and_si128(tmp, b2));
				__dp_code_block3;
#endif
			}
		} else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
			__m128i *pr = p + (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (t = st_; t <= en_; ++t) {
				__m128i d, z, a, b, a2, b2, ut, tmp;
				__dp_code_block1;
#ifdef __SSE4_1__
				d = ksw_mask01(_mm_cmpgt_epi8(a, z));                            // d = a  > z? 1 : 0
				z = _mm_max_epi8(z, a);
				d = ksw_bitsel(d, _mm_set1_epi8(2), _mm_cmpgt_epi8(b,  z)); // d = b  > z? 2 : d
				z = _mm_max_epi8(z, b);
				d = ksw_bitsel(d, _mm_set1_epi8(3), _mm_cmpgt_epi8(a2, z)); // d = a2 > z? 3 : d
				z = _mm_max_epi8(z, a2);
				d = ksw_bitsel(d, _mm_set1_epi8(4), _mm_cmpgt_epi8(b2, z)); // d = b2 > z? 4 : d
				z = _mm_max_epi8(z, b2);
				z = _mm_min_epi8(z, sc_mch_);
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
				tmp = _mm_cmpgt_epi8(a,  z);
				d = _mm_and_si128(tmp, _mm_set1_epi8(1));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(b,  z);
				d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, _mm_set1_epi8(2)));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(a2, z);
				d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, _mm_set1_epi8(3)));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a2));
				tmp = _mm_cmpgt_epi8(b2, z);
				d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, _mm_set1_epi8(4)));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b2));
				tmp = _mm_cmplt_epi8(sc_mch_, z);
				z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
#endif
				__dp_code_block2;
				tmp = _mm_cmpgt_epi8(a, zero_);
				__dp_store_x(_mm_and_si128(tmp, a));
				d = ksw_bitins_flag(d, tmp, _mm_set1_epi8(0x08)); // d = a > 0? 1<<3 : 0
				tmp = _mm_cmpgt_epi8(b, zero_);
				__dp_store_y(_mm_and_si128(tmp, b));
				d = ksw_bitins_flag(d, tmp, _mm_set1_epi8(0x10)); // d = b > 0? 1<<4 : 0
				tmp = _mm_cmpgt_epi8(a2, zero_);
				__dp_store_x2(_mm_and_si128(tmp, a2));
				d = ksw_bitins_flag(d, tmp, _mm_set1_epi8(0x20)); // d = a2 > 0? 1<<5 : 0
				tmp = _mm_cmpgt_epi8(b2, zero_);
				__dp_store_y2(_mm_and_si128(tmp, b2));
				d = ksw_bitins_flag(d, tmp, _mm_set1_epi8(0x40)); // d = b2 > 0? 1<<6 : 0
				_mm_store_si128(&pr[t], d);
				__dp_code_block3;
			}
		} else { // gap right-alignment
			__m128i *pr = p + (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (t = st_; t <= en_; ++t) {
				__m128i d, z, a, b, a2, b2, ut, tmp;
				__dp_code_block1;
#ifdef __SSE4_1__
				d = _mm_andnot_si128(_mm_cmpgt_epi8(z, a), _mm_set1_epi8(1));    // d = z > a?  0 : 1
				z = _mm_max_epi8(z, a);
				d = _mm_blendv_epi8(_mm_set1_epi8(2), d, _mm_cmpgt_epi8(z, b));  // d = z > b?  d : 2
				z = _mm_max_epi8(z, b);
				d = _mm_blendv_epi8(_mm_set1_epi8(3), d, _mm_cmpgt_epi8(z, a2)); // d = z > a2? d : 3
				z = _mm_max_epi8(z, a2);
				d = _mm_blendv_epi8(_mm_set1_epi8(4), d, _mm_cmpgt_epi8(z, b2)); // d = z > b2? d : 4
				z = _mm_max_epi8(z, b2);
				z = _mm_min_epi8(z, sc_mch_);
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
				tmp = _mm_cmpgt_epi8(z, a);
				d = _mm_andnot_si128(tmp, _mm_set1_epi8(1));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(z, b);
				d = _mm_or_si128(_mm_and_si128(tmp, d), _mm_andnot_si128(tmp, _mm_set1_epi8(2)));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(z, a2);
				d = _mm_or_si128(_mm_and_si128(tmp, d), _mm_andnot_si128(tmp, _mm_set1_epi8(3)));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, a2));
				tmp = _mm_cmpgt_epi8(z, b2);
				d = _mm_or_si128(_mm_and_si128(tmp, d), _mm_andnot_si128(tmp, _mm_set1_epi8(4)));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, b2));
				tmp = _mm_cmplt_epi8(sc_mch_, z);
				z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
#endif
				__dp_code_block2;
#if KSW_EXTD2_NEON_KERNEL
				// the flag masks are a >= 0 (not a > 0 as on the left), computed directly
				// so the store is one and and the flag one insert, as on the left
				tmp = _mm_cmpgt_epi8(a, m1_);
				__dp_store_x(_mm_and_si128(tmp, a));
				d = ksw_bitins_flag(d, tmp, _mm_set1_epi8(0x08)); // d = a >= 0? 1<<3 : 0
				tmp = _mm_cmpgt_epi8(b, m1_);
				__dp_store_y(_mm_and_si128(tmp, b));
				d = ksw_bitins_flag(d, tmp, _mm_set1_epi8(0x10)); // d = b >= 0? 1<<4 : 0
				tmp = _mm_cmpgt_epi8(a2, m1_);
				__dp_store_x2(_mm_and_si128(tmp, a2));
				d = ksw_bitins_flag(d, tmp, _mm_set1_epi8(0x20)); // d = a2 >= 0? 1<<5 : 0
				tmp = _mm_cmpgt_epi8(b2, m1_);
				__dp_store_y2(_mm_and_si128(tmp, b2));
				d = ksw_bitins_flag(d, tmp, _mm_set1_epi8(0x40)); // d = b2 >= 0? 1<<6 : 0
#else
				tmp = _mm_cmpgt_epi8(zero_, a);
				__dp_store_x(_mm_andnot_si128(tmp, a));
				d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x08))); // d = a > 0? 1<<3 : 0
				tmp = _mm_cmpgt_epi8(zero_, b);
				__dp_store_y(_mm_andnot_si128(tmp, b));
				d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x10))); // d = b > 0? 1<<4 : 0
				tmp = _mm_cmpgt_epi8(zero_, a2);
				__dp_store_x2(_mm_andnot_si128(tmp, a2));
				d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x20))); // d = a > 0? 1<<5 : 0
				tmp = _mm_cmpgt_epi8(zero_, b2);
				__dp_store_y2(_mm_andnot_si128(tmp, b2));
				d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x40))); // d = b > 0? 1<<6 : 0
#endif
				_mm_store_si128(&pr[t], d);
				__dp_code_block3;
			}
		}
		if (!approx_max) { // find the exact max with a 32-bit score array
			int32_t max_H, max_t;
			// compute H[], max_H and max_t
			if (r > 0) {
				int32_t HH[4], tt[4], en1 = st0 + (en0 - st0) / 4 * 4, i;
				__m128i max_H_, max_t_, t_, t4_ = _mm_set1_epi32(4);
				max_H = H[en0] = en0 > 0? H[en0-1] + __dp_u8(en0) : H[en0] + __dp_v8(en0); // special casing the last element
				max_t = en0;
				max_H_ = _mm_set1_epi32(max_H);
				max_t_ = _mm_set1_epi32(max_t);
				t_ = _mm_set1_epi32(st0);
				// the four-lane accumulator and the update order are what keep max_t
				// reproducible; only the widening changes here
#define __ksw_hmax_step(HT, W) do { \
					__m128i H1 = _mm_loadu_si128((__m128i*)(HT)), msk; \
					H1 = _mm_add_epi32(H1, (W)); \
					_mm_storeu_si128((__m128i*)(HT), H1); \
					msk = _mm_cmpgt_epi32(H1, max_H_); \
					_ksw_hmax_sel(H1, msk); \
					t_ = _mm_add_epi32(t_, t4_); \
				} while (0)
#ifdef __SSE4_1__
#define _ksw_hmax_sel(H1, msk) do { \
					max_H_ = _mm_max_epi32(max_H_, (H1)); \
					max_t_ = _mm_blendv_epi8(max_t_, t_, (msk)); \
				} while (0)
#else
#define _ksw_hmax_sel(H1, msk) do { \
					max_H_ = _mm_or_si128(_mm_and_si128((msk), (H1)), _mm_andnot_si128((msk), max_H_)); \
					max_t_ = _mm_or_si128(_mm_and_si128((msk), t_), _mm_andnot_si128((msk), max_t_)); \
				} while (0)
#endif
				{
					int en16 = st0 + (en0 - st0) / 16 * 16;
					for (t = st0; t < en16; t += 16) {
						__m128i w[4];
						__dp_v8_i32x16(&v8[t], w);
						__ksw_hmax_step(&H[t],      w[0]);
						__ksw_hmax_step(&H[t +  4], w[1]);
						__ksw_hmax_step(&H[t +  8], w[2]);
						__ksw_hmax_step(&H[t + 12], w[3]);
					}
				}
				for (; t < en1; t += 4) { // this implements: H[t]+=v8[t]; if(H[t]>max_H) max_H=H[t],max_t=t;
					__ksw_hmax_step(&H[t], __dp_v8_i32x4(&v8[t]));
				}
#undef __ksw_hmax_step
#undef _ksw_hmax_sel
				_mm_storeu_si128((__m128i*)HH, max_H_);
				_mm_storeu_si128((__m128i*)tt, max_t_);
				for (i = 0; i < 4; ++i)
					if (max_H < HH[i]) max_H = HH[i], max_t = tt[i] + i;
				for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE
					H[t] += (int32_t)__dp_v8(t);
					if (H[t] > max_H)
						max_H = H[t], max_t = t;
				}
			} else H[0] = __dp_v8(0) - qe0, max_H = H[0], max_t = 0; // special casing r==0
			// update ez
			if (en0 == tlen - 1 && H[en0] > ez->mte)
				ez->mte = H[en0], ez->mte_q = r - en0;
			if (r - st0 == qlen - 1 && H[st0] > ez->mqe)
				ez->mqe = H[st0], ez->mqe_t = st0;
			if (ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e2)) break;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H[tlen - 1];
			// early stopping in the extension-only mode; this block should not change the alignment
			if (flag & KSW_EZ_EXTZ_ONLY) {
				int32_t rH = last_max_H > max_H? last_max_H : max_H; // NB: if diagonal r is on the optimal path, r-1 or r+1 is not on the optimal path
				int32_t rq = qlen - (r - st0), rt = tlen - en0;
				// 3 conditions: hitting bottom of the DP matrix, hitting the right boundary of the matrix, and bracketing the diagonal
				int32_t rm = rq >= tlen - st0? tlen - st0 : rt >= qlen - (r - en0)? qlen - (r - en0) : tlen + qlen - 1 - r;
				if (rH + rm * max_sc + end_bonus < ez->max) break;
			}
			last_max_H = max_H; // need this for calculating rH
		} else { // find approximate max; Z-drop might be inaccurate, too.
			if (r > 0) {
				if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
					int32_t d0 = __dp_v8(last_H0_t);
					int32_t d1 = __dp_u8(last_H0_t + 1);
					if (d0 > d1) H0 += d0;
					else H0 += d1, ++last_H0_t;
				} else if (last_H0_t >= st0 && last_H0_t <= en0) {
					H0 += __dp_v8(last_H0_t);
				} else {
					++last_H0_t, H0 += __dp_u8(last_H0_t);
				}
			} else H0 = __dp_v8(0) - qe0, last_H0_t = 0;
			if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e2)) break;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H0;
		}
		last_st = st, last_en = en;
		//for (t = st0; t <= en0; ++t) printf("(%d,%d)\t(%d,%d,%d,%d)\t%d\n", r, t, ((int8_t*)u)[t], ((int8_t*)v)[t], ((int8_t*)x)[t], ((int8_t*)y)[t], H[t]); // for debugging
	}
	kfree(km, mem);
	if (!approx_max) kfree(km, H);
	if (with_cigar) { // backtrack
		int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
		if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if ((flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
			ez->reach_end = 1;
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (ez->max_t >= 0 && ez->max_q >= 0) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		}
		kfree(km, mem2); kfree(km, off);
	}
}
