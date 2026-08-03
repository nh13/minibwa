#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "mbpriv.h"
#include "kalloc.h"
#include "ksort.h"

#define key_sai0(a) ((a).x[0])
KRADIX_SORT_INIT(mb_sai0, mb_sai_t, key_sai0, 8)

#define key_sais(a) ((a).size)
KRADIX_SORT_INIT(mb_sais, mb_sai_t, key_sais, 8)

#define key_saii(a) ((a).info)
KRADIX_SORT_INIT(mb_saii, mb_sai_t, key_saii, 8)

#define key_anchor(a) ((a).tpos)
KRADIX_SORT_INIT(mb_anchor, mb_anchor_t, key_anchor, 8)

/***********
 * Seeding *
 ***********/

void mb_seed_intv(void *km, const mb_bwt_t *bwt, int32_t len, const uint8_t *seq, int32_t min_len, int32_t max_sub_occ, int32_t min_sub_occ, mb_sai_v *v)
{
	int64_t x = 0, i, n_a0;
	mb_sai_t p;

	v->n = 0;
	do { // pass 1: standard SMEMs
		x = mb_bwt_smem(bwt, len, seq, x, min_len, 1, &p);
		if (p.size > 0) {
			Kgrow(km, mb_sai_t, v->a, v->n, v->m);
			v->a[v->n++] = p;
		}
	} while (x < len);

	n_a0 = v->n;
	for (i = 0; i < n_a0; ++i) { // pass 2: sub-SMEMs
		int32_t sub_min_len;
		uint32_t st = v->a[i].info>>32, en = (uint32_t)v->a[i].info;
		if (en - st < min_len * 2 || v->a[i].size > max_sub_occ || v->a[i].size < min_sub_occ)
			continue;
		x = st;
		sub_min_len = (en - st) / 2 > min_len? (en - st) / 2 : min_len;
		do { // if two SMEMs have large overlaps, we may find the same sub intervals in both. A rare case not worth optimizing
			x = mb_bwt_smem(bwt, en, seq, x, sub_min_len, v->a[i].size + 1, &p);
			if (p.size > v->a[i].size) {
				Kgrow(km, mb_sai_t, v->a, v->n, v->m);
				v->a[v->n++] = p;
			}
		} while (x < en);
	}
}

void mb_seed_intv_batch(void *km, const mb_bwt_t *bwt, int32_t n_seq, const int32_t *len, uint8_t *const* seq, int32_t min_len, int32_t max_sub_occ, int32_t min_sub_occ, mb_sai_v *v)
{ // identical to mb_seed_intv() though the order of intervals is often different
	const int max_batch_size = 50;
	mb_smem_entry_t *s;
	int32_t i, j, n_s, *nv;

	// first pass: standard SMEMs
	s = Kcalloc(km, mb_smem_entry_t, max_batch_size);
	nv = Kcalloc(km, int32_t, n_seq);
	for (i = 0; i < n_seq; ++i) v[i].n = 0;
	for (i = 0; i < n_seq; i += max_batch_size) {
		int32_t en = i + max_batch_size < n_seq? i + max_batch_size : n_seq;
		for (j = i; j < en; ++j) {
			mb_smem_entry_t *t = &s[j - i];
			t->min_len = min_len;
			t->min_occ = 1;
			t->st = 0, t->en = len[j];
			t->q = seq[j];
			t->v = &v[j];
		}
		mb_bwt_smem_batch(km, bwt, en - i, s);
	}

	// second pass; sub-SMEMs
	for (i = 0; i < n_seq; ++i) nv[i] = v[i].n;
	for (i = n_s = 0; i < n_seq; ++i) {
		for (j = 0; j < nv[i]; ++j) {
			mb_smem_entry_t *t;
			uint32_t st = v[i].a[j].info>>32, en = (uint32_t)v[i].a[j].info;
			if (en - st < min_len * 2 || v[i].a[j].size > max_sub_occ || v[i].a[j].size < min_sub_occ)
				continue;
			t = &s[n_s++];
			t->min_len = (en - st) / 2 > min_len? (en - st) / 2 : min_len;
			t->min_occ = v[i].a[j].size + 1;
			t->st = st, t->en = en;
			t->q = seq[i];
			t->v = &v[i];
			if (n_s == max_batch_size) {
				mb_bwt_smem_batch(km, bwt, n_s, s);
				n_s = 0;
			}
		}
	}
	if (n_s > 0)
		mb_bwt_smem_batch(km, bwt, n_s, s);
	kfree(km, nv);
	kfree(km, s);
}

/*****************************
 * Seed/anchor deduplication *
 *****************************/

static void mb_seed_sort_dedup(mb_sai_v *u)
{
	int64_t i, i0, j;
	if (u->n <= 1) return;
	// sort by ::{x[0],size,info}
	radix_sort_mb_sai0(u->a, u->a + u->n); // sort by ::x[0]
	for (i = 1, i0 = 0; i <= u->n; ++i) {
		if (i == u->n || u->a[i].x[0] != u->a[i0].x[0]) {
			if (i - i0 > 1) {
				int64_t k, k0, n1 = i - i0;
				mb_sai_t *a1 = &u->a[i0];
				radix_sort_mb_sais(&u->a[i0], &u->a[i]); // sort by ::size
				kom_reverse(mb_sai_t, n1, a1);
				for (k = i0 + 1, k0 = i0; k <= i; ++k) {
					if (k == i || u->a[k0].size != u->a[k].size) {
						if (k - k0 > 1)
							radix_sort_mb_saii(&u->a[k0], &u->a[k]); // sort by ::info
						k0 = k;
					}
				}
			}
			i0 = i;
		}
	}
	// dedup
	for (i = 1, j = 0; i < u->n; ++i)
		if (!(u->a[i].x[0] == u->a[j].x[0] && u->a[i].size == u->a[j].size && u->a[i].info == u->a[j].info))
			u->a[++j] = u->a[i];
	u->n = j + 1;
}

/* Remove duplicated anchors. The two-round seeding algorithm may lead an
 * anchor precisely contained in a longer anchor. This routine filters out the
 * shorter anchor. This wouldn't happen to minimap2. */
static void mb_anchor_dedup(mb_anchor_v *v) // NB: assuming sorted by tpos
{
	const int max_back = 100; // to avoid quadratic behavior in the worst case
	int64_t i, j, k;
	for (i = 1; i < v->n; ++i) {
		mb_anchor_t *ai = &v->a[i];
		int64_t tsj, tsi = ai->tpos + 1 - ai->len;
		int32_t qsj, qsi = ai->qpos + 1 - ai->len;
		for (j = i - 1, k = 0; j >= 0 && k < max_back; --j, ++k) {
			mb_anchor_t *aj = &v->a[j];
			if (aj->sid != ai->sid) break;
			if (aj->tpos < tsi) break;
			tsj = aj->tpos + 1 - aj->len;
			qsj = aj->qpos + 1 - aj->len;
			if (tsj >= tsi) { // then j is contained in i
				if (tsj - tsi == qsj - qsi) aj->flt = 1;
			} else if (ai->tpos == aj->tpos) { // then i is contained in j
				if (ai->qpos == aj->qpos) ai->flt = 1;
			}
		}
	}
	for (i = j = 0; i < v->n; ++i)
		if (!v->a[i].flt) v->a[j++] = v->a[i];
	v->n = j;
}

/****************************************
 * ALT-seed -> primary anchor projection *
 ****************************************/

/* Maximum number of projected primary anchors injected per (pri_tid, folded
 * strand) bucket.  Segduplicated loci can have hundreds of ALT/paralog copies;
 * without a cap, near-duplicate paralog projections would blow up mb_lchain_dp's
 * O(n*max_iter) inner loop on exactly the repeat-heavy reads.  A small cap is
 * sufficient: we only need ONE surviving primary anchor at the lifted locus to
 * seed a primary candidate (the rest are near-duplicates the following
 * mb_anchor_dedup collapses), so capping here never under-recovers. */
#define MB_PROJ_CAP_PER_LOCUS 4

/* Maximum number of DISTINCT projected primary loci tracked per read.  Sized so
 * it effectively never fires on real data: a read projecting to more distinct
 * loci than this is in a massive repeat family and will be MAPQ 0 regardless, so
 * stopping early is correctness-neutral.  The bound exists only to keep the
 * per-locus bookkeeping table on the stack.  Both caps are silent in production
 * but counted and emitted under --dbg-alt-proj so truncation is observable. */
#define MB_PROJ_MAX_LOCI 256

/* Project ALT-contig anchors onto the primary assembly so a segduplicated
 * primary locus gets a candidate even when max_occ subsampling drops the
 * primary's OWN seed (DRAGEN's mechanism: use ALT-contig seed matches to obtain
 * the corresponding primary alignment).
 *
 * Called from mb_anchor() AFTER process_batch() has filled v->a[] and BEFORE the
 * radix_sort + tpos-rebase + mb_anchor_dedup at the tail of mb_anchor().  All
 * coordinates here are in the CONCATENATED frame (mirroring seed.c:207); the
 * existing sort/rebase/dedup then handle ordering and exact-duplicate removal of
 * the injected anchors for free.
 *
 * Each injected anchor is written as a NATIVE-equivalent primary anchor (no
 * provenance bit any chaining code reads): sid = pri_tid<<1 | folded_rev, with
 * tpos/qpos = last base in the strand-FOLDED concatenated frame and len = seed
 * length.  It then chains normally under comput_sc() as a genuine primary anchor
 * -- this is NOT a chainer-merge of an ALT anchor into an ALT chain. */
static void mb_anchor_project_alt(void *km, const l2b_t *l2b, int32_t qlen, mb_anchor_v *v)
{
	int64_t i, n0 = v->n;
	/* --dbg-no-alt-proj ablates projection so the segdup regression tests can
	 * compare WITH vs WITHOUT projection from a single binary; production leaves
	 * it off (projection is unconditional, gated only on per-anchor is_alt). */
	if (kom_dbg_flag & MB_DBG_NO_ALT_PROJ) return;
	/* --dbg-alt-proj traces each projected primary anchor (see the trace below). */
	int proj_trace = (kom_dbg_flag & MB_DBG_ALT_PROJ) != 0;
	/* Per-locus cap bookkeeping: a tiny rolling table keyed by the projected
	 * sid (pri_tid<<1|folded_rev) AND the projected forward last base
	 * (fold_last), so the cap is per distinct projected LOCUS rather than per
	 * (contig,strand) -- two paralog seeds that lift to different positions on
	 * the same contig+strand must NOT share a cap slot.  Loci are few per read
	 * in practice, so a linear scan is fine; this also dedups
	 * projected-vs-projected at the same locus. */
	int32_t cap_sid[MB_PROJ_MAX_LOCI];
	int64_t cap_pos[MB_PROJ_MAX_LOCI];
	int32_t cap_cnt[MB_PROJ_MAX_LOCI];
	int32_t n_cap = 0;
	int64_t n_drop_locuscap = 0, n_drop_tablefull = 0; /* observability (trace only) */

	if (n0 == 0) return;

	for (i = 0; i < n0; ++i) {
		const mb_anchor_t *q = &v->a[i];
		int64_t alt_tid = q->sid >> 1;
		int32_t alt_rev = q->sid & 1;
		const l2b_ctg_t *alt_ctg;
		int64_t alt_cst, alt_clast;       /* ALT contig-local forward span [cst, clast] inclusive */
		int64_t pt_lo, pt_hi;             /* lifted primary tids of the two endpoints */
		uint64_t pp_lo, pp_hi;            /* lifted primary positions (forward, contig-local) */
		uint8_t rv_lo, rv_hi;
		int64_t pri_tid, pri_st, pri_en;  /* primary forward span [st, en] inclusive */
		uint8_t blk_rev, folded_rev;
		const l2b_ctg_t *pri_ctg;
		int64_t qf_s;                     /* query forward start of the seed */
		int64_t new_qpos, fold_last, new_tpos;
		int32_t new_sid, j, c;
		mb_anchor_t *p;

		/* Gate: only ALT contigs with lift blocks project. */
		if (alt_tid < 0 || alt_tid >= (int64_t)l2b->n_ctg) continue;
		alt_ctg = &l2b->ctg[alt_tid];
		if (!alt_ctg->is_alt || alt_ctg->n_lift == 0) continue;

		/* Recover the ALT contig-local FORWARD span from the concatenated tpos
		 * (inverse of process_batch's q->tpos = off*2 + len*rev + cst + len-1).
		 * The position component is the FOLDED last base for the seed strand;
		 * recover that fold first, then unfold per strand to a forward span.
		 * For a reverse seed the folded last base is the forward FIRST base of
		 * the span, so cst = len-1-fold_last; for forward it is the last base. */
		int64_t alt_fold_last = q->tpos - alt_ctg->off * 2 - alt_ctg->len * alt_rev;
		if (alt_fold_last < 0 || alt_fold_last >= (int64_t)alt_ctg->len) continue;
		alt_cst = alt_rev ? (int64_t)alt_ctg->len - 1 - alt_fold_last : alt_fold_last - (q->len - 1);
		alt_clast = alt_cst + q->len - 1; /* inclusive last forward ALT base */
		if (alt_cst < 0 || alt_clast >= (int64_t)alt_ctg->len) continue;

		/* Lift both inclusive endpoints; a reverse block maps low ALT -> high
		 * primary, so take min/max over the two lifted outputs.  Require both to
		 * lift, to the same primary tid and the same block strand (a hole or a
		 * cross-block seed yields no clean primary anchor -> skip). */
		if (!l2b_lift(l2b, alt_tid, (uint64_t)alt_cst,   &pt_lo, &pp_lo, &rv_lo)) continue;
		if (!l2b_lift(l2b, alt_tid, (uint64_t)alt_clast, &pt_hi, &pp_hi, &rv_hi)) continue;
		if (pt_lo != pt_hi || rv_lo != rv_hi) continue;
		pri_tid = pt_lo;
		blk_rev = rv_lo;
		pri_st = (int64_t)pp_lo < (int64_t)pp_hi ? (int64_t)pp_lo : (int64_t)pp_hi;
		pri_en = (int64_t)pp_lo > (int64_t)pp_hi ? (int64_t)pp_lo : (int64_t)pp_hi;
		if (pri_tid < 0 || pri_tid >= (int64_t)l2b->n_ctg) continue;
		pri_ctg = &l2b->ctg[pri_tid];

		/* Reject length-changing lifts: a seed spanning an indel or two adjacent
		 * lift blocks maps to a primary span whose length differs from the seed
		 * length.  Injecting it as a contiguous len-bp anchor would corrupt the
		 * chain coordinates, so skip it. */
		if (pri_en - pri_st + 1 != q->len) continue;

		/* Strand fold: the projected primary strand is the .alt block strand
		 * XOR the seed's strand on the ALT contig (mirrors mb_hit_place's
		 * blk_rev ^ h->rev). */
		folded_rev = (uint8_t)(blk_rev ^ alt_rev);

		/* Recover the query FORWARD start of the seed, then re-fold qpos for the
		 * projected strand (qpos = last base in the folded query frame). */
		qf_s = alt_rev ? (int64_t)qlen - 1 - q->qpos : q->qpos - (q->len - 1);
		if (qf_s < 0 || qf_s + q->len > qlen) continue;
		new_qpos = folded_rev ? (int64_t)qlen - 1 - qf_s : qf_s + q->len - 1;

		/* Forward contig-local LAST base of the primary span.  process_batch
		 * stores tpos's position component in the FORWARD contig frame for BOTH
		 * strands (strand lives in sid&1 plus the len*rev half-frame shift; the
		 * consumer mb_hit_set_coor does ts = tpos+1-len with no reverse-unfold).
		 * The min/max over the two lifted endpoints already handled the
		 * reverse-block low-alt -> high-primary inversion, so use pri_en
		 * unconditionally -- do NOT re-fold for folded_rev. */
		fold_last = pri_en;
		if (fold_last < 0 || fold_last >= (int64_t)pri_ctg->len) continue;

		/* Concatenated-frame tpos, mirroring seed.c:207. */
		new_tpos = pri_ctg->off * 2 + pri_ctg->len * folded_rev + fold_last;
		new_sid = (int32_t)(pri_tid << 1 | folded_rev);

		/* --dbg-alt-proj emits one line per projected anchor giving the primary
		 * contig, 1-based POS, and projected strand -- the load-bearing
		 * coordinates produced by the reverse-span recovery and forward-frame
		 * fold.  Observable even when the resulting alignment is masked at SAM
		 * level by identical-scoring paralog collapse (reverse RC repeats), so a
		 * fixture can assert the projected locus directly.  Diagnostics only. */
		if (proj_trace) {
			int64_t pri_pos1 = fold_last - q->len + 2; /* 1-based POS = (ts 0-based)+1 = (fold_last+1-len)+1 */
			fprintf(stderr, "MB_PROJ\t%s\t%lld\t%c\tlen=%d\n",
				pri_ctg->name, (long long)pri_pos1, folded_rev ? '-' : '+', q->len);
		}

		/* Per-locus cap + projected-vs-projected dedup at the same projected
		 * locus (sid + forward last base).  (Projected-vs-native exact
		 * duplicates are removed by the mb_anchor_dedup that runs right after
		 * this; here we only bound volume and squash redundant paralog
		 * projections to the same coordinate.) */
		c = -1;
		for (j = 0; j < n_cap; ++j)
			if (cap_sid[j] == new_sid && cap_pos[j] == fold_last) { c = j; break; }
		if (c < 0) {
			if (n_cap < (int32_t)(sizeof(cap_sid) / sizeof(cap_sid[0]))) {
				c = n_cap++;
				cap_sid[c] = new_sid;
				cap_pos[c] = fold_last;
				cap_cnt[c] = 0;
			} else { ++n_drop_tablefull; continue; } /* table full: stop projecting new loci */
		}
		if (cap_cnt[c] >= MB_PROJ_CAP_PER_LOCUS) { ++n_drop_locuscap; continue; }
		++cap_cnt[c];

		/* Inject the native-equivalent primary anchor (flag/flt = 0 via memset). */
		Kgrow(km, mb_anchor_t, v->a, v->n, v->m);
		p = &v->a[v->n++];
		memset(p, 0, sizeof(*p));
		p->sid  = new_sid;
		p->len  = q->len;
		p->qpos = (int32_t)new_qpos;
		p->tpos = new_tpos;
	}
	/* Make cap-driven truncation observable (default-off diagnostic seam): in
	 * production both caps are correctness-neutral, but a non-zero drop count on
	 * a repeat-heavy read is worth seeing when investigating recovery. */
	if (proj_trace && (n_drop_locuscap || n_drop_tablefull))
		fprintf(stderr, "MB_PROJ_CAP\tdropped_locuscap=%lld\tdropped_tablefull=%lld\tn_loci=%d\n",
			(long long)n_drop_locuscap, (long long)n_drop_tablefull, n_cap);
}

/************************
 * Get contig positions *
 ************************/

typedef struct { int64_t st, en; } anchor_aux_t;
typedef struct { int64_t a, i; } sa_aux_t;

static void process_batch(void *km, const mb_idx_t *idx, const anchor_aux_t *aux, int32_t m, const sa_aux_t *b, uint64_t *a, int32_t qlen, l2b_meth_t mt, const mb_sai_v *u, mb_anchor_v *v)
{
	int64_t j, k;
	for (k = 0; k < m; ++k) a[k] = b[k].a;
	mb_bwt_sa_batch(km, idx->bwt, m, a);
	for (k = 0; k < m; ++k) {
		const anchor_aux_t *p = &aux[b[k].i];
		for (j = p->st; j < p->en; ++j) {
			int32_t qs = u->a[j].info>>32, qe = (int32_t)u->a[j].info;
			int32_t rev, len = qe - qs;
			int64_t tid, cst;
			l2b_meth_t mt_anchor;
			const l2b_ctg_t *ctg;
			mb_anchor_t *q;
			if (mt != L2B_METH_NONE) {
				tid = l2b_intv2cid_meth(idx->l2b, a[k], a[k] + len, &mt_anchor, &cst, &rev);
				if (tid < 0) continue;
				// R1(C2T): keep c2t_f(copy0) and g2a_r(copy2); R2(G2A): keep g2a_f(copy1) and c2t_r(copy3)
				if ((mt_anchor == mt) != (rev == 0)) continue; // filter
			} else {
				tid = l2b_intv2cid(idx->l2b, a[k], a[k] + len, &cst, &rev);
				if (tid < 0) continue;
			}
			rev = !!rev;
			ctg = &idx->l2b->ctg[tid];
			Kgrow(km, mb_anchor_t, v->a, v->n, v->m);
			q = &v->a[v->n++];
			memset(q, 0, sizeof(*q));
			q->sid = tid << 1 | rev;
			q->len = len;
			q->qpos = rev? qlen - 1 - qs : qs + len - 1;
			q->tpos = ctg->off * 2 + ctg->len * rev + cst + len - 1; // for sorting; will be adjusted later
		}
	}
}

static void mb_anchor_split_meth(void *km, const l2b_t *l2b, int32_t min_len, int32_t qlen, const uint8_t *qseq0, l2b_meth_t mt0, mb_anchor_v *v)
{
	int64_t i, m_a = v->n * 2, n_a = 0;
	int32_t max_len = 0;
	mb_anchor_t *a;
	uint8_t *tseq, *qseq2[2];

	for (i = 0; i < v->n; ++i)
		if (max_len < v->a[i].len)
			max_len = v->a[i].len;
	tseq = Kmalloc(km, uint8_t, max_len + qlen * 2);
	qseq2[0] = tseq + max_len;
	qseq2[1] = qseq2[0] + qlen;
	memcpy(qseq2[0], qseq0, qlen);
	for (i = 0; i < qlen; ++i)
		qseq2[1][qlen - 1 - i] = qseq0[i] > 3? 4 : 3 - qseq0[i];

	a = Kmalloc(km, mb_anchor_t, m_a);
	for (i = 0; i < v->n; ++i) {
		mb_anchor_t *p, *q = &v->a[i];
		const l2b_ctg_t *ctg = &l2b->ctg[q->sid>>1];
		int32_t rev = q->sid&1;
		int64_t tpos = q->tpos - (ctg->off * 2 + ctg->len * rev); // NB: requiring concatenated ::tpos!!
		int64_t ts = tpos + 1 - q->len;
		int32_t qs = q->qpos + 1 - q->len, j, j0;
		uint8_t t_allow, q_allow;
		l2b_meth_t mt;
		const uint8_t *qseq = qseq2[rev] + qs;
		l2b_getseq(l2b, q->sid>>1, ts, ts + q->len, tseq);
		mt = q->sid&1? l2b_meth_rev(mt0) : mt0;
		t_allow = mt == L2B_METH_C2T? 1 : 2;
		q_allow = mt == L2B_METH_C2T? 3 : 0;
		for (j0 = j = 0; j <= q->len; ++j) {
			if (j == q->len || tseq[j] == 4 || qseq[j] == 4 || (tseq[j] != qseq[j] && !(tseq[j] == t_allow && qseq[j] == q_allow))) {
				if (j - j0 >= min_len) {
					Kgrow(km, mb_anchor_t, a, n_a, m_a);
					p = &a[n_a++];
					*p = *q;
					p->len = j - j0;
					p->qpos = q->qpos - (q->len - j);
					p->tpos = q->tpos - (q->len - j);
				}
				j0 = j + 1;
			}
		}
	}
	kfree(km, tseq);
	Kgrow(km, mb_anchor_t, v->a, n_a, v->m);
	memcpy(v->a, a, n_a * sizeof(mb_anchor_t));
	v->n = n_a;
	kfree(km, a);
}

/* Converting seed intervals to anchors. This function batches small SA
 * intervals and calls mb_bwt_sa_batch() in process_batch(). With prefetch, the
 * strategy noticeably improves the performance. */
double mb_anchor(void *km, const mb_idx_t *idx, mb_sai_v *u, int32_t min_len, int32_t qlen, const uint8_t *qseq, l2b_meth_t mt, int32_t max_occ, mb_anchor_v *v)
{
	const int batch_size = 20;
	int32_t n_aux, m, m_a;
	int64_t i, i0, j, k;
	uint64_t *a;
	double seed_ratio = 1.0;
	sa_aux_t *b;
	anchor_aux_t *aux;

	v->n = 0;
	if (u->n == 0) return seed_ratio; // no anchors
	mb_seed_sort_dedup(u);

	for (i = 0, k = 0; i < u->n; ++i) // pre-calculate the size of v->a
		k += u->a[i].size < max_occ? u->a[i].size : max_occ;
	Kgrow(km, mb_anchor_t, v->a, k - 1, v->m); // preallocate

	for (i = 1, i0 = 0, n_aux = 0; i <= u->n; ++i) // pre-compute n_aux
		if (i == u->n || u->a[i].x[0] != u->a[i0].x[0] || u->a[i].size != u->a[i0].size)
			++n_aux, i0 = i;
	aux = Kmalloc(km, anchor_aux_t, n_aux);
	for (i = 1, i0 = 0, n_aux = 0; i <= u->n; ++i) // populate aux[]
		if (i == u->n || u->a[i].x[0] != u->a[i0].x[0] || u->a[i].size != u->a[i0].size)
			aux[n_aux].st = i0, aux[n_aux++].en = i, i0 = i;

	m_a = max_occ > batch_size? max_occ : batch_size; // max size of a[] and b[]
	a = Kmalloc(km, uint64_t, m_a);
	b = Kmalloc(km, sa_aux_t, m_a);
	for (i = 0, m = 0; i < n_aux; ++i) {
		const anchor_aux_t *p = &aux[i];
		const mb_sai_t *q = &u->a[p->st];
		if (q->size + m > batch_size) {
			process_batch(km, idx, aux, m, b, a, qlen, mt, u, v);
			m = 0;
		}
		if (q->size <= max_occ) { // get SA for all of them
			for (j = 0; j < q->size; ++j)
				b[m].a = q->x[0] + j, b[m++].i = i;
		} else { // sample up to max_occ
			int32_t n = 0;
			for (j = 0; j < q->size && n < max_occ; ++n) {
				int32_t step = (q->size - j) / (max_occ - n);
				if (step < 1) step = 1;
				b[m].a = q->x[0] + j, b[m++].i = i;
				j += step;
			}
		}
		assert(m <= m_a); // shouldn't happen!
	}
	process_batch(km, idx, aux, m, b, a, qlen, mt, u, v);
	kfree(km, b);
	kfree(km, a);
	kfree(km, aux);

	if (mt != L2B_METH_NONE && v->n > 0) {
		int64_t t0, t1;
		for (i = 0, t0 = 0; i < v->n; ++i) t0 += v->a[i].len;
		mb_anchor_split_meth(km, idx->l2b, min_len, qlen, qseq, mt, v);
		for (i = 0, t1 = 0; i < v->n; ++i) t1 += v->a[i].len;
		seed_ratio = (double)t1 / t0;
	}

	/* ALT-seed -> primary anchor projection (segdup recovery).  Inject in the
	 * concatenated frame so the radix_sort + tpos-rebase + mb_anchor_dedup below
	 * order and dedup the injected anchors for free.  No-op unless a seed landed
	 * on an ALT contig (so non-ALT references are byte-identical). */
	mb_anchor_project_alt(km, idx->l2b, qlen, v);

	radix_sort_mb_anchor(v->a, v->a + v->n);
	for (i = 0; i < v->n; ++i) { // adjust mb_anchor_t::tpos
		mb_anchor_t *q = &v->a[i];
		const l2b_ctg_t *ctg = &idx->l2b->ctg[q->sid>>1];
		q->tpos -= ctg->off * 2 + ctg->len * (q->sid&1);
	}
	mb_anchor_dedup(v);
	return seed_ratio;
}

void mb_anchor_sort(const l2b_t *l2b, int64_t n_a, mb_anchor_t *a)
{
	int64_t i;
	if (n_a <= 1) return;
	for (i = 0; i < n_a; ++i) {
		const l2b_ctg_t *ctg = &l2b->ctg[a[i].sid>>1];
		a[i].tpos += ctg->off * 2 + ctg->len * (a[i].sid&1);
	}
	radix_sort_mb_anchor(a, a + n_a);
	for (i = 0; i < n_a; ++i) {
		const l2b_ctg_t *ctg = &l2b->ctg[a[i].sid>>1];
		a[i].tpos -= ctg->off * 2 + ctg->len * (a[i].sid&1);
	}
}
