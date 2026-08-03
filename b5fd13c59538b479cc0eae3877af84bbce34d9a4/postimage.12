#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "mbpriv.h"
#include "kommon.h"

static char mb_rg_id[256];

/******************
 * string helpers *
 ******************/

static inline void str_enlarge(kstring_t *s, size_t l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kom_roundup64(s->m);
		s->s = kom_realloc(char, s->s, s->m);
	}
}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

// inline equivalents of the kom_sprintf_lite "%c"/"%s"/"%d"/"%u" conversions,
// without the per-field format-string parsing
static inline void str_putc(kstring_t *s, char c)
{
	str_enlarge(s, 1);
	s->s[s->l++] = c;
}

static inline void str_puts(kstring_t *s, const char *t)
{
	str_copy(s, t, t + strlen(t));
}

static inline void str_puti(kstring_t *s, long c) // matches "%d" and "%ld"
{
	char buf[24];
	int i, l = 0;
	unsigned long x = c >= 0? (unsigned long)c : -(unsigned long)c;
	do { buf[l++] = x % 10 + '0'; x /= 10; } while (x > 0);
	if (c < 0) buf[l++] = '-';
	str_enlarge(s, l);
	for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
}

static inline void str_putu(kstring_t *s, unsigned x) // matches "%u"
{
	char buf[16];
	int i, l = 0;
	do { buf[l++] = x % 10 + '0'; x /= 10; } while (x > 0);
	str_enlarge(s, l);
	for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
}

/**************
 * PAF output *
 **************/

static inline void write_tags(kstring_t *s, const mb_hit_t *p)
{
	int32_t nm = p->blen - p->mlen + p->p->n_ambi;
	str_puts(s, "\tNM:i:"); str_puti(s, nm);
	str_puts(s, "\tAS:i:"); str_puti(s, p->p->dp_score);
	str_puts(s, "\tms:i:"); str_puti(s, p->p->dp_max0);
	str_puts(s, "\tmd:i:"); str_puti(s, p->p->dp_max - p->p->dp_max2);
}

/* Build the Bismark XM:Z per-base methylation-call string, a NUL-terminated
 * string of length qlen written in SAM SEQ orientation. SAM SEQ is always in
 * top-strand-ref orientation (sam_write_sq revcomps when r->rev), so a
 * CIGAR-ordered walk over SEQ positions produces a string that already aligns
 * with SEQ in display order; no end-of-walk reversal is needed.
 *
 * Alphabet: . non-cytosine; z/Z CpG un/methylated; x/X CHG; h/H CHH;
 * u/U cytosine in unknown context (off the reference window or N base).
 *
 * Returns a malloc()ed string; the caller frees it. This is the expensive half
 * of the methylation tags -- a reference window fetch plus a per-base walk, and
 * a result that is as long as the read -- so it lives in its own function that
 * write_meth_tags only calls when XM is actually selected. */
static char *build_meth_xm(const l2b_t *l2b, const mb_bseq1_t *t, int xg_is_ga, const mb_hit_t *r)
{
	int qlen = t->l_seq;
	int64_t tlen = l2b->ctg[r->tid].len;
	int64_t fetch_st = r->ts >= 2 ? r->ts - 2 : 0;
	int64_t fetch_en = r->te + 2 <= tlen ? r->te + 2 : tlen;
	uint8_t *ref = (uint8_t*)kom_malloc(uint8_t, fetch_en - fetch_st);
	char *xm = (char*)kom_malloc(char, qlen + 1);
	/* XM is written in SAM SEQ orientation, so the first aligned base lands
	 * at the SAM leading soft-clip offset. That offset is r->qs for forward
	 * reads but qlen - r->qe for reverse reads (see clip_len[0] in mb_fmt_sam);
	 * r->qs/r->qe are in original-read coordinates. Using r->qs unconditionally
	 * mis-frames reverse reads with asymmetric clipping. */
	int seq_pos = r->rev ? (int)(qlen - r->qe) : r->qs;
	int64_t ref_pos = r->ts;
	uint32_t k;
	int j;

	l2b_getseq(l2b, r->tid, fetch_st, fetch_en, ref);
	for (j = 0; j < qlen; ++j) xm[j] = '.';
	xm[qlen] = 0;

	for (k = 0; k < r->p->n_cigar; ++k) {
		int op = r->p->cigar[k] & 0xf;
		int len = r->p->cigar[k] >> 4;
		if (op == 0 || op == 7 || op == 8) { /* M/=/X */
			for (j = 0; j < len; ++j) {
				int read_b, off, c1, c2;
				char meth = '.', unmeth = '.', call = '.';
				if (ref_pos < fetch_st || ref_pos >= fetch_en) goto next;
				off = (int)(ref_pos - fetch_st);
				/* Top-strand-ref-frame read base. SEQ at seq_pos is
				 * t->seq[seq_pos] for r->rev=0 and comp(t->seq[qlen-1-seq_pos])
				 * for r->rev=1; either way it aligns to top_strand_ref[ref_pos]. */
				if (r->rev) {
					int b = kom_nt4_table[(uint8_t)t->seq[qlen - 1 - seq_pos]];
					read_b = b < 4 ? 3 - b : 4;
				} else {
					read_b = kom_nt4_table[(uint8_t)t->seq[seq_pos]];
				}
				if (!xg_is_ga) {
					/* XG=CT: top-strand C; context on top strand. */
					if (ref[off] != 1) goto next;
					c1 = (off + 1 < fetch_en - fetch_st) ? ref[off + 1] : 4;
					c2 = (off + 2 < fetch_en - fetch_st) ? ref[off + 2] : 4;
					if (c1 == 2)        { meth = 'Z'; unmeth = 'z'; }
					else if (c1 == 4)   { meth = 'U'; unmeth = 'u'; }
					else if (c2 == 2)   { meth = 'X'; unmeth = 'x'; }
					else if (c2 == 4)   { meth = 'U'; unmeth = 'u'; }
					else                { meth = 'H'; unmeth = 'h'; }
					if      (read_b == 1) call = meth;
					else if (read_b == 3) call = unmeth;
				} else {
					/* XG=GA: bottom-strand C, top-strand G. Bottom-strand
					 * 5'->3' context base is comp(top[ref_pos-1]); CpG iff
					 * top[ref_pos-1]=C, etc. */
					if (ref[off] != 2) goto next;
					c1 = (ref_pos - 1 >= fetch_st) ? ref[off - 1] : 4;
					c2 = (ref_pos - 2 >= fetch_st) ? ref[off - 2] : 4;
					if (c1 == 1)        { meth = 'Z'; unmeth = 'z'; }
					else if (c1 == 4)   { meth = 'U'; unmeth = 'u'; }
					else if (c2 == 1)   { meth = 'X'; unmeth = 'x'; }
					else if (c2 == 4)   { meth = 'U'; unmeth = 'u'; }
					else                { meth = 'H'; unmeth = 'h'; }
					/* meth bottom C preserved -> comp = G (read_b=2);
					 * unmeth bottom C->T -> comp = A (read_b=0). */
					if      (read_b == 2) call = meth;
					else if (read_b == 0) call = unmeth;
				}
next:
				xm[seq_pos] = call;
				++seq_pos;
				++ref_pos;
			}
		} else if (op == 1) { /* I */
			seq_pos += len;
		} else if (op == 2 || op == 3) { /* D, N */
			ref_pos += len;
		}
	}

	free(ref);
	return xm;
}

/* Bismark XR/XG/XM. Under directional --meth: XR follows mate (R1=CT,
 * R2=GA, SE=CT); XG = "GA" iff (is_R2 XOR r->rev), reproducing Bismark's
 * (XR,XG) -> {OT,OB,CTOT,CTOB} encoding.
 *
 * `tags` is an MB_METH_TAG_* mask (see --meth-tags). It only ever removes
 * tags: with the default MB_METH_TAG_ALL the emitted bytes are exactly what
 * this function wrote before the mask existed, in the same XR/XG/XM order. */
static void write_meth_tags(kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t,
							int n_seg, int seg_idx, const mb_hit_t *r, int32_t tags)
{
	int is_r2 = (n_seg == 2 && seg_idx == 1);
	int xg_is_ga = is_r2 ^ r->rev;
	if (tags & MB_METH_TAG_XR) kom_sprintf_lite(s, "\tXR:Z:%s", is_r2 ? "GA" : "CT");
	if (tags & MB_METH_TAG_XG) kom_sprintf_lite(s, "\tXG:Z:%s", xg_is_ga ? "GA" : "CT");
	if (tags & MB_METH_TAG_XM) {
		char *xm = build_meth_xm(l2b, t, xg_is_ga, r);
		kom_sprintf_lite(s, "\tXM:Z:%s", xm);
		free(xm);
	}
}

void mb_fmt_paf(kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, const mb_hit_t *p, uint64_t opt_flag, int n_seg, int seg_idx)
{
	str_puts(s, t->name);
	if (n_seg > 1 && seg_idx >= 0) {
		str_putc(s, '/'); str_puti(s, seg_idx + 1);
	}
	str_putc(s, '\t'); str_puti(s, t->l_seq);
	if (p == 0) { // for unmapped reads
		str_puts(s, "\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\n");
		s->s[s->l] = 0; // keep the always-room-for-an-extra-byte NUL invariant
		return;
	}
	str_putc(s, '\t'); str_puti(s, p->qs);
	str_putc(s, '\t'); str_puti(s, p->qe);
	str_putc(s, '\t'); str_putc(s, p->rev? '-' : '+');
	str_putc(s, '\t'); str_puts(s, l2b->ctg[p->tid].name);
	str_putc(s, '\t'); str_puti(s, l2b->ctg[p->tid].len);
	str_putc(s, '\t'); str_puti(s, p->ts);
	str_putc(s, '\t'); str_puti(s, p->te);
	str_putc(s, '\t'); str_puti(s, p->mlen);
	str_putc(s, '\t'); str_puti(s, p->blen);
	str_putc(s, '\t'); str_puti(s, p->mapq);
	str_puts(s, "\ttp:A:"); str_putc(s, p->parent == p->id? 'P' : 'S');
	str_puts(s, "\ts1:i:"); str_puti(s, p->score);
	str_puts(s, "\tcm:i:"); str_puti(s, p->cnt);
	if (p->parent == p->id) { str_puts(s, "\ts2:i:"); str_puti(s, p->subsc >= 0? p->subsc : 0); }
	if (p->p) {
		write_tags(s, p);
		if (p->p->n_cigar > 0) {
			int32_t i;
			str_puts(s, "\tcg:Z:");
			for (i = 0; i < p->p->n_cigar; ++i) {
				str_puti(s, p->p->cigar[i]>>4); str_putc(s, MB_CIGAR_STR[p->p->cigar[i]&0xf]);
			}
		}
		if (p->p->cs) { str_putc(s, '\t'); str_puts(s, (char*)&p->p->cigar[p->p->n_cigar]); }
	}
	if ((opt_flag & MB_F_COPY_COMMENT) && t->comment) {
		str_putc(s, '\t'); str_puts(s, t->comment);
	}
	str_putc(s, '\n');
	s->s[s->l] = 0; // keep the always-room-for-an-extra-byte NUL invariant (matches mb_fmt_sam)
}

/**************
 * SAM header *
 **************/

char *mb_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == 'n') *q++ = '\n';
			else if (*p == 'r') *q++ = '\r';
			else if (*p == '\\') *q++ = '\\';
			else if (*p == '\0') break;
		} else *q++ = *p;
	}
	*q = '\0';
	return s;
}

static int sam_write_rg_line(kstring_t *str, const char *s)
{
	char *p, *q, *r, *rg_line = 0;
	memset(mb_rg_id, 0, 256);
	if (s == 0) return 0;
	if (strstr(s, "@RG") != s) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] the read group line is not started with @RG\n");
		goto err_set_rg;
	}
	if (strstr(s, "\t") != NULL) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] the read group line contained literal <tab> characters -- replace with escaped tabs: \\t\n");
		goto err_set_rg;
	}
	rg_line = kom_strdup(s);
	mb_escape(rg_line);
	if ((p = strstr(rg_line, "\tID:")) == 0) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] no ID within the read group line\n");
		goto err_set_rg;
	}
	p += 4;
	for (q = p; *q && *q != '\t' && *q != '\n'; ++q);
	if (q - p + 1 > 256) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] @RG:ID is longer than 255 characters\n");
		goto err_set_rg;
	}
	for (q = p, r = mb_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
		*r++ = *q;
	kom_sprintf_lite(str, "%s\n", rg_line);
	free(rg_line);
	return 0;

err_set_rg:
	free(rg_line);
	return -1;
}

int mb_fmt_sam_hdr(kstring_t *str, const l2b_t *idx, const char *rg, const char *ver, int argc, char *argv[])
{
	int i, ret = 0;
	str->l = 0;
	kom_sprintf_lite(str, "@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
	if (idx)
		for (i = 0; i < idx->n_ctg; ++i)
			kom_sprintf_lite(str, "@SQ\tSN:%s\tLN:%ld\n", idx->ctg[i].name, idx->ctg[i].len);
	if (rg) ret = sam_write_rg_line(str, rg);
	kom_sprintf_lite(str, "@PG\tID:minibwa\tPN:minibwa");
	if (ver) kom_sprintf_lite(str, "\tVN:%s", ver);
	if (argc > 1) {
		kom_sprintf_lite(str, "\tCL:minibwa");
		for (i = 0; i < argc; ++i)
			kom_sprintf_lite(str, " %s", argv[i]);
	}
	kom_sprintf_lite(str, "\n");
	return ret;
}

/**************
 * SAM output *
 **************/

static void sam_write_sq(kstring_t *s, char *seq, int l, int rev, int comp)
{
	if (rev) {
		int i;
		str_enlarge(s, l);
		for (i = 0; i < l; ++i) {
			int c = seq[l - 1 - i];
			s->s[s->l + i] = c < 128 && comp? kom_comp_table[c] : c;
		}
		s->l += l;
	} else str_copy(s, seq, seq + l);
}

static inline const mb_hit_t *get_sam_pri(int n_hit, const mb_hit_t *hit)
{
	int i;
	for (i = 0; i < n_hit; ++i)
		if (hit[i].sam_pri)
			return &hit[i];
	assert(n_hit == 0);
	return NULL;
}

static void write_sam_cigar(kstring_t *s, int sam_flag, int in_tag, int qlen, const mb_hit_t *r, int64_t opt_flag)
{
	if (r->p == 0) {
		str_putc(s, '*');
	} else {
		uint32_t k, clip_len[2];
		clip_len[0] = r->rev? qlen - r->qe : r->qs;
		clip_len[1] = r->rev? r->qs : qlen - r->qe;
		if (in_tag) {
			int clip_char = (((sam_flag&0x800) || ((sam_flag&0x100) && (opt_flag&MB_F_2ND_SEQ))) &&
							 !(opt_flag&MB_F_SUPP_SOFT)) ? 5 : 4;
			str_puts(s, "\tCG:B:I");
			if (clip_len[0]) { str_putc(s, ','); str_putu(s, clip_len[0]<<4|clip_char); }
			for (k = 0; k < r->p->n_cigar; ++k) {
				str_putc(s, ','); str_putu(s, r->p->cigar[k]);
			}
			if (clip_len[1]) { str_putc(s, ','); str_putu(s, clip_len[1]<<4|clip_char); }
		} else {
			int clip_char = (((sam_flag&0x800) || ((sam_flag&0x100) && (opt_flag&MB_F_2ND_SEQ))) &&
							 !(opt_flag&MB_F_SUPP_SOFT)) ? 'H' : 'S';
			assert(clip_len[0] < qlen && clip_len[1] < qlen);
			if (clip_len[0]) { str_puti(s, clip_len[0]); str_putc(s, clip_char); }
			for (k = 0; k < r->p->n_cigar; ++k) {
				str_puti(s, r->p->cigar[k]>>4); str_putc(s, MB_CIGAR_STR[r->p->cigar[k]&0xf]);
			}
			if (clip_len[1]) { str_puti(s, clip_len[1]); str_putc(s, clip_char); }
		}
	}
}

void mb_fmt_sam(void *km, kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, const mb_opt_t *opt, int seg_idx, int32_t mate_qlen)
{
	int flag, n_h = n_hit[seg_idx];
	int this_tid = -1, this_pos = -1;
	const mb_hit_t *h = hit[seg_idx], *r_prev = NULL, *r_next;
	const mb_hit_t *r = n_h > 0 && hit_idx < n_h && hit_idx >= 0? &h[hit_idx] : NULL;

	assert(n_seg == 1 || n_seg == 2);

	// find the primary of the previous and the next segments, if they are mapped
	if (n_seg > 1) {
		int next_sid = (seg_idx + 1) % n_seg;
		r_prev = r_next = get_sam_pri(n_hit[next_sid], hit[next_sid]);
	} else r_prev = r_next = NULL;

	// write QNAME and FLAG
	str_puts(s, t->name);
	flag = n_seg > 1? 0x1 : 0x0;
	if (r == 0) {
		flag |= 0x4;
	} else {
		if (r->rev) flag |= 0x10;
		if (r->parent != r->id) flag |= 0x100;
		else if (!r->sam_pri) flag |= 0x800;
	}
	if (n_seg > 1) {
		if (r && r->proper_pair) flag |= 0x2;
		if (seg_idx == 0) flag |= 0x40;
		else if (seg_idx == n_seg - 1) flag |= 0x80;
		if (r_next == NULL) flag |= 0x8;
		else if (r_next->rev) flag |= 0x20;
	}
	str_putc(s, '\t'); str_puti(s, flag);

	// write coordinate, MAPQ and CIGAR
	if (r == 0) {
		if (r_prev) {
			this_tid = r_prev->tid, this_pos = r_prev->ts;
			str_putc(s, '\t'); str_puts(s, l2b->ctg[this_tid].name); str_putc(s, '\t'); str_puti(s, this_pos+1); str_puts(s, "\t0\t*");
		} else str_puts(s, "\t*\t0\t0\t*");
	} else {
		this_tid = r->tid, this_pos = r->ts;
		str_putc(s, '\t'); str_puts(s, l2b->ctg[r->tid].name); str_putc(s, '\t'); str_puti(s, r->ts+1); str_putc(s, '\t'); str_puti(s, r->mapq); str_putc(s, '\t');
		write_sam_cigar(s, flag, 0, t->l_seq, r, opt->flag);
	}

	// write mate positions
	if (n_seg > 1) {
		int tlen = 0;
		if (this_tid >= 0 && r_next) {
			if (this_tid == r_next->tid) {
				if (r) {
					int this_pos5 = r->rev? r->te - 1 : this_pos;
					int next_pos5 = r_next->rev? r_next->te - 1 : r_next->ts;
					tlen = next_pos5 - this_pos5;
				}
				str_puts(s, "\t=\t");
			} else { str_putc(s, '\t'); str_puts(s, l2b->ctg[r_next->tid].name); str_putc(s, '\t'); }
			str_puti(s, r_next->ts + 1); str_putc(s, '\t');
		} else if (r_next) { // && this_tid < 0
			str_putc(s, '\t'); str_puts(s, l2b->ctg[r_next->tid].name); str_putc(s, '\t'); str_puti(s, r_next->ts + 1); str_putc(s, '\t');
		} else if (this_tid >= 0) { // && r_next == NULL
			str_puts(s, "\t=\t"); str_puti(s, this_pos + 1); str_putc(s, '\t'); // next segment will take r's coordinate
		} else str_puts(s, "\t*\t0\t"); // neither has coordinates
		if (tlen > 0) ++tlen;
		else if (tlen < 0) --tlen;
		str_puti(s, tlen); str_putc(s, '\t');
	} else str_puts(s, "\t*\t0\t0\t");

	// write SEQ and QUAL
	if (r == 0) {
		sam_write_sq(s, t->seq, t->l_seq, 0, 0);
		str_putc(s, '\t');
		if (t->qual) sam_write_sq(s, t->qual, t->l_seq, 0, 0);
		else str_putc(s, '*');
	} else {
		if ((flag & 0x900) == 0 || (opt->flag & MB_F_SUPP_SOFT)
		    || ((opt->flag & MB_F_ALT_RECORDS) && r->is_alt)) {
			sam_write_sq(s, t->seq, t->l_seq, r->rev, r->rev);
			str_putc(s, '\t');
			if (t->qual) sam_write_sq(s, t->qual, t->l_seq, r->rev, 0);
			else str_putc(s, '*');
		} else if ((flag & 0x100) && !(opt->flag & MB_F_2ND_SEQ)){
			str_puts(s, "*\t*");
		} else {
			sam_write_sq(s, t->seq + r->qs, r->qe - r->qs, r->rev, r->rev);
			str_putc(s, '\t');
			if (t->qual) sam_write_sq(s, t->qual + r->qs, r->qe - r->qs, r->rev, 0);
			else str_putc(s, '*');
		}
	}

	// write tags
	if (mb_rg_id[0]) { str_puts(s, "\tRG:Z:"); str_puts(s, mb_rg_id); }
	if (n_seg > 2) { str_puts(s, "\tFI:i:"); str_puti(s, seg_idx); }
	if (r) {
		write_tags(s, r);
		if ((opt->flag & MB_F_METH) && opt->meth_tags)
			write_meth_tags(s, l2b, t, n_seg, seg_idx, r, opt->meth_tags);
		// MC:Z mate CIGAR and MQ:i mate MAPQ; r_next is the mate's primary (see above).
		if (n_seg > 1 && r_next && r_next->p && r_next->p->n_cigar > 0 && mate_qlen > 0) {
			str_puts(s, "\tMC:Z:");
			write_sam_cigar(s, 0, 0, mate_qlen, r_next, opt->flag);
			str_puts(s, "\tMQ:i:"); str_puti(s, r_next->mapq);
		}
		if (r->p->cs) { str_putc(s, '\t'); str_puts(s, (char*)&r->p->cigar[r->p->n_cigar]); }
		if (r->parent == r->id && r->p && n_h > 1 && h && r >= h && r - h < n_h) { // supplementary aln may exist
			int i, n_sa = 0; // n_sa: number of SA fields
			for (i = 0; i < n_h; ++i)
				if (i != r - h && h[i].parent == h[i].id && h[i].p)
					++n_sa;
			if (n_sa > 0) {
				str_puts(s, "\tSA:Z:");
				for (i = 0; i < n_h; ++i) {
					const mb_hit_t *q = &h[i];
					int l_M, l_I = 0, l_D = 0, clip5 = 0, clip3 = 0;
					if (r == q || q->parent != q->id || q->p == 0) continue;
					if (q->qe - q->qs < q->te - q->ts) l_M = q->qe - q->qs, l_D = (q->te - q->ts) - l_M;
					else l_M = q->te - q->ts, l_I = (q->qe - q->qs) - l_M;
					clip5 = q->rev? t->l_seq - q->qe : q->qs;
					clip3 = q->rev? q->qs : t->l_seq - q->qe;
					str_puts(s, l2b->ctg[q->tid].name); str_putc(s, ','); str_puti(s, q->ts+1); str_putc(s, ','); str_putc(s, "+-"[q->rev]); str_putc(s, ',');
					if (clip5) { str_puti(s, clip5); str_putc(s, 'S'); }
					if (l_M) { str_puti(s, l_M); str_putc(s, 'M'); }
					if (l_I) { str_puti(s, l_I); str_putc(s, 'I'); }
					if (l_D) { str_puti(s, l_D); str_putc(s, 'D'); }
					if (clip3) { str_puti(s, clip3); str_putc(s, 'S'); }
					str_putc(s, ','); str_puti(s, q->mapq); str_putc(s, ','); str_puti(s, q->blen - q->mlen + q->p->n_ambi); str_putc(s, ';');
				}
			}
			if (opt->xa_max > 0) {
				int i, n_xa = 0;
				for (i = 0; i < n_h; ++i)
					if (i != r - h && h[i].parent == r - h && h[i].p->dp_max >= (double)opt->out_s * r->p->dp_max)
						++n_xa;
				if (n_xa > 0) kom_sprintf_lite(s, "\tn2:i:%d", n_xa);
				if (n_xa > 0 && n_xa <= opt->xa_max) {
					kom_sprintf_lite(s, "\tXA:Z:");
					for (i = 0; i < n_h; ++i) {
						const mb_hit_t *q = &h[i];
						if (i != r - h && q->parent == r - h && q->p->dp_max >= (double)opt->out_s * r->p->dp_max) {
							kom_sprintf_lite(s, "%s,%c%d,", l2b->ctg[q->tid].name, "+-"[q->rev], q->ts+1);
							write_sam_cigar(s, 0, 0, t->l_seq, q, opt->flag);
							kom_sprintf_lite(s, ",%d;", q->blen - q->mlen + q->p->n_ambi);
						}
					}
				}
			}
		}
	}

	if ((opt->flag & MB_F_COPY_COMMENT) && t->comment) {
		str_putc(s, '\t'); str_puts(s, t->comment);
	}
	str_putc(s, '\n');
	s->s[s->l] = 0; // we always have room for an extra byte
}

void mb_format(void *km, kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, const mb_opt_t *opt, int seg_idx, int32_t mate_qlen)
{
	if (!(opt->flag & MB_F_PAF))
		mb_fmt_sam(km, s, l2b, t, n_seg, n_hit, hit, hit_idx, opt, seg_idx, mate_qlen);
	else
		mb_fmt_paf(s, l2b, t, hit_idx >= 0? &hit[seg_idx][hit_idx] : 0, opt->flag, n_seg, seg_idx);
}
