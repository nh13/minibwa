#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define __STDC_LIMIT_MACROS
#include "bseq.h"
#include "fr_fastq.h"
#include <zlib.h>

/* Byte source for the fr_fastq parser: stock zlib gzread. gzread transparently
 * passes plaintext through when the input is not gzip, so the same path serves
 * an in-process .gz file and a decompressed plaintext pipe on stdin. */
static int mb_src_read(void *ctx, unsigned char *buf, int len) { return gzread((gzFile)ctx, buf, len); }

#define kvec_t(type) struct { size_t n, m; type *a; }

#define kv_resize(type, v, s) do { \
		if ((v).m < (s)) { \
			(v).m = (s); \
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
		} \
	} while (0)

#define kv_push(type, v, x) do { \
		if ((v).n == (v).m) { \
			(v).m += ((v).m>>1) + 16; \
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
		} \
		(v).a[(v).n++] = (x); \
	} while (0)

#define kv_pushp(type, v, p) do { \
		if ((v).n == (v).m) { \
			(v).m += ((v).m>>1) + 16; \
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
		} \
		*(p) = &(v).a[(v).n++]; \
	} while (0)

#define CHECK_PAIR_THRES 1000000

struct mb_bseq_file_s {
	gzFile fp;
	fr_fastq_t *parser;
	mb_bseq1_t s;
};

mb_bseq_file_t *mb_bseq_open(const char *fn)
{
	mb_bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (f == 0) return 0;
	fp = (mb_bseq_file_t*)calloc(1, sizeof(mb_bseq_file_t));
	fp->fp = f;
	fp->parser = fr_fastq_init(mb_src_read, fp->fp);
	return fp;
}

void mb_bseq_close(mb_bseq_file_t *fp)
{
	fr_fastq_destroy(fp->parser);
	gzclose(fp->fp);
	free(fp->s.name); free(fp->s.seq); free(fp->s.qual); free(fp->s.comment);
	free(fp);
}

/* Length-aware dup: malloc len+1 and NUL-terminate (fr_fastq slices are not
 * NUL-terminated). Mirrors the storage kstrdup produced from a kstring_t. */
static inline char *fr_dup(const char *s, size_t l)
{
	char *t = (char*)malloc(l + 1);
	memcpy(t, s, l);
	t[l] = 0;
	return t;
}

/* Copy one parsed record into mb_bseq1_t. Byte-for-byte the same as the old
 * kseq2bseq: warn on empty name, own all strings, convert U->T via the same
 * `--c` decrement ('U'->'T', 'u'->'t'), and store comment/qual only when
 * present and requested. */
static inline void fr_rec_to_bseq1(const fr_fastq_rec_t *r, mb_bseq1_t *s, int with_qual, int with_comment)
{
	int i;
	if (r->name_l == 0)
		fprintf(stderr, "[WARNING]\033[1;31m empty sequence name in the input.\033[0m\n");
	s->name = fr_dup(r->name, r->name_l);
	s->seq = fr_dup(r->seq, r->seq_l);
	for (i = 0; i < (int)r->seq_l; ++i) // convert U to T
		if (s->seq[i] == 'u' || s->seq[i] == 'U')
			--s->seq[i];
	s->qual = with_qual && r->qual_l? fr_dup(r->qual, r->qual_l) : 0;
	s->comment = with_comment && r->comment_l? fr_dup(r->comment, r->comment_l) : 0;
	s->l_seq = r->seq_l;
}

mb_bseq1_t *mb_bseq_read(mb_bseq_file_t *fp, int64_t chunk_size, int with_qual, int with_comment, int frag_mode, int min_cnt, int64_t max_chunk_size, int *n_)
{
	int64_t size = 0;
	int ret;
	fr_fastq_rec_t rec;
	kvec_t(mb_bseq1_t) a = {0,0,0};
	fr_fastq_t *parser = fp->parser;
	*n_ = 0;
	if (fp->s.seq) {
		kv_resize(mb_bseq1_t, a, 256);
		kv_push(mb_bseq1_t, a, fp->s);
		size = fp->s.l_seq;
		memset(&fp->s, 0, sizeof(mb_bseq1_t));
	}
	if (max_chunk_size < chunk_size)
		max_chunk_size = chunk_size;
	for (;;) {
		ret = fr_fastq_next(parser, &rec);
		if (ret <= 0) break;
		int32_t to_stop = 0;
		mb_bseq1_t *s;
		assert(rec.seq_l <= INT32_MAX);
		if (a.m == 0) kv_resize(mb_bseq1_t, a, 256);
		kv_pushp(mb_bseq1_t, a, &s);
		fr_rec_to_bseq1(&rec, s, with_qual, with_comment);
		size += rec.seq_l;
		if (chunk_size <= 0 || max_chunk_size <= 0) to_stop = 1;
		else if (size >= max_chunk_size) to_stop = 1;
		else if (size >= chunk_size && a.n >= min_cnt) to_stop = 1;
		if (to_stop) {
			if (frag_mode && a.a[a.n-1].l_seq < CHECK_PAIR_THRES) {
				while ((ret = fr_fastq_next(parser, &rec)) == 1) {
					fr_rec_to_bseq1(&rec, &fp->s, with_qual, with_comment);
					size += rec.seq_l;
					if (mb_qname_same(fp->s.name, a.a[a.n-1].name)) {
						kv_push(mb_bseq1_t, a, fp->s);
						memset(&fp->s, 0, sizeof(mb_bseq1_t));
					} else break;
				}
			}
			break;
		}
	}
	if (ret == -2) {
		if (a.n) fprintf(stderr, "[WARNING]\033[1;31m failed to parse the FASTA/FASTQ record next to '%s'. Continue anyway.\033[0m\n", a.a[a.n-1].name);
		else fprintf(stderr, "[WARNING]\033[1;31m failed to parse the first FASTA/FASTQ record. Continue anyway.\033[0m\n");
	}
	*n_ = a.n;
	return a.a;
}

mb_bseq1_t *mb_bseq_read_frag(int n_fp, mb_bseq_file_t **fp, int64_t chunk_size, int with_qual, int with_comment, int *n_)
{
	int i;
	int64_t size = 0;
	fr_fastq_rec_t *rec;
	kvec_t(mb_bseq1_t) a = {0,0,0};
	*n_ = 0;
	if (n_fp < 1) return 0;
	/* One record slot per file: each file has its own parser, so the slices
	 * stay valid across the per-file reads until we copy them out below. */
	rec = (fr_fastq_rec_t*)calloc(n_fp, sizeof(fr_fastq_rec_t));
	while (1) {
		int n_read = 0;
		for (i = 0; i < n_fp; ++i)
			if (fr_fastq_next(fp[i]->parser, &rec[i]) == 1)
				++n_read;
		if (n_read < n_fp) {
			if (n_read > 0)
				fprintf(stderr, "[W::%s]\033[1;31m query files have different number of records; extra records skipped.\033[0m\n", __func__);
			break; // some file reaches the end
		}
		if (a.m == 0) kv_resize(mb_bseq1_t, a, 256);
		for (i = 0; i < n_fp; ++i) {
			mb_bseq1_t *s;
			kv_pushp(mb_bseq1_t, a, &s);
			fr_rec_to_bseq1(&rec[i], s, with_qual, with_comment);
			size += s->l_seq;
		}
		if (size >= chunk_size) break;
	}
	free(rec);
	*n_ = a.n;
	return a.a;
}

int mb_bseq_eof(mb_bseq_file_t *fp)
{
	return (fr_fastq_eof(fp->parser) && fp->s.seq == 0);
}
