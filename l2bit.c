#include <zlib.h>
#include <stdio.h>
#include <assert.h>
#include "kommon.h"
#include "l2bit.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static int64_t l2b_pos2cid(const l2b_t *l2b, int64_t s, int64_t len, int64_t *cst)
{
	int64_t lo = 0, hi = l2b->n_ctg, mid;
	while (lo < hi) {
		const l2b_ctg_t *ctg;
		mid = (lo + hi) / 2;
		ctg = &l2b->ctg[mid];
		if (ctg->off <= s && s < ctg->off + ctg->len) {
			*cst = s - ctg->off;
			return s + len <= ctg->off + ctg->len? mid : -1;
		} else if (s < ctg->off) hi = mid;
		else lo = mid + 1;
	}
	return -1; // s is in no contig
}

int64_t l2b_intv2cid(const l2b_t *l2b, uint64_t st, uint64_t en, int64_t *cst, int *rev)
{
	int64_t s;
	assert(st < en);
	if (en > l2b->tot_len * 2) return -3;
	if (st < l2b->tot_len && l2b->tot_len < en) return -2;
	*rev = (st >= l2b->tot_len);
	s = st < l2b->tot_len? st : l2b->tot_len * 2 - en;
	return l2b_pos2cid(l2b, s, en - st, cst);
}

int64_t l2b_intv2cid_meth(const l2b_t *l2b, uint64_t st, uint64_t en, l2b_meth_t *mt, int64_t *cst, int *rev)
{
	int64_t s, len = en - st;
	int32_t copy;
	uint64_t tot_len = l2b->tot_len;

	assert(st < en);
	if (en > tot_len * 4) return -3;
	copy = st / tot_len; // 0: c2t_f, 1: g2a_f, 2: g2a_r, 3: c2t_r
	*mt = (copy == 0 || copy == 3)? L2B_METH_C2T : L2B_METH_G2A;
	*rev = (copy >= 2);
	s = st - tot_len * copy;
	if (copy >= 2) s = tot_len - len - s; // flip for reverse copies
	return s < 0? -2 : l2b_pos2cid(l2b, s, len, cst);
}

// Maps a packed byte (4 x 2-bit bases, base 0 in the low bits) to its 4
// unpacked bases, one output byte each. This is the bulk-unpack equivalent of
// l2b_get0() (the canonical per-base accessor): entry b holds the four bases
// l2b_get0 would return for the four positions packed in byte b, ordered so a
// single little-endian uint32 store emits them in ascending position order.
// Compile-time constant -- no runtime constructor. The uint32 byte order is
// load-bearing (a memcpy of the entry writes base 0 to the lowest address), so
// the little-endian assumption is enforced at compile time rather than left to
// produce silently wrong output on a big-endian host:
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
#error "l2b_getseq LUT (g_l2b_unpack4) assumes a little-endian host"
#endif
static const uint32_t g_l2b_unpack4[256] = {
	0x00000000u, 0x00000001u, 0x00000002u, 0x00000003u, 0x00000100u, 0x00000101u, 0x00000102u, 0x00000103u,
	0x00000200u, 0x00000201u, 0x00000202u, 0x00000203u, 0x00000300u, 0x00000301u, 0x00000302u, 0x00000303u,
	0x00010000u, 0x00010001u, 0x00010002u, 0x00010003u, 0x00010100u, 0x00010101u, 0x00010102u, 0x00010103u,
	0x00010200u, 0x00010201u, 0x00010202u, 0x00010203u, 0x00010300u, 0x00010301u, 0x00010302u, 0x00010303u,
	0x00020000u, 0x00020001u, 0x00020002u, 0x00020003u, 0x00020100u, 0x00020101u, 0x00020102u, 0x00020103u,
	0x00020200u, 0x00020201u, 0x00020202u, 0x00020203u, 0x00020300u, 0x00020301u, 0x00020302u, 0x00020303u,
	0x00030000u, 0x00030001u, 0x00030002u, 0x00030003u, 0x00030100u, 0x00030101u, 0x00030102u, 0x00030103u,
	0x00030200u, 0x00030201u, 0x00030202u, 0x00030203u, 0x00030300u, 0x00030301u, 0x00030302u, 0x00030303u,
	0x01000000u, 0x01000001u, 0x01000002u, 0x01000003u, 0x01000100u, 0x01000101u, 0x01000102u, 0x01000103u,
	0x01000200u, 0x01000201u, 0x01000202u, 0x01000203u, 0x01000300u, 0x01000301u, 0x01000302u, 0x01000303u,
	0x01010000u, 0x01010001u, 0x01010002u, 0x01010003u, 0x01010100u, 0x01010101u, 0x01010102u, 0x01010103u,
	0x01010200u, 0x01010201u, 0x01010202u, 0x01010203u, 0x01010300u, 0x01010301u, 0x01010302u, 0x01010303u,
	0x01020000u, 0x01020001u, 0x01020002u, 0x01020003u, 0x01020100u, 0x01020101u, 0x01020102u, 0x01020103u,
	0x01020200u, 0x01020201u, 0x01020202u, 0x01020203u, 0x01020300u, 0x01020301u, 0x01020302u, 0x01020303u,
	0x01030000u, 0x01030001u, 0x01030002u, 0x01030003u, 0x01030100u, 0x01030101u, 0x01030102u, 0x01030103u,
	0x01030200u, 0x01030201u, 0x01030202u, 0x01030203u, 0x01030300u, 0x01030301u, 0x01030302u, 0x01030303u,
	0x02000000u, 0x02000001u, 0x02000002u, 0x02000003u, 0x02000100u, 0x02000101u, 0x02000102u, 0x02000103u,
	0x02000200u, 0x02000201u, 0x02000202u, 0x02000203u, 0x02000300u, 0x02000301u, 0x02000302u, 0x02000303u,
	0x02010000u, 0x02010001u, 0x02010002u, 0x02010003u, 0x02010100u, 0x02010101u, 0x02010102u, 0x02010103u,
	0x02010200u, 0x02010201u, 0x02010202u, 0x02010203u, 0x02010300u, 0x02010301u, 0x02010302u, 0x02010303u,
	0x02020000u, 0x02020001u, 0x02020002u, 0x02020003u, 0x02020100u, 0x02020101u, 0x02020102u, 0x02020103u,
	0x02020200u, 0x02020201u, 0x02020202u, 0x02020203u, 0x02020300u, 0x02020301u, 0x02020302u, 0x02020303u,
	0x02030000u, 0x02030001u, 0x02030002u, 0x02030003u, 0x02030100u, 0x02030101u, 0x02030102u, 0x02030103u,
	0x02030200u, 0x02030201u, 0x02030202u, 0x02030203u, 0x02030300u, 0x02030301u, 0x02030302u, 0x02030303u,
	0x03000000u, 0x03000001u, 0x03000002u, 0x03000003u, 0x03000100u, 0x03000101u, 0x03000102u, 0x03000103u,
	0x03000200u, 0x03000201u, 0x03000202u, 0x03000203u, 0x03000300u, 0x03000301u, 0x03000302u, 0x03000303u,
	0x03010000u, 0x03010001u, 0x03010002u, 0x03010003u, 0x03010100u, 0x03010101u, 0x03010102u, 0x03010103u,
	0x03010200u, 0x03010201u, 0x03010202u, 0x03010203u, 0x03010300u, 0x03010301u, 0x03010302u, 0x03010303u,
	0x03020000u, 0x03020001u, 0x03020002u, 0x03020003u, 0x03020100u, 0x03020101u, 0x03020102u, 0x03020103u,
	0x03020200u, 0x03020201u, 0x03020202u, 0x03020203u, 0x03020300u, 0x03020301u, 0x03020302u, 0x03020303u,
	0x03030000u, 0x03030001u, 0x03030002u, 0x03030003u, 0x03030100u, 0x03030101u, 0x03030102u, 0x03030103u,
	0x03030200u, 0x03030201u, 0x03030202u, 0x03030203u, 0x03030300u, 0x03030301u, 0x03030302u, 0x03030303u,
};

int64_t l2b_getseq(const l2b_t *l2b, int64_t tid, int64_t st, int64_t en, uint8_t *seq)
{
	const l2b_ctg_t *ctg;
	int64_t i, aid;
	int32_t n_ambi;
	if (tid < 0 || tid >= l2b->n_ctg) return -1;
	ctg = &l2b->ctg[tid];
	if (st < 0) st = 0;
	if (en > ctg->len) en = ctg->len;
	st += ctg->off;
	en += ctg->off;
	// head: unpack scalar until i is 4-base (byte) aligned
	for (i = st; i < en && (i & 3); ++i)
		seq[i - st] = l2b_get0(l2b, i);
	// middle: unpack a full packed byte (4 bases) per step via the LUT
	for (; i + 4 <= en; i += 4) {
		uint64_t w = l2b->pac[i >> 5];
		unsigned byte = (unsigned)(w >> (((i >> 2) & 7) << 3)) & 0xff;
		memcpy(&seq[i - st], &g_l2b_unpack4[byte], 4);
	}
	// tail: remaining scalar bases
	for (; i < en; ++i)
		seq[i - st] = l2b_get0(l2b, i);
	// retrieve ambiguous bases
	aid = l2b_getambi(l2b, tid, st - ctg->off, en - ctg->off, &n_ambi);
	if (aid >= 0) {
		for (i = 0; i < n_ambi; ++i) {
			const l2b_intv_t *iv = &l2b->ambi[aid + i];
			int64_t s = iv->st, e = iv->en;
			if (s < st) s = st;
			if (e > en) e = en;
			if (s < e) memset(&seq[s - st], 4, e - s);
		}
	}
	return en - st;
}

void l2b_meth_convert(l2b_meth_t mt, int64_t len, uint8_t *seq)
{
	int64_t i;
	if (mt == L2B_METH_C2T) {
		for (i = 0; i < len; ++i)
			if (seq[i] == 1) seq[i] = 3;
	} else if (mt == L2B_METH_G2A) {
		for (i = 0; i < len; ++i)
			if (seq[i] == 2) seq[i] = 0;
	}
}

int64_t l2b_getambi(const l2b_t *l2b, int64_t tid, int64_t st, int64_t en, int32_t *n_ambi)
{
	int64_t g_beg, g_end, lo, hi, mid, i_st, i_en;
	*n_ambi = 0;
	if (tid < 0 || tid >= l2b->n_ctg) return -1;
	if (st < 0) st = 0;
	if (en > l2b->ctg[tid].len) en = l2b->ctg[tid].len;
	if (st >= en) return -1;
	g_beg = l2b->ctg[tid].off + st;
	g_end = l2b->ctg[tid].off + en;

	lo = 0, hi = l2b->n_ambi;
	while (lo < hi) {
		mid = (lo + hi) / 2;
		if (l2b->ambi[mid].en > g_beg) hi = mid;
		else lo = mid + 1;
	}
	i_st = lo;

	lo = i_st, hi = l2b->n_ambi;
	while (lo < hi) {
		mid = (lo + hi) / 2;
		if (l2b->ambi[mid].st >= g_end) hi = mid;
		else lo = mid + 1;
	}
	i_en = lo;

	*n_ambi = i_en - i_st;
	if (*n_ambi == 0) return -1;
	return i_st;
}

static void l2b_format_seq(uint64_t len, char *seq, uint64_t *rng)
{
	uint64_t i;
	for (i = 0; i < len; ++i) {
		int b = (uint8_t)seq[i];
		uint8_t c = kom_nt4_table[b];
		if (c > 4) c = 4;
		if (c == 4) c |= kom_splitmix64(rng) & 3;
		if (b < 'A' || b > 'Z') c |= 1<<3;
		seq[i] = c;
	}
}

static void l2b_add_seq(l2b_t *l2b, uint64_t len, const char *seq, const char *name, const char *comm, uint64_t *rng)
{
	uint64_t i, ambi_len, mask_len, off, m_pac_old;
	l2b_ctg_t *ctg;

	off = l2b->tot_len;
	kom_grow(l2b_ctg_t, l2b->ctg, l2b->n_ctg, l2b->m_ctg);
	ctg = &l2b->ctg[l2b->n_ctg++];
	ctg->name = kom_strdup(name);
	ctg->comm = comm? kom_strdup(comm) : 0;
	ctg->len = len;
	ctg->off = l2b->tot_len;
	l2b->tot_len += len;

	m_pac_old = l2b->m_pac;
	l2b->n_pac = (l2b->tot_len + 31) / 32;
	kom_grow(uint64_t, l2b->pac, l2b->n_pac, l2b->m_pac);
	if (m_pac_old < l2b->m_pac) // zero out newly allocated part
		memset(&l2b->pac[m_pac_old], 0, (l2b->m_pac - m_pac_old) * 8);

	for (i = 0, ambi_len = mask_len = 0; i < len; ++i) {
		uint64_t c = (uint8_t)seq[i], x = off + i;
		if (c & 1<<3) { // soft-masked base
			++mask_len;
		} else if (mask_len > 0) {
			kom_grow(l2b_intv_t, l2b->mask, l2b->n_mask, l2b->m_mask);
			l2b->mask[l2b->n_mask].st = x - mask_len;
			l2b->mask[l2b->n_mask].en = x;
			l2b->n_mask++;
			mask_len = 0;
		}
		if (c & 1<<2) { // ambiguous base
			++ambi_len;
		} else if (ambi_len > 0) {
			kom_grow(l2b_intv_t, l2b->ambi, l2b->n_ambi, l2b->m_ambi);
			l2b->ambi[l2b->n_ambi].st = x - ambi_len;
			l2b->ambi[l2b->n_ambi].en = x;
			l2b->n_ambi++;
			ambi_len = 0;
		}
		l2b->pac[x>>5] |= (c&3) << (x&0x1f)*2;
	}
}

static void l2b_collate_str(l2b_t *l2b)
{
	uint64_t i, tot_name = 0, tot_comm = 0;
	char *p_name, *p_comm;
	if (l2b->cat_name || l2b->cat_comm) return;
	for (i = 0; i < l2b->n_ctg; ++i) {
		l2b_ctg_t *ctg = &l2b->ctg[i];
		tot_name += strlen(ctg->name) + 1;
		tot_comm += ctg->comm? strlen(ctg->comm) + 1 : 1;
	}
	p_name = l2b->cat_name = kom_calloc(char, tot_name);
	p_comm = l2b->cat_comm = kom_calloc(char, tot_comm);
	for (i = 0; i < l2b->n_ctg; ++i) {
		l2b_ctg_t *ctg = &l2b->ctg[i];
		uint64_t len;
		len = strlen(ctg->name);
		memcpy(p_name, ctg->name, len + 1);
		free(ctg->name);
		ctg->name = p_name;
		p_name += len + 1;
		if (ctg->comm) {
			len = strlen(ctg->comm);
			memcpy(p_comm, ctg->comm, len + 1);
			free(ctg->comm);
			ctg->comm = p_comm;
		} else len = 0;
		p_comm += len + 1;
	}
}

l2b_t *l2b_import(const char *fn, uint64_t seed)
{
	gzFile fp;
	kseq_t *ks;
	l2b_t *l2b;
	uint64_t rng = seed;

	fp = fn == 0 || strcmp(fn, "-") == 0? gzdopen(0, "r") : gzopen(fn, "r");
	if (fp == 0) return 0;
	ks = kseq_init(fp);
	l2b = kom_calloc(l2b_t, 1);
	while (kseq_read(ks) >= 0) {
		l2b_format_seq(ks->seq.l, ks->seq.s, &rng);
		l2b_add_seq(l2b, ks->seq.l, ks->seq.s, ks->name.s, ks->comment.l? ks->comment.s : 0, &rng);
	}
	kseq_destroy(ks);
	gzclose(fp);
	l2b_collate_str(l2b);
	return l2b;
}

void l2b_destroy(l2b_t *l2b)
{
	if (l2b->mmap) { // ambi/mask/pac/cat_name/cat_comm point into the mapped file
		free(l2b->ctg); // ctg[] is always heap-allocated
		kom_munmap(l2b->mmap, l2b->mmap_len);
	} else {
		free(l2b->cat_name); free(l2b->cat_comm);
		free(l2b->pac); free(l2b->ambi); free(l2b->mask); free(l2b->ctg);
	}
	free(l2b);
}

int l2b_save(const char *fn, const l2b_t *l2b)
{
	FILE *fp;
	uint64_t i, len_name = 0, len_comm = 0;
	uint32_t dummy = 0;
	fp = fn == 0 || strcmp(fn, "-") == 0? stdout : fopen(fn, "wb");
	if (fp == 0) return -1;
	for (i = 0; i < l2b->n_ctg; ++i) {
		const l2b_ctg_t *ctg = &l2b->ctg[i];
		len_name += strlen(ctg->name) + 1;
		len_comm += ctg->comm? strlen(ctg->comm) + 1 : 1;
	}
	fwrite(L2B_MAGIC, 1, 4, fp);
	fwrite(&dummy, 4, 1, fp);
	fwrite(&l2b->n_ctg, 8, 1, fp);
	fwrite(&l2b->tot_len, 8, 1, fp);
	fwrite(&l2b->n_ambi, 8, 1, fp);
	fwrite(&l2b->n_mask, 8, 1, fp);
	fwrite(&len_name, 8, 1, fp);
	fwrite(&len_comm, 8, 1, fp);
	fwrite(&l2b->n_pac, 8, 1, fp);
	for (i = 0; i < l2b->n_ctg; ++i)
		fwrite(&l2b->ctg[i].len, 8, 1, fp);
	fwrite(l2b->ambi, 16, l2b->n_ambi, fp);
	fwrite(l2b->mask, 16, l2b->n_mask, fp);
	fwrite(l2b->pac, 8, l2b->n_pac, fp);
	fwrite(l2b->cat_name, 1, len_name, fp); // put strings at the end to make sure uint64_t are all aligned
	fwrite(l2b->cat_comm, 1, len_comm, fp);
	fclose(fp);
	return 0;
}

l2b_t *l2b_load(const char *fn)
{
	FILE *fp;
	char magic[4], *p_name, *p_comm;
	uint32_t dummy;
	uint64_t off, i, len_name, len_comm;
	l2b_t *l2b;
	fp = fn == 0 || strcmp(fn, "-") == 0? stdin : fopen(fn, "rb");
	if (fp == 0) return 0;
	fread(magic, 1, 4, fp);
	if (strncmp(magic, L2B_MAGIC, 4) != 0) {
		if (fp != stdin) fclose(fp);
		return 0;
	}
	l2b = kom_calloc(l2b_t, 1);
	fread(&dummy, 4, 1, fp);
	fread(&l2b->n_ctg, 8, 1, fp);
	fread(&l2b->tot_len, 8, 1, fp);
	fread(&l2b->n_ambi, 8, 1, fp);
	fread(&l2b->n_mask, 8, 1, fp);
	fread(&len_name, 8, 1, fp);
	fread(&len_comm, 8, 1, fp);
	fread(&l2b->n_pac, 8, 1, fp);
	l2b->ctg = kom_calloc(l2b_ctg_t, l2b->n_ctg);
	for (i = 0, off = 0; i < l2b->n_ctg; ++i) { // read contig lengths
		fread(&l2b->ctg[i].len, 8, 1, fp);
		l2b->ctg[i].off = off;
		off += l2b->ctg[i].len;
	}
	if (off != l2b->tot_len) goto load_failure;
	l2b->ambi = kom_malloc(l2b_intv_t, l2b->n_ambi);
	l2b->mask = kom_malloc(l2b_intv_t, l2b->n_mask);
	l2b->pac = kom_malloc(uint64_t, l2b->n_pac);
	l2b->cat_name = kom_malloc(char, len_name);
	l2b->cat_comm = kom_malloc(char, len_comm);
	fread(l2b->ambi, 16, l2b->n_ambi, fp);
	fread(l2b->mask, 16, l2b->n_mask, fp);
	fread(l2b->pac, 8, l2b->n_pac, fp);
	fread(l2b->cat_name, 1, len_name, fp);
	fread(l2b->cat_comm, 1, len_comm, fp);
	p_name = l2b->cat_name, p_comm = l2b->cat_comm;
	for (i = 0; i < l2b->n_ctg; ++i) { // synchronize contig names and comments
		l2b_ctg_t *ctg = &l2b->ctg[i];
		ctg->name = p_name;
		p_name += strlen(p_name) + 1;
		ctg->comm = *p_comm? p_comm : 0;
		p_comm += *p_comm? strlen(p_comm) + 1 : 1;
	}
	if (p_name - l2b->cat_name != len_name || p_comm - l2b->cat_comm != len_comm) goto load_failure;
	if (fp != stdin) fclose(fp);
	return l2b;
load_failure:
	if (fp != stdin) fclose(fp);
	l2b_destroy(l2b);
	return 0;
}

l2b_t *l2b_load_mmap(const char *fn, int preload)
{
	uint8_t *base;
	const uint64_t *hdr, *lens;
	uint64_t off, i, len_name, len_comm, foff;
	size_t map_len;
	char *p_name, *p_comm;
	l2b_t *l2b;

	base = (uint8_t*)kom_mmap_file(fn, &map_len, preload);
	if (base == 0) return 0;
	if (map_len < 64 || strncmp((const char*)base, L2B_MAGIC, 4) != 0) { kom_munmap(base, map_len); return 0; }

	l2b = kom_calloc(l2b_t, 1);
	l2b->mmap = base;
	l2b->mmap_len = map_len;
	hdr = (const uint64_t*)base; // hdr[0] is magic+dummy; fields start at hdr[1]
	l2b->n_ctg   = hdr[1];
	l2b->tot_len = hdr[2];
	l2b->n_ambi  = hdr[3];
	l2b->n_mask  = hdr[4];
	len_name     = hdr[5];
	len_comm     = hdr[6];
	l2b->n_pac   = hdr[7];

	l2b->ctg = kom_calloc(l2b_ctg_t, l2b->n_ctg); // ctg[] is not stored in the file
	lens = &hdr[8]; // contig lengths follow the 64-byte header
	for (i = 0, off = 0; i < l2b->n_ctg; ++i) {
		l2b->ctg[i].len = lens[i];
		l2b->ctg[i].off = off;
		off += lens[i];
	}
	if (off != l2b->tot_len) goto mmap_failure;

	foff = 64 + l2b->n_ctg * 8; // point large arrays into the mapped file
	l2b->ambi = (l2b_intv_t*)(base + foff); foff += l2b->n_ambi * 16;
	l2b->mask = (l2b_intv_t*)(base + foff); foff += l2b->n_mask * 16;
	l2b->pac  = (uint64_t*)  (base + foff); foff += l2b->n_pac  * 8;
	l2b->cat_name = (char*)(base + foff);   foff += len_name;
	l2b->cat_comm = (char*)(base + foff);   foff += len_comm;
	if (foff != map_len) goto mmap_failure;

	p_name = l2b->cat_name, p_comm = l2b->cat_comm;
	for (i = 0; i < l2b->n_ctg; ++i) { // synchronize contig names and comments
		l2b_ctg_t *ctg = &l2b->ctg[i];
		ctg->name = p_name;
		p_name += strlen(p_name) + 1;
		ctg->comm = *p_comm? p_comm : 0;
		p_comm += *p_comm? strlen(p_comm) + 1 : 1;
	}
	if (p_name - l2b->cat_name != len_name || p_comm - l2b->cat_comm != len_comm) goto mmap_failure;
	return l2b;
mmap_failure:
	l2b_destroy(l2b);
	return 0;
}

/******************************
 * Save .paf files for bwtgen *
 ******************************/

int l2b_save_pac(const char *fn, const l2b_t *l2b, int both_strand)
{
	FILE *fp;
	uint64_t n_pac, x;
	int64_t i;
	uint8_t *pac, ct;

	fp = fn == 0 || strcmp(fn, "-") == 0? stdout : fopen(fn, "wb");
	if (fp == 0) return -1;

	// fill pac[]
	n_pac = ((both_strand? l2b->tot_len * 2 : l2b->tot_len) + 3) / 4;
	pac = kom_calloc(uint8_t, n_pac);
	for (i = 0, x = 0; i < l2b->tot_len; ++i, ++x)
		pac[x>>2] |= l2b_get0(l2b, i) << (~x&3) * 2;
	if (both_strand)
		for (i = l2b->tot_len - 1; i >= 0; --i, ++x)
			pac[x>>2] |= (3 - l2b_get0(l2b, i)) << (~x&3) * 2;

	// write pac
	fwrite(pac, 1, (x>>2) + ((x&3) == 0? 0 : 1), fp);
	// the following codes make the pac file size always (x/4+1+1)
	if (x % 4 == 0) {
		ct = 0;
		fwrite(&ct, 1, 1, fp);
	}
	ct = x % 4;
	fwrite(&ct, 1, 1, fp);

	fclose(fp);
	free(pac);
	return 0;
}

static inline uint8_t l2b_c2t(uint8_t b) { return b == 1? 3 : b; } // C(1) -> T(3)
static inline uint8_t l2b_g2a(uint8_t b) { return b == 2? 0 : b; } // G(2) -> A(0)

int l2b_save_pac_meth(const char *fn, const l2b_t *l2b, int both_strand)
{
	FILE *fp;
	uint64_t n_pac, x, len;
	int64_t i;
	uint8_t *pac, ct, b;

	fp = fn == 0 || strcmp(fn, "-") == 0? stdout : fopen(fn, "wb");
	if (fp == 0) return -1;

	len = l2b->tot_len * 2; // c2t + g2a
	if (both_strand) len *= 2; // forward + reverse
	n_pac = (len + 3) / 4;
	pac = kom_calloc(uint8_t, n_pac);

	// c2t forward
	for (i = 0, x = 0; i < l2b->tot_len; ++i, ++x) {
		b = l2b_c2t(l2b_get0(l2b, i));
		pac[x>>2] |= b << ((~x&3) * 2);
	}
	// g2a forward
	for (i = 0; i < l2b->tot_len; ++i, ++x) {
		b = l2b_g2a(l2b_get0(l2b, i));
		pac[x>>2] |= b << ((~x&3) * 2);
	}
	if (both_strand) {
		// g2a reverse (reverse complement of g2a converted sequence)
		for (i = l2b->tot_len - 1; i >= 0; --i, ++x) {
			b = 3 - l2b_g2a(l2b_get0(l2b, i));
			pac[x>>2] |= b << ((~x&3) * 2);
		}
		// c2t reverse (reverse complement of c2t converted sequence)
		for (i = l2b->tot_len - 1; i >= 0; --i, ++x) {
			b = 3 - l2b_c2t(l2b_get0(l2b, i));
			pac[x>>2] |= b << ((~x&3) * 2);
		}
	}

	fwrite(pac, 1, (x>>2) + ((x&3) == 0? 0 : 1), fp);
	if (x % 4 == 0) {
		ct = 0;
		fwrite(&ct, 1, 1, fp);
	}
	ct = x % 4;
	fwrite(&ct, 1, 1, fp);

	fclose(fp);
	free(pac);
	return 0;
}
