#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include "minibwa.h"
#include "mbpriv.h"   /* mb_idx_s -> l2b for mb_hit_place */
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

/* Probe for Task 2 (mb_hit_place).  For every hit of every read it prints:
 *   <qname> <ctg> ts=<ts> te=<te> rev=<rev> is_alt=<is_alt> cg=<cigar|.> ||
 *     place: pri=<pri_ctg> lst=<lifted_st> len=<lifted_en> prev=<rev> liftable=<liftable>
 * so the test script can assert on the placement (pri contig name, lifted_st,
 * folded rev, liftable) without depending on internal tid numbering.
 *
 * Optional 3rd arg <hitspec> = "all" (default) or "coarse": when "coarse" the
 * probe NULLs h->p before calling mb_hit_place to exercise the pre-DP [ts,te]
 * path against the very same hit (so coarse vs exact can be compared directly). */
int main(int argc, char *argv[])
{
	mb_opt_t opt;
	int coarse = 0;
	mb_opt_init(&opt);

	if (argc < 3) {
		fprintf(stderr, "Usage: ex-place-check <idxPrefix> <query.fa> [all|coarse]\n");
		return 1;
	}
	if (argc >= 4 && strcmp(argv[3], "coarse") == 0) coarse = 1;

	gzFile f = gzopen(argv[2], "r");
	assert(f);
	kseq_t *ks = kseq_init(f);

	mb_idx_t *idx = mb_idx_load(argv[1], 0);
	assert(idx);
	/* mb_idx_load auto-detects <prefix>.alt; reload explicitly to be safe. */
	{
		char buf[1024];
		snprintf(buf, sizeof buf, "%s.alt", argv[1]);
		mb_idx_set_alt(idx, buf);
	}

	while (kseq_read(ks) >= 0) {
		mb_hit_t *hit;
		int32_t i, j, n_hit;
		hit = mb_map(&opt, idx, ks->seq.l, ks->seq.s, 0, &n_hit, 0, ks->name.s);
		for (j = 0; j < n_hit; ++j) {
			mb_hit_t *h = &hit[j];
			mb_extra_t *saved = h->p;
			mb_place_t pl;
			const char *pri_name;

			printf("%s\t%s\tts=%ld\tte=%ld\trev=%d\tis_alt=%d\tcg=",
				ks->name.s, mb_idx_ctg_name(idx, h->tid),
				(long)h->ts, (long)h->te, h->rev, h->is_alt);
			if (h->p && h->p->n_cigar > 0) {
				for (i = 0; i < h->p->n_cigar; ++i)
					printf("%d%c", h->p->cigar[i]>>4, MB_CIGAR_STR[h->p->cigar[i]&0xf]);
			} else putchar('.');

			if (coarse) h->p = 0;           /* force the pre-DP [ts,te] path */
			pl = mb_hit_place(idx->l2b, h);
			/* Restore before printing so ctg name lookup via idx still works if
			 * needed.  Free `saved` directly -- it is the one allocation we own;
			 * h->p after the restore is the same pointer, but freeing via `saved`
			 * makes ownership unambiguous and avoids any confusion with the
			 * subsequent free(hit) that releases only the flat hit array. */
			h->p = saved;
			free(saved);
			saved = h->p = NULL;            /* prevent any accidental re-use */

			pri_name = (pl.pri_tid >= 0) ? mb_idx_ctg_name(idx, pl.pri_tid) : "?";
			printf("\t||\tplace:\tpri=%s\tlst=%ld\tlen=%ld\tprev=%d\tliftable=%d\n",
				pri_name, (long)pl.lifted_st, (long)pl.lifted_en, pl.rev, pl.liftable);
		}
		free(hit);
	}
	mb_idx_destroy(idx);
	kseq_destroy(ks);
	gzclose(f);
	return 0;
}
