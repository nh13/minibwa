#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "minibwa.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

/* Probe for Task 4 (mb_reconcile_alt grouping).  mb_map() runs the full SE
 * pipeline, including the gated mb_reconcile_alt pass, so the per-hit `parent`
 * field already reflects the FINAL liftover grouping: a hit is its own group
 * representative iff parent == id, otherwise it is a subordinate whose parent is
 * the id of its group's representative.  For every hit of every read it prints:
 *
 *   <qname> ctg=<ctg> qs=<qs> qe=<qe> rev=<rev> is_alt=<is_alt> \
 *           id=<id> parent=<parent> subsc=<subsc> dp_max2=<dp_max2|.>
 *
 * so the test can assert directly on the grouping (parent == id => independent
 * representative; parent == <other id> => merged subordinate) without depending
 * on SAM-flag heuristics that are insensitive to the merge decision. */
int main(int argc, char *argv[])
{
	mb_opt_t opt;
	mb_opt_init(&opt);

	if (argc < 3) {
		fprintf(stderr, "Usage: ex-group-check <idxPrefix> <query.fa>\n");
		return 1;
	}

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
		int32_t j, n_hit;
		hit = mb_map(&opt, idx, ks->seq.l, ks->seq.s, 0, &n_hit, 0, ks->name.s);
		for (j = 0; j < n_hit; ++j) {
			mb_hit_t *h = &hit[j];
			printf("%s\tctg=%s\tqs=%d\tqe=%d\trev=%d\tis_alt=%d\tid=%d\tparent=%d\tsubsc=%d\tdp_max2=",
				ks->name.s, mb_idx_ctg_name(idx, h->tid),
				h->qs, h->qe, h->rev, h->is_alt, h->id, h->parent, h->subsc);
			if (h->p) printf("%d\n", h->p->dp_max2);
			else      printf(".\n");
			free(h->p);
		}
		free(hit);
	}
	mb_idx_destroy(idx);
	kseq_destroy(ks);
	gzclose(f);
	return 0;
}
