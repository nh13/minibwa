#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "minibwa.h"
#include "l2bit.h"
/* Loads <idx>.l2b, sets .alt, then either:
 *   (default, no extra args) prints lift of (chrP_alt, pos) -> (chrP, pos):
 *     for pos<300 expect identity; for pos>=305 expect pos-5 (the 5bp insertion).
 *   (check mode) with trailing "<ctg> <pos> <expected>" triples, asserts each
 *     lift's primary position equals <expected> (or expected<0 means "in a hole /
 *     not liftable"); prints one line per triple and exits nonzero on any
 *     mismatch.  Also verifies the contig's lift[] is sorted by alt_st. */
static int64_t find_tid(const l2b_t *l2b, const char *name){
    int64_t i; for(i=0;i<(int64_t)l2b->n_ctg;i++) if(!strcmp(l2b->ctg[i].name,name)) return i; return -1;
}
int main(int argc, char **argv){
    if (argc<2){fprintf(stderr,"usage: ex-lift-check <idxprefix> [<ctg> <pos> <expected> ...]\n");return 2;}
    char buf[1024]; snprintf(buf,sizeof buf,"%s.l2b",argv[1]);
    l2b_t *l2b=l2b_load(buf); if(!l2b){fprintf(stderr,"no l2b\n");return 2;}
    snprintf(buf,sizeof buf,"%s.alt",argv[1]); l2b_set_alt(l2b,buf);
    if (argc>2){ /* check mode: (ctg,pos,expected) triples */
        int rc=0, a;
        if ((argc-2)%3!=0){fprintf(stderr,"check mode needs triples\n");return 2;}
        for(a=2;a+2<argc;a+=3){
            const char *ctg=argv[a]; uint64_t pos=strtoull(argv[a+1],0,10); long exp=strtol(argv[a+2],0,10);
            int64_t at=find_tid(l2b,ctg);
            if(at<0||!l2b->ctg[at].is_alt){fprintf(stderr,"%s not flagged ALT\n",ctg);return 1;}
            unsigned k; int sorted=1;
            for(k=1;k<l2b->ctg[at].n_lift;k++) if(l2b->ctg[at].lift[k].alt_st < l2b->ctg[at].lift[k-1].alt_st) sorted=0;
            int64_t pt; uint64_t pp=0; uint8_t rev;
            int ok=l2b_lift(l2b,at,pos,&pt,&pp,&rev);
            long got = ok? (long)pp : -1;
            int pass = (got==exp) && sorted;
            printf("lift(%s,%llu)=%ld rev=%d sorted=%d expected=%ld %s\n",
                   ctg,(unsigned long long)pos,got,ok?rev:-1,sorted,exp,pass?"OK":"MISMATCH");
            if(!pass) rc=1;
        }
        l2b_destroy(l2b);
        return rc;
    }
    /* default (legacy) mode: chrP_alt at 100 and 400 */
    int64_t at=find_tid(l2b,"chrP_alt");
    if(at<0||!l2b->ctg[at].is_alt){fprintf(stderr,"chrP_alt not flagged ALT\n");return 1;}
    int64_t pt1, pt2; uint64_t pp1, pp2; uint8_t rev1, rev2;
    int ok1=l2b_lift(l2b,at,100,&pt1,&pp1,&rev1);
    int ok2=l2b_lift(l2b,at,400,&pt2,&pp2,&rev2);
    printf("lift(100)=%d pri=%llu  lift(400)=%d pri=%llu\n",ok1,(unsigned long long)pp1,ok2,(unsigned long long)pp2);
    return 0;
}
