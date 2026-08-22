#include <stdio.h>
#include <string.h>
#include "minibwa.h"
#include "l2bit.h"
/* Loads <idx>.l2b, sets .alt, prints lift of (chrP_alt, pos) -> (chrP, pos).
 * For pos<300 expect identity; for pos>=305 expect pos-5 (the 5bp insertion). */
int main(int argc, char **argv){
    if (argc<2){fprintf(stderr,"usage: ex-lift-check <idxprefix>\n");return 2;}
    char buf[1024]; snprintf(buf,sizeof buf,"%s.l2b",argv[1]);
    l2b_t *l2b=l2b_load(buf); if(!l2b){fprintf(stderr,"no l2b\n");return 2;}
    snprintf(buf,sizeof buf,"%s.alt",argv[1]); l2b_set_alt(l2b,buf);
    /* find chrP_alt tid */
    int64_t at=-1,i; for(i=0;i<(int64_t)l2b->n_ctg;i++) if(!strcmp(l2b->ctg[i].name,"chrP_alt")) at=i;
    if(at<0||!l2b->ctg[at].is_alt){fprintf(stderr,"chrP_alt not flagged ALT\n");return 1;}
    /* lift two positions via the block list */
    int64_t pt1, pt2; uint64_t pp1, pp2; uint8_t rev1, rev2;
    int ok1=l2b_lift(l2b,at,100,&pt1,&pp1,&rev1);
    int ok2=l2b_lift(l2b,at,400,&pt2,&pp2,&rev2);
    printf("lift(100)=%d pri=%llu  lift(400)=%d pri=%llu\n",ok1,(unsigned long long)pp1,ok2,(unsigned long long)pp2);
    return 0;
}
