/*
  topology.c : routines for determining topological properties of RNA
  energy landscapes

  Last changed Time-stamp: <2016-08-31 16:20:51 mtw>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <lila.h>
#include <moves.h>

/*
  check whether a given structure (provided as pairtable) is a local
  minimum of the energy landscape. return 1 if it is a minimum, 0 if
  it is not a minimum and -1 if it is a degenerate minimum (i.e. there
  are neighbour structures at the same energy level).
 */
int
lila_is_minimum_or_shoulder_pt(const char *seq, short int *pt)
{
  move_str r,*mvs=NULL;
  int rv,i,count,emove;
  
  count = lila_construct_moves((const char *)seq,pt,1,&mvs);
  for(i=0;i<count;i++) {
    emove = vrna_eval_move_pt(vc,pt,mvs[i].left,mvs[i].right);
    if (emove < 0){ /* lower energy neighbor found */
      free(mvs);
      return 0;
    }
    if (emove == 0){ /* degenerate neighbor found */
      free(mvs);
      return -1;
    }
  }
  free(mvs);
  return 1;
}

/*
  Construct the connected component that contains a structure x, whose
  pt is provided as argument. Populates a list of type Lila2seT
  containing all dot bracket structures of the cc and their respective
  energies. The argument 'ismin' indicates whether the component is
  minimal, ie. whether it has lower energy neighbors (0 or 1).
 */
GList *
lila_get_cc_pt(const char *seq,
	       short int *pt,
	       int *ismin)
{
  int min,e;
  short int *minpt = NULL;
  char *v,*vS,*vT;
  GQueue *TODO = g_queue_new(); /* the TODO list */
  GQueue *SEEN = g_queue_new(); /* list of structures already processed */
  GList *cc = NULL;
  LilaDBE *element =NULL;
  
  min = 1; /* initially mark this component minimal */
  v = lila_db_from_pt(pt);
  vS = strdup(v);
  vT = strdup(v);
  
  /* add first structure to TODO & SEEN queue */
  g_queue_push_tail(TODO,vT);
  g_queue_push_tail(SEEN,vS);

  /* add first element to cc list */
  element = (LilaDBE*)calloc(1,sizeof(LilaDBE));
  element->structure = strdup(v);
  element->energy = vrna_eval_structure(vc,v);
  /* fprintf(stderr, "\n%s (%6.2f) added to cc list I\n",v,element->energy); */
  cc = g_list_append(cc,element);
  /* fprintf(stderr, "%s PUSH\n",v); */
  free(v);
 
  while(g_queue_is_empty(TODO) != TRUE){
    char *v2=NULL;
    move_str *mvs = NULL;
    int i,count,emove;
    short int *p = NULL;

    v2 = (char*)g_queue_pop_head(TODO);
    p = vrna_ptable(v2);
    e = vrna_eval_structure_pt(vc,p);
    /* fprintf(stderr,"%s (%6.2f) POP queuelen=%i\n", */
    /*    v2,(float)e/100, g_queue_get_length(TODO)); */ 
    
    /* generate neighbors of v2  */
    count = lila_construct_moves((const char *)seq,p,1,&mvs);
    /* fprintf(stderr,">>> %i moves possible\n", count); */

    for(i=0;i<count;i++) {
      emove = vrna_eval_move_pt(vc,p,mvs[i].left,mvs[i].right);
      /* fprintf(stderr," move #%i l:%3i|r:%3i  ",i,mvs[i].left,mvs[i].right); */
      /* fprintf(stderr," d(%6.2f) ",(float)emove/100); */
      
      if(emove == 0){ 	/* degenerate neighbor */
	/* fprintf(stderr, "degenerate  \n"); */
	short int *ptbak = vrna_ptable_copy(p);
	lila_apply_move_pt(ptbak,mvs[i]);
	char *w,*wT,*wS;
	w = lila_db_from_pt(ptbak);
	wT = strdup(w);
	wS = strdup(w);
	
	/* insert into queue if not yet present */
	if (g_queue_find_custom(SEEN,w,(GCompareFunc)lila_cmp_db) == NULL){
	  g_queue_push_tail(TODO,wT);
	  g_queue_push_tail(SEEN,wS);
	  /* fprintf(stderr, "%s d(%6.2f) PUSH\n", w, (float)emove/100); */
	  
	  /* add structure and energy to connected component list */
	  LilaDBE *e = NULL;
	  e = (LilaDBE*)calloc(1,sizeof(LilaDBE));
	  e->structure = strdup(w);
	  e->energy = vrna_eval_structure(vc,w);
	  cc = g_list_append(cc,e);
	  /* fprintf(stderr, "%s (%6.2f) added to cc list\n",e->structure,e->energy); */
	  /* fprintf(stderr,"\n");  */
	}
	else{
	  free(wT);
	  free(wS);
	}
	free(w);
	free(ptbak);
      } /* end if */
      else if (emove <0){ /* shoulder */
	min = 0; /* mark component as non-minimal */
	/* fprintf(stderr, "SHOULDER\n");  */
      }
      else {
	/* fprintf(stderr, "\n");  */
      }
    } /* end for */
    free(mvs);
    free(v2);
    free(p);
  } /* end while */

  g_queue_free_full(TODO,(GDestroyNotify)lila_free_string);
  g_queue_free_full(SEEN,(GDestroyNotify)lila_free_string);

  *ismin = min;
  return cc;
}
