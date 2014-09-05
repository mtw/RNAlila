/*
  moves.c : move-set related routines for RNAlila
  Last changed Time-stamp: <2014-09-05 15:12:28 mtw>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <lila.h>
#include <moves.h>

#define MINGAP 3
  
static int lila_construct_moves(const char*, const short*, int , move_str **);
static int lila_RNAlexicographicalOrder(const void *, const void *);
inline int try_insert_seq2(const char*, int, int);
inline int compat(const char, const char);
void lila_dump_pt(const short*);
void lila_structure2hash(gpointer,gpointer);
  
typedef struct _nb {
  char *struc;
  int ediff;
  short n;
} nbT;


/*
  Generate all neighbors of a secondary structure (provided via a pair
  table). Returns a pointer to an energy-lexicographically sorted
  GQueue of LilaDBE-type elements .
*/
GQueue *
lila_generate_neighbors_pt(const char *seq,
			   short int *pt)
{
  int i,e,count;
  move_str r,*mvs=NULL;

  GQueue *N = g_queue_new();
  e = vrna_eval_structure_pt(seq,pt,P);
  count = lila_construct_moves(seq,pt,1,&mvs);

  for(i=0;i<count;i++) {
    LilaDBE *nb = NULL;
    nb = (LilaDBE*)calloc(1,sizeof(LilaDBE));
    move_str m = mvs[i];
    int emove =  vrna_eval_move_pt(pt,s0,s1,m.left,m.right,P);
    
    short int *ptbak = vrna_pt_copy(pt);
    lila_apply_move_pt(ptbak,m);
    char *struc =  lila_db_from_pt(ptbak);
    nb->structure = strdup(struc);
    nb->energy = (float)(e + emove)/100;
    g_queue_push_tail(N,nb);
    free(ptbak);
    free(struc);
  }
  free(mvs);
  return N;
  // call lila_get_cc_pt for each neighbor
}

/*
  get random move operation on a pair table
  returns move operations to be applied to pt in order to perform the move
*/
move_str
lila_random_move_pt(const char *seq,
		    short int *pt)
{
  move_str r,*mvs=NULL;
  int count;
  
  count = lila_construct_moves((const char *)seq,pt,1,&mvs);
  r.left  = mvs[0].left;
  r.right = mvs[0].right;
  free(mvs);
  return r;
}

/*
  Get gradient (steepest descent) move operation to lexicographically
  smallest neighbor on a pair table. Returns move operations to be
  applied to pt in order to perform the move
*/
move_str
lila_gradient_move_pt(const char *seq, short int *pt)
{
  move_str r,*mvs=NULL;
  int emove,count,i,k=0,mindiff=100000;
  int max_neighbours;
  nbT *neighbours = NULL;
  
  count = lila_construct_moves((const char *)seq,pt,1,&mvs);
  neighbours = (nbT*)calloc(count+1,sizeof(nbT));
  
  for(i=0;i<count;i++) {
    neighbours[i].ediff = vrna_eval_move_pt(pt,s0,s1,mvs[i].left,mvs[i].right,P);
    /* seems we're storing the same pt here sinde the move is not applied */
    neighbours[i].struc = vrna_pt_to_db(pt); 
    neighbours[i].n = i;
  }
  qsort(neighbours, count, sizeof(nbT), lila_RNAlexicographicalOrder);
  
  r.left  = mvs[neighbours[0].n].left;
  r.right = mvs[neighbours[0].n].right;

  for(i=0;i<count;i++) 
    free(neighbours[i].struc);
  free(neighbours);
  free(mvs);
  return r;
}

/*
  get one adaptive move operation on a pair table
  returns move operations to be applied to pt in order to perform the move
*/
move_str
lila_adaptive_move_pt(const char *seq,
		      short int *pt)
{
  move_str r,*mvs=NULL;
  int i,count,have_aw_nb=0;

  count = lila_construct_moves((const char *)seq,pt,1,&mvs);
  for(i=0;i<count;i++) {
    if(vrna_eval_move_pt(pt,s0,s1,mvs[i].left,mvs[i].right,P) < 0){
      /* here we should also accept degenerate moves IFF the target
	 structure is lexicographically smaller that the source
	 structure - check that !*/
	 
      r.left = mvs[i].left;
      r.right = mvs[i].right;
      have_aw_nb=1;
      /* fprintf(stderr,"applying adaptive move operation [l:%i r:%i]\n",
	 r.left,r.right); */
      break;
    }
    /* check case when we get a degenerate neighbor at the same energy
       level: call lila_is_minimum() and if this returns 0
       ensure that it is lexicographically smaller AND ensure
       that no further back-and-forth moves are possible in follow-up
       moves. this can eventually turn out to be a bit tricky */
  }
  if(have_aw_nb == 0){
    r.left = 0;
    r.right = 0;
    /* fprintf (stderr,"no adaptive neighbor found, hence moves are
       [l:%i r: %i]\n", r.left,r.right);*/
  }
  free(mvs);
  return r;
}

/*
  get all adaptive move operations on a pair table. Returns an array
  of individual move operations to be applied to pt or NULL if no
  moves are possible. The 'want_degenerate' argument indicates whether
  degenerate neighbors are considered adaptive.
*/
move_str* 
lila_all_adaptive_moves_pt(const char *seq,
			   short int *pt,
			   int *ct,
			   int want_degenerate)
{
  move_str r,*mvs=NULL, *allmvs=NULL;
  int i,count,e,j=0;
  
  count = lila_construct_moves((const char *)seq,pt,1,&mvs);
  allmvs = (move_str*)calloc(count,sizeof(move_str));

#ifdef PARANOID
  { /* begin paranoid energy evaluation  REMOVE ME */
     e = vrna_eval_structure_pt(seq,pt,P);
  } /* end paranoid energy evaluation */
#endif
  
  for(i=0;i<count;i++){
    int ediff = vrna_eval_move_pt(pt,s0,s1,mvs[i].left,mvs[i].right,P);
    if (want_degenerate == 1){
      if(ediff <= 0){
	allmvs[j] = mvs[i];
	j++;
      }
    }
    else{
      if(ediff < 0){
	allmvs[j] = mvs[i];
	j++;
	
#ifdef PARANOID      
	{ /* begin paranoid energy evaluation REMOVE ME */
	  int emove,enew,e2;
	  short int *ptbak=NULL;
	  move_str m;
	  m.left  = mvs[i].left;
	  m.right = mvs[i].right;
	  /* copy pair table, operate on this copy */
	  ptbak = vrna_pt_copy(pt);
	  /* compute energy difference for this move */
	  emove = vrna_eval_move_pt(ptbak,s0,s1,m.left,m.right,P);
	  /* evaluate energy of the new structure */
	  enew = e + emove;
	  /* do the move */
	  lila_apply_move_pt(ptbak,m);
	  /* eval energy of the new structure */
	  e2 =  vrna_eval_structure_pt(seq,ptbak,P);
	  
	  if (e2 != enew){
	    fprintf(stderr, "energy evaluation against vrna_eval_structure_pt() mismatch... HAVE %6.2f != %6.2f (SHOULD BE)\n",(float)enew/100, (float)e2/100);
	    
	    fprintf(stderr,"INPUT pt:\n");
	    lila_dump_pt(pt);
	    fprintf(stderr, "AW NEIGHBOUR pt:\n");
	    lila_dump_pt(ptbak);
	  }
	  free(ptbak);
	} /* end paranoid energy evaluation */
#endif      
	
      } /* end if */
    } /* end else */
  } /* end for */
  *ct = j;
  free(mvs);
  if(j>0)
    return allmvs;
  else{
    free(allmvs);
    return NULL;
  }
}


/*
  apply one move operation on a pair table
*/
void
lila_apply_move_pt(short int *pt,
		   move_str m)
{
  if(m.left < 0){
    pt[(int)(fabs(m.left))] = 0;
    pt[(int)(fabs(m.right))] = 0;
  }
  else {
    pt[m.left] = m.right;
    pt[m.right] = m.left;
  }
  //print_str(stdout,pt);printf("\n");
}

static int
lila_construct_moves(const char *seq,
		     const short *structure,
		     int permute,
		     move_str **array)
{
  /* generate all possible moves (less than n^2)*/
  int i;
  int size = 4;
  int count = 0;
  move_str *res = (move_str*) malloc(sizeof(move_str)*(size+1));
  
  for (i=1; i<=structure[0]; i++) {
    if (structure[i]!=0) {
      if (structure[i]<i) continue;
      count ++;
      // need to reallocate the array?
      if (count>size) {
      	size *= 2;
      	res = realloc(res, sizeof(move_str)*(size));
      }
      res[count-1].left = -i;
      res[count-1].right = -structure[i];
      //fprintf(stderr, "add  d(%d, %d)\n", i, structure[i]);
    } else {
      int j;
      for (j=i+1; j<=structure[0]; j++) {
        /* fprintf(stderr, "check (%d, %d)\n", i, j); */
        if (structure[j]==0) {
          if (try_insert_seq2(seq,i,j)) {
            count ++;
	    // need to reallocate the array?
	    if (count>size) {
	      size *= 2;
	      res = realloc(res, sizeof(move_str)*(size));
	    }
	    res[count-1].left = i;
	    res[count-1].right = j;
            //fprintf(stderr, "add  i(%d, %d)\n", i, j);
            continue;
          }
        } else if (structure[j]>j) { /*  '(' */
          j = structure[j];
        } else break;
      }
    }
  }
  
  res = realloc(res, sizeof(move_str)*(count));
  
  /* permute them */
  if (permute) {
    for (i=0; i<count; i++) {
      int rnd = rand(); 
      rnd = rnd % (count-i) + i;
      move_str mv;
      mv = res[i];
      res[i] = res[rnd];
      res[rnd] = mv;
    }
  }
  *array = res;
  return count;
}

/*  try insert base pair (i,j) */
inline int
try_insert_seq2(const char *seq,
	       int i,
	       int j)
{
  if (i<=0 || j<=0) return 0;
  return (j-i>MINGAP && compat(seq[i-1], seq[j-1]));
}

/* compatible base pair?*/
inline int
compat(const char a,
       const char b)
{
  if (a=='A' && b=='U') return 1;
  if (a=='C' && b=='G') return 1;
  if (a=='G' && b=='U') return 1;
  if (a=='U' && b=='A') return 1;
  if (a=='G' && b=='C') return 1;
  if (a=='U' && b=='G') return 1;
  /* and with T's*/
  if (a=='A' && b=='T') return 1;
  if (a=='T' && b=='A') return 1;
  if (a=='G' && b=='T') return 1;
  if (a=='T' && b=='G') return 1;
  return 0;
}


void
lila_dump_pt(const short *pairtable)
{
  int i;
  printf("> ");
  for (i=0;i<(*pairtable)+1;i++){
    printf("%i ",*(pairtable+i));
  }
  printf("\n");
}

static int 
lila_RNAlexicographicalOrder(const void *a, const void *b)
{
  const nbT *ma = a;
  const nbT *mb = b;
  
  int comp = ma->ediff - mb->ediff;

  if (comp < 0)
    return -1;
  
  if (comp > 0)
    return 1;
  
  comp = strcmp(ma->struc, mb->struc);
  return comp;
}

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
    emove = vrna_eval_move_pt(pt,s0,s1,mvs[i].left,mvs[i].right,P);
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
  convert pair table to dot bracket notation
*/
char *
lila_db_from_pt(short int *pt)
{
  int i;
  char *db=NULL;
  //fprintf(stderr,"[[lila_db_from_pt]]:\n");
  //lila_dump_pt(pt);
  if(pt){
    db = (char*)calloc(pt[0]+1,sizeof(char));
    for(i=1;i<=pt[0];i++){
      if(*(pt+i)==0)
	db[i-1] = '.';
      else if (*(pt+i)<i)
	db[i-1] = ')';
      else
	db[i-1] = '(';
    }
    db[i-1] = '\0';
  }
  //fprintf(stderr,"%s\n",db);
  return db;
}

/*
  Get the minimal (lexicographically smallest) representative of a
  degenerate connected component, starting from a pair table of one
  representative.
*/

/*
  Construct the connected component that contains a structure x, whose
  pt is provided as argument. Populates a list of type Lila2seT
  containing all dot bracket structures of the cc and their respective
  energies. The return value indicates whether the component is
  minimal, ie. whether it has lower energy neighbors (0 or 1).
 */
GList *
lila_get_cc_pt(const char *seq,
	       short int *pt,
	       int *ismin)
{
  int min,e;
  short int *minpt = NULL;
  char *v = NULL;
  GQueue *TODO = g_queue_new(); /* the TODO list */
  GQueue *SEEN = g_queue_new(); /* list of structures already processed */
  GList *cc = NULL;
  Lila2seT *element =NULL;
  
  min = 1; /* initially mark this component minimal */
  v = lila_db_from_pt(pt);
  
  /* add first structure to TODO & SEEN queue */
  g_queue_push_tail(TODO,v);
  g_queue_push_tail(SEEN,v);

  /* add first element to cc list */
  element = (Lila2seT*)calloc(1,sizeof(Lila2seT));
  element->structure = strdup(v);
  element->energy = vrna_eval_structure(seq,v,P);
  /* fprintf(stderr, "\n%s (%6.4f) added to cc list\n",v,element->energy);*/
  cc = g_list_append(cc,element);
  
  /* fprintf(stderr, ">> PUSH %s\n",v); */
 
  while(g_queue_is_empty(TODO) != TRUE){
    move_str *mvs = NULL;
    int i,count,emove;
    short int *p;

    v = (char*)g_queue_pop_head(TODO);
    p = vrna_pt_get(v);
    e = vrna_eval_structure_pt(seq,p,P);
    /* fprintf(stderr,"%s (%6.2f) POP queuelen=%i\n",
       v,(float)e/100, g_queue_get_length(TODO)); */
    
    /* generate neighbors of v  */
    count = lila_construct_moves((const char *)seq,p,1,&mvs);
    /* fprintf(stderr,">>> %i moves possible\n", count); */

    for(i=0;i<count;i++) {
      emove = vrna_eval_move_pt(p,s0,s1,mvs[i].left,mvs[i].right,P);
      if(emove == 0){ 	/* degenerate neighbor */
	/* fprintf(stderr," move #%i l:%3i|r:%3i  ",i,mvs[i].left,mvs[i].right);
	fprintf(stderr," d(%6.2f) ",(float)emove/100);
	fprintf(stderr, "degenerate  \n");*/
	
	short int *ptbak = vrna_pt_copy(p);
	lila_apply_move_pt(ptbak,mvs[i]);
	char *w = lila_db_from_pt(ptbak);
	
	/* insert into queue if not yet present */
	if (g_queue_find_custom(SEEN,w,(GCompareFunc)lila_cmp_db) == NULL){
	 
	  g_queue_push_tail(TODO,w);
	  g_queue_push_tail(SEEN,w);
	  /* fprintf(stderr, "%s d(%6.2f) PUSH", w, (float)emove/100); */
	  
	  /* add structure and energy to connected component list */
	  Lila2seT *e = NULL;
	  e = (Lila2seT*)calloc(1,sizeof(Lila2seT));
	  e->structure = strdup(w);
	  e->energy = vrna_eval_structure(seq,w,P);
	  cc = g_list_append(cc,e);
	  /* fprintf(stderr, "%s (%6.4f) added to cc list\n",w,e->energy);*/
	  /* fprintf(stderr,"\n"); */
	}
	free(ptbak);
      } /* end if */
      else if (emove <0){ /* shoulder */
	min = 0; /* mark component as non-minimal */
	/* fprintf(stderr," move #%i l:%3i|r:%3i  ",i,mvs[i].left,mvs[i].right);
	fprintf(stderr," d(%6.2f) ",(float)emove/100);
	fprintf(stderr, "we have a shoulder\n"); */
      }
      else {
	//fprintf(stderr,"\n");
      }
    } /* end for */
    free(mvs);
    //free(v);
  } /* end while */
  g_queue_free(TODO);
  g_queue_free(SEEN);

  *ismin = min;
  return cc;
}
/*
  Insert all elements of a GList to a GHash
*/
void
lila_cc2hash(GHashTable *S,
	     GList *cc)
{
  fprintf(stderr,"[[lila_cc2hash]]::BEFORE %i\n", (int)g_hash_table_size(S));
  g_list_foreach(cc,(GFunc)lila_structure2hash,S);
  fprintf(stderr,"[[lila_cc2hash]]::AFTER %i\n", (int)g_hash_table_size(S));
}

/*
  Add a secondary structure from Lila2seT into a hash
*/
void
lila_structure2hash(gpointer data,
		    gpointer user_data)
{
  Lila2seT *foo = (Lila2seT *)data;
  GHashTable *S = (GHashTable *)user_data;
  char *form = strdup(foo->structure);
  if (g_hash_table_contains(S,form))
    fprintf(stderr, "%s ALREADY in hash ",form);
  else{
    g_hash_table_add(S, form);
    fprintf(stderr, "%s INSERTED to hash ",form);
  }
  fprintf(stderr, " [hashsize(S): %i]\n",(int)g_hash_table_size(S));
  
}

/*
  Return the lexicographically lowest member of a connected component
*/
GList *
lila_lexmin_cc(GList *c)
{
  //qsort(cc,cc_entries,sizeof(Lila2seT),lila_cmp_sse_lex);
  fprintf(stderr,"[[lila_lexmin_cc]]::GList length is %i\n", (int)g_list_length(c));
  c = g_list_sort(c,(GCompareFunc)lila_cmp_sse_lex);
  GList *first = g_list_first(c);
  return first;
}

/*
  Dump connected component array to stderr
  SEE lila_output_dbe_queue in lila.c
*/
void
lila_dump_cc(GList *c)
{
  /* fprintf(stderr,"[[lila_dump_cc]]::START\n"); */
  g_list_foreach(c,(GFunc)lila_print_dbe,NULL);
  /* fprintf(stderr,"[[lila_dump_cc]]::END\n");	  */ 
}

/*
  Cleanup connected component array
*/
void
lila_cleanup_cc(GList *c)
{
  g_list_free_full(c,(GDestroyNotify)lila_free_cc);
}
