/*
  moves.c : move-set related routines for RNAlila
  Last changed Time-stamp: <2014-08-14 23:19:08 mtw>
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

typedef struct _nb {
  char *struc;
  int ediff;
  short n;
} nbT;

/*
  get random move operation on a pair table
  returns move operations to be applied to pt in order to perform the move
*/
move_str
lila_random_move_pt(const char *seq, short int *pt)
{
  move_str r,*mvs=NULL;
  int count;
  
  count = lila_construct_moves((const char *)seq,pt,1,&mvs);
  /*
    {
    for (int i = 0; i<count; i++) {  
    printf("%d %d\n", mvs[i].left, mvs[i].right);
    }
    }
  */
  
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
  get adaptive move operation on a pair table
  returns move operations to be applied to pt in order to perform the move
*/
move_str
lila_adaptive_move_pt(const char *seq, short int *pt)
{
  move_str r,*mvs=NULL;
  int i,count;

  count = lila_construct_moves((const char *)seq,pt,1,&mvs);
  fprintf(stderr," count is %i\n",count);
  for(i=0;i<count;i++) {
    if(vrna_eval_move_pt(pt,s0,s1,mvs[i].left,mvs[i].right,P) < 0){
      r.left = mvs[i].left;
      r.right = mvs[i].right;
      fprintf(stderr,"applying adaptive move operation [l:%i r:%i]\n",
	      r.left,r.right);
    }
    else {
      r.left = 0;
      r.right = 0;
      fprintf (stderr,"RNAlila::moves.c::lila_adaptive_move_pt no adaptive neighbor found, hence moves are [l:%i r: %i]\n", r.left,r.right);
    }
  }
  
  return r;
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
  for (i=0;i<*pairtable;i++){
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

