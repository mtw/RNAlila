/* ranstrucS.c */
/* count the number of secondary structures compatible with an RNA sequence */
/* and produce random structures compatible with the sequence */
/* Last changed Time-stamp: <2015-02-08 23:24:11 mtw> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <lila.h>

#include <gsl/gsl_rng.h>
#ifdef __MACH__
#include <mach/mach_time.h>
#define CLOCK_REALTIME 0
#define CLOCK_MONOTONIC 0
int clock_gettime(int clk_id, struct timespec *t){
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);
    uint64_t time;
    time = mach_absolute_time();
    double nseconds = ((double)time * (double)timebase.numer)/((double)timebase.denom);
    double seconds = ((double)time * (double)timebase.numer)/((double)timebase.denom * 1e9);
    t->tv_sec = seconds;
    t->tv_nsec = nseconds;
    return 0;
}
#else
#include <time.h>
#endif


#define MAXLENGTH 500
#define TURN 3

int L=1;           /* minimum helix length */
double p=1.;       /* stickiness */

void       count(void);
extern double      drand48();

short *S = NULL;
double **Q = NULL;
double *PP = NULL;

int Nmax;
unsigned short xsubi[3];

char *seq, *s;
int length, n_of_structs=1;

static unsigned long seed;    /* random seed */
static gsl_rng *rng = NULL;     /* GSL random number generator */
static struct timespec ts;    /* timespec struct for random seed */
static double rnum;           /* random number */

/**/
void
count()
{
  double pp;
  int i,j,n,k;

  for (i=1; i<=length; i++)
    Q[i][i-1] = 1;

  for (i=length; i>0; i--)
    for (j=i; j<=length; j++) {
      Q[i][j] = Q[i][j-1];
      for (k=i; k<=j-TURN-1; k++)
	if (md.pair[s0[k]][s0[j]])
	  Q[i][j] += Q[i][k-1]*Q[k+1][j-1];
    }
}

/**/
char *
random_struc(int n)
{
  char *seq, *s;
  int n_of_structs=1;

  int i,j,k,l,ns,psi=0;
  double r,SS,pk;
  char *struc;
  struct {
    int i,j;
  } sector[MAXLENGTH/2+1];

  if (n>Nmax) {Nmax=n; count();}

  struc = (char *) calloc(sizeof(char),(n+2));

  sector[0].i=i=1;
  sector[0].j=j=n;
  ns=1;

  do {
    if (j-i<TURN+1) {
      for (l=i-1; l<j; l++) struc[l]='.';
      i=sector[--ns].i;
      j=sector[ns].j;
      psi=0;
      continue;
    }
    rnum =  gsl_rng_uniform (rng);
    // printf(" >> rnum = %6.4g\n",rnum);
    r = (rnum*Q[i][j]);
    if (r<Q[i+1][j]) {
      struc[i-1]='.';
      i++;
    } else {
      struc[i-1]='(';
      SS=Q[i+1][j];
      for (k=i+TURN; k+1<j; k++) { /* pair i with k+1 */
	if (md.pair[s0[i]][s0[k+1]]) {
	  SS += Q[i+1][k]*Q[k+2][j];
	  if (r<SS) break;
	}
      }
      struc[k]=')';
      if (k+2<=j) {
	sector[ns].i = k+2;
	sector[ns++].j= j;
      }
      j = k;
      i = i+1;
    }
  } while(ns>0);

  return struc;
}

/**/
char *
lila_random_structureS(char *seq)
{
  int i, display = 0;
  int count_only=0;
  char *line;
  
  length = strlen(seq);
  
  Q = (double **) calloc(sizeof(double *), length+1);
  for (i=0; i<=length; i++) {
    Q[i] = (double *) calloc(sizeof(double), length+1);
    if (!Q[i]){
      fprintf(stderr, "Out of memory\n");
      exit (EXIT_FAILURE);
    }
  }

  /* prepare gsl random-number generation */
  (void) clock_gettime(CLOCK_REALTIME, &ts);
  //if(wanglandau_opt.seed_given){
  //  seed = wanglandau_opt.seed;
  //}
  //else {
  seed =   ts.tv_sec ^ ts.tv_nsec;
    //}
  //fprintf(stderr, "initializing random seed: %d\n",seed);
  gsl_rng_env_setup();
  rng = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set( rng, seed );
  /* end gsl */

  s=random_struc(length);

  for (i=0;i<=length;i++)
    free (Q[i]);
  free(Q);

  gsl_rng_free(rng);
  return s;  
}



/* End of file */
