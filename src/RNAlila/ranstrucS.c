/* ranstrucS.c */
/* count the number of secondary structures compatible with an RNA sequence */
/* and produce random structures compatible with the sequence */
/* Last changed Time-stamp: <2015-02-04 18:21:52 mtw> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <errno.h>
#include "ViennaRNA/pair_mat.h"

/* RNA package variables used in pair_mat.h */

#define PRIVATE static
#define MAXLENGTH 500
#define TURN 3

int L=1;           /* minimum helix length */
double p=1.;       /* stickiness */

PRIVATE void       count(void);
extern double      drand48();

short *S = NULL;
double **P = NULL;
double *PP = NULL;

int Nmax;
unsigned short xsubi[3];


char *seq, *s;
int length, n_of_structs=1;

PRIVATE void count()
{
  double pp;
  int i,j,n,k;
  
  for (i=1; i<=length; i++) 
    P[i][i-1] = 1;
  
  for (i=length; i>0; i--)
    for (j=i; j<=length; j++) {
      P[i][j] = P[i][j-1];
      for (k=i; k<=j-TURN-1; k++)
	if (pair[S[k]][S[j]])
	  P[i][j] += P[i][k-1]*P[k+1][j-1];
    }
}

PRIVATE char *random_struc(int n)
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
    r = (drand48()*P[i][j]);
    if (r<P[i+1][j]) {
      struc[i-1]='.';
      i++;
    } else {
      struc[i-1]='(';
      SS=P[i+1][j];
      for (k=i+TURN; k+1<j; k++) { /* pair i with k+1 */
	if (pair[S[i]][S[k+1]]) {
	  SS += P[i+1][k]*P[k+2][j];
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


void b2d(char *s)
{
   int i,j,o;
   char ch='A';

   for (i=0; i<strlen(s); i++) {
      if (s[i]=='.') s[i]=' ';
      if (s[i]=='(') {
	 for (o=1, j=i+1; o>0; j++) {
	    if (s[j]==')') o--;
	    if (s[j]=='(')o++;
	 }
	 if ((s[i-1]!=ch)||(s[j]!=ch)) ch++;
	 s[--j]=s[i]=ch;
      }
   }
}

int lila_random_structureS(char *seq)
{
  int i, display = 0;
  int count_only=0;
  char *line;
  // make_pair_matrix();
  
  length = strlen(seq);
  S = (short *) calloc(sizeof(short),length+1);
    
  for (i=1; i<=length; i++)  //make numerical encoding of sequence 
    S[i] = ENCODE(seq[i-1]);

  P = (double **) calloc(sizeof(double *), length+1);
  for (i=0; i<=length; i++) {
    P[i] = (double *) calloc(sizeof(double), length+1);
    if (!P[i]){
      fprintf(stderr, "Out of memory\n");
      exit (EXIT_FAILURE);
    }
  }
  count();

  if (count_only) {
    s=random_struc(length);
    for (i=1; i<=length; i++) {
      printf("%3d %g\n", i, P[1][i]);
    }
    free(s);
    return 0;
  } 
  fprintf(stderr, "%g structures\n", P[1][length]);
  
  srand48(time(NULL));
  
  for (i=0; i<n_of_structs; i++) {
    s=random_struc(length);
    if (display) b2d(s);
    printf("%s\n",s);
    free(s);
  }
  
  
}



/* End of file */
