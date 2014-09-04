/*
  lila.c: common routines for RNAlila
  Last changed Time-stamp: <2014-09-04 12:28:48 mtw>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <lila.h>

/* ==== */
void 
lila_ini_vcd_options(void)
{
  vcd.temperature         = VRNA_MODEL_DEFAULT_TEMPERATURE;
  vcd.betaScale           = VRNA_MODEL_DEFAULT_BETA_SCALE;
  vcd.dangles             = VRNA_MODEL_DEFAULT_DANGLES;
  vcd.noLP                = VRNA_MODEL_DEFAULT_NO_LP;
}

/* ==== */
void
lila_ini_vRNA (const char *seq)
{
  P  = vrna_get_energy_contributions(md);
  vc = vrna_get_fold_compound(seq, &md,VRNA_OPTION_MFE);
  s0 = vrna_seq_encode_simple(seq,&(P->model_details));
  s1 = vrna_seq_encode(seq,&(P->model_details));
}

/* === */
void
lila_vRNA_cleanup(void)
{
  free(s0);
  free(s1);
  free(P);
  vrna_free_fold_compound(vc);
}


/* ==== */
/*  set temperature, dangles, betaScale, noLP as vRNA modelDetails*/
void
lila_set_vcd_options(const unsigned int temp_given,
		     const unsigned int betaScale_given,
		     const unsigned int dangles_given,
		     const unsigned int noLP_given,
		     const double temp_arg,
		     const double betaScale_arg,
		     const int dangles_arg,
		     const int noLP_flag)
{
  /* temperature */
  if(temp_given){
    if( (vcd.temperature = temp_arg) < -K0 ){
      fprintf(stderr, "Value of --temp must be > -273.15\n");
      exit (EXIT_FAILURE);
    }
    vrna_md_set_temperature(&md,vcd.temperature);
  }
  /* dangles */
  if(dangles_given){
    if( (vcd.dangles = dangles_arg) < 0 ){
      fprintf(stderr, "Value of --dangles must be > 0\n");
      exit (EXIT_FAILURE);
    }
    else if ( (vcd.dangles = dangles_arg) > 3 ){
      fprintf(stderr, "Value of --dangles must be <= 3\n");
      exit (EXIT_FAILURE);
    }
    vrna_md_set_dangles(&md,vcd.dangles);
  }
  /* betaScale */
  if(betaScale_given){
    vcd.betaScale = betaScale_arg;
    vrna_md_set_betascale(&md,vcd.betaScale);
  }
  /* noLP */
  if(noLP_given){
    vcd.noLP = noLP_flag; /* 0 or 1 */
    vrna_md_set_nolp(&md,vcd.noLP);
  }
}


/**/
/* TODO: implement lila_parse_seq_struc(FILE *fp) in rnalila/io.c */
void
lila_parse_seq_struc(FILE *fp)
{
  char *line=NULL;
  
  line = get_line(fp);
  /* skip comment lines */
  while ((*line == '*')||(*line == '\0')||(*line == '>')) {
    free(line);
    line = get_line(fp);
  }
  lilass.sequence  = (char *) calloc (strlen(line)+1, sizeof(char));
  lilass.structure = (char *) calloc (strlen(line)+1, sizeof(char));
  assert(lilass.sequence != NULL); assert(lilass.structure != NULL);
  sscanf(line, "%s", lilass.sequence);
  free (line);
  line = get_line(fp);
  sscanf(line, "%s", lilass.structure);
  free (line);
  lilass.length = strlen(lilass.sequence);
}

/**/
char *
lila_basename(char *arg)
{
  int len;
  char *s=NULL, *t=NULL;
  
  s = strdup(arg);
  len = strlen(s);
  t = rindex(s, '/');
  if (t != NULL) memmove(s, t+1, (len-(t-s))*sizeof(char));
  t = NULL; t = index(s, '.');
  if (t != NULL) *t = '\0';
  return s;
}

/**/
int
lila_cmp_db(void *a,
	    void *b)
{
  return strcmp(a,b);	
}

/*
  Compare function for sseT (see lila.h); Compares (dot bracket)
  structures lexicographically.
*/
int
lila_cmp_sse_lex(const void *a, const void *b)
{
  const Lila2seT *ma = a;
  const Lila2seT *mb = b;

  return strcmp(ma->structure,mb->structure);
}

/*
  Compare function for sseT (see lila.h); Compares energies.
*/
int
lila_cmp_sse_en(const void *a,
		const void *b)
{
  const Lila2seT *ma = a;
  const Lila2seT *mb = b;

  if (ma->energy > mb->energy)
    return 1;
  else if (ma->energy < mb->energy)
    return -1;
  else
    return 0;
}

/*
  Compare function for sseT (see lila.h); Compares first by energy,
  then by lexicographical order of structure.
*/
int
lila_cmp_sse_lexen(const void *a,
		   const void *b)
{
  const Lila2seT *ma = a;
  const Lila2seT *mb = b;

  int comp = ma->energy - mb->energy;

  if (comp < 0)
    return -1;
  
  if (comp > 0)
    return 1;
  
  comp = strcmp(ma->structure, mb->structure);
  return comp;
}

/*
  Dump one element of a Lila2seT list.
*/
void
lila_print_2se(gpointer data,
	       gpointer user_data)
{
  Lila2seT *foo = (Lila2seT *)data;
  fprintf(stderr,"%s (%6.4f) CC\n",foo->structure,foo->energy);
}

/*
  Free space allocated for a Lils2seT
*/
void
lila_free_cc(gpointer data)
{
  Lila2seT *foo = (Lila2seT *)data;
  free(foo->structure);
}
