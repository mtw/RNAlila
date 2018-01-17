/*
  RNAwalk.c
  Last changed Time-stamp: <2017-06-11 10:09:10 michl>
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "RNAlila/lila.h"
#include "RNAlila/moves.h"
#include "RNAlila/topology.h"
#include "RNAwalk_cmdl.h"

static void parse_infile(FILE *fp);
static void RNAwalk_memoryCleanup(void);
static void process_app_options(int args, char **argv);
static void ini_AWmin(void);
static void AWmin(short int *);
static void AWmin_memoryCleanup(void);
static void dump_items(gpointer data,  gpointer user_data);
static void free_key(gpointer data);
static void walk(int);

typedef struct _rnawalk {
  FILE *input;            /* file pointer to input */
  char *basename;         /* base name of input */
  char  walktype;         /* walk type (A|G|R) */
  int   walklen;          /* walk length */
  bool  awmin;            /* flag for awmin */
} local_optT;

GHashTable *S, *M;
static local_optT local_opt;
struct RNAwalk_args_info args_info;

int
main(int argc, char **argv)
{
  int energy;
  float mfe;
  double erange;
  
  /* initialize ViennaRNA common options */
  lila_ini_vcd_options();

  /* process all application-specific options */
  process_app_options(argc,argv);

  /* TODO: use ViennaRNA routines to parse input file */
  lila_parse_seq_struc(local_opt.input);

  /* set default vRNA model details */
  vrna_md_set_default(&md);      

  /* adjust common vRNA options */
  lila_set_vcd_options(args_info.temp_given,
		       args_info.betaScale_given,
		       args_info.dangles_given,
		       args_info.noLP_given,
		       args_info.temp_arg,
		       args_info.betaScale_arg,
		       args_info.dangles_arg,
		       args_info.noLP_flag); 
  lila_ini_vRNA(lilass.sequence);
  srand(time(NULL));

  /*
    printf ("%s\n",lilass.sequence);
    printf ("MFE is:\n");
    { // compute mfe
    char *ss=(char*)calloc((strlen(lilass.sequence)+1),sizeof(char));
    mfe = vrna_fold(vc,ss);
    vrna_free_fold_compound(vc);
    printf ("%s %6.2f\n",ss,mfe);
    free(ss);
    }
  */
  
  pt = vrna_ptable(lilass.structure);
  energy = vrna_eval_structure_pt(vc,pt); 
  // this could be done directly via vrna_eval_structure_simple
  //float e2 = vrna_eval_structure(lilass.sequence,lilass.structure,P);
  
  /* print start structure */
  printf ("%s\n",lilass.sequence);
  printf ("%s %6.2f S\n",lilass.structure, (float)energy/100);
  
  walk(energy);

  RNAwalk_memoryCleanup();
  
  /* free ViennaRNA data */
  lila_vRNA_cleanup();
  
  free(pt);
  
  return 0;
}
  

static void
walk (int e)
{
  int enew,emove,ismin,len = 0;
  move_str m;
  
  if(local_opt.walktype == 'R'){
    /* printf ("performing random walk\n");*/
    while(len<local_opt.walklen){
      /* make a random move */
      m = lila_random_move_pt(lilass.sequence,pt);
      /* compute energy difference for this move */
      emove = vrna_eval_move_pt(vc,pt,m.left,m.right);
      /* evaluate energy of the new structure */
      enew = e + emove;
      /* do the move */
      lila_apply_move_pt(pt,m);
      print_str(stdout,pt);
      printf(" %6.2f\n", (float)enew/100);
      e = enew;
      len++;
    }
  }
  else if (local_opt.walktype == 'G'){
    /* printf ("performing gradient walk\n");*/
    while(len<local_opt.walklen){
      char status[] = "x";
      /* backup current structure */
      short int *ptbak = vrna_ptable_copy(pt);
      /* return gradient move */
      m = lila_gradient_move_pt(lilass.sequence,pt);
      /* compute energy difference for this move */
      emove = vrna_eval_move_pt(vc,pt,m.left,m.right);
      /* evaluate energy of the new structure */
      enew = e + emove;
      /* do the move */
      lila_apply_move_pt(pt,m);
      if (emove > 0.){
	break;
      }
      else {
	if (memcmp(ptbak,pt,sizeof(short int)*(lilass.length+1)) == 0){
	  fprintf(stderr,"gradient walk ended in same structure, finishing\n");
	  break;
	}
      }
      
      {
	/* validate topological status: transient or minimum?  NOTE:
	   this is expensive since ALL neighbors are re-generated and
	   re-evaluated */
	ismin = lila_is_minimum_or_shoulder_pt(lilass.sequence,pt);
	if (ismin == 1)
	  status[0]='*';
	else if (ismin == 0)
	  status[0]='T';
	else
	  status[0]='D';
      }
      
      print_str(stdout,pt);
      printf(" %6.2f %s\n", (float)enew/100, status);
      e = enew;
      len++;
    }
  }
  else if (local_opt.walktype == 'A'){
    
    /* if(local_opt.awmin==true) { */
    /*   GList *L=NULL; */
    /*   /\* printf ("performing AWmin\n");*\/ */
    /*   ini_AWmin(); */
    /*   AWmin(pt); */
    /*   /\* fprintf(stdout,"Minima:\n");*\/ */
    /*   /\* L = g_hash_table_get_keys(M); *\/ */
    /*   /\* g_list_foreach(L, dump_items, "struc %s\n"); *\/ */
    /*   /\* fprintf(stderr,"detroying hash M\n");*\/ */
    /*   g_hash_table_destroy(M); */
    /*   /\* fprintf(stderr,"destroying hash S\n");*\/ */
    /*   g_hash_table_destroy(S); */
    /*   g_list_free(L); */
    /*   //AWmin_memoryCleanup(); */
    /* } */
    /* else{ */
     
    /* printf ("performing adaptive walk\n");*/
    while(len<local_opt.walklen){
      char status[] = "x";
      /* make a adaptive move */
      m = lila_adaptive_move_pt(lilass.sequence,pt);
      /* no further moves possible */
      if (m.left == 0 && m.right == 0){ 
	fprintf(stderr, "no further moves posible after %i moves\n",len);
	break;
      }
      /* compute energy difference for this move */
      emove = vrna_eval_move_pt(vc,pt,m.left,m.right);
      /* evaluate energy of the new structure */
      enew = e + emove;
      /* do the move */
      lila_apply_move_pt(pt,m);
      
      {
	/* validate topological status: transient or minimum?  NOTE:
	   this is expensive since ALL neighbors are re-generated and
	   re-evaluated */
	ismin = lila_is_minimum_or_shoulder_pt(lilass.sequence,pt);
	if (ismin == 1)
	  status[0]='*';
	else if (ismin == 0)
	  status[0]='T';
	else
	  status[0]='D';
      }
      
      print_str(stdout,pt);
      printf(" %6.2f %s\n", (float)enew/100, status);
      e = enew;
      len++;
    }
  }
  else if (local_opt.walktype == 'N'){
    /* printf ("computing neighbors only\n");*/
    GQueue *nb = NULL;
    nb = lila_generate_neighbors_pt((const char *)lilass.sequence,pt);
    lila_output_dbe_gqueue(nb, NULL);
    fprintf(stderr,"=============\n");
 
    { 
      LilaDBE *bar = NULL;
      while ( (bar = (LilaDBE*)g_queue_pop_head(nb)) != NULL){
	int min;
	GList *gliste;
	short int *pt = vrna_ptable(bar->structure);
	gliste = lila_get_cc_pt(lilass.sequence,pt,&min);
	lila_output_dbe_glist(gliste, "CC");
	lila_dealloc_dbe_glist(gliste);
	free(pt);
	fprintf(stderr,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
	lila_free_dbe(bar);
      }  
    }
    lila_dealloc_dbe_gqueue(nb);
  }
  else {
    printf ("unknown walktype, exiting ...\n");
    exit(EXIT_FAILURE);
  }
}

static void
dump_items(gpointer data, gpointer user_data)
{
  float energy;
  char *structure = (char*)data;
  energy = vrna_eval_structure(vc,structure);
  fprintf(stdout,"%s %6.2f\n", structure, energy);
}

static void
free_key(gpointer data)
{
  //fprintf(stderr, "[[free_key]] freeing %s\n", data);
  free(data);
}


/**/
/* static void */
/* ini_AWmin(void) */
/* { */
/*   /\* initialize S as glib hash. S is a key==value type hash containing */
/*      just secondary structures in dot bracket notation *\/ */
/*   S = g_hash_table_new_full(g_str_hash,g_str_equal,free_key,NULL); */
/*   //S = g_hash_table_new(g_str_hash,g_str_equal); */
/*   /\* initialize list of local minima M *\/ */
/*   M = g_hash_table_new(g_str_hash,g_str_equal);   */
/* } */

/*
  Compute all minima reachable from a single structure via adaptive
  walks via recursive depth first search.
*/
/* static void */
/* AWmin(short int *pt) */
/* { */
/*   int i=0,size,count; */
/*   move_str *moves=NULL; */
/*   char *struc; */
/*   short int *ptbak=NULL; */
/*   static int counter=0; */

/*   /\* fprintf(stderr,"\n[[AWmin]] [[%i]]\n",counter++); *\/ */
/*   struc = lila_db_from_pt(pt); */
/*   if(struc == NULL) */
/*     fprintf(stderr,"struc is a NULl pointer; this shouldn't happen\n"); */
/*   /\* fprintf(stderr,"  %s INSERTED ",struc); *\/ */
/*   if(g_hash_table_contains(S,struc)) */
/*     return; */
/*   else */
/*     g_hash_table_add(S,struc); /\* key == value *\/ */
/*   /\* fprintf(stderr,"[hashsize=%i]\n", g_hash_table_size(S); *\/ */
  
/*   /\* get move operations for all AW neighbors, INCUDING DEGENERATE */
/*      NEIGHBORS (!) *\/ */
/*   moves = lila_all_adaptive_moves_pt(lilass.sequence,pt,&count,1); */
/*   /\* //count number of adaptive move operation */
/*      { */
/*      int k; */
/*      fprintf(stderr,"%i adaptive move operations possible\n",count); */
/*      for (k=0;k<count;k++){ */
/*      fprintf(stderr," move #%i l:%3i|r:%3i\n",k,moves[k].left,moves[k].right); */
/*      }    */
/*      } */
/*   *\/ */
 
/*   if(count == 0){ /\* no adaptive walks neighbors found *\/ */
/*     if (lila_is_minimum_or_shoulder_pt(lilass.sequence,pt) == 1){ */
/*       /\* it's a true minimun *\/ */
/*       int e; */
/*       e = vrna_eval_structure_pt(vc,pt); */
/*       /\* add struc to the list of minima *\/ */
/*       /\* fprintf(stderr,"M %s %6.2f ADDED TO MINIMA\n",struc,(float)e/100); *\/ */
/*       g_hash_table_add(M,struc); */
/*     } */
/*     else{ */
/*       fprintf(stderr, "ERROR: no AW neighbors found, but structure isn't a minimum either\n"); */
/*       fprintf(stderr, "This shouldn't happen. Exiting ...\n"); */
/*       fprintf(stderr, "%s\n",struc); */
/*     } */
/*   } */
/*   else if (lila_is_minimum_or_shoulder_pt(lilass.sequence,pt) == -1) { */
/*     /\* it's a degenerate minimum or a shoulder *\/ */
/*     GList *f,*conncomp; */
/*     LilaDBE *l = NULL; */
/*     char *lexmin = NULL; */
/*     int ismin; */
  
/*     fprintf(stderr, "%s D\n",struc); */
/*     conncomp = lila_get_cc_pt(lilass.sequence,pt,&ismin); */
/*     lila_output_dbe_glist(conncomp, "CC"); */
/*     if (ismin == 1) */
/*       fprintf(stderr, "MINIMUM COMPONENT\n"); */
/*     else */
/*       fprintf(stderr, "SHOULDER COMPONENT\n"); */
/*     f = lila_lexmin_dbe_glist(conncomp); */
/*     l = (LilaDBE *)f; */
/*     lexmin = strdup(l->structure); */
/*     g_hash_table_add(M,lexmin); */
/*     g_list_foreach(conncomp,(GFunc)lila_dbe_structure2ghashtable,S); */
/*     lila_dealloc_dbe_glist(conncomp); */
/*     fprintf(stderr,"---\n"); */
/*   } */
/*   else{ */
/*     for(i=0;i<count;i++){ /\* loop over all move operations *\/ */
/*       ptbak = vrna_ptable_copy(pt); */
/*       char *v; */
/*       /\* */
/* 	v = lila_db_from_pt(pt); */
/* 	fprintf(stderr,"S %s | move #%i l:%3i|r:%3i\n",v,i,moves[i].left,moves[i].right); */
/* 	free(v); */
/*       *\/ */
/*       lila_apply_move_pt(pt,moves[i]); */
/*       v = lila_db_from_pt(pt); */
/*       /\* */
/* 	fprintf(stderr,"T %s\n",v); */
/*       *\/ */
/*       if(g_hash_table_lookup(S,v)){ */
/* 	/\* */
/* 	  fprintf(stderr,"  %s already processed\n",v); */
/* 	*\/ */
/* 	free(v); */
/*       } */
/*       else{ */
/* 	free(v); */
/* 	/\* */
/* 	  fprintf(stderr,"calling AWmin()\n"); */
/* 	*\/ */
/* 	AWmin(pt); */
/*       } */

/*       /\* fprintf(stderr, "RETURNING\n"); *\/ */
/*       /\* reset to the version of pt this function was called with *\/ */
/*       /\* */
/* 	v = lila_db_from_pt(pt); */
/* 	fprintf(stderr, "B %s [before resetting]\n",v); */
/* 	free(v); */
/*       *\/ */
/*       pt = ptbak; */
/*       /\* */
/* 	v = lila_db_from_pt(pt); */
/* 	fprintf(stderr, "B %s [after resetting]\n",v); */
/* 	free(v); */
/*       *\/ */
/*     } /\* end for *\/ */
/*   } /\* end else *\/ */
/*   free(ptbak); */

/*   if(moves != NULL) free(moves); */
/* } */

/**/
/* static void */
/* AWmin_memoryCleanup(void) */
/* { */
/*   // g_hash_table_destroy(S); */
/*   //g_hash_table_destroy(M); */
/* } */

/**/
static void
process_app_options(int argc, char **argv)
{
  /* initialize local options */
  local_opt.input      = NULL;
  local_opt.walktype   = 'G';
  local_opt.walklen    = 100;
  local_opt.awmin      = false;

  /* parse command line, overweite local options */
  if (RNAwalk_cmdline_parser (argc, argv, &args_info) != 0){
    fprintf(stderr, "error while parsing command-line options\n");
    exit(EXIT_FAILURE);
  } 

  /* walktype */
  if(args_info.walktype_given){
    if( strncmp(args_info.walktype_arg, "A", 1)==0 )
      local_opt.walktype = 'A';
    else if(  strncmp(args_info.walktype_arg, "G", 1)==0 )
      local_opt.walktype = 'G';
    else if(  strncmp(args_info.walktype_arg, "R", 1)==0 )
      local_opt.walktype = 'R';
    else if(  strncmp(args_info.walktype_arg, "N", 1)==0 )
      local_opt.walktype = 'N';
     else {
       fprintf(stderr, "argument of --walktype must be A, G, R or N\n");
       exit(EXIT_FAILURE);
     }
  }

  /* walk length */
  if(args_info.walklength_given){
    if( (local_opt.walklen = args_info.walklength_arg)  < 1) {
      fprintf(stderr, "argument of --walklength must >= 1\n");
      exit(EXIT_FAILURE);
    }
  }

  /* compute awmin */  
  /* if(args_info.awmin_given){ */
  /*   if(args_info.awmin_given){ */
  /*     local_opt.awmin = true; */
  /*   } */
  /* } */

  /* input file */
  if (args_info.inputs_num){
    char *infile=NULL;
    infile = strdup(args_info.inputs[0]);
    local_opt.basename = lila_basename(args_info.inputs[0]);
    local_opt.input = fopen(infile, "r");
    free(infile);
  }
  else
    local_opt.input = stdin;
}

/**/
static void
RNAwalk_memoryCleanup (void)
{
  RNAwalk_cmdline_parser_free(&args_info);
  if (args_info.inputs_num){
    free(local_opt.basename);
  }
  fclose(local_opt.input);
  free(lilass.sequence);
  free(lilass.structure);
}
