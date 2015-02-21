/*
  RNArandstruc.c
  Last changed Time-stamp: <2015-02-08 22:50:08 mtw>
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
#include "RNArandstruc_cmdl.h"

static void parse_infile(FILE *fp);
static void RNArandstruc_memoryCleanup(void);
static void process_app_options(int args, char **argv);
static void ini_AWmin(void);
static void AWmin(short int *);
static void AWmin_memoryCleanup(void);
static void dump_items(gpointer data,  gpointer user_data);
static void free_key(gpointer data);
static void walk(void);

typedef struct _rnawalk {
  FILE *input;            /* file pointer to input */
  char *basename;         /* base name of input */
  int   number;           /* number of random structures */
} local_optT;

GHashTable *S, *M;
static local_optT local_opt;
struct RNArandstruc_args_info args_info;

int
main(int argc, char **argv)
{
  int i;
  char *form = NULL;
  float e;
  
  /* initialize ViennaRNA common options */
  lila_ini_vcd_options();

  /* process all application-specific options */
  process_app_options(argc,argv);

  /* TODO: use ViennaRNA routines to parse input file */
  lila_parse_sequence(local_opt.input);

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
  //srand(time(NULL));

  printf ("%s\n",lilass.sequence);
  
  for (i=0;i<local_opt.number;i++){
    form = lila_random_structureS(lilass.sequence);
    e=vrna_eval_structure(lilass.sequence,form,P);
    printf ("%s %6.2f\n",form, e);
  }

  
  RNArandstruc_memoryCleanup();
  
  /* free ViennaRNA data */
  lila_vRNA_cleanup();
  
  free(pt);
  free(form);
  return 0;
}
  

static void
free_key(gpointer data)
{
  //fprintf(stderr, "[[free_key]] freeing %s\n", data);
  free(data);
}


/**/
static void
process_app_options(int argc, char **argv)
{
  /* initialize local options */
  local_opt.input      = NULL;
  local_opt.number     = 1;

  /* parse command line, overweite local options */
  if (RNArandstruc_cmdline_parser (argc, argv, &args_info) != 0){
    fprintf(stderr, "error while parsing command-line options\n");
    exit(EXIT_FAILURE);
  } 

  /* number of structures */
  if(args_info.number_given){
    if( (local_opt.number = args_info.number_arg)  < 1) {
      fprintf(stderr, "argument of --number must >= 1\n");
      exit(EXIT_FAILURE);
    }
  }
  
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
RNArandstruc_memoryCleanup (void)
{
  RNArandstruc_cmdline_parser_free(&args_info);
  if (args_info.inputs_num){
    free(local_opt.basename);
  }
  fclose(local_opt.input);
  free(lilass.sequence);
  free(lilass.structure);
}
