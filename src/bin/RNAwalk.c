#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "RNAlila/lila.h"
#include "RNAlila/moves.h"
#include "RNAwalk_cmdl.h"


/* functions */
static void parse_infile(FILE *fp);
static void RNAwalk_memoryCleanup(void);
static void process_app_options(int args, char **argv);

/* structures */
typedef struct _rnawalk {
  FILE *input;            /* file pointer to input */
  char *basename;         /* base name of input */
  char walktype;          /* walk type (A|G|R) */
  int walklen;            /* walk length */
} local_optT;

static local_optT local_opt;
struct RNAwalk_args_info args_info;

int
main(int argc, char **argv)
{
  int e,enew,emove,len=0;
  float mfe;
  double erange;
  move_str m;
  char *newstruc=NULL;
  
  /* initialize ViennaRNA common options */
  lila_ini_vcd_options();

  /* process all application-specific options */
  process_app_options(argc,argv);

  /* TODO: get sequence from input file */
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

  printf ("input sequence:\n");
  printf ("%s\n",lilass.sequence);
  printf ("MFE is:\n");
  { // compute mfe
    char *ss=(char*)calloc(strlen((lilass.sequence)+1),sizeof(char));
    mfe = vrna_fold(vc,ss);
    vrna_free_fold_compound(vc);
    //   ss = vrna_pt_to_db(pt);
    printf ("%s (%6.2f)\n",ss,mfe);
    free(ss);
  }
  
  pt = vrna_pt_get(lilass.structure);
  e = vrna_eval_structure_pt(lilass.sequence,pt,P);
  
  printf ("Start structure:\n");
  printf ("%s (%6.2f)\n",lilass.structure, (float)e/100);
 
  // lila_dump_pt(pt);

  if(local_opt.walktype == 'R'){
    printf ("performing random walk\n");
    while(len<local_opt.walklen){
      /* make a random move */
      m = lila_random_move_pt(lilass.sequence,pt);
      /* compute energy difference for this move */
      emove = vrna_eval_move_pt(pt,s0,s1,m.left,m.right,P);
      /* evaluate energy of the new structure */
      enew = e + emove;
      /* do the move */
      lila_apply_move_pt(pt,m);
      /*
	newstruc = vrna_pt_to_db(pt);
      printf("%s (%6.2f)\n", newstruc, (float)enew/100);
      */
      print_str(stdout,pt);
      printf(" (%6.2f)\n", (float)enew/100);
      e = enew;
      len++;
    }
    free(newstruc);
  }
  else if (local_opt.walktype == 'G'){
    printf ("performing gradient walk\n");
    while(len<local_opt.walklen){
      /* make a random move */
      m = lila_gradient_move_pt(lilass.sequence,pt);
      /* compute energy difference for this move */
      emove = vrna_eval_move_pt(pt,s0,s1,m.left,m.right,P);
      /* evaluate energy of the new structure */
      enew = e + emove;
      /* do the move */
      lila_apply_move_pt(pt,m);
      /*
	newstruc = vrna_pt_to_db(pt);
      printf("%s (%6.2f)\n", newstruc, (float)enew/100);
      */
      print_str(stdout,pt);
      printf(" (%6.2f)\n", (float)enew/100);
      e = enew;
      len++;
    }
    
  }
  else if (local_opt.walktype == 'A'){
    printf ("performing adaptive walk\n");
  }
  else {
    printf ("unknown walktype, exiting ...\n");
    exit(EXIT_FAILURE);
  }
  //char *str = vrna_pt_to_db(pt);
  //printf(">%s<\n",str);

  /*    
    m = get_random_move_pt(pt);
    emove = vrna_eval_move_pt(pt,s0,s1,m.left,m.right,P);
    lila_apply_move_pt(pt,m);
    //mtw_dump_pt(pt);
    enew = e + emove;
    //printf ("performed move l:%4d r:%4d\t Energy +/- %6.2f\n",m.left,m.right,(float)emove/100);
    print_str(stdout,pt);printf(" %6.2f\n",(float)enew/100);
    //e = vrna_eval_structure_pt(move_opt.sequence,pt,P);
    //print_str(stdout,pt);printf(" %6.2f\n",(float)e/100);
    */
  RNAwalk_memoryCleanup();
  free(pt);
  free(s0);
  free(s1);

  return 0;
}

/**/
static void
process_app_options(int argc, char **argv)
{
  /* initialize local options */
  local_opt.input      = NULL;
  local_opt.walktype   = 'G';

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
     else {
       fprintf(stderr, "argument of --walktype must be A, G or R\n");
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
  free(lilass.sequence);
  free(lilass.structure);
}
