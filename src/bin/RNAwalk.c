#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "rnalila/lila.h"
#include "rnalila/moves.h"
#include "RNAwalk_cmdl.h"


/* functions */
static void parse_infile(FILE *fp);
int get_list(struct_en*, struct_en*);
static void RNAwalk_memoryCleanup(void);
static void process_app_options(void);

/* structures */
typedef struct _rnawalk {
  FILE *INFILE;
  int len;                /* walk length */
} rnawalk_optT;

static rnawalk_optT rnawalk_opt;
struct RNAwalk_args_info args_info;

int
main(int argc, char **argv)
{
  int e,enew,emove;
  float mfe;
  double erange;
  move_str m;

  lila_ini_vcd_options();
  process_app_options(argc,argv); /* app-specific options*/
  rnawalk_opt.INFILE = stdin;
  /* TODO: get sequence from input file */
  lila_parse_seq_struc(rnawalk_opt.INFILE);   /* process input */
  vrna_md_set_default(&md);       /* set default vRNA model details */
  lila_set_vcd_options(args_info.temp_given, /* adjust common vRNA options */
		       args_info.betaScale_given,
		       args_info.dangles_given,
		       args_info.noLP_given,
		       args_info.temp_arg,
		       args_info.betaScale_arg,
		       args_info.dangles_arg,
		       args_info.noLP_flag); 
  lila_ini_vRNA(lilass.startseq);
  srand(time(NULL));
  
  { // compute mfe
    mfe = vrna_fold(vc,NULL);
    destroy_fold_compound(vc);
    printf ("mfe = %6.2f\n",mfe);
  }
  
  pt = vrna_pt_get(lilass.startstruc);
  /*
    s0 = get_sequence_encoding(rnawalk_opt.my_seq,0,&(P->model_details));
    s1 = get_sequence_encoding(rnawalk_opt.my_seq,1,&(P->model_details));
  */
		 
  //mtw_dump_pt(pt);
  //char *str = vrna_pt_to_db(pt);
  //printf(">%s<\n",str);

  /*
    e = vrna_eval_structure_pt(rnawalk_opt.my_seq,pt,P);
    printf("%s\n", rnawalk_opt.my_seq);
    print_str(stdout,pt);printf(" %6.2f\n",(float)e/100);
    
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
process_app_options(void)
{
  rnawalk_opt.INFILE      = NULL;
  if (RNAwalk_cmdline_parser (argc, argv, &args_info) != 0){
    fprintf(stderr, "error while parsing command-line options\n");
    exit(EXIT_FAILURE);
  } 

  /* TODO non vcd options */

  RNAwalk_cmdline_parser_free(&args_info);
}

/**/
static void
RNAwalk_memoryCleanup (void)
{
  free(lilass.startseq);
  free(lilass.startstruc);
}
