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
void ini_ViennaRNA(const char*);
static void RNAwalk_memoryCleanup(void);
static void ini_RNAwalk(void);
static void ini_globs(void);

/* structures */
typedef struct _rnawalk {
  FILE *INFILE;
  char *seq;
  char *struc;
  int len;
  double temp;
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

 
  process_options(argc,argv);
  rnawalk_opt.INFILE = stdin;
  parse_infile(rnawalk_opt.INFILE);
  ini_ViennaRNA();
  srand(time(NULL));
  
  { // compute mfe
    mfe = vrna_fold(vc,NULL);
    destroy_fold_compound(vc);
    printf ("mfe = %6.2f\n",mfe);
  }
  
  pt = vrna_pt_get(rnawalk_opt.my_struc);
  s0 = get_sequence_encoding(rnawalk_opt.my_seq,0,&(P->model_details));
  s1 = get_sequence_encoding(rnawalk_opt.my_seq,1,&(P->model_details));
  		 
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
process_options(void)
{
  ini_globs();
  if (RNAwalk_cmdline_parser (argc, argv, &args_info) != 0){
    fprintf(stderr, "error while parsing command-line options\n");
    exit(EXIT_FAILURE);
  }

  process_vcd_options(); /* temperature, dangles, noLP, betascale */
 
  
}


/**/
void
ini_ViennaRNA (const char *seq)
{
  set_model_details(&md); /* use current global model */
  P = vrna_get_energy_contributions(md);
  vc = vrna_get_fold_compound(seq, &md,VRNA_OPTION_MFE);
}

/**/
static void
ini_globs(void)
{
  rnawalk_opt.INFILE      = NULL;
  vcd.temperature         = VRNA_TEMP_DEFAULT;
}


/**/
static void
parse_infile(FILE *fp)
{
  char *line=NULL;
  
  line = get_line(fp);
  /* skip comment lines */
  while ((*line == '*')||(*line == '\0')||(*line == '>')) {
    free(line);
    line = get_line(fp);
  }
  my_seq  = (char *) calloc (strlen(line)+1, sizeof(char));
  my_struc = (char *) calloc (strlen(line)+1, sizeof(char));
  assert(my_seq != NULL); assert(my_struc != NULL);
  sscanf(line, "%s", my_seq);
  free (line);
  line = get_line(fp);
  sscanf(line, "%s", my_struc);
  free (line);
  my_len = strlen(my_seq);
}

/**/
static void
RNAwalk_memoryCleanup (void)
{
  free(rnawalk_opt.my_seq);
  free(rnawalk_opt.my_struc);
}
