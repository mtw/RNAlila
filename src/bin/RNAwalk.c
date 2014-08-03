#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "rnalila/lila.h"


static void parse_infile(FILE *fp);
int get_list(struct_en*, struct_en*);
void ini_RNA(const char*);


model_detailsT md;
paramT *P = NULL;
vrna_fold_compound *vc = NULL;

int main() {
  int e,enew,emove;
  float mfe;
  short int *pt,*s0,*s1;
  double erange;
  
  move_opt.INFILE = stdin;
  parse_infile(move_opt.INFILE);
  ini_RNA(move_opt.sequence);
  srand(time(NULL));
  
  { // compute mfe
    mfe = vrna_fold(vc,NULL);
    destroy_fold_compound(vc);
    printf ("mfe = %6.2f\n",mfe);
  }
  
  pt = vrna_pt_get(move_opt.structure);
  s0 = get_sequence_encoding(move_opt.sequence,0,&(P->model_details));
  s1 = get_sequence_encoding(move_opt.sequence,1,&(P->model_details));
  		 
  //mtw_dump_pt(pt);
  //char *str = vrna_pt_to_db(pt);
  //printf(">%s<\n",str);
  e = vrna_eval_structure_pt(move_opt.sequence,pt,P);
  printf("%s\n", move_opt.sequence);
  print_str(stdout,pt);printf(" %6.2f\n",(float)e/100);

  m = get_random_move_pt(pt);
  emove = vrna_eval_move_pt(pt,s0,s1,m.left,m.right,P);
  apply_move_pt(pt,m);
  //mtw_dump_pt(pt);
  enew = e + emove;
  //printf ("performed move l:%4d r:%4d\t Energy +/- %6.2f\n",m.left,m.right,(float)emove/100);
  print_str(stdout,pt);printf(" %6.2f\n",(float)enew/100);
  //e = vrna_eval_structure_pt(move_opt.sequence,pt,P);
  //print_str(stdout,pt);printf(" %6.2f\n",(float)e/100);
  
  free(list);
  free(pt);
  free(s0);
  free(s1);
  free(move_opt.sequence);
  free(move_opt.structure);
  return 0;
}

void
ini_RNA (const char *seq)
{
  set_model_details(&md); /* use current global model */
  P = vrna_get_energy_contributions(md);
  vc = vrna_get_fold_compound(seq, &md,VRNA_OPTION_MFE);
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
  move_opt.sequence  = (char *) calloc (strlen(line)+1, sizeof(char));
  move_opt.structure = (char *) calloc (strlen(line)+1, sizeof(char));
  assert(move_opt.sequence != NULL); assert(move_opt.structure != NULL);
  sscanf(line, "%s", move_opt.sequence);
  free (line);
  line = get_line(fp);
  sscanf(line, "%s", move_opt.structure);
  free (line);
  move_opt.len = strlen(move_opt.sequence);

}

