/*
  lila.h : common definitions for RNAlila
  Last changed Time-stamp: <2014-09-03 23:44:07 mtw>
*/

#ifndef __RNA_LILA_H__
#define __RNA_LILA_H__

#include "config.h"
#include <glib.h>
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/structure_utils.h"
#include "ViennaRNA/move_set.h"
#include "ViennaRNA/subopt.h"
#include "ViennaRNA/energy_const.h"

/* structs */
typedef struct _vRNAcommon {
  double temperature;
  double betaScale;
  int dangles;
  int noLP;
} vRNAcommonT;

typedef struct _lilass {
  char *sequence;
  char *structure;
  int length;
} lilassT;                 /* sequence, secondary structure and length */       

typedef struct _lila2se {
  char *structure;
  float energy;
} Lila2seT;                 /* secondary structure and energy */

/* variables */
paramT *P;
model_detailsT md;
vrna_fold_compound *vc;
short int *pt,*s0,*s1;
vRNAcommonT vcd;         /* ViennaRNA common (model) details */
lilassT lilass;          /* start sequence && start structure */

/* functions */
char *lila_basename(char *);
void lila_parse_seq_struc(FILE *fp);
void lila_ini_vRNA(const char *seq);
void lila_vRNA_cleanup(void);
void lila_ini_vcd_options(void);
void lila_set_vcd_options(const unsigned int temp_given,
			  const unsigned int betaScale_given,
			  const unsigned int dangles_given,
			  const unsigned int noLP_given,
			  const double temp_arg,
			  const double betaScale_arg,
			  const int dangles_arg,
			  const int noLP_flag);
char *lila_db_from_pt(short int*);
int   lila_cmp_db(void *,void *);
int   lila_cmp_sse_lex(const void *,const void *);
int   lila_cmp_sse_en(const void *,const void *);
int   lila_cmp_sse_lexen(const void *,const void *);
void  lila_print_2se(gpointer,gpointer);
void  lila_free_cc(gpointer);

#endif
