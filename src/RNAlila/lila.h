/*
  lila.h : common definitions for RNAlila
  Last changed Time-stamp: <2014-08-07 00:30:03 mtw>
*/

#ifndef __RNA_LILA_H__
#define __RNA_LILA_H__

#include "config.h"
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
  char *startseq;
  char *startstruc;
  int length;
} lilassT;

/* variables */
paramT *P;
model_detailsT md;
vrna_fold_compound *vc;
short int *pt,*s0,*s1;
vRNAcommonT vcd;         /* ViennaRNA common (model) details */
lilassT lilass;          /* start sequence && start structure */

/* functions */
void lila_parse_seq_struc(FILE *fp);
void lila_ini_vRNA(const char *seq);
void lila_ini_vcd_options(void);
void lila_set_vcd_options(const unsigned int temp_given,
			  const unsigned int betaScale_given,
			  const unsigned int dangles_given,
			  const unsigned int noLP_given,
			  const double temp_arg,
			  const double betaScale_arg,
			  const int dangles_arg,
			  const int noLP_flag);

#endif
