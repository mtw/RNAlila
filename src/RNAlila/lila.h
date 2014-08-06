/*
  lila.h : common definitions for RNAlila
  Last changed Time-stamp: <2014-08-06 17:03:42 mtw>
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
typedef _vRNAcommon {
  double temperature;
  double betaScale;
  int dangles;
  int noLP;
} vRNAcommonT;

/* variables */
paramT *P;
model_detailsT md;
vrna_fold_compound *vc;
short int *pt,*s0,*s1;
vRNAcommonT vcd; /* ViennaRNA common (model) details */

/* functions */
void lila_initialize_vRNA (const char *seq);
void process_vcd_options(const short temp_given,
		    const double temp_arg,
		    const short dangles_given,
			 const int dangles_arg);
{
#endif
