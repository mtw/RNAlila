/*
  lila.h : common definitions for RNAlila
  Last changed Time-stamp: <2014-08-05 23:05:13 mtw>
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

/* functions */
void lila_initialize_vRNA (const char *seq);

/* variables */
vrna_fold_compound *vc;
paramT *P;
short int *pt,*s0,*s1;

#endif
