/*
  lila.h : common definitions for RNAlila
  Last changed Time-stamp: <2014-09-07 23:32:02 mtw>
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

typedef struct _lilass {  /* looks like this one is only used for */
  char *sequence;         /* passing input sequence & structure */
  char *structure;        /* throughout the library */
  int length;             /* TODO: consolidate me !!! */
} lilassT;                /* sequence, secondary structure and length */       

typedef struct _lila2se {
  char *structure;
  float energy;
} Lila2seT;                 /* secondary structure and energy */

typedef Lila2seT LilaDBE;   /* dot bracket & energy */

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
void lila_dump_pt(const short *);

/* functions defined in ds_utils.c */
int   lila_cmp_db(const void *,const void *);
int   lila_cmp_dbe_lex(const void *,const void *);
int   lila_cmp_dbe_en(const void *,const void *);
int   lila_cmp_dbe_lexen(const void *,const void *);
GList *lila_lexmin_dbe_glist(GList *);
void  lila_output_dbe_gqueue(GQueue *);
void  lila_output_dbe_glist(GList *);
void  lila_print_dbe(void *,void *);
void  lila_dealloc_dbe_gqueue(GQueue *);
void  lila_dealloc_dbe_glist(GList *);
void  lila_dbe_structure2ghashtable(void *,void *);

#endif
