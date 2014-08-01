/*
  lila.c: common routines for RNAlila
  Last changed Time-stamp: <2014-08-01 23:06:47 mtw>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <lila.h>

/* ==== */
void
lila_initialize_vRNA (const char *seq)
{
model_detailsT md;
set_model_details(&md); /* use current global model */
P = vrna_get_energy_contributions(md);
vc = vrna_get_fold_compound(seq, &md,VRNA_OPTION_MFE);
s0 = vrna_seq_encode_simple(seq,&(P->model_details));
s1 = vrna_seq_encode(seq,&(P->model_details));
}
