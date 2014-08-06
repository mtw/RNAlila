/*
  lila.c: common routines for RNAlila
  Last changed Time-stamp: <2014-08-06 17:02:38 mtw>
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
  vrna_md_set_default(&md); /* set default model details */
  P = vrna_get_energy_contributions(md);
  vc = vrna_get_fold_compound(seq, &md,VRNA_OPTION_MFE);
  s0 = vrna_seq_encode_simple(seq,&(P->model_details));
  s1 = vrna_seq_encode(seq,&(P->model_details));
}

/* ==== */
/*  set temperature, dangles, noLP, betascale */
void
process_vcd_options(const short temp_given,
		    const double temp_arg,
		    const short dangles_given,
		    const int dangles_arg)
{
  /* temperature */
  if(temp_given){
    if( (vcd.temp = temp_arg) < -K0 ){
      fprintf(stderr, "Value of --temp must be > -273.15\n");
      exit (EXIT_FAILURE);
    }
    vrna_md_set_temperature(&md,vcdtemp);
  }
  /*etc ...*/
  /* dangles */

  /* noLP */

  /* betaScale */
}
