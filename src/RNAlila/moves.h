/*  Last changed Time-stamp: <2014-08-28 17:25:02 mtw> */

#ifndef __MOVES__
#define __MOVES__

typedef struct move_str {
  int left;
  int right;
} move_str;

move_str   lila_random_move_pt(const char *,short int*);
move_str   lila_gradient_move_pt(const char *,short int*);
move_str   lila_adaptive_move_pt(const char *,short int*);
void       lila_apply_move_pt(short int *,const move_str);
int        lila_is_minimum_pt(const char *,short int*);
move_str  *lila_all_adaptive_moves_pt(const char *, short int*,int*);
short int *lila_degenerate_cc_min_pt(const char *, short int *);

#endif
