/*  Last changed Time-stamp: <2014-09-02 17:56:57 mtw> */

#ifndef __MOVES__
#define __MOVES__

typedef struct move_str {
  int left;
  int right;
} move_str;

Lila2seT *cc;  /* connected component */

move_str   lila_random_move_pt(const char *,short int*);
move_str   lila_gradient_move_pt(const char *,short int*);
move_str   lila_adaptive_move_pt(const char *,short int*);
void       lila_apply_move_pt(short int *,const move_str);
int        lila_is_minimum_or_shoulder_pt(const char *,short int*);
move_str  *lila_all_adaptive_moves_pt(const char *, short int*,int*,int);
short int *lila_degenerate_cc_min_pt(const char *, short int *);

#endif
