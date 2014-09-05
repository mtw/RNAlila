/*  Last changed Time-stamp: <2014-09-06 00:36:24 mtw> */

#ifndef __MOVES__
#define __MOVES__

typedef struct move_str {
  int left;
  int right;
} move_str;

GQueue     *lila_generate_neighbors_pt(const char *,short int *);
int        lila_construct_moves(const char*, const short*, int , move_str **);
move_str   lila_random_move_pt(const char*,short int*);
move_str   lila_gradient_move_pt(const char*,short int*);
move_str   lila_adaptive_move_pt(const char*,short int*);
void       lila_apply_move_pt(short int*,const move_str);
move_str  *lila_all_adaptive_moves_pt(const char*,short int*,int*,int);

#endif
