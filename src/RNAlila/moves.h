/*  Last changed Time-stamp: <2014-08-01 14:42:31 mtw> */

#ifndef __MOVES__
#define __MOVES__

typedef struct move_str {
  int left;
  int right;
} move_str;

move_str lila_random_move_pt(const char *,const short int*);
move_str lila_gradient_move_pt(const char *,const short int*);
move_str lila_adaptive_move_pt(const char *,const short int*);
void     lila_apply_move_pt(short int *,const move_str);

#endif
