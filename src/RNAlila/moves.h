/*  Last changed Time-stamp: <2014-09-04 14:23:13 mtw> */

#ifndef __MOVES__
#define __MOVES__

typedef struct move_str {
  int left;
  int right;
} move_str;

move_str   lila_random_move_pt(const char*,short int*);
move_str   lila_gradient_move_pt(const char*,short int*);
move_str   lila_adaptive_move_pt(const char*,short int*);
void       lila_apply_move_pt(short int*,const move_str);
int        lila_is_minimum_or_shoulder_pt(const char*,short int*);
move_str  *lila_all_adaptive_moves_pt(const char*,short int*,int*,int);
GList     *lila_get_cc_pt(const char*,short int*,int*);
void       lila_cleanup_cc(GList*);
void       lila_dump_cc(GList*);
GList     *lila_lexmin_cc(GList*);
void       lila_cc2hash(GHashTable*,GList*);
#endif
