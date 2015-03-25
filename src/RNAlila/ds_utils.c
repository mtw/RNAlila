/*
  ds_utils.c: RNAlila data structure utility functions
  Last changed Time-stamp: <2014-09-08 17:08:20 mtw>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <lila.h>

/*
  Compare function for two strings.
*/
int
lila_cmp_db(const void *a,
	      const void *b)
{
return strcmp( (const char *)a,(const char *)b );	
}

/*
  Compare function for LilaDBE elements: compares (dot bracket)
  structures lexicographically.
*/
int
lila_cmp_dbe_lex(const void *a,
		   const void *b)
{
const LilaDBE *ma = a;
const LilaDBE *mb = b;

return strcmp(ma->structure,mb->structure);
}

/*
  Compare function for LilaDBE elements: compares energies.
*/
int
lila_cmp_dbe_en(const void *a,
		  const void *b)
{
const LilaDBE *ma = a;
  const LilaDBE *mb = b;

  if (ma->energy > mb->energy)
    return 1;
  else if (ma->energy < mb->energy)
    return -1;
  else
    return 0;
}

/*
  Compare function for LilaDBE elements: compares first by energy,
  then by lexicographical order of structure.
*/
int
lila_cmp_dbe_lexen(const void *a,
		   const void *b)
{
  const LilaDBE *ma = a;
  const LilaDBE *mb = b;

  int comp = ma->energy - mb->energy;

  if (comp < 0)
    return -1;
  
  if (comp > 0)
    return 1;
  
  comp = strcmp(ma->structure, mb->structure);
  return comp;
}

/*
  Return the lexicographically lowest member of a connected component
*/
GList *
lila_lexmin_dbe_glist(GList *c)
{
  fprintf(stderr,"[[lila_lexmin_cc]]::GList length is %i\n", (int)g_list_length(c));
  c = g_list_sort(c,(GCompareFunc)lila_cmp_dbe_lex);
  GList *first = g_list_first(c);
  return first;
}

/*
  Dump one LilaDBE element.
*/
void
lila_print_dbe(void *data,
	       void *user_data)
{
  LilaDBE *foo = (LilaDBE *)data;
  fprintf(stderr,"%s %6.2f ",foo->structure,foo->energy);
  if(user_data != NULL)
    fprintf(stderr,"%s ",(char*)user_data);
  fprintf(stderr,"\n");
}

/*
  Output entire GQueue of LilaDBE elements.
*/
void
lila_output_dbe_gqueue(GQueue *Q, const char *userdata)
{
  g_queue_foreach(Q,(GFunc)lila_print_dbe, (char*)userdata);
}

/*
  Output entire GList of LilaDBE elements.
*/
void
lila_output_dbe_glist(GList *L, const char* userdata)
{
  g_list_foreach(L,(GFunc)lila_print_dbe, (char*)userdata);
}

/*
  Deallocate GQueue of LilaDBE elements
 */
void
lila_dealloc_dbe_gqueue(GQueue *Q)
{
  g_queue_free_full(Q,(GDestroyNotify)lila_free_dbe);
}

/*
  Deallocate GList of LilaDBE elements
*/
void
lila_dealloc_dbe_glist(GList *L)
{
  g_list_free_full(L,(GDestroyNotify)lila_free_dbe);
}

/*
  Free one element of type LilaDBE
*/
void
lila_free_dbe(void *data)
{
  LilaDBE *foo = (LilaDBE*)data;
  free(foo->structure);
  free(foo);
}

/*
  Free one string
*/
void
lila_free_string(void *data)
{
  char *foo = (char*)data;
  free(foo);
}

/*
  Add a LilaDBE 'structure' element into a GHashTable
*/
void
lila_dbe_structure2ghashtable(void *data,
			      void *user_data)
{
  LilaDBE *foo = (LilaDBE *)data;
  GHashTable *S = (GHashTable *)user_data;
  char *form = strdup(foo->structure);
  if (g_hash_table_contains(S,form))
    fprintf(stderr, "%s ALREADY in hash ",form);
  else{
    g_hash_table_add(S, form);
    fprintf(stderr, "%s INSERTED to hash ",form);
  }
  fprintf(stderr, " [hashsize(S): %i]\n",(int)g_hash_table_size(S));
}
